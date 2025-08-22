#!/bin/bash
#SBATCH --job-name=kmer14
#SBATCH --output=kmer14_%j.out
#SBATCH --error=kmer14_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tn2220@nyu.edu

set -euo pipefail

#####################
# EDIT THESE INPUTS #
#####################
REF="/scratch/tn2220/Mitonuclear/REF/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
PARENT_VCF="parents_diff.vcf.gz"     # vcf with BOTH parents; biallelic sites where they differ
PAT="ECA1493"                        # paternal sample ID (must match VCF sample name)
MAT="ECA1229"                        # maternal sample ID (must match VCF sample name)

# Your 6 sample BAMs (will be converted to FASTQ then counted)
BAMS=(ECA1229.bam ECA1493.bam QG4448.bam QG4454.bam QG4499.bam QG4505.bam)

# Jellyfish hash sizes
HASH_REF="500M"   # ref k-mers
HASH_SAMP="4G"    # per-sample reads

#########################
# MODULES (ADJUST VERS) #
#########################
module purge
module load samtools/intel/1.14
module load bcftools/intel/1.14
module load bedtools/intel/2.29.2
module load jellyfish/2.3.0
# prefer pigz; fall back to zcat if missing
if module avail pigz 2>/dev/null | grep -q pigz; then module load pigz; fi

#########################
# FOLDERS & SANITY CHECK#
#########################
mkdir -p work kmers counts plots
test -s "$REF"        || { echo "ERROR: REF not found: $REF"; exit 1; }
test -s "$PARENT_VCF" || { echo "ERROR: PARENT_VCF not found: $PARENT_VCF"; exit 1; }

echo "[INFO] Indexing reference if needed ..."
samtools faidx "$REF"

#############################################
# STEP 0a: QUICK DIAGNOSTICS ON NAMES/BUILD #
#############################################
echo "[CHECK] Comparing contig names between REF and parental VCF ..."
awk '{print $1}' "${REF}.fai" | sort > work/ref.contigs.txt
bcftools index -f "$PARENT_VCF" >/dev/null 2>&1 || true
bcftools view -h "$PARENT_VCF" | grep -E '^##contig=<ID=' | sed -E 's/.*ID=([^,>]+).*/\1/' | sort > work/vcf.contigs.txt || true

if ! comm -12 work/ref.contigs.txt work/vcf.contigs.txt | head -1 | grep -q . ; then
  echo "[WARN] No shared contig names between REF and VCF. Likely 'chr' vs no 'chr' or build mismatch."
  echo "[INFO] Attempting a simple rename guess (add/remove 'chr' prefix)."
  # build a rename map guessing 'I' <-> 'chrI', etc.
  awk '{print $1"\tchr"$1}' "${REF}.fai" > work/ref_to_chr.map
  awk '{print "chr"$1"\t"$1}' "${REF}.fai" > work/chr_to_ref.map

  if grep -q '^chr' work/vcf.contigs.txt; then
    echo "[INFO] VCF has 'chr*', REF has no 'chr'. Renaming VCF -> REF names."
    bcftools annotate --rename-chrs work/chr_to_ref.map -Oz -o work/parents_diff.renamed.vcf.gz "$PARENT_VCF"
    tabix -f -p vcf work/parents_diff.renamed.vcf.gz
    PARENT_VCF="work/parents_diff.renamed.vcf.gz"
  else
    echo "[INFO] VCF has no 'chr', REF maybe has 'chr'. Renaming VCF -> chr*."
    bcftools annotate --rename-chrs work/ref_to_chr.map -Oz -o work/parents_diff.renamed.vcf.gz "$PARENT_VCF" || true
    if [ -s work/parents_diff.renamed.vcf.gz ]; then
      tabix -f -p vcf work/parents_diff.renamed.vcf.gz
      PARENT_VCF="work/parents_diff.renamed.vcf.gz"
    fi
  fi
fi

echo "[CHECK] Ensuring PAT/MAT sample IDs exist in VCF ..."
if ! bcftools query -l "$PARENT_VCF" | grep -qx "$PAT"; then
  echo "ERROR: Sample '$PAT' not found in $PARENT_VCF"; exit 1;
fi
if ! bcftools query -l "$PARENT_VCF" | grep -qx "$MAT"; then
  echo "ERROR: Sample '$MAT' not found in $PARENT_VCF"; exit 1;
fi

#############################################
# STEP 1: Parental SNPs -> 25-mer quartets  #
#############################################
echo "[STEP 1] Extract parental SNPs and build 25-mer quartets ..."
# SNPs (biallelic) CHROM POS REF ALT
bcftools view -m2 -M2 -v snps "$PARENT_VCF" \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > work/parental_snps.tsv

# if empty, stop early
if [ ! -s work/parental_snps.tsv ]; then
  echo "ERROR: No SNPs extracted from $PARENT_VCF after contig-name harmonization."; exit 1;
fi

# 25bp windows (BED: 0-based, half-open): [POS-13, POS+12)
awk -v OFS='\t' '{print $1, $2-13, $2+12, $1":"$2}' work/parental_snps.tsv > work/snps.25bp.bed

# fetch sequences name<TAB>sequence
bedtools getfasta -fi "$REF" -bed work/snps.25bp.bed -name -tab > work/snps.25bp.tab

# build quartets: id chrom pos center kA kC kG kT
awk -v OFS='\t' '
  {
    split($1,a,":"); chrom=a[1]; pos=a[2];
    seq=toupper($2);
    if(length(seq)!=25) next;
    pre=substr(seq,1,12); cen=substr(seq,13,1); post=substr(seq,14,12);
    print chrom":"pos, chrom, pos, cen, pre "A" post, pre "C" post, pre "G" post, pre "T" post
  }' work/snps.25bp.tab \
| awk 'BEGIN{print "id\tchrom\tpos\tcenter\tkA\tkC\tkG\tkT"}1' > kmers/quartets.tsv

# flat list of all quartet kmers
awk 'NR>1{print $5"\n"$6"\n"$7"\n"$8}' kmers/quartets.tsv > kmers/quartets.list

#############################################################
# STEP 2: Reference single-copy & only-one-match filter     #
#############################################################
echo "[STEP 2] Count reference 25-mers and filter to single-copy unique sites ..."
jellyfish count -m 25 -s "$HASH_REF" -t "$SLURM_CPUS_PER_TASK" -C -o kmers/ref25.jf "$REF"
jellyfish dump -c kmers/ref25.jf > kmers/ref25.dump

# attach ref counts for A/C/G/T
awk 'NR==FNR{cnt[$1]=$2; next}
     NR>1{
       cA=(($5 in cnt)?cnt[$5]:0); cC=(($6 in cnt)?cnt[$6]:0);
       cG=(($7 in cnt)?cnt[$7]:0); cT=(($8 in cnt)?cnt[$8]:0);
       print $0"\t"cA"\t"cC"\t"cG"\t"cT
     }' kmers/ref25.dump kmers/quartets.tsv > kmers/quartets.refcnt.tsv

# keep where exactly ONE quartet member present in reference AND that count==1
awk 'BEGIN{OFS="\t"}
     NR==1{print $0,"\tkeep"; next}
     {
       a=$(NF-3); c=$(NF-2); g=$(NF-1); t=$(NF);
       n=(a>0)+(c>0)+(g>0)+(t>0);
       single=(a==1)+(c==1)+(g==1)+(t==1);
       keep=(n==1 && single==1)?"YES":"NO";
       print $0, keep
     }' kmers/quartets.refcnt.tsv | awk '$NF=="YES"{sub(/\tYES$/,""); print}' > kmers/quartets.filtered.tsv

sites_step2=$(($(wc -l < kmers/quartets.filtered.tsv)-1))
echo "[STEP 2] Retained $sites_step2 candidate sites after reference uniqueness filter."
if [ "$sites_step2" -lt 1 ]; then
  echo "ERROR: No sites passed the reference uniqueness filter. Aborting."; exit 1;
fi

##########################################################
# STEP 3: Assign PAT/MAT bases from parental VCF         #
#        (robust: use REF/ALT + GT -> bases; no %TGT)    #
##########################################################
echo "[STEP 3] Pull PAT/MAT bases from parental VCF and build diagnostic_sites.tsv ..."

# Pull per-site REF, ALT and each parent's GT (numeric alleles)
# Keep only biallelic SNPs in the query (m2/M2/v snps) to ensure ALT is a single base.
bcftools view -m2 -M2 -v snps "$PARENT_VCF" \
| bcftools query -s "$PAT","$MAT" -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' \
> work/parent_refalt_gt.tsv

# Translate GT -> base for each parent, drop het/missing/ambiguous; keep PAT_BASE != MAT_BASE and both A/C/G/T
# Output: chrom pos PAT_BASE MAT_BASE
awk -v OFS='\t' -v pat="$PAT" -v mat="$MAT" '
  function gt2base(gt, ref, alt,   a,b) {
    # Return a single base for homozygous (0/0 or 1/1). Return "" for het/missing/complex.
    if (gt == "0/0" || gt == "0|0") return ref;
    if (gt == "1/1" || gt == "1|1") return alt;
    return ""; # het (0/1, 1/0), missing (./.), or multi-allelic: skip
  }
  {
    chrom=$1; pos=$2; ref=toupper($3); alt=toupper($4);
    # parse sample=GT pairs (two samples expected)
    patgt = ""; matgt = "";
    for (i=5; i<=NF; i++) {
      split($i, p, "="); if (p[1]==pat) patgt=p[2]; else if (p[1]==mat) matgt=p[2];
    }
    if (patgt=="" || matgt=="") next;
    patb = toupper(gt2base(patgt, ref, alt));
    matb = toupper(gt2base(matgt, ref, alt));
    if (patb ~ /^[ACGT]$/ && matb ~ /^[ACGT]$/ && patb != matb)
      print chrom, pos, patb, matb;
  }
' work/parent_refalt_gt.tsv > work/parent_bases.clean.tsv

# Join to quartets.filtered.tsv -> diagnostic_sites.tsv
# quartets.filtered.tsv cols: id chrom pos center kA kC kG kT refA refC refG refT
awk -v OFS='\t' '
  BEGIN{
    while((getline < "work/parent_bases.clean.tsv")>0){
      key=$1":"$2; PB[key]=$3"\t"$4;
    }
    print "id","chrom","pos","PAT_BASE","MAT_BASE","kA","kC","kG","kT";
  }
  FNR==1 {next}  # skip header
  {
    id=$1; chrom=$2; pos=$3; key=chrom":"pos;
    if(!(key in PB)) next;
    split(PB[key],b,"\t"); patb=b[1]; matb=b[2];
    center=toupper($4); if(center !~ /^[ACGT]$/) next;
    print id,chrom,pos,patb,matb,$5,$6,$7,$8;
  }
' kmers/quartets.filtered.tsv > kmers/diagnostic_sites.tsv

sites_final=$(($(wc -l < kmers/diagnostic_sites.tsv)-1))
echo "[STEP 3] diagnostic_sites.tsv ready with $sites_final sites."
if [ "$sites_final" -lt 1 ]; then
  echo "ERROR: 0 diagnostic sites after joining with parental bases."
  echo "DEBUG HINTS:"
  echo "  - Inspect a few rows with GTs:"
  echo "      head -5 work/parent_refalt_gt.tsv"
  echo "  - Confirm parents are homozygous and different at some sites:"
  echo "      awk '\''$3!=$4{print}'\'' work/parent_bases.clean.tsv | head"
  echo "  - Check a key overlap:"
  echo "      head -3 kmers/quartets.filtered.tsv | awk '\''NR>1{print $2\":\"$3}'\''"
  echo "      head -3 work/parent_bases.clean.tsv | awk '\''{print $1\":\"$2}'\''"
  exit 1
fi

# k-mers to query in samples
awk 'NR>1{print $6"\n"$7"\n"$8"\n"$9}' kmers/diagnostic_sites.tsv > kmers/diag_kmers.list

##########################################################
# STEP 4: Sample k-mer counting & per-site PAT/MAT/OTHER #
##########################################################
# clean + dedup k-mer list
tr -d '\r' < kmers/diag_kmers.list | awk '!seen[$0]++' > kmers/diag_kmers.list.tmp && mv kmers/diag_kmers.list.tmp kmers/diag_kmers.list

# choose decompressor
if command -v pigz >/dev/null 2>&1; then DECOMP="pigz -dc"; else DECOMP="zcat -f"; fi

count_with_fifos() {
  local s="$1"
  local r1="${s}_1.fq.gz"
  local r2="${s}_2.fq.gz"
  local jf="kmers/${s}.jf"

  # if jf exists and is non-empty/valid, skip
  if [[ -s "$jf" ]]; then
    if jellyfish stats "$jf" > "kmers/${s}.stats.txt" 2>&1 && ! grep -q "Total = 0" "kmers/${s}.stats.txt"; then
      echo "[INFO] Using existing $jf"
      return 0
    fi
    echo "[WARN] Rebuilding $jf because stats failed or zero total."
  fi

  # ensure FASTQs exist
  if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    echo "[INFO] Converting BAM -> FASTQ for $s"
    samtools fastq -@ "$SLURM_CPUS_PER_TASK" -1 "$r1" -2 "$r2" "${s}.bam"
  fi

  # make FIFOs
  tmpd=$(mktemp -d)
  f1="$tmpd/${s}_R1.fq"; f2="$tmpd/${s}_R2.fq"
  mkfifo "$f1" "$f2"

  # background decompress into FIFOs
  (set -o pipefail; $DECOMP "$r1" > "$f1") &
  pid1=$!
  (set -o pipefail; $DECOMP "$r2" > "$f2") &
  pid2=$!

  # count reading from FIFOs (no stdin, no gz dependency)
  jellyfish count \
    -m 25 \
    -s "$HASH_SAMP" \
    -t "$SLURM_CPUS_PER_TASK" \
    -C \
    -o "$jf" \
    "$f1" "$f2"

  # wait and clean
  wait $pid1
  wait $pid2
  rm -rf "$tmpd"

  # verify DB
  if ! jellyfish stats "$jf" > "kmers/${s}.stats.txt" 2>&1; then
    echo "[ERROR] jellyfish stats failed for $jf"; exit 1
  fi
  if grep -q "Total = 0" "kmers/${s}.stats.txt"; then
    echo "[ERROR] $jf contains 0 kmers — check FASTQs / HASH_SAMP / memory"; exit 1
  fi
}

for bam in "${BAMS[@]}"; do
  s=$(basename "$bam" .bam)

  # Count via FIFOs (or reuse valid DB)
  echo "[INFO] Counting kmers for $s ..."
  count_with_fifos "$s"

  echo "[INFO] Dumping and filtering counts for $s ..."
  jellyfish dump -c "kmers/${s}.jf" > "kmers/${s}.dump"

  # filter to diagnostic kmers (write query-like output: "<kmer> <count>")
  awk 'NR==FNR{wanted[$1]=1; next} ($1 in wanted){print $1, $2}' \
    kmers/diag_kmers.list "kmers/${s}.dump" > "kmers/${s}.rawq"

  if [[ ! -s kmers/${s}.rawq ]]; then
    echo "[ERROR] No matches for diagnostic kmers in kmers/${s}.dump"; exit 1
  fi

  echo "[INFO] Summarizing pat/mat/other for $s ..."
  awk -v OFS='\t' -v RAW="kmers/${s}.rawq" '
    BEGIN{
      while((getline<RAW)>0){ split($0,a," "); Q[a[1]]=a[2]+0 }
      print "sample","chrom","pos","pat_reads","mat_reads","other_reads","total"
    }
    NR==1{next}
    {
      chrom=$2; pos=$3; PATB=$4; MATB=$5; kA=$6; kC=$7; kG=$8; kT=$9;
      cA=Q[kA]+0; cC=Q[kC]+0; cG=Q[kG]+0; cT=Q[kT]+0;
      pr=(PATB=="A"?cA:(PATB=="C"?cC:(PATB=="G"?cG:cT)));
      mr=(MATB=="A"?cA:(MATB=="C"?cC:(MATB=="G"?cG:cT)));
      other=(cA+cC+cG+cT) - pr - mr; if(other<0) other=0;
      tot=pr+mr+other;
      print "'$s'",chrom,pos,pr,mr,other,tot
    }
  ' kmers/diagnostic_sites.tsv > "counts/${s}.patmat.tsv"
done

# combine all samples
head -n1 "counts/$(basename "${BAMS[0]}" .bam).patmat.tsv" > counts/all.patmat.tsv
for bam in "${BAMS[@]}"; do
  s=$(basename "$bam" .bam)
  tail -n +2 "counts/${s}.patmat.tsv" >> counts/all.patmat.tsv
done

echo "[DONE] Steps 1–4 complete."
echo "Outputs:"
echo "  - kmers/quartets.tsv"
echo "  - kmers/quartets.filtered.tsv  (single-copy unique)"
echo "  - kmers/diagnostic_sites.tsv   (with PAT/MAT bases)"
echo "  - kmers/diag_kmers.list"
echo "  - counts/*.patmat.tsv and counts/all.patmat.tsv"
echo ""
echo "Next: plot counts/all.patmat.tsv with your R script."

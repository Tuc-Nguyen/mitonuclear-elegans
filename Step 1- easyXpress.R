library(data.table)
library(easyXpress)
require(stringi)
require(stringr)
require(plyr)
require(ggplot2)
require(tibble)

ori.wd = getwd()

datafiles = list.files(path ="./cp_data", pattern = ".RData", full.names = FALSE, recursive = FALSE)

for (datafile in datafiles){
  raw <- easyXpress::readXpress(filedir = ori.wd, rdafile = datafile, design = FALSE)
  raw = subset(raw, Worm_Length > 100)
  consolidate = vector("list",nrow(raw) )
  c <- 1
  for (i in c(1:nrow(raw))){
    tempdf= raw[i,]
    tempdf$Metadata_Group = NULL
    condition = (unlist(strsplit(tempdf$Metadata_Plate,"_")))[1]
    well = (unlist(strsplit(tempdf$FileName_RawBF,"_")))[2]
    tempdf$Metadata_Well = (unlist(strsplit(well,"[.]")))[1]
    tempdf$Metadata_Plate = condition
    consolidate[[c]]<-tempdf
    c <- c + 1
  }
  raw = as.data.frame(data.table::rbindlist(consolidate))
  raw$model = ifelse(raw$model == "L1_N2_HB101_100w_NonOverlappingWorms.model.outputs","L1", 
                     ifelse(raw$model == "L2L3_N2_HB101_100w_NonOverlappingWorms.model.outputs","L2L3",
                            ifelse(raw$model == "L4_N2_HB101_100w_NonOverlappingWorms.model.outputs","L4",
                                   ifelse(raw$model == "MDHD_NonOverlappingWorms.model.outputs", "MDHD", NA))))
  model_selected <- easyXpress::modelSelection(raw)
  edge_flagged <- easyXpress::edgeFlag(model_selected, radius=825, center_x=1024, center_y=1024)
  raw_flagged <- easyXpress::setFlags(edge_flagged, cluster_flag = TRUE, well_edge_flag = TRUE)
  
  processed <- easyXpress::process(raw_flagged, Metadata_Plate, Metadata_Well)
  name <- unlist(strsplit(datafile,"[.]"))[1]
  assign(name, processed)
}


raw.data <- rbind(`HeavyMetals-1A`[[1]], `HeavyMetals-1B`[[1]],
                  `HeavyMetals-2A`[[1]], `HeavyMetals-2B`[[1]],
                  `HeavyMetals-3A`[[1]], `HeavyMetals-3B`[[1]],
                  `HeavyMetals-4A`[[1]], `HeavyMetals-4B`[[1]])
processed.data <- rbind(`HeavyMetals-1A`[[2]], `HeavyMetals-1B`[[2]],
                        `HeavyMetals-2A`[[2]], `HeavyMetals-2B`[[2]],
                        `HeavyMetals-3A`[[2]], `HeavyMetals-3B`[[2]],
                        `HeavyMetals-4A`[[2]], `HeavyMetals-4B`[[2]])
raw.summary <- rbind(`HeavyMetals-1A`[[3]], `HeavyMetals-1B`[[3]],
                     `HeavyMetals-2A`[[3]], `HeavyMetals-2B`[[3]],
                     `HeavyMetals-3A`[[3]], `HeavyMetals-3B`[[3]],
                     `HeavyMetals-4A`[[3]], `HeavyMetals-4B`[[3]])
processed.summary <- rbind(`HeavyMetals-1A`[[4]], `HeavyMetals-1B`[[4]],
                           `HeavyMetals-2A`[[4]], `HeavyMetals-2B`[[4]],
                           `HeavyMetals-3A`[[4]], `HeavyMetals-3B`[[4]],
                           `HeavyMetals-4A`[[4]], `HeavyMetals-4B`[[4]])

df = join(processed.data,Plate_Design)
write.table(df, "Processed.txt", row.names = FALSE)

df = plyr::join(processed.summary,Plate_Design)
write.table(df, "Processed Summary.txt", row.names = FALSE)

df = join(raw.data,Plate_Design)
write.table(df, "Raw.txt", row.names = FALSE)

df = join(raw.summary,Plate_Design)
write.table(df, "Raw Summary.txt", row.names = FALSE)




rm(list = ls())
rstudioapi::restartSession()



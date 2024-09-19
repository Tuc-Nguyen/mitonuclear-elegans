require(ggplot2)
require(reshape2)
require(lme4)
require(plyr)
require(Rmisc)
require(ggbeeswarm)
require(gridExtra)
require(grid)
require(tidyverse)



df = read.table("FinalDataframe.txt", header = T, na.strings = NA)


require(lme4)
treatment=factor(df$Condition)
nuclear=factor(df$Nuclear)
mito=factor(df$Mito)
batch = factor(df$Batch)

FM=lmer(df$PHENO~(1|nuclear)+(1|mito)+treatment+
          (1|nuclear:treatment)+(1|mito:treatment)+(1|nuclear:mito)+
          (1|nuclear:mito:treatment))
summary(FM)

#no Nuclear
FM1=lmer(df$PHENO~(1|mito)+treatment+
           (1|nuclear:treatment)+(1|mito:treatment)+(1|nuclear:mito)+
           (1|nuclear:mito:treatment))
a = anova(FM,FM1)
chisq1 = a$Chisq[2]
degreefreedom1 = a$Df[2]
p1 = a$`Pr(>Chisq)`[2]

#no Mito
FM2=lmer(df$PHENO~(1|nuclear)+treatment+
           (1|nuclear:treatment)+(1|mito:treatment)+(1|nuclear:mito)+
           (1|nuclear:mito:treatment))
a = anova(FM,FM2)
chisq2 = a$Chisq[2]
degreefreedom2 = a$Df[2]
p2 = a$`Pr(>Chisq)`[2]

#no Treatment
FM3=lmer(df$PHENO~(1|nuclear)+(1|mito)+
           (1|nuclear:treatment)+(1|mito:treatment)+(1|nuclear:mito)+
           (1|nuclear:mito:treatment))
a = anova(FM,FM3)
chisq3 = a$Chisq[2]
degreefreedom3 = a$Df[2]
p3 = a$`Pr(>Chisq)`[2]

#no Nuclear-Treatment
FM4=lmer(df$PHENO~(1|nuclear)+(1|mito)+treatment+
           (1|mito:treatment)+(1|nuclear:mito)+
           (1|nuclear:mito:treatment))
a=anova(FM,FM4)
chisq4 = a$Chisq[2]
degreefreedom4 = a$Df[2]
p4 = a$`Pr(>Chisq)`[2]

#no Mito-treatment
FM5=lmer(df$PHENO~(1|nuclear)+(1|mito)+treatment+
           (1|nuclear:treatment)+(1|nuclear:mito)+
           (1|nuclear:mito:treatment))
a=anova(FM,FM5)
chisq5 = a$Chisq[2]
degreefreedom5 = a$Df[2]
p5 = a$`Pr(>Chisq)`[2]

#no Mitonuclear
FM6=lmer(df$PHENO~(1|nuclear)+(1|mito)+treatment+
           (1|nuclear:treatment)+(1|mito:treatment)+
           (1|nuclear:mito:treatment))
a=anova(FM,FM6)
chisq6 = a$Chisq[2]
degreefreedom6 = a$Df[2]
p6 = a$`Pr(>Chisq)`[2]

#no GXGXE
FM7=lmer(df$PHENO~(1|nuclear)+(1|mito)+treatment+
           (1|nuclear:treatment)+(1|mito:treatment)+(1|nuclear:mito))
a=anova(FM,FM7)
chisq7 = a$Chisq[2]
degreefreedom7 = a$Df[2]
p7 = a$`Pr(>Chisq)`[2]



tables9 = data.frame()
for(condition in conditions){
  subcondition = subset(df, Condition == condition) 
  nuclear=factor(subcondition$Nuclear)
  mito=factor(subcondition$Mito)
  
  FM=lmer(subcondition$PHENO~(1|nuclear)+(1|mito)+(1|nuclear:mito))
  #no Nuclear
  FM1=lmer(subcondition$PHENO~(1|mito)+(1|nuclear:mito))
  a = anova(FM,FM1)
  p1 = a$`Pr(>Chisq)`[2]
  d1 = a$Df[2]
  print(d1)
  #no Mito
  FM2=lmer(subcondition$PHENO~(1|nuclear)+(1|nuclear:mito))
  b = anova(FM,FM2)
  p2 = b$`Pr(>Chisq)`[2]
  d2 = b$Df[2]
  print(d2)
  #no Mitonuclear
  FM3=lmer(subcondition$PHENO~(1|nuclear)+(1|mito))
  c=anova(FM,FM3)
  p3 = c$`Pr(>Chisq)`[2]
  d3 = c$Df[2]
  print(d3)
  tempdf = data.frame(Condition = condition, Nuclear = p1, Mito = p2, Mitonuclear = p3)
  tables7 = rbind(tables7, tempdf)
}


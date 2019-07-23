setwd("/Users/lindsayvm/Documents/Avillach_Lab/projects/PheWAS")
library(dplyr)
library(data.table)
library(tidyr)

##### CREATE DATAFRAME

source("src/functions.R")
# Go to https://copdgene.hms.harvard.edu/transmart/login/auth, utilities, user profile, IRCT Token
token = "eyJraWQiOiJSa05CUWpFNU9VTkVOelkzTmpJd04wVkNNVGd3TmpFM01EVXdSVEpETVVaRk5EZzROa0ZFUmciLCJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJhdWQiOiJHSjhmWWlOZ3VsbEhyZTZOYWptV25QblQ1NTRlNW4xdSIsInN1YiI6InNhbWxwfExpbmRzYXlfTGVla0BobXMuaGFydmFyZC5lZHUiLCJpc3MiOiJodHRwczovL2F2aWxsYWNobGFiLmF1dGgwLmNvbS8iLCJleHAiOjE1NjExMDU2MzcsImlhdCI6MTU2MTA2OTYzNywiZW1haWwiOiJMaW5kc2F5X0xlZWtAaG1zLmhhcnZhcmQuZWR1In0.UrKOcV4_Bwv0XIaa5xO7bf4PtfUkaW5IXhmIEimNXVw"
IRCT_REST_BASE_URL = "https://copdgene.hms.harvard.edu"

# Phenotypes of interest for all patients in freeze 5 (n = 10.000)
pheno.df1 = reloadDfFormFile("data/ids.txt", token, IRCT_REST_BASE_URL)
# make a raw copy
raw.df = pheno.df1

rm(list=setdiff(ls(), c("raw.df", "pheno.names1")))
pheno.df1 = raw.df



#summary(as.factor(pheno.df1$Affection_status))
# COPD cases and controls as in Parket et al 2019 and Yi et al 2018
#p heno.df1 = pheno.df1[pheno.df1$Affection_status == "Control" | pheno.df1$Affection_status =="Case", ] 
# only COPD cases
pheno.df1 = pheno.df1[pheno.df1$Affection_status =="Case", ] 
summary(as.factor(pheno.df1$Affection_status))

# Remove phenotypes
pheno.df1 = subset(pheno.df1, select = -c(Lung_cancer,              # lung cancer values are off????
                                          COPD, 
                                          Affection_status,
                                          Chest_injuries,
                                          Chest_operations,
                                          Other_chest_illnesses))  # all samples are positive for COPD or affection_status

dim(pheno.df1)
#Change phenotypic scores to 1s, 0s and NAs 
pheno.df1[pheno.df1 == "Yes"] = 1
pheno.df1[pheno.df1 == "No"]  = 0 
pheno.df1[pheno.df1 == "Do not know"] = NA
#Change phenotypic scores from character to numeric, keep patient ID as character
pheno.df1[ ,-1] = lapply(pheno.df1[ ,-1], function(x) {
  if(is.character(x)) as.numeric((x))
})
#Check class:
#str(pheno.df1)


# GT.df has GT calls of selected variants for samples of interest: freeze 4 (n = 2000)
GT.df = read.csv("Parker/Parker_GThits.csv",       
                 stringsAsFactors = FALSE, 
                 header = TRUE,  
                 sep = ",")          
length(unique(GT.df$dbGaP_Subject_ID))

#data/NHW_top6.csv         
#data/NHW_Sign7e-06.csv     
#data/combi_pval5e-06.csv   
#data/combi_noprunepval5e-06.csv   #
#data/combi_FEV1FVC_pval5e-06.csv  
#data/combi_FEV1_pval5e-06.csv     

#data/Parker_GThits.csv
#data/Yi_GThits.csv

#head(GT.df)
#dim(GT.df)
#str(GT.df)

#CLUMPING

# Assign 1 if GT call is [1,0] or [0,1] and 0 if [0,0]
GT.df$GT.alleles[GT.df$GT.alleles != "[0, 0]"] = 1
GT.df$GT.alleles[GT.df$GT.alleles == "[0, 0]"] = 0
GT.df$locus = paste0(GT.df$locus.contig, ":", GT.df$locus.position)

# values are duplicated, and we need the second set. But if you do duplicated(), it will remove the second one
# and if you want the opposite, you will remove all the other data. 
GT.df = GT.df[!GT.df$locus == "chr2:228705203", ]

dim(GT.df)
##BOOSDOENER "chr5:157502070"
locus.names = unique(GT.df$locus)
for(l in 1:length(locus.names)){
  tmp = GT.df[GT.df$locus == locus.names[l], ]
  if(sum(as.numeric(tmp$GT.alleles)) < 10){
    GT.df = GT.df[!GT.df$locus == locus.names[l], ]
  }
}
dim(GT.df)

# MERGE data
# pheno.df and GT.df are both from freeze 5, but pheno.df has more samples
# therefore remove those in GT.df that are not present in GT.df
final.df1 = pheno.df1[pheno.df1$Patient.Id %in% GT.df$dbGaP_Subject_ID, ]
# Next, pheno.df1 is filtered on cases, and GT.df not yet.
GT.df = GT.df[GT.df$dbGaP_Subject_ID %in% final.df1$Patient.Id, ]

#PHEWAS PER VARIANT
pheno.names1 = names(pheno.df1[ ,c(2:ncol(pheno.df1))]) #all phenotypes
locus.names = unique(GT.df$locus) # all loci
signlv = 0.05
for(l in 1:length(locus.names)){
  # select for variant
  var.df = GT.df[GT.df$locus == locus.names[l], ]
  # remove duplicates
  var.df = var.df[!duplicated(var.df$dbGaP_Subject_ID), ]
  #annotate GT calls of selected variant
  final.df1$GTalleles = var.df$GT.alleles    #GTalleles is character
  for(i in 1:length(pheno.names1)){
    #general linear regression on phenotype of interest(POI)
    POI.glm1 = glm(formula = final.df1[ ,pheno.names1[i]]~GTalleles, family = binomial(), data = final.df1, na.action = na.omit)
    #summary(POI.glm1)
    if(summary(POI.glm1)$coefficients[2,4] < signlv){    #/n_variants
      print(c(locus.names[l] , pheno.names1[i], summary(POI.glm1)$coefficients[2,4]))
    }
  }
}

#[1] "chr6:32183666"        "Pneumothorax"         "0.000219912890125386"

n_va = length(locus.names)
n_ph = length(pheno.names1)
print(c("va = ", signlv/n_va, "ph =", signlv/n_ph, "va&ph =", signlv/n_ph/n_va))


######################### FOR SINGLE SNP

l = 10 #"chr6:32183666" 
var.df = GT.df[GT.df$locus == locus.names[l], ]
final.df1$GTalleles = var.df$GT.alleles
phewasOutput.df1 = as.data.frame(matrix(ncol = 7)) #df to store output in
for(i in 1:length(pheno.names1)){
  #general linear regression on phenotype of interest(POI)
  POI.glm1 = glm(formula = final.df1[ ,pheno.names1[i]]~GTalleles, family = binomial(), data = final.df1, na.action = na.omit)
  #summary(POI.glm1)
  
  # confidence interval
  ci1 = exp(summary(POI.glm1)$coefficients["GTalleles1", 1] + qnorm(c(0.025, 0.975)) * summary(POI.glm1)$coefficients["GTalleles1", 2])
  
  # caco.df1: for COPD cases, which patient has pheno (case (1)), and which do not (ctrl (0)) 
  caco.df1        = final.df1[ ,c((pheno.names1)[i], "Patient.Id", "GTalleles")]
  caco.df1        = na.omit(caco.df1)
  # Count number of cases with disease (CaseDisease1), without etc etc 
  CaseDisease1   = length(unique(caco.df1[caco.df1[ ,1] == 1 & caco.df1$GTalleles == "1", "Patient.Id"])) # pheno yes; variant yes
  CaseNoDisease1 = length(unique(caco.df1[caco.df1[ ,1] == 0 & caco.df1$GTalleles == "1", "Patient.Id"])) # pheno no; variant yes
  CtrlDisease1   = length(unique(caco.df1[caco.df1[ ,1] == 1 & caco.df1$GTalleles == "0", "Patient.Id"])) # pheno yes; variant no
  CtrlNoDisease1 = length(unique(caco.df1[caco.df1[ ,1] == 0 & caco.df1$GTalleles == "0", "Patient.Id"])) # pheno no; variant no
  newRow1        = c(pheno.names1[i], 
                     summary(POI.glm1)$coefficients[2], 
                     exp(summary(POI.glm1)$coefficients[2]), 
                     paste0("[", round(ci1[1][1],3), ", ", round(ci1[2][1],3), "]"),
                     summary(POI.glm1)$coefficients[2,4], 
                     paste0(CaseDisease1   + CtrlDisease1,   "(", CaseDisease1,  "/", CtrlDisease1,  ")"), 
                     paste0(CaseNoDisease1 + CtrlNoDisease1, "(", CaseNoDisease1,"/", CtrlNoDisease1,")"))
  phewasOutput.df1 = rbind(newRow1, phewasOutput.df1)
}
colnames(phewasOutput.df1) = c("Phenotype", "Coefficient", "OR", "Confidence_interval", "Pvalue", "Pheno_present(variant/novariant)", "Pheno_absent(variant/novariant)")
#Remove all rows with only NAs
phewasOutput.df1 = phewasOutput.df1[complete.cases(phewasOutput.df1), ]
phewasOutput.df1$adjustPvalue <- p.adjust( as.numeric(phewasOutput.df1$Pvalue), method = "BH")
phewasOutput.df1

# Significant phenotypes based onf  adj
#adjPvalue_Sign.df1 = phewasOutput.df1[as.numeric(phewasOutput.df1$adjPvalue) < 0.05, ]
#OR_Sign.df1   = adjPvalue_Sign.df1[as.numeric(adjPvalue_Sign.df1$OR) > 2, ] 

#OR_Sign.df1   = phewasOutput.df1[as.numeric(phewasOutput.df1$OR) > 1.5, ] 
#OR_Sign.df1

l = 2

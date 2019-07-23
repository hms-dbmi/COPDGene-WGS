setwd("/Users/lindsayvm/Documents/Avillach_Lab/scripts/PheWAS")
library(dplyr)
library(data.table)
library(tidyr)

##### CREATE DATAFRAME

source("src/functions.R")
# Go to https://copdgene.hms.harvard.edu/transmart/login/auth, utilities, user profile, IRCT Token
token = "eyJraWQiOiJSa05CUWpFNU9VTkVOelkzTmpJd04wVkNNVGd3TmpFM01EVXdSVEpETVVaRk5EZzROa0ZFUmciLCJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJhdWQiOiJHSjhmWWlOZ3VsbEhyZTZOYWptV25QblQ1NTRlNW4xdSIsInN1YiI6InNhbWxwfExpbmRzYXlfTGVla0BobXMuaGFydmFyZC5lZHUiLCJpc3MiOiJodHRwczovL2F2aWxsYWNobGFiLmF1dGgwLmNvbS8iLCJleHAiOjE1NjM3NjI5NzQsImlhdCI6MTU2MzcyNjk3NCwiZW1haWwiOiJMaW5kc2F5X0xlZWtAaG1zLmhhcnZhcmQuZWR1In0.SFIaAb8LrCtg9wJ1qpx3iVtKBRyYIaECA_3Jqtx11ws"
IRCT_REST_BASE_URL = "https://copdgene.hms.harvard.edu"

# Phenotypes of interest for all patients in freeze 5 (n = 10.000)
pheno.df = reloadDfFormFile("data/ids.txt", token, IRCT_REST_BASE_URL)
# make a raw copy
raw.df = pheno.df

rm(list=setdiff(ls(), "raw.df"))
pheno.df = raw.df
# Remove phenotypes that are not used
#pheno.df = subset(pheno.df, select = -c(Gender, Race))
# Column lung cancer is class factor, change to character like others
# pheno.df$Lung_cancer = as.character(pheno.df$Lung_cancer)
# ??? lung cancer still weird, for now remove it
#pheno.df = pheno.df[ ,0:27]  # depending if you want to look at lungfucntions
pheno.df = subset(pheno.df, select = -c(COPD, Lung_cancer))  #Lung_cancer

# Remove "other" and "exclusionary disease" from Affection status such that you only have binary value
pheno.df = pheno.df[pheno.df$Affection_status == "Control" | pheno.df$Affection_status =="Case", ] 
# Check factors:
summary(as.factor(pheno.df$Affection_status))

#Change phenotypic scores to 1s, 0s and NAs 
pheno.df[pheno.df == "Yes"] = 1
pheno.df[pheno.df == "No"]  = 0 
pheno.df[pheno.df == "Do not know"] = NA
#Change affection status to 1s, 0s, but keep as class character
pheno.df$Affection_status[pheno.df$Affection_status == "Case"]    = 1
pheno.df$Affection_status[pheno.df$Affection_status == "Control"] = 0
#Change phenotypic scores from character to numeric
pheno.df[ ,3:ncol(pheno.df)] = lapply(pheno.df[ ,3:ncol(pheno.df)], function(x) {
  if(is.character(x)) as.numeric((x))
})
#Check class:
str(pheno.df)

# s.df has samples of interest: freeze 4 (n = 2000)
s.df = read.csv("data/COPDannotations.txt", 
                stringsAsFactors = FALSE, 
                header = TRUE,  
                sep = "\t")
str(s.df)


#Subset the samples in pheno.df based on samples in s.df
final.df = pheno.df[pheno.df$Patient.Id %in% s.df$dbGaP_Subject_ID, ]
dim(final.df)


##### BUILDING OUTPUT
#df of possible phenotypes
pheno.names = names(final.df[ ,c(3:ncol(final.df))])
#create df to store results in
phewasOutput.df = as.data.frame(matrix(ncol = 7)) 

for(i in 1:length(pheno.names)){
  #general linear regression on phenotype of interest(POI)
  POI.glm = glm(formula = final.df[ ,pheno.names[i]]~Affection_status, family = binomial(), data = final.df, na.action = na.omit)
  #summary(POI.glm)
  
  # confidence interval
  ci = exp(summary(POI.glm)$coefficients["Affection_status1", 1] + qnorm(c(0.025, 0.975)) * summary(POI.glm)$coefficients["Affection_status1", 2])
  #ci = confint(POI.glm)
  
  # caco.df is df with specific phenotype, patient ID, case/ctrl
  caco.df        = final.df[ ,c((pheno.names)[i], "Patient.Id", "Affection_status")]
  caco.df        = na.omit(caco.df)
  # Count number of cases with disease (CaseDisease), without etc etc 
  CaseDisease   = length(unique(caco.df[caco.df[ ,1] == 1 & caco.df$Affection_status == "1", "Patient.Id"]))
  CaseNoDisease = length(unique(caco.df[caco.df[ ,1] == 0 & caco.df$Affection_status == "1", "Patient.Id"]))
  CtrlDisease   = length(unique(caco.df[caco.df[ ,1] == 1 & caco.df$Affection_status == "0", "Patient.Id"]))
  CtrlNoDisease = length(unique(caco.df[caco.df[ ,1] == 0 & caco.df$Affection_status == "0", "Patient.Id"]))
  newRow        = c(pheno.names[i], 
                    round(summary(POI.glm)$coefficients[2] ,2), 
                    round(exp(summary(POI.glm)$coefficients[2]), 2), 
                    paste0("[", round(ci[1][1],3), ", ", round(ci[2][1],3), "]"),
                    if(summary(POI.glm)$coefficients[2,4] > 10e-150){
                      summary(POI.glm)$coefficients[2,4]
                    }else{
                      10e-150
                    }, 
                    paste0(CaseDisease   + CtrlDisease,   "(", CaseDisease,  "/", CtrlDisease,  ")"), 
                    paste0(CaseNoDisease + CtrlNoDisease, "(", CaseNoDisease,"/", CtrlNoDisease,")")
  )
  phewasOutput.df = rbind(newRow, phewasOutput.df)
}

colnames(phewasOutput.df) = c("Phenotype", "Coefficient", "OR", "Confidence_interval", "Pvalue", "Phenotype_present", "Phenotype_absent")

#HOW DO THEY DO BONFERRONI????
phewasOutput.df$adjPvalue = p.adjust(as.numeric(phewasOutput.df$Pvalue), method = "bonferroni")
#Remove all rows with only NAs
phewasOutput.df = phewasOutput.df[complete.cases(phewasOutput.df), ]

phewasOutput.df$Phenotype = gsub("_", " ", phewasOutput.df$Phenotype)


# Significant phenotypes based onf  adj
adjPvalue_Sign.df = phewasOutput.df[as.numeric(phewasOutput.df$adjPvalue) < 0.05, ]
OR_Sign.df   = adjPvalue_Sign.df[as.numeric(adjPvalue_Sign.df$OR) > 2, ] 


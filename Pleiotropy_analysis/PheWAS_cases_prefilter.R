setwd("/Users/lindsayvm/Documents/Avillach_Lab/projects/PheWAS")
library(dplyr)
library(data.table)
library(tidyr)
library(Hmisc)


##### CREATE DATAFRAME

source("src/functions.R")
# Go to https://copdgene.hms.harvard.edu/transmart/login/auth, utilities, user profile, IRCT Token
token = "eyJraWQiOiJSa05CUWpFNU9VTkVOelkzTmpJd04wVkNNVGd3TmpFM01EVXdSVEpETVVaRk5EZzROa0ZFUmciLCJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJhdWQiOiJHSjhmWWlOZ3VsbEhyZTZOYWptV25QblQ1NTRlNW4xdSIsInN1YiI6InNhbWxwfExpbmRzYXlfTGVla0BobXMuaGFydmFyZC5lZHUiLCJpc3MiOiJodHRwczovL2F2aWxsYWNobGFiLmF1dGgwLmNvbS8iLCJleHAiOjE1NjM3NzQ5ODYsImlhdCI6MTU2MzczODk4NiwiZW1haWwiOiJMaW5kc2F5X0xlZWtAaG1zLmhhcnZhcmQuZWR1In0.mUEzJIyt-2cyD9_pbtmwHYmiRC_e9h21P_T0AtGT3f0"
IRCT_REST_BASE_URL = "https://copdgene.hms.harvard.edu"

# Phenotypes of interest for all patients in freeze 5 (n = 10.000)
pheno.df1 = reloadDfFormFile("data/ids.txt", token, IRCT_REST_BASE_URL)
# make a raw copy
raw.df = pheno.df1

rm(list=setdiff(ls(), "raw.df"))
pheno.df1 = raw.df

#summary(as.factor(pheno.df1$Affection_status))
# only COPD cases
pheno.df1 = pheno.df1[pheno.df1$Affection_status =="Case", ] 
summary(as.factor(pheno.df1$Affection_status))

# Remove phenotypes
pheno.df1 = subset(pheno.df1, select = -c(Lung_cancer,              # lung cancer values are all "no"
                                          COPD,                     # all samples are positive for COPD or affection_status
                                          Affection_status,        
                                          Chest_injuries,           #vague description/ no disease
                                          Chest_operations,
                                          Other_chest_illnesses,
                                          Prostate_cancer,
                                          High_cholesterol,
                                          High_blood_pressure
                                          ))

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
GT.df = read.csv("data/Combi/combi_FEV1_pval5e-06.csv",       
                 stringsAsFactors = FALSE, 
                 header = TRUE,  
                 sep = ",")           
dim(GT.df) #155928      5
# Assign 1 if GT call is [1,0] or [0,1] and 0 if [0,0]
GT.df$GT.alleles[GT.df$GT.alleles != "[0, 0]"] = 1
GT.df$GT.alleles[GT.df$GT.alleles == "[0, 0]"] = 0
GT.df$locus = paste0(GT.df$locus.contig, ":", GT.df$locus.position)

#only keep GT calls of cases 
GT.df = GT.df[GT.df$dbGaP_Subject_ID %in% pheno.df1$Patient.Id, ]
dim(GT.df) #66305     6

#Subset the samples in pheno.df1 (freeze 5) based on samples in GT.df (freeze 4)
pheno.df1 = pheno.df1[pheno.df1$Patient.Id %in% GT.df$dbGaP_Subject_ID, ]
dim(pheno.df1) #745  29




####### FILTER

#absolute count
patientID = pheno.df1$Patient.Id
pheno.df1 = pheno.df1[ ,-1]
# Check absolute counts
colSums(pheno.df1, na.rm = TRUE)[order(colSums(pheno.df1, na.rm = TRUE))]
x = colSums(pheno.df1, na.rm = TRUE) > 25 # 0.05*dim(pheno.df1)[1]    #IF you do 25 you can keep stroke; and then remove angina at correlation
pheno.df1 = pheno.df1[ ,x]
# check
colSums(pheno.df1, na.rm = TRUE)

#2. Check correlation between pheno's
correlation = rcorr(as.matrix(pheno.df1), type="pearson")
# Represent data in clearer way
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
corflat = flattenCorrMatrix(correlation$r, correlation$P)
#corflat$Padjust = p.adjust(as.numeric(corflat$p), method = "bonferroni") 
### filter on adj pval?? The fact that it is non-significant (not below < 0.05) is explained by small sample size. 
#cor_pfilter = na.omit(corflat[corflat$Padjust>0.05, ])   ##Rheumatoid_arthritis , Macular_degeneration, Breast_cancer, Transient_Ischemic_Attack, Prostate_cancer, 
#unique(cor_pfilter$row)
#pheno.df1 = subset(pheno.df1, select = unique(cor_pfilter$row))

# Keep those combinations that show low correlation ==> WRONG, it is possible something has high correlation with one pheno has a low corrleation with another, so then you will still not filter it out. 
#tmp = na.omit(corflat[corflat$cor<0.7, ])    # cor 1 means identical, the more identical the more you want them out
#tmp = droplevels(tmp)
#unique(tmp$row)

# Show which combinations have high correlation
na.omit(corflat[corflat$cor>0.3, ])    
pheno.df1 = subset(pheno.df1, select = -c(Bronchitis, Heart_attack, Angina)) #Bronchitis, Angina, Heart_attack, high_cholesterol, 

pheno.df1$Patient.Id = patientID

########PHEWAS PER VARIANT
final.df1 = pheno.df1
pheno.names1 = names(final.df1[ ,-ncol(final.df1)]) #all phenotypes
locus.names = unique(GT.df$locus) # all loci

#Bonferroni correction
signlv = 0.05 # significance level
va = signlv/length(locus.names) #  by number of variants
ph = signlv/length(pheno.names1) # ,, ,, phenotypes
vaph = signlv/length(locus.names)/length(pheno.names1) # variants & phenotypes
print(c("corrected by variants = ", round(va, 4), 
        "corrected by phenotypes = ", round(ph, 4), 
        "corrected by variants & phenotypes = ", round(vaph,8)))

ph = round(ph, 3)

for(l in 1:length(locus.names)){
  # select for variant
  var.df = GT.df[GT.df$locus == locus.names[l], ]
  #annotate GT calls of selected variant
  if(sum(as.numeric(var.df$GT.alleles)) == 0) next # skip iteration
  final.df1$GTalleles = as.character(var.df$GT.alleles)    
  for(i in 1:length(pheno.names1)){
    #general linear regression on phenotype of interest(POI)
    #????? COVARIATES: AGE
    POI.glm1 = glm(formula = final.df1[ ,pheno.names1[i]]~GTalleles, 
                   family = binomial(), 
                   data = final.df1, 
                   na.action = na.omit)
    #summary(POI.glm1)
    if(summary(POI.glm1)$coefficients[2,4] < 0.05){
      print(c(locus.names[l] , pheno.names1[i], summary(POI.glm1)$coefficients[2,4], round(exp(summary(POI.glm1)$coefficients[2]) ,3)))
    }
  }
}

#[5] "3:71096494"          "Stroke"              "0.00280072593954425"
#[7] "3:71098330"          "Stroke"              "0.00280072593954425"
#[54] "6:126189301"         "Pneumonia"          "0.0024143954398921"
#[61] "8:4695411"           "Compression_fractures_in_your_back" "0.00296559921837648"  

# [19] "4:145448056" most sign GWAS hit
# [10] "4:145448056"same but then posisiton in NHW

######################### FOR SINGLE SNP

l = 61
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
                     round(summary(POI.glm1)$coefficients[2] , 3), 
                     round(exp(summary(POI.glm1)$coefficients[2]) ,3), 
                     paste0("[", round(ci1[1][1],3), ", ", round(ci1[2][1],3), "]"),
                     round(summary(POI.glm1)$coefficients[2,4] , 4), 
                     paste0(CaseDisease1   + CtrlDisease1,   "(", CaseDisease1,  "/", CtrlDisease1,  ")"), 
                     paste0(CaseNoDisease1 + CtrlNoDisease1, "(", CaseNoDisease1,"/", CtrlNoDisease1,")"))
  phewasOutput.df1 = rbind(newRow1, phewasOutput.df1)
}
colnames(phewasOutput.df1) = c("Phenotype", "Coefficient", "OR", "Confidence_interval", "Pvalue", "Pheno_present(variant/novariant)", "Pheno_absent(variant/novariant)")
#Remove all rows with only NAs
phewasOutput.df1 = phewasOutput.df1[complete.cases(phewasOutput.df1), ]

phewasOutput.df1$Phenotype = gsub("_", " ", phewasOutput.df1$Phenotype)


phewasOutput.df1$adjustPvalue <- p.adjust( as.numeric(phewasOutput.df1$Pvalue), method = "BH")
#phewasOutput.df1

# Significant phenotypes based onf  adj
#adjPvalue_Sign.df1 = phewasOutput.df1[as.numeric(phewasOutput.df1$adjPvalue) < 0.05, ]
#OR_Sign.df1   = adjPvalue_Sign.df1[as.numeric(adjPvalue_Sign.df1$OR) > 2, ] 

#OR_Sign.df1   = phewasOutput.df1[as.numeric(phewasOutput.df1$OR) > 1.5, ] 
#OR_Sign.df1


#install.packages("dplyr")
setwd("/Users/lindsayvm/Documents/Avillach_Lab/projects/PheWAS")

source("functions.R")

token <- "eyJraWQiOiJSa05CUWpFNU9VTkVOelkzTmpJd04wVkNNVGd3TmpFM01EVXdSVEpETVVaRk5EZzROa0ZFUmciLCJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJhdWQiOiJHSjhmWWlOZ3VsbEhyZTZOYWptV25QblQ1NTRlNW4xdSIsInN1YiI6InNhbWxwfExpbmRzYXlfTGVla0BobXMuaGFydmFyZC5lZHUiLCJpc3MiOiJodHRwczovL2F2aWxsYWNobGFiLmF1dGgwLmNvbS8iLCJleHAiOjE1NjEwODAxMDQsImlhdCI6MTU2MTA0NDEwNCwiZW1haWwiOiJMaW5kc2F5X0xlZWtAaG1zLmhhcnZhcmQuZWR1In0.ue2Y6Hl1glA65aUgpYfRA5QyHGeNo6VDjpWDj__LTZE"

IRCT_REST_BASE_URL <- "https://copdgene.hms.harvard.edu"

transmartVars <-c()
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/00 Affection status/00 Affection status/", targetName="Affection_status"))

transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/01 Demographics/01 Demographics/Race",targetName="Race"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/01 Demographics/01 Demographics/Gender",targetName="Gender"))

transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Angina",targetName="Angina"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Blood clots in legs or lungs",targetName="Blood_clots_in_legs_or_lungs"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Congestive heart failure",targetName="Congestive_heart_failure"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Coronary artery disease",targetName="Coronary_artery_disease"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Diabetes",targetName="Diabetes"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Heart attack",targetName="Heart_attack"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/High blood pressure",targetName="High_blood_pressure"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/High cholesterol",targetName="High_cholesterol"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Peripheral vascular disease",targetName="Peripheral_vascular_disease"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Cardiology/Pneumothorax",targetName="Pneumothorax"))

transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Gastrointestinal/Gastroesophageal reflux",targetName="Gastroesophageal_reflux"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Gastrointestinal/Stomach ulcers",targetName="Stomach_ulcers"))

transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Musculoskeletal/Compression fractures in your back",targetName="Compression_fractures_in_your_back"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Musculoskeletal/Hip fracture",targetName="Hip_fracture"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Musculoskeletal/Osteoarthritis",targetName="Osteoarthritis"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Musculoskeletal/Osteoporosis thin bones",targetName="Osteoporosis_thin_bones"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Musculoskeletal/Rheumatoid arthritis",targetName="Rheumatoid_arthritis"))

transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Neurology/Macular degeneration",targetName="Macular_degeneration"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Neurology/Stroke",targetName="Stroke"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Neurology/Transient Ischemic Attck",targetName="Transient_Ischemic_Attack"))

transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Oncology/Bladder cancer",targetName="Bladder_cancer"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Oncology/Breast cancer",targetName="Breast_cancer"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Oncology/Colon cancer",targetName="Colon_cancer"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Oncology/Lung cancer",targetName="Lung_cancer"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Oncology/Prostate cancer",targetName="Prostate_cancer"))

transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/01 Asthma/01 Have you ever had asthma",targetName="Asthma"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/02 Hay fever/01 Have you ever had hay fever (allergy involving nose or eyes)",targetName="Hay_fever"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/03 Bronchitis/01 Have you ever had an attack of bronchitis",targetName="Bronchitis"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/04 Pneumonia/01 Have you ever had pneumonia or bronchopneumonia",targetName="Pneumonia"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/05 Chronic bronchitis/01 Have you ever had chronic bronchitis",targetName="Chronic_bronchitis"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/06 Emphysema/01 Have you ever had emphysema",targetName="Emphysema"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/07 COPD/01 Have you ever had COPD",targetName="COPD"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/08 Sleep apnea/01 Have you ever had sleep apnea",targetName="Sleep_apnea"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/09 Chest/Any chest injuries",targetName="Chest_injuries"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/09 Chest/Any chest operations",targetName="Chest_operations"))
transmartVars <- c(transmartVars,transmartVar(path="/i2b2-wildfly-default/Demo/03 Clinical data/03 Clinical data/Respiratory disease form/04 Respiratory Conditions/09 Chest/Any other chest illnesses",targetName="Other_chest_illnesses"))
#No encoding supplied: defaulting to UTF-8.
#Show Traceback

#Rerun with Debug
#Error in as.data.frame.default(x[[i]], optional = TRUE) : 
#  cannot coerce class ‘"externalptr"’ to a data.frame 
#In addition: There were 50 or more warnings (use warnings() to see the first 50)



df <- buildDf(transmartVars, token, IRCT_REST_BASE_URL)
df2 <- reloadDfFormFile("COPDtest.txt",token,IRCT_REST_BASE_URL)




###########deprecated

#req <- getData4Variable("/i2b2-wildfly-default/Demo/01 Demographics/01 Demographics/Race","race",T)
#dfRace<-getDataFromReq(req)
#req <- getData4Variable("/i2b2-wildfly-default/Demo/01 Demographics/01 Demographics/Gender","gender",T)
#dfgender<-getDataFromReq(req)
#req <- getData4Variable("/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/01 Health status","health_status",T)
#dfHStatus<-getDataFromReq(req)
#req <- getData4Variable("/i2b2-wildfly-default/Demo/00 Affection status/00 Affection status/","affection_status",T)
#dfafftatus<-getDataFromReq(req)
#req <- getData4Variable("/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Gastrointestinal/Stomach ulcers","history_stomach_ulcers",T)
#dfulc<-getDataFromReq(req)
#req <- getData4Variable("/i2b2-wildfly-default/Demo/02 Medical history/02 Medical history/Medical history form/02 Disease history/Gastrointestinal/Gastroesophageal reflux","history_gastrointestinal_refluc",T)
#dfreflux<-getDataFromReq(req)
#df <- merge(df,df,all = T, by = "Patient.Id")
#df <- merge(df,dfHStatus,all = T)
#df <- merge(df,dfafftatus,all = T)
#df <- merge(df,dfreflux,all = T)
#df <- merge(df,dfulc,all = T)
#Collapse
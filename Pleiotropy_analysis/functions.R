#install.packages("devtools")
#install_github( "hms-dbmi/Rcupcake" )
library(httr)
library(dplyr)
library(devtools)
#library( Rcupcake )


transmartVar <- setClass("transmartVar", slots=list(path="character",targetName="character"))
getRessources <- function (IRCT_REST_BASE_URL,token,path){
  url <- paste0(IRCT_REST_BASE_URL,"/rest/v1/resourceService/path",path)
  getdata<-httr::GET(url=url, httr::add_headers(Authorization=paste0("bearer ",token)))
  manageHttrRequest(getdata)
  
  ressource <- data.frame(unlist(httr::content(getdata)))
  
  return(httr::content(getdata,encoding = "UTF-8"))
}
getDataFromReq <- function (req,token,IRCT_REST_BASE_URL,saveIds,saveIdsFile){
  print(req)
  IRCT_REST_BASE_URL <- "https://copdgene.hms.harvard.edu"
  url <- paste0(IRCT_REST_BASE_URL,"/rest/v1/queryService/runQuery")
  res<-httr::POST(url=url, httr::add_headers(Authorization=paste0("bearer ",token)),body = req)
  manageHttrRequest(res)
  reqId <- httr::content(res,encoding = "UTF-8")$resultId
  
  if(saveIds){
    write(reqId,file = saveIdsFile, append = T)
  }
  print(reqId)
  
  mes <- "`Result` has been initialized."
  while (mes == "`Result` has been initialized.") {
    url <- paste0(IRCT_REST_BASE_URL,paste0("/rest/v1/resultService/resultStatus/",reqId))
    getdata<-httr::GET(url=url, httr::add_headers(Authorization=paste0("bearer ",token)))
    mes <- httr::content(getdata,encoding = "UTF-8")$message    
    print(mes)
    Sys.sleep(3)
  }
  
  dfResult <- retrieveDataFromId(reqId,token,IRCT_REST_BASE_URL)
  return(dfResult)
}
retrieveDataFromId <- function(reqId,token,IRCT_REST_BASE_URL){
  url <- paste0(IRCT_REST_BASE_URL,paste0("/rest/v1/resultService/result/",reqId,"/JSON?download=yes"))
  print(paste0("Retrieve data for ",reqId))
  getdata<-httr::GET(url=url, httr::add_headers(Authorization=paste0("bearer ",token)))
  manageHttrRequest(getdata)
  jsonResult <- jsonlite::fromJSON(httr::content(getdata,type="text",encoding = "UTF-8"),simplifyDataFrame = FALSE)
  
  listResult <- jsonResult[[2]]
  
  dfResult <- data.frame(t(unlist(listResult[[1]])))
  for (i in 2:length(listResult)){
    dfResult<-bind_rows(dfResult,data.frame(t(unlist(listResult[[i]]))))
  }  
  
  return(dfResult)
}
getData4Variable <- function(pui,alias,modality = T,token,IRCT_REST_BASE_URL) {
  if(modality){
    values <- get.children(
      url = IRCT_REST_BASE_URL,
      fieldname = pui,
      verbose = F
    )
    if (is.null(values)){
      values <- c(pui)
    }
    
    print(values)  
    print(length(values))  
    
    if (length(values < 6)){
      req <- paste0(
        '{',
        '"select": ['
      )
      for(i in 1:length(values)){
        req <- paste0(req,
                      '{',
                      '"field": {',
                      '"pui": "',values[i],'",',
                      '"dataType": "STRING"',
                      '},',
                      '"alias": "',alias,'"',
                      '}'
        )
        if (i < length(values)){
          req <- paste0(req,',')
        }
      }
      req <- paste0(req,
                    '],',
                    '"where": [',
                    '{',
                    '"field": {',
                    '"pui": "/i2b2-wildfly-default/Demo/00 Affection status/00 Affection status/",',
                    '"dataType": "STRING"',
                    '},',
                    '"predicate": "CONTAINS",',
                    '"fields": {"ENCOUNTER": "YES"}',
                    '}',
                    ']',
                    '}'
      )
    }
    print(req)
    print(jsonlite::prettify(req))
    return (req)
  }
}
buildDf <- function(transmartVars,token,IRCT_REST_BASE_URL,saveIds = T,saveIdsFile = "ids.txt"){
  if(saveIds){
    write("",file = saveIdsFile, append = F)
  }
  for(i in 1:length(transmartVars)){
    print("======Querying data ======") 
    print(transmartVars[[i]]@path) 
    print(transmartVars[[i]]@targetName)
    print("======Querying data ======") 
    req <- getData4Variable(transmartVars[[i]]@path,transmartVars[[i]]@targetName,T,token,IRCT_REST_BASE_URL)
    dfT<-getDataFromReq(req,token,IRCT_REST_BASE_URL,saveIds,saveIdsFile)
    print("======Querying end ======") 
    
    
    if (i == 1){
      df <- dfT
    }else{
      df <- merge(df,dfT,all = T,by = "Patient.Id")
    }
  }
  
  return(df)
}
reloadDfFormFile <- function(file,token,IRCT_REST_BASE_URL){
  print(paste0("=============== Building from File ==> ", file, "===================="))
  con = file(file)
  line = readLines(con)
  j=0
  for( i in 1:length(line)){
    if (line[i] !=""){
      j= j+1
      dfT <- retrieveDataFromId(line[i],token,IRCT_REST_BASE_URL)
      if (j == 1){
        df <- dfT
      }else{
        df <- merge(df,dfT,all = T,by = "Patient.Id")
      }
    }
  }
  
  return(df)
}
printRessourceList<-function(ressourceList){
  for (i in 1:length(ressourceList)){
    print(paste0("========== Ressource => [", i ,"] ==============="))
    print(paste0("Path ==> ", ressourceList[[i]]@path))
    print(paste0("targetName ==> ", ressourceList[[i]]@targetName))
}}
manageHttrRequest <- function (httrRequest){
  if(httr::http_error(httrRequest)){
    print("PIC-SURE - Request ERROR ******************")
    print(paste("Status : ", httr::http_status(httrRequest)))
    print(paste("message : ", httr::content(httrRequest,as="text",encoding = "UTF-8")))
  }
}


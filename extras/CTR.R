# setting environments 

.libPaths("/home/aj/env/")
setwd("/home/aj/sooj/epigeneral/")

library(evidnet)
library(DBI)
library(SqlRender)

oracleTempSchema <- NULL
tempEmulationSchema = oracleTempSchema
options(sqlRenderTempEmulationSchema = oracleTempSchema)
options(andromedaTempFolder = "./androtemp")
maxCores <- 10
# Specify where the temporary files (used by the ff package) will be created:
options(fftempdir = "./temp")

cohortTable <- "ContKtAki_epigeneral"

URL <- Sys.getenv("MEDICAL_RECORDS_URL")
DB_HOST <- strsplit(URL, ":")[[1]][1]
options(fftempdir = "./temp")
dbms <- "postgresql"
user <- Sys.getenv("MEDICAL_RECORDS_USER")
pw <- Sys.getenv("MEDICAL_RECORDS_PW")
server <- paste0(DB_HOST)
port <- strsplit(URL, ":")[[1]][2]
maxCores <- 1
outputFolder <- "/data/results"
# downloadJdbcDrivers('postgresql','.')

connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = dbms,
  server = paste0(server, "/", "omop"),
  user = user,
  password = pw,
  port = port,
  pathToDriver = "/home/aj/sooj/jdbc/"
)
cat("DBconnector library load OK!\n")

cdmDatabaseSchema <- Sys.getenv("MEDICAL_RECORDS_SCHEMA")
SCHEMA <- Sys.getenv("MEDICAL_RECORDS_SCHEMA")
oracleTempSchema <- NULL
options(sqlRenderTempEmulationSchema = NULL)
cohortDatabaseSchema <- paste0(cdmDatabaseSchema, "_results_dq_v276")
execute(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = paste0(cdmDatabaseSchema, "_results_dq_v276"),
  cohortTable = cohortTable,
  outputFolder = outputFolder,
  databaseName = "omop",
  verifyDependencies = FALSE,
  createCohorts = TRUE,
  synthesizePositiveControls = FALSE,
  runAnalyses = F,
  packageResults = FALSE,
  maxCores = maxCores
)

#name <- paste0(cohortDatabaseSchema, ".", cohortTable)
#connect <- DatabaseConnector::connect(connectionDetails)
#cohortoutcome <- dbReadTable(conn = connect, name = name)
#cohortoutput <- paste0(outputFolder, "/", "cohortorigin.csv")
#write.csv(cohortoutcome, file = cohortoutput)

connection <- connect(connectionDetails)
sql <- readSql("/home/aj/sooj/epigeneral/extras/covariateCohorts.sql")
sql <- SqlRender::render(sql,
                         cohort_database_schema = cohortDatabaseSchema
)
sql <- SqlRender::translate(sql, targetDialect = attr(connection, "dbms"), tempEmulationSchema = tempEmulationSchema)
DatabaseConnector::executeSql(connection, sql)


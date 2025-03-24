library(mdrd)

outputFolder <- "/home/gansujin/synology/storage/10/output/2025/ContKT/mdrd/T01"

# Specifmdrd# Specify where the temporary files (used by the ff package) will be created:
options(fftempdir = "./temp")

dbms <- "sql server"
user <- 'gansujin'
pw <- 'sujin30401@'
server <- '10.5.99.50'
port <- '1433'
DATABASECONNECTOR_JAR_FOLDER<- '~/Users/gansujin/path/mssql'
Sys.setenv("DATABASECONNECTOR_JAR_FOLDER" = '~/Users/gansujin/path/mssql')

connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = dbms,
                                                                server = server,
                                                                user = user,
                                                                password = pw,
                                                                port = port,
                                                                pathToDriver = DATABASECONNECTOR_JAR_FOLDER)

# Add the database containing the OMOP CDM data
cdmDatabaseSchema <- 'CDMPv536.dbo'
cohortDatabaseSchema <- 'cohortDb.dbo'
cdmDatabaseName <- 'AUSOM'

# first only sql query 수정 to preserve all the records at the first cohort 
cohortTable <- "sooj_mdrd_250321_T01"

oracleTempSchema <- NULL
options(sqlRenderTempEmulationSchema = NULL)
options(andromedaTempFolder = "./androtemp")

maxCores <- 1

execute(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  outputFolder = outputFolder,
  verifyDependencies = FALSE,
  createCohorts = F,
  synthesizePositiveControls = FALSE,
  runAnalyses = T,
  packageResults = T,
  maxCores = maxCores
)

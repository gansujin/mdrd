library(SqlRender)
library(DatabaseConnector)

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


cdmDatabaseSchema <- 'CDMPv536.dbo'
cohortDatabaseSchema <- 'cohortDb.dbo'
cdmDatabaseName <- 'AUSOM'
cohortTable <- "sooj_mdrd_250321_T01"

sql <- readSql("~/synology/storage/10/study/2025/ContKT/mdrd/package/extras/eGFR.sql")

connection <- connect(connectionDetails)
renderTranslateExecuteSql(
  connection = connection,
  sql = sql,
  target_database_schema = cohortDatabaseSchema,
  cohort_database_schema = cohortDatabaseSchema,
  cohort_table = cohortTable
)


sql <- paste(
  "SELECT cohort_definition_id,
COUNT(*) AS count",
  "FROM @cohort_database_schema.@cohort_table",
  "GROUP BY cohort_definition_id"
)

renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  cohort_database_schema = cohortDatabaseSchema,
  cohort_table = cohortTable
)


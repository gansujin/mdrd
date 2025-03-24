--- DROP TABLE IF EXISTS @cohort_database_schema.@cohort_table;
DROP TABLE IF EXISTS #gfr_cohort;
DROP TABLE IF EXISTS #gfr_temp_cohort;

select person_id, start_date, start_date as cohort_end_date, gfr
into #gfr_cohort
from @target_database_schema.gfr_comparator_mdrd
union all
select person_id, start_date, start_date as cohort_end_date, gfr
from @target_database_schema.gfr_target_mdrd;


CREATE TABLE #gfr_temp_cohort (
    cohort_definition_id INT,
    subject_id BIGINT,
    cohort_start_date DATE,
    cohort_end_date DATE
);

-- Assign cohort_definition_id based on GFR range
INSERT INTO #gfr_temp_cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
SELECT 
    CASE 
        WHEN ROUND(gfr, 0) BETWEEN 0 AND 15 THEN 1011
        WHEN ROUND(gfr, 0) BETWEEN 15 AND 30 THEN 1012
        WHEN ROUND(gfr, 0) BETWEEN 30 AND 45 THEN 1013
        WHEN ROUND(gfr, 0) BETWEEN 45 AND 60 THEN 1014
        WHEN ROUND(gfr, 0) > 60 THEN 1015
    END AS cohort_definition_id,
    person_id,
    start_date, start_date
FROM #gfr_cohort;

INSERT INTO @cohort_database_schema.@cohort_table (
  cohort_definition_id,
  subject_id,
  cohort_start_date,
  cohort_end_date
)
select cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
FROM #gfr_temp_cohort;

DELETE FROM @cohort_database_schema.@cohort_table
WHERE Cohort_definition_id is null;
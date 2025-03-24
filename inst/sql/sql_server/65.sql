drop table if exists #cr;
drop table if exists #ct_cohort;
drop table if exists #final_cohort;

select B.person_id, B.measurement_datetime, B.measurement_date, B.value_as_number into #cr 
from @target_database_schema.@target_cohort_table A, @cdm_database_schema.measurement B
where A.subject_id = B.person_id
and b.measurement_date between dateadd(day, -2, a.cohort_start_date) and dateadd(day, 7, a.cohort_start_date) --temporary 3
and measurement_concept_id IN (3016723,3051825)
and value_as_number is not null; -- serum creatinine

select subject_id as person_id, cohort_start_date as ct_date, 
       avg(value_as_number) as base_cr -- median value of CT -2~0 days creatinin
       into #ct_cohort
from @target_database_schema.@target_cohort_table A, #cr B -- using target, compartor cohorts
where a.subject_id = B.person_id and b.measurement_date between dateadd(day, -2, a.cohort_start_date) and a.cohort_start_date
group by subject_id, cohort_start_date;

select distinct A.person_id, b.measurement_date as start_date, dateadd(day, 1, b.measurement_date) as end_date
into #final_cohort
from #ct_cohort A, #cr B
where A.person_id = B.person_id
  and (((A.base_cr + 0.3 < B.value_as_number and b.measurement_date between a.ct_date and  dateadd(day, 2, a.ct_date))) -- >base_cr + 0.3 in 1~3 days
 or (A.base_cr * 1.5 < B.value_as_number and b.measurement_date between a.ct_date and  dateadd(day, 7, a.ct_date))); -- >base_cr *1.5 in 1~7 days

DELETE FROM @target_database_schema.@target_cohort_table where cohort_definition_id = @target_cohort_id; -- target_cohort_id
INSERT INTO @target_database_schema.@target_cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
select @target_cohort_id as cohort_definition_id, person_id, start_date, end_date
FROM #final_cohort CO;

drop table if exists #final_cohort;

SET search_path TO mimiciii;

DROP MATERIALIZED VIEW IF EXISTS Glen_Multivariate_CPMs_Cohort CASCADE;

CREATE MATERIALIZED VIEW Glen_Multivariate_CPMs_Cohort AS
WITH cohort AS (
SELECT ie.subject_id, ie.hadm_id, ie.icustay_id

, ie.intime
, ie.outtime
, ie.los AS icu_los_days
, DENSE_RANK() OVER (PARTITION BY ie.hadm_id ORDER BY ie.intime) AS icustay_seq --calculate a count of how many ICU admissions each patient had within a given hospitalisation

-- Patient level factors
, pat.dob
, ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) AS age
, pat.gender

-- Hospital level factors
, adm.admittime
, adm.dischtime
, adm.deathtime
, adm.hospital_expire_flag
, DENSE_RANK() OVER (PARTITION BY adm.subject_id ORDER BY adm.admittime) AS hospstay_seq --define the number of hostpial admissions for each patient
, adm.admission_type
, CASE WHEN adm.ethnicity IN
  (
       'WHITE'
     , 'WHITE - RUSSIAN' 
     , 'WHITE - OTHER EUROPEAN' 
     , 'WHITE - BRAZILIAN'
     , 'WHITE - EASTERN EUROPEAN' 
  ) THEN 'white' --define a white ethnic group
  WHEN adm.ethnicity IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) THEN 'black' --define a 'black' ethnic group
  WHEN adm.ethnicity IN
    (
      'HISPANIC OR LATINO' 
    , 'HISPANIC/LATINO - PUERTO RICAN' 
    , 'HISPANIC/LATINO - DOMINICAN' 
    , 'HISPANIC/LATINO - GUATEMALAN' 
    , 'HISPANIC/LATINO - CUBAN' 
    , 'HISPANIC/LATINO - SALVADORAN' 
    , 'HISPANIC/LATINO - CENTRAL AMERICAN (OTHER)' 
    , 'HISPANIC/LATINO - MEXICAN' 
    , 'HISPANIC/LATINO - COLOMBIAN' 
    , 'HISPANIC/LATINO - HONDURAN' 
  ) THEN 'hispanic' --define an 'hispanic' ethnic group
  WHEN adm.ethnicity IN
  (
      'ASIAN'
    , 'ASIAN - CHINESE' 
    , 'ASIAN - ASIAN INDIAN'
    , 'ASIAN - VIETNAMESE' 
    , 'ASIAN - FILIPINO' 
    , 'ASIAN - CAMBODIAN' 
    , 'ASIAN - OTHER' 
    , 'ASIAN - KOREAN' 
    , 'ASIAN - JAPANESE' 
    , 'ASIAN - THAI' 
  ) THEN 'asian' --define an 'asian' ethnic group
  WHEN adm.ethnicity IN
  (
      'UNKNOWN/NOT SPECIFIED' 
    , 'UNABLE TO OBTAIN' 
    , 'PATIENT DECLINED TO ANSWER' 
  ) THEN 'unknown' --define an unknown group
  ELSE 'other' END AS ethnicity_grouped

-- First 24 hour lab results (lfd prefix = "lab first day")
, lfd.BICARBONATE_Min
, lfd.BICARBONATE_Mean
, lfd.BICARBONATE_Max
, lfd.CREATININE_Min
, lfd.CREATININE_Mean
, lfd.CREATININE_Max
, lfd.CHLORIDE_Min
, lfd.CHLORIDE_Mean
, lfd.CHLORIDE_Max
, lfd.HEMOGLOBIN_Min
, lfd.HEMOGLOBIN_Mean
, lfd.HEMOGLOBIN_Max
, lfd.PLATELET_Min
, lfd.PLATELET_Mean
, lfd.PLATELET_Max
, lfd.POTASSIUM_Min
, lfd.POTASSIUM_Mean
, lfd.POTASSIUM_Max
, lfd.PTT_Min
, lfd.PTT_Mean
, lfd.PTT_Max
, lfd.INR_Min
, lfd.INR_Mean
, lfd.INR_Max
, lfd.PT_Min
, lfd.PT_Mean
, lfd.PT_Max
, lfd.BUN_Min
, lfd.BUN_Mean
, lfd.BUN_Max
, lfd.WBC_Min
, lfd.WBC_Mean
, lfd.WBC_Max
-- Calculate baseline eGFR using the modification of diet in renal disease equation (Levey et al):
, CASE 
   WHEN (ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2)) > 89 --if age over 89 then need to assume 89 due to additions that MIMIC add to ages over 89 (i.e. adding 300 years)
   THEN 
   CASE WHEN adm.ethnicity IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender = 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * (89^(-0.203)) * 0.742 * 1.212
  WHEN adm.ethnicity NOT IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender = 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * (89^(-0.203)) * 0.742
  WHEN adm.ethnicity IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender <> 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * (89^(-0.203)) * 1.212
  WHEN adm.ethnicity NOT IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender <> 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * (89^(-0.203))
  ELSE null --eGFR not defined for age = 0 years: i.e. if ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) = 0
  END
  
 ELSE  
  CASE WHEN adm.ethnicity IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender = 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * ((ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2))^(-0.203)) * 0.742 * 1.212
  WHEN adm.ethnicity NOT IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender = 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * ((ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2))^(-0.203)) * 0.742
  WHEN adm.ethnicity IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender <> 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * ((ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2))^(-0.203)) * 1.212
  WHEN adm.ethnicity NOT IN
  (
      'BLACK/AFRICAN AMERICAN' 
    , 'BLACK/CAPE VERDEAN' 
    , 'BLACK/HAITIAN' 
    , 'BLACK/AFRICAN' 
    , 'CARIBBEAN ISLAND' 
  ) AND pat.gender <> 'F' 
    AND ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) > 0
  THEN 175 * (lfd.CREATININE_Min^(-1.154)) * ((ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2))^(-0.203))
  ELSE null --eGFR not defined for age = 0 years: i.e. if ROUND((cast(ie.intime as date) - cast(pat.dob as date))/365.242, 2) = 0
  END 
END AS baseline_eGFR
  


-- First 24 hour vital sign results (vfd prefix = "vital first day")
, vfd.HeartRate_Min
, vfd.HeartRate_Mean
, vfd.HeartRate_Max
, vfd.SysBP_Min
, vfd.SysBP_Mean
, vfd.SysBP_Max
, vfd.DiasBP_Min
, vfd.DiasBP_Mean
, vfd.DiasBP_Max
, vfd.MeanBP_Min
, vfd.MeanBP_Mean
, vfd.MeanBP_Max
, vfd.RespRate_Min
, vfd.RespRate_Mean
, vfd.RespRate_Max
, vfd.TempC_Min
, vfd.TempC_Mean
, vfd.TempC_Max
, vfd.SpO2_Min
, vfd.SpO2_Mean
, vfd.SpO2_Max
, vfd.Glucose_Min
, vfd.Glucose_Mean
, vfd.Glucose_Max

-- Select the variables relating to the creatinine values on day 2/3 of each ICU admission
, followup_labs.CREATININE_Min AS CREATININE_Min_Day2_3
, followup_labs.CREATININE_Mean AS CREATININE_Mean_Day2_3
, followup_labs.CREATININE_Max AS CREATININE_Max_Day2_3


FROM icustays ie


INNER JOIN admissions adm
    ON ie.hadm_id = adm.hadm_id


INNER JOIN patients pat
    ON ie.subject_id = pat.subject_id


INNER JOIN
(
SELECT
  pvt.subject_id, pvt.hadm_id, pvt.icustay_id
  --Summarise each patient's lab tests over the first 24 hours for a given ICU admission
  , min(CASE WHEN lab_label = 'BICARBONATE' THEN valuenum ELSE null END) AS BICARBONATE_Min
  , avg(CASE WHEN lab_label = 'BICARBONATE' THEN valuenum ELSE null END) AS BICARBONATE_Mean
  , max(CASE WHEN lab_label = 'BICARBONATE' THEN valuenum ELSE null END) AS BICARBONATE_Max
  , min(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Min
  , avg(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Mean
  , max(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Max
  , min(CASE WHEN lab_label = 'CHLORIDE' THEN valuenum ELSE null END) AS CHLORIDE_Min
  , avg(CASE WHEN lab_label = 'CHLORIDE' THEN valuenum ELSE null END) AS CHLORIDE_Mean
  , max(CASE WHEN lab_label = 'CHLORIDE' THEN valuenum ELSE null END) AS CHLORIDE_Max
  , min(CASE WHEN lab_label = 'HEMOGLOBIN' THEN valuenum ELSE null END) AS HEMOGLOBIN_Min
  , avg(CASE WHEN lab_label = 'HEMOGLOBIN' THEN valuenum ELSE null END) AS HEMOGLOBIN_Mean
  , max(CASE WHEN lab_label = 'HEMOGLOBIN' THEN valuenum ELSE null END) AS HEMOGLOBIN_Max
  , min(CASE WHEN lab_label = 'PLATELET' THEN valuenum ELSE null END) AS PLATELET_Min
  , avg(CASE WHEN lab_label = 'PLATELET' THEN valuenum ELSE null END) AS PLATELET_Mean
  , max(CASE WHEN lab_label = 'PLATELET' THEN valuenum ELSE null END) AS PLATELET_Max
  , min(CASE WHEN lab_label = 'POTASSIUM' THEN valuenum ELSE null END) AS POTASSIUM_Min
  , avg(CASE WHEN lab_label = 'POTASSIUM' THEN valuenum ELSE null END) AS POTASSIUM_Mean
  , max(CASE WHEN lab_label = 'POTASSIUM' THEN valuenum ELSE null END) AS POTASSIUM_Max
  , min(CASE WHEN lab_label = 'PTT' THEN valuenum ELSE null END) AS PTT_Min
  , avg(CASE WHEN lab_label = 'PTT' THEN valuenum ELSE null END) AS PTT_Mean
  , max(CASE WHEN lab_label = 'PTT' THEN valuenum ELSE null END) AS PTT_Max
  , min(CASE WHEN lab_label = 'INR' THEN valuenum ELSE null END) AS INR_Min
  , avg(CASE WHEN lab_label = 'INR' THEN valuenum ELSE null END) AS INR_Mean
  , max(CASE WHEN lab_label = 'INR' THEN valuenum ELSE null END) AS INR_Max
  , min(CASE WHEN lab_label = 'PT' THEN valuenum ELSE null END) AS PT_Min
  , avg(CASE WHEN lab_label = 'PT' THEN valuenum ELSE null END) AS PT_Mean
  , max(CASE WHEN lab_label = 'PT' THEN valuenum ELSE null END) AS PT_Max
  , min(CASE WHEN lab_label = 'BUN' THEN valuenum ELSE null end) AS BUN_Min
  , avg(CASE WHEN lab_label = 'BUN' THEN valuenum ELSE null end) AS BUN_Mean
  , max(CASE WHEN lab_label = 'BUN' THEN valuenum ELSE null end) AS BUN_Max
  , min(CASE WHEN lab_label = 'WBC' THEN valuenum ELSE null end) AS WBC_Min
  , avg(CASE WHEN lab_label = 'WBC' THEN valuenum ELSE null end) AS WBC_Mean
  , max(CASE WHEN lab_label = 'WBC' THEN valuenum ELSE null end) AS WBC_Max
FROM
(
  SELECT ie.subject_id, ie.hadm_id, ie.icustay_id
  -- combine multiple ITEMIDs containing the same data:
  , CASE
        WHEN itemid = 50882 THEN 'BICARBONATE'
        WHEN itemid = 50912 THEN 'CREATININE'
        WHEN itemid = 50806 THEN 'CHLORIDE'
        WHEN itemid = 50902 THEN 'CHLORIDE'
        WHEN itemid = 50811 THEN 'HEMOGLOBIN'
        WHEN itemid = 51222 THEN 'HEMOGLOBIN'
        WHEN itemid = 51265 THEN 'PLATELET'
        WHEN itemid = 50822 THEN 'POTASSIUM'
        WHEN itemid = 50971 THEN 'POTASSIUM'
        WHEN itemid = 51275 THEN 'PTT'
        WHEN itemid = 51237 THEN 'INR'
        WHEN itemid = 51274 THEN 'PT'
        WHEN itemid = 51006 THEN 'BUN'
        WHEN itemid = 51300 THEN 'WBC'
        WHEN itemid = 51301 THEN 'WBC'
      ELSE null
    END AS lab_label
    
  -- Conduct sanity checks on the values
  ,  CASE
      WHEN itemid = 50882 and valuenum > 10000 THEN null -- mEq/L 'BICARBONATE'
      WHEN itemid = 50806 and valuenum > 10000 THEN null -- mEq/L 'CHLORIDE'
      WHEN itemid = 50902 and valuenum > 10000 THEN null -- mEq/L 'CHLORIDE'
      WHEN itemid = 50912 and valuenum >   150 THEN null -- mg/dL 'CREATININE'
      WHEN itemid = 50811 and valuenum >    50 THEN null -- g/dL 'HEMOGLOBIN'
      WHEN itemid = 51222 and valuenum >    50 THEN null -- g/dL 'HEMOGLOBIN'
      WHEN itemid = 51265 and valuenum > 10000 THEN null -- K/uL 'PLATELET'
      WHEN itemid = 50822 and valuenum >    30 THEN null -- mEq/L 'POTASSIUM'
      WHEN itemid = 50971 and valuenum >    30 THEN null -- mEq/L 'POTASSIUM'
      WHEN itemid = 51275 and valuenum >   150 THEN null -- sec 'PTT'
      WHEN itemid = 51237 and valuenum >    50 THEN null -- 'INR'
      WHEN itemid = 51274 and valuenum >   150 THEN null -- sec 'PT'
      WHEN itemid = 51006 and valuenum >   300 THEN null -- 'BUN'
      WHEN itemid = 51300 and valuenum >  1000 THEN null -- 'WBC'
      WHEN itemid = 51301 and valuenum >  1000 THEN null -- 'WBC'
    ELSE le.valuenum
    END AS valuenum

  FROM icustays ie
  LEFT JOIN labevents le
    ON le.subject_id = ie.subject_id 
    AND le.hadm_id = ie.hadm_id
    AND le.charttime BETWEEN (ie.intime - interval '6' hour) 
                     AND (ie.intime + interval '1' day) --extract the lab events occuring within 1 day (minus a small buffer to capture any 'pre ICU lab tests')
    AND le.ITEMID in --Only extract the lab events we are interested in for this study:
    (
      50882, -- BICARBONATE 
      50912, -- CREATININE 
      50902, -- CHLORIDE 
      50806, -- CHLORIDE, WHOLE BLOOD 
      51222, -- HEMOGLOBIN 
      50811, -- HEMOGLOBIN
      51265, -- PLATELET COUNT 
      50971, -- POTASSIUM
      50822, -- POTASSIUM, WHOLE BLOOD
      51275, -- PARTIAL THROMBOPLASTIN TIME (PTT) 
      51237, -- INTERNATIONAL NORMALIZED RATIO (INR) 
      51274, -- PROTHROMBIN TIME (PT) 
      51006, -- BLOOD UREA NITROGEN (BUN)
      51301, -- WHITE BLOOD CELLS (WBC)
      51300  -- WBC COUNT 
    )
    AND valuenum IS NOT null AND valuenum > 0 -- lab values cannot be 0 and cannot be negative
 ) pvt
GROUP BY pvt.subject_id, pvt.hadm_id, pvt.icustay_id --Group by statement to ensure summarises of lab values are taken for each patient/hospitalisation/ICU event
) lfd --prefix with lfd (lab first day)
    ON (ie.subject_id = lfd.subject_id AND 
        ie.hadm_id = lfd.hadm_id AND 
        ie.icustay_id = lfd.icustay_id)


INNER JOIN --Here, inner join will cause some ICU admissions to be dropped (any pt. with only vital signs that were collected with error)
(
SELECT 
  pvt.subject_id, pvt.hadm_id, pvt.icustay_id
  --Summarise each patient's vital signs over the first 24 hours for a given ICU admission
  , min(CASE WHEN vital_label = 'HR' THEN valuenum ELSE null END) AS HeartRate_Min
  , avg(CASE WHEN vital_label = 'HR' THEN valuenum ELSE null END) AS HeartRate_Mean
  , max(CASE WHEN vital_label = 'HR' THEN valuenum ELSE null END) AS HeartRate_Max
  , min(CASE WHEN vital_label = 'SBP' THEN valuenum ELSE null END) AS SysBP_Min
  , avg(CASE WHEN vital_label = 'SBP' THEN valuenum ELSE null END) AS SysBP_Mean
  , max(CASE WHEN vital_label = 'SBP' THEN valuenum ELSE null END) AS SysBP_Max
  , min(CASE WHEN vital_label = 'DBP' THEN valuenum ELSE null END) AS DiasBP_Min
  , avg(CASE WHEN vital_label = 'DBP' THEN valuenum ELSE null END) AS DiasBP_Mean
  , max(CASE WHEN vital_label = 'DBP' THEN valuenum ELSE null END) AS DiasBP_Max
  , min(CASE WHEN vital_label = 'BP' THEN valuenum ELSE null END) AS MeanBP_Min
  , avg(CASE WHEN vital_label = 'BP' THEN valuenum ELSE null END) AS MeanBP_Mean
  , max(CASE WHEN vital_label = 'BP' THEN valuenum ELSE null END) AS MeanBP_Max
  , min(CASE WHEN vital_label = 'RR' THEN valuenum ELSE null END) AS RespRate_Min
  , avg(CASE WHEN vital_label = 'RR' THEN valuenum ELSE null END) AS RespRate_Mean
  , max(CASE WHEN vital_label = 'RR' THEN valuenum ELSE null END) AS RespRate_Max
  , min(CASE WHEN vital_label = 'Temp' THEN valuenum ELSE null END) AS TempC_Min
  , avg(CASE WHEN vital_label = 'Temp' THEN valuenum ELSE null END) AS TempC_Mean
  , max(CASE WHEN vital_label = 'Temp' THEN valuenum ELSE null END) AS TempC_Max
  , min(CASE WHEN vital_label = 'SpO2' THEN valuenum ELSE null END) AS SpO2_Min
  , avg(CASE WHEN vital_label = 'SpO2' THEN valuenum ELSE null END) AS SpO2_Mean
  , max(CASE WHEN vital_label = 'SpO2' THEN valuenum ELSE null END) AS SpO2_Max
  , min(CASE WHEN vital_label = 'Glucose' THEN valuenum ELSE null END) AS Glucose_Min
  , avg(CASE WHEN vital_label = 'Glucose' THEN valuenum ELSE null END) AS Glucose_Mean
  , max(CASE WHEN vital_label = 'Glucose' THEN valuenum ELSE null END) AS Glucose_Max
  
FROM  (
  SELECT ie.subject_id, ie.hadm_id, ie.icustay_id
  -- combine multiple ITEMIDs containing the same data:
  , CASE
	WHEN itemid IN (211,220045) AND valuenum > 0 AND valuenum < 300 THEN 'HR' -- Heart Rate
	WHEN itemid IN (51,442,455,6701,220179,220050) AND valuenum > 0 AND valuenum < 400 THEN 'SBP' -- Systolic Blood Pressure
	WHEN itemid IN (8368,8440,8441,8555,220180,220051) AND valuenum > 0 AND valuenum < 300 THEN 'DBP' -- Diastolic Blood Pressure
	WHEN itemid IN (456,52,6702,443,220052,220181,225312) AND valuenum > 0 AND valuenum < 300 THEN 'BP' -- Mean Arterial Blood Pressure
	WHEN itemid IN (615,618,220210,224690) AND valuenum > 0 AND valuenum < 70 THEN 'RR' -- Respiration Rate
	WHEN itemid IN (223761,678) AND valuenum > 70 AND valuenum < 120  THEN 'Temp' 
	WHEN itemid IN (223762,676) AND valuenum > 10 AND valuenum < 50  THEN 'Temp' 
	WHEN itemid IN (646,220277) AND valuenum > 0 AND valuenum <= 100 THEN 'SpO2' 
	WHEN itemid IN (807,811,1529,3745,3744,225664,220621,226537) AND valuenum > 0 THEN 'Glucose'
    ELSE null END AS vital_label

  , CASE WHEN itemid IN (223761,678) THEN (valuenum-32)/1.8 ELSE valuenum END AS valuenum --convert temperature to C

  FROM icustays ie
  LEFT JOIN chartevents ce
    ON ie.subject_id = ce.subject_id 
    AND ie.hadm_id = ce.hadm_id 
    AND ie.icustay_id = ce.icustay_id --extract only the vital measures relating to each ICU admission
    AND ce.charttime BETWEEN ie.intime AND (ie.intime + interval '1' day) --extract the vital events occuring within 1 day after admission
    AND ce.error IS DISTINCT FROM 1 -- exclude rows marked as error for the vital signs
    AND ce.itemid IN --only link the chart events we are interested in for this study:
    (
    211,    --Heart Rate
    220045, --Heart Rate
    51,     --Arterial BP [Systolic]
    442,    --Manual BP [Systolic]
    455,    --NBP [Systolic]
    6701,   --Arterial BP #2 [Systolic]
    220179, --Non Invasive Blood Pressure systolic
    220050, --Arterial Blood Pressure systolic
    8368,   --Arterial BP [Diastolic]
    8440,   --Manual BP [Diastolic]
    8441,   --NBP [Diastolic]
    8555,   --Arterial BP #2 [Diastolic]
    220180, --Non Invasive Blood Pressure diastolic
    220051, --Arterial Blood Pressure diastolic
    456,    --NBP Mean
    52,     --Arterial BP Mean
    6702,   --Arterial BP Mean #2
    443,    --Manual BP Mean(calc)
    220052, --Arterial Blood Pressure mean
    220181, --Non Invasive Blood Pressure mean
    225312, --ART BP mean
    618,    --Respiratory Rate
    615,    --Resp Rate (Total)
    220210, --Respiratory Rate
    224690, --Respiratory Rate (Total)
    646,    --SPO2, peripheral
    220277, --SPO2, peripheral
    807,    --Fingerstick Glucose
    811,    --Glucose (70-105)
    1529,   --Glucose
    3745,   --BloodGlucose
    3744,   --Blood Glucose
    225664, --Glucose finger stick
    220621, --Glucose (serum)
    226537, --Glucose (whole blood)
    223762, --Temperature Celsius
    676,    --Temperature C
    223761, --Temperature Fahrenheit
    678     --Temperature F

  )
) pvt
GROUP BY pvt.subject_id, pvt.hadm_id, pvt.icustay_id --Group by statement to ensure summarises of vital signs are taken for each patient/hospitalisation/ICU event
) vfd --prefix with vfd (vital first day)
    ON (ie.subject_id = vfd.subject_id AND 
        ie.hadm_id = vfd.hadm_id AND 
        ie.icustay_id = vfd.icustay_id)

-- Extract the creatinine values on day 2 and day 3 of each ICU admission (summarised as min, mean, max for each patient/ICU admission):
LEFT JOIN
(
SELECT
  pvt.subject_id, pvt.hadm_id, pvt.icustay_id
  --Summarise each patient's lab tests over day 2 for a given ICU admission
  , min(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Min
  , avg(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Mean
  , max(CASE WHEN lab_label = 'CREATININE' THEN valuenum ELSE null END) AS CREATININE_Max
FROM
(
  SELECT ie.subject_id, ie.hadm_id, ie.icustay_id
  -- combine multiple ITEMIDs containing the same data:
  , CASE
        WHEN itemid = 50912 THEN 'CREATININE'
      ELSE null
    END AS lab_label
    
  , -- Conduct sanity checks on the values
    CASE
      WHEN itemid = 50912 and valuenum >   150 THEN null -- mg/dL 'CREATININE'
    ELSE le.valuenum
    END AS valuenum

  FROM icustays ie
  LEFT JOIN labevents le
    ON le.subject_id = ie.subject_id 
    AND le.hadm_id = ie.hadm_id
    AND le.charttime BETWEEN (ie.intime + interval '1' day) AND (ie.intime + interval '3' day) --extract the lab events occuring between day 1 and day 3
    AND le.ITEMID IN (50912) -- CREATININE
    AND valuenum IS NOT null AND valuenum > 0 -- lab values cannot be 0 and cannot be negative
 ) pvt
GROUP BY pvt.subject_id, pvt.hadm_id, pvt.icustay_id --Group by statement to ensure summarises of lab values are taken for each patient/hospitalisation/ICU event
) followup_labs
ON (ie.subject_id = followup_labs.subject_id AND 
    ie.hadm_id = followup_labs.hadm_id AND 
    ie.icustay_id = followup_labs.icustay_id)

ORDER BY ie.subject_id, adm.admittime, ie.intime
)
SELECT *
FROM cohort
-- Apply inclusion/exclusion criteria to the extracted cohort:
WHERE age >= 18 --only extract patients aged over 18 years (i.e. exclude neonatal)
  AND icu_los_days > 1 --exclude any ICU admission with a LOS less than 24 hours
  AND icustay_seq = 1 --only include a patient's first ICU admission within a given hospitalisation
  AND hospital_expire_flag = 0 --only include those who did not die within the hospitalisation  
  AND baseline_eGFR IS NOT null AND baseline_eGFR >= 60 --exclude anyone with baseline eGFR <60 mL/min/1.73m2 (or who was missing baseline creatinine)
;
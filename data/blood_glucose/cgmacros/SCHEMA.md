# CGMacros Schema (Working Draft)

This schema is based on the downloaded dictionary files:

- `DataDictionary_Bio.csv`
- `DataDictionary_CGMacros-00X.csv`
- `DataDictionary_Gut_Health_Test.csv`
- `DataDictionary_Microbes.csv`

The full archive `CGMacros_dateshifted365.zip` is still downloading, so this is a dictionary-driven schema draft.

## Core Tables

### 1) `bio` (participant baseline)
- Grain: one row per participant
- Purpose: baseline demographics and fasting labs
- Fields (23 total; key examples):
  - `Age` (numeric)
  - `Gender` (`F`/`M`)
  - `BMI`
  - `Body weight` (lb)
  - `Height` (inches)
  - `Self-identify` (race/ethnicity category)
  - `A1c PDL (Lab)`
  - `Fasting GLU - PDL (Lab)` (mg/dL)
  - `Insulin`
  - Lipids: `Triglycerides`, `Cholesterol`, `HDL`, `LDL (Cal)`, etc.
  - Fingerstick fields and collection times

### 2) `cgm_minute` (minute-level time series)
- Grain: one row per participant per minute
- Purpose: glucose, activity, and meal events
- Fields (14 total):
  - `Timestamp` (date-shifted)
  - `Libre GL` (mg/dL)
  - `Dexcom GL` (mg/dL)
  - `HR`
  - `Calories (Activity)`
  - `Mets`
  - Meal/event fields:
    - `Meal Type` (`Breakfast`/`Lunch`/`Dinner`)
    - `Calories`, `Carbs`, `Protein`, `Fat`, `Fiber`
    - `Amount Consumed` (%)
    - `Image Path`

### 3) `gut_health_scores`
- Grain: one row per participant (or per test date if repeated)
- Purpose: summarized gut-health scores
- Fields: 22 ordinal score fields
- Value encoding: mostly `1=Not Optimal`, `2=Average`, `3=Good`
- Examples:
  - `Gut Lining Health`
  - `Metabolic Fitness`
  - `Gut Microbiome Health`
  - `Digestive Efficiency`
  - `Inflammatory Activity`

### 4) `microbe_presence`
- Grain: one row per participant (or per sample)
- Purpose: microbiome taxon presence indicators
- Fields: 1,979 indicator columns
- Encoding: binary `0=No`, `1=Yes`
- Each field name is a taxon/species label.

## Keys And Joins

The dictionaries do not explicitly list a canonical participant ID column, but in practice the dataset is likely split by participant (e.g., files named like `CGMacros-00X`).

Use this join strategy once full files are extracted:

1. Identify participant key column(s) in each table/file (`participant_id`, `subject`, filename token, etc.).
2. Build `participant_id` as the conformed key across `bio`, `cgm_minute`, `gut_health_scores`, `microbe_presence`.
3. Join:
   - `bio` 1:many `cgm_minute`
   - `bio` 1:1 (or 1:many by sample date) `gut_health_scores`
   - `bio` 1:1 (or 1:many by sample date) `microbe_presence`

## Data Types (Recommended)

- `Timestamp`: `datetime`
- Glucose/activity/macros: `float`
- Binary indicators: `int8` or `boolean`
- Ordinal gut scores: `int8`
- Categorical labels (`Gender`, `Meal Type`, `Self-identify`): `string/category`

## Notes For Analysis

- `Timestamp` is date-shifted for privacy; use within-participant temporal patterns, not calendar dates.
- Prefer one glucose source at a time (`Libre GL` vs `Dexcom GL`) or harmonize with a sensor-source flag.
- Meal nutrient fields appear event-annotated; minutes without meals likely have null/empty meal columns.
- Treat microbiome table as high-dimensional sparse features (feature selection/regularization recommended).

## Immediate Next Validation (After ZIP Completes)

1. Extract archive and inventory all files.
2. Confirm actual table names and participant key columns.
3. Update this document from draft to final schema with exact primary/foreign keys.

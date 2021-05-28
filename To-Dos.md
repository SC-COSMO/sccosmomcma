# To-Do's for SC-COSMO

## Inputs
1. Design a template for setting-specific inputs
    - Demographic data
        - Population by years of age
        - Number of deaths by years of age
        - Proportion of population leaving in urban setting
        - Proportion of land that is urban
2. Handling of dates
    - First case (when becomes symptomatic)
    - Interventions
    - End of projections
    - Calibration 
        - Lag of infectiousness
        - Lag of confirmation

## Model
1. Severity levels
2. Inflow and non-fatal outflows
3. Time-varying NPIs
    - How to specify timing, intensity and age variations
4. Partial reopening interventions
5. Computing R0 using next generation matrix
6. Compartments to model PIs (i.e., vaccines and treatments)
7. Compartments to model screening 
    - For infectious
    - For recovered
8. Omega, possibility of waning immunity
9. High performance computing (HPC)

## Calibration
1. Multi-parameter calibration
    - Epi parameters with informed priors (?)
    - Beta
    - Intervention effect
    - time-varying CDR
    - hospitalization parameters
        - Probability of hospitalization
        - ICU | hosp
        - All by severity classes
2. Multi-target calibration
    - Cases
    - hospitaliations
    - deaths
3. Hierarchical calibration
    - Simultaneoulsy calibrate all states/counties assuming hyperpriors
4. Identifiability analysis

## Validation
1. Summarize uncertainty quantification (UQ) outputs

## Outputs
1. Calibration outputs
    - Comparison to targets
    - Describing marginal and joint distributions of parameters
2. Summarize uncertainty quantification (UQ) outputs
    - Posterior means
    - Uncertainty bounds (e.g., 75%CI and 95%CI)
3. Web app or Shiny
    - Display
    - Fully interactive for exploration
        - Pull previoulsy run scenarios
        - Live running of the model
    

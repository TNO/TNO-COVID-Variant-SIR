# TNO-COVID-delta2omicron-SIR

Running

for further information on this code and referencing to this code:

(not peer reviewed)
MS ID#: MEDRXIV/2021/268394
MS TITLE: SIR model for assessing the impact of the advent of Omicron and mitigating measures on infection pressure and hospitalization needs


Jan-Diederik van Wees,  Martijn van der Kuip,Sander Osinga, Bart Keijser, David van Westerloo, Maurice Hanegraaf
Maarten Pluymaekers, Olwijn Leeuwenburgh,  Logan Brunner,  Marceline Tutu van Furth




The program has been tested under python 3.9


To run a  forecast:

1. Create a .json configuration file.  example files are provided in this repository in the configs directory
    voc_input_omicronNL.json
    voc_input_omicronSA.json


2. run  voc_main.py and afterwards you can optionally run voc_plot.py (not supported yet)
```json
For voc_main.py
  script path: {your git omicron dir} voc_main.py
  parameters: voc_input_omicronNL.json
  working directory {your git omicron dir}
For voc_plot.py (not supported yet)
  script path: {your git omicron dir} voc_plot.py
  parameters: voc_input_omicronNL.json
  working directory {your git omicron dir}

```

3. formats of input and meaning of parameters

voc_input_omicronNL.json :
```json
{
  "base_dir": "C:\\Users\\weesjdamv\\PycharmProjects\\omicron", # base directory of the run (choose equal to git dir)
    "calcmethod":1,                                     # calculation method 1) SIR,using ts  2) SEIR using tlatent, tinf
  "run_name": "voc_runNL",          # output plotfiles names header
  "n_samples": 100,                 # number of  monte carlo samples
   "startdate": "2021-12-10",      # starting date of the simulation
   "enddate": "2022-4-10",          # end date of the simulation
   "boosterdayx":  [0,  50],        # booster days,
   "boosterdayn":  [0,   17.7],       # mln of  people boostered at the booster days
   "demage": [90, 75, 60, 45, 20, 0] ,  # age cohort definition starting age (from oldest)
  "demn": [ 0,    1.46,    4.51,    8.23, 13.63, 17.7],  # cumulative number of people in age cohorts
  "agegroups_p_vac": [ 0.9,   0.9,    0.8,    0.7, 0.7],  #  vaccination percentage in the age cohorts (one less that demn)
  "w_hosp_age": [8000, 8500, 6000, 2500, 750], # hospitalizations in each age cohort (prior to vaccination)
 "r_lockdowndayx":  [0,  10, 15],     # days at which effective reproduction number (for delta) is specified afterwards exptrapolated from last value 
  "r_lockdown":  [0.88,  0.88, 0.7],   # effective reproduction number values (for delta),
   "comment": "For all variables below, you can either enter a scalar value, or a description of a distribution: [Uniform, min_value, max_value], [(Log)Normal, mean/median, std_dev], [Triangular, min, peak, max]",
  "r_lockdownscale" : 1,  # scaling of the effective reproduction number
  "k_voc" : ["Uniform",0.26, 0.30],  # k-value (natural logarithm of growth rate of variant  cases f(t) )
   "f0_voc": ["Uniform", -2.5, -1.5],  # initial value of log10 value of f(t), so -2 is 1%
 "vac_ratio": ["Uniform",0.8,1.2],   # scaling of booster campaign (boosterdayn)
  "p_vac_start": 0.75,               # fraction of population vaccinated
  "p_booster_start": 0.05,           # fraction of population boostered
  "startinfperc": ["Uniform",0.5,1.5],  # starting value of infected as % of population
  "startsusceptible":  ["Uniform",0.6,0.7], # starting value of susceptible fraction
    "ts": 4.8,    #  ts (days) = 1/gamma for the SIR model (calculation method 1)
  "seir_tlatent": 2,  # latency time (days), for the SEIR model (calculation method 2)
  "seir_tinf": 5,# infectious time (days), for the SEIR model (calculation method 2)
  "ve_trans": 0.4, # vaccine efficacy for transmission 
  "ve_immune_o": 0.2, # natural infected efficacy against immune escape from omicron 
  "ve_immune_hosp_o": ["Uniform",0.4,0.7], # non vaccinated efficacy against hospitalisation from omicron 
  "ve_vac_d": 0.6,# vaccine efficacy against infection with delta 
  "ve_vac_o": 0.1,,#  vaccine efficacy against infection with omicron 
  "ve_booster_d": 0.9, # booster efficacy against infection with delta 
  "ve_booster_o": 0.6,#  booster efficacy against infection with omicron
  "ve_vac_hosp_d": 0.90,, #  vaccine efficacy against hospitalisation for delta
  "ve_vac_hosp_o": ["Uniform",0.85,0.95],#  vaccine efficacy against hospitalisation for omicron
   "ve_booster_hosp_d": 0.975,#  booster efficacy against hospitalisation for delta
  "ve_booster_hosp_o": ["Uniform",0.95, 0.975],#  booster efficacy against hospitalisation for omicron
    "re0_voc": -1,# not supported
  "plot": {
     "plotnames": ["inf", "hosp", "ft", "Rt", "Rtd", "Rtratio", "pu", "pv", "pb"], # plot files created
     "plottitles": ["infected index (1=start)","daily hospitalized index (1=start)", "VOC prevalence", "Reproduction number", "Rtd", "Rtratio", "pu", "pv", "pb"],
      "plotcolors": ["mistyrose","powderblue", "orange","plum","orange","orange","orange","orange","orange"],
      "plothistograms": false,
     "scales": [1.0, 1.0 ,1.0,1.0,1.0,1.0,1.0, 1.0, 1.0, 1.0] # scaling to particular value
  }
}
  
```

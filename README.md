# multifaceted_FD_multifaceted_YLD
Data and code for the study "Multifaceted functional diversity for multifaceted crop yield: towards ecological assembly rules for varietal mixtures"

Two data files are available:

- "CWM_FD.csv" contains one row per experimental plot with community-weighted mean (CWM) and functional diversity (FD) indices computed on the 19 functional traits. 

The first two columns ("genotype_1" & "genotype_2") are the identity of the two genotypes in the plot, which are identical in single-variety plots. The third column ("assoc") refers to the plot type: single-variety ("M") or mixed-variety ("P") plot. Then, all trait CWMs and FDs are reported as "CWM_trait_name" and the "FD_trait_name" , respectively. For single-variety plots, only CWMs are reported but in this case they correspond to unweighted-averaged trait values across replicated measurements within plots. Root trait names are followed by "sem" or "adv" depending if they were measured on seminal or adventitious roots. Reported traits are: "Angle_aer" (Aerial angle, °), "Angle_root" (Root angle, °), "Diam_sem/adv" (mean root diameter, mm), "SRL_sem/adv" (specific root length, m/g), "RTD_sem/adv" (root tissue density, g/cm3), "RBI_sem/adv" (root branching intensity, nb of root tips/cm), "RLD_sem/adv" (root length density, cm root/cm3 soil), "Till_nb" (tiller number per capita), "Ear_bio" (early biomass per capita, g), "SLA" (specific leaf area, m²/kg), "Leaf N" (leaf nitrogen content, %), "Height" (plant height, cm), "Heading" (heading date, Growing Degree Days), and "Maturity" (maturity date, Growing Degree Days).

- "RAW_RYT.csv" contains one row per experimental plot with absolute and relative measures of performance on several agronomic variables.

The first two columns ("genotype_1" & "genotype_2") are the identity of the two genotypes in the plot, which are identical in single-variety plots. The third column ("assoc") refers to the plot type: single-variety ("M") or mixed-variety ("P") plot. Then, all absolute and relative measures of agronomic performance are reported as "RAW_performance_variable_name" and "RYT_performance_variable_name", respectively. For single-variety plots, only absolute performances are reported. Reported performance variable are "GY" (Grain yield, g/m²), "GNb" (Grain number per m²), "SY" (Spike yield, g/m²), "SNb" (Spike number per m²), "BY" (Biomass yield, g/m²), "PY" (Protein yield, g/m²), "TKW" (Thousand kernel weight, g), "SeY" (Semolina yield, %), "GPC" (Grain protein concentration, %), "TW" (Test weight, kg/hL), "RLVA" (Rate of loss of vitreous aspect, %), "YI" (Yellowness index), "GPD" (Grain protein deviation, %).

One R code file is available:

- "Mu_FD_Mu_CY_Analysis.R" contains all statistical analysis performed to produce the results presented in the main text and in the Supplementary Information of the study. It uses "CWM_FD.csv" and "RAW_RYT.csv" files as inputs. 
 

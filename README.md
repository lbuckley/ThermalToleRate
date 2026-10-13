# ThermalToleRate
ThermalToleRate working group aphid project integrating TDT survival estimates with constant temperature TPCs.

This repository adapts the analysis in this repo: <https://github.com/lbuckley/ThermalHistory/>

Citation for corresponding paper: Buckley LB, Huey RB, and Ma C-S. 2025.
How damage, recovery, and repair alter the fitness impacts of thermal stress.
Integrative and Comparative Biology 65:1061–1075. <https://doi.org/10.1093/icb/icaf019>

Synopsis: Integrating estimates of organisms' thermal sensitivity of performance across time series of body temperatures is a common strategy for estimating the fitness implications of climate variability and change.
The approach is an important first step in accounting for organismal sensitivity when estimating climate responses, but the approach omits the influence of past exposure to thermal stress.
We fit data from experiments exposing aphids to stressful, fluctuating temperatures to a function accounting for thermal history to estimate rates of damage, repair, and carryover effects in response to thermal stress.
Our analysis indicates that thermal stress and damage initiate at temperatures beyond thermal optima and considerably below critical thermal limits; that thermal stress can accumulate multiplicatively with the magnitude and duration of exposure; and that considerable and rapid rates of recovery and repair are possible but that rates decline quickly as temperatures diverge from the thermal optima.
Accounting for damage and repair across thermal history is important to accurately assessing the fitness outcomes of climate change and variability.

# Work flow
AphidExpts_TempsFitness: Assembles the temperature and fitness data from the initial experiments. Included survival data that wasn't used in paper.

FecundityEstimatesTDT: Exploration of using TDT from Li et al. 2023 to estimate survival and combining estimates with fecundity data.
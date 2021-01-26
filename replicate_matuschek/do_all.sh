#!/bin/bash -eu

NMC=$1

R --vanilla --args 50 20 $NMC < sim_lmer.R
R --vanilla --args 50 20 $NMC < sim_anova.R

R --vanilla --args 30 10 $NMC < sim_lmer.R
R --vanilla --args 30 10 $NMC < sim_anova.R

R --vanilla < results.R

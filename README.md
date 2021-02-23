# Modeling Indels

This project contains partial (not all) scripts used to build a custom empirical model of HIV-1 genetic mutations. 
To perform this, I coded a custom Markov Chain Monte Carlo sampler used to sample the posterior distribution of chose parameter values.
We needed a custom MCMC sampler in order to incorporate genetic sequence reconstruction steps into the sampler itself. 
Essentially, this code aims to simulate specific mutations (insertions & deletions) that were observed within patient-derived HIV-1 sequence data under a specific model of generation.
In its current state, it is able to reconstruct these insertion/deletion mutations to a moderate/high level of accuracy.  

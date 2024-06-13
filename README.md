# IntegrationLCStim
Integration analysis for Grimm et al. 

Load in the correlation matrices for activity in pre and during stimulation.

Within the code, the community Louvain algorithm is run, creating a clustering of regions. The Participation coefficient is calculated for each cortical ROI from this. The algorithm is repeated multiple times due to the approach's stochasticity. 

Finally, boxplots and statistical testing are outputted!

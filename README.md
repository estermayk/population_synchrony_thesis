# Spatial synchrony and population dynamics in a Scottish blue tit metapopulation
**Code and figures relating to a thesis submitted in partial fulfilment of a masters of science in ecology, evolution and biodiversity**

### Abstract
Spatial synchrony, when the population growth rates of distinct populations connected by dispersal fluctuate simultaneously, can be induced by spatially auto-correlated weather, dispersal, and synchrony of related taxa. It has important consequences for metapopulation persistence, especially where systems are vulnerable to synchronous collapse, and reduces capacity for rescue in bad seasons which are expected to become more frequent under climate change. The processes which determine population growth rate (survival, reproduction etc) are known as demographic rates and can also exhibit spatial synchrony as well as influencing overall population synchrony. However, limited research has been dedicated to these higher resolution influences, and the extent to which they vary synchronously and impact population synchrony is largely unknown. Here I show, in a Scottish blue tit (*Cyanistes caeruleus*) metapopulation, (the synchrony of) population size and demographic rates by constructing an integrated population model (IPM) fitted to 12-years of capture mark recapture, nestbox occupancy and productivity data from a 220km transect. Both population size and growth rate were fairly stable, while population synchrony was low-moderate (intraclass correlation coefficient (ICC) = 0.40). Adult survival, on the other hand, showed exceptionally high synchrony (ICC = 0.90) and was positively correlated with population growth rate, suggesting demographic buffering. While more moderately synchronous, productivity (fledglings per female) also showed a positive relationship here, as expected for a fast pace of life species. These results represent the second attempt to address this question under an integrated Bayesian framework, highlighting a methodology with full propagation of uncertainty and estimation of parameters with no direct data source. These results motivate further investigation into this system to investigate the sensitivity and elasticity of population growth rate synchrony to synchrony of demographic rates, and suggest further expansion on this methodology such as investigating spatial decay within the IPM. 

## FAIR Principles
### Findable
All code and data is available at this repository (data under data, code under R_files and stan_models, and figures under figs) which is publicly visible and searchable. 

## Accessible 
All data and metadata used is stored as excel sheets and csv files contained in the folder titled ‘data’ which is accessible through this repository. 

## Interpolatable
All parameters are named using the convention from the literature (eg phi for survival). In addition, code sections relating to a given parameter are labelled with the full English title for the parameter. Code is annotated and dataframes are intuitively named. 

## Reusable
All code is annotated throughout for understanding and signposting, including breaking plotting down to sections for easy navigation. The first iteration of each phrase is annotated, though annotation may become less thorough where chunks are repeated for a different parameter.


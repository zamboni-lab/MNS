# MNS Toolbox
The Metabolic Network Segmentation (MNS) toolbox contains algorithms employing Markov random fields to identify sites and sequential order of metabolic regulations from large-scale metabolomics datasets and genome-scale metabolic network reconstructions. In version 2 we added the possibility to extract active metabolic modules, i.e. subgraphs, from metabolomics or mixed metabolomics, transcriptomics data. The toolbox is implemented in Matlab and calls specific C++ functions of the OpenGM toolbox (http://hciweb2.iwr.uni-heidelberg.de/opengm/). The toolbox can be run on Windows and Mac operating systems.

## Getting started
### Installation of the toolbox
Copy/unpack the MNS toolbox on your hard drive. **IMPORTANT: Make sure that the path in which you install the MNS toolbox has no blank character otherwise the toolbox won’t run in Matlab.** Start Matlab, go to the folder of the MNS toolbox and initialize the toolbox by executing

### Documentation
Detailed documentation for MNS toolbox v1 can be found in the *doc* subfolder

### Examples
Application examples for different applications are in the *example* subfolder

## Citations
Kuehne A, Mayr U, Sévin DC, Claassen M, Zamboni N (2017) Metabolic network segmentation: A probabilistic graphical modeling approach to identify the sites and sequential order of metabolic regulation from non-targeted metabolomics data. PLoS Comput Biol 13(6): e1005577. [https://doi.org/10.1371/journal.pcbi.1005577](https://doi.org/10.1371/journal.pcbi.1005577)

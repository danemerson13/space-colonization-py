## Note: This readme file and documentation still needs updating (10/1/24)
Results must be added, and some of the recently added code needs to be better commented.

## Overview
This repository contains code to generate 3D vascular models of the human liver. We then simulate CPA loading/perfusion through the vasculature and test out a range of pressure and flow rate boundary conditions to examine how the loading is affected. There are many classes central to the model's data structure, namely colony.py, branch.py, superlobule.py, node.py, attractor.py and segment.py. Models can be created in the manner outlined in main.py, as well as in many of the study specific files: nSLConvergence.py, PQGrid.py, tConvergence.py, and consistency.py. Their respective analysis files for generating plots and results are contained within the analysis folder. 

## Data Availability
Where available, results are made available in the results folder, but many of the files are too large to be hosted on Github. We provide all code necessary to generate these results, but files can also be made available upon request to <danielem@andrew.cmu.edu> or <lkara@cmu.edu>.


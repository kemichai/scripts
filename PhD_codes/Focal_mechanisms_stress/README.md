# Stress inversions for the central Southern Alps, NZ

## Description ##
Number of codes in Python, Matlab and R to reproduce Figure 2
from Michailos et al., 2020 in [Tectonophysics](https://www.sciencedirect.com/science/article/abs/pii/S0040195119303208?via%3Dihub).

### Requirements:
- Python
- Matlab
- R

### Order of codes to use with brief description
1) `quadtree_prep.py`

    Read quakeml catalog of earthquakes with focal mechanisms and 
    prepares input file for the Matlab quadtree code.

2) `quadtree_NZ.m`

    Run quadtree. 

3) `read_quadtree_output.py`

    Reads quadtree output and creates input for rstress codes.


4) `drive_inversion.q`

    Run rstress inversions.  

5) `read_rstress_outputs.py`

    Read rstress outputs and store a file with each bin's
    details as well as a file to plot the bowties with GMT.


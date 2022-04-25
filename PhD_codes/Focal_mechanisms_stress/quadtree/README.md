# Stress inversions for the central Southern Alps, NZ

## Description ##
Number of codes in Python, Matlab and R to reproduce Figure 2
from Michailos et al., 2020 (Tectonophysics).

### Requirements:
- Python
- Matlab
- R

### Order of codes to use with brief description
1) _quadtree_prep.py_

    Read quakeml catalog of earthquakes with focal mechanisms and 
    prepares input file for the Matlab quadtree code.

2) _quadtree_NZ.m_

    Run quadtree. 

3) _read_quadtree_output.py_

    Reads quadtree output and creates input for rstress codes.


4) _drive_inversion.q_

    Run rstress inversions.  

5) _read_rstress_outputs.py_

    Read rstress outputs and store a file with each bin's
    details as well as a file to plot the bowties with GMT.


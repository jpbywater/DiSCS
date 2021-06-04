# DiSCS
The Differential Segmentation of Categorical Sequences (DiSCS) algorithm 

DiSCS is a machine-learning algorithm that takes a categorial sequence and splits it into segments so that each of the segments is as different from its immediate neighbors as possible. In effect, it finds the most distinct segments within a sequence of categories.

### Using DiSCS
The DiSCS functions are written in Python. To use them:
 - copy the file `DiSCSfunctions.py` into your project folder
 - import the DiSCS functions with:
```
 from DiSCSfunctions import *
```
 - check the python packages `numpy`, `scipy`, and `sklearn` are installed.
 - segment a sequence with the function:
 ```
find_best_cut_positions(sequence)
```

For example, the code:
```
sequence = ['A','A','A','A','B','B','B','C','C','C','C','C']
cut_positions, stat = find_best_cut_positions(sequence)
print(cut_postions)
```
will output:
```
[0, 4, 7, 12]
```
This indicates that the DiSCS algorithm found the most distinct segments when the sequence was cuts at position 4 (i.e. between the As and the Bs) and position 7 (i.e. between the Bs and the Cs).

### DiSCS Performance

This [peer-reviewed paper (AIED2021)](https://github.com/jpbywater/DiSCS/blob/main/DiSCS%20_final_camera_ready.pdf) reports detailed DiSCS performance data from a simulation study. The script used for the simulation study is:
```
AIED2021_simulation_script.py
```
Within the `DiSCSfunctions.py` file there are several parameters that can be changed to adjust DiSCS performance (see paper for details). Many of these parameters involve trade-offs. For example, the more segments that DiSCS searches over, the longer DiSCS will take to find a result.  

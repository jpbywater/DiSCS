## The Differential Segmentation of Categorical Sequences (DiSCS) algorithm 
DiSCS is an algorithm that takes a categorial sequence and splits it into segments so that each of the segments is as different from its immediate neighbors as possible. In effect, it finds the most distinct segments within a sequence of categories.

For example:
![graphical example of segmentation](https://github.com/jpbywater/DiSCS/blob/main/example_img.png)

### DiSCS Demonstration site
Click [here](https://jpbywater.pythonanywhere.com/) for an interactive visual demo of DiSCS segmenting sequences.

### Using DiSCS
The DiSCS functions are written in Python. To use them:
 - copy the file `DiSCSfunctions.py` into your project folder
 - check the python packages `numpy`, `scipy`, and `sklearn` are installed in your virtual environment
 - in your python script, import the DiSCS functions with:
```
 from DiSCSfunctions import *
```
 - in your python script, segment a sequence with the function:
 ```
find_best_cut_positions(sequence)
```

For example, the code:
```
sequence = ['A','A','A','A','B','B','B','A','A','A','A','A']
cut_positions, stat = find_best_cut_positions(sequence)
print(cut_postions)
```
will output:
```
[0, 4, 7, 12]
```
This indicates that the DiSCS algorithm found the most distinct segments when the sequence was cut at positions 4 and 7 (i.e. at the boundaries of the As and the Bs).

### DiSCS Performance

This [peer-reviewed paper (AIED2021)](https://link.springer.com/chapter/10.1007/978-3-030-78292-4_8) reports detailed DiSCS performance data from a simulation study. The script used for the simulation study is:
```
AIED2021_simulation_script.py
```
Within the `DiSCSfunctions.py` file there are several parameters that can be changed to adjust DiSCS performance (see paper for details). Many of these parameters involve trade-offs. For example, the more segments that DiSCS searches over, the longer DiSCS will take to find a result.  

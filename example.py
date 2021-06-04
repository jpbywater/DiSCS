from DiSCS.DiSCSfunctions import *

# example from AIED 2021 paper
action_list = ['R','R','R','R','G','G','G','Y','R','R','R','R','R','R','R','G','G','G','Y','R','R','R','R']

# specify maximum number of segments to search for
max_sections = 10

# do DiSCS segmentation
print("Working on DiSCS segmentation. Please wait... this might take a few minutes")
cut_positions, stat = find_best_cut_positions(action_list, max_sections)

# report positions in sequence where cut into segments
print("Cut Positions:", cut_positions)

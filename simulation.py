from simulation_study_AIED2021.DiSCSfunctions import *
from simulation_study_AIED2021.MarkovModels import *
from simulation_study_AIED2021.additionalfunctions import *
import numpy as np
import datetime

### PROCESS ###
# 0) Specify default parameters for sequences
number_of_states_values = [2, 3, 4]
number_of_actions_in_each_state_values = [2, 4, 6]
p_between_state_values = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
sequence_length_values = [25, 50, 75, 100, 150, 200, 300, 400]
grainsize_values = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
p_common_action_values = [0.05, 0.1, 0.15, 0.2, 0.25]

number_of_states = number_of_states_values[1]
number_of_actions_in_each_state = number_of_actions_in_each_state_values[0]
p_between_state = p_between_state_values[0]
sequence_length = sequence_length_values[5]
grainsize=grainsize_values[8]
p_common_action = 0

# 1) Specify what to loop over
for p_between_state in [0.05, 0.1]:
    max_sections = find_max_sections(sequence_length, p_between_state)
    for x in range(20): # repeat 20 times at each value
        start_time = datetime.datetime.now()
        # 2) Generate sequences
        states, actions = generate_sequence(number_of_states, number_of_actions_in_each_state, p_between_state, sequence_length, p_common_action)
        print("State names:", states)
        print("Action names:", actions)

        # 3) Run DiSCS algorithm
        cut_positions, stat = find_best_cut_positions(actions, max_sections, grainsize)
        clustering_stat, number_of_clusters, clusterIDs, cluster_centers = cluster_blocks(actions, cut_positions)
        D_states = DiSCS_states(cut_positions, clusterIDs)
        #print("DiSCS states:", D_states)

        # 4) Find association between original states and DiSCS states using Cramers V
        V = cramers_v(states, D_states)
        print("Cramer's V =", V)

        # 5) Save output.
        end_time = datetime.datetime.now()
        duration = end_time - start_time
        print("Duration", duration)
        file_object = open('datalog.csv', 'a')
        new_text = str(number_of_states) + "," + str(number_of_actions_in_each_state) + "," + str(p_between_state) + "," + str(sequence_length) + "," + str(V) + "," + str(max_sections) + "," + str(duration) + "," + str(grainsize) + "," + str(p_common_action) + "\n"
        file_object.write(new_text)
        file_object.close()

### END PROCESS ###
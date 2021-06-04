from DiSCS.DiSCSfunctions import *
from scipy import stats
import pandas as pd
import numpy as np
import random
import datetime


# FUNCTIONS #
def generate_sequence(number_of_states=2,
                      number_of_actions_in_each_state=3,
                      p_between_state=0.1,
                      sequence_length=20,
                      p_common_action = 0):

    state_names = ["A", "B", "C", "D", "E", "F", "G"][:number_of_states]
    state_id = 0
    action_id = 0
    output_state_names = []
    output_action_names = []

    for i in range(sequence_length):
        # save current state and action to output lists
        output_state_names.append(state_names[state_id])
        if random.random() < p_common_action:
            output_action_names.append('common')
        else:
            output_action_names.append(state_names[state_id].lower() + str(action_id))
        # transition state
        if random.random() < p_between_state:
            # transition to another state randomly (equally likely excluding current).
            state_ids = list(range(number_of_states))
            state_ids.remove(state_id)
            state_id = random.choice(state_ids)
        # whether transition state or not, randomly select new action id (equally likely including current)
        action_ids=list(range(number_of_actions_in_each_state))
        action_id=random.choice(action_ids)
    return output_state_names, output_action_names


def DiSCS_states(cut_positions, clusterIDs, labels=["Z", "Y", "X", "W", "V", "U", "T", "S"]):
    output=[]
    for i in range(len(cut_positions)-1):
        output.extend(labels[clusterIDs[i]]*(cut_positions[i+1]-cut_positions[i]))
    return output


def cramers_v(states, DiSCS_states):
    df = pd.DataFrame(columns=['states', 'DiSCS_states'])
    df['states'] = states
    df['DiSCS_states'] = DiSCS_states
    contingency_table = pd.crosstab(df['states'], df['DiSCS_states'])
    print(contingency_table)
    chi2, p, dof, expected = stats.chi2_contingency(contingency_table, correction=False)
    n = contingency_table.sum().sum()
    phi2 = chi2/n
    r,k = contingency_table.shape
    V = np.sqrt(phi2/min((k-1),(r-1)))
    return V


def find_max_sections(n, p):
    max = int((n*p) + (3 * np.sqrt(n*p*(1-p)))) #expected number of extions plus 3 stdevs
    if max > n:
        max = n
    return max
# END FUNCTIONS #


# SIMULATION #
# 0) Specify default parameters for sequences.
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

# 1) Specify what to loop over. Change the for loop variables to test the different parameters above
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

# END SIMULATION #
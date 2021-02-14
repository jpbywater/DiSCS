from scipy import stats
import pandas as pd
import numpy as np


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
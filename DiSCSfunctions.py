import numpy as np
from scipy import stats, special
import random
from operator import itemgetter
import itertools
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def find_best_cut_positions(input_list, max_sections=7, grainsize=1):
    max_sections = min(max_sections, len(input_list)-1)
    grainsize = min(max_sections/2, grainsize)
    cut_positions_gen, stat_gen = find_best_cut_positions_genetic(input_list, max_sections, grainsize)
    cut_positions_dyn, stat_dyn = find_best_cut_positions_dynamic(input_list, max_sections, grainsize)
    if stat_dyn>stat_gen:
        cut_positions = cut_positions_dyn
        stat = stat_dyn
    else:
        cut_positions = cut_positions_gen
        stat = stat_gen
    return cut_positions, stat


def get_array_comparison_stat(array1, array2, fn="two_props"):
    if fn == "sum_sq":
        diff=array1 - array2
        return np.square(diff).sum()
    if fn == "sum_abs":
        diff = array1 - array2
        return np.absolute(diff).sum()
    if fn == "chi2_homo":
        table = np.vstack((array1, array2))
        chi2, p, dof, ex = stats.chi2_contingency(table+1)
        return chi2
    if fn == "two_props":
        array1 = array1 + 1
        array2 = array2 + 1
        props1 = array1 / array1.sum()
        props2 = array2 / array2.sum()
        pooledp = (array1 + array2) / (array1.sum() + array2.sum())
        variance = pooledp * (1 - pooledp) * ((1 / array1.sum()) + (1 / array2.sum()))
        se = np.sqrt(variance)
        abs_zscore = abs(props1 - props2) / se
        stat = abs_zscore.sum() / len(abs_zscore)
        return stat


def find_best_cut_positions_brute_force(input_list, max_sections=3, grainsize=1):
    n = len(input_list)
    counts = make_counts_array(input_list)
    ft = []
    cpt = []
    for c in range(2, max_sections+1):
        # print("working on", c, "sections...")
        if n>=c:
            fs = []
            cps = []
            for combos in itertools.combinations(list(range(1, n)), c - 1):
                cut_positions = [0] + list(combos) + [n]
                f=0
                for i in range (len(cut_positions)-2):
                    array1 = counts[cut_positions[i], cut_positions[i+1]]
                    array2 = counts[cut_positions[i+1], cut_positions[i+2]]
                    g = get_array_comparison_stat(array1, array2)
                    f = combine_g(f, i, g, grainsize)
                fs.append(f)
                cps.append(cut_positions)
            ft.append(np.array(fs).max())
            cpt.append(cps[np.array(fs).argmax()])
    best_f = np.array(ft).max()
    best_cuts = cpt[np.array(ft).argmax()]
    sections = []
    for x in range (len(best_cuts)-1):
        start=best_cuts[x]
        end=best_cuts[x+1]
        sections.append(input_list[start:end])
    return best_cuts, best_f


def find_best_cut_positions_genetic(input_list, max_sections=5, grainsize=1, mate_prob=0.5, mut_prob=0.3):
    def get_fitness(cuts, counts, grainsize):
        f = 0
        for i in range(len(cuts) - 2):
            array1 = counts[cuts[i], cuts[i + 1]]
            array2 = counts[cuts[i + 1], cuts[i + 2]]
            g = get_array_comparison_stat(array1, array2)
            f = combine_g(f, i, g, grainsize)
        return f
    n = len(input_list)
    counts = make_counts_array(input_list)
    best_cuts = []
    for c in range(2, max_sections+1):
        # print("working on", c, "sections...")
        if n>=c:
            # MAKE ORIGINAL POPULATION
            # set population size (reduce if small data set)
            pop_size = 5000
            cut_combos = special.comb(n - 1, c - 1, exact=True)
            if round(cut_combos/2) < pop_size:
                pop_size = round(cut_combos/2)
            # create original parents
            population_cuts = []
            while len(population_cuts) < pop_size:
                cut_positions = [0] + sorted(random.sample(range(1,n), k=c-1)) + [n]
                if cut_positions not in population_cuts:
                    population_cuts.append(cut_positions)
            population = []
            for cuts in population_cuts:
                f = get_fitness(cuts, counts, grainsize)
                population.append([cuts,f])
            population = sorted(population, key=itemgetter(1), reverse=True)
            # START EVOLUTION
            for generation in range(500):
                # IDENTIFY MATES
                mates = []
                i=0
                while len(mates) < 2:
                    if random.random() < mate_prob:
                        mates.append(population[i])
                        if i < len(population) - 1:
                            i = i + 1
                    else:
                        if i < len(population) - 1:
                            i = i + 1
                        else:
                            mates.append(population[i])
                # MAKE OFFSPRINGS
                # combine mates
                number_of_offspring = 2
                offspring = []
                while len(offspring) < number_of_offspring:
                    splice_position = random.randint(2, c)
                    # need to figure out while parent to come first based on which has lower number at splice point
                    if mates[0][0][splice_position] < mates[1][0][splice_position]:
                        cuts = mates[0][0][:splice_position] + mates[1][0][splice_position:]
                    else:
                        cuts = mates[1][0][:splice_position] + mates[0][0][splice_position:]
                    offspring.append(cuts)
                # mutate offspring cuts
                for cuts in offspring:
                    if random.random() < mut_prob:
                        for i in range(1, c):
                            nudge_down_max = int((cuts[i] - cuts[i-1]) / 2)
                            nudge_up_max = int((cuts[i+1] - cuts[i]) / 2)
                            cuts[i] = cuts[i] + random.randint(-nudge_down_max, nudge_up_max)
                # check offspring not already in population
                population_cuts = []
                for member in population:
                    population_cuts.append(member[0])
                offspring_to_add = []
                for cuts in offspring:
                    if cuts not in population_cuts:
                        offspring_to_add.append(cuts)
                # REPLACE UNFIT POPULATION WITH OFFSPRING
                population = population[:pop_size-len(offspring_to_add)]
                for cuts in offspring_to_add:
                    f = get_fitness(cuts, counts, grainsize)
                    population.append([cuts, f])
                population = sorted(population, key=itemgetter(1), reverse=True)
            best_cuts.append(population[0]+[c])
    bestist = sorted(best_cuts, key=itemgetter(1), reverse=True)[0]
    sections = []
    for x in range(len(bestist[0])-1):
        start = bestist[0][x]
        end = bestist[0][x+1]
        sections.append(input_list[start:end])
    return bestist[0], bestist[1]


def find_best_cut_positions_dynamic(input_list, max_sections=5, grainsize=1):
    n = len(input_list)
    counts = make_counts_array(input_list)
    # initialize arrays
    # A stores best value of the optimization function for [first d rows, with c cuts]
    A = np.zeros((n, max_sections+1), dtype=object)
    # P stores previous cut position that A used (e.g. 0,1,2 v 3,4 => 3)
    P = np.zeros((n, max_sections+1), dtype=object)
    # first two columns are zero/empty since 0 cuts/sections is meaningless and 1 cut/section has no comparison.
    for c in range(2, max_sections+1):
        # print("working on", c, "sections...")
        for d in range(n):
            if d+1>=c:
                tmpf = []
                tmpi = []
                for i in range(c-2, d):
                    array1 = counts[P[i, c-1], i+1]
                    array2 = counts[i+1, d+1]
                    g = get_array_comparison_stat(array1, array2)
                    f = combine_g(A[i, c-1], c-2, g, grainsize)
                    tmpf.append(f)
                    tmpi.append(i)
                A[d,c] = np.array(tmpf).max()
                P[d,c] = tmpi[np.array(tmpf).argmax()]+1
    best_f = A[n-1,].max()
    cut_positions = [n]
    row = n-1
    section = A[row,].argmax()
    while section > 0:
        cut_pos = P[row, section]
        cut_positions.append(cut_pos)
        row = cut_pos-1
        section = section -1
    cut_positions_asc = list(reversed(cut_positions))
    sections = []
    for x in range (len(cut_positions_asc)-1):
        start = cut_positions_asc[x]
        end = cut_positions_asc[x+1]
        sections.append(input_list[start:end])
    return cut_positions_asc, best_f


def cluster_blocks(input_list, cut_positions):
    all_block_proportions_list = []
    counts = make_counts_array(input_list)
    for i in range(len(cut_positions)-1):
        block_counts = counts[cut_positions[i], cut_positions[i+1]]
        block_proportions = block_counts/block_counts.sum()
        all_block_proportions_list.append(block_proportions)
    all_block_proportions_array = np.array(all_block_proportions_list)
    max_allowed_clusters = len(cut_positions)
    clustering_output = []
    for j in range(2, max_allowed_clusters-1):
        clusterer = KMeans(n_clusters=j)
        for i in range(100):
            clusterIDs = clusterer.fit_predict(np.array(all_block_proportions_array))
            cluster_centers = clusterer.cluster_centers_
            clustering_score = silhouette_score(all_block_proportions_array, clusterIDs, metric='euclidean')
            clustering_output.append([clustering_score, j, clusterIDs, cluster_centers])
    best_clustering_output = sorted(clustering_output, key=itemgetter(0), reverse=True)[0]
    clustering_stat = best_clustering_output[0]
    number_of_clusters = best_clustering_output[1]
    clusterIDs = best_clustering_output[2]
    cluster_centers = best_clustering_output[3]
    return clustering_stat, number_of_clusters, clusterIDs, cluster_centers


def cluster_blocks_multi(list_of_pairs_of_actions_and_cut_positions, max_allowed_clusters = 7):
    # get all categories
    all_actions = []
    for pair in list_of_pairs_of_actions_and_cut_positions:
        all_actions.extend(pair[0])
    all_categories = get_category_names_from_list(all_actions)
    print("ALL:", all_actions)
    print("lenALL:", len(all_actions))
    print("col headings", all_categories)

    all_block_proportions_list = []
    for pair in list_of_pairs_of_actions_and_cut_positions:
        input_list = pair[0]
        cut_positions = pair[1]
        sparse_array = np.zeros((len(input_list), len(all_categories)), dtype=int)
        for index, item in enumerate(input_list):
            sparse_array[index, all_categories.index(item)]=1
        for i in range(len(cut_positions)-1):
            segment = sparse_array[cut_positions[i]:cut_positions[i+1], :]
            block_proportions = segment.sum(axis=0)/segment.sum()
            all_block_proportions_list.append(block_proportions)
    all_block_proportions_array = np.array(all_block_proportions_list)
    clustering_output = []
    for j in range(2, max_allowed_clusters-1):
        clusterer = KMeans(n_clusters=j)
        for i in range(100):
            clusterIDs = clusterer.fit_predict(all_block_proportions_array)
            cluster_centers = clusterer.cluster_centers_
            clustering_score = silhouette_score(all_block_proportions_array, clusterIDs, metric='euclidean')
            clustering_output.append([clustering_score, j, clusterIDs, cluster_centers])
    best_clustering_output = sorted(clustering_output, key=itemgetter(0), reverse=True)[0]
    clustering_stat = best_clustering_output[0]
    number_of_clusters = best_clustering_output[1]
    clusterIDs = best_clustering_output[2]
    cluster_centers = best_clustering_output[3]
    return clustering_stat, number_of_clusters, clusterIDs, cluster_centers, all_categories


def combine_g(f, i, g, grainsize=0):
    if grainsize < 0:
        grainsize = 0
    output = ((f * (i + grainsize)) + g) / (i + 1 + grainsize)
    return output


def get_category_names_from_list(input_list):
    return list(sorted(set(input_list)))


def make_counts_array(input_list):
    # make a 2D array with each entry containing an 1D array with the total counts for each category
    # where the totals are [from row, to (but not including) row]
    n = len(input_list)
    column_headings = get_category_names_from_list(input_list)
    sparse_array = np.zeros((len(input_list), len(column_headings)), dtype=int)
    for index, item in enumerate(input_list):
        sparse_array[index, column_headings.index(item)]=1
    counts_array = np.zeros((n, n + 1), dtype=object)
    for start in range(0, n):
        for stop in range(0, n + 1):
            if stop > start:  # only do this if it makes sense
                counts_array[start, stop] = counts_array[start, stop - 1] + sparse_array[stop - 1]
    return counts_array

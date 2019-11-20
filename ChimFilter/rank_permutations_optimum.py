from ChimFilter.code_as_binary import code_as_binary
from ChimFilter.functions import bidirectional
from collections import Counter

def rank_permutations_optimum(aln):
    # aln: the output of the code_alignment function
    coded = aln[0]
    coverage = aln[3]
    names = aln[2]

    # Identify the permutation that maximises the difference between parent sequences:
    permutations, distances, maxBreaks, biBreaks, stats, parentcov, parentNames = [], [], [], [], [], [], []

    for k, ref in enumerate(coded):
        this_permutation = [tuple([0 if s[i] == y else 1 for i, y in enumerate(ref)]) for s in coded]

        dist = [sum(x) for x in this_permutation]
        parent1 = [i for i, d in enumerate(dist) if d == max(dist)][0]  # id of the second parent in the list
        parent0_name = names[k]
        parent1_name = names[parent1]

        verbose = False
        if verbose:
            print("Maximum distance: " + str(max([sum(x) for x in this_permutation])))
            print("coverage of parent 0: " + str(coverage[k]))
            print("coverage of parent 1: " + str(coverage[parent1]))
            print("coverage summed: " + str(coverage[k] + coverage[parent1]))
            print("-------------------")

        cov = coverage[k] + coverage[parent1]
        ds = max([sum(x) for x in this_permutation])
        mb = max([sum(x) for x in code_as_binary(this_permutation)])
        bb = bidirectional(this_permutation)[0]

        permutations.append(this_permutation)
        stats.append((ds, cov, mb, bb))
        distances.append(ds)
        maxBreaks.append(mb)
        biBreaks.append(bb)
        parentcov.append(cov)
        parentNames.append([parent0_name, parent1_name])

    sorted_permutations = [x for _, x in sorted(zip(stats, permutations), key=lambda x: (-x[0][0], -x[0][1], x[0][2], x[0][3]))]
    sorted_stats = sorted(stats, key=lambda x: (-x[0], -x[1], x[2], x[3]))
    sorted_names = [x for _, x in sorted(zip(stats, parentNames), key=lambda x: (-x[0][0], -x[0][1], x[0][2], x[0][3]))]

    # We can remove nearly half the solutions, since each "good" solution should appear twice:
    solution_id = []
    for i, s in enumerate(sorted_stats):
        solution_id.append(tuple(sorted(sorted_names[i]) + [x for x in s]))
    redundancy = Counter(solution_id)

    keptPermutations, keptStats, keptNames, isRemoved = [], [], [], []
    for i, s in enumerate(solution_id):
        if redundancy[s] > 1 and s not in isRemoved:
            isRemoved.append(s)
            keptPermutations.append(sorted_permutations[i])
            keptStats.append(sorted_stats[i])
            keptNames.append(sorted_names[i])
        elif redundancy[s] == 1:
            keptPermutations.append(sorted_permutations[i])
            keptStats.append(sorted_stats[i])
            keptNames.append(sorted_names[i])


    # Check if best distance is also best coverage:
    if sorted_stats[0][1] == max(parentcov):
        print("Best distance is best coverage")

    return([keptPermutations, keptStats, keptNames])
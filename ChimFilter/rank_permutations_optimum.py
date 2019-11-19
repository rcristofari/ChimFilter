def rank_permutations_optimum(coded, coverage, leaveOneOut=False):
    # Identify the permutation that maximises the difference between parent sequences:
    permutations, distances, maxBreaks, biBreaks, stats, parentcov = [], [], [], [], [], []

    for k, ref in enumerate(coded):
        this_permutation = [tuple([0 if s[i] == y else 1 for i, y in enumerate(ref)]) for s in coded]

        dist = [sum(x) for x in this_permutation]
        parent1 = [i for i, d in enumerate(dist) if d == max(dist)][0]

        # print("Maximum distance: " + str(max([sum(x) for x in this_permutation])))
        # print("coverage of parent 0: " + str(coverage[k]))
        # print("coverage of parent 1: " + str(coverage[parent1]))
        # print("coverage summed: " + str(coverage[k] + coverage[parent1]))
        # print("-------------------")

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

    sorted_permutations = [x for _, x in
                           sorted(zip(stats, permutations), key=lambda x: (-x[0][0], -x[0][1], x[0][2], x[0][3]))]
    sorted_stats = sorted(stats, key=lambda x: (-x[0], -x[1], x[2], x[3]))

    # Check if best distance is also best coverage:
    if sorted_stats[0][1] == max(parentcov):
        print("Best distance is best coverage")

    return ([sorted_permutations, sorted_stats])
def splitlist(x, by):
    levels = [x for x in set(by)]
    levels.sort()
    res = []
    for l in levels:
        this_level = []
        for i, val in enumerate(x):
            if by[i] == l:
                this_level.append(val)
        res.append(this_level)
    return([levels, res])


# How many breakpoints are bidirectional?
def bidirectional(permutation):
    c = code_as_012(permutation)
    ambiguous = 0
    bidir = []
    for i, x in enumerate(c[0]):
        this_position = [x[i] for x in c]
        if 1 in this_position and 2 in this_position:
            # print(this_position)
            ambiguous += 1
            bidir.append(i)
    return([ambiguous, bidir])
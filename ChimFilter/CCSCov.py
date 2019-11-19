def CCSCov(file):
    reads, assignment = [], []
    with open(file) as handle:
        names = next(handle)
        lines = []
        for h in handle:
            data = h.strip('\n').split(',')
            reads.append(data[0])
            assignment.append([int(round(float(x))) for x in data[1:]])

    names = names.strip('\n').split(',')[1:]
    ndict = dict(zip(names, [0] * len(names)))
    belldict = dict(zip(names, [[] for x in range(len(names))]))

    # Assign each subread's SMRTbell to its cluster:
    for i, r in enumerate(reads):
        SMRTbell = r.split('/')[1]
        clust = [index for index, x in enumerate(assignment[i]) if x == 1]
        if clust:
            ndict[names[clust[0]]] += 1
            if SMRTbell not in belldict[names[clust[0]]]:
                belldict[names[clust[0]]].append(SMRTbell)
    trueCov = dict([(x, len(belldict[x])) for x in belldict.keys()])

    return (trueCov)

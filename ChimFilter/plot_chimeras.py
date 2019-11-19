def plot_chimeras(permutation, bidirectional=None, w=20, h=10, s=60):
    plt.rcParams['figure.figsize'] = [w, h]

    # Plot the permutation:
    spermutation = sorted(permutation)
    for i, c in enumerate(spermutation):
        for p in range(len(c)):
            if c[p] == 1:
                col = 'r'
            else:
                col = 'g'
            plt.scatter(p, len(spermutation) - i, c=col, marker='s', s=s)

    # Identify all state changes:
    known_breaks = []
    for i, c in enumerate(permutation):
        for p in range(len(c) - 1):
            if c[p] != c[p + 1]:
                if p not in known_breaks:
                    known_breaks.append(p)

    known_breaks.sort()

    for i, k in enumerate(known_breaks):
        if bidirectional:
            if i not in bidirectional:
                plt.axvline(k + .5, c='k', ls='--', lw=3)
            else:
                plt.axvline(k + .5, c='r', ls='--', lw=3)
        else:
            plt.axvline(k + .5, c='k', ls='--', lw=3)

    plt.show()
# This function takes an alignment, and returns a list of SNPs, a coded version (0/1), and a windows averaging.
# The SNP calling does not obey to any model. As these are LAA sequences, I consider them "correct" SNP wise.
# Of course this is not true, and would need to be handled somehow using the error scores returned by LAA.
from collections import Counter

def code_alignment(sequences, names, covdict, minCov=0, ignoreSingle=True):
    # sequences is a list of DNA sequences as str
    # names is a list of sequence names, as in the coverage dictionary
    # covdict is a dictionary containing sequence name and coverage

    # Distribution of unique sequences:
    counter = Counter(sequences)
    counts = []
    for x in set(sequences):
        counts.append(counter[x])
    counts.sort()
    counts.reverse()

    ## ENCODE AS 0/1 COMPARED TO A REFERENCE SEQUENCE:
    # The reference is arbitrarily the most common sequence:
    ref = counter.most_common(1)[0][0]

    # Length of the sequence:
    L = len(ref)

    # A list that will hold L lists, each with the group membership of each sequence at that position:
    belonging = []

    # A list of the SNP positions
    SNPs = []

    # We exclude:
    # - singleton SNPs (documented in only one sequence) {<<- this should be from CCS...? but alignment problem}
    # Singletons are only considered as such if there are at least 5 sequences (otherwise sampling error)
    # - InDels (too hard to score)
    # - positions with more than 2 alleles

    # Filter out sequences that do not pass the coverage threshold:
    keptSequences, keptNames = [], []
    for i, s in enumerate(sequences):
        if covdict[names[i]] >= minCov:
            keptSequences.append(s)
            keptNames.append(names[i])

    # Scanning along the sequence:
    for l in range(L):

        # We examine the position "vertically":
        this_position = [s[l] for s in keptSequences]

        # We ignore positions with an indel:
        if '-' not in this_position:
            # We ignore positions with more than 2 alleles:
            if len(set(this_position)) == 2:
                # We ignore positions with singleton and >= 5 sequences: ###THIS HERE CAN BE REFINED A LOT
                c = Counter(this_position)
                if len(this_position) >= 5 and 1 not in [x[1] for x in c.most_common()]:
                    # If it passed all the criteria:
                    belonging.append([0 if s[l] == ref[l] else 1 for s in keptSequences])
                    SNPs.append(l)

    if not SNPs:
        return(None)

    nSNPs = len(SNPs)
    coded = [[belonging[p][s] for p in range(nSNPs)] for s in range(len(keptSequences))]

    return (coded, SNPs, keptNames, [covdict[n] for n in keptNames])
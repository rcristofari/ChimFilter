# A function based on code_as_binary,
# but that distinguishes between moving *away* from reference or *towards* reference

def code_as_012(permutation, wSize=1, ofile=None):
    # If desired, we perform window-averaging:
    if wSize > 1:
        wavg = []
        # wSize = 2 # This is now passed as an argument in the function call
        nWindows = len(permutation[0]) // wSize

        for s in range(len(permutation)):
            this_seq = []
            for w in range(nWindows):
                avg = sum(permutation[s][w * wSize:w * wSize + wSize]) // wSize
                this_seq.append(avg)
            wavg.append(this_seq)

        # Clean out the full identity windows (everywhere there is no polymorphism after window averaging)
        for i, c in enumerate(wavg[0]):
            if all(x == c for x in [y[i] for y in wavg]):
                for x in wavg:
                    x[i] = 2
        permutation = wavg[:]

    #####################################################################################################

    # Identify all state changes:
    breakpoints = []
    for i, c in enumerate(permutation):
        for p in range(len(c) - 1):
            if c[p] != c[p + 1]:
                if p not in breakpoints:
                    breakpoints.append(p)
    breakpoints.sort()

    #####################################################################################################

    # Recode the sequences as binary breakpoint sequences (1 if they break away from ref, 2 towards ref,
    # 0 if they don't)

    chimseq = []
    this_ref = permutation[0]

    for p in permutation:
        this_chimseq = []
        for b in breakpoints:
            if p[b] != p[b + 1] and p[b] == this_ref[b]:
                this_chimseq.append(1)
            elif p[b] != p[b + 1] and p[b] != this_ref[b]:
                this_chimseq.append(2)
            elif p[b] == p[b + 1]:
                this_chimseq.append(0)
            else:
                warning("Unplanned case")
        chimseq.append(this_chimseq)

    # For writing out the characters, everything needs to be coded as "0" or "1".
    # Non ambigious positons go from 0/2 to 0/1
    # Ambiguous positions are split as two separate 0/1 positions (the 0/1 state first, the 0/2 state next)
    all_positions = []

    for i, x in enumerate(chimseq[0]):
        this_position = [x[i] for x in chimseq]
        if 1 not in this_position and 2 in this_position:
            k = [1 if x == 2 else 0 for x in this_position]
            all_positions.append(k)
        elif 1 in this_position and 2 in this_position:
            k1 = [1 if x == 1 else 0 for x in this_position]
            k2 = [1 if x == 2 else 0 for x in this_position]
            all_positions.append(k1)
            all_positions.append(k2)
        else:
            all_positions.append(this_position)

    haplotypes = []
    # Convert back to haplotype format:
    for i, j in enumerate(all_positions[0]):
        this_haplotype = [x[i] for x in all_positions]
        haplotypes.append(this_haplotype)

    if ofile:
        nSeq = len(haplotypes)
        nChar = len(haplotypes[0])
        with open(ofile, 'w') as handle:
            handle.write(str(nSeq) + ' ' + str(nChar) + '\n')
            for i, s in enumerate(haplotypes[:len(haplotypes)]):
                handle.write('chimera_' + str("%02d" % i) + '\t' + ''.join([str(x) for x in s]) + '\n')

    return (chimseq)

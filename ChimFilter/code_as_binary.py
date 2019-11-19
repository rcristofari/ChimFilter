# A function coding the sequence along its breakpoints (1 if the breakpoint is realised, 0 otherwise)

def code_as_binary(permutation, wSize=1, ofile=None):
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

    # Identify all state changes:
    breakpoints = []
    for i, c in enumerate(permutation):
        for p in range(len(c) - 1):
            if c[p] != c[p + 1]:
                if p not in breakpoints:
                    breakpoints.append(p)
    breakpoints.sort()

    # Recode the sequences as binary breakpoint sequences (1 if they break, 0 if they don't)
    chimseq = []
    for p in permutation:
        this_chimseq = []
        for b in breakpoints:
            if p[b] != p[b + 1]:
                this_chimseq.append(1)
            else:
                this_chimseq.append(0)
        chimseq.append(this_chimseq)

    if ofile:
        nSeq = len(chimseq)
        nChar = len(chimseq[0])
        with open(ofile, 'w') as handle:
            handle.write(str(nSeq) + ' ' + str(nChar) + '\n')
            for i, s in enumerate(chimseq[:len(chimseq)]):
                handle.write('chimera_' + str("%02d" % i) + '\t' + ''.join([str(x) for x in s]) + '\n')

    return (chimseq)
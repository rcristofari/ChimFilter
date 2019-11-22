# Log-likelihood function for the probability of a chimera arising p:
import numpy as np

def lnlp(p, d, q):
    freqs = [x / sum(d) for x in d]
    lnls = [(np.log(p) - np.log(di))**(1 - q[i]) * (np.log(di - p) - np.log(di))**q[i] for i, di in enumerate(freqs)]
    lnl = sum(lnls)
    print(lnls)
    print(lnl)
    return(lnl)

import galois as gl
import numpy as np
import math
from typing import Iterable,SupportsIndex

syms = {0:'a',1:'b',2:'c',3:'d'}


 # Is the LUT supposed to keep track of probabilities or symbol combinations that give each symbol?
def generate_LUT(initial_probs: Iterable[SupportsIndex]):
    num_symnodes,q = np.shape(initial_probs)
    GF = gl.GF(q)
    LUT = []
    for sym in range(num_symnodes): # Symbol node to calculate probabilities for
        probs = [0 for _ in range(q)]
        other_vectors = [(x+sym+1)%num_symnodes for x in list(range(num_symnodes-1))] # Get input vectors
        for sym_prob in range(q): # Symbol to calculate probability for
            for element1 in range(q):
                for element2 in range(q):
                    # Check if GF elements equal the symbol of interest
                    if GF(element1) + GF(element2) == GF(sym_prob):
                        probs[sym_prob] += initial_probs[other_vectors[0]][element1] * initial_probs[other_vectors[1]][element2]
        LUT.append(probs)
    return LUT

class Symnode:
    def __init__(self,initial_probs):
        self.to_check = initial_probs
        self.dv,self.q = np.shape(initial_probs)
        self.GF = gl.GF(self.q)

    def update(self,from_check):
        # dv,q = np.shape(from_check) # degree of symbol node and size of GF
        # to_check = np.matrix([[0 for _ in range(q)] for _ in range(dv)]) # initialize empty matrix
        for i in range(self.dv): # select output array
            input_idxs = [(x+i+1)%self.dv for x in list(range(self.dv-1))] # Get input vectors
            for sym_out in range(self.q): # select symbol probability to calculate
                # Select input symbol probabilities
                print(f"Calculating probability for symbol: {[syms[sym_out]]}")
                p_sym_out = from_check[input_idxs[0]][sym_out]*from_check[input_idxs[1]][sym_out]
                self.to_check[i][sym_out] = p_sym_out
        return self.to_check

class Checknode:
    def __init__(self,LUT):
        self.LUT = LUT
        self.dc,self.q = np.shape(LUT)
        self.from_check = np.matrix([[0 for _ in range(q)] for _ in range(dc)])

    def update(self,to_check, LUT):
        # dc,q = np.shape(to_check) # degree of check node and size of GF
         # initialize empty matrix
        for i in range(self.dc): # select output array
            for sym in range(self.q): # select symbol probability to calculate

                # I am guessing that the LUT should not be calculated with symbol probabilities
                # otherwise the input from the symbol nodes is never used and it is just a matrix copy
                self.from_check[i][sym] = self.LUT[i][sym]

        return self.from_check


LUT = generate_LUT([[.2,.4,.3,.1],[.5,.1,.1,.3],[.2,.2,.4,.2]])
# r1 = checknode([[.2,.4,.3,.1],[.5,.1,.1,.3],[.2,.2,.4,.2]], LUT)
# print(r1)
# r2 = checknode([[.2,.7,0,.1],[.6,0,.1,.3],[1,0,0,0]])
# r3 = checknode([[.2,.7,0,.1],[.6,0,.1,.3],[0,0,1,0]])
# s = symnode([r1[0],r2[0],r3[0]])
# print(s)
# print([sum(x) for x in r])




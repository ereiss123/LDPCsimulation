import galois as gl
import numpy as np
import math

syms = {0:'a',1:'b',2:'c',3:'d'}

def symnode(from_check):
    dv = len(from_check) # degree of symbol node
    q = len(from_check[0]) # Size of GF
    to_check = []
    for i in range(dv): # select output array
        input_idxs = [(x+i+1)%dv for x in list(range(dv-1))]
        print(i, input_idxs)
        sym_vector = []
        for sym_out in range(q): # select symbol probability to calculate
            # Select input symbol probabilities
            print(f"Calculating probability for symbol: {[syms[sym_out]]}")
            p_sym_out = from_check[input_idxs[0]][sym_out]*from_check[input_idxs[1]][sym_out]
            sym_vector.append(p_sym_out)
        to_check.append(sym_vector.index(max(sym_vector)))
    return to_check

def checknode(to_check):
    dc = len(to_check) # degree of check node
    q = len(to_check[0]) # Size of GF
    Gf = gl.GF(q)
    from_check = []
    for i in range(dc): # select output array
        input_idxs = [(x+i+1)%dc for x in list(range(dc-1))]
        print(i, input_idxs)
        sym_vector = []
        for sym_out in range(q): # select symbol probability to calculate
            # Select input symbol probabilities
            p_sym_out = 0
            print(f"Calculating probability for symbol: {[syms[sym_out]]}")
            for sym1 in range(q):
                for sym2 in range(q):
                    if (Gf(sym1) + Gf(sym2)) == Gf(sym_out):
                        print(f"P({syms[sym1]})+P({syms[sym2]})")
                        print(f"{to_check[input_idxs[0]][sym1]}*{to_check[input_idxs[1]][sym2]}")
                        p_sym_out += to_check[input_idxs[0]][sym1]*to_check[input_idxs[1]][sym2]
            sym_vector.append(p_sym_out)
        from_check.append(sym_vector)
    return from_check


r1 = checknode([[.2,.4,.3,.1],[.5,.1,.1,.3],[.2,.2,.4,.2]])
r2 = checknode([[.2,.7,0,.1],[.6,0,.1,.3],[1,0,0,0]])
r3 = checknode([[.2,.7,0,.1],[.6,0,.1,.3],[0,0,1,0]])
s = symnode([r1[0],r2[0],r3[0]])
print(s)
# print([sum(x) for x in r])




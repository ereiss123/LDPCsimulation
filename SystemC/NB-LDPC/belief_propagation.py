import galois as gl
import numpy as np
import math

"""
This file contains an implemntation of a generalized belief propagation algorithm for NB-LDPC decoding. 
The assumption is made that the all-zero codeword is transmitted.    
"""

def build_H(N:int,K:int,weight:int, q:int) -> gl.FieldArray:
    """
    Builds a random NB-LDPC parity check matrix.
    
    Args:
        N (int): transmitted block length
        K (int): source block length
        weight (int): mean weight of each column, must be greater than 2
        q (int): order of the Galois field
    """
    assert N > K
    assert weight > 2
    M = N - K # number of parity check equations
    GF = gl.GF(q)
    H = GF(np.zeros((M,N),dtype=int)) # inirialize parity check matrix
    for i in range(N):
        for _ in range(weight):
            idx = np.random.randint(0,M) # select a random row
            while H[idx,i] != GF(0):
                idx = np.random.randint(0,M)
            H[idx,i] = GF(np.random.randint(1,q-1)) # select a random non-zero element
    return H 

def simulate_channel(N: int, sigma: float, mu:int, q: int) -> gl.FieldArray:
    """
    Generates a random noise sample vector.
    
    Args:
        N (int): transmitted block length
        sigma (float): standard deviation of the Gaussian distribution
        mu (int): mean of the Gaussian distribution
        q (int): order of the Galois field
    """
    GF = gl.GF(q)
    noise = np.sign(np.random.normal(0,sigma,int(N*math.log2(q)))) # Generate noise bits and make decision
    print(noise)
    noise = [int(x) for x in np.where(noise == -1, 0, noise)] # Convert to binary
    samples = GF(np.zeros(N,dtype=int))
    idx = 0
    for i in range(0,int(N*math.log2(q)),int(math.log2(q))):
        sample = noise[i:i+int(math.log2(q))]
        sample = int(''.join([str(x) for x in sample]),2)
        samples[idx] = GF(sample)
        idx += 1
    return samples

def belief_propagation(H: gl.FieldArray, samples: gl.FieldArray, max_iter: int, q: int) -> gl.FieldArray:
    """
    Performs belief propagation decoding on a NB-LDPC code assuming the all-zero codeword is transmitted.
    
    Args:
        H (galois.FieldArray): parity check matrix
        samples (galois.FieldArray): received samples
        max_iter (int): maximum number of iterations
        q (int): order of the Galois field
    """
    M,N = H.shape
    GF = gl.GF(q)
    
    # Initialize messages
    q = GF(np.zeros((M,N),dtype=int))
    r = GF(np.zeros((M,N),dtype=int))
    
    # Initialize q
    for m in range(M):
        for n in range(N):
            if H[m,n] != GF(0):
                sample = int(samples[n]) # pertinent sample
                sample_bin = bin(sample)[2:].zfill(int(math.log2(q))) # convert to binary
                
    return decision
    
    
    
        


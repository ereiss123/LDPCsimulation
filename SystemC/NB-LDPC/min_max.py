"""
This file contains an implemntation of a generalized min/max decoding algorithm for NB-LDPC decoding.
"""
import numpy as np
import galois as gl
import math

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
        
    Returns:
        gl.FieldArray: noise sample vector
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

def likelihood(y: float, sigma: float, bit: int) -> float:
    """_summary_
    
    Calculates the likelihood of a bit being 1 or 0 given the received sample. Assumes s = 1

    Args:
        y (float): received sample
        sigma (float): standard deviation of the Gaussian distribution
        bit (int): bit to calculate the likelihood for

    Returns:
        float: likelihood of the bit being 1 or 0
    """
    assert bit in range(2)
    g = 1/(1+math.exp(2*abs(y)/sigma**2 ))
    return g if bit == 1 else 1-g    

def decode():
    """
    
    """
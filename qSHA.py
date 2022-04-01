import cirq
from sage.all import *
from itertools import combinations, permutations
from sbox_gates import *
import numpy as np

initial_reg = [
            cirq.NamedQubit("initial_reg" + str(i)) for i in range(1600)
        ]
iota_reg = [
            cirq.NamedQubit("iota_reg" + str(i)) for i in range(1600)
        ]
chi_reg = [
            cirq.NamedQubit("chi_reg" + str(i)) for i in range(1600)
        ]
pi_reg = [
            cirq.NamedQubit("pi_reg" + str(i)) for i in range(1600)
        ]
rho_reg = [
            cirq.NamedQubit("rho_reg" + str(i)) for i in range(1600)
       ]
theta_reg = [ cirq.NamedQubit("theta_reg" + str(i)) for i in range(1600)]

circuit = cirq.Circuit()

def main():
    print("created vars...")
    msg = [0]*1600 # CLASSICAL CONSTRUCTION
    SHA3(msg)
    print(f"First count: {count_gates(circuit)}")


def configure_circuit_input(configuration, orig_circuit=None, input_bit_list = [0]*1600):
    moment_ops = []
    for k in range(len(input_bit_list)):  # initialize Ubox to given ival
        if input_bit_list[k]:
            # X | U[k]
            # circuit.append()
            moment_ops.append(cirq.ops.X.on(configuration.U[k]))
    # TODO: Doubts
    circuit = cirq.Moment(moment_ops) + orig_circuit
    return circuit

import numpy as np
import random
l = 6  # value of l = {0, 1, 2, 3, 4, 5, 6}
b = 25*(2**l)  # b = state size (value of b = {25, 50, 100, 200, 400, 800, 1600} )
# For SHA-3 the value of ‘l’ was 6 and the
# so rounds turn out to be
rounds = 12 + 2*l  # 24 rounds
print(str(rounds)+' Rounds in SHA-3')
# So SHA-3 has state size of 1600 bits and the number of rounds of computations will be 24

def quantum_parity_addition(summands, target):
    for i in summands:
        circuit.append(cirq.CNOT(i,target))

# 1600 bits(1 dimensional array) to 3 dimensional array of 5x5x64
def _1Dto3D(A):
    A = np.array(A).reshape(5,5,64)
    A_out = np.array(initial_reg).reshape(5,5,64) # Initialize empty 5x5x64 array
    for i in range(5):
        for j in range(5):
            for k in range(64):
                if A[i][j][k]:
                    circuit.append(cirq.X(A_out[i][j][k]))
    return A_out


def theta(A):
        A_out = np.array(theta_reg).reshape(5,5,64)  # Initialize empty 5x5x64 array
       #A_out = [[[0 for _ in range(64)] for _ in range(5)] for _ in range(5)] #without numpy
        for i in range(5):
                for j in range(5):
                        for k in range(64):
                            target = A_out[i][j][k]
                            C=[A[(i-1)%5][ji][k] for ji in range(5)] # 5 bit column "to the left" of the original bit
                            D=[A[((i+1) % 5)][ji][(k-1)%64] for ji in range(5)] #5 bit column "to the right"  and one position "to the front" of the original bit
                            quantum_parity_addition(C, target)
                            quantum_parity_addition(D, target)
                            circuit.append(cirq.CNOT(A[i][j][k],A_out[i][j][k]))

        return A_out

#Rho : Each word is rotated by a fixed number of position according to table.
def rho(A):
    rhomatrix=[[0,36,3,41,18],[1,44,10,45,2],[62,6,43,15,61],[28,55,25,21,56],[27,20,39,8,14]]
    rhom = np.array(rhomatrix, dtype=int)  # Initialize empty 5x5x64 array
    A_out = np.array(rho_reg).reshape(5,5,64)
    for i in range(5):
        for j in range(5):
            for k in range(64):
                A_out[i][j][k] = A[i][j][k - rhom[i][j]] #  A[i][j][k − (t + 1)(t + 2)/2] so here rhom[i][j] Use lookup table to "calculate" (t + 1)(t + 2)/2
    return A_out

#Pi: Permutate the 64 bit words
def pi(A):
    A_out = np.array(pi_reg).reshape(5,5,64) # Initialize empty 5x5x64 array
    for i in range(5):
        for j in range(5):
            for k in range(64):
                A_out[j][(2*i+3*j)%5][k] = A[i][j][k]
    return A_out

# A_out [i][j][k] = A[i][j][k] XOR ( (A[i + 1][j][k] XOR 1) AND (ain[i + 2][j][k]) )
def chi(A):
    A_out = np.array(chi_reg).reshape(5,5,64) # Initialize empty 5x5x64 array
    for i in range(5):
        for j in range(5):
            for k in range(64):
                circuit.append(cirq.X(A[(i + 1)%5][j][k]))
                circuit.append(cirq.TOFFOLI(A[(i + 1)%5][j][k],A[(i + 2)%5][j][k], A_out[i][j][k]))
                circuit.append(cirq.X(A[(i + 1)%5][j][k]))
                circuit.append(cirq.CNOT(A[i][j][k], A_out[i][j][k]))
    return A_out

#iota: add constants  to word (0,0)
# aout[i][j][k] = ain[i][j][k] ⊕ bit[i][j][k]
# for 0 ≤ ℓ ≤ 6, we have bit[0][0][2ℓ − 1] = rc[ℓ + 7ir]
def iota(A, round):
    # Initialize empty arrays
    A_out = A.copy()
    bit = np.zeros((5,5,64), dtype=int)
    rc = np.zeros((168), dtype=int)

    #generation of rc as Linear Feedback Shift Register
    w = np.array([1,0,0,0,0,0,0,0], dtype = int)
    rc[0] = w[0]
    for i in range(1, 168): #7*24
        w = [w[1],w[2],w[3],w[4],w[5],w[6],w[7], (w[0]+w[4]+w[5]+w[6]) % 2]
        rc[i] = w[0]

    # Calculate A_out
    for l in range(7):
        if rc[l + 7*round]:
            circuit.append(cirq.X(A_out[0][0][2**l - 1] ))

# 5x5x64 (three-dimensional array) into 1600 bits(one-dimensional array)
def _3Dto1D(A):
    A_out = np.zeros(1600, dtype = int) # Initialize empty array of size 1600
    for i in range(5):
        for j in range(5):
            for k in range(64):
                A_out[64*(5*j+i)+k] = A[i][j][k]
    return A_out

# 24 X (ι ◦ χ ◦ π ◦ ρ ◦ θ)
def SHA3(SHA_in):
    length=len(SHA_in)
    A_3D = _1Dto3D(SHA_in)
    for r in range(24):
        SHA_out_3D = iota(chi(pi(rho(theta(A_3D)))), r)
    print("complete!")
    return 0

main()
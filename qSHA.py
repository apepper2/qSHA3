import cirq
from sage.all import *
from itertools import combinations, permutations
from sbox_gates import *
import numpy as np
from ClassicalSimulator import ClassicalSimulator

import classical


circuit = cirq.Circuit()

def main():
    simulator = ClassicalSimulator()
    print("simulator is on ...")
    print("created vars...")
    msg = [1] + [0]*1598 + [1] # CLASSICAL CONSTRUCTION
    print("Rendering Circuit...")
    output_register = SHA3(msg)
    print(f"First count: {count_gates(circuit)}")
    print("simulating circuit...")
    result = simulator.run(circuit)
    print("simulation complete!")
    quantum_values = read_meas_results(output_register, result)
    print(quantum_values)


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
# For SHA-3 the value of ālā was 6 and the
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
    initial_reg = [cirq.NamedQubit("initial_reg" + str(i)) for i in range(1600)]
    A_out = np.array(initial_reg).reshape(5,5,64) # Initialize empty 5x5x64 array
    for i in range(5):
        for j in range(5):
            for k in range(64):
                if A[i][j][k]:
                    circuit.append(cirq.X(A_out[i][j][k]))
    return A_out


def theta(A):
        theta_reg = [ cirq.NamedQubit("theta_reg" + str(i)) for i in range(1600)]
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
    A_out = A.copy()
    for i in range(5):
        for j in range(5):
            for k in range(64):
                A_out[i][j][k] = A[i][j][k - rhom[i][j]] #  A[i][j][k ā (t + 1)(t + 2)/2] so here rhom[i][j] Use lookup table to "calculate" (t + 1)(t + 2)/2
    return A_out

#Pi: Permutate the 64 bit words
def pi(A):
    A_out = A.copy() # Initialize empty 5x5x64 array
    for i in range(5):
        for j in range(5):
            for k in range(64):
                A_out[i,j,k] = A[(i + 3*j)%5,i,k]
    return A_out

# A_out [i][j][k] = A[i][j][k] XOR ( (A[i + 1][j][k] XOR 1) AND (ain[i + 2][j][k]) )
def chi(A):
    for j in range(5):
        for k in range(64):
            temp_reg = [cirq.NamedQubit("temp_reg" + str(i)) for i in range(2)]
            temp_reg[0] = A[0,j,k]
            temp_reg[1] = A[1,j,k]
            for i in range(3):
                circuit.append(cirq.X(A[(i + 1),j,k]))
                circuit.append(cirq.TOFFOLI(A[(i + 1),j,k],A[(i + 2),j,k], A[i,j,k]))
                circuit.append(cirq.X(A[(i + 1),j,k]))

            circuit.append(cirq.X(A[4,j,k]))
            circuit.append(cirq.TOFFOLI(A[4,j,k], temp_reg[0], A[3,j,k]))
            circuit.append(cirq.X(A[4,j,k]))
            circuit.append(cirq.X(temp_reg[0]))
            circuit.append(cirq.TOFFOLI(temp_reg[1], temp_reg[0], A[4,j,k]))
            circuit.append(cirq.X(temp_reg[0]))
    return A

#iota: add constants  to word (0,0)
# aout[i][j][k] = ain[i][j][k] ā bit[i][j][k]
# for 0 ā¤ ā ā¤ 6, we have bit[0][0][2ā ā 1] = rc[ā + 7ir]

def rc(t):
    top = t % 255
    if top == 0:
        return 1
    else:
        R= [1,0,0,0,0,0,0,0]
        for i in range(1,top+1):
            R = [0]+R
            R[0] = R[0] + R[8] % 2
            R[4] = R[4] + R[8] % 2
            R[5] = R[5] + R[8] % 2
            R[6] = R[6] + R[8] % 2
            R = R[:8]
        return R[0]

def iota(A, round):
    # Initialize empty arrays
    RC = np.zeros((64), dtype=int)
    for l in range(7):
        RC[2**l - 1] = rc(l+7*round)
    for z in range(len(RC)):
        if RC[z]:
            circuit.append(cirq.X(A[0][0][z]))
    return A

# 5x5x64 (three-dimensional array) into 1600 bits(one-dimensional array)
def _3Dto1D(A):
    A_out = [0]*1600# Initialize empty array of size 1600
    for i in range(5):
        for j in range(5):
            for k in range(64):
                circuit.append(cirq.ops.measure(A[i][j][k]))
                A_out[64*(5*j+i)+k] = A[i][j][k]   
    return A_out

def read_meas_results(register, result):
    quantum_u_values = []
    for qubit in register:
        quantum_u_values += [int(result.measurements[qubit.name])]
    return quantum_u_values


# 24 X (Ī¹ ā¦ Ļ ā¦ Ļ ā¦ Ļ ā¦ Īø)
def SHA3(SHA_in):
    length=len(SHA_in)
    A_3D = _1Dto3D(SHA_in)
    for r in range(24):
        SHA_out_3D = iota(chi(pi(rho(theta(A_3D)))), r)
    sha_out_1D = _3Dto1D(SHA_out_3D)
    print("complete!")
    return(sha_out_1D)

main()
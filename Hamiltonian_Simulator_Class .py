import matplotlib.pyplot as plt
import numpy as np

# qutip
from qutip import *

# Numerical
from scipy.linalg import expm
from scipy.sparse.linalg import expm as sparse_expm

# Trotter
from numpy.linalg import matrix_power

# Time performance
import time

# 
from scipy import sparse


class Hamiltonian_Simulator:

    def __init__(self):
        pass
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    ## Methods---------------------------------------------------------
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    def benchmark(self, qubits_nums, initial_states, Js, h, trotter_nums, times, options):
        
        print("Starting Numerical Simulation")
        numerical_times = []
        numerical_states_at_diff_qubits = []
        
        # Time performance
        begin_time = time.time()

        for i in range(len(Js)):
            states = numerical_sim(initial_states[i], Js[i], h, times)
            numerical_times.append(states[1])
            numerical_states_at_diff_qubits.append(states[0])
            print(str(i + 1) + " qubits: done")
            
        numerical_entire_sim_time = time.time() - begin_time
        print("Numerical simulation done")
        print("Total time to complete Numerical simulation is " + str(numerical_entire_sim_time) + ".")
        
        # Time performance
        begin_time = time.time()
        
        trotter_times = []
        trotter_states_at_diff_qubits_at_diff_trotter_num = []# a list of a list(states for n qubits at diff trotter nums) 
                                                              # of a list(states for spe)
        for i in range(len(Js)):
            t = []
            trotter_states_for_diff_trotter_num = []
            for j in range(len(trotter_nums)):
                states = trotter_sim(initial_states[i], Js[i], h, trotter_nums[j], times)
                trotter_states_for_diff_trotter_num.append(states[0])
                t.append(states[1])
            trotter_states_at_diff_qubits_at_diff_trotter_num.append(trotter_states_for_diff_trotter_num)
            trotter_times.append(t)
            print(str(i + 1) + " qubits: done")
        trotter_entire_sim_time = time.time() - begin_time
        print("Qutip simulation done")
        print("Total time to complete qutip simulation is " + str(trotter_entire_sim_time) + ".")
        
        
        
    # perfrom qutip simulation for an n-qubits system where n = len(J). 
    # initial_state must be a ket qutip Qobj representing the initial_state of the system with shape(2^n, 1).
    # J must be an nxn matrix and h must be a real number.
    # times is a list of time points for which it is desired to know the state of the system.
    # options is of type Options() to control the qutip integrator.
    # returns a list of Qobj of type ket representing the state of the system at every time point in times.
    # returns the time needed to simulate
    def qutip_simulate(self, initial_state, J, h, times, options):
        
        # Time performance
        begin_time = time.time()
        
        N = len(J)
        H = self.construct_hamiltonian(J, h)
        result = sesolve(H, initial_state, times, options = options)
        return result.states, time.time() - begin_time


    # Construct the Trapped Ion Hamiltonian (Top section) for an n-qubits system where n = len(J) for a given J and h
    # J must be an nxn matrix and h must be a real number
    # returns a qutip Qobj representing the Trapped Ion Hamiltonian for the n-qubits system
    def construct_hamiltonian(self, J, h):
        N = len(J)
        
        sum_sx = []
        # build h times sum sigma_j^x
        for i in range(N):
            # make a list of N qeye matrices for one qubit
            op_list = []
            for _ in range(N):
                op_list.append(qeye(2))  
            # replace the jth element with sigmax, tensor, then append sum_sx
            op_list[i] = sigmax()
            sum_sx.append(h*tensor(op_list))
            
        sum_sz_sz = []
        # bild sum J_{i,j} sigma_i^z sigma_j^z, we ignore the case where i = j
        for i in range(N):      
            for j in range(N):
                if i == j:
                    continue
                # make a list of N qeye matrices for one qubit
                op_list = []
                for _ in range(N):
                    op_list.append(qeye(2))  
                # replace the ith and jth elements with sigma_z, tensor, then append to sum_sz_sz
                op_list[i] = sigmaz()
                op_list[j] = sigmaz()
                sum_sz_sz.append(J[i][j]*tensor(op_list))
        H = 0
        for op in sum_sz_sz:
            H -= op
        for op in sum_sx:
            H -=op
        return H
    
    
    # perfrom Trotter formula simulation for an n-qubits system where n = len(J). 
    # initial_state must be a ket qutip Qobj representing the initial_state of the system with shape(2^n, 1).
    # J must be an nxn matrix and h must be a real number.
    # Trotter_num is a integer
    # times is a list of time points for which it is desired to know the state of the system.
    # returns a list representing the state of the system at every time point in times.
    # returns the time needed to simulate

    def trotter_sim(self, initial_state, J, h, trotter_num, times):
        # Time performance
        begin_time = time.time()
        
        initial_state = np.asarray(initial_state)
        
        N = len(J)

        sum_sx = []
        # build h times sum sigma_j^x
        for i in range(N):
            # make a list of N qeye matrices for one qubit
            op_list = []
            for _ in range(N):
                op_list.append(qeye(2))  
            # replace the jth element with sigmax, tensor, then append sum_sx
            op_list[i] = sigmax()
            sum_sx.append(h*tensor(op_list))
            
        B = 0
        for op in sum_sx:
            B -= op
        B = np.asarray(B)
        
        sum_sz_sz = []
        # bild sum J_{i,j} sigma_i^z sigma_j^z, we ignore the case where i = j
        for i in range(N):      
            for j in range(N):
                if i == j:
                    continue
                # make a list of N qeye matrices for one qubit
                op_list = []
                for _ in range(N):
                    op_list.append(qeye(2))  
                # replace the ith and jth elements with sigma_z, tensor, then append to sum_sz_sz
                op_list[i] = sigmaz()
                op_list[j] = sigmaz()
                sum_sz_sz.append(J[i][j]*tensor(op_list))
                
        A = 0        
        for op in sum_sz_sz:
            A -= op
        A = np.asarray(A)
        
        # IMPORTANT: if we are simulating a 1-qubit system, there won't be a sum_sz_sz term
        # so we don't have to perform the matrix multiplication of np.matmul(expm(-1j*A*1/trotter_num), expm(-1j*B*1/trotter_num))
        states = []
        if N == 1:
            for t in times:
                U = expm(-1j*B*t/trotter_num)
                U = matrix_power(U, trotter_num)
                states.append(np.matmul(U, initial_state))
        else:
            for t in times:
                U = np.matmul(expm(-1j*A*t/trotter_num), expm(-1j*B*t/trotter_num))
                U = matrix_power(U, trotter_num)
                states.append(np.matmul(U, initial_state))            
        return states, time.time() - begin_time

    # This algorithm perform the same things as trotter_sim but with sparse matrcies rather than normal matricies.
    def sparse_trotter_sim(self, initial_state, J, h, trotter_num, times):
        # Time performance
        begin_time = time.time()
        
        initial_state = np.asarray(initial_state)
        
        N = len(J)

        sum_sx = []
        # build h times sum sigma_j^x
        for i in range(N):
            # make a list of N qeye matrices for one qubit
            op_list = []
            for _ in range(N):
                op_list.append(qeye(2))  
            # replace the jth element with sigmax, tensor, then append sum_sx
            op_list[i] = sigmax()
            sum_sx.append(h*tensor(op_list))
            
        B = 0
        for op in sum_sx:
            B -= op
        B = np.asarray(B)
        B = sparse.csr_matrix(B)
        
        sum_sz_sz = []
        # bild sum J_{i,j} sigma_i^z sigma_j^z, we ignore the case where i = j
        for i in range(N):      
            for j in range(N):
                if i == j:
                    continue
                # make a list of N qeye matrices for one qubit
                op_list = []
                for _ in range(N):
                    op_list.append(qeye(2))  
                # replace the ith and jth elements with sigma_z, tensor, then append to sum_sz_sz
                op_list[i] = sigmaz()
                op_list[j] = sigmaz()
                sum_sz_sz.append(J[i][j]*tensor(op_list))
                
        A = 0        
        for op in sum_sz_sz:
            A -= op
        A = np.asarray(A)
        A = sparse.csr_matrix(A)
        
        # IMPORTANT: if we are simulating a 1-qubit system, there won't be a sum_sz_sz term
        # so we don't have to perform the matrix multiplication of np.matmul(expm(-1j*A*1/trotter_num), expm(-1j*B*1/trotter_num))
        states = []
        if N == 1:
            for t in times:
                U = sparse_expm(-1j*B*t/trotter_num)
                U = U**trotter_num
                states.append(U @ initial_state)
        else:
            for t in times:
                U = sparse_expm(-1j*A*t/trotter_num) @ sparse_expm(-1j*B*t/trotter_num)
                U = U**trotter_num
                states.append(U @ initial_state)        
        return states, time.time() - begin_time
    
    
    
    def numerical_sim(self, initial_state, J, h, times):
        # Time performance
        begin_time = time.time()
        
        H_numerical = np.asarray(self.construct_hamiltonian(J, h))

        initial_state_numerical = np.asarray(initial_state)

        numerical_states = []
        for t in times:
            U = expm(-1j*H_numerical*t)
            mult = np.matmul(U, initial_state_numerical)
            numerical_states.append(mult)
        return numerical_states, time.time() - begin_time
h = 1
J = np.array([[1, 1], [1, 1]])
initial_state = tensor(basis(2, 0), basis(2, 0))
trotter_num = 100
times = np.linspace(0.0, 10, 100)
s = Hamiltonian_Simulator()
r = s.numerical_sim(initial_state, J, h, times)
print(r)
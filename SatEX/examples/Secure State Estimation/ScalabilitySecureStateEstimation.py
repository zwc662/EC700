#  * --------------------------------------------------------------------------
#  * File: ScalabilityBoolean.py
#  * ---------------------------------------------------------------------------
#  * Copyright (c) 2016 The Regents of the University of California.
#  * All rights reserved.
#  *
#  * Redistribution and use in source and binary forms, with or without
#  * modification, are permitted provided that the following conditions
#  * are met:
#  * 1. Redistributions of source code must retain the above copyright
#  *    notice, this list of conditions and the following disclaimer.
#  * 2. Redistributions in binary form must reproduce the above
#  *    copyright notice, this list of conditions and the following
#  *    disclaimer in the documentation and/or other materials provided
#  *    with the distribution.
#  * 3. All advertising materials mentioning features or use of this
#  *    software must display the following acknowledgement:
#  *       This product includes software developed by Cyber-Physical &
#  *       Systems Lab at UCLA.
#  * 4. Neither the name of the University nor that of the Laboratory
#  *    may be used to endorse or promote products derived from this
#  *    software without specific prior written permission.
#  *
#  * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS''
#  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#  * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#  * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS
#  * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#  * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
#  * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#  * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
#  * SUCH DAMAGE.
#  *
#  * Developed by: Yasser Shoukry
#  */

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(dir_path+'/../../solver/')
from SMConvexSolver import *

sys.path.append(dir_path+'/../../solver/z3/z3-4.4.1-x64-osx-10.11/bin/')
import z3 as z3



import scipy as sci
from scipy.io import loadmat
from scipy.linalg import expm

import timeit
import math

import matplotlib.pyplot as plt

import scipy as sci
from scipy.io import loadmat
from scipy.linalg import expm

import sys

import numpy as np

import cplex as cplex
from cplex.exceptions import CplexError

import random as rand

import signal


TIMEOUTE_LIMIT              = 600

# Register an handler for the timeout
#def timeout_handler(signum, frame):
#    print "------> Timeout <------"
#    raise Exception("end of time")


# Register an handler for the keyboard inteeruption
def keyboard_handler(signal, frame):
    print 'You pressed Ctrl+C! --- Keep pressing more until you exit'
    sys.exit(0)

#***************************************************************************************************
#***************************************************************************************************
#
#         TEST1: INCREASE NUMBER OF BOOL CONSTRAINTS
#
#***************************************************************************************************
#***************************************************************************************************

def scalability_main():
    signal.signal(signal.SIGINT, keyboard_handler)
    
    max_sensors_under_attack                        = 5
    attack_power                                    = 1.0
    n                                               = 5
    number_of_sensors                               = [25, 50, 75, 100, 125, 150, 175, 200]



    # ------------------------------------------------------------------------
    #       INITIALIZE DATASTRUCTURES
    # ------------------------------------------------------------------------
    time_SMC                = list()
    time_SMC_IIS            = list()
    time_MILP_1core         = list()
    time_MILP_2core         = list()
    time_MILP_3core         = list()
    time_MILP_4core         = list()
    time_SMT                = list()
    
    
    for p in number_of_sensors:
        print '\n********************************************************************'
        print '************* New SSE Test Case: #Sensors = ', p,' *****************'
        print '********************************************************************'
        
        # ------------------------------------------------------------------------
        #       GENERATE DYNAMICAL SYSTEM
        # ------------------------------------------------------------------------
        n,p,A,B,C,Y,x_init,attacked_sensor_index,safe_sensors = SSE_Random(n,p,max_sensors_under_attack,attack_power)
        attacked_sensor_index                       = sorted(attacked_sensor_index)
        
        print 'attacked_sensor_index', attacked_sensor_index, 'state',x_init
        
        
        # ------------------------------------------------------------------------
        #       GENERATE CONSTRAINTS
        # ------------------------------------------------------------------------
        constraints             = list()
        numberOfBoolVars        = p
        numberOfRealVars        = n
        numberOfConvexConstraints = 0
        for counter in range(0, p):
            O_i                 = obsv(A, C[counter,:])
            Y_i                 = Y[:,counter]
            
            Q_i                 = np.dot(np.transpose(O_i),O_i)
            c_i                 = -2*np.dot(np.transpose(Y_i),O_i)
            c_i                 = [c_i[0,count] for count in range(0,n)]
            
            b_i                 = -1*np.dot(np.transpose(Y_i),Y_i)


            convexConstraint    = {'Q': Q_i, 'c':c_i, 'b':b_i, 'sense':'L'}
            boolConstraint      = [counter]

            constraints.append({'bool':boolConstraint, 'convex':convexConstraint})
            numberOfConvexConstraints = numberOfConvexConstraints + 1

        # ------------------------------------------------------------------------
        #       LAUNCH SMC SOLVER - Sum Of Slacks
        # ------------------------------------------------------------------------
        
        time_SMC.append(scalability_SMC(numberOfBoolVars, numberOfRealVars, numberOfConvexConstraints, constraints, max_sensors_under_attack, 'SMC'))
        
        # Write SMC results
        print '--> Writing Execution Time Results to File'
        SMCfilename             = 'scalbality_SMC_results.txt'
        SMCfile                 = open(SMCfilename, 'w')
        for item in time_SMC:
            SMCfile.write("%s\n" % item)
        SMCfile.close()
        
        # ------------------------------------------------------------------------
        #       LAUNCH SMC SOLVER - IIS
        # ------------------------------------------------------------------------
        
        time_SMC_IIS.append(scalability_SMC(numberOfBoolVars, numberOfRealVars, numberOfConvexConstraints, constraints, max_sensors_under_attack, 'IIS'))
        
        # Write SMC results
        print '--> Writing Execution Time Results to File'
        SMCfilename             = 'scalbality_SSE_SMC_IIS_results.txt'
        SMCfile                 = open(SMCfilename, 'w')
        for item in time_SMC_IIS:
            if item > TIMEOUTE_LIMIT:
                SMCfile.write("TIME OUT\n")
            else:
                SMCfile.write("%s\n" % item)
        SMCfile.close()


        # ------------------------------------------------------------------------
        #       LAUNCH MILP SOLVER - 1 core
        # ------------------------------------------------------------------------
        time_MILP_1core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, max_sensors_under_attack, 1))

        # Write CPLEX results - 1 core
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_SSE_MILP1_results.txt'
        MILPfile                = open(MILPfilename, 'w')
        for item in time_MILP_1core:
            if item > TIMEOUTE_LIMIT:
                MILPfile.write("TIME OUT\n")
            else:
                MILPfile.write("%s\n" % item)
        MILPfile.close()


        # ------------------------------------------------------------------------
        #       LAUNCH MILP SOLVER - 2 core2
        # ------------------------------------------------------------------------
        time_MILP_2core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, max_sensors_under_attack, 2))

        # Write CPLEX results - 2 core2
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_SSE_MILP2_results.txt'
        MILPfile                = open(MILPfilename, 'w')
        for item in time_MILP_2core:
            if item > TIMEOUTE_LIMIT:
                MILPfile.write("TIME OUT\n")
            else:
                MILPfile.write("%s\n" % item)
        MILPfile.close()

        # ------------------------------------------------------------------------
        #       LAUNCH MILP SOLVER - 3 cores
        # ------------------------------------------------------------------------
        time_MILP_3core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, max_sensors_under_attack, 3))

        # Write CPLEX results - 3 cores
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_SSE_MILP3_results.txt'
        MILPfile                = open(MILPfilename, 'w')
        for item in time_MILP_3core:
            if item > TIMEOUTE_LIMIT:
                MILPfile.write("TIME OUT\n")
            else:
                MILPfile.write("%s\n" % item)
        MILPfile.close()


        # ------------------------------------------------------------------------
        #       LAUNCH MILP SOLVER - 4 cores
        # ------------------------------------------------------------------------
        time_MILP_4core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, max_sensors_under_attack, 4))

        # Write CPLEX results - 4 cores
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_SSE_MILP4_results.txt'
        MILPfile                = open(MILPfilename, 'w')
        for item in time_MILP_4core:
            if item > TIMEOUTE_LIMIT:
                MILPfile.write("TIME OUT\n")
            else:
                MILPfile.write("%s\n" % item)
        MILPfile.close()
        
        
        # ------------------------------------------------------------------------
        #       LAUNCH SMT SOLVER
        # ------------------------------------------------------------------------
        '''
        time_SMT.append(scalability_Z3(numberOfBoolVars, numberOfRealVars, constraints, max_sensors_under_attack))

        # Write SMT results
        print '--> Writing Execution Time Results to File'
        Z3filename             = 'scalbality_SSE_Z3_results.txt'
        Z3file                 = open(Z3filename, 'w')
        for item in time_SMT:
            if item > TIMEOUTE_LIMIT:
                Z3file.write("TIME OUT\n")
            else:
                Z3file.write("%s\n" % item)
        Z3file.close()
        '''
        
        # ------------------------------------------------------------------------
        #       PRINT RESULTS
        # ------------------------------------------------------------------------

        print '-----> Execution Time <------'
        print 'time_SMC',time_SMC
        print 'time_SMC_IIS',time_SMC_IIS
        print 'time_MILP_1core',time_MILP_1core
        print 'time_MILP_2core',time_MILP_2core
        print 'time_MILP_3core',time_MILP_3core
        print 'time_MILP_4core',time_MILP_4core
        print 'time_SMT',time_SMT
        


    # ------------------------------------------------------------------------
    #       PLOT THE FINAL RESULTS
    # ------------------------------------------------------------------------
    
    plt.semilogy(number_of_sensors, time_SMC, 'b',
                 number_of_sensors, time_SMC_IIS, 'b--',
                 number_of_sensors, time_MILP_1core, 'r',
                 number_of_sensors, time_MILP_2core, 'r--',
                 number_of_sensors, time_MILP_3core, 'r:',
                 number_of_sensors, time_MILP_4core, 'r.-',
                 number_of_sensors, time_SMT, 'g'
                 )
    plt.legend(['SatEX','CPLEX 1 core', 'CPLEX 2 cores', 'CPLEX 3 cores', 'CPLEX 4 cores', 'Z3'])
    plt.xlabel('number of sensors')
    plt.ylabel('Execution Time')
    plt.draw()
    plt.show()

#***************************************************************************************************
#***************************************************************************************************
#
#         SCALABILITY Test - MILP
#
#***************************************************************************************************
#***************************************************************************************************


def scalability_CPLEX(numberOfBoolVars, numOfRealVars, constraints, max_sensors_under_attack, numberOfCores):
    print '\n======================================================='
    print '      TEST CPLEX'
    print '            #Real Varabiles = ',     numOfRealVars
    print '            #Bool Varabiles = ',     numberOfBoolVars
    print '            #Bool Constraints = ',   len(constraints)
    print '            #Cores = ',   numberOfCores
    print '=======================================================\n'


    mipSolver = cplex.Cplex()
    rVars                               = ['x'+str(i) for i in range(0, numOfRealVars)]
    bVars                               = ['b'+str(i) for i in range(0, numberOfBoolVars)]
    mipSolver.variables.add(    names = bVars,
                                types = [mipSolver.variables.type.binary] * len(bVars)
                            )
    
    mipSolver.variables.add(names       =   rVars,
                            ub          =   [cplex.infinity] * len(rVars),
                            lb          =   [-cplex.infinity] * len(rVars)
                            )
    
    mipSolver.parameters.threads.set(numberOfCores)
    mipSolver.parameters.timelimit.set(TIMEOUTE_LIMIT)
    bigM                                = 1E6
    
    #----------------------------------------------------------------------------
    print '--> CPLEX: FEEDING CONSTRAINTS'
    #----------------------------------------------------------------------------
    convexIFCounter                     = 0
    for constraint in constraints:
        boolConstraint                  = constraint['bool']
        #vars                            = ['b'+str(abs(i)-1)    for i in boolConstraint]
        #signs                           = [np.sign(i)           for i in boolConstraint]
        #numOfNegativeSigns              = sum(x < 0             for x in signs)
        #mipSolver.linear_constraints.add(lin_expr   = [cplex.SparsePair(ind = vars,val =  signs)],
        #                                rhs         = [1 -numOfNegativeSigns],
        #                                senses      = ['G'])
        convexConstraint                = constraint['convex']

        clause                          = QPClause(convexConstraint['Q'], convexConstraint['c'], convexConstraint['b'], rVars)
        my_ind                          = ['x'+str(i) for i in range(0, numOfRealVars)] + ['b'+str(abs(boolConstraint[0]))]
        my_val                          = convexConstraint['c'] + [-1*bigM]
        
        mip_lin_expr                    = cplex.SparsePair(ind = my_ind, val = my_val)
                                                                                 
        mipSolver.quadratic_constraints.add(quad_expr = clause['quad_expr'], lin_expr = mip_lin_expr, rhs = clause['rhs'])
        convexIFCounter             = convexIFCounter + 1

    # add constraint on maximum number of sensors being under attack
    mipSolver.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = bVars, val =  [1.0]*len(bVars))],
                                     rhs      = [max_sensors_under_attack + 1],
                                     senses   = 'L')

    #----------------------------------------------------------------------------
    print '--> CPLEX: SEARCHING FOR A SOLUTION'
    #----------------------------------------------------------------------------
    start                               = timeit.default_timer()
    mipSolver.solve()
    end                                 = timeit.default_timer()
    time_smt_MILP                       = end - start

    print(mipSolver.solution.get_status(), mipSolver.solution.get_status_string())
    
    print 'CPLEX time = ', time_smt_MILP
    
    return time_smt_MILP






#***************************************************************************************************
#***************************************************************************************************
#
#         SCALABILITY Test - SMC
#
#***************************************************************************************************
#***************************************************************************************************


def scalability_SMC(numberOfBoolVars, numOfRealVars, numOfConvIFClauses, constraints, max_sensors_under_attack, strategy):
    print '\n======================================================='
    print '      TEST SMC'
    print '            #Real Varabiles = ',     numOfRealVars
    print '            #Bool Varabiles = ',     numberOfBoolVars
    print '            #Bool Constraints = ',   len(constraints)
    print '            UNSAT Certificate = ',   strategy
    print '=======================================================\n'
    
    solver                              = SMConvexSolver(numberOfBoolVars, numOfRealVars, numOfConvIFClauses,
                                                     maxNumberOfIterations = 10000,
                                                     counterExampleStrategy = strategy,
                                                     verbose = 'OFF',
                                                     profiling = 'false',
                                                     numberOfCores = 1,
                                                     slackTolerance = 1E-03)

    #----------------------------------------------------------------------------
    print '--> SMC: FEEDING CONSTRAINTS'
    #----------------------------------------------------------------------------
    convexIFCounter                     = 0
    
    for constraint in constraints:
        boolConstraint                  = constraint['bool']
        convexConstraint                = constraint['convex']
        
        firstVar                    = boolConstraint[0]
        b                           = convexConstraint['b']
        solver.setConvIFClause(QPClause(np.array(convexConstraint['Q']),
                                            np.array(convexConstraint['c']), b,
                                            solver.rVars,
                                            sense = convexConstraint['sense']),
                                    convexIFCounter)
        solver.addBoolConstraint(IMPLIES(NOT(solver.bVars[firstVar]), solver.convIFClauses[convexIFCounter]))
        convexIFCounter             = convexIFCounter + 1

    solver.addBoolConstraint(sum([BoolVar2Int(solver.bVars[i]) for i in range (0, numberOfBoolVars)]) == max_sensors_under_attack)

    #----------------------------------------------------------------------------
    print '--> SMC: SEARCHING FOR A SOLUTION'
    #----------------------------------------------------------------------------
    start                           = timeit.default_timer()
    x,b                             = solver.solve()
    end                             = timeit.default_timer()
    time_smt_SMC                    = end - start

    print('SMC time = ', time_smt_SMC)
    return time_smt_SMC



#***************************************************************************************************
#***************************************************************************************************
#
#         SCALABILITY Test - Z3
#
#***************************************************************************************************
#***************************************************************************************************

def scalability_Z3(numberOfBoolVars, numOfRealVars, constraints, max_sensors_under_attack):
    print '\n======================================================='
    print '      TEST Z3'
    print '            #Real Varabiles = ',     numOfRealVars
    print '            #Bool Varabiles = ',     numberOfBoolVars
    print '            #Bool Constraints = ',   len(constraints)
    print '=======================================================\n'

    SMTsolver                  = z3.Solver()
    bVars                      = z3.BoolVector('b', numberOfBoolVars)
    rVars                      = z3.RealVector('y', numOfRealVars)

    #----------------------------------------------------------------------------
    print '--> Z3: FEEDING CONSTRAINTS'
    #----------------------------------------------------------------------------
    for constraint in constraints:
        boolConstraint                  = constraint['bool']
        convexConstraint                = constraint['convex']

        firstVar                    = boolConstraint[0]
        Q                           = convexConstraint['Q']
        b                           = convexConstraint['b']
        c                           = convexConstraint['c']
        quadConstraint              = 0
        for counter1 in range(0,numOfRealVars):
            for counter2 in range(0,numOfRealVars):
                quadConstraint      = quadConstraint + rVars[counter1]*rVars[counter2]*Q[counter1,counter2]

        linearConstraint            = 0
        for counter1 in range(0,numOfRealVars):
            linearConstraint        = linearConstraint + rVars[counter1] * c[counter1]
        
        #quadConstraint                    = sum([i * j for i,j in zip(rVars,A)] + [-1*b[0]])
        SMTsolver.add(z3.Implies(bVars[firstVar], quadConstraint + linearConstraint < b))

    SMTsolver.add(sum([BoolVar2Int(bVars[i]) for i in range (0, numberOfBoolVars)]) == max_sensors_under_attack)
    #----------------------------------------------------------------------------
    print '--> Z3: SEARCHING FOR A SOLUTION'
    #----------------------------------------------------------------------------
    SMTsolver.set("timeout", TIMEOUTE_LIMIT*1000)
    start                               = timeit.default_timer()
    SMTsolver.check()
    end                                 = timeit.default_timer()
    time_smt_SMT                        = end - start

    print('Z3 time = ', time_smt_SMT)
    return time_smt_SMT



#***************************************************************************************************
#***************************************************************************************************
#
#         Generate Random Secure State Estimation Problem
#
#***************************************************************************************************
#***************************************************************************************************

def SSE_Random(n,p,s_bar,attackpower):
    print("Secure State Estimation - Random")
    np.random.seed(0)
    
    A                               = sci.sparse.rand(n, n, density=1.0)
    B                               = [0.0] * n
    C                               = 10*sci.sparse.rand(p, n, density=1.0)
    
    A                               = A.toarray()
    U, s, V                         = np.linalg.svd(A)
    max_eig                         = max(s)
    A                               = A/(max_eig + 0.1)
    
    C                               = C.toarray()
    
    max_sensors_under_attack        = s_bar;           # maximum sensors under attack
    
    
    per                             = np.random.permutation(p-1); # the attacker randmoizes the sensor indecies
    attacked_sensor_index           = per[0:max_sensors_under_attack]; # pick the first 14 sensors to attack them
    safe_sensors                    = []
    
    # simulate the system bus network
    x                               = np.random.rand(n,1);      # unknown initial condition
    x_init                          = x
    simulation_time                 = n;                      #simulation time (number of samples)
    Y                               = np.transpose(np.dot(C, x))
    for t in range(0,simulation_time-1):
        x                                       = np.dot(A, x)
        # Generate a random attack vector
        attack_signal                           = np.zeros((p,1))
        attack_signal[attacked_sensor_index]    = attackpower*np.random.rand(len(attacked_sensor_index),1);
        
        y                                       = np.dot(C, x) + attack_signal
        Y                                       = np.append(Y, np.transpose(y), axis = 0)


    return n,p,A,B,C,Y,x_init,attacked_sensor_index,safe_sensors


#***************************************************************************************************
#***************************************************************************************************
#
#         Helper Functions
#
#***************************************************************************************************
#***************************************************************************************************



def obsv(A, C):
    """Observability matrix
        Parameters
        ----------
        A, C: array_like or string
        Dynamics and output matrix of the system
        Returns
        -------
        O: matrix
        Observability matrix
        Examples
        --------
        >>> O = obsv(A, C)
        """
    
    # Convert input parameters to matrices (if they aren't already)
    amat = np.mat(A)
    cmat = np.mat(C)
    n = np.shape(amat)[0]
    
    # Construct the controllability matrix
    obsv = cmat
    for i in range(1, n):
        obsv = np.vstack((obsv, cmat*amat**i))
    return obsv


if __name__ == "__main__":
    np.random.seed(0)
    scalability_main()



     
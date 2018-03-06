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


TIMEOUTE_LIMIT              = 100

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
    
    
    numberOfRealVars_testCases  = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
    numberOfBooleanConstraints  = 7000

    print '===================== TEST: SCALABILITY_REAL ======================='
    
    
    # ------------------------------------------------------------------------
    #       READ THE BOOLEAN CONSTRIANTS
    # ------------------------------------------------------------------------
    infile = '003-23-80.cnf'
    with open(infile) as f:
        content = f.readlines()
    
    # ---- skip over all comments
    documentStart       = 0
    for line in content:
        if line[0] == 'c':
            documentStart = documentStart  +1
        else:
            break

    # ---- get number of constraints and number of variables
    probelmDefLine          = content[documentStart].split(" ")
    numberOfBoolVars        = int(probelmDefLine[2])
    maxNumberOfConstraints  = int(probelmDefLine[3])
    print 'maxNumberOfConstraints', maxNumberOfConstraints


    # ------------------------------------------------------------------------
    #       INITIALIZE DATASTRUCTURES
    # ------------------------------------------------------------------------
    probConvex              = 1.0
    time_SMC                = list()
    time_MILP_1core         = list()
    time_MILP_2core         = list()
    time_MILP_3core         = list()
    time_MILP_4core         = list()
    time_SMT                = list()
    
    
    for numberOfRealVars in numberOfRealVars_testCases:
        numberOfConstraints = numberOfBooleanConstraints
        print '\n********************************************************************'
        print '************* New Test Case: #Real Variabels = ', numberOfConstraints,' *****************'
        print '********************************************************************'
        # ------------------------------------------------------------------------
        #       GENERATE CONVEX CONSTRAINTS
        # ------------------------------------------------------------------------
        convexChoice            = np.random.choice([0,1], numberOfConstraints, p=[1- probConvex, probConvex])
        senseChoice             = np.random.choice([0,1], numberOfConstraints, p=[0.999, 0.001])


        numberOfConvexConstraints = int(numberOfRealVars/2)
        A                = sci.sparse.rand(numberOfConvexConstraints,i numberOfRealVars, density=0.9)
        A                = 10 * A - A.ceil()
        A                = A.toarray()
        b                = sci.sparse.rand(numberOfConvexConstraints, 1, density=1.0)
        b                = 100  * b - b.ceil()
        b                = b.toarray()

        constraints             = list()
        convexConstaintCounter   = 0

        # ------------------------------------------------------------------------
        #       CONNECT CONVEX AND BOOLEAN CONSTRAINTS
        # ------------------------------------------------------------------------
        for counter in range(0, numberOfConstraints):
            constraintLine      = content[documentStart + counter + 1].split(" ")
            # get rid of the trailing 0
            constraintLine      = constraintLine[0:-1]
            boolConstraint      = [int(i) for i in constraintLine]
            convexConstraint    = None
            if convexConstaintCounter < numberOfConvexConstraints:
                if convexChoice[counter] == 1:
                    sense           = 'L'
                    convexConstraint    = {'A': A[convexConstaintCounter,:], 'b':b[convexConstaintCounter], 'sense':sense}
                    convexConstaintCounter = convexConstaintCounter + 1

            constraints.append({'bool':boolConstraint, 'convex':convexConstraint})

        # ------------------------------------------------------------------------
        #       LAUNCH SMC SOLVER
        # ------------------------------------------------------------------------
        time_SMC.append(scalability_SMC(numberOfBoolVars, numberOfRealVars, numberOfConvexConstraints, constraints, 'IIS'))
        
        # Write SMC results
        print '--> Writing Execution Time Results to File'
        SMCfilename             = 'scalbality_reals_SMC_results.txt'
        SMCfile                 = open(SMCfilename, 'w')
        for item in time_SMC:
            if item > TIMEOUTE_LIMIT:
                SMCfile.write("TIME OUT\n")
            else:
                SMCfile.write("%s\n" % item)
        SMCfile.close()


        # ------------------------------------------------------------------------
        #       LAUNCH MILP SOLVER - 1 core
        # ------------------------------------------------------------------------
        time_MILP_1core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, 1))

        # Write CPLEX results - 1 core
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_reals_MILP1_results.txt'
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
        time_MILP_2core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, 2))

        # Write CPLEX results - 2 core2
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_reals_MILP2_results.txt'
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
        time_MILP_3core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, 3))

        # Write CPLEX results - 3 cores
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_reals_MILP3_results.txt'
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
        time_MILP_4core.append(scalability_CPLEX(numberOfBoolVars, numberOfRealVars, constraints, 4))

        # Write CPLEX results - 4 cores
        print '--> Writing Execution Time Results to File'
        MILPfilename            = 'scalbality_reals_MILP4_results.txt'
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
        time_SMT.append(scalability_Z3(numberOfBoolVars, numberOfRealVars, constraints))

        # Write SMT results
        print '--> Writing Execution Time Results to File'
        Z3filename             = 'scalbality_reals_Z3_results.txt'
        Z3file                 = open(Z3filename, 'w')
        for item in time_SMT:
            if item > TIMEOUTE_LIMIT:
                Z3file.write("TIME OUT\n")
            else:
                Z3file.write("%s\n" % item)
        Z3file.close()


        # ------------------------------------------------------------------------
        #       PRINT RESULTS
        # ------------------------------------------------------------------------

        print '-----> Execution Time <------'
        print 'time_SMC',time_SMC
        print 'time_MILP_1core',time_MILP_1core
        print 'time_MILP_2core',time_MILP_2core
        print 'time_MILP_3core',time_MILP_3core
        print 'time_MILP_4core',time_MILP_4core
        print 'time_SMT',time_SMT



    # ------------------------------------------------------------------------
    #       PLOT THE FINAL RESULTS
    # ------------------------------------------------------------------------
    
    plt.semilogy(numberOfRealVars_testCases, time_SMC, 'b',
                 numberOfRealVars_testCases, time_MILP_1core, 'r',
                 numberOfRealVars_testCases, time_MILP_2core, 'r--',
                 numberOfRealVars_testCases, time_MILP_3core, 'r:',
                 numberOfRealVars_testCases, time_MILP_4core, 'r.-',
                 numberOfRealVars_testCases, time_SMT, 'g'
                 )
    plt.legend(['SatEX','CPLEX 1 core', 'CPLEX 2 cores', 'CPLEX 3 cores', 'CPLEX 4 cores', 'Z3'])
    plt.xlabel('number of Real vars')
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


def scalability_CPLEX(numberOfBoolVars, numOfRealVars, constraints, numberOfCores):
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
    bigM                                = 1E16
    
    #----------------------------------------------------------------------------
    print '--> CPLEX: FEEDING CONSTRAINTS'
    #----------------------------------------------------------------------------
    convexIFCounter                     = 0
    for constraint in constraints:
        boolConstraint                  = constraint['bool']
        vars                            = ['b'+str(abs(i)-1)    for i in boolConstraint]
        signs                           = [np.sign(i)           for i in boolConstraint]
        numOfNegativeSigns              = sum(x < 0             for x in signs)
        mipSolver.linear_constraints.add(lin_expr   = [cplex.SparsePair(ind = vars,val =  signs)],
                                        rhs         = [1 -numOfNegativeSigns],
                                        senses      = ['G'])
        convexConstraint                = constraint['convex']
        if convexConstraint is not None:
            # the mip part takes the form of  a_i => (Ax \le b)
            # which is equal to Ax \le b + (1 - a_i)M
            # which is equal to (A + M) [x a_i]^T \le b + M
            A                           = np.append(convexConstraint['A'], bigM)
            b                           = convexConstraint['b'] + bigM
            cvars                       = rVars + [vars[-1]]
            mipSolver.linear_constraints.add(lin_expr   = [cplex.SparsePair(ind = cvars,val =  A)],
                                            rhs         = b,
                                            senses      = convexConstraint['sense']
                                        )
            convexIFCounter             = convexIFCounter + 1

    
    #----------------------------------------------------------------------------
    print '--> CPLEX: SEARCHING FOR A SOLUTION'
    #----------------------------------------------------------------------------
    start                               = timeit.default_timer()
    mipSolver.solve()
    end                                 = timeit.default_timer()
    time_smt_MILP                       = end - start
    
    print 'CPLEX time = ', time_smt_MILP
    
    return time_smt_MILP






#***************************************************************************************************
#***************************************************************************************************
#
#         SCALABILITY Test - SMC
#
#***************************************************************************************************
#***************************************************************************************************


def scalability_SMC(numberOfBoolVars, numOfRealVars, numOfConvIFClauses, constraints, strategy):
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
                                                     numberOfCores = 1)

    #----------------------------------------------------------------------------
    print '--> SMC: FEEDING CONSTRAINTS'
    #----------------------------------------------------------------------------
    convexIFCounter                     = 0
    
    for constraint in constraints:
        boolConstraint                  = constraint['bool']
        positiveVars                    = [abs(i)-1 for i in boolConstraint if i > 0]
        negativeVars                    = [abs(i)-1 for i in boolConstraint if i < 0]
        Z3Constraint                    = [solver.bVars[i] for i in positiveVars] + [NOT(solver.bVars[i]) for i in negativeVars]
        solver.addBoolConstraint(OR(*Z3Constraint))
        
        convexConstraint                = constraint['convex']
        if convexConstraint is not None:
            firstVar                    = abs(boolConstraint[-1])-1
            solver.setConvIFClause(LPClause(np.array([convexConstraint['A']]), convexConstraint['b'],
                                            solver.rVars,
                                            sense = convexConstraint['sense']),
                                    convexIFCounter)
            solver.addBoolConstraint(IMPLIES(solver.bVars[firstVar], solver.convIFClauses[convexIFCounter]))
            convexIFCounter             = convexIFCounter + 1


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

def scalability_Z3(numberOfBoolVars, numOfRealVars, constraints):
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
        positiveVars                    = [abs(i)-1 for i in boolConstraint if i > 0]
        negativeVars                    = [abs(i)-1 for i in boolConstraint if i < 0]
        Z3Constraint                    = [bVars[i] for i in positiveVars] + [NOT(bVars[i]) for i in negativeVars]
        
        SMTsolver.add(z3.Or(*Z3Constraint))
        convexConstraint                = constraint['convex']
        if convexConstraint is not None:
            firstVar                    = abs(boolConstraint[-1])-1
            A                           = convexConstraint['A']
            b                           = convexConstraint['b']
            linConst                    = sum([i * j for i,j in zip(rVars,A)] + [-1*b[0]])
            SMTsolver.add(z3.Implies(bVars[firstVar], linConst < 0))

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



if __name__ == "__main__":
    np.random.seed(0)
    scalability_main()



     

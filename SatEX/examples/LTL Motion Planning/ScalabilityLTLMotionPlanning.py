#  * --------------------------------------------------------------------------
#  * File: ScalabilityLTLMotionPlanning.py
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

import scipy as sci
from scipy.io import loadmat
from scipy.linalg import expm

import timeit
import math

import matplotlib.pyplot as plt


import sys

import numpy as np

import cplex as cplex
from cplex.exceptions import CplexError

#from collections import defaultdict



import csv

from itertools import izip

import signal


TIMEOUTE_LIMIT              = 600

# Register an handler for the timeout
def timeout_handler(signum, frame):
    print "------> Timeout <------"
    raise Exception("end of time")

# Register an handler for the keyboard inteeruption
def keyboard_handler(signal, frame):
    print 'You pressed Ctrl+C! --- Keep pressing more until you exit'
    sys.exit(0)

#***************************************************************************************************
#***************************************************************************************************
#
#         Robotic LTL Motion Planning Test
#
#***************************************************************************************************
#***************************************************************************************************


def robotic_LTL_main():
    signal.signal(signal.SIGINT, keyboard_handler)
    signal.signal(signal.SIGALRM, timeout_handler)
    
    x           = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]#, [3, 4,  4, 3]#,  4, 5, 5,  5, 5]
    y           = [4,4,4,4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]# [2, 2,  3, 3]#,  3, 3, 4, 5, 5]
    z           = [4,4,4,4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]#[2, 2,  2, 3]#,  3, 3, 3, 4, 5]


    time_IIS                = list()
    time_SSF                = list()
    time_PREFIX             = list()
    time_MILP_1core         = list()

    for counter in range(len(x)):
        xMaxTick                        = x[counter]
        yMaxTick                        = y[counter]
        zMaxTick                        = z[counter]
        # ---------- LOAD WORKSPACE -------------------------------------------------------------------
        partitions, obstaclPartitions, tag2PartitionCounterDic, wkspDimensions  = robotic_LTL_workspace(xMaxTick, yMaxTick, zMaxTick)
        numOfPartitions                 = len(partitions)
        
        #print 'obstaclPartitions', obstaclPartitions, 'partitions',partitions
    
        #return
        # ---------- PHYSICAL CONSTRAINTS ----------------------------------------------------------------
        x0                              = 0.0
        y0                              = 0.0
        z0                              = 0.0
        p0                              = 0.0
    
        X_0                             = [x0, 0.0, 0.0]
        Y_0                             = [y0, 0.0, 0.0]
        Z_0                             = [z0, 0.0, 0.0]
        P_0                             = [p0, 0.0]
    
        XDot_f                          = [0.0, 0.0]
        YDot_f                          = [0.0, 0.0]
        ZDot_f                          = [0.0, 0.0]
        PDot_f                          = [0.0]
    
        xf                              = None
        yf                              = None
        zf                              = None #0.4
        pf                              = 0.0
    
    
        startRegion                     = tag2PartitionCounterDic['111']
        #endRegion                       = tag2PartitionCounterDic['113']
        endRegion                       = tag2PartitionCounterDic[str(xMaxTick-1) + str(yMaxTick-1) + str(zMaxTick-1)]
    

        accBound                        = 1.0
        smoothBound                     = 1.0

        xMin                            = wkspDimensions[0]
        xMax                            = wkspDimensions[1]
        yMin                            = wkspDimensions[2]
        yMax                            = wkspDimensions[3]
        zMin                            = wkspDimensions[4]
        zMax                            = wkspDimensions[5]
    
    
        Ts                              = 1.0/1.0
        maxHorizon                      = 1000
        final_time                      = 0.0
    
    
    
    
        # ---------- MILP ----------------------------------------------------------------
        
        print '\n======================================================='
        print '      TEST CPLEX'
        print '            #Real Varabiles = \n',
        print '            #Bool Varabiles = \n',
        print '            #Bool Constraints = \n',
        print '            #Cores = \n',
        print '=======================================================\n'

        counterExamples                 = list()
        start                           = timeit.default_timer()
        for horizon in range(xMaxTick, maxHorizon):
            dwellTime                       = [1] + [1]*horizon                # dwell time starts with one (for initial condition) then dwell
                                                                        # for each step in horizon

            traj_1Hz_MILP               = robotic_LTL_segment_MILP(
                                                partitions, obstaclPartitions, numOfPartitions,
                                                horizon, counterExamples, dwellTime,
                                                accBound, smoothBound, Ts, startRegion, endRegion,
                                                X_0, Y_0, Z_0, P_0, XDot_f, YDot_f, ZDot_f, PDot_f,
                                                x_f = xf, y_f = yf, z_f = zf, p_f = pf
                                        )
                        
            if traj_1Hz_MILP['xTraj']:     # if no more counter examples then return
                break
        end                             = timeit.default_timer()
        time_smt_MILP                   = end - start
        print 'Time CPLEX = ', time_smt_MILP
        time_MILP_1core.append(time_smt_MILP)
        
        # Write execution time results
        print '--> Writing Execution Time Results to File'
        MILPfilename             = 'scalbality_motionplanning_MILP_results.txt'
        MILPfilename                 = open(MILPfilename, 'w')
        for item in time_MILP_1core:
            if item > TIMEOUTE_LIMIT:
                MILPfilename.write("TIME OUT\n")
            else:
                MILPfilename.write("%s\n" % item)
        MILPfilename.close()
        # ---------- IIS STRATEGY ----------------------------------------------------------------
        
        print '\n======================================================='
        print '      TEST SMC'
        print '            UNSAT Certificate = IIS \n',
        print '=======================================================\n'


        strategy                        = 'IIS'
        counterExamples                 = list()

        signal.alarm(TIMEOUTE_LIMIT)
        try:
            start                           = timeit.default_timer()
            for horizon in range(xMaxTick, maxHorizon):
                print 'horizon', horizon
                dwellTime                       = [1] + [1]*horizon                # dwell time starts with one (for initial condition) then dwell
                                                                        # for each step in horizon

                traj_1Hz = robotic_LTL_segment(strategy, partitions, obstaclPartitions, numOfPartitions, horizon, counterExamples, dwellTime, accBound, smoothBound, Ts, startRegion, endRegion, X_0, Y_0, Z_0, P_0, XDot_f, YDot_f, ZDot_f, PDot_f, x_f = xf, y_f = yf, z_f = zf, p_f = pf
                )


                #counterExamples = counterExamples + traj_1Hz['counterExamples']
                if traj_1Hz['xTraj']:     # if no more counter examples then return
                    break
        except Exception, exc:
            print 'Time out exception'
        finally:

            end                             = timeit.default_timer()
            time_smc_IIS                   = end - start
            #print 'Time PREFIX = ', time_smt_PREFIX, '#Counter Exmples = ', len(traj_1Hz['counterExamples']), traj_1Hz['xTraj']
            time_IIS.append(time_smc_IIS)
            print 'Time IIS = ', time_smc_IIS
            
            # Write execution time results
            print '--> Writing Execution Time Results to File'
            SMCfilename             = 'scalbality_motionplanning_IIS_results.txt'
            SMCfilename                 = open(SMCfilename, 'w')
            for item in time_IIS:
                if item > TIMEOUTE_LIMIT:
                    SMCfilename.write("TIME OUT\n")
                else:
                    SMCfilename.write("%s\n" % item)
            SMCfilename.close()
            
        
        # ---------- SSF STRATEGY ----------------------------------------------------------------
        '''
        print '\n======================================================='
        print '      TEST SMC'
        print '            UNSAT Certificate = SSF \n',
        print '=======================================================\n'


        strategy                        = 'SMC'
        counterExamples                 = list()


        start                           = timeit.default_timer()
        for horizon in range(xMaxTick, maxHorizon):
            print 'horizon', horizon
            dwellTime                       = [1] + [1]*horizon                # dwell time starts with one (for initial condition) then dwell
                                                                        # for each step in horizon

            traj_1Hz = robotic_LTL_segment(strategy, partitions, obstaclPartitions, numOfPartitions, horizon, counterExamples, dwellTime, accBound, smoothBound, Ts, startRegion, endRegion, X_0, Y_0, Z_0, P_0, XDot_f, YDot_f, ZDot_f, PDot_f, x_f = xf, y_f = yf, z_f = zf, p_f = pf
            )


            #counterExamples = counterExamples + traj_1Hz['counterExamples']
            if traj_1Hz['xTraj']:     # if no more counter examples then return
                break

        end                             = timeit.default_timer()
        time_smc_SSF                   = end - start
        #print 'Time PREFIX = ', time_smt_PREFIX, '#Counter Exmples = ', len(traj_1Hz['counterExamples']), traj_1Hz['xTraj']
        time_SSF.append(time_smc_SSF)

        '''

        # ---------- PREFIX STRATEGY ----------------------------------------------------------------
        
        print '\n======================================================='
        print '      TEST SMC'
        print '            UNSAT Certificate = PREFIX \n',
        print '=======================================================\n'


        strategy                        = 'PREFIX'
        counterExamples                 = list()

        signal.alarm(TIMEOUTE_LIMIT)
        try:
            start                           = timeit.default_timer()
            for horizon in range(xMaxTick, maxHorizon):
                print 'horizon', horizon
                dwellTime                       = [1] + [1]*horizon                # dwell time starts with one (for initial condition) then dwell
                                                                        # for each step in horizon

                traj_1Hz = robotic_LTL_segment(strategy, partitions, obstaclPartitions, numOfPartitions, horizon, counterExamples, dwellTime, accBound, smoothBound, Ts, startRegion, endRegion, X_0, Y_0, Z_0, P_0, XDot_f, YDot_f, ZDot_f, PDot_f, x_f = xf, y_f = yf, z_f = zf, p_f = pf
                )


                #counterExamples = counterExamples + traj_1Hz['counterExamples']
                if traj_1Hz['xTraj']:     # if no more counter examples then return
                    break

        except Exception, exc:
            print 'Time out exception'
        finally:
            end                             = timeit.default_timer()
            time_smc_PREFIX                 = end - start
            #print 'Time PREFIX = ', time_smt_PREFIX, '#Counter Exmples = ', len(traj_1Hz['counterExamples']), traj_1Hz['xTraj']
            time_PREFIX.append(time_smc_PREFIX)
            print 'Time PREFIX = ', time_smc_PREFIX
        
            # Write execution time results
            print '--> Writing Execution Time Results to File'
            SMCfilename             = 'scalbality_motionplanning_PREFIX_results.txt'
            SMCfilename                 = open(SMCfilename, 'w')
            for item in time_PREFIX:
                if item > TIMEOUTE_LIMIT:
                    SMCfilename.write("TIME OUT\n")
                else:
                    SMCfilename.write("%s\n" % item)
            SMCfilename.close()
        

        print '-----> Execution Time <------'
        print 'time_IIS',time_IIS
        print 'time_SSF',time_SSF
        print 'time_PREFIX',time_PREFIX
        print 'time_MILP_1core',time_MILP_1core
        #print 'time_MILP_2core',time_MILP_2core
        #print 'time_MILP_3core',time_MILP_3core
        #print 'time_MILP_4core',time_MILP_4core
        #print 'time_SMT',time_SMT



#***************************************************************************************************
#***************************************************************************************************
#
#         Robotic LTL Motion Planning Test - MILP
#
#***************************************************************************************************
#***************************************************************************************************

def robotic_LTL_segment_MILP(partitions, obstaclPartitions, numOfPartitions, horizon, counterExamples, dwellTime,
                        accBound, smoothBound, Ts, startRegion, endRegion,
                        X_0, Y_0, Z_0, P_0, XDot_f, YDot_f, ZDot_f, PDot_f, x_f = None, y_f = None, z_f = None, p_f = None
                        ):

    # ---------- STATE VAR INDEX ----------------------------------------------------------------
    indexInputX                     = 0
    indexInputY                     = 1
    indexInputZ                     = 2
    indexInputPsi                   = 3
    U                               = [indexInputX, indexInputX, indexInputY, indexInputPsi]

    indexX                          = 4
    indexXDot                       = 5
    indexXDDot                      = 6
    indexXDDDot                     = 7
    X                               = [indexX, indexXDot, indexXDDot, indexXDDDot]

    indexY                          = 8
    indexYDot                       = 9
    indexYDDot                      = 10
    indexYDDDot                     = 11
    Y                               = [indexY, indexYDot, indexYDDot, indexYDDDot]

    indexZ                          = 12
    indexZDot                       = 13
    indexZDDot                      = 14
    indexZDDDot                     = 15
    Z                               = [indexZ, indexZDot, indexZDDot, indexZDDDot]

    indexPsi                        = 16
    indexPsiDot                     = 17
    Psi                             = [indexPsi, indexPsiDot]

    state                           = range(0,18)


    # ---------- INSTANTIATE SOLVER ----------------------------------------------------------------
    mipSolver = cplex.Cplex()
    mipSolver.parameters.threads.set(1)
    mipSolver.parameters.timelimit.set(TIMEOUTE_LIMIT)
    
    totalDwell                      = sum([i for i in dwellTime])
    
    numOfBooleanVars                = numOfPartitions * horizon * (totalDwell * 4) # 4 states
    numberOfStates                  = 3 * 4 + 2         # x-y-z are 4th order integrators, psi 2nd order integrator
    numberOfStates                  = numberOfStates * totalDwell
    
    
    numberOfInputs                  = 4                 # one input for each integrator chain
    numberOfInputs                  = numberOfInputs  * totalDwell

    numOfRealVars                   = numberOfStates + numberOfInputs
    numOfConvIFClauses              = numOfPartitions * horizon



    rVars                               = ['x'+str(i) for i in range(0, numOfRealVars)]
    bVars                               = ['b'+str(i) for i in range(0, numOfBooleanVars)]
    
    mipSolver.variables.add(names = bVars, types = [mipSolver.variables.type.binary] * len(bVars))
    
    mipSolver.variables.add(names       =   rVars,
                            ub          =   [cplex.infinity] * len(rVars),
                            lb          =   [-cplex.infinity] * len(rVars)
    )
    
    
    bigM                                = 1E15
    
    #mipSolver.parameters.mip.tolerances.absmipgap.set(1E-16)
    
    
    # ---------- ADD ADJACENY MATRIX ----------------------------------------------------------------
    # for the $k$th time step , generate clause ensuring that only one partition is active
    # P_{1,k} + P_{2,k} + P_{p,k} = 1
    # where $p$ = number of partitions
    for horizonCounter in range(0, horizon):
        indexShift                      = horizonCounter*numOfPartitions
        indecies                        = [bVars[i+indexShift] for i in range (0, numOfPartitions)]
        #print 'indecies', indecies
        mipSolver.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = indecies, val =  [1]*len(indecies))],
                                        rhs      = [1],
                                        senses   = ['E'],
                                        names    = ['selOne_'+str(horizonCounter)]
                                        )
        #solver.addBoolConstraint(sum([solver.bVars[i+indexShift] for i in range (0, numOfPartitions)]) == 1)

    # unroll the state machine representing the workspace adjacency
    # For the $i$ th partition, generate the follwoing clause
    # P_{i,k} => P_{adj1,k+1} OR P_{adj2,k+1} OR ... P_{adjn,k+1}
    # NOT P_{i,k} OR P_{adj1,k+1} OR P_{adj2,k+1} OR ... P_{adjn,k+1}
    
    # to encode b_1 OR b_2 as ILP:  b_1 + b_2 \ge 1
    # hence to encode NOT b_1 OR b_2 <==>  (1 - b_1) + b_2 \ge 1 <==> - b_1 + b_2 \ge 0
    for horizonCounter in range(0, horizon-1):
        indexShift                  = horizonCounter*numOfPartitions
        indexShift_plus             = (horizonCounter+1)*numOfPartitions
        for counter in range(0,numOfPartitions):
            adjacent                = partitions[counter]['adjacent']
            anticedent              = [bVars[i + indexShift_plus] for i in adjacent]
            sourceRegion            = [bVars[counter+indexShift]]
            #print ['adj_'+str(counter) +'_'+str(horizonCounter)], 'anticedent of', sourceRegion, 'is ', anticedent
            mipSolver.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = sourceRegion + anticedent,
                                                                            val =  [-1] + [1]*len(anticedent))],
                                            rhs      = [0],
                                            senses   = ['G'],
                                            names    = ['adj_'+str(counter) +'_'+str(horizonCounter)]
                                            )

            #solver.addBoolConstraint(
            #    IMPLIES(    solver.bVars[counter+indexShift], OR(*anticedent)
            #                # the * used to unpack the python list to arguments
            #    )
            #)
    
    # ---------- SPECIFICATIONS -------------------------------------------------------------------
    
    #mipSolver.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [bVars[startRegion]],val =  [1])],
    #                                        rhs      = [1],
    #                                        senses   = ['E'],
    #                                        names   = ['start_region'])
    #print 'start = ', [bVars[startRegion]]
    
    #solver.addBoolConstraint(solver.bVars[startRegion])
    
    goalIndecies                    = [i*numOfPartitions + endRegion for i in range(0,horizon)]
    goals                           = [bVars[i] for i in goalIndecies]
    #goals                           = [goals[-1]]
    
    mipSolver.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = goals,val =  [1]*len(goals))],
                                            rhs      = [1],
                                            senses   = ['E'],
                                            names    = ['goal_region'])

    #solver.addBoolConstraint(OR(*goals))

    # avoid obstacles regions
    
    badIndecies                     = [bVars[i*numOfPartitions + j] for i in range(0,horizon-1) for j in obstaclPartitions]
    #print 'bad index = ', badIndecies
    for badIndex in badIndecies:
        mipSolver.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [badIndex],val =  [1.0])],
                                            rhs      = [0],
                                            senses   = ['E'],
                                            names    = ['avoid_obst_'+str(badIndex)])
    
        #solver.addBoolConstraint(NOT(solver.bVars[badIndex]))


    # ---------- ADD INPUT AND DYNAMICS CONSTRAINTS ------------------------------------------------
    #input constraints (can be added here if needed)


    # discrete dynamics of 4th order chain of integrators
    A4_cont                         = np.array([[0,1,0,0],[0,0,1,0], [0,0,0,1], [0,0,0,0] ])       #integrators
    A4_disc                         = expm(A4_cont*Ts)

    B4_cont                         = np.array([[0,0,0,1]])
    B4_disc                         = np.array([[0,0,0,Ts]])

    I4                              = np.array([[1,0,0,0], [0,1,0,0],[0,0,1,0], [0,0,0,1] ])       #identity

    A4_LP                           = np.concatenate((I4, -1*A4_disc, -1*B4_disc.T), axis = 1)
    b4_LP                           = [0]*4
    
    # discrete dynamics of 2nd order chain of integrators
    A2_cont                         = np.array([[0,1],[0,0] ])       #integrators
    A2_disc                         = expm(A2_cont*Ts)

    B2_cont                         = np.array([[0,1]])
    B2_disc                         = np.array([[0,Ts]])

    I2                              = np.array([[1,0], [0,1] ])       #identity

    A2_LP                           = np.concatenate((I2, -1*A2_disc, -1*B2_disc.T), axis = 1)
    b2_LP                           = [0]*2

    # initial state constraints
    vars                            = [rVars[i] for i in [indexX, indexXDot, indexXDDot]]
    initialStateConstraint          = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), X_0, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = initialStateConstraint['lin_expr'],
                                    senses      = initialStateConstraint['senses'],
                                    rhs         = initialStateConstraint['rhs'],
                                    names       = ['initState_X_'+str(i) for i in range(0,len(initialStateConstraint['rhs']))]
                          )
    #solver.addConvConstraint(initialStateConstraint)

    vars                            = [rVars[i] for i in [indexY, indexYDot, indexYDDot]]
    initialStateConstraint          = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), Y_0, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = initialStateConstraint['lin_expr'],
                                    senses      = initialStateConstraint['senses'],
                                    rhs         = initialStateConstraint['rhs'],
                                    names       = ['initState_Y_'+str(i) for i in range(0,len(initialStateConstraint['rhs']))]
                          )
#    solver.addConvConstraint(initialStateConstraint)

    vars                            = [rVars[i] for i in [indexZ, indexZDot, indexZDDot]]
    initialStateConstraint          = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), Z_0, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = initialStateConstraint['lin_expr'],
                                    senses      = initialStateConstraint['senses'],
                                    rhs         = initialStateConstraint['rhs'],
                                    names       = ['initState_Z_'+str(i) for i in range(0,len(initialStateConstraint['rhs']))]
                          )
#    solver.addConvConstraint(initialStateConstraint)

    vars                            = [rVars[i] for i in Psi]
    initialStateConstraint          = LPClause(np.array([[1,0], [0,1] ]), P_0, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = initialStateConstraint['lin_expr'],
                                    senses      = initialStateConstraint['senses'],
                                    rhs         = initialStateConstraint['rhs'],
                                    names       = ['initState_P_'+str(i) for i in range(0,len(initialStateConstraint['rhs']))]
                          )
#    solver.addConvConstraint(initialStateConstraint)
    

    
    #for horizonCounter in range(0, (horizon * totalDwell) - 1):
    for horizonCounter in range(0, totalDwell - 1):
        # Dynamics Constraints, x-axis
        indexShift                  = (horizonCounter)* len(state)
        indexShiftplus              = (horizonCounter+1)* len(state)
        
        xVars                       = [rVars[i + indexShift] for i in X]
        xplusVars                   = [rVars[i + indexShiftplus] for i in X]
        vars                        = xplusVars + xVars + [rVars[indexInputX + indexShift]]
        
        dynamicConstraint           = LPClause(A4_LP, b4_LP, vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = dynamicConstraint['lin_expr'],
                                    senses      = dynamicConstraint['senses'],
                                    rhs         = dynamicConstraint['rhs'],
                                    names       = ['dyn_X_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(dynamicConstraint['rhs']))]
                          )
        
#        solver.addConvConstraint(dynamicConstraint)
        
        # Dynamics Constraints, y-axis
        yVars                       = [rVars[i + indexShift] for i in Y]
        yplusVars                   = [rVars[i + indexShiftplus] for i in Y]
        vars                        = yplusVars + yVars + [rVars[indexInputY + indexShift]]
        
        
        dynamicConstraint = LPClause(A4_LP, b4_LP, vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = dynamicConstraint['lin_expr'],
                                    senses      = dynamicConstraint['senses'],
                                    rhs         = dynamicConstraint['rhs'],
                                    names       = ['dyn_Y_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(dynamicConstraint['rhs']))]
                          )
        
#        solver.addConvConstraint(dynamicConstraint)


        # Dynamics Constraints, z-axis
        zVars                       = [rVars[i + indexShift] for i in Z]
        zplusVars                   = [rVars[i + indexShiftplus] for i in Z]
        vars                        = zplusVars + zVars + [rVars[indexInputZ + indexShift]]
        
        dynamicConstraint = LPClause(A4_LP, b4_LP, vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = dynamicConstraint['lin_expr'],
                                    senses      = dynamicConstraint['senses'],
                                    rhs         = dynamicConstraint['rhs'],
                                    names       = ['dyn_Z_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(dynamicConstraint['rhs']))]
                          )
        
#        solver.addConvConstraint(dynamicConstraint)
        

        # Dynamics Constraints, psi-axis
        pVars                       = [rVars[i + indexShift] for i in Psi]
        pplusVars                   = [rVars[i + indexShiftplus] for i in Psi]
        vars                        = pplusVars + pVars + [rVars[indexInputPsi + indexShift]]

        dynamicConstraint = LPClause(A2_LP, b2_LP, vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = dynamicConstraint['lin_expr'],
                                    senses      = dynamicConstraint['senses'],
                                    rhs         = dynamicConstraint['rhs'],
                                    names       = ['dyn_P_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(dynamicConstraint['rhs']))]
                          )
        
#        solver.addConvConstraint(dynamicConstraint)



        # acceleration must be between -3 m/s^2 and 3 m/s^2
        
        mipSolver.variables.set_upper_bounds(rVars[indexXDDot + indexShift], 0.5 * accBound)
        mipSolver.variables.set_lower_bounds(rVars[indexXDDot + indexShift], -0.5 * accBound)
        
        
        
        mipSolver.variables.set_upper_bounds(rVars[indexYDDot + indexShift], 0.5 * accBound)
        mipSolver.variables.set_lower_bounds(rVars[indexYDDot + indexShift], -0.5 * accBound)
        
        
        mipSolver.variables.set_upper_bounds(rVars[indexZDDot + indexShift], 0.5 * accBound)
        mipSolver.variables.set_lower_bounds(rVars[indexZDDot + indexShift], -0.5 * accBound)
        
        
        #C4 smootheness constraint
        


        
        vars                        = [rVars[i] for i in [indexXDDDot + indexShiftplus, indexXDDDot + indexShift]]
        smoothConstraint            = LPClause(np.array([[1.0, -1.0]]), [smoothBound], vars)
        mipSolver.linear_constraints.add(
                                    lin_expr    = smoothConstraint['lin_expr'],
                                    senses      = smoothConstraint['senses'],
                                    rhs         = smoothConstraint['rhs'],
                                    names       = ['smooth_X1_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(smoothConstraint['rhs']))]
                          )
#        solver.addConvConstraint(smoothConstraint)
        
        smoothConstraint            = LPClause(np.array([[-1.0, 1.0]]), [smoothBound], vars)
        mipSolver.linear_constraints.add(
                                    lin_expr    = smoothConstraint['lin_expr'],
                                    senses      = smoothConstraint['senses'],
                                    rhs         = smoothConstraint['rhs'],
                                    names       = ['smooth_X2_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(smoothConstraint['rhs']))]
                          )
#        solver.addConvConstraint(smoothConstraint)

        vars                        = [rVars[i] for i in [indexYDDDot + indexShiftplus, indexYDDDot + indexShift]]

        smoothConstraint            = LPClause(np.array([[1.0, -1.0]]), [smoothBound], vars)
        mipSolver.linear_constraints.add(
                                    lin_expr    = smoothConstraint['lin_expr'],
                                    senses      = smoothConstraint['senses'],
                                    rhs         = smoothConstraint['rhs'],
                                    names       = ['smooth_Y1_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(smoothConstraint['rhs']))]
                          )
        
#        solver.addConvConstraint(smoothConstraint)
        
        smoothConstraint            = LPClause(np.array([[-1.0, 1.0]]), [smoothBound], vars)
        mipSolver.linear_constraints.add(
                                    lin_expr    = smoothConstraint['lin_expr'],
                                    senses      = smoothConstraint['senses'],
                                    rhs         = smoothConstraint['rhs'],
                                    names       = ['smooth_Y2_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(smoothConstraint['rhs']))]
                          )
        
#        solver.addConvConstraint(smoothConstraint)
        
        vars                        = [rVars[i] for i in [indexZDDDot + indexShiftplus, indexZDDDot + indexShift]]
    
        smoothConstraint            = LPClause(np.array([[1.0, -1.0]]), [smoothBound], vars)
        mipSolver.linear_constraints.add(
                                    lin_expr    = smoothConstraint['lin_expr'],
                                    senses      = smoothConstraint['senses'],
                                    rhs         = smoothConstraint['rhs'],
                                    names       = ['smooth_Z1_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(smoothConstraint['rhs']))]
                          )
        
#        solver.addConvConstraint(smoothConstraint)
        
        smoothConstraint            = LPClause(np.array([[-1.0, 1.0]]), [smoothBound], vars)
        mipSolver.linear_constraints.add(
                                    lin_expr    = smoothConstraint['lin_expr'],
                                    senses      = smoothConstraint['senses'],
                                    rhs         = smoothConstraint['rhs'],
                                    names       = ['smooth_Z2_'+ str(horizonCounter) + '_' + str(i) for i in range(0,len(smoothConstraint['rhs']))]
                          )
    
#        solver.addConvConstraint(smoothConstraint)



    # velocity, acceleration, jerk should be equal to zero at the end.
    #indexShift                      = ((horizon * totalDwell) - 1)* len(state)
    indexShift                      = (totalDwell - 1)* len(state)
    
    
    vars                            = [rVars[i] for i in [indexXDot + indexShift, indexXDDot + indexShift, indexXDDDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), XDot_f, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )
#    solver.addConvConstraint(finalStateConstraint)

    vars                            = [rVars[i] for i in [indexYDot + indexShift, indexYDDot + indexShift, indexYDDDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), YDot_f, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )
#    solver.addConvConstraint(finalStateConstraint)

    vars                            = [rVars[i] for i in [indexZDot + indexShift, indexZDDot + indexShift, indexZDDDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), ZDot_f, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )
#    solver.addConvConstraint(finalStateConstraint)

    vars                            = [rVars[i] for i in [indexPsiDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1] ]), PDot_f, vars,sense = "E")
    mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )
#    solver.addConvConstraint(finalStateConstraint)




    # if assigned, we add constraint on the final position
    if x_f is not None:
        vars                            = [rVars[indexX + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [x_f], vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )
#        solver.addConvConstraint(finalStateConstraint)
        
    if y_f is not None:
        vars                            = [rVars[indexY + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [y_f], vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )

#        solver.addConvConstraint(finalStateConstraint)
        

    if z_f is not None:
        vars                            = [rVars[indexZ + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [z_f], vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )

#        solver.addConvConstraint(finalStateConstraint)

    if p_f is not None:
        vars                            = [rVars[indexPsi + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [p_f], vars,sense = "E")
        mipSolver.linear_constraints.add(
                                    lin_expr    = finalStateConstraint['lin_expr'],
                                    senses      = finalStateConstraint['senses'],
                                    rhs         = finalStateConstraint['rhs'],
                          )

#        solver.addConvConstraint(finalStateConstraint)


    # ---------- REGIONS CONSTRAINTS ----------------------------------------------------------------
    for horizonCounter in range(0, horizon):
        indexShift                  = horizonCounter * numOfPartitions
        for partitionCounter in range(0,numOfPartitions):
            #solver.addBoolConstraint(IMPLIES(solver.bVars[partitionCounter + indexShift], solver.convIFClauses[partitionCounter + indexShift]))
            # skip initial condition and starts first step in the plan

            #stateIndexShift             = (horizonCounter * totalDwell * len(state))
            stateIndexShift             = sum([dwellTime[i] for i in range(0,horizonCounter+1)]) * len(state)
            #print 'stateIndexShift', stateIndexShift, len(rVars)
            
            vars                        = [rVars[i + stateIndexShift] for i in [indexX, indexY, indexZ]]
            
            partition_A                 = partitions[partitionCounter]['A']
            partition_b                 = partitions[partitionCounter]['b']
            
            
            for dwellCounter in range(1, dwellTime[horizonCounter+1]):
                #stateIndexShift         = (horizonCounter * totalDwell * len(state))  + (dwellCounter * len(state) )
                stateIndexShift         = (sum([dwellTime[i] for i in range(0,horizonCounter+1)]) * len(state))  + (dwellCounter * len(state) )
                
                partition_A_new_rows    = np.concatenate(
                                            (np.zeros((   np.shape(partitions[partitionCounter]['A'])[0], np.shape(partition_A)[1] )),
                                            partitions[partitionCounter]['A']),
                                            axis = 1)
                
                partition_A             = np.concatenate((partition_A,
                                                np.zeros((np.shape(partition_A)[0], np.shape(partitions[partitionCounter]['A'])[1]))),
                                                axis = 1)
                partition_A             = np.concatenate((partition_A, partition_A_new_rows))
                partition_b             = np.concatenate((partition_b,partitions[partitionCounter]['b']))
                
                vars                    = vars + [rVars[i + stateIndexShift] for i in [indexX, indexY, indexZ]]
            
            regionconstraint            = LPClause(partition_A, partition_b, vars)
            # use big M to slack the constraint whenever the Boolean variable is assigned to zero
            
            numOfRows                   = len(partition_b)
            slacked_regionconstraint    = list()
            slacked_rhs                 = list()
            slackVariable               = bVars[partitionCounter + indexShift]
            AAA                         = partition_A
            BBB                         = partition_b
            #print 'vars',vars, 'slackVariable', slackVariable
            #print 'partition_A', partition_A
            for counter in range(0, numOfRows):
                A_row                   = AAA[counter,:]
                slacked_regionconstraint.append(cplex.SparsePair(
                            ind                     = vars + [slackVariable],
                            val                     = np.append(A_row, bigM)
                            ))
                slacked_rhs.append(BBB[counter] + bigM)
            
            #print ['regConstr_'+str(partitionCounter) +'_' + str(horizonCounter)], slacked_regionconstraint, slacked_rhs
            #if partitionCounter in [21, 63, 108, 156]:
            #    print '----> added'
            mipSolver.linear_constraints.add(
                                    lin_expr    = slacked_regionconstraint,
                                    senses      = regionconstraint['senses'],
                                    rhs         = slacked_rhs,
                                    names       = ['regConstr_'+str(partitionCounter) +'_' + str(horizonCounter)+ '_'+str(i) for i in range(0, numOfRows)]
                          )
            
#            solver.setConvIFClause(regionconstraint, partitionCounter + indexShift)


    # ---------- RUN AND GET RESULTS ----------------------------------------------------------------
    mipSolver.solve()

    print(mipSolver.solution.get_status(), mipSolver.solution.get_status_string())


    xTraj = []
    if mipSolver.solution.get_status() == 101:
        traj = [i for i, x in enumerate(mipSolver.solution.get_values(bVars)) if x == 1.0]
        horizonShifts = [i * numOfPartitions for i in range(0,horizon)]
        #print traj
        
        traj = [a_i - b_i for a_i, b_i in zip(traj, horizonShifts)]
        #print traj
        
        stateTraj = mipSolver.solution.get_values(rVars)
        #print xTraj
        regionTraj = mipSolver.solution.get_values(bVars)

        inputXTraj                      = [stateTraj[i + indexInputX]   for i in  range(0, len(state) * totalDwell, len(state))]
        inputYTraj                      = [stateTraj[i + indexInputY]   for i in  range(0, len(state) * totalDwell, len(state))]
        inputZTraj                      = [stateTraj[i + indexInputZ]   for i in  range(0, len(state) * totalDwell, len(state))]
        inputPTraj                      = [stateTraj[i + indexInputPsi] for i in  range(0, len(state) * totalDwell, len(state))]

        xTraj                           = [stateTraj[i + indexX]        for i in  range(0, len(state) * totalDwell, len(state))]
        xDotTraj                        = [stateTraj[i + indexXDot]     for i in  range(0, len(state) * totalDwell, len(state))]
        xDDotTraj                       = [stateTraj[i + indexXDDot]    for i in  range(0, len(state) * totalDwell, len(state))]
        xDDDotTraj                      = [stateTraj[i + indexXDDDot]   for i in  range(0, len(state) * totalDwell, len(state))]

        yTraj                           = [stateTraj[i + indexY]        for i in  range(0, len(state) * totalDwell, len(state))]
        yDotTraj                        = [stateTraj[i + indexYDot]     for i in  range(0, len(state) * totalDwell, len(state))]
        yDDotTraj                       = [stateTraj[i + indexYDDot]    for i in  range(0, len(state) * totalDwell, len(state))]
        yDDDotTraj                      = [stateTraj[i + indexYDDDot]   for i in  range(0, len(state) * totalDwell, len(state))]

        zTraj                           = [stateTraj[i + indexZ]        for i in  range(0, len(state) * totalDwell, len(state))]
        zDotTraj                        = [stateTraj[i + indexZDot]     for i in  range(0, len(state) * totalDwell, len(state))]
        zDDotTraj                       = [stateTraj[i + indexZDDot]    for i in  range(0, len(state) * totalDwell, len(state))]
        zDDDotTraj                      = [stateTraj[i + indexZDDDot]   for i in  range(0, len(state) * totalDwell, len(state))]


        pTraj                           = [stateTraj[i + indexPsi]      for i in  range(0, len(state) * totalDwell, len(state))]
        pDotTraj                        = [stateTraj[i + indexPsiDot]   for i in  range(0, len(state) * totalDwell, len(state))]


        regionSeq_augmented             = regionTraj
        seqTraj                         = [y for i in range(0,horizon+1) for y in [regionSeq_augmented[i]] * dwellTime[i]]


    traj                            = { 'inputXTraj': [],
                                        'inputYTraj': [],
                                        'inputZTraj': [],
                                        'inputPTraj': [],
                                        'xTraj'     : xTraj,
                                        'xDotTraj'  : [],
                                        'xDDotTraj' : [],
                                        'xDDDotTraj': [],
                                        'yTraj'     : [],
                                        'yDotTraj'  : [],
                                        'yDDotTraj' : [],
                                        'yDDDotTraj': [],
                                        'zTraj'     : [],
                                        'zDotTraj'  : [],
                                        'zDDotTraj' : [],
                                        'zDDDotTraj': [],
                                        'pTraj'     : [],
                                        'pDotTraj'  : [],
                                        'seqTraj'   : [],
                                        'counterExamples' : []
                                        }
    return  traj









#***************************************************************************************************
#***************************************************************************************************
#
#         Robotic LTL Motion Planning Test - SMC
#
#***************************************************************************************************
#***************************************************************************************************

def robotic_LTL_segment(strategy, partitions, obstaclPartitions, numOfPartitions, horizon, counterExamples, dwellTime,
                        accBound, smoothBound, Ts, startRegion, endRegion,
                        X_0, Y_0, Z_0, P_0, XDot_f, YDot_f, ZDot_f, PDot_f, x_f = None, y_f = None, z_f = None, p_f = None
                        ):
    

    # ---------- STATE VAR INDEX ----------------------------------------------------------------
    indexInputX                     = 0
    indexInputY                     = 1
    indexInputZ                     = 2
    indexInputPsi                   = 3
    U                               = [indexInputX, indexInputX, indexInputY, indexInputPsi]

    indexX                          = 4
    indexXDot                       = 5
    indexXDDot                      = 6
    indexXDDDot                     = 7
    X                               = [indexX, indexXDot, indexXDDot, indexXDDDot]

    indexY                          = 8
    indexYDot                       = 9
    indexYDDot                      = 10
    indexYDDDot                     = 11
    Y                               = [indexY, indexYDot, indexYDDot, indexYDDDot]

    indexZ                          = 12
    indexZDot                       = 13
    indexZDDot                      = 14
    indexZDDDot                     = 15
    Z                               = [indexZ, indexZDot, indexZDDot, indexZDDDot]

    indexPsi                        = 16
    indexPsiDot                     = 17
    Psi                             = [indexPsi, indexPsiDot]

    state                           = range(0,18)

    # ---------- INSTANTIATE SOLVER ----------------------------------------------------------------
    totalDwell                      = sum([i for i in dwellTime])
    
    
    numOfBooleanVars                = numOfPartitions * horizon * (totalDwell * 4) # 4 states
    numberOfStates                  = 3 * 4 + 2         # x-y-z are 4th order integrators, psi 2nd order integrator
    numberOfStates                  = numberOfStates * totalDwell
    
    
    numberOfInputs                  = 4                 # one input for each integrator chain
    numberOfInputs                  = numberOfInputs  * totalDwell

    numOfRealVars                   = numberOfStates + numberOfInputs
    numOfConvIFClauses              = (numOfPartitions * horizon) + (horizon * (3 * 4 + 2))
                    # we add dynamics as interface clauses to be slacked later on. We have 3 integrators with 4 states and one integrator with w states.
    
    cores                           = 1
    
    tolerance                       = 1E-3
    solver                          = SMConvexSolver(numOfBooleanVars, numOfRealVars, numOfConvIFClauses,
                                                     maxNumberOfIterations = 10000,
                                                     counterExampleStrategy = strategy,
                                                     verbose = 'OFF',
                                                     profiling = 'false',
                                                     numberOfCores = cores,
                                                     slackTolerance = tolerance)
                                                     
                                                     
    # ---------- ADD ADJACENY MATRIX ----------------------------------------------------------------
    # for the $k$th time step , generate clause ensuring that only one partition is active
    # P_{1,k} + P_{2,k} + P_{p,k} = 1
    # where $p$ = number of partitions
    for horizonCounter in range(0, horizon):
        indexShift                   = horizonCounter*numOfPartitions
        solver.addBoolConstraint(sum([BoolVar2Int(solver.bVars[i+indexShift]) for i in range (0, numOfPartitions)]) == 1)



    # unroll the state machine representing the workspace adjacency
    # For the $i$ th partition, generate the follwoing clause
    # P_{i,k} => P_{adj1,k+1} OR P_{adj2,k+1} OR ... P_{adjn,k+1}
    for horizonCounter in range(0, horizon-1):
        indexShift                  = horizonCounter*numOfPartitions
        indexShift_plus             = (horizonCounter+1)*numOfPartitions
        for counter in range(0,numOfPartitions):
            adjacent                = partitions[counter]['adjacent']
            anticedent              = [solver.bVars[i + indexShift_plus] for i in adjacent]
            #print 'anticedent of', counter+indexShift, 'is ', anticedent
            solver.addBoolConstraint(
                IMPLIES(    solver.bVars[counter+indexShift], OR(*anticedent)
                            # the * used to unpack the python list to arguments
                )
            )

    # ---------- ADD PREVIOUS COUNTER EXAMPLES ----------------------------------------------------
    for counterExample in counterExamples:
        #print 'ADD CE', counterExample
        if not counterExample:
            continue
        constraint                  = [ NOT(solver.convIFClauses[counter]) for counter in counterExample ]
        solver.addBoolConstraint(OR(*constraint))

    # ---------- SPECIFICATIONS -------------------------------------------------------------------

    solver.addBoolConstraint(solver.bVars[startRegion])
    
    goalIndecies                    = [i*numOfPartitions + endRegion for i in range(0,horizon)]
    goals                           = [solver.bVars[i] for i in goalIndecies]
    solver.addBoolConstraint(OR(*goals))

    # avoid obstacles regions
    
    badIndecies                     = [i*numOfPartitions + j for i in range(0,horizon-1) for j in obstaclPartitions]
    for badIndex in badIndecies:
        solver.addBoolConstraint(NOT(solver.bVars[badIndex]))


    # ---------- ADD INPUT AND DYNAMICS CONSTRAINTS ------------------------------------------------
    #input constraints (can be added here if needed)



    # discrete dynamics of 4th order chain of integrators
    A4_cont                         = np.array([[0,1,0,0],[0,0,1,0], [0,0,0,1], [0,0,0,0] ])       #integrators
    A4_disc                         = expm(A4_cont*Ts)

    B4_cont                         = np.array([[0,0,0,1]])
    B4_disc                         = np.array([[0,0,0,Ts]])

    I4                              = np.array([[1,0,0,0], [0,1,0,0],[0,0,1,0], [0,0,0,1] ])       #identity

    #print "A4_cont", A4_cont
    #print "A4_disc", A4_disc
    A4_LP                           = np.concatenate((I4, -1*A4_disc, -1*B4_disc.T), axis = 1)
    b4_LP                           = [0]*4
    
    # discrete dynamics of 2nd order chain of integrators
    A2_cont                         = np.array([[0,1],[0,0] ])       #integrators
    A2_disc                         = expm(A2_cont*Ts)

    B2_cont                         = np.array([[0,1]])
    B2_disc                         = np.array([[0,Ts]])

    I2                              = np.array([[1,0], [0,1] ])       #identity

    A2_LP                           = np.concatenate((I2, -1*A2_disc, -1*B2_disc.T), axis = 1)
    b2_LP                           = [0]*2

    # initial state constraints
    vars                            = [solver.rVars[i] for i in [indexX, indexXDot, indexXDDot]]
    initialStateConstraint          = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), X_0, vars,sense = "E")
    solver.addConvConstraint(initialStateConstraint)

    vars                            = [solver.rVars[i] for i in [indexY, indexYDot, indexYDDot]]
    initialStateConstraint          = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), Y_0, vars,sense = "E")
    solver.addConvConstraint(initialStateConstraint)

    vars                            = [solver.rVars[i] for i in [indexZ, indexZDot, indexZDDot]]
    initialStateConstraint          = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), Z_0, vars,sense = "E")
    solver.addConvConstraint(initialStateConstraint)

    vars                            = [solver.rVars[i] for i in Psi]
    initialStateConstraint          = LPClause(np.array([[1,0], [0,1] ]), P_0, vars,sense = "E")
    solver.addConvConstraint(initialStateConstraint)
    

    # TBD, works for Dwell == 1 only <------- FIX LATER
    slackIFdynamicsOffset            = horizon * numOfPartitions
    #for horizonCounter in range(0, (horizon * totalDwell) - 1):
    for horizonCounter in range(0, totalDwell - 1):
        slackIFdynamicsShift        = slackIFdynamicsOffset + (3 * 4 + 2) * horizonCounter
        # Dynamics Constraints, x-axis
        indexShift                  = (horizonCounter)* len(state)
        indexShiftplus              = (horizonCounter+1)* len(state)
        
        xVars                       = [solver.rVars[i + indexShift] for i in X]
        xplusVars                   = [solver.rVars[i + indexShiftplus] for i in X]
        vars                        = xplusVars + xVars + [solver.rVars[indexInputX + indexShift]]
        #print "vars = ", vars
        #dynamicConstraint           = LPClause(A4_LP, b4_LP, vars,sense = "E")
        #solver.addConvConstraint(dynamicConstraint)
        for stateCounter in range(0,4):
            dynamicConstraint           = LPClause(np.array([A4_LP[stateCounter,:]]), [b4_LP[stateCounter]], vars,sense = "E")
            solver.addBoolConstraint(solver.convIFClauses[slackIFdynamicsShift + stateCounter])
            solver.setConvIFClause(dynamicConstraint, slackIFdynamicsShift + stateCounter)
        
        # Dynamics Constraints, y-axis
        yVars                       = [solver.rVars[i + indexShift] for i in Y]
        yplusVars                   = [solver.rVars[i + indexShiftplus] for i in Y]
        vars                        = yplusVars + yVars + [solver.rVars[indexInputY + indexShift]]
        
        
        #dynamicConstraint = LPClause(A4_LP, b4_LP, vars,sense = "E")
        #solver.addConvConstraint(dynamicConstraint)
        for stateCounter in range(0,4):
            dynamicConstraint           = LPClause(np.array([A4_LP[stateCounter,:]]), [b4_LP[stateCounter]], vars,sense = "E")
            solver.addBoolConstraint(solver.convIFClauses[slackIFdynamicsShift + 4 + stateCounter])
            solver.setConvIFClause(dynamicConstraint, slackIFdynamicsShift + 4 + stateCounter)


        # Dynamics Constraints, z-axis
        zVars                       = [solver.rVars[i + indexShift] for i in Z]
        zplusVars                   = [solver.rVars[i + indexShiftplus] for i in Z]
        vars                        = zplusVars + zVars + [solver.rVars[indexInputZ + indexShift]]
        
        #dynamicConstraint = LPClause(A4_LP, b4_LP, vars,sense = "E")
        #solver.addConvConstraint(dynamicConstraint)
        for stateCounter in range(0,4):
            dynamicConstraint           = LPClause(np.array([A4_LP[stateCounter,:]]), [b4_LP[stateCounter]], vars,sense = "E")
            solver.addBoolConstraint(solver.convIFClauses[slackIFdynamicsShift + 8 + stateCounter])
            solver.setConvIFClause(dynamicConstraint, slackIFdynamicsShift + 8 + stateCounter)
        

        # Dynamics Constraints, psi-axis
        pVars                       = [solver.rVars[i + indexShift] for i in Psi]
        pplusVars                   = [solver.rVars[i + indexShiftplus] for i in Psi]
        vars                        = pplusVars + pVars + [solver.rVars[indexInputPsi + indexShift]]

        #dynamicConstraint = LPClause(A2_LP, b2_LP, vars,sense = "E")
        #solver.addConvConstraint(dynamicConstraint)
        for stateCounter in range(0,2):
            dynamicConstraint           = LPClause(np.array([A2_LP[stateCounter,:]]), [b2_LP[stateCounter]], vars,sense = "E")
            solver.addBoolConstraint(solver.convIFClauses[slackIFdynamicsShift + 12 + stateCounter])
            solver.setConvIFClause(dynamicConstraint, slackIFdynamicsShift + 12 + stateCounter)


        # acceleration must be between -3 m/s^2 and 3 m/s^2
        
        solver.setUpperBound(indexXDDot + indexShift, 0.5 * accBound)
        solver.setLowerBound(indexXDDot + indexShift, -0.5 * accBound)
        
        
        
        solver.setUpperBound(indexYDDot + indexShift, 0.5 * accBound)
        solver.setLowerBound(indexYDDot + indexShift, -0.5 * accBound)
        
        
        solver.setUpperBound(indexZDDot + indexShift, 0.5 * accBound)
        solver.setLowerBound(indexZDDot + indexShift, -0.5 * accBound)
        
        
                
        #C4 smootheness constraint
        


        
        vars                        = [solver.rVars[i] for i in [indexXDDDot + indexShiftplus, indexXDDDot + indexShift]]
        smoothConstraint            = LPClause(np.array([[1.0, -1.0]]), [smoothBound], vars)
        solver.addConvConstraint(smoothConstraint)
        
        smoothConstraint            = LPClause(np.array([[-1.0, 1.0]]), [smoothBound], vars)
        solver.addConvConstraint(smoothConstraint)

        vars                        = [solver.rVars[i] for i in [indexYDDDot + indexShiftplus, indexYDDDot + indexShift]]

        smoothConstraint            = LPClause(np.array([[1.0, -1.0]]), [smoothBound], vars)
        solver.addConvConstraint(smoothConstraint)
        
        smoothConstraint            = LPClause(np.array([[-1.0, 1.0]]), [smoothBound], vars)
        solver.addConvConstraint(smoothConstraint)
        
        vars                        = [solver.rVars[i] for i in [indexZDDDot + indexShiftplus, indexZDDDot + indexShift]]
    
        smoothConstraint            = LPClause(np.array([[1.0, -1.0]]), [smoothBound], vars)
        solver.addConvConstraint(smoothConstraint)
        
        smoothConstraint            = LPClause(np.array([[-1.0, 1.0]]), [smoothBound], vars)
        solver.addConvConstraint(smoothConstraint)
        
        

    # velocity, acceleration, jerk should be equal to zero at the end.
    #indexShift                      = ((horizon * totalDwell) - 1)* len(state)
    indexShift                      = (totalDwell - 1)* len(state)
    

    vars                            = [solver.rVars[i] for i in [indexXDot + indexShift, indexXDDot + indexShift, indexXDDDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), XDot_f, vars,sense = "E")
    solver.addConvConstraint(finalStateConstraint)

    vars                            = [solver.rVars[i] for i in [indexYDot + indexShift, indexYDDot + indexShift, indexYDDDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), YDot_f, vars,sense = "E")
    solver.addConvConstraint(finalStateConstraint)

    vars                            = [solver.rVars[i] for i in [indexZDot + indexShift, indexZDDot + indexShift, indexZDDDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1,0,0], [0,1,0], [0,0,1] ]), ZDot_f, vars,sense = "E")
    solver.addConvConstraint(finalStateConstraint)

    vars                            = [solver.rVars[i] for i in [indexPsiDot + indexShift]]
    finalStateConstraint            = LPClause(np.array([[1] ]), PDot_f, vars,sense = "E")
    solver.addConvConstraint(finalStateConstraint)
    
    
    # if assigned, we add constraint on the final position
    if x_f is not None:
        vars                            = [solver.rVars[indexX + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [x_f], vars,sense = "E")
        solver.addConvConstraint(finalStateConstraint)
        
    if y_f is not None:
        vars                            = [solver.rVars[indexY + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [y_f], vars,sense = "E")
        solver.addConvConstraint(finalStateConstraint)
        

    if z_f is not None:
        vars                            = [solver.rVars[indexZ + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [z_f], vars,sense = "E")
        solver.addConvConstraint(finalStateConstraint)

    if p_f is not None:
        vars                            = [solver.rVars[indexPsi + indexShift]]
        finalStateConstraint            = LPClause(np.array([[1]]), [p_f], vars,sense = "E")
        solver.addConvConstraint(finalStateConstraint)

    
    # ---------- REGIONS CONSTRAINTS ----------------------------------------------------------------
    for horizonCounter in range(0, horizon):
        indexShift                  = horizonCounter * numOfPartitions
        for partitionCounter in range(0,numOfPartitions):
            solver.addBoolConstraint(IMPLIES(solver.bVars[partitionCounter + indexShift], solver.convIFClauses[partitionCounter + indexShift]))
            # skip initial condition and starts first step in the plan

            #stateIndexShift             = (horizonCounter * totalDwell * len(state))
            stateIndexShift             = sum([dwellTime[i] for i in range(0,horizonCounter+1)]) * len(state)
            
            vars                        = [solver.rVars[i + stateIndexShift] for i in [indexX, indexY, indexZ]]
            
            partition_A                 = partitions[partitionCounter]['A']
            partition_b                 = partitions[partitionCounter]['b']
            
            
            for dwellCounter in range(1, dwellTime[horizonCounter+1]):
                #stateIndexShift         = (horizonCounter * totalDwell * len(state))  + (dwellCounter * len(state) )
                stateIndexShift         = (sum([dwellTime[i] for i in range(0,horizonCounter+1)]) * len(state))  + (dwellCounter * len(state) )
                
                partition_A_new_rows    = np.concatenate(
                                            (np.zeros((   np.shape(partitions[partitionCounter]['A'])[0], np.shape(partition_A)[1] )),
                                            partitions[partitionCounter]['A']),
                                            axis = 1)
                
                partition_A             = np.concatenate((partition_A,
                                                np.zeros((np.shape(partition_A)[0], np.shape(partitions[partitionCounter]['A'])[1]))),
                                                axis = 1)
                partition_A             = np.concatenate((partition_A, partition_A_new_rows))
                partition_b             = np.concatenate((partition_b,partitions[partitionCounter]['b']))
                
                vars                    = vars + [solver.rVars[i + stateIndexShift] for i in [indexX, indexY, indexZ]]
            
            
            regionconstraint            = LPClause(partition_A, partition_b, vars)
            solver.setConvIFClause(regionconstraint, partitionCounter + indexShift)

    #print solver.ConvSolver.linear_constraints.get_senses()
    #print solver.ConvSolver.linear_constraints.get_rhs()
    #print solver.ConvSolver.linear_constraints.get_rows()

    start                           = timeit.default_timer()
    stateTraj, regionTraj           = solver.solve()
    end                             = timeit.default_timer()
    time_smt                        = end - start

    counter_examples                = solver.counterExamples
    number_counter_examples         = len(counter_examples)
    #print '#CE = ', number_counter_examples
    #print 'Time = ', time_smt

    #print 'counterExamples', counter_examples
    if not stateTraj:
        traj                            = { 'inputXTraj': [],
                                        'inputYTraj': [],
                                        'inputZTraj': [],
                                        'inputPTraj': [],
                                        'xTraj'     : [],
                                        'xDotTraj'  : [],
                                        'xDDotTraj' : [],
                                        'xDDDotTraj': [],
                                        'yTraj'     : [],
                                        'yDotTraj'  : [],
                                        'yDDotTraj' : [],
                                        'yDDDotTraj': [],
                                        'zTraj'     : [],
                                        'zDotTraj'  : [],
                                        'zDDotTraj' : [],
                                        'zDDDotTraj': [],
                                        'pTraj'     : [],
                                        'pDotTraj'  : [],
                                        'seqTraj'   : [],
                                        'counterExamples' : counter_examples
                                        }
        return  traj
    
    
    #inputXTraj                      = [stateTraj[i + indexInputX]   for i in  range(0, len(state) * horizon * totalDwell, len(state))]
    inputXTraj                      = [stateTraj[i + indexInputX]   for i in  range(0, len(state) * totalDwell, len(state))]
    inputYTraj                      = [stateTraj[i + indexInputY]   for i in  range(0, len(state) * totalDwell, len(state))]
    inputZTraj                      = [stateTraj[i + indexInputZ]   for i in  range(0, len(state) * totalDwell, len(state))]
    inputPTraj                      = [stateTraj[i + indexInputPsi] for i in  range(0, len(state) * totalDwell, len(state))]

    xTraj                           = [stateTraj[i + indexX]        for i in  range(0, len(state) * totalDwell, len(state))]
    xDotTraj                        = [stateTraj[i + indexXDot]     for i in  range(0, len(state) * totalDwell, len(state))]
    xDDotTraj                       = [stateTraj[i + indexXDDot]    for i in  range(0, len(state) * totalDwell, len(state))]
    xDDDotTraj                      = [stateTraj[i + indexXDDDot]   for i in  range(0, len(state) * totalDwell, len(state))]

    yTraj                           = [stateTraj[i + indexY]        for i in  range(0, len(state) * totalDwell, len(state))]
    yDotTraj                        = [stateTraj[i + indexYDot]     for i in  range(0, len(state) * totalDwell, len(state))]
    yDDotTraj                       = [stateTraj[i + indexYDDot]    for i in  range(0, len(state) * totalDwell, len(state))]
    yDDDotTraj                      = [stateTraj[i + indexYDDDot]   for i in  range(0, len(state) * totalDwell, len(state))]

    zTraj                           = [stateTraj[i + indexZ]        for i in  range(0, len(state) * totalDwell, len(state))]
    zDotTraj                        = [stateTraj[i + indexZDot]     for i in  range(0, len(state) * totalDwell, len(state))]
    zDDotTraj                       = [stateTraj[i + indexZDDot]    for i in  range(0, len(state) * totalDwell, len(state))]
    zDDDotTraj                      = [stateTraj[i + indexZDDDot]   for i in  range(0, len(state) * totalDwell, len(state))]


    pTraj                           = [stateTraj[i + indexPsi]      for i in  range(0, len(state) * totalDwell, len(state))]
    pDotTraj                        = [stateTraj[i + indexPsiDot]   for i in  range(0, len(state) * totalDwell, len(state))]


    regionSeq_augmented             = regionTraj
    seqTraj                         = [y for i in range(0,horizon+1) for y in [regionSeq_augmented[i]] * dwellTime[i]]

    traj                            = { 'inputXTraj': inputXTraj,
                                        'inputYTraj': inputYTraj,
                                        'inputZTraj': inputZTraj,
                                        'inputPTraj': inputPTraj,
                                        'xTraj'     : xTraj,
                                        'xDotTraj'  : xDotTraj,
                                        'xDDotTraj' : xDDotTraj,
                                        'xDDDotTraj': xDDDotTraj,
                                        'yTraj'     : yTraj,
                                        'yDotTraj'  : yDotTraj,
                                        'yDDotTraj' : yDDotTraj,
                                        'yDDDotTraj': yDDDotTraj,
                                        'zTraj'     : zTraj,
                                        'zDotTraj'  : zDotTraj,
                                        'zDDotTraj' : zDDotTraj,
                                        'zDDDotTraj': zDDDotTraj,
                                        'pTraj'     : pTraj,
                                        'pDotTraj'  : pDotTraj,
                                        'seqTraj'   : seqTraj,
                                        'counterExamples' : counter_examples
                                        }


    '''

    print 'inputXtraj =',   inputXTraj
    print 'inputYtraj =',   inputYTraj
    print 'inputZtraj =',   inputZTraj
    print ''


    print 'xTraj =',        xTraj
    print 'xDotTraj =',     xDotTraj
    print 'xDDotTraj =',    xDDotTraj
    print 'xDDDotTraj =',   xDDDotTraj
    print ''


    print 'yTraj =',        yTraj
    print 'yDotTraj =',     yDotTraj
    print 'yDDotTraj =',    yDDotTraj
    print 'yDDDotTraj =',   yDDDotTraj
    print ''


    print 'zTraj =',        zTraj
    print 'zDotTraj =',     zDotTraj
    print 'zDDotTraj =',    zDDotTraj
    print 'zDDDotTraj =',   zDDDotTraj


    print 'pTraj =',        pTraj
    print 'pDotTraj =',     pDotTraj


    print 'seqTraj = ',     seqTraj

    
    '''


    #solver.setConvIFClause(QPClause(Q_i, c_i, b_i, x_var), counter)
    #solver.addInterfaceConstraint(IMPLIES(NOT(solver.bVars[counter]), solver.convIFClauses[counter]))


    #solver.addInterfaceConstraint(sum([BoolVar2Int(solver.bVars[i]) for i in range (0, p)]) <= max_sensors_under_attack)

    return traj



#***************************************************************************************************
#***************************************************************************************************
#
#         Robotic LTL Motion Planning Test
#
#***************************************************************************************************
#***************************************************************************************************
def robotic_LTL_workspace(xMaxTic, yMaxTick, zMaxTick):
    ticksX                      = range(0,xMaxTic)       #[-6.0, -5.0, -2.0, 2.0, 5.0, 6.0]
    ticksY                      = range(0,yMaxTick)        #[-1.0, -0.5, 0.5, 1.0]
    ticksZ                      = range(0,zMaxTick)        #[0.0,   0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]

    wkspDimensions              = [ticksX[0], ticksX[-1], ticksY[0], ticksY[-1], ticksZ[0], ticksZ[-1]]
    obstaclesXTag               = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
    obstaclesZTag               = [1,3]
    
    partitions                  = list()
    obstaclPartitions           = list()
    
    tag2PartitionCounter        = list()
    partitionCounter            = 0
    for counterX in range(1, len(ticksX)):
        for counterY in range(1, len(ticksY)):
            for counterZ in range(1, len(ticksZ)):
                tag             = ''.join(str(i) for i in [counterX, counterY, counterZ])
                tag2PartitionCounter.append((tag, partitionCounter))
                partitionCounter= partitionCounter + 1

    tag2PartitionCounterDic     = dict(tag2PartitionCounter)
    
    for counterX in range(1, len(ticksX)):
        for counterY in range(1, len(ticksY)):
            for counterZ in range(1, len(ticksZ)):
                A = np.array([  [-1.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [0.0, -1.0, 0.0],
                                [0.0, 1.0, 0.0],
                                [0.0, 0.0, -1.0],
                                [0.0, 0.0, 1.0]])
                b = np.array ([ -1 * ticksX[counterX - 1], ticksX[counterX],
                                -1 * ticksY[counterY - 1], ticksY[counterY],
                                -1 * ticksZ[counterZ - 1], ticksZ[counterZ]
                                ])
                '''
                point1          = [ticksX[counterX - 1], ticksX[counterY - 1], ticksX[counterZ - 1]]
                point2          = [ticksX[counterX], ticksX[counterY - 1], ticksX[counterZ - 1]]
                point3          = [ticksX[counterX], ticksX[counterY], ticksX[counterZ - 1]]
                point4          = [ticksX[counterX - 1], ticksX[counterY], ticksX[counterZ - 1]]
                point5          = [ticksX[counterX - 1], ticksX[counterY - 1], ticksX[counterZ]]
                point6          = [ticksX[counterX], ticksX[counterY - 1], ticksX[counterZ]]
                point7          = [ticksX[counterX], ticksX[counterY], ticksX[counterZ]]
                point8          = [ticksX[counterX - 1], ticksX[counterY], ticksX[counterZ]]
                points          = [point1,point2,point3,point4,point5,point6,point7,point8]
                polyhedra       = Vrep (points)
                '''
                tag             = ''.join(str(i) for i in [counterX, counterY, counterZ])
                
                
                isObstacle      = 0
                if counterX in obstaclesXTag and counterZ in obstaclesZTag:
                    isObstacle  = 1
                    obstaclPartitions.append(tag2PartitionCounterDic[tag])
                
                adjacents       = list()
                adjacentTag     = ''.join(str(i) for i in [counterX, counterY, counterZ])        #self is always adjacent
                adjacents.append(tag2PartitionCounterDic[adjacentTag])
                
                if counterX-1 >= 1:
                    adjacentTag = ''.join(str(i) for i in [counterX-1, counterY, counterZ])
                    adjacents.append(tag2PartitionCounterDic[adjacentTag])
                
                if counterX+1 <= len(ticksX)-1:
                    adjacentTag = ''.join(str(i) for i in [counterX+1, counterY, counterZ])
                    adjacents.append(tag2PartitionCounterDic[adjacentTag])
                
                if counterY-1 >= 1:
                    adjacentTag = ''.join(str(i) for i in [counterX, counterY-1, counterZ])
                    adjacents.append(tag2PartitionCounterDic[adjacentTag])
                    
                if counterY+1 <= len(ticksY)-1:
                    adjacentTag = ''.join(str(i) for i in [counterX, counterY+1, counterZ])
                    adjacents.append(tag2PartitionCounterDic[adjacentTag])
                
                if counterZ-1 >= 1:
                    adjacentTag = ''.join(str(i) for i in [counterX, counterY, counterZ-1])
                    adjacents.append(tag2PartitionCounterDic[adjacentTag])

                if counterZ+1 <= len(ticksZ)-1:
                    adjacentTag = ''.join(str(i) for i in [counterX, counterY, counterZ+1])
                    adjacents.append(tag2PartitionCounterDic[adjacentTag])

                partition       = {'A':A, 'b':b, 'adjacent': adjacents, 'tag':tag, 'isObstacle':isObstacle}
                partitions.append(partition)
                


    return partitions, obstaclPartitions, tag2PartitionCounterDic, wkspDimensions


if __name__ == "__main__":
    np.random.seed(0)
    #secureStateEstimation_main()
    #scalability_main()
    robotic_LTL_main()


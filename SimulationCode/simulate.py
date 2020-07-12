#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Wright-Fisher simulations with mutation for two-locus system.
	@author: Saumil Shah
	@email: saumil.shah@students.iiserpune.ac.in '''

''' Genotypes are of the form 'ab', where 'a' is the allele at locus A and 'b' 
is the allele at locus B. There are three possible alleles at both loci that 
correspond to low expression, intermediate expression, and high expression. We 
will denote these alleles as L, O, and H. Hence the set of all possible 
genotypes is {LL, LO, LH, OL, OO, OH, HL, HO, HH}. An allele at a locus may 
change due to mutations. '''

import gc
import os
import time
import random
import numpy as np
import pandas as pd
import argparse as ap
from sys import stdout

# take simulation parameters as command line arguments
PARSER = ap.ArgumentParser()
PARSER.add_argument('POPULATION_SIZE')
PARSER.add_argument('TOTAL_GENERATIONS')
PARSER.add_argument('TOTAL_REPLICATES')
PARSER.add_argument('BLOCK_ID')
ARGUMENTS = PARSER.parse_args()
POPULATION_SIZE = int(ARGUMENTS.POPULATION_SIZE)
TOTAL_GENERATIONS = int(ARGUMENTS.TOTAL_GENERATIONS)
TOTAL_REPLICATES = int(ARGUMENTS.TOTAL_REPLICATES)
BLOCK_ID = int(ARGUMENTS.BLOCK_ID)

# define mutation rates
c1 = 1 * 1e-3 # very common
c2 = 5 * 1e-4 # common
r1 = 5 * 1e-5 # rare
r2 = 5 * 1e-8 # very rare
r3 = 1 * 1e-9 # very very rare

# mutation rates at locus A
A_MU_RATES = np.array([[1-r1-r2, r1, r2],
                 [c1, 1-c1-r1, r1],
                 [c1, c1, 1-c1-c1]])

# mutation rates at locus B
B_MU_RATES = np.array([[1-r2-r3, r2, r3],
                 [c2, 1-c2-r2, r2],
                 [c2, c2, 1-c1-c2]])

''' We use following digits for computations
0 - locus A
1 - locus B

0 - low
1 - intermediate 
2 - high

Hence, the genotypes are

[ 00 01 02       [ LL LO LH
  10 11 12    =    OL OO OH
  20 21 22 ]       HL HO HH ] '''

GENOTYPES = np.array([0, 1, 2, 10, 11, 12, 20, 21, 22], dtype = np.int8)

# fitness landscape
FITNESS = np.array([[1, 2, 9],
                 [3, 4, 8],
                 [5, 6, 7]], dtype = np.int8)

# change of alleles in digit basis
CHANGE = np.array([[0, 1, 2],
                 [-1, 0, 1],
                 [-2, -1, 0]], dtype = np.int8)

# changes on locus A follow 10*CHANGE
# (allele a is in second pace from right in digit basis)
A_CHANGE = 10*np.copy(CHANGE)

# changes on locus B follow CHANGE
B_CHANGE = np.copy(CHANGE)

# converts genotype 'ab' into list of alleles (a, b)
# can also process a list of genotypes as an input, i.e. a population 'p'
def index(ab):
    return np.divmod(ab, 10)

# returns count of an allele 'a' on a locus 'l' for a population 'p'
def freq_allele(p, l, a):
    return np.flatnonzero(index(p)[l] == a)

# returns counts of all alleles on a locus 'l' for a population 'p' 
def freq_all(p, l):
    return np.array([freq_allele(p, l, 0).size, 
        freq_allele(p, l, 1).size, 
        freq_allele(p, l, 2).size]).reshape((3,1))

# chooses n individuals out of p randomly without replacement
def pick(p, n):
    return p[ random.sample(range(p.size), min(n,p.size)) ]

# returns counts of each genotype in the population 'p'
def count(p):
    return np.bincount(p, minlength = 23)[GENOTYPES]

# start time benchmark counter
T1 = time.perf_counter()

# initiate an empty array
df = np.zeros((11))

for r in range(TOTAL_REPLICATES):
    replicate_id = 1000*BLOCK_ID + r

    # initiate wild type population
    population = 11*np.ones(POPULATION_SIZE, dtype = np.int8)
    df = np.vstack((df, np.append([replicate_id, 0],count(population))))

    for g in range(TOTAL_GENERATIONS):
        gc.collect() # clear unused memory

        # generate offspring population
        offsprings = np.repeat(population, FITNESS[index(population)])

        # find mutation events
        nmoa, nmob = np.zeros(A_MU_RATES.shape), np.zeros(B_MU_RATES.shape)
        nmoa = np.random.poisson(freq_all(offsprings, 0)*A_MU_RATES)
        nmob = np.random.poisson(freq_all(offsprings, 1)*B_MU_RATES)
        np.fill_diagonal(nmoa,0)
        np.fill_diagonal(nmob,0)
        nmoa, nmob = np.int64([nmoa, nmob])
        
        # mutate offspring individuals
        for ma in np.flatnonzero(nmoa):
            f1,t1 = np.divmod(ma, 3)
            offsprings[ pick(freq_allele(offsprings, 0, f1), nmoa[f1, t1]) ] += A_CHANGE[f1, t1]
            
        for mb in np.flatnonzero(nmob):
            f2,t2 = np.divmod(mb, 3)
            offsprings[ pick(freq_allele(offsprings, 1, f2), nmob[f2, t2]) ] += B_CHANGE[f2, t2]
        
        # Wright-Fisher step
        population = pick(offsprings, POPULATION_SIZE)
        df = np.vstack((df, np.append([replicate_id, g+1],count(population))))
        
        # update progressbar
        pbar = int( (g+1)*100 / TOTAL_GENERATIONS )
        stdout.write('\r [{:-<50}] gen:{:d}% rep:{:d}%\t'.format('#'*int(pbar/2), pbar, int((r)*100/TOTAL_REPLICATES)))
        stdout.flush()

# stop time benchmark counter
T2 = time.perf_counter()

# dataframe column names
names=['rep_id','g','LL', 'LO', 'LH', 'OL', 'OO', 'OH', 'HL', 'HO', 'HH']
df = pd.DataFrame(df,columns=names,dtype=np.int).drop(0,0)

# file name containing population size and generations
fn = 'p{:.0e}_g{}.txt'.format(POPULATION_SIZE, TOTAL_GENERATIONS)

if not os.path.isfile(fn): # datafile doesnt exist
   df.to_csv(fn, '\t', index=False)
else: # datafile exists so append without the header
   df.to_csv(fn, '\t', index=False, mode='a', header=False)
   
stdout.write('\n Done! Time elapsed: {:.1f} Hours'.format((T2-T1) / 3600))
stdout.flush()
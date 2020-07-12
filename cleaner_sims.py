#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Wright-Fisher simulations 
    for qu_switch with mutations.
    @author: Saumil Shah'''
    
import gc
import random
import numpy as np
import pandas as pd
import argparse as ap
from time import time
from sys import stdout
from datetime import datetime
from matplotlib import pyplot as plt

#maximum individuals to simulate
prsr = ap.ArgumentParser()
prsr.add_argument('N')
prsr.add_argument('G')
prsr.add_argument('R')
args = prsr.parse_args()
ccap = int(args.N)
mgen = int(args.G)
mrep = int(args.R)

#mutation rates per generation per individual per locus
c1 = 1 * 1e-3 #very common
c2 = 5 * 1e-4 #common
r1 = 5 * 1e-5 #rare
r2 = 5 * 1e-8 #very rare
r3 = 1 * 1e-9 #very very rare

#mutation rates at locus A
amut = np.array([[1-r1-r2, r1, r2],
                 [c1, 1-c1-r1, r1],
                 [c1, c1, 1-c1-c1]])

#mutation rates at locus B
bmut = np.array([[1-r2-r3, r2, r3],
                 [c2, 1-c2-r2, r2],
                 [c2, c2, 1-c1-c2]])

#genotypes
gref = np.array([0,1,2,10,11,12,20,21,22], dtype=np.int8)

#fitness landscape
fitl = np.array([[1, 2, 9],
                 [3, 4, 8],
                 [5, 6, 7]], dtype=np.int8)

#changes on locus A
aswp = np.array([[0, 10, 20],
                 [-10, 0, 10],
                 [-20, -10, 0]], dtype=np.int8)

#changes on locus B
bswp = np.array([[0, 1, 2],
                 [-1, 0, 1],
                 [-2, -1, 0]], dtype=np.int8)

#converts genotype to index
def to_index(genotype):
    return np.divmod(genotype,10)

#returns counts of allele a on locus A
def find_index(p, locus, allele):
    return np.flatnonzero(np.divmod(p,10)[locus]==allele)

def count_freq(p):
    return np.bincount(p, minlength=23)[dn.gref]

#updates population with offsprings
def reproduce(p,fitness_landscape):
    return np.repeat(p, fitness_landscape[to_index(p)])

def 
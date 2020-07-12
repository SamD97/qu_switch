#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Basic parameters for simulations.
    @author: Saumil Shah'''

import numpy as np

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

#mutation rates at locus A
cmut = np.array([[1-r1-r2, r1, r2],
                 [c1, 1-c1-r1, r1],
                 [c1, c1, 1-c1-c1]])

#genotypes
gref = np.array([  0,   1,   2,
                  10,  11,  12,
                  20,  21,  22,
                 100, 101, 102,
                 110, 111, 112,
                 120, 121, 122,
                 200, 201, 202,
                 210, 211, 212,
                 220, 221, 222], dtype=np.int16)

#fitness landscape
fitl = np.array([[[ 2,  3, 10],
                  [ 5,  6, 10],
                  [14, 13, 11]],

                 [[ 4,  5, 12],
                  [ 7,  8, 12],
                  [13, 13, 12]],
        
                 [[ 6,  7, 13],
                  [ 9, 10, 14],
                  [12, 13, 13]]], dtype=np.int16)

#changes on locus A
aswp = np.array([[   0,  100,  200],
                 [-100,    0,  100],
                 [-200, -100,    0]], dtype=np.int16)

#changes on locus B
bswp = np.array([[  0,  10,  20],
                 [-10,   0,  10],
                 [-20, -10,   0]], dtype=np.int16)

#changes on locus A
cswp = np.array([[0, 1, 2],
                 [-1, 0, 1],
                 [-2, -1, 0]], dtype=np.int16)
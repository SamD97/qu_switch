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
bmut = np.array([[1-r1-r2, r1, r2],
                 [c1, 1-c1-r1, r1],
                 [c1, c1, 1-c1-c1]])

#genotypes
gref = np.array([0,1,2,10,11,12,20,21,22], dtype=np.int8)

#fitness landscape
fitl = np.array([[1, 2, 7],
                 [3, 4, 8],
                 [5, 6, 9]], dtype=np.int8)

#changes on locus A
aswp = np.array([[0, 10, 20],
                 [-10, 0, 10],
                 [-20, -10, 0]], dtype=np.int8)

#changes on locus B
bswp = np.array([[0, 1, 2],
                 [-1, 0, 1],
                 [-2, -1, 0]], dtype=np.int8)
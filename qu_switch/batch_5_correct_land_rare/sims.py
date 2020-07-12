#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Wright-Fisher simulations 
    for qu_switch with mutations.
    @author: Saumil Shah'''


'''
Genotypes looks like 'ab';    a,b = 0,1,2
 where position 1 ~ A allele,
  and position 2 ~ B allele,

0 -  low expression
1 -  mid expression
2 - high expression

possible genotypes

[ 00 01 02       [ LL LO LH
  10 11 12    ~    OL OO OH
  20 21 22 ]       HL HO HH ]

'''

import gc
import random
from time import time
from datetime import datetime
import numpy as np
import argparse as ap
from sys import stdout
from matplotlib import pyplot as plt

#maximum individuals to simulate
prsr = ap.ArgumentParser()
prsr.add_argument('ccap')
prsr.add_argument('mgen')
prsr.add_argument('mrep')
args = prsr.parse_args()
ccap = int(args.ccap)
mgen = int(args.mgen)
mrep = int(args.mrep)

import defn as dn

#returns indices of given genotype
def indx(x):
    return np.divmod(x,10)

#returns counts of allele a on locus A
def fnda(p, a):
    return np.flatnonzero(np.divmod(p,10)[0]==a)

#returns counts of allele b on locus B
def fndb(p, b):
    return np.flatnonzero(np.divmod(p,10)[1]==b)

#returns n random indices of p without replacements
def pick(p, n):
    return np.random.choice(p,min(n,p.size), replace=False)

#returns counts of individuals with given allele on locus A
def frea(p):
    return np.array([fnda(p,0).size, fnda(p,1).size, fnda(p,2).size]).reshape((3,1))

#returns counts of individuals with given allele on locus B
def freb(p):
    return np.array([fndb(p,0).size, fndb(p,1).size, fndb(p,2).size]).reshape((3,1))

t1 = time()

df=[]

for r in range(mrep):
    pop0 = 11*np.ones(ccap, dtype=np.int8)
    popi = np.copy(pop0)

    for g in range(mgen):
        gc.collect()
        popo = np.repeat(popi,dn.fitl[indx(popi)])
        nmoa,nmob = np.zeros(dn.amut.shape),np.zeros(dn.bmut.shape)
        nmoa = np.random.poisson(frea(popo)*dn.amut)
        nmob = np.random.poisson(freb(popo)*dn.bmut)
        np.fill_diagonal(nmoa,0)
        np.fill_diagonal(nmob,0)
        nmoa, nmob = np.int64([nmoa,nmob])
        
        for ma in np.flatnonzero(nmoa):
            f1,t1 = np.divmod(ma,3)
            popo[ pick(fnda(popo,f1), nmoa[f1,t1]) ] += dn.aswp[f1,t1]
            
        for mb in np.flatnonzero(nmob):
            f2,t2 = np.divmod(mb,3)
            popo[ pick(fndb(popo,f2), nmob[f2,t2]) ] += dn.bswp[f2,t2]
        
        popi = pick(popo,ccap)
        genc = np.bincount(popi, minlength=23)[dn.gref]
        
        pbar = int( (g+1)*100 / mgen )
        stdout.write('\r [{:-<50}] gen:{:d}% rep:{:d}%\t'.format('#'*int(pbar/2), pbar, int((r)*100/mrep)))
        stdout.flush()

    df += [genc]

t2 = time()
nw = datetime.now()
print('\n {:.1f} secs'.format(t2-t1))
np.savetxt('{}_{}_{}_{:.1e}.txt'.format(nw.strftime('%y%m%d_%H%M%S'), mgen, mrep, ccap),df,'%10d','\t')
frem = np.average(df, axis=0).round(2)
fred = np.    std(df, axis=0).round(2)
plt.figure()
plt.xticks(np.arange(9),['LL', 'LO', 'LH', 'OL', 'OO', 'OH', 'HL', 'HO', 'HH'])
plt.ylim((0,1))
plt.bar(np.arange(9), frem/ccap, yerr=fred/ccap)
plt.xlabel('genotypes')
plt.ylabel('fraction')
plt.savefig('plots/{}_{}_{}_{:.1e}.svg'.format(nw.strftime('%y%m%d_%H%M%S'), mgen, mrep, ccap))
plt.close()
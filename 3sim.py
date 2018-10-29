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
import os
import time
import random
import numpy as np
import pandas as pd
import argparse as ap
from sys import stdout

#maximum individuals to simulate
prsr = ap.ArgumentParser()
prsr.add_argument('ccap')
prsr.add_argument('mgen')
prsr.add_argument('mrep')
args = prsr.parse_args()
ccap = int(args.ccap)
mgen = int(args.mgen)
mrep = int(args.mrep)

import def3 as dn

#returns indices of given genotype
def indx(x):
    a, d = np.divmod(x, 100)
    b, c = np.divmod(d, 10)
    return a,b,c

#returns counts of allele a on locus A
def fnda(p,l,a):
    return np.flatnonzero(indx(p)[l]==a)

#reproduction step skipping offsprings
def rprd(p,f,n):
    r = np.random.multinomial(n,f/sum(f))
    return np.repeat(p,r)

#returns n random indices of p without replacements
def pick(p, n):
    return p[ random.sample(range(p.size), min(n,p.size)) ]

#returns counts of individuals with given allele on given locus 
def frea(p,l):
    return np.array([fnda(p,l,0).size, fnda(p,l,1).size, fnda(p,l,2).size]).reshape((3,1))

def cunt(p):
    return np.bincount(p, minlength=223)[dn.gref]

t1 = time.time()

df=np.zeros((29))

for r in range(mrep):
    popi = 111*np.ones(ccap, dtype=np.int64)
    larn = np.random.randint(1000)
    df = np.vstack((df,np.append([larn,0],cunt(popi))))

    for g in range(mgen):
        gc.collect()
        popo = np.repeat(popi,dn.fitl[indx(popi)])
        nmoa,nmob = np.zeros(dn.amut.shape),np.zeros(dn.bmut.shape)
        nmoa = np.random.poisson(frea(popo,0)*dn.amut)
        nmob = np.random.poisson(frea(popo,1)*dn.bmut)        
        nmoc = np.random.poisson(frea(popo,2)*dn.cmut)
        np.fill_diagonal(nmoa,0)
        np.fill_diagonal(nmob,0)
        np.fill_diagonal(nmoc,0)
        nmoa, nmob, nmoc = np.int64([nmoa,nmob,nmoc])
        
        for ma in np.flatnonzero(nmoa):
            f1,t1 = np.divmod(ma,3)
            popo[ pick(fnda(popo,0,f1), nmoa[f1,t1]) ] += dn.aswp[f1,t1]
            
        for mb in np.flatnonzero(nmob):
            f2,t2 = np.divmod(mb,3)
            popo[ pick(fnda(popo,1,f2), nmob[f2,t2]) ] += dn.bswp[f2,t2]
            
        for mc in np.flatnonzero(nmoc):
            f3,t3 = np.divmod(mc,3)
            popo[ pick(fnda(popo,2,f3), nmob[f3,t3]) ] += dn.cswp[f3,t3]
        
        popi = pick(popo,ccap)
        genc = cunt(popi)
        df = np.vstack((df,np.append([larn,g+1],cunt(popi))))
        
        pbar = int( (g+1)*100 / mgen )
        stdout.write('\r [{:-<50}] gen:{:d}% rep:{:d}%\t'.format('#'*int(pbar/2), pbar, int((r)*100/mrep)))
        stdout.flush()

t2 = time.time()

names=['r','g',
       'LLL', 'LLO', 'LLH',
       'LOL', 'LOO', 'LOH',
       'LHL', 'LHO', 'LHH',
       'OLL', 'OLO', 'OLH',
       'OOL', 'OOO', 'OOH',
       'OHL', 'OHO', 'OHH',
       'HLL', 'HLO', 'HLH',
       'HOL', 'HOO', 'HOH',
       'HHL', 'HHO', 'HHH']

df = pd.DataFrame(df,columns=names,dtype=np.int).drop(0,0)
fn = '{}_{:.0e}.txt'.format(mgen, ccap)
if not os.path.isfile(fn):
   df.to_csv(fn, '\t', index=False)
else: # else it exists so append without writing the header
   df.to_csv(fn, '\t', index=False, mode='a', header=False)
   
stdout.write('\n Done! {:.1f}'.format((t2-t1)/3600))
stdout.flush()
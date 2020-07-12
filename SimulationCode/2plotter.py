#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Plotter for simulations.
    @author: Saumil Shah'''

import os
import gc
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from defn import fitl

path = input(' Please enter path to folder: ')
if path == '': path = os.getcwd()
os.chdir(path)

files=[]
for f in os.listdir(path):
        if f.endswith('pop'):
            files.append(f)
files.sort()
no_f = len(files)

print ()
print (' {}'.format('='*70))
print (' {:^70s}'.format('File-list'))
print (' {}'.format('='*70))
for f in files:
    print (' {:^70s}'.format(f))
print (' {}'.format('='*70))
print (' \t\tTotal {} .pop files were found.'.format(no_f))
print (' {}'.format('='*70))
print ()

names=['r','g','LL', 'LO', 'LH', 'OL', 'OO', 'OH', 'HL', 'HO', 'HH']

pops = np.zeros(no_f)
mfit = np.zeros(no_f)
sfit = np.zeros(no_f)

for i in range(no_f):
    gc.collect()
    
    mgen = int(files[i].split('_')[0])
    ccap = float(files[i].split('_')[1].split('.')[0])
    df = pd.read_csv(files[i], '\t')
    frem = df.groupby('g').mean().drop('r',1).iloc[mgen].values
    fred = df.groupby('g').std().drop('r',1).iloc[mgen].values / np.sqrt(ccap)
    pops[i] = ccap
    mfit[i] = np.average(fitl.flatten(), weights=frem)
    sfit[i] = np.sqrt(np.average( (fitl.flatten()-mfit[i])**2, weights=frem))
    
    plt.figure()
    plt.xticks(np.arange(9), names[2:], rotation=75)
    plt.ylim((0,1))
    plt.bar(np.arange(9), frem/ccap, yerr=fred/ccap)
    plt.xlabel('genotypes')
    plt.ylabel('fraction')
    plt.tight_layout()
    plt.savefig('plots/{}_{:.0e}.svg'.format(mgen, ccap))
    plt.close()
    
    sys.stdout.write('\r[{}{}] {}% {}\t'.format('#'*int(((i+1)*50)/no_f),
                     '-'*int(((no_f-i-1)*50)/no_f), int(((i+1)*100)/no_f), files[i][:-4]))
    sys.stdout.flush()

plt.figure()
plt.xscale('log')
plt.errorbar(pops,mfit,yerr=sfit,fmt='ko')
plt.xlabel('population size')
plt.ylabel('average fitness')
plt.savefig('plots/avg_fitness.svg')
plt.tight_layout()
plt.close()

sys.stdout.write('\n Done!')
sys.stdout.flush()
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

for i in range(no_f):
    gc.collect()
    
    mgen = int(files[i].split('_')[0])
    ccap = float(files[i].split('_')[1].split('.')[0])
    if not os.path.isdir('plots/{:.0e}'.format(ccap)):
        os.mkdir('plots/{:.0e}'.format(ccap))
    df = pd.read_csv(files[i], '\t')
    for rep in pd.unique(df.r):
        frem = df[df.r==rep].drop('r',1)
        mfit = fitl.flatten() * frem[frem.columns[1:]] / ccap
        sfit = np.sqrt( np.sum((fitl.flatten()-mfit)**2 * frem[frem.columns[1:]],axis=1) / ccap )
    
        plt.figure()
        plt.ylim((0,12))
        plt.errorbar(frem.g, np.sum(mfit,axis=1), yerr=sfit, fmt='k--', ecolor='0.7')
        plt.xlabel('generations')
        plt.ylabel('mean fitness')
        plt.tight_layout()
        plt.savefig('plots/{:.0e}/{}_{:.0e}.svg'.format(ccap, rep, ccap))
        plt.close()
        
    sys.stdout.write('\r[{}{}] {}% {}\t'.format('#'*int(((i+1)*50)/no_f),
                     '-'*int(((no_f-i-1)*50)/no_f), int(((i+1)*100)/no_f), files[i][:-4]))
    sys.stdout.flush()
sys.stdout.write('\n Done!')
sys.stdout.flush()
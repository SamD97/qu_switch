#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Plotter for simulations.
    @author: Saumil Shah'''

import os
import gc
import sys
import numpy as np
import pandas as pd
import seaborn as sb
from matplotlib import rc
from matplotlib import pyplot as plt

rc('font', size=16)
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)

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
repsn=['rep 0', 'rep 1', 'rep 2', 'rep 3', 'rep 4',
       'rep 5', 'rep 6', 'rep 7', 'rep 8', 'rep 9']

for i in range(no_f):
    gc.collect()
    
    mgen = int(files[i].split('_')[0])
    ccap = float(files[i].split('_')[1].split('.')[0])
    #if not os.path.isdir('plots/{:.0e}'.format(ccap)):
        #os.mkdir('plots/{:.0e}'.format(ccap))
    df = pd.read_csv(files[i], '\t')
    
    plt.figure(figsize=(8,16))
    plt.subplot(311)
    plt.title('Fitness vs Generations',fontsize=18)
    plt.ylim((3,10))
    for rep in pd.unique(df.r):
        frem = df[df.r==rep].drop('r',1)
        mfit = np.sum(fitl.flatten() * frem[frem.columns[1:]] / ccap,1)
    
        plt.plot(frem.g, mfit, label='rep {}'.format(rep))
    plt.xlabel('Generations')
    plt.ylabel('Fitness')
    plt.legend()
    
    frem = df.groupby('g').mean().drop('r',1).iloc[mgen].values
    fred = df.groupby('g').std().drop('r',1).iloc[mgen].values / np.sqrt(pd.unique(df.r).size)
    plt.subplot(312)
    plt.title('Relative Fraction vs Genotypes',fontsize=18)
    plt.xticks(rotation=15)
    plt.ylim((0,1))
    uplm = ((frem+fred)/ccap > 0.99)
    lwlm = ((frem-fred)/ccap < 0.01)
    fuer = np.copy(fred/ccap)
    fler = np.copy(fred/ccap)
    fuer[uplm] = 0
    fler[lwlm] = 0
    plt.bar(names[2:], frem/ccap, yerr=[fler, fuer])
    plt.xlabel('Genotypes')
    plt.ylabel('Relative Fraction')
    
    frem = df[df.g==mgen].drop('g',1)
    frem = frem[frem.columns[1:]]/ccap
    plt.subplot(313)
    plt.title('Genotypes vs Replicates',fontsize=18)
    sb.heatmap(frem.T,0,1,cmap='magma',linewidths=0.1)
    plt.yticks(rotation=15)
    plt.xticks(np.arange(10)+0.5,repsn,rotation=15)
    plt.ylabel('Genotypes')
    plt.xlabel('Replicates')
    plt.suptitle('Population size = $10^{:.0f}$'.format(len(str(ccap)[1:-2])), y=0.99)
    plt.tight_layout(rect=(0,0,1,0.97))
    plt.savefig('plots/combi_{}_{:.0e}.svg'.format(mgen, ccap))
    plt.close()  
        
    sys.stdout.write('\r[{}{}] {}% {}\t'.format('#'*int(((i+1)*50)/no_f),
                     '-'*int(((no_f-i-1)*50)/no_f), int(((i+1)*100)/no_f), files[i][:-4]))
    sys.stdout.flush()
    
sys.stdout.write('\n Done!')
sys.stdout.flush()
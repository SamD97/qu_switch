#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Combiner for simulation data.
    @author: Saumil Shah'''

import os
import gc
import sys
import pandas as pd
from numpy.random import randint

path = input(' Please enter path to folder: ')
if path == '': path = os.getcwd()
os.chdir(path)

files=[]
for f in os.listdir(path):
        if f.endswith('txt'):
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
print (' \t\tTotal {} .txt files were found.'.format(no_f))
print (' {}'.format('='*70))
print ()

for i in range(no_f):
    gc.collect()
    
    mgen = int(files[i].split('_')[2])
    ccap = float(files[i].split('_')[4][:-4])
    fn = '{}_{:.0e}.pop'.format(mgen, ccap)
    
    df = pd.read_csv(files[i], '\t')
    df.r = df.r.replace(pd.unique(df.r), randint(0,1000,len(pd.unique(df.r))))
    
    if not os.path.isfile(fn):
       df.to_csv(fn, '\t', index=False)
    else: 
       df.to_csv(fn, '\t', index=False, mode='a', header=False)    
    
    sys.stdout.write('\r[{}{}] {}% {}\t'.format('#'*int(((i+1)*50)/no_f),
                     '-'*int(((no_f-i-1)*50)/no_f), int(((i+1)*100)/no_f), files[i][:-4]))
    sys.stdout.flush()

sys.stdout.write('\n Done!')
sys.stdout.flush()
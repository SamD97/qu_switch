#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Code to plot fitness landscape 
    @author: Saumil Shah'''

import numpy as np
import defn as dn
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xticks([0,1,2])
ax.set_xticklabels(['B0','B1','B2'])
ax.set_yticks([-0,-1,-2])
ax.set_yticklabels(['A0','A1','A2'])
X, Y = np.meshgrid(np.arange(3),-np.arange(3))
cset = ax.plot_wireframe(X, Y, dn.fitl, rstride=0)
plt.savefig('plots/fitness_landscape.svg')
plt.close()
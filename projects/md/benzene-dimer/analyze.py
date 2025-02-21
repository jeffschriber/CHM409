#! /usr/bin/env python

import numpy as np
import pandas as pd
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

rgb_colors = [( 23, 71,158), # dark blue 
              (189, 32, 46), # red 
              (  0,161, 75), # green
              ( 76,150,209), # light blue
              (242,101, 34), # orange
              #(159,110,175), # lavender
              (109,110,175), # lavender
             ]

colors = []
for c in rgb_colors:
    colors.append(tuple([float(p) / 255.0 for p in c]))


plt.rcParams['xtick.major.pad']='2'
plt.rcParams['ytick.major.pad']='2'
mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size']=11
plt.rcParams['figure.figsize'] = 3.375*2, 2.25*2

fig, ax= plt.subplots()


time = []
stre = []
bend = []
tors = []
nonb = []
with open('file.txt','r') as infile:
    for t, line in enumerate(infile):
        line = line.split('kJ/mol')
        print(line)
        try:
            stre.append(float(line[0]))
            bend.append(float(line[1]))
            tors.append(float(line[2]))
            nonb.append(float(line[3]))
            time.append(float(t))
        except:
            pass
stre = np.asarray(stre)
bend = np.asarray(bend)
tors = np.asarray(tors)
nonb = np.asarray(nonb)

print(stre)

tots = stre + bend + tors + nonb

idx = []
for n, item in enumerate(tots):
    if item <= 20:
        idx.append(n)
        print(f"Frame {n} energy {item}")

ax.plot(time, stre, marker='', linestyle='-', color=colors[0], label='Stretch')
ax.plot(time, bend, marker='', linestyle='-', color=colors[1], label='Bend')
ax.plot(time, tors, marker='', linestyle='-', color=colors[2], label='Torsion')
ax.plot(time, nonb, marker='', linestyle='-', color=colors[3], label='Non-Bonded')
ax.plot(time, tots, marker='', linestyle='-', color='k', label='Total')
plt.legend()
plt.savefig('ens.pdf')


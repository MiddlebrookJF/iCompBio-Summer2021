# File: ScanAnalysis.py
# Project: The local topological free energy of SARS-CoV-2

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean

def plot_vas(proteinName, scanlength, section, vas, maxVas=False):
    knot = pd.read_csv(f'Coordinates/{proteinName}.csv').values.tolist()

    x, y, z = [], [], []
    section = section.split(':')
    section[1] = section[1][0:-1]
    for i in range (int(section[0]), int(section[1])):
        x.append(knot[i][0])
        y.append(knot[i][1])
        z.append(knot[i][2])

    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z)

    section = section[0] + '-' + section[1]
    if maxVas:
        plt.title(f'{proteinName} with Max Vas for {section} of {vas}')
    else:
        plt.title(f'{proteinName} at {section} with Vas of {vas}')
    if not os.path.exists(f'Vas-Data/Scans/{proteinName}'):
        os.mkdir(f'Vas-Data/Scans/{proteinName}')
    if not os.path.exists(f'Vas-Data/Scans/{proteinName}/{scanlength}'):
        os.mkdir(f'Vas-Data/Scans/{proteinName}/{scanlength}')
    plt.savefig(f'Vas-Data/Scans/{proteinName}/{scanlength}/{section}.png')
    plt.clf()

print('\n\n\n')
with open('Vas-Data/1000Scan0.txt') as scanFile:
    lines = scanFile.readlines()

proteinName = ''
max_loc = 0
scanlength = 200
for line in lines:
    words = line.split(' ')
    if len(words) < 2: continue
    if words[1] == 'length':
        proteinName = words[3]
        scanlength = 200
    elif words[1] == 'at':          #Plot section at this line along with vas
        plot_vas(proteinName, scanlength, words[2], words[0])
    elif words[1] == 'maximum':
        max_loc = words[-2][1:-1] + ':' + words[-1][0:-1]
        plot_vas(proteinName, int(words[5]), max_loc, words[7][0:-1], maxVas=True)
        scanlength += 200

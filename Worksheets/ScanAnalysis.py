# File: ScanAnalysis.py
# Project: The local topological free energy of SARS-CoV-2

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean

#Used by plot_all to plot an indiviual protein section and then save the figure as a png
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

def plot_all(lines):
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

def plot_vas_change(lines):
    proteinName = ''
    scanlength = 200
    vas_list = []
    start_list = []

    for line in lines:
        words = line.split(' ')
        if len(words) < 2: continue
        if words[1] == 'length':
            proteinName = words[3]
            scanlength = 200
        if words[1] == 'at':
            vas_list.append(float(words[0]))
            start_list.append(words[2].split(':')[0])
        if words[1] == 'maximum':
            plt.plot(start_list, vas_list)
            plt.title(f'Change in V2 for {proteinName} at {scanlength} scan length')
            plt.xlabel('Starting Point')
            plt.ylabel('Local Second Vassiliev Measure')

            if not os.path.exists(f'Vas-Data/Change/{proteinName}'):
                os.mkdir(f'Vas-Data/Change/{proteinName}')
            plt.savefig(f'Vas-Data/Change/{proteinName}/{proteinName}{scanlength}.png')
            plt.clf()
            start_list.clear()
            vas_list.clear()
            scanlength += 200

def vas_matrix(lines):
    proteinName = ''
    scanlength = 1
    matrix = []

    for line in lines:
        words = line.split(' ')
        if len(words) < 2: break
        if words[1] == 'at':        #Add the section and its V2 to the matrix
            section = words[2].split(':')
            section[0], section[1] = int(section[0]), int(section[1])
            matrix[section[0]][section[1]-1] = float(words[0]) * 100        #Puts v2 into the spot
            continue
        if words[1] == 'length':
            proteinName = words[3]
            proteinLength = int(words[5][0:-2])
            matrix = np.zeros((proteinLength, proteinLength), dtype=float)
        if words[1] == 'maximum':
            with open(f"Vas-Data/Matrices/{proteinName}.csv", mode='w', newline='') as f:
                pd.DataFrame(matrix).to_csv(f, header = f.tell()==0)
            scanlength += 1

with open('Vas-Data/All-50inter-one.txt') as scanFile:
    lines = scanFile.readlines()

plot_vas_change(lines)

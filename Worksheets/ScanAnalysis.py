# File: ScanAnalysis.py
# Project: The local topological free energy of SARS-CoV-2

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean

from pandas.io.pytables import dropna_doc

#Used by plot_all to plot an indiviual protein section and then save the figure as a png
def plot_vas(proteinName, scanlength, section, vas, maxVas=False):
    knot = pd.read_csv(f'Coordinates/{proteinName}.csv').values.tolist()

    x, y, z = [], [], []
    section = section.split(':')
    section[1] = section[1]
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
    if not os.path.exists(f'Vas-Data/Scan-Graphs/{proteinName}'):
        os.mkdir(f'Vas-Data/Scan-Graphs/{proteinName}')
    if not os.path.exists(f'Vas-Data/Scan-Graphs/{proteinName}/{scanlength}'):
        os.mkdir(f'Vas-Data/Scan-Graphs/{proteinName}/{scanlength}')
    plt.savefig(f'Vas-Data/Scan-Graphs/{proteinName}/{scanlength}/{section}.png')
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
            #Change start and end to the appropriate indices based on the dictionary of coord_index to amino_index
            amino_index_cols = ['Coord_Index', '6acd', '7krq', '7lwt', '7lws', '7lww', '7lyn', '6zgh',
                    '6zge', '6zgi', '6xkl', '7lyl', '7kdk', '6zgg', '7m8k', '7mjg', '6xra']
            col_index = [i for i in range(len(amino_index_cols)) if amino_index_cols[i] == proteinName]
            amino_indices = pd.read_csv('Vas-Data/AminoAcid-Indices.csv', usecols=col_index)
            amino_indices = amino_indices.dropna().values.tolist()

            old_section = words[2].split(':')
            section = str(int(amino_indices[int(old_section[0])][0])) + ':'
            section += str(int(amino_indices[int(old_section[1][0:-2])][0]))
            plot_vas(proteinName, scanlength, section, words[0])
        elif words[1] == 'maximum':
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

            if not os.path.exists(f'Vas-Data/Change-Graphs/{proteinName}'):
                os.mkdir(f'Vas-Data/Change-Graphs/{proteinName}')
            plt.savefig(f'Vas-Data/Change-Graphs/{proteinName}/{proteinName}{scanlength}.png')
            plt.clf()
            start_list.clear()
            vas_list.clear()
            scanlength += 200

def plot_all_change(lines, scanlengths=[200,400,600], title='Change in V2'):
    proteinName = ''
    vas_list = []
    start_list = []
    proteinNames = {
        '6zge': 'Uncleaved Closed',
        '6zgi': 'Cleaved Closed',
        '6zgg': 'Wild Cleaved',
        '6zgh': 'Cleaved Intermediate',
        '7lww': 'Brazil 3-Mutant',
        '7lyn': 'South African',
        '7lwt': 'United Kingdom',
        '6xkl': 'Hexapro',
        '6acd': '2003 SARS-CoV',
        '7lws': 'United Kingdom',
        '7lyl': 'South African',
        '7kdk': 'D614G-Mutated'
    }

    for scanlength in scanlengths:
        currentScanLength = 0
        colors = ['black', 'darkmagenta', 'sandybrown', 'cornflowerblue', 'red', 'forestgreen']
        plt.figure(figsize=(10, 6)) 
        for line in lines:
            words = line.split(' ')
            if len(words) < 2: continue
            if words[1] == 'length':
                proteinName = proteinNames[words[3]]
                proteinLabel = words[3]
            if words[0] == 'Scanlength':
                currentScanLength = int(words[2][0:-1])
                continue
            if currentScanLength == scanlength:
                if words[1] == 'at':
                    #If current graph being created is for the scanlength at the lines being read, then add those
                    vas_list.append(float(words[0]))
                    start_list.append(words[2].split(':')[0])
                if words[1] == 'maximum':

                    #Change start_list to the appropriate indices based on the dictionary of coord_index to amino_index
                    amino_index_cols = ['Coord_Index', '6acd', '7krq', '7lwt', '7lws', '7lww', '7lyn', '6zgh',
                            '6zge', '6zgi', '6xkl', '7lyl', '7kdk', '6zgg', '7m8k', '7mjg', '6xra']
                    col_index = [i for i in range(len(amino_index_cols)) if amino_index_cols[i] == proteinLabel]
                    amino_indices = pd.read_csv('Vas-Data/AminoAcid-Indices.csv', usecols=col_index)
                    amino_indices = amino_indices.dropna().values.tolist()

                    amino_start_list = []
                    for start in start_list:
                        amino_start_list.append(int(amino_indices[int(start)][0]))
                    
                    plt.plot(amino_start_list, vas_list, 'o-', label=proteinName, color=colors.pop())
                    start_list.clear()
                    vas_list.clear()

        plt.title(title, y=1.035)
        plt.suptitle(f'Scan Length {scanlength}', y=0.915, fontsize=10)
        plt.xlabel('Amino Acid Starting Point')
        plt.ylabel('Local Second Vassiliev Measure')
        plt.legend()
        # plt.savefig(rf'Vas-Data\Change-Graphs\Groups\6zg{scanlength}.png')
        plt.savefig(rf'Vas-Data\Change-Graphs\Groups\RBD-Up{scanlength}.png')
        # plt.savefig(rf'Vas-Data\Change-Graphs\Groups\RBD-Down{scanlength}.png')
        plt.clf()           #Clears the graph


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

with open('Vas-Data/Scan-Text/50inter-RBD_Up.txt') as scanFile:
    lines = scanFile.readlines()
# plot_all_change(lines, title='Change in V2 for Wild SARS-CoV-2 S Proteins in Various Conformations')
plot_all_change(lines, title='Change in V2 for SARS-CoV-2 Variant S Proteins in Open RBD-Up Conformation')
# plot_all_change(lines, title='Change in V2 for SARS-CoV-2 Variant S Proteins in Closed RBD-Down Conformation')

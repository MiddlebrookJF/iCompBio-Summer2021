#!/usr/bin/python
import os
import matplotlib as mpl
mpl.use('Agg')
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle as pl
from operator import itemgetter
from matplotlib.colors import BoundaryNorm

def matrix_vis(matrixDF):

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_aspect('equal')
    cmap = plt.get_cmap('PuOr')
    # define the colormap
    cmap = plt.get_cmap('PuOr')
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    print("\n\nMatrix:")
    print(matrixDF)
    matrix = np.matrix(np.array((matrixDF)))
    print(matrix)

    vm = max(abs(np.min(matrix)), abs(np.max(matrix)))
    norm = mpl.colors.Normalize(vmin=-vm, vmax=vm, clip=False)

    plt.imshow(matrix, interpolation='none', cmap=cmap, norm=norm) #cmap=plt.cm.ocean)
    plt.colorbar()
    #plt.show()
    plt.savefig(f"draft.png")

def read_matrix(path):
    with open(path, 'r') as scanFile:
        lines = scanFile.readlines()
    protein_matrix = []
    for line in lines[1: len(lines)]:
        values = line.split(',')
        values[-1] = values[-1][0:-2]
        for val in values[1: len(values)]:
            val = float(val)
        protein_matrix.append(values)

    return pd.DataFrame(protein_matrix, columns=lines[0][1: len(lines[0])])

#Executed stuff below here
protein_matrix = read_matrix("Vas-Data/Matrices/6zge_700-1100.csv")
matrix_vis(protein_matrix)

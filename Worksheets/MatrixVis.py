#!/usr/bin/python
import os
import matplotlib as mpl
# mpl.use('Agg')
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle as pl
from operator import itemgetter
from matplotlib.colors import BoundaryNorm

def matrix_vis(matrixDF):
    print("\n\nMatrix:")
    print(matrixDF)

    matrix_range = (int(matrixDF.columns[0]), int(matrixDF.columns[-1]), int(matrixDF.columns[0]), int(matrixDF.columns[-1]))

    fig = plt.figure()
    ax = fig.add_subplot(111)                               # means only one plot
    ax.set_aspect('equal')
    cmap = plt.get_cmap('PuOr')                             # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]             # extract all colors from the .jet map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)  # create the new map

    matrix = np.matrix(matrixDF)
    print(matrix)

    vm = max(abs(np.min(matrix)), abs(np.max(matrix)))
    norm = mpl.colors.Normalize(vmin=-vm, vmax=vm, clip=False)

    plt.imshow(matrix, interpolation='none', cmap=cmap, norm=norm, extent=matrix_range) #cmap=plt.cm.ocean)
    plt.colorbar()
    ax.set_title('Uncleaved Closed S: V2 Values in Matrix Plot')
    plt.savefig(f"Vas-Data/Matrices/Matrix_{matrix_range[0]}-{matrix_range[1]}.png")
    plt.show()

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
# protein_matrix = read_matrix('Vas-Data/Matrices/6zge_700-1100.csv')
proteinDF = pd.read_csv('Vas-Data/Matrices/6zge_700-1100.csv', index_col=0)
matrix_vis(proteinDF)

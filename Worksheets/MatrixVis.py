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

    start = int(matrixDF.columns[0])
    end = int(matrixDF.columns[-1])
    matrix_range = (start, end, start, end)

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
    ax.set_title('Uncleaved Closed S: V2 Values in Fingerprint Matrix')
    plt.savefig(f"Vas-Data/Matrices/NewMatrix_{matrix_range[0]}-{matrix_range[1]}.png")
    plt.show()

#Executed stuff below here
proteinDF = pd.read_csv('Vas-Data/Matrices/6zge_0-1100.csv', index_col=0)
matrix_vis(proteinDF)

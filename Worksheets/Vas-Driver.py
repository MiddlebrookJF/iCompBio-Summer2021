import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from .Vassiliev import Vas

def vas_closed(knot):
    vas = Vas.vas_closed(knot)
    print("\nThe Vassiliev Measure of the closed knot is " + str(vas) + "\n")

    x, y, z = [], [], []
    for i in range (0, len(knot)):
        x.append(knot[i][0])
        y.append(knot[i][1])
        z.append(knot[i][2])
    x.append(knot[0][0])
    y.append(knot[0][1])
    z.append(knot[0][2])

    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z)
    plt.show()

def vas_open(knot):
    vas = Vas.vas_open(knot)
    print("\nThe Vassiliev Measure of the open 3d curve is " + str(vas) + "\n")

    x, y, z = [], [], []
    for i in range (0, len(knot)):
        x.append(knot[i][0])
        y.append(knot[i][1])
        z.append(knot[i][2])

    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z)
    plt.show()

#########################

spikeDF = pd.read_csv(r'Coordinates\6zge.csv')
print(spikeDF.head())
spikeList = spikeDF.values.tolist()                 #change df to a list of atoms' coordinates

vas_open(spikeList)
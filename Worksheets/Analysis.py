# File: Analysis.py
# Project: The local topological free energy of SARS-CoV-2

import numpy as np
import pandas as pd
from statistics import mean

print('\n\n\n\n')
proteinsDF = pd.read_csv('Vas-Data/ProteinVas.csv')
milleDF = proteinsDF[proteinsDF['NumProjections'] == 1000]
centDF = proteinsDF[proteinsDF['NumProjections'] == 100]

# milleSTD = dict()
# #For every protein name, get the standard deviation for the Vas of that protein's rows
# for pName in milleDF['Name'].unique():
#     thisVasDF = milleDF[milleDF['Name'].isin([pName])]['Vassiliev']
#     milleSTD[pName] = np.std(thisVasDF, ddof=1)

def avgVasDF(df):
    dfAVG = []
    #For every protein name, find average Vassiliev
    for pName in df['Name'].unique():
        #Get the average for the Vas of that protein's rows
        thisProteinDF = df[df['Name'].isin([pName])]
        thisMean = mean(thisProteinDF['Vassiliev'])
        dfAVG.append((pName, thisMean))
    return pd.DataFrame(dfAVG, columns=['Name', 'AvgVas'])

milleAVG = avgVasDF(milleDF)
centAVG = avgVasDF(centDF)

allSTD = []
for pName in milleAVG['Name']:
    #thisVasList is the list containing the two AvgVas values for protein pName
    print(pName)
    print('\t' + milleAVG[milleAVG['Name'].isin([pName])]['AvgVas'].to_string())
    print('\t' + centAVG[centAVG['Name'].isin([pName])]['AvgVas'].to_string())
    thisVasList = milleAVG[milleAVG['Name'].isin([pName])]['AvgVas'].append(centAVG[centAVG['Name'].isin([pName])]['AvgVas'])
    allSTD.append((pName, np.std(thisVasList, ddof=1)))
allSTD = pd.DataFrame(allSTD, columns=['Name', 'VasSTD'])

print(allSTD)
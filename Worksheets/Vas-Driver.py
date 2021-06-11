import time
import pandas as pd
import Read_File as read
import Vas_Defs_JW as Vas

# proteins = read.readAll_pdb(fr"Proteins-PDB")
# proteins = sorted(proteins.items(), key=lambda x: len(x[1])) #returns tuples of (protein, chain) sorted by chain length

# #estimate vas of each point
# vasValues = {}
# for tuple in proteins:
#     # value = None
#     # while value == None: #ensure an actual return result before moving on
#     startTime = time.time()
#     value = Vas.vas_open_parallel(tuple[1], 100, size=20)
#     execTime = Vas.runtime(startTime)
#     if(value!=None):
#         print (tuple[0], ':' , len(tuple[1]))
#         print('Vas: %f' %(value))
#         print('Runtime: %f \n'%(execTime))
#         vasValues.update({tuple[0]:value})

# protein = read.read_pdb("Proteins-PDB/6zge.pdb")

# value = Vas.vas_open_parallel(protein)
# print(value)

# proteinList = [[1, 0, 0],
#             [4, 0, 0],
#             [1, 6, 2],
#             [0, 2, -5],
#             [5, 2, 5],
#             [4, 6, -2]]

# proteinName = "6acd"
# proteinDF = pd.read_csv(fr'Coordinates\{proteinName}.csv')
# proteinList = proteinDF.values.tolist()

# measure = Vas.vas_measure(proteinList, closed=True)
# print(measure)
# #vas is 1.0 for closed trefoil


### Trefoil for testing ###

# proteinList = [[1, 0, 0],             #trefoil
#             [4, 0, 0],
#             [1, 6, 2],
#             [0, 2, -5],
#             [5, 2, 5],
#             [4, 6, -2]]

# startTime = time.time()
# value = vas_measure(proteinList, closed=True)
# execTime = runtime(startTime)
# if(value!=None):
#     print (proteinList, ':' , len(proteinList))
#     print(f'Vas: {value}')
#     print(f'Runtime: {execTime} seconds or {execTime/60} minutes\n')
#vas is 1.0 for closed trefoil

    # trefoil = [[1, 0, 0],
    #             [4, 0, 0],
    #             [1, 6, 2],
    #             [0, 2, -5],
    #             [5, 2, 5],
    #             [4, 6, -2],
    #             [0.5, 0.5, 0.5]]

    # plot_vas(trefoil)
    # #vas is 0.995 for open trefoil
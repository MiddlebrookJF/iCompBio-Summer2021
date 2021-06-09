
import glob

proteinList = glob.glob("Coordinates/*.csv")
for proteinPath in proteinList:
    print(proteinPath[-8:-4])
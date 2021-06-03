import glob
import os

#read pdb file and find coordinates of CA of first model
def read_pdb(f):
    with open(f, 'r') as reader:
        data = []
        for line in reader:
            data.append(line.split())

    #search for those that contain CA and list there coordinates
    alphaCCoords = []
    for i in range(len(data)):
        if(data[i][0] == "ATOM"):
            if(data[i][2] == "CA"): #split to prevent out of bounds errors
                if(data[i][4].isalpha()): #with AA>1000 a space is deleted and there is a shift in the table
                    alphaCCoords.append([float(data[i][6]), float(data[i][7]), float(data[i][8])])
                else:
                    alphaCCoords.append([float(data[i][5]), float(data[i][6]), float(data[i][7])])  
     
        #break at the end of the first model or segment
        if(data[i][0] == "ENDMDL" or data[i][0] == "TER"):
            break;

    return alphaCCoords

#in current file location, read all pdb files and return dictionary of name: coordinate list
def readAll_pdb():
    pdbFiles = glob.glob('*.pdb') #find all files with extension .pdb
    dictNameCoords = {}
    for file in pdbFiles:
        fname, fext = os.path.splitext(file) #split file name from extension
        dictNameCoords.append({fname:read_pdb(file)}) #add to dict name and its coords
    return dictNameCoords

#read all pdb files in a specific file location
def readAll_pdb(fileLocation):
    cwd = os.getcwd
    os.chdir(fileLocation)
    pdbFiles = glob.glob('*.pdb') #find all files with extension .pdb
    dictNameCoords = {}
    for file in pdbFiles:
        fname, fext = os.path.splitext(file) #split file name from extension
        dictNameCoords.update({fname:read_pdb(file)}) #add to dict name and its coords
    return dictNameCoords
    os.chdir(cwd)

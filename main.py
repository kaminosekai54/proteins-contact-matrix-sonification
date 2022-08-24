import numpy as np
import pandas as pd
import Bio.PDB
from conkit.io.pdb import PdbParser
import conkit.plot
from matplotlib import pylab
import pygame.midi as pm
from midiutil import MIDIFile
from midi2audio import FluidSynth
import os
import random




################################################################################
# global variable

Octave={
"0":[0,1,2,3,4,5,6,7,8,9,10,11],
"1":[12,13,14,15,16,17,18,19,20,21,22,23],
"2":[24,25,26,27,28,29,30,31,32,33,34,35],
"3":[36,37,38,39,40,41,42,43,44,45,46,47],
"4":[48,49,50,51,52,53,54,55,56,57,58,59],
"5":[60,61,62,63,64,65,66,67,68,69,70,71],
"6":[72,73,74,75,76,77,78,79,80,81,82,83],
"7":[84,85,86,87,88,89,90,91,92,93,94,95],
"8":[96,97,98,99,100,101,102,103,104,105,106,107],
"9":[108,109,110,111,112,113,114,115,116,117,118,119],
"10":[120,121,122,123,124,125,126,127, 127, 127, 127, ]
}

noteSymbToNoteNumber = {
"C": [noteNumber[0] for noteNumber in Octave.values()],
"C#": [noteNumber[1] for noteNumber in Octave.values()],
"D": [noteNumber[2] for noteNumber in Octave.values()],
"D#": [noteNumber[3] for noteNumber in Octave.values()],
"E": [noteNumber[4] for noteNumber in Octave.values()],
"F": [noteNumber[5] for noteNumber in Octave.values()],
"F#": [noteNumber[6] for noteNumber in Octave.values()],
"G": [noteNumber[7] for noteNumber in Octave.values()],
"G#": [noteNumber[8] for noteNumber in Octave.values()],
"A": [noteNumber[9] for noteNumber in Octave.values()],
# "A#": [noteNumber[10] for noteNumber in Octave.values()],
"B": [noteNumber[10] for noteNumber in Octave.values()],
}

AAMap = {
"ALA": {"note" : "C", "instrument": ""},
"ARG": {"note": "C", "instrument" : ""},
"ASN": {"note": "C", "instrument": ""},
"ASP": {"note": "C", "instrument": ""},
"CYS": {"note": "D", "instrument": ""},
"GLU": {"note": "E", "instrument": ""},
"GLN": {"note": "E", "instrument": ""},
"GLY": {"note": "E", "instrument": ""},
"HIS": {"note": "F", "instrument": ""},
"ILE": {"note": "F", "instrument": ""},
"LEU": {"note": "G", "instrument": ""},
"LYS": {"note": "G", "instrument": ""},
"MET": {"note": "G", "instrument": ""},
"PHE": {"note": "A", "instrument": ""},
"PRO": {"note": "A", "instrument": ""},
"SER": {"note": "A", "instrument": ""},
"THR": {"note": "B", "instrument": ""},
"TRP": {"note": "B", "instrument": ""},
"TYR": {"note": "B", "instrument": ""},
"VAL": {"note": "D", "instrument": ""},
}

################################################################################
# using conkit

def test(list):
    for i in range(max(list)):
        if i not in list: print("not in list : ", i)

def checkMissingValue(dict):
    if len(dict) < max(dict.keys()):
        for i in range(1, max(dict.keys())):
            if i not in dict.keys() :
                print("new value to add : ", i)
                dict[i] = "Unknown"

        return dict
    else: return dict

def getContactMatrixWithConkit(pdbFile, pdbID):
    p = PdbParser()
    with open(pdbFile, "r") as pdb_fhandle:
        pdb = p.read(pdb_fhandle, f_id=pdbID, atom_type="CA")
        for cmap in pdb:
            print(cmap.get_unique_distances())
            rowId = {}
            columnId = {}
            idDict = {}

        # fig = conkit.plot.ContactMapMatrixFigure(cmap)
            # fig2 = conkit.plot.ContactMapFigure(cmap)
            # fig.savefig("testMatrixMap.png", overwrite=True)
            # fig2.savefig("testContactMap.png", overwrite=True)
            for index in range(len(cmap)):
                dist = cmap[index]
                if dist.id[0] not in rowId .keys(): 
                    rowId[dist.id[0]] = dist.res1
                if dist.id[1] not in columnId.keys(): 
                    columnId[dist.id[1]] = dist.res2

            print(len(rowId))
            print(len(columnId))
            rowId= checkMissingValue(rowId)
            columnId= checkMissingValue(columnId)
            print(len(rowId))
            print(len(columnId))
            matrix = np.zeros((len(rowId), len(columnId)), np.float)
            for dist in cmap:
                row_id = dist.id[0]-1
                column_id = dist.id[1]-1
                matrix[[row_id], column_id] = dist.raw_score
            np.save("distance_matrix.txt", matrix)
            df = pd.DataFrame(matrix, columns=columnId.values(), index=rowId.values())
            df.to_csv("matrix_df.csv")



################################################################################
# using biopython


# function computeResidueDist
# Returns the C-alpha distance between two residues
# @param,
# @residue1, the first residue
# @residue2, the second residue
def computeResidueDist(residue1, residue2) :
    vectorDiff = residue1.get_coord() - residue2.get_coord() 
    return np.sqrt(np.sum(vectorDiff  * vectorDiff ))

# computeDistMatrix
# Returns a matrix of C-alpha distances 
def computeDistMatrix(structure) :
    CAList = [atom for atom in structure.get_atoms() if atom.name == "CA"]
    residuesList = structure.get_residues()
    residuesID= []
    for residue in residuesList:
        if residue.get_resname() not in residuesID: residuesID.append(residue.get_resname())

    tmp = [x for x in residuesID if x not in AAMap .keys()]
    if len(tmp) > 0 : print(tmp)
    for i in range(len(list(residuesList))):
        curRes = list(residuesList)[i]
        curCA=[a for a in curRes.get_atoms() if a.name == "CA"]
        if len(curCA) > 1 : print("should have only one atome")
        for ca in curCA:
            if ca != CAList[i] : print("soucis")


    distMatrix= np.zeros((len(CAList), len(CAList)), np.float)
    for row, ca1 in enumerate(CAList) :
        for col, ca2 in enumerate(CAList) :
            distMatrix[[row], col] = computeResidueDist(ca1, ca2)
            # distMatrix[row, col] = computeResidueDist(ca1, ca2)
    return distMatrix




def getPDBStruct(pdbFile, pdbCode):
    structure = Bio.PDB.PDBParser().get_structure(pdbCode, pdbFile)
    model = structure[0]
    return structure 


def mapValue(value, min_value, max_value, min_result, max_result):
    result = min_result + (value - min_value)/(max_value - min_value) * (max_result - min_result)
    return result

def getAANote(AA, id, length):
    note= AAMap[AA]["note"]
    octave = int(mapValue(id, 0, length, 0, 127) /12)
    noteNumber =noteSymbToNoteNumber[note][octave]
    return noteNumber

def getContactMatrixWithBiopython(pdbFile, pdbID, contactThreshold):
    structure = getPDBStruct(pdbFile, pdbID)
    distMatrix = computeDistMatrix(structure)
    CAList = [atom for atom in structure.get_atoms() if atom.name == "CA"]
    residuesList = list(structure.get_residues())
    contactMatrix = distMatrix < contactThreshold
    # showContactMap(distMatrix , contactMatrix )
    noteList = mapping(distMatrix, residuesList, CAList, contactThreshold)
    creatMidiFile(noteList)

    print("Maximum distance", np.max(distMatrix))
    print("minimum distance", np.min(distMatrix))

def showContactMap(distMatrix , contactMatrix):
    pylab.matshow(np.transpose(distMatrix))
    pylab.colorbar()
    pylab.show()
    pylab.savefig("biopython_distance_map.png")
    pylab.autumn()
    pylab.imshow(np.transpose(contactMatrix))
    pylab.show()
    pylab.savefig("biopython_contact_map.png")


def mapping(distMatrix, residuesList, CAList, contactThreshold):
    noteList = []
    print(len(distMatrix))
    time = 0
    duration = 1 # longueur d'une noire
    for line, dist1 in enumerate(distMatrix):
        aaName = residuesList[line].get_resname()
        # aaName = residuesList[line]
        note = getAANote(aaName, line, len(distMatrix))
            #  element in this order :
            # track, channel, note, time, duration, volume
        # noteList.append((0, 0, note, time, duration, 50))
        # time +=1
        # isAdded = False

        for col, dist2 in enumerate(distMatrix):
            if (distMatrix[[line], col] > 0 and distMatrix[[line], col]  <= contactThreshold) and (line != col):
                aaInContact= residuesList[col].get_resname()
                # aaInContact= residuesList[col]
                contactNote = getAANote(aaInContact, col, len(distMatrix))
                volume = int(mapValue(distMatrix[line, col], 0, contactThreshold, 127, 20))
                # if not isAdded : 
                    # noteList.append((0, 0, note, time, duration, 50))
                    # time +=1
                    # isAdded = True
                noteList.append((0, 0, note, time, duration, volume))
                noteList.append((0, 0, contactNote, time, duration, volume))
                time +=1

    return noteList

def creatMidiFile(noteList):
    time = 0 # in bits
    duration = 1 # in bits
    tempo = 100 # in bpm (bits per minutes)
    print(len(noteList))
    MyMIDI = MIDIFile(2)
    MyMIDI.addTempo(0, time, tempo)
    for i in range(len(noteList)):
        track, channel, note, time, duration, volume = noteList[i]
        # print(track, channel, note, time, duration, volume)
        
        if note < 0 : note = 0
        MyMIDI.addNote(track, channel, note, time, duration, volume)

    if os.path.isfile("test.mid") : os.remove("test.mid")
    with open("test.mid", "wb") as output_file:
        MyMIDI.writeFile(output_file)

    # fs = FluidSynth()
    # fs.midi_to_audio(os.path.abspath('test.mid'), os.path.abspath('output.wav'))
    # fs.midi_to_audio(os.path.abspath('output.mid'), os.path.abspath('output.wav'))


    print("finish")
    


def generateTestData(matrixSize, nbValueToChange, contactThreshold):
    testMatrix= np.zeros((matrixSize, matrixSize), np.float)
    testResList = ['ALA', 'SER', 'ALA', 'CYS', 'LYS', 'LEU', 'HIS', 'LEU', 'ALA', 'PRO', 'GLY', 'ARG', 'PHE', 'MET', 'THR', 'PHE', 'GLN', 'TYR', 'ASN', 'THR']
    # for  i in range(20):
        # testResList .append(list(AAMap.keys())[random.randint(0, len(AAMap.keys())-1)])


    testMatrix[[0], matrixSize-2] = 5.0
    testMatrix[[matrixSize-2], matrixSize-1] = random.randint(1, contactThreshold)
    nbModif = 0
    modified = []
    while nbModif < nbValueToChange:
        line = random.randint(0, matrixSize-1)
        col = random.randint(0, matrixSize-1)
        if not (line, col) in modified and testMatrix[[line], col] == 0:
            testMatrix[[line], col] = random.randint(1, contactThreshold)
            modified.append((line, col))
            nbModif +=1
    print(testMatrix)
    print(testResList)
    np.save("testMatrix.txt", testMatrix)

    contactMatrix = testMatrix < contactThreshold
    # showContactMap(testMatrix, contactMatrix )
    noteList = mapping(testMatrix, testResList, [], contactThreshold)
    creatMidiFile(noteList)

 






################################################################################
# main function
def main():
    print("matrix with conkit")
    getContactMatrixWithConkit("AF-F4HVG8-F1-model_v2.pdb", "AF-F4HVG8")
    print("matrix with biopython")
    getContactMatrixWithBiopython("AF-F4HVG8-F1-model_v2.pdb", "AF-F4HVG8", 12.0)

    # to generate sample test data
    # commant the line above, and uncommant this line

    # generateTestData(20, 100, 12) # generateTestData(matrixSize, nbValToChange, contactThreshold)


if __name__ == '__main__':
    main()
from conkit.io.pdb import PdbParser
import conkit.plot
import numpy as np
import pandas as pd
import Bio.PDB
from matplotlib import pylab

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
    # for residue in residuesList:
        # print(list(residue.get_atoms()))
        # break
    # print(dir(structure))

    answer = np.zeros((len(CAList), len(CAList)), np.float)
    for row, ca1 in enumerate(CAList) :
        for col, ca2 in enumerate(CAList) :
            answer[row, col] = computeResidueDist(ca1, ca2)
    return answer
    return answer


def getPDBStruct(pdbFile, pdbCode):
    structure = Bio.PDB.PDBParser().get_structure(pdbCode, pdbFile)
    model = structure[0]
    return structure 


def map_value(value, min_value, max_value, min_result, max_result):
    result = min_result + (value - min_value)/(max_value - min_value) * (max_result - min_result)
    return result

def getContactMatrixWithBiopython(pdbFile, pdbID, contactThreshold):
    structure = getPDBStruct(pdbFile, pdbID)
    distMatrix = computeDistMatrix(structure)
    contactMatrix = distMatrix < contactThreshold
    showContactMap(distMatrix , contactMatrix )

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






################################################################################
# main function
def main():
    print("matrix with conkit")
    getContactMatrixWithConkit("AF-F4HVG8-F1-model_v2.pdb", "AF-F4HVG8")
    print("matrix with biopython")
    getContactMatrixWithBiopython("AF-F4HVG8-F1-model_v2.pdb", "AF-F4HVG8", 12.0)


if __name__ == '__main__':
    main()
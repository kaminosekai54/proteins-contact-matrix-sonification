from conkit.io.pdb import PdbParser
import conkit.plot
import numpy as np
import pandas as pd


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

p = PdbParser()
# with open("1qhw.pdb", "r") as pdb_fhandle:
with open("AF-F4HVG8-F1-model_v2.pdb", "r") as pdb_fhandle:
    pdb = p.read(pdb_fhandle, f_id="1QHW", atom_type="CA")
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
            # print(row_id, column_id)
            matrix[[row_id], column_id] = dist.raw_score
            # matrix[[dist.id[0]-1, dist.id[1]-1]] = dist.raw_score
        # print(matrix)
        np.save("distance_matrix.txt", matrix)
        df = pd.DataFrame(matrix, columns=columnId.values(), index=rowId.values())
        df.to_csv("matrix_df.csv")






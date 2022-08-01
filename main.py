from conkit.io.pdb import PdbParser
import conkit.plot
import numpy as np

p = PdbParser()
with open("AF-F4HVG8-F1-model_v2.pdb", "r") as pdb_fhandle:
    pdb = p.read(pdb_fhandle, f_id="1QHW", atom_type="CA")
    for cmap in pdb:
        print(cmap)

        fig = conkit.plot.ContactMapMatrixFigure(cmap)
        fig2 = conkit.plot.ContactMapFigure(cmap)
        fig.savefig("testMatrixMap.png", overwrite=True)
        fig2.savefig("testContactMap.png", overwrite=True)

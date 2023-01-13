import numpy as np
import openpmd_api as io
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import sys


def main(pathToShadowgramWithPath):
    
    series = io.Series(pathToShadowgramWithPath, io.Access.read_only)
    i = series.iterations[[i for i in series.iterations][0]]

    chunkdata = i.meshes["shadowgram"][io.Mesh_Record_Component.SCALAR].load_chunk()
    unit = i.meshes["shadowgram"].get_attribute("unitSI")
    series.flush()

    openpmddata = chunkdata * unit
    del chunkdata
    
    #dx = i.meshes["shadowgram"].get_attribute("gridSpacing")[0] * i.meshes["shadowgram"].get_attribute("gridUnitSI") 
    #dy = i.meshes["shadowgram"].get_attribute("gridSpacing")[1] * i.meshes["shadowgram"].get_attribute("gridUnitSI") 
    
    print(f"value of interest: {np.sum(openpmddata):.3e}")
    
    plt.imshow(openpmddata)
    plt.savefig("shadowgram.png")


if __name__ == "__main__":
    main(sys.argv[1])

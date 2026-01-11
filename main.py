from read_data import *
from MES_Classes import *
import pandas as pd
def main():
    npc = 2
    glob = GlobalData(npc)
    g = Grid()
    # filename = "Test2_4_4_MixGrid.txt"
    #filename = "Test1_4_4.txt"
    filename = "Test3_31_31_kwadrat.txt"
    if not reading_data(filename, glob, g):
        return
    deltaT = glob.simulationTimeStep
    time = deltaT
    
    T_init = glob.initialTemp
    t0 = np.full(g.nNode, T_init)

    printStructure = False #zmienic na true jak chce wypisac co wczytal (za duzo linijek przy 31x31 zeby na kazdym wypisywac)
    if printStructure:
        print("File Name ", filename)
        print(f"Number of nodes: {glob.nNode}")
        print(f"Number of elements: {g.nElem}\n")
    
        for i, elem in enumerate(g.elems, start=1):
            print(f"Element {i}:")
            print("  Node IDs:", " ".join(map(str, elem.ID)))
            print("  Node coordinates:")
            for node_id in elem.ID:
                n = g.nodes[node_id - 1]
                print(f"    Node {node_id}: ({n.x}, {n.y})")
            print()
        print("Boundary condition nodes:", glob.BC)

    while(time <= glob.simulationTime):
        obliczeniaWKroku(glob,g)
        temperatura = obliczTkroku(glob.MatrixH, glob.MatrixC, -glob.P, t0, deltaT)
        print("Temperatura w czasie t = ", time)
        print("MIN:", np.min(temperatura), "MAX", np.max(temperatura))
        t0 = temperatura
        time += deltaT


if __name__ == "__main__":
    main()

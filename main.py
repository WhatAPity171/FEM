from read_data import *
from MES_Classes import *
import pandas as pd
def main():
    data = GlobalData()
    g = Grid()

    #s = Surface(3,25 ,2)
    


    if not reading_data("Test1_4_4.txt", data, g):
        return

    #print(f"Number of nodes: {data.nNode}")
    #print(f"Number of elements: {g.nElem}\n")
    
    #for i, elem in enumerate(g.elems, start=1):
    #    print(f"Element {i}:")
    #    print("  Node IDs:", " ".join(map(str, elem.ID)))
    #    print("  Node coordinates:")
    #    for node_id in elem.ID:
    #        n = g.nodes[node_id - 1]  # adjust for 1-based indexing
    #        print(f"    Node {node_id}: ({n.x}, {n.y})")
    #    print()

    #print("Boundary condition nodes:", data.BC)
    np.set_printoptions(precision=5, suppress=True, linewidth=200)
    print(data.MatrixH)
    print(data.P)
    #zapytać o to na zaj!!!!
    t = obliczT(data.MatrixH, data.P)
    print(t)
    #print(np.allclose(data.MatrixH @ t + data.P, 0))



if __name__ == "__main__":
    main()

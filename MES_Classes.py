from GaussLegendre import *
import numpy as np
class Node:
    def __init__(self,i = 0, x=0.0, y=0.0):
        self.i = i
        self.x = x
        self.y = y
import numpy as np

class GlobalData:
    def __init__(self):
        self.simulationTime = 0
        self.simulationTimeStep = 0
        self.conductivity = 0.0
        self.alpha = 0.0
        self.Tot = 0.0
        self.initialTemp = 0.0
        self.specificHeat = 0.0
        self.density = 0.0
        self.nNode = 0
        self.nElem = 0
        self.npc = 4
        self.MatrixH = np.zeros((self.nNode,self.nNode))


class ElemUniv:
    def __init__(self):
        self.npc = 4
        self.sqrt = int(np.sqrt(self.npc))
        self.params = GaussLegendreParams(self.sqrt)

        self.points = np.array(
            [(xi, yi) for xi in self.params.points for yi in self.params.points],
            dtype=float
        )

        self.dN_dE = np.zeros((self.npc, self.npc))
        self.dN_dn = np.zeros((self.npc, self.npc))

        for i, (x, y) in enumerate(self.points):
            self.dN_dE[i, 0] = -0.25 * (1 - y)
            self.dN_dn[i, 0] = -0.25 * (1 - x)

            self.dN_dE[i, 1] = 0.25 * (1 - y)
            self.dN_dn[i, 1] = -0.25 * (1 + x)

            self.dN_dE[i, 2] = 0.25 * (1 + y)
            self.dN_dn[i, 2] = 0.25 * (1 + x)

            self.dN_dE[i, 3] = -0.25 * (1 + y)
            self.dN_dn[i, 3] = 0.25 * (1 - x)

class Grid:
    def __init__(self):
        self.nNode = 0
        self.nElem = 0
        self.nodes = []
        self.elems = []
        self.BC = []

import numpy as np

class Jakobian:
    def __init__(self):
        self.J = np.zeros((2, 2))
        self.J1 = np.zeros((2, 2))
        self.detJ = 0.0

    def obliczDet(self):
        self.detJ = np.linalg.det(self.J)

    def obliczJ1(self):
        if abs(self.detJ) < 1e-12:
            raise ValueError("Determinant Jacobian is too close to zero.")
        self.J1 = np.linalg.inv(self.J)


class Element:
    def __init__(self):
        self.npc = 4
        self.ID = []
        self.Jakobian = [Jakobian() for _ in range(self.npc)]
        self.H = np.zeros((4, 4))

    def obliczJakobiany(self, grid):
        mapping = {n.i: n for n in grid.nodes}
        nodes = [mapping[i] for i in self.ID]
        univ = ElemUniv()

        for i in range(self.npc):
            J = np.zeros((2, 2))
            for j in range(4):
                J[0, 0] += univ.dN_dE[i][j] * nodes[j].x  # dx/dksi
                J[0, 1] += univ.dN_dE[i][j] * nodes[j].y  # dy/dksi
                J[1, 0] += univ.dN_dn[i][j] * nodes[j].x  # dx/deta
                J[1, 1] += univ.dN_dn[i][j] * nodes[j].y  # dy/deta

            self.Jakobian[i].J = J
            print(f"Jakobian {i}:\n{J}\n")

            self.Jakobian[i].obliczDet()
            self.Jakobian[i].obliczJ1()

    def obliczH(self, conductivity, glob: GlobalData):
        univ = ElemUniv()
        dN_dx = np.zeros_like(univ.dN_dE)
        dN_dy = np.zeros_like(univ.dN_dn)
        for i, jak in enumerate(self.Jakobian):
            grad_local = np.vstack((univ.dN_dE[i], univ.dN_dn[i]))  # shape (2,4)
            grad_global = jak.J1 @ grad_local              # shape (2,4)
            dN_dx[i, :] = grad_global[1, :]
            dN_dy[i, :] = grad_global[0, :]
            tempY = grad_global[0, : ].reshape(-1,1)
            tempX = grad_global[1, : ].reshape(-1,1)
            self.H += (tempX @ tempX.T + tempY @ tempY.T)*jak.detJ
        #print(dN_dx)
        #print(dN_dy)
        self.H *= conductivity
        for (i,j), value in np.ndenumerate(self.H):
            globI = self.ID[i] - 1
            globJ = self.ID[j] - 1
            glob.MatrixH[globI,globJ] += value
        print(self.H)












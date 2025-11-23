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
        self.BC = []

#note to self, pozmieniac npc na 2 zamiast 4 i wczytywanie z pliku, teraz pokazuje sie w 2d zamiast w 1d i mnozenie przez siebie smh
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
class Surface:
    def __init__(self, wall, glob: GlobalData, npc = 2):#npc w formie 1d, jeśli podawane w obecnej wersji kodu trzeba podac sqrt
        self.N = np.zeros((npc,4))
        params = GaussLegendreParams(npc)
        points = []
        self.Hbc = np.zeros((4,4))
        if wall == 0: # 0 - sciana dolna
            for i in range(npc):
                point = [params.points[i], -1]
                points.append(point)
        if wall == 1: # 1 - sciana prawa
            for i in range(npc):
                point = [1, params.points[i]]
                points.append(point)        
        if wall == 2: # 2 - sciana gorna
            for i in range(npc):
                point = [params.points[i], 1]
                points.append(point)
        if wall == 3: # 3 - sciana lewa
            for i in range(npc):
                point = [-1, params.points[i]]
                points.append(point)
        #print(points)
        for i, element in enumerate(points):
            ksi = element[0]
            eta = element[1]
            self.N[i][0] = 0.25*(1-ksi)*(1-eta) #N1
            self.N[i][1] = 0.25*(1+ksi)*(1-eta) #N2
            self.N[i][2] = 0.25*(1+ksi)*(1+eta) #N3
            self.N[i][3] = 0.25*(1-ksi)*(1+eta) #N4
        for i, element in enumerate(self.N):
            tempHbc = element.reshape(-1,1) @ element.reshape(1,-1)
            #print(tempHbc)
            tempHbc *= glob.alpha
            tempHbc *= params.weights[i]
            self.Hbc += tempHbc

        
class SurfaceUniv:
    def __init__(self, npc = 2):
        self.walls = []
        for i in range(4):
            pass

class Grid:
    def __init__(self):
        self.nNode = 0
        self.nElem = 0
        self.nodes = []
        self.elems = []
        self.BC = []


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
            #print(f"Jakobian {i}:\n{J}\n")

            self.Jakobian[i].obliczDet()
            self.Jakobian[i].obliczJ1()
    def obliczHbc(self, glob: GlobalData, grid: Grid):
        Hbc = np.zeros((4,4))
        for i in range(4):
            if self.ID[i] in glob.BC and self.ID[(i+1) % 4] in glob.BC:
                wall = Surface(i,glob,int(self.npc**0.5))
                tempHbc = wall.Hbc
                x = grid.nodes[self.ID[(i+1)%4]-1].x - grid.nodes[self.ID[i]-1].x
                y = grid.nodes[self.ID[(i+1)%4]-1].y - grid.nodes[self.ID[i]-1].y
                c = x*x + y*y
                c = c**0.5
                detJ = c/2
                tempHbc *= detJ 
                Hbc += tempHbc   
        #print(Hbc)        
        return Hbc
    def obliczH(self, conductivity, glob: GlobalData, grid: Grid):
        univ = ElemUniv()
        dN_dx = np.zeros_like(univ.dN_dE)
        dN_dy = np.zeros_like(univ.dN_dn)
        for i, jak in enumerate(self.Jakobian):
            grad_local = np.vstack((univ.dN_dE[i], univ.dN_dn[i]))  # shape (2,4)
            grad_global = jak.J1 @ grad_local              # shape (2,4)
            dN_dx[i, :] = grad_global[0, :]
            dN_dy[i, :] = grad_global[1, :]
            tempY = grad_global[0, : ].reshape(-1,1)
            tempX = grad_global[1, : ].reshape(-1,1)
            self.H += (tempX @ tempX.T + tempY @ tempY.T)*jak.detJ #trzeba dodać wagi do 3 i 4 wymiarowego gaussa
        #print(dN_dx)
        #print(dN_dy)
        self.H *= conductivity
        #print("1")
        #print(self.H)
        self.H += self.obliczHbc(glob, grid)
        #print("2")
        #print(self.H)
        for (i,j), value in np.ndenumerate(self.H):
            globI = self.ID[i] - 1
            globJ = self.ID[j] - 1
            glob.MatrixH[globI,globJ] += value
        #print(self.H)












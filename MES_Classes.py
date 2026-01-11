from GaussLegendre import *
import numpy as np

class Node:
    def __init__(self,i = 0, x=0.0, y=0.0):
        self.i = i
        self.x = x
        self.y = y
import numpy as np

class GlobalData: #Zawiera wszystko co potrzebujemy, w teorii mozna tez wrzucic grid tutaj ale jest rzadko uzywany
    def __init__(self, npc):
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
        self.npc = npc
        self.MatrixH = np.zeros((self.nNode,self.nNode))
        self.MatrixC = np.zeros((self.nNode,self.nNode))
        self.BC = []
        self.P = np.zeros(1)
        self.elemUniv = ElemUniv(self.npc)

class ElemUniv: #Element uniwersalny w ukladzie lokalnym -1 1 do wyliczania funkcji ksztaltu
    def __init__(self, npc):
        self.npc = npc**2
        self.sqrt = npc
        self.params = GaussLegendreParams(self.sqrt)

        self.points = np.array(
            [(xi, yi) for xi in self.params.points for yi in self.params.points],
            dtype=float
        )
        self.N = np.zeros((self.npc,4))
        self.dN_dE = np.zeros((self.npc, 4))
        self.dN_dn = np.zeros((self.npc, 4))

        for i, (x, y) in enumerate(self.points):
            self.N[i,0] = 0.25*(1-x)*(1-y)
            self.N[i,1] = 0.25*(1+x)*(1-y)
            self.N[i,2] = 0.25*(1+x)*(1+y)
            self.N[i,3] = 0.25*(1-x)*(1+y)

            self.dN_dE[i, 0] = -0.25 * (1 - y)
            self.dN_dn[i, 0] = -0.25 * (1 - x)

            self.dN_dE[i, 1] = 0.25 * (1 - y)
            self.dN_dn[i, 1] = -0.25 * (1 + x)

            self.dN_dE[i, 2] = 0.25 * (1 + y)
            self.dN_dn[i, 2] = 0.25 * (1 + x)

            self.dN_dE[i, 3] = -0.25 * (1 + y)
            self.dN_dn[i, 3] = 0.25 * (1 - x)
class Surface: #do obliczania warunków brzegowych
    def __init__(self, wall, glob: GlobalData):
        self.npc = glob.npc
        self.N = np.zeros((self.npc,4))
        params = GaussLegendreParams(self.npc)
        points = []
        self.Hbc = np.zeros((4,4))
        self.P = 0.0
        if wall == 0: # 0 - sciana dolna
            for i in range(self.npc):
                point = [params.points[i], -1]
                points.append(point)
        if wall == 1: # 1 - sciana prawa
            for i in range(self.npc):
                point = [1, params.points[i]]
                points.append(point)        
        if wall == 2: # 2 - sciana gorna
            for i in range(self.npc):
                point = [params.points[i], 1]
                points.append(point)
        if wall == 3: # 3 - sciana lewa
            for i in range(self.npc):
                point = [-1, params.points[i]]
                points.append(point)
        for i, element in enumerate(points):
            ksi = element[0]
            eta = element[1]
            self.N[i][0] = 0.25*(1-ksi)*(1-eta) #N1
            self.N[i][1] = 0.25*(1+ksi)*(1-eta) #N2
            self.N[i][2] = 0.25*(1+ksi)*(1+eta) #N3
            self.N[i][3] = 0.25*(1-ksi)*(1+eta) #N4
        for i, element in enumerate(self.N):
            tempHbc = element.reshape(-1,1) @ element.reshape(1,-1)
            tempHbc *= glob.alpha
            tempHbc *= params.weights[i]
            self.Hbc += tempHbc
            tempP = element * glob.alpha * glob.Tot * params.weights[i]
            self.P += tempP
        


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
            raise ValueError("DetJ za blisko zera.")
        self.J1 = np.linalg.inv(self.J)


class Element:
    def __init__(self, glob: GlobalData):
        self.npc = glob.npc**2
        self.ID = []
        self.Jakobian = [Jakobian() for _ in range(self.npc)]
        self.H = np.zeros((4, 4))
        self.P = np.zeros(4)
        self.C = np.zeros((4,4))
        self.glob = glob

    def obliczJakobiany(self, grid):
        mapping = {n.i: n for n in grid.nodes}
        nodes = [mapping[i] for i in self.ID]
        univ = self.glob.elemUniv

        for i in range(self.npc):
            J = np.zeros((2, 2))
            for j in range(4):
                J[0, 0] += univ.dN_dE[i][j] * nodes[j].x  # dx/dksi
                J[0, 1] += univ.dN_dE[i][j] * nodes[j].y  # dy/dksi
                J[1, 0] += univ.dN_dn[i][j] * nodes[j].x  # dx/deta
                J[1, 1] += univ.dN_dn[i][j] * nodes[j].y  # dy/deta

            self.Jakobian[i].J = J

            self.Jakobian[i].obliczDet()
            self.Jakobian[i].obliczJ1()

    def obliczHbc(self, glob: GlobalData, grid: Grid):
        Hbc = np.zeros((4,4))
        for i in range(4):
            if self.ID[i] in glob.BC and self.ID[(i+1) % 4] in glob.BC:
                wall = Surface(i,glob)
                tempHbc = wall.Hbc
                x = grid.nodes[self.ID[(i+1)%4]-1].x - grid.nodes[self.ID[i]-1].x
                y = grid.nodes[self.ID[(i+1)%4]-1].y - grid.nodes[self.ID[i]-1].y
                c = x*x + y*y
                c = c**0.5
                detJ = c/2
                tempHbc *= detJ 
                Hbc += tempHbc

                tempP = wall.P
                tempP *= detJ
                self.P += tempP           
        return Hbc
    
    def obliczH(self, glob: GlobalData, grid: Grid):
        conductivity = glob.conductivity
        univ = glob.elemUniv
        dN_dx = np.zeros_like(univ.dN_dE)
        dN_dy = np.zeros_like(univ.dN_dn)
        gauss = GaussLegendreParams(glob.npc)
        weights = np.outer(gauss.weights, gauss.weights).ravel()

        for i, jak in enumerate(self.Jakobian):
            grad_local = np.vstack((univ.dN_dE[i], univ.dN_dn[i])) 
            grad_global = jak.J1 @ grad_local              
            dN_dx[i, :] = grad_global[0, :]
            dN_dy[i, :] = grad_global[1, :]
            tempY = grad_global[0, : ].reshape(-1,1)
            tempX = grad_global[1, : ].reshape(-1,1)
            B = tempX @ tempX.T + tempY @ tempY.T
            self.H += B*jak.detJ*weights[i]*conductivity
            
            N = univ.N[i].reshape(4,1)     # N1..N4 w punkcie Gaussa
            NtN = N @ N.T                 # (4×4)
            self.C += glob.density * glob.specificHeat * NtN * jak.detJ * weights[i]


        self.H += self.obliczHbc(glob, grid)

        for (i,j), value in np.ndenumerate(self.H):
            globI = self.ID[i] - 1
            globJ = self.ID[j] - 1
            glob.MatrixH[globI,globJ] += value
        for (i,j), value in np.ndenumerate(self.C):
            globI = self.ID[i] - 1
            globJ = self.ID[j] - 1
            glob.MatrixC[globI,globJ] += value

    def dodajP(self, glob: GlobalData):
        for i in range(4):
            glob.P[self.ID[i]-1] += self.P[i]
def obliczT(matrixH, vectorP):
    return -np.linalg.solve(matrixH, -vectorP)
    #całość jest na minusie bo w równaniu ogólnym był minus który pominęliśmy etc etc

def obliczTkroku(H,C,P,t0,dt):
    A = H + C / dt
    b = (C / dt)@t0 - P
    T1 = np.linalg.solve(A, b)
    return T1

def obliczeniaWKroku(glob: GlobalData, g: Grid):
    for elem in g.elems:
        elem.obliczJakobiany(g)
        elem.obliczH(glob,g)
        elem.dodajP(glob)








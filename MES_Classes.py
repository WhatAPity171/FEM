from GaussLegendre import *
import numpy as np
class Node:
    def __init__(self,i = 0, x=0.0, y=0.0):
        self.i = i
        self.x = x
        self.y = y
class ElemUniv:
    def __init__(self):
        self.npc = 4
        self.sqrt = int(self.npc**0.5)
        self.params = GaussLegendreParams(self.sqrt)
        self.points = [[0.0 for _ in range(self.sqrt)] for _ in range(self.sqrt)]
        for i,xi in enumerate(self.params.points):
            for j,yi in enumerate(self.params.points):
                self.points[i][j] = (xi,yi)
        self.points = [p for row in self.points for p in row]        
            
        self.dN_dE = [[0.0 for _ in range(self.npc)] for _ in range(4)]
        self.dN_dn = [[0.0 for _ in range(self.npc)] for _ in range(4)]
        rows = len(self.dN_dE)
        cols = len(self.dN_dE[0])
        for i, (x,y)  in enumerate(self.points):
            for j in range(cols):
                if j == 0:
                    self.dN_dE[i][j] = -1/4*(1 - x)
                    self.dN_dn[i][j] = -1/4*(1 - y)
                if j == 1:
                    self.dN_dE[i][j] = 1/4*(1 - x)
                    self.dN_dn[i][j] = -1/4*(1 + y)
                if j == 2:
                    self.dN_dE[i][j] = 1/4*(1 + x)
                    self.dN_dn[i][j] = 1/4*(1 + y)
                if j == 3:
                    self.dN_dE[i][j] = -1/4*(1 + x)
                    self.dN_dn[i][j] = 1/4*(1 - y)
                #print(f"{self.dN_dn[i][j]} ")
            #print("\n")        

class Grid:
    def __init__(self):
        self.nNode = 0
        self.nElem = 0
        self.nodes = []
        self.elems = []
        self.BC = []

class Jakobian:
    def __init__(self):
        self.J = [[0.0 for _ in range(2)] for _ in range(2)]
        self.J1 = [[0.0 for _ in range(2)] for _ in range(2)]
        self.detJ = 0.0
    def obliczDet(self):
        self.detJ = self.J[0][0]*self.J[1][1] - self.J[0][1]*self.J[1][0]
    def obliczJ1(self):
        self.J1[0][0] = self.J[1][1] / self.detJ
        self.J1[1][1] = self.J[0][0] / self.detJ
        self.J1[0][1] = -self.J[1][0] / self.detJ
        self.J1[1][0] = -self.J[0][1] / self.detJ
   
class Element:
    def __init__(self):
        self.npc = 4
        self.ID = []
        self.Jakobian = [Jakobian() for _ in range(self.npc)]
    def obliczJakobiany(self, grid : Grid):
        mapping = {n.i: n for n in grid.nodes}
        nodes = [mapping[i] for i in self.ID]
        univ = ElemUniv()
        for i, x in enumerate(self.Jakobian):
            for j in range(4):
                self.Jakobian[i].J[0][0] += univ.dN_dE[i][j] * nodes[j].x  #x po ksi
                self.Jakobian[i].J[1][1] += univ.dN_dn[i][j] * nodes[j].y  #y po eta
                self.Jakobian[i].J[1][0] += univ.dN_dn[i][j] * nodes[j].x  #x po eta
                self.Jakobian[i].J[0][1] += univ.dN_dE[i][j] * nodes[j].y  #y po ksi
            print(f"Jakobian{i}\n{self.Jakobian[i].J[0][0]}\t{self.Jakobian[i].J[0][1]}\n{self.Jakobian[i].J[1][0]}\t{self.Jakobian[i].J[1][1]}\n")
            self.Jakobian[i].obliczDet()
            self.Jakobian[i].obliczJ1()
            print("\n")


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









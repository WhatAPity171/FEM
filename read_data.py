from MES_Classes import *

def reading_data(filename, data, g):
    try:
        with open(filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: Cannot open file {filename}")
        return False
    i : int
    i = 0
    while i < len(lines):
        line = lines[i]

        if line.startswith("SimulationTime"):
            data.simulationTime = int(line.split()[1])
        elif line.startswith("SimulationStepTime"):
            data.simulationTimeStep = int(line.split()[1])
        elif line.startswith("Conductivity"):
            data.conductivity = float(line.split()[1])
        elif line.startswith("Alfa"):
            data.alpha = float(line.split()[1])
        elif line.startswith("Tot"):
            data.Tot = float(line.split()[1])
        elif line.startswith("InitialTemp"):
            data.initialTemp = float(line.split()[1])
        elif line.startswith("SpecificHeat"):
            data.specificHeat = float(line.split()[1])
        elif line.startswith("Density"):
            data.density = float(line.split()[1])
        elif line.startswith("Nodes number"):
            data.nNode = int(line.split()[2])
        elif line.startswith("Elements number"):
            data.nElem = int(line.split()[2])
        elif line.startswith("*Node"):
            g.nodes = []
            for j in range(data.nNode):
                i = i + 1
                parts = [p.strip() for p in lines[i].split(",")]
                id, x, y = parts
                g.nodes.append(Node(int(id), float(x), float(y)))
        elif line.startswith("*Element"):
            g.elems = []
            for j in range(data.nElem):
                i = i + 1
                parts = [p.strip() for p in lines[i].split(",")]
                # example: 1, 1, 2, 6, 5
                _, n1, n2, n3, n4 = parts
                elem = Element()
                elem.ID = [int(n1), int(n2), int(n3), int(n4)]
                elem.obliczJakobiany(g)
                g.elems.append(elem)
        elif line.startswith("*BC"):
            i = i +1
            parts = [p.strip() for p in lines[i].split(",")]
            g.BC = [int(p) for p in parts]

        i = i + 1

    g.nNode = data.nNode
    g.nElem = data.nElem
    return True
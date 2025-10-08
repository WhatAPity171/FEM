class Node:
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y


class Element:
    def __init__(self):
        self.ID = []


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


class Grid:
    def __init__(self):
        self.nNode = 0
        self.nElem = 0
        self.nodes = []
        self.elems = []
        self.BC = []


def read_data(filename, data, g):
    try:
        with open(filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: Cannot open file {filename}")
        return False

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
                i += 1
                parts = [p.strip() for p in lines[i].split(",")]
                _, x, y = parts
                g.nodes.append(Node(float(x), float(y)))
        elif line.startswith("*Element"):
            g.elems = []
            for j in range(data.nElem):
                i += 1
                parts = [p.strip() for p in lines[i].split(",")]
                # example: 1, 1, 2, 6, 5
                _, n1, n2, n3, n4 = parts
                elem = Element()
                elem.ID = [int(n1), int(n2), int(n3), int(n4)]
                g.elems.append(elem)
        elif line.startswith("*BC"):
            i += 1
            parts = [p.strip() for p in lines[i].split(",")]
            g.BC = [int(p) for p in parts]

        i += 1

    g.nNode = data.nNode
    g.nElem = data.nElem
    return True


def main():
    data = GlobalData()
    g = Grid()

    if not read_data("Test1_4_4.txt", data, g):
        return

    print("File loaded successfully!")
    print(f"Number of nodes: {g.nNode}")
    print(f"Number of elements: {g.nElem}\n")

    for i, elem in enumerate(g.elems, start=1):
        print(f"Element {i}:")
        print("  Node IDs:", " ".join(map(str, elem.ID)))
        print("  Node coordinates:")
        for node_id in elem.ID:
            n = g.nodes[node_id - 1]  # adjust for 1-based indexing
            print(f"    Node {node_id}: ({n.x}, {n.y})")
        print()

    print("Boundary condition nodes:", g.BC)


if __name__ == "__main__":
    main()

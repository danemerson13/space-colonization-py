class Branch:
    # Properties of a branch
    # proximal: type = Node
    # distal: type = Node
    # parent: type = Branch
    # children: type = list(Branch)
    # radius: type = double
    # resistance: type = double
    # flowrate: type = double
    # generation: type = int

    def __init__(self, prox, dist):
        self.proximal = prox
        self.distal = dist
        self.parent = None
        self.children = list()

        self.radius = None
        self.resistance = None
        self.flowrate = None
        self.generation = None

    def getProximal(self):
        return self.proximal

    def setDistal(self, dist):
        self.distal = dist

    def getDistal(self):
        return self.distal

    def setParent(self, parent):
        self.parent = parent

    def getParent(self):
        return self.parent

    def setChild(self, child):
        self.children.append(child)

    def getChildren(self):
        return self.children

    def removeChild(self, child):
        self.children.remove(child)

    def isRoot(self):
        return self.getProximal().isRoot()

    def isTerminal(self):
        return self.getDistal().isTerminal()

    def setRadius(self, rad):
        self.radius = rad
    
    def getRadius(self):
        return self.radius

    def setResistance(self, R):
        self.resistance = R

    def getResistance(self):
        return self.resistance

    def setFlowRate(self, Q):
        self.flowrate = Q

    def getFlowRate(self):
        return self.flowrate

    def setGeneration(self, gen):
        self.generation = gen
    
    def getGeneration(self):
        return self.generation
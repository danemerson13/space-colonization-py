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
        self.volume = None
        self.fluid = 0
        self.buffer = 0

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

    def removeParent(self):
        self.parent = None

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

    def setVolume(self, V):
        self.volume = V

    def getVolume(self):
        return self.volume

    def addFluid(self, V):
        if self.fluid + V > self.getVolume():
            self.setBuffer(self.fluid + V - self.getVolume())
            self.fluid = self.volume
        else:
            self.fluid += V

    def getFluid(self):
        return self.fluid

    def percentFull(self):
        return self.getFluid()/self.getVolume()

    def setBuffer(self, V):
        self.buffer = V

    def getBuffer(self):
        return self.buffer
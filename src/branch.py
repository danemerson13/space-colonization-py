class Branch:
    # Properties of a branch
    # proximal: type = Node
    # distal: type = Node
    # parents: type = list(Branch)
    # children: type = list(Branch)
    # type: type = string. "Root", "Terminal", "Inlet", "Outlet", None
    # radius: type = double
    # resistance: type = double
    # flowrate: type = double
    # volume: type = double
    # fluid: type = double
    # buffer: type = double


    def __init__(self, prox, dist):
        self.proximal = prox
        self.distal = dist
        self.parents = list()
        self.children = list()
        self.type = None

        self.setType()

        self.length = None
        self.radius = None
        self.resistance = None
        self.flowrate = None
        self.volume = None
        self.concentration = 0
        self.updateFlag = None

    def setProximal(self, prox):
        self.proximal = prox
        self.setType()

    def getProximal(self):
        return self.proximal

    def setDistal(self, dist):
        self.distal = dist
        self.setType()

    def getDistal(self):
        return self.distal

    def setParent(self, parent):
        self.parents.append(parent)

    def getParents(self):
        return self.parents

    def removeParent(self, parent):
        self.parents.remove(parent)

    def setChild(self, child):
        self.children.append(child)

    def getChildren(self):
        return self.children

    def removeChild(self, child):
        self.children.remove(child)

    def setType(self):
        if self.getProximal().getType() == "Inlet" or self.getDistal().getType() == "Inlet":
            self.type = "Inlet"
        elif self.getProximal().getType() == "Outlet" or self.getDistal().getType() == "Outlet":
            self.type = "Outlet"
        else:
            self.type = None

    def getType(self):
        return self.type

    def setLength(self, L):
        self.length = L

    def getLength(self):
        return self.length

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
    
    def getConcentration(self):
        return self.concentration

    def updateConcentration(self, Cin, Vin):
        Vol = self.getVolume()
        if Vin < Vol:
            self.concentration = (Cin * Vin + (Vol - Vin) * self.getConcentration())/Vol
        else:
            self.concentration = Cin
        self.setUpdateFlag(True)

    def setUpdateFlag(self, val):
        self.updateFlag = val

    def isUpdated(self):
        return self.updateFlag

##### UTILITY FUNCTIONS #####

    def isRoot(self):
        if self.getProximal().getType() == "Root" or self.getDistal().getType() == "Root":
            return True
        else:
            return False

    def isTerminal(self):
        if self.getProximal().getType() == "Terminal" or self.getDistal().getType() == "Terminal":
            return True
        else:
            return False

    def isSL(self):
        if self.getProximal().isSL() or self.getDistal().isSL():
            return True
        else:
            return False

    def percentFull(self):
        return self.getFluid()/self.getVolume()
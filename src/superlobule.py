class SuperLobule:
    # Properties of a super lobule:
    # location: type = 3x1 double array
    # prox: type = Node
    # dist: type = Node
    # parents: type = list(branch) (from inlet tree)
    # children: type = list(branch) (from outlet tree)
    # resistance: type = double
    # flowrate: type = double
    
    def __init__(self, loc, prox, dist, parent, child, vol = 0, R = 0):
        self.location = loc
        self.proximal = prox
        self.distal = dist
        self.parents = list([parent])
        self.children = list([child])
        self.resistance = R
        self.flowrate = None
        self.volume = vol
        self.concentration = 0
        self.concentrationList = list([self.concentration])
        self.updateFlag = None

    def getLocation(self):
        return self.location

    def getProximal(self):
        return self.proximal

    def getDistal(self):
        return self.distal
    
    def getParents(self):
        return self.parents
    
    def getChildren(self):
        return self.children

    def setResistance(self, R):
        self.resistance = R

    def getResistance(self):
        return self.resistance

    def setFlowRate(self, Q):
        self.flowrate = Q

    def getFlowRate(self):
        return self.flowrate

    def getVolume(self):
        return self.volume
    
    def setVolume(self, vol):
        self.volume = vol

    def getConcentration(self):
        return self.concentration

    def updateConcentration(self, Cin, Vin):
        Vol = self.getVolume()
        if Vin < Vol:
            self.concentration = (Cin * Vin + (Vol - Vin) * self.getConcentration())/Vol
        else:
            self.concentration = Cin
        # Append to concentrationList
        self.concentrationList.append(self.concentration)
        self.setUpdateFlag(True)

    def setUpdateFlag(self, val):
        self.updateFlag = val

    def isUpdated(self):
        return self.updateFlag
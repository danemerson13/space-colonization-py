class SuperLobule:
    # Properties of a super lobule:
    # location: type = 3x1 double array
    # parent: type = Node (from inlet tree)
    # child: type = Node (from outlet tree)
    # resistance: type = double
    # flowrate: type = double
    
    def __init__(self, loc, prox, dist, vol = 0, R = 0):
        self.location = loc
        self.proximal = prox
        self.distal = dist
        self.resistance = R
        self.flowrate = None
        self.volume = vol

    def getLocation(self):
        return self.location

    def getProximal(self):
        return self.proximal

    def getDistal(self):
        return self.distal

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
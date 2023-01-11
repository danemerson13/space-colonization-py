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
        self.fluid = 0
        self.buffer = 0

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

    def addFluid(self, V):
        if self.fluid + V > self.getVolume():
            self.setBuffer(self.fluid + V - self.getVolume())
            self.fluid = self.volume
        else:
            self.fluid += V
    
    def getFluid(self):
        return self.fluid

    def setBuffer(self, V):
        self.buffer = V

    def getBuffer(self):
        return self.buffer
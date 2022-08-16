class SuperLobule:
    # Properties of a super lobule:
    # location: type = 3x1 double array
    # parent: type = Node (from inlet tree)
    # child: type = Node (from outlet tree)
    # resistance: type = double
    # current: type = double
    
    def __init__(self, loc, parent, child, R = 0):
        self.location = loc
        self.parent = parent
        self.child = child
        self.resistance = R
        self.flowrate = None

    def getLocation(self):
        return self.location

    def getParent(self):
        return self.parent

    def getChild(self):
        return self.child

    def setResistance(self, R):
        self.resistance = R

    def getResistance(self):
        return self.resistance

    def setFlowRate(self, Q):
        self.flowrate = Q

    def getFlowRate(self):
        return self.flowrate
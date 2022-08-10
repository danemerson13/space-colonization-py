class Node:
    # Properties of a node:
    # location: type = 3x1 double array
    # type: type = {"Inlet", "Outlet"}
    # parents: type = list(Node). If none, then root
    # children: type = list(Node). If none - terminal, if more than one - furcation
    # pressure: type = double
    # current: type = double

    def __init__(self, loc, type, parent = None, terminal = False):
        self.location = loc
        self.type = type
        self.terminal = terminal

        self.children = list()
        self.parents = list()
        if parent != None:
            self.parents.append(parent)

        self.pressure = None
        self.flowrate = None

    def getLocation(self):
        return self.location

    def setParent(self, parent):
        self.parents.append(parent)

    def getParents(self):
        return self.parents
    
    def setChild(self, child):
        self.children.append(child)

    def getChildren(self):
        return self.children

    def removeChild(self, child):
        self.children.remove(child)

    def isRoot(self):
        if len(self.parents) < 1:
            return True
        else:
            return False

    def setTerminal(self, flag):
        self.terminal = flag

    def isTerminal(self):
        return self.terminal

    def isFurcation(self):
        if self.type == "Inlet":
            if len(self.children) > 1:
                return True
            else:
                return False
        elif self.type == "Outlet":
            if len(self.parents) > 1:
                return True
            else:
                return False
        else:
            print("Invalid node type")
            return None

    def setPressure(self, P):
        self.pressure = P

    def getPressure(self):
        return self.pressure

    def setFlowRate(self, Q):
        self.flowrate = Q

    def getFlowRate(self):
        return self.flowrate

    
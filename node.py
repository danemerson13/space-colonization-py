class Node:
    # Properties of a node:
    # location: type = 3x1 double array
    # parent: type = Node. If none, then root
    # children: type = list(Node). If none - terminal, if more than one - furcation
    # pressure: type = double
    # current: type = double

    def __init__(self, loc, parent = None, terminal = False):
        self.location = loc
        self.parent = parent
        self.terminal = terminal
        self.children = list()
        self.pressure = None
        self.current = None

    def getLocation(self):
        return self.location

    def getParent(self):
        return self.parent

    def setChild(self, child):
        self.children.append(child)

    def getChildren(self):
        return self.children

    def removeChild(self, child):
        self.children.remove(child)

    def isRoot(self):
        if self.parent == None:
            return True
        else:
            return False

    def setTerminal(self, flag):
        self.terminal = flag

    def isTerminal(self):
        return self.terminal

    def isFurcation(self):
        if len(self.children) > 1:
            return True
        else:
            return False

    def setPressure(self, P):
        self.pressure = P

    def getPressure(self):
        return self.pressure

    def setCurrent(self, I):
        self.current = I

    def getCurrent(self):
        return self.current
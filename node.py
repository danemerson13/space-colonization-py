class Node:
    # Properties of a node:
    # location: type = 3x1 double array
    # parents: type = list(Node).
    # children: type = list(Node).
    # type: type = string. "Interior", "Inlet", "Outlet", "Root", "Terminal", None
    # sl: type = boolean. True if superlobule exists at this node
    # pressure: type = double
    # current: type = double

    def __init__(self, loc):
        self.location = loc
        self.parents = list()
        self.children = list()
        self.type = None
        self.sl = False
        self.pressure = None
        self.current = None

    def getLocation(self):
        return self.location

    def setParent(self, parent):
        self.parents.append(parent)
        self.setType()

    def getParents(self):
        return self.parents

    def removeParent(self, parent):
        self.parents.remove(parent)
        self.setType()

    def setChild(self, child):
        self.children.append(child)
        self.setType()
            
    def getChildren(self):
        return self.children

    def removeChild(self, child):
        self.children.remove(child)
        self.setType()

    def setType(self):
        if len(self.parents) == 1 and len(self.children) == 1:
            self.type = "Interior"
        elif len(self.parents) == 1 and len(self.children) > 1:
            self.type = "Inlet"
        elif len(self.parents)  > 1 and len(self.children) == 1:
            self.type = "Outlet"
        elif len(self.parents) == 0 and len(self.children) > 0:
            self.type = "Root"
        elif len(self.parents) > 0 and len(self.children) == 0:
            self.type = "Terminal"
        else:
            self.type = None

    def getType(self):
        return self.type

    def setSL(self, bool):
        self.sl = bool
    
    def isSL(self):
        return self.sl

    def setPressure(self, P):
        self.pressure = P

    def getPressure(self):
        return self.pressure

    def setCurrent(self, I):
        self.current = I

    def getCurrent(self):
        return self.current

##### UTILITY FUNCTIONS #####

    def isRoot(self):
        return self.type == "Root"

    def isFurcation(self):
        return self.type == "Terminal"

    
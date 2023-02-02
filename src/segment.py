import numpy as np
from src import branch, superlobule

class Segment:
    # Properties of a segment
    # proximal: type = 3x1 np.array
    # distal: type = 3x1 np.array
    # parents: type = list(Segment)
    # children: type = list(Segment)
    # ancestor: type = either branch or superlobule class
    # root: type = bool 
    # SL: type = bool
    def __init__(self, prox, dist, ancestor):
        self.proximal = prox
        self.distal = dist
        self.parents = list()
        self.children = list()
        self.ancestor = ancestor
        self.root = False
        self.SL = False

        self.volume = None
        self.setVolume()
        self.concentration = 0

    def setProximal(self, prox):
        self.proximal = prox

    def getProximal(self):
        return self.proximal

    def setDistal(self, dist):
        self.distal = dist

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

    def getAncestor(self):
        return self.ancestor
    
    def getRadius(self):
        if self.ancestor == branch.Branch:
            return self.ancestor.getRadius()
        else: # sl ancestor, compute rad from eq of sphere volume
            return np.power((3/4*self.getVolume()/np.pi),1/3)

    def getFlowRate(self):
        return self.ancestor.getFlowRate()

    def setVolume(self):
        if self.ancestor == branch.Branch:
            self.volume =  np.pi * self.getRadius()**2 * np.linalg.norm(self.getProximal() - self.getDistal())
        else:
            self.volume = self.ancestor.getVolume()

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

    def setRoot(self):
        self.root = True

    def isRoot(self):
        return self.root

    def setSL(self):
        self.SL = True

    def isSL(self):
        return self.SL

    def getType(self):
        return type(self.ancestor)

    def isSLAdjacent(self):
        # Return True if the segment lies adjacent to an SL or is an SL
        # If it does we dont want to connect it to the other segment 
        # that shares a prox/dist node, only want to connect it to the SL
        if self.isSL():
            return True
        for child in self.children:
            if child.isSL():
                return True
        for parent in self.parents:
            if parent.isSL():
                return True
        # Otherwise
        return False
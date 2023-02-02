import numpy as np
from src import branch, superlobule

class Segment:
    # Properties of a segment
    # proximal: type = 3x1 np.array
    # distal: type = 3x1 np.array
    # parents: type = list(Segment)
    # children: type = list(Segment)
    # ancestor: type = either branch or superlobule class
    def __init__(self, prox, dist, ancestor):
        self.proximal = prox
        self.distal = dist
        self.parents = list()
        self.children = list()
        self.ancestor = ancestor
        self.root = False

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
    
    def getRadius(self):
        if self.ancestor == branch.Branch:
            return self.ancestor.getRadius()
        else:
            return None

    def getFlowRate(self):
        return self.ancestor.getFlowRate

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
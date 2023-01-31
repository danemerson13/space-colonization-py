from src import attractor, node, branch, superlobule, segment
import numpy as np
from util import plotter
import imageio, os
import time
import pickle
from matplotlib import pyplot as plt

class Colony:
    def __init__(self, D, dk, di, ):
        
        self.D = D
        self.dk = dk
        self.di = di

        self.attractorList = list()
        self.nodeList = list()
        self.slList = list()
        self.branchList = list()
        self.segList = list()

        self.fillTime = None
        self.overfillTime = None

##### HIGH LEVEL FUNCTIONS #####

    def initTree(self, targets, attractors, root):
        pass

    def createTree(self, targets, attractors, initNodeList, Rinitial, mu):
        pass

    def mergeTrees(self, secondTree):
        pass

    def solveResistanceNetwork(self, Pin, Pout, ):
        pass

    def solveRSL(self, Qactual, Pin, Pout, iterations, tol, verbosity, a = 0, b = 1e-4):
        pass

##### PRIMARY FUNCTIONS #####

    def runSpaceColonization(self):
        pass

    def reverseTree(self):
        pass

    def fillTree(self, dt):
        # Note do so with segment methodology
        pass

    def saveModel(self):
        pass

##### HELPER FUNCTIONS #####

    def addAttractors(self, targets, attractors):
        pass

    def getLocationArray(self):
        pass

    def getAttractorArray(self):
        pass

    def getNodeArray(self):
        pass

    def getNumTargets(self):
        pass

    def computeDistanceVec(self, arr, pt):
        pass

    def findClosestNode(self):
        pass

    def attractorNormal(self, attractorIdx, nodePtr):
        pass

    def repeatedStepCheck(self, loc, growthNode):
        pass

    def addNode(self, growthNode, dir, activeAttractors):
        pass

    def addSLNode(self, growthNode, dir, D_exact):
        pass

    def checkTargetDistance(self, attractorIdx, nodePtr):
        pass

    def killAttractors(self):
        pass

    def trimNodes(self):
        pass

    def removeNode(self, node):
        pass

    def createBranches(self):
        pass

    def traverseBranch(self, endNode):
        pass

    def findParentBranch(self, childBranch):
        pass

    def findCorrespondingBranch(self, parentBranch):
        # formerly findChildBranchCrossTree
        pass

    def countTerminals(self):
        pass

    def setRadii(self, initial):
        pass

    def branchLength(self, branch):
        pass

    def setResistance(self):
        pass

    def setBranchVolume(self):
        pass

    def setRSL(self, RSL):
        pass

    def queryQin(self, Pin, Pout, RSL):
        pass

    def setSLVolume(self, bloodVolume):
        totalSLVolume = bloodVolume - self.totalVolume()
        SLVolume = totalSLVolume / len(self.slList)
        print(SLVolume)
        for sl in self.slList:
            sl.setVolume(SLVolume)

    def totalVolume(self):
        pass

    def getBranchNodes(self, branch):
        pass

    def createSegments(self, lmax):
        pass

    def connectSegments(self):
        pass

    def getTotalConcentration(self):
        pass

    def findRootSegment(self):
        pass

    

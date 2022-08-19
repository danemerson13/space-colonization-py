import attractor, node, branch, superlobule
import numpy as np
import plotter
import imageio, os

class Colony:
    def __init__(self, D, dk, di, ):
        
        self.D = D
        self.dk = dk
        self.di = di

        self.attractorList = list()
        self.nodeList = list()
        self.slList = list()
        self.branchList = list()

##### HIGH LEVEL FUNCTIONS #####

    def initTree(self, targets, attractors, root):
        # Add the attractors and target attractors
        self.addAttractors(targets, attractors)
        # Add the root node
        self.nodeList.append(node.Node(root))

        # Run SC Algorithm
        self.runSpaceColonization(initflag = True)

        # Trim the node list
        self.trimNodes()

    def createTree(self, targets, attractors, initNodeList, Rinitial, mu):
        # Add the attractors and target attractors
        self.addAttractors(targets, attractors)
        # Add the initialized nodes
        for inNode in initNodeList:
            self.nodeList.append(inNode)
            self.nodeList[-1].setSL(False)

        # Run SC Algorithm
        self.runSpaceColonization(initflag = False)

        # Trim nodes and create branches
        self.trimNodes()
        self.createBranches()

        # Assign radii and resistances
        self.setRadii(Rinitial)
        self.setResistance(mu)

    def mergeTrees(self, outTree):
        # First reverse the outlet tree
        outTree = outTree.reverseTree()

        # Find matching SL nodes between trees and place a SL to join them
        for branch in self.branchList:
            if branch.isTerminal():
                self.slList.append(superlobule.SuperLobule(branch.getDistal().getLocation(), branch.getDistal(), outTree.findNodeByLoc(branch.getDistal().getLocation())))

        # Add the nodes and branches from the outlet tree
        for node in outTree.nodeList:
            self.nodeList.append(node)
        for branch in outTree.branchList:
            self.branchList.append(branch)

    def solveResistanceNetwork(self, Pin, Pout):
        # Build square A matrix to solve linear system Ax = b
        # x is the all of the branchwise flow rates and nodal pressures, concatenated into one vector
        # A will be comprised of three types of equations
        # (1) Branch eqns: Vdist - Vprox - IR = 0 (Ohms law)
        # (2) Nodal eqns: Iin - sum(Iout) = 0 (Kirchoffs current law)
        # (3) BC eqns: V = const (Direct assignment)

        # Need to assign boundary conditions and count number of nodes, BCs
        nBranch = len(self.branchList)
        nSL = len(self.slList)
        nNode = 0
        nBC = 0
        VNodeList = list()
        for node in self.nodeList:
            if node.isSL():
                VNodeList.append(node)
                nNode += 1
            else:
                if node.getType() == "Root":
                    node.setPressure(Pin)
                    VNodeList.append(node)
                    nNode += 1
                    nBC += 1
                elif node.getType() == "Terminal":
                    node.setPressure(Pout)
                    VNodeList.append(node)
                    nNode += 1
                    nBC += 1
                elif node.getType() == "Inlet" or node.getType() == "Outlet": # Furcation node
                    VNodeList.append(node)
                    nNode += 1

        A = np.zeros((nBranch+nNode+nSL, nBranch+nNode+nSL))
        b = np.zeros(nBranch+nNode+nSL)
        # x = [branchList slList, nodeList]

        # (1.1) Branch eqns: vdist - vprox - IR = 0 (Ohms law)
        for i, branch in enumerate(self.branchList):
            A[i,i] = -branch.getResistance()
            A[i, VNodeList.index(branch.getProximal()) + nBranch + nSL] = 1
            A[i, VNodeList.index(branch.getDistal()) + nBranch + nSL] = -1

        # (1.2) SL eqns: vdist - vprox - IR = 0 (Ohms law)
        j = nBranch
        for sl in self.slList:
            A[j,j] = -sl.getResistance()
            A[j, VNodeList.index(sl.getProximal()) + nBranch + nSL] = 1
            A[j, VNodeList.index(sl.getDistal()) + nBranch + nSL] = -1
            j += 1

        # (2) Nodal eqns: Iin - Iout1 - Iout2 = 0 (Kirchoffs current law)
        k = nBranch + nSL
        for branch in self.branchList:
            if branch.isSL():
                if branch.getType() == "Inlet":
                    A[j, self.branchList.index(branch)] = 1
                    # Find the matching SL
                    sl = self.findSLByLoc(branch.getDistal().getLocation())
                    A[j, nBranch + self.slList.index(sl)] = -1
                else: # branch.getType() == "Outlet"
                    A[j, self.branchList.index(branch)] = -1
                    # Find the matching SL
                    sl = self.findSLByLoc(branch.getProximal().getLocation())
                    A[j, nBranch + self.slList.index(sl)] = 1
                j += 1
            else: # branch not SL
                if branch.getType() == "Inlet":
                    A[j, self.branchList.index(branch)] = 1
                    for child in branch.getChildren():
                        A[j, self.branchList.index(child)] = -1
                else: # branch.getType() == "Outlet"
                    A[j, self.branchList.index(branch)] = -1
                    for parent in branch.getParents():
                        A[j, self.branchList.index(parent)] = 1
                j += 1

        # (3) BC eqns
        k = len(A) - nBC
        for node in VNodeList:
            if node.getType() == "Root" or node.getType() == "Terminal":
                if not node.isSL():
                    A[k, VNodeList.index(node) + nBranch + nSL] = 1
                    b[k] = node.getPressure()
                    k += 1

        # Solve Ax = b
        x = np.linalg.solve(A, b)

        # Assign computed I and V values
        for i, branch in enumerate(self.branchList):
            branch.setFlowRate(x[i])

        j = nBranch
        for sl in self.slList:
            sl.setFlowRate(x[j])
            j += 1

        k = nBranch + nSL
        for node in VNodeList:
            node.setPressure(x[k])
            k += 1

##### PRIMARY FUNCTIONS #####

    def runSpaceColonization(self, initflag):
        prevNodeCount = len(self.nodeList)
        while self.numTargets() > 0:
            # First find the node closest to each attractor
            # Note: this is returned as a list with pointers to the closest nodes, 
            # with None returned where there are no nodes within di
            influencers = self.findClosestNode(initflag)
            # Take the unique entries of the influencer list
            growthNodes = list(set(influencers))
            # If there is a None entry, remove it
            if None in growthNodes:
                growthNodes.remove(None)
            # Now step through all of the nodes that the tree will be grown from
            for growthNode in growthNodes:
                # Get indices of attractors for the specific growthNode
                activeAttractors = [i for i, node in enumerate(influencers) if node == growthNode]
                # Check for case where there is only on attractor node, it is a target attractor, and within D distance
                if len(activeAttractors) == 1 and self.attractorList[activeAttractors[0]].isTarget() and self.checkTargetDistance(activeAttractors[0], growthNode) < self.D:
                    # Compute the direction to step in
                    dir = self.attractorNormal(activeAttractors, growthNode)
                    # What is the distance < D that we need to step to exactly reach the target
                    D_exact = self.checkTargetDistance(activeAttractors[0], growthNode)
                    # Add a *terminal* node in this direction at distance D_exact
                    self.addSLNode(dir, D_exact, growthNode)
                    # Flag this target attractor as killed
                    self.attractorList[activeAttractors[0]].setKill()
                else:
                    # Compute direction to step in
                    dir = self.attractorNormal(activeAttractors, growthNode)
                    # Add node in this direction at distance D
                    self.addNode(dir, growthNode, activeAttractors)
            # Once all of the nodes have been added, we will kill any attractors within dk of
            # existing nodes. Note: we cannot kill target attractors, unless they have been reached
            self.killAttractors()
            # Check to see if no steps were taken
            nodeCount = len(self.nodeList)
            if nodeCount == prevNodeCount:
                output = 'Broken with ' + str(self.numTargets()) + ' targets not reached.'
                print(output)
                break
            prevNodeCount = nodeCount

    def reverseTree(self):
        new = Colony(self.D, self.dk, self.di)
        # Copy the nodes and branches over
        for oldNode in self.nodeList:
            new.nodeList.append(node.Node(oldNode.getLocation()))
            new.nodeList[-1].setSL(oldNode.isSL())
        for oldBranch in self.branchList:
            # Make sure to assign the (reversed) prox/dist nodes from the new tree
            new.branchList.append(branch.Branch(new.findNodeByLoc(oldBranch.getDistal().getLocation()), new.findNodeByLoc(oldBranch.getProximal().getLocation())))
            # Carry over the branch radius, and resistance
            new.branchList[-1].setRadius(oldBranch.getRadius())
            new.branchList[-1].setResistance(oldBranch.getResistance())
        # Add parent child relationships for the new nodes
        for i, oldNode in enumerate(self.nodeList):
            for parent in oldNode.getParents():
                new.nodeList[i].setChild(new.findNodeByLoc(parent.getLocation()))
            for child in oldNode.getChildren():
                new.nodeList[i].setParent(new.findNodeByLoc(child.getLocation()))
        # Make sure to set the branch type again after the nodes have had their parents and children assigned
        for bran in new.branchList:
            bran.setType()
        # Add parent child relationships for new branches
        for i, oldBranch in enumerate(self.branchList):
            for parent in oldBranch.getParents():
                new.branchList[i].setChild(new.findBranchByLoc(parent.getDistal(), parent.getProximal()))
            for child in oldBranch.getChildren():
                new.branchList[i].setParent(new.findBranchByLoc(child.getDistal(), child.getProximal()))
        # Return the reversed tree
        return new

##### HELPER FUNCTIONS #####

    def addAttractors(self, targets, attractors):
        # Clear the list just to be safe
        self.attractorList = list()
        # Add all the target attractors to the attractorList
        for i in range(len(targets)):
            self.attractorList.append(attractor.Attractor(targets[i,:], True))
        # Add all the regular attractors to the attractorList
        for i in range(len(attractors)):
            self.attractorList.append(attractor.Attractor(attractors[i,:], False))

    def getLocationArray(self, inputList):
        # Returns np array with Euclidean coordinates of input list (attractors, nodes, etc.)
        arr = np.zeros((len(inputList),3))
        for i in range(len(inputList)):
            arr[i,:] = inputList[i].getLocation()
        return arr

    def getAttractorArray(self):
        # Returns np array with locations of the attractors
        return self.getLocationArray(self.attractorList)

    def getNodeArray(self, initflag):
        # Returns np array with locations of the nodes
        # Fill terminal nodes with inf since we dont want to grow from them
        # If we are growing the tree from initialized tree, we also dont want to grow from the root node (initflag serves this purpose)
        arr = self.getLocationArray(self.nodeList)
        for i in range(len(self.nodeList)):
            if self.nodeList[i].isSL():
                arr[i,:] = np.full((3,), np.inf)
            elif not initflag and self.nodeList[i].isRoot():
                arr[i,:] = np.full((3,), np.inf)
        return arr

    def numTargets(self):
        # Counts the number or target attractors in attractorList
        ct = 0
        for att in self.attractorList:
            if att.isTarget():
                ct += 1
        return ct

    def computeDistanceVec(self,arr,pt):
        # Computes the Euclidean distance between pt and every point in arr
        return np.linalg.norm(arr - pt, axis = 1)

    def findClosestNode(self, initflag):
        # For each of the attractors, return a list with a pointer to the closest node, 
        # permitting that the attractor is within the sphere of influence, di. Otherwise return None
        nodeArr = self.getNodeArray(initflag)
        closestNode = list()
        for att in self.attractorList:
            dist = self.computeDistanceVec(nodeArr, att.getLocation())
            if min(dist) < self.di:
                closestNode.append(self.nodeList[np.argmin(dist)])
            else:
                closestNode.append(None)
        return closestNode

    def attractorNormal(self, attractorIdx, nodePtr):
        attractorLocs = self.getLocationArray(self.attractorList)[attractorIdx,:]
        nodeLoc = nodePtr.getLocation()
        attractorNormals = attractorLocs - nodeLoc
        # Normalize attractor normals
        attractorNormals = (attractorNormals.T / np.linalg.norm(attractorNormals, axis = 1)).T
        # Compute mean of all unit normals
        meanNormal = np.mean(attractorNormals, axis = 0)
        # If the computed step normal is zero, randomly step towards one of the normals
        # This is edge case 1 where the algorithm would get stuck between two attractors 
        # and the computed normal would be zero and effectively no step was taken
        if np.linalg.norm(meanNormal) < 1e-6:
            return attractorNormals[np.random.choice(range(len(attractorNormals)))]
        else:
            return meanNormal/np.linalg.norm(meanNormal)

    def repeatedStepCheck(self, loc, growthNode):
        # Checks for existing children of the growthNode that have the same location as the prospective new node
        # Returns True if there is a existing child with the same location, False otherwise
        eps = 1e-6
        for child in growthNode.getChildren():
            if np.linalg.norm(loc - child.getLocation()) < eps:
                return True
        return False

    def addNode(self, dir, growthNode, activeAttractors):
        # Add a node in direction dir, D distance from node
        loc = growthNode.getLocation() + dir * self.D
        # Check for edge case 2 where the step taken is identical to previous and too far past the attractors
        if self.repeatedStepCheck(loc, growthNode):
            # print("REPEATED STEP, Num Attractors: ", len(activeAttractors))
            modifiedLoc = np.zeros((3,))
            for att in activeAttractors:
                modifiedLoc += self.attractorList[att].getLocation()
            modifiedLoc /= len(activeAttractors)
            if np.linalg.norm(growthNode.getLocation() - modifiedLoc) < self.D:
                loc = modifiedLoc
        self.nodeList.append(node.Node(loc))
        self.nodeList[-1].setParent(growthNode)
        # Also update parent node to have this new node as its child
        growthNode.setChild(self.nodeList[-1])

    def addSLNode(self, dir, D_exact, growthNode):
        # Add a SL node in direction dir, D_exact distance from node
        loc = growthNode.getLocation() + dir * D_exact
        self.nodeList.append(node.Node(loc))
        self.nodeList[-1].setParent(growthNode)
        self.nodeList[-1].setSL(True)
        # Also update parent node to have this new node as its child
        growthNode.setChild(self.nodeList[-1])

    def checkTargetDistance(self, attractorIdx, nodePtr):
        attractorLoc = self.attractorList[attractorIdx].getLocation()
        nodeLoc = nodePtr.getLocation()
        return np.linalg.norm(attractorLoc - nodeLoc)

    def killAttractors(self):
        nodeLocs = self.getLocationArray(self.nodeList)
        # If you remove items of a list while iterating through we will mess things up, instead keep a list
        killList = list()
        for att in self.attractorList:
            # First check for killed target attractors
            if att.isKilled():
                killList.append(att)
            # Now check for attractors within dk of any nodes
            elif min(self.computeDistanceVec(nodeLocs, att.getLocation())) < self.dk:
                if not att.isTarget():
                    killList.append(att)
        for att in killList:
            self.attractorList.remove(att)

    def trimNodes(self):
        criticalNodes = list()
        for curNode in self.nodeList:
            if curNode.isSL():
                criticalNodes.append(curNode)
                while curNode.getParents():
                    curNode = curNode.getParents()[0]
                    criticalNodes.append(curNode)
        # Remove the nodes not in criticalNodes list
        for node in reversed(self.nodeList):
            if criticalNodes.count(node) < 1:
                self.removeNode(node)

    def removeNode(self, node):
        # Remove the node
        self.nodeList.remove(node)
        # Also remove node from parents child list
        node.getParents()[0].removeChild(node)

    def createBranches(self):
        # Will create branches by taking all of the distal nodes in the tree
        # These will either be leaf nodes, or furcation ndoes
        distList = list()
        for node in self.nodeList:
            if node.getType() == "Terminal" or node.getType() == "Inlet":
                distList.append(node)
        # From the distal list we can traverse each branch backwards until we reach 
        # a node with multiple children, or no parent in the case of the root branch
        for node in distList:
            prox, dist = self.traverseBranch(node)
            self.branchList.append(branch.Branch(prox, dist))
        # Need to assign parent - child relationships between branches
        for bran in self.branchList:
            parent = self.findParentBranch(bran)
            # Assign to both
            if parent != None:
                bran.setParent(parent)
                parent.setChild(bran)
            elif not bran.isRoot():
                raise RuntimeError("Cannot find parent branch for non-root branch")

    def traverseBranch(self, endNode):
        dist = endNode
        curNode = endNode.getParents()[0]
        while curNode.getType() == "Interior":
            curNode = curNode.getParents()[0]
        prox = curNode
        return prox, dist

    def findParentBranch(self, childBranch):
        # Find the branch in branchList where the proximal node of the child
        # branch matches the distal node of the parent branch
        prox = childBranch.getProximal()
        for bran in self.branchList:
            if bran.getDistal() == prox:
                return bran

    def countTerminals(self, branch):
        count = 0
        for curBranch in self.branchList:
            if curBranch.isTerminal():
                while not curBranch.isRoot() and curBranch != branch:
                    curBranch = curBranch.getParents()[0]
                if curBranch == branch:
                    count += 1
        return count

    def setRadii(self, initial):
        # To assign the radii, count the # of terminals below branch, and the number of terminals below it's parent branch
        # This will give the ratio of area the branch should take from its parent branch according to our formula
        # rnew = sqrt(ratio * rold^2)
        # Set initial radius
        self.branchList[0].setRadius(initial)
        # Set the subsequent radii
        for i in range(1,len(self.branchList)):
            if not self.branchList[i].getParents():
                raise RuntimeError('Non root branch has no parents')
            ratio = self.countTerminals(self.branchList[i])/self.countTerminals(self.branchList[i].getParents()[0])
            self.branchList[i].setRadius(np.sqrt(ratio * self.branchList[i].getParents()[0].getRadius()**2))

    def branchLength(self, branch):
        # Walk from distal to proximal node
        curNode = branch.getDistal()
        length = 0
        while not(curNode == branch.getProximal()):
            length += np.linalg.norm(curNode.getLocation() - curNode.getParents()[0].getLocation())
            curNode = curNode.getParents()[0]
        return length

    def setResistance(self, mu):
        # Convert mu from centipoise to dyn-s/cm^2
        mu = mu * 0.01
        for branch in self.branchList:
            L = self.branchLength(branch)
            R = (8 * mu * L) / (np.pi * branch.getRadius()**4)
            branch.setResistance(R)

    def setBranchVolume(self):
        for branch in self.branchList:
            rad = branch.getRadius()
            length = self.branchLength(branch)
            branch.setVolume(np.pi * rad**2 * length)

    def percentFull(self):
        fluid = 0
        vol = 0
        for branch in self.branchList:
            fluid += branch.getFluid()
            vol += branch.getVolume()
        return fluid/vol

    def totalVolume(self):
        self.setBranchVolume()
        vol = 0.
        for branch in self.branchList:
            vol += branch.getVolume()
        return vol

    def totalLength(self):
        length = 0.
        for branch in self.branchList:
            length += self.branchLength(branch)
        return length

    def findNodeByLoc(self, loc):
        eps = 1e-6
        for node in self.nodeList:
            if np.linalg.norm(node.getLocation() - loc) < eps:
                return node

    def findBranchByLoc(self, prox, dist):
        eps = 1e-6
        proxLoc = prox.getLocation()
        distLoc = dist.getLocation()
        for branch in self.branchList:
            if np.linalg.norm(branch.getProximal().getLocation() - proxLoc) < eps and np.linalg.norm(branch.getDistal().getLocation() - distLoc) < eps:
                return branch

    def findSLByLoc(self, loc):
        eps = 1e-6
        for sl in self.slList:
            if np.linalg.norm(sl.getLocation() - loc) < eps:
                return sl

    def setRSL(self, value):
        for sl in self.slList:
            sl.setResistance(value)
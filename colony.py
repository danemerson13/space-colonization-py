import attractor, node, superlobule, branch
import numpy as np
import plotter
import imageio, os

class Colony:
    def __init__(self, D, dk, di, type, targets = np.array([]), attractors = np.array([]), root = None, nodes = None):

        # Expectation is that either an empty tree with a root node is passed
        # OR an initialized tree with node and branchLists is passed
        self.D = D
        self.dk = dk
        self.di = di

        # Assign type - either "Inlet" or "Outlet"
        self.type = type

        self.attractorList = list()
        self.branchList = list()
        self.nodeList = list()
        if nodes:
            for inNode in nodes:
                self.nodeList.append(inNode)
                self.nodeList[-1].setTerminal(False)

        # Add all the target attractors to the attractorList
        for i in range(len(targets)):
            self.attractorList.append(attractor.Attractor(targets[i,:], True))
        # Add all the regular attractors to the attractorList
        for i in range(len(attractors)):
            self.attractorList.append(attractor.Attractor(attractors[i,:], False))
        # Add the root node if there is one
        if root is not None:
            self.nodeList.append(node.Node(root, type))

    def getLocationArray(self, inputList):
        # Returns np array with euclidean coordinates of input list (attractors, nodes, etc.)
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
            if self.nodeList[i].isTerminal():
                arr[i,:] = np.full((3,), np.inf)
            elif not initflag and self.nodeList[i].isRoot():
                arr[i,:] = np.full((3,), np.inf)
        return arr

    def numTargets(self):
        # Counts the number or target attractors in attractorList
        ct = 0
        for i in range(len(self.attractorList)):
            if self.attractorList[i].isTarget():
                ct += 1
        return ct

    def computeDistanceVec(self,arr,pt):
        # Computes the euclidean distance between pt and every point in arr
        return np.linalg.norm(arr - pt, axis = 1)

    def findClosestNode(self, initflag):
        # For each of the attractors, return a list with a pointer to the closest node, 
        # permitting that the attractor is within the sphere of influence, di. Otherwise return None
        nodeArr = self.getNodeArray(initflag)
        closestNode = list()
        for i in range(len(self.attractorList)):
            dist = self.computeDistanceVec(nodeArr, self.attractorList[i].getLocation())
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
        # This was an edge case where the algorithm would get stuck between two attractors 
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

    def addNode(self, dir, growthNode):
        # Add a node in direction dir, D distance from node
        loc = growthNode.getLocation() + dir * self.D
        while self.repeatedStepCheck(loc, growthNode):
            print("REPEATED STEP")
            dir = dir * 0.5
            loc = growthNode.getLocation() + dir * self.D
        self.nodeList.append(node.Node(loc, self.type, parent = growthNode))
        # Also update parent node to have this new node as its child
        growthNode.setChild(self.nodeList[-1])

    def addTerminalNode(self, dir, D_exact, growthNode):
        # Add a *terminal* node in direction dir, D_exact distance from node
        loc = growthNode.getLocation() + dir * D_exact
        self.nodeList.append(node.Node(loc, self.type, parent = growthNode, terminal = True))
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

    def createTree(self, initflag):
        prevNodeCount = len(self.nodeList)
        loop = 0
        while self.numTargets() > 0:
            # if loop == 2500:
            #     print("Pause")
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
                    self.addTerminalNode(dir, D_exact, growthNode)
                    # Flag this target attractor as killed
                    self.attractorList[activeAttractors[0]].setKill()
                else:
                    # Compute direction to step in
                    dir = self.attractorNormal(activeAttractors, growthNode)
                    # Add node in this direction at distance D
                    self.addNode(dir, growthNode)
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
            loop += 1

    def createBranches(self):
        # Will create branches by taking all of the distal nodes in the tree
        # These will either be leaf nodes, or furcation ndoes
        distList = list()
        # Leaf nodes have 0 children, furcation nodes have > 1 child
        for node in self.nodeList:
            nChildren = len(node.getChildren())
            if nChildren == 0 or nChildren > 1 and len(node.getParents()) > 0:
                distList.append(node)
        # From the distal list we can traverse each branch backwards until we reach 
        # a node with multiple children, or no parent in the case of the root branch
        for node in distList:
            prox, dist = self.traverseBranch(node)
            self.branchList.append(branch.Branch(self.type, prox, dist))
        # Now assign parent - child relationships
        for bran in self.branchList:
            parent = self.findParentBranch(bran)
            # Assign to both
            if parent != None:
                bran.setParent(parent)
                parent.setChild(bran)
            elif not bran.isRoot():
                print('flag')

    def traverseBranch(self, endNode):
        # From the distal list we can traverse each branch backwards until we reach 
        # a node with multiple children, or no parent in the case of the root branch
        dist = endNode
        curNode = endNode.getParents()[0]
        while len(curNode.getChildren()) < 2 and len(curNode.getParents()) > 0:
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

    def trimNodes(self):
        criticalNodes = list()
        for curNode in self.nodeList:
            if curNode.isTerminal():
                criticalNodes.append(curNode)
                while len(curNode.getParents()) > 0:
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
            if len(self.branchList[i].getParents()) == 0:
                print('Err')
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

    def setBranchVolume(self):
        for branch in self.branchList:
            rad = branch.getRadius()
            length = self.branchLength(branch)
            branch.setVolume(np.pi * rad**2 * length)

    def setResistance(self, mu):
        # Convert mu from centipoise to dyn-s/cm^2
        mu = mu * 0.01
        for branch in self.branchList:
            L = self.branchLength(branch)
            R = (8 * mu * L) / (np.pi * branch.getRadius()**4)
            branch.setResistance(R)

    def solveResistanceNetwork(self, Pin, Pout):
        # Build square A matrix to solve linear system Ax = b
        # x is the all of the branchwise flow rates and nodal pressures, concatenated into one vector
        # A will be comprised of three types of equations
        # (1) Branch eqns: Vdist - Vprox - IR = 0 (Ohms law)
        # (2) Nodal eqns: Iin - sum(Iout) = 0 (Kirchoffs current law)
        # (3) BC eqns: V = const (Direct assignment)
        
        # Need to assign necessary boundary conditions. Input pressure, output pressure
        nBranch = len(self.branchList)
        nNode = 0
        nBC = 0
        VNodeList = list()
        for node in self.nodeList:
            if node.isRoot():
                node.setPressure(Pin)
                VNodeList.append(node)
                nNode += 1
                nBC += 1
            elif node.isTerminal():
                node.setPressure(Pout)
                VNodeList.append(node)
                nNode += 1
                nBC += 1
            elif node.isFurcation():
                VNodeList.append(node)
                nNode += 1

        A = np.zeros((nBranch+nNode, nBranch+nNode))
        b = np.zeros(nBranch+nNode)
        # x = [branchList, nodeList]

        # (1) Branch eqns: vdist - vprox - IR = 0 (Ohms law)
        for i, branch in enumerate(self.branchList):
            A[i,i] = -branch.getResistance()
            A[i, VNodeList.index(branch.getProximal()) + nBranch] = 1
            A[i, VNodeList.index(branch.getDistal()) + nBranch] = -1

        # (2) Nodal eqns: Iin - Iout1 - Iout2 = 0 (Kirchoffs current law)
        j = nBranch
        for branch in self.branchList:
            if not branch.isTerminal():
                A[j, self.branchList.index(branch)] = 1
                for child in branch.getChildren():
                    A[j, self.branchList.index(child)] = -1
                j += 1

        # (3) BC eqns
        k = len(A) - nBC
        for node in VNodeList:
            if node.isRoot() or node.isTerminal():
                A[k, VNodeList.index(node) + nBranch] = 1
                b[k] = node.getPressure()
                k += 1

        # Solve Ax = b
        x = np.linalg.solve(A, b)

        # Assign computed I and V values
        for i, branch in enumerate(self.branchList):
            branch.setFlowRate(x[i])

        j = nBranch
        for node in VNodeList:
            node.setPressure(x[j])
            j += 1

    def findNodeByLoc(self, loc):
        eps = 1e-10
        loc = loc.getLocation()
        for curnode in self.nodeList:
            if np.linalg.norm(loc- curnode.getLocation()) < eps:
                return curnode

    def findMatchingNode(self, node):
        eps = 1e-6
        for branch in self.branchList:
            if branch.getDistal().isTerminal():
                if np.linalg.norm(node.getLocation() - branch.getDistal().getLocation()) < eps:
                    return branch.getDistal()

    def findMatchingBranch(self, branch):
        eps = 1e-6
        for bran in self.branchList:
            if bran.getDistal().isTerminal():
                if np.linalg.norm(bran.getDistal().getLocation() - branch.getDistal().getLocation()) < eps:
                    return bran

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

    def fillStep(self, dt):
        # Start with root branch
        V = self.branchList[0].getFlowRate()*dt
        self.branchList[0].addFluid(V)
        # Now step through all the branches and move around the excess fluid
        for branch in self.branchList:
            # If the branch has excess fluid and isnt terminal, transfer to its children
            if branch.getBuffer() > 0 and not branch.isTerminal():
                Q = branch.getFlowRate()
                # For each child give its corresponding amount of excess fluid based on proportion of total flow rate
                for child in branch.getChildren():
                    ratio = child.getFlowRate() / Q
                    child.addFluid(branch.getBuffer() * ratio)
                # Zero out the buffer from the parent branch
                branch.setBuffer(0)
            else:
                # Otherwise can just zero the buffer on the terminal nodes
                branch.setBuffer(0)
        
    def fillNetwork(self, dt):
        # Assign volumes
        self.setBranchVolume()
        with imageio.get_writer('filling.gif', mode  = 'I') as writer:
            # Begin iteration
            time = 0
            plotter.plotFilling(self, time)
            image = imageio.imread('fill_%.2f.png' %time)
            writer.append_data(image)
            os.remove('fill_%.2f.png' %time)
            while self.percentFull() < 1:
                self.fillStep(dt)
                time += dt
                plotter.plotFilling(self, time)
                image = imageio.imread('fill_%.2f.png' %time)
                writer.append_data(image)
                os.remove('fill_%.2f.png' %time)
  
def solveResistanceNetwork2x(inlet, outlet, Pin, Pout):
        nBranchIn = len(inlet.branchList); nBranchOut = len(outlet.branchList)
        nNode = 0
        nBC = 0
        VNodeList = list()
        for node in inlet.nodeList:
            if node.isRoot():
                node.setPressure(Pin)
                VNodeList.append(node)
                nNode += 1
                nBC += 1
            elif node.isFurcation() or node.isTerminal():
                VNodeList.append(node)
                nNode += 1
        for node in outlet.nodeList:
            if node.isRoot():
                node.setPressure(Pout)
                VNodeList.append(node)
                nNode += 1
                nBC += 1
            elif node.isFurcation():
                VNodeList.append(node)
                nNode += 1

        nBranch = nBranchIn + nBranchOut
        
        A = np.zeros((nBranch+nNode, nBranch+nNode))
        b = np.zeros(nBranch+nNode)
        # x = [branchList, nodeList]

        # (1) Branch eqns: vdist - vprox - IR = 0 (Ohms law)
        for i, branch in enumerate(inlet.branchList):
            A[i,i] = -branch.getResistance()
            A[i, VNodeList.index(branch.getProximal()) + nBranch] = 1
            A[i, VNodeList.index(branch.getDistal()) + nBranch] = -1

        j = nBranchIn
        for branch in outlet.branchList:
            A[j,j] = -branch.getResistance()
            # Reverse signs for outlet tree since direction is reversed
            A[j, VNodeList.index(branch.getProximal()) + nBranch] = -1
            if branch.isTerminal():
                node = inlet.findMatchingNode(branch.getDistal())
                A[j, VNodeList.index(node) + nBranch] = 1
            else:
                A[j, VNodeList.index(branch.getDistal()) + nBranch] = 1
            j += 1

        # (2) Nodal eqns: sum(Iin) - sum(Iout) = 0 (Kirchoffs current law)
        k = nBranch
        for branch in inlet.branchList:
            # Normal branches
            if not branch.isTerminal():
                A[k, inlet.branchList.index(branch)] = 1
                for child in branch.getChildren():
                    A[k, inlet.branchList.index(child)] = -1
                k += 1
            # Terminal branches
            if branch.isTerminal():
                # Find branch in outlet that matches distal node
                A[k, inlet.branchList.index(branch)] = 1
                bran = outlet.findMatchingBranch(branch)
                A[k, outlet.branchList.index(bran) + nBranchIn] = -1
                k += 1
        for branch in outlet.branchList:
            if not branch.isTerminal():
                # Reverse signs for outlet tree since direction is reversed
                A[k, outlet.branchList.index(branch) + nBranchIn] = -1
                for child in branch.getChildren():
                    A[k, outlet.branchList.index(child) + nBranchIn] = 1
                k += 1

        # (3) BC eqns
        l = len(A) - nBC
        for node in VNodeList:
            if node.isRoot():
                A[l, VNodeList.index(node) + nBranch] = 1
                b[l] = node.getPressure()
                l += 1

        # Solve Ax = b
        x = np.linalg.solve(A, b)

        # Assign computed I and V values
        for i, branch in enumerate(inlet.branchList):
            branch.setFlowRate(x[i])
        
        j = nBranchIn
        for branch in outlet.branchList:
            branch.setFlowRate(x[j])
            j += 1

        k = nBranch
        for node in VNodeList:
            node.setPressure(x[k])
            k += 1

def reverseColony(col):
    # Takes inlet colony and reverses to be outlet tree
    new = Colony(col.D, col.dk, col.di, "Outlet")
    # Step through all the nodes and reverse their parent and child lists
    for curnode in col.nodeList:
        new.nodeList.append(node.Node(curnode.getLocation(), "Outlet"))
    for i, curnode in enumerate(new.nodeList):

        # NEED TO SET CHILDREN AND PARENTS TO CORRESPONDING NODE IN NEW DATA STRUCTURE, CURRENTLY POINTING TO THE OLD NODES...
        
        for parent in col.nodeList[i].getParents():
            curnode.setChild(parent)
        for child in col.nodeList[i].getChildren():
            curnode.setParent(child)
    # Step through all the branches and reverse their prox/distal nodes, and parent child lists
    for curbranch in col.branchList:
        new.branchList.append(branch.Branch("Outlet", prox = new.findNodeByLoc(curbranch.getDistal()), dist = new.findNodeByLoc(curbranch.getProximal())))
    for i, curbranch in enumerate(new.branchList):
        for parent in col.branchList[i].getParents():
            curbranch.setChild(parent)
        for child in col.branchList[i].getChildren():
            curbranch.setParent(child)
    return new

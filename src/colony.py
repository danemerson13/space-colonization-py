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
        # Add the attractors and target attractors
        self.addAttractors(targets, attractors)
        # Add the root node
        self.nodeList.append(node.Node(root))

        # Run SC Algorithm
        self.runSpaceColonization()

        # Trim the node list
        self.trimNodes()

    def createTree(self, targets, attractors, initNodeList, Rinitial, mu):
        # Add the attractors and target attractors
        self.addAttractors(targets, attractors)
        # Add the initialized nodes
        if type(initNodeList) == np.ndarray:
            self.nodeList.append(node.Node(initNodeList))
        else:
            for inNode in initNodeList:
                self.nodeList.append(inNode)
                self.nodeList[-1].setSL(False)

        # Run SC Algorithm
        self.runSpaceColonization()

        # Trim nodes and create branches
        self.trimNodes()
        self.createBranches()

        # Assign radii and resistances
        self.setRadii(Rinitial)
        self.setResistance(mu)

        # Create and link segments
        # self.createSegments()

    def createTreeNoTrim(self, targets, attractors, initNodeList, Rinitial, mu):
        # Add the attractors and target attractors
        self.addAttractors(targets, attractors)
        # Add the initialized nodes
        if type(initNodeList) == np.ndarray:
            self.nodeList.append(node.Node(initNodeList))
        else:
            for inNode in initNodeList:
                self.nodeList.append(inNode)
                self.nodeList[-1].setSL(False)

        # Run SC Algorithm
        self.runSpaceColonization()

        # Trim nodes and create branches
        # self.trimNodes()
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
                prox = branch.getDistal()
                dist = outTree.findNodeByLoc(branch.getDistal().getLocation())
                self.slList.append(superlobule.SuperLobule(branch.getDistal().getLocation(), prox, dist))
                # Pair the nodes from the inlet and outlet trees together
                prox.setChild(dist)
                dist.setParent(prox)
                # Pair the branches from the inlet and outlet trees together
                distBranch = outTree.findChildBranchCrossTree(branch)
                branch.setChild(distBranch)
                distBranch.setParent(branch)

        # Add the nodes and branches from the outlet tree
        for node in outTree.nodeList:
            self.nodeList.append(node)
        for branch in outTree.branchList:
            self.branchList.append(branch)

    def solveResistanceNetwork(self, Pin, Pout, flag = False):
        # Convert Pressure from mmHg to dyn/mm^2
        Pin = Pin * 13.3322
        Pout = Pout * 13.3322

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

        # (1.1) Branch eqns: vprox - vdist - IR = 0 (Ohms law)
        for i, branch in enumerate(self.branchList):
            A[i,i] = -branch.getResistance()
            A[i, VNodeList.index(branch.getProximal()) + nBranch + nSL] = 1
            A[i, VNodeList.index(branch.getDistal()) + nBranch + nSL] = -1

        # (1.2) SL eqns: vprox - vdist - IR = 0 (Ohms law)
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

        # Print statistics about matrix if flag True
        if flag:
            print("2 Norm: ", np.linalg.norm(A, 2))
            print("Frobenius Norm: ", np.linalg.norm(A, 'fro'))
            print("Infinity Norm: ", np.linalg.norm(A, np.inf))
            print("-Infinity Norm: ", np.linalg.norm(A, -np.inf))
            print("Condition Number: ", np.linalg.cond(A))

    def solveRSL(self, Qactual, Pin, Pout, n_iter, tol, verbal, a = 0, b = 1e-4):
        # Solves for RSL value via bisection method with brackets a and b initialized such that f(a) * f(b) < 0, 
        # terminating when Qi is within tol of Qactual, 
        # or after iter iterations

        # verbal = 0: print nothing
        # verbal = 1: print final RSL, iter, tol
        # verbal = 2: print 1, and starting bounds
        # verbal = 3: print 2, and RSL and tol at each iter

        start = time.time()

        # Convert Qactual from mL/min to mm^3/s (1000/60)
        Qactual = Qactual * 50/3 # Same as Qact * 16.6666

        # Convert Pressure from mmHg to dyn/mm^2
        Pin = Pin * 13.3322
        Pout = Pout * 13.3322

        # First need to find brackets a,b for which Q(a) and Q(b) lie on opposite sides of Qactual
        a = a; Qa, _ = self.queryQin(a, Pin, Pout)
        b = b; Qb, _ = self.queryQin(b, Pin, Pout)
        while (Qa - Qactual) * (Qb - Qactual) >= 0:
            b *= 10
            Qb, _ = self.queryQin(b, Pin, Pout)
        if verbal >= 2:
            print("Starting bisection method with initial bounds {%d, %d}" %(a, b))

        # Save the Req for RSL = 0
        Req = (Pin - Pout)/Qa

        # Now iterate until tolerance is met or n_iter is broken
        t_list = list()
        c = (a+b)/2
        Qc, t = self.queryQin(c, Pin, Pout); t_list.append(t)
        iter = 0
        while np.abs(Qactual - Qc)/Qactual > tol and iter < n_iter:
            if verbal >= 3:
                print("Iter: %d, RSL = %.2f, Tolerance: %.2E" %(iter, c, np.abs(Qactual - Qc)/Qactual))
            # Bisection Method logic
            if (Qa - Qactual) * (Qc - Qactual) < 0:
                b = c
                Qb = Qc
            else:
                a = c
                Qa = Qc
            # Update c
            c = (a+b)/2
            Qc, t = self.queryQin(c, Pin, Pout); t_list.append(t)
            # Update counter
            iter += 1

        end = time.time()
        elapsed = end - start
        avg_mat_solve = np.mean(t_list)

        if iter >= n_iter:
            print("DID NOT CONVERGE - Method terminated after %d iterations with RSL = %.2f and a tolerance of %.2E" %(iter, c, np.abs(Qactual - Qc)/Qactual))
        else:
            if verbal >= 1:
                nBranch, nNode, nSL = self.getRNetworkDims()
                # print("Solution converged to RSL = %.2f dyn-s/mm^5, Q(RSL) = %.2f mL/min, after %d iterations to a tolerance of %.2E" %(c, Qc * 3/50, iter, np.abs(Qactual - Qc)/Qactual))
                print("Solution for %d SL converged to RSL = %.2f dyn-s/mm^5, after %d iterations and %.2f secs to a tolerance of %.2E, Matrix Dim: %d x %d, nBranch: %d, nNode: %d, nSL: %d, Avg Matrix Solve %.2E" %(len(self.slList), c, iter, elapsed, np.abs(Qactual - Qc)/Qactual, nBranch+nNode+nSL, nBranch+nNode+nSL, nBranch, nNode, nSL, avg_mat_solve))
        # Check for backflow
        self.flowCheck()
        # Return the RSL from bisection method
        return c

##### PRIMARY FUNCTIONS #####

    def runSpaceColonization(self):
        prevNodeCount = len(self.nodeList)
        while self.numTargets() > 0:
            # First find the node closest to each attractor
            # Note: this is returned as a list with pointers to the closest nodes, 
            # with None returned where there are no nodes within di
            influencers = self.findClosestNode()
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
            # Check to see if no steps were taken, in this case: increase di (sphere of influence) to find growthNodes
            nodeCount = len(self.nodeList)
            if nodeCount == prevNodeCount:
                self.di = self.di*2
                print("Updated di to %d" %(self.di))
                # output = 'Broken with ' + str(self.numTargets()) + ' targets not reached.'
                # print(output)
                # break
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

    def fillNetwork(self, dt):
        # Assign volumes
        self.setBranchVolume()
        with imageio.get_writer('filling.gif', mode  = 'I') as writer:
            # Begin iteration
            time = 0
            plotter.plotFilling(self, time, excess = 0)
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

    def overfillNetwork(self, dt, path):
        # Assign volumes and print total vessel volume
        # (Volume assignment done within self.totalVolume() function)
        totalVol = self.totalVolume()
        print("Total Vessel Volume: ", totalVol)
        overfill = 100000 - totalVol
        with imageio.get_writer('filling.gif', mode  = 'I') as writer:
            # Begin iteration
            time = 0
            excess = 0.
            plotter.plotFilling(self, time, excess)
            image = imageio.imread('fill_%.2f.png' %time)
            writer.append_data(image)
            os.remove('fill_%.2f.png' %time)
            trip = False
            while excess < overfill:
                excess += self.fillStep(dt)
                time += dt
                percFull = self.percentFull()
                if percFull > 0.9 and trip == False:
                    plotter.plotFilling(self, time, excess, path)
                else:
                    plotter.plotFilling(self, time, excess)
                if self.getFillTime() == None:
                    if percFull >= 1.:
                        self.setFillTime(time)
                image = imageio.imread('fill_%.2f.png' %time)
                writer.append_data(image)
                os.remove('fill_%.2f.png' %time)

        self.setOverfillTime(time)

    def saveModel(self, path):
        A = self.createResistanceMatrix()
        self.overfillNetwork(0.1, os.getcwd() + '/' + path + '/fill90.png')
        with open(os.getcwd() + '/' + path + '/model.pkl', 'wb') as activeFile:
            pickle.dump(self, activeFile)
        np.save(os.getcwd() + '/' + path + '/matrix.npy', A)
        
    def saveTopo(self, path):
        with open(os.getcwd() + '/' + path + '/model.pkl', 'wb') as activeFile:
            pickle.dump(self, activeFile)
        plotter.plotRadius(self, path)

##### HELPER FUNCTIONS #####

    def findRootBranch(self):
        for branch in self.branchList:
            if branch.getType() == "Inlet" and branch.isRoot():
                return branch
        raise RuntimeError('No root branch found')

    def fillStep(self, dt):
        excess = 0.
        rootBranch = self.findRootBranch()
        V = rootBranch.getFlowRate() * dt
        rootBranch.addFluid(V)
        while self.checkBuffers() > 0:
            excess += self.transferBuffers()
        return excess
    
    def transferBuffers(self):
        # Moves buffers onto child or excess
        excess = 0.
        for branch in self.branchList:
            if branch.getBuffer() > 0:
                # If there are multiple children, transfer fluid according to ratio of flowrates
                if len(branch.getChildren()) > 1:
                    Q = branch.getFlowRate()
                    V = branch.getBuffer()
                    for child in branch.getChildren():
                        ratio = child.getFlowRate() / Q
                        child.addFluid(V * ratio)
                # If there is a single child, transfer all fluid to child
                elif len(branch.getChildren()) == 1:
                    branch.getChildren()[0].addFluid(branch.getBuffer())
                # If there are no children transfer the fluid to excess (this should only be triggered once per transferBuffers() call)
                else:
                    excess += branch.getBuffer()
                # Clear the buffer from the parent branch
                branch.setBuffer(0)
        return excess

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

    def getNodeArray(self):
        # Returns np array with locations of the nodes
        # Fill terminal nodes with inf since we dont want to grow from them
        # Never want to grow from the "root" node (two branches from root mess up the radius splitting) 
        # Note: node is not "Root" until it has a single child, so initializing the tree from a root is not a problem
        arr = self.getLocationArray(self.nodeList)
        for i in range(len(self.nodeList)):
            if self.nodeList[i].isSL():
                arr[i,:] = np.full((3,), np.inf)
            elif self.nodeList[i].isRoot():
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

    def findClosestNode(self):
        # For each of the attractors, return a list with a pointer to the closest node, 
        # permitting that the attractor is within the sphere of influence, di. Otherwise return None
        nodeArr = self.getNodeArray()
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
            # Edge case 3, similar to edge case 2, BUT, modified loc is more than D away from growth node, so the repeated step was being taken
            # Instead, step D in the direction of the modified loc
            else:
                dir = (modifiedLoc - growthNode.getLocation())/np.linalg.norm(modifiedLoc - growthNode.getLocation())
                loc = growthNode.getLocation() + dir * self.D

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

    def findChildBranchCrossTree(self, parentBranch):
        # Find the branch in branchList where the distal node of the parent 
        # branch matches the proximal node of the child branch from another tree
        dist = self.findNodeByLoc(parentBranch.getDistal().getLocation())
        for bran in self.branchList:
            if bran.getProximal() == dist:
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
        # If outlet branch:
        if branch.getType() == "Outlet":
            # Walk from proximal to distal node
            curNode = branch.getProximal()
            length = 0
            while not(curNode == branch.getDistal()):
                length += np.linalg.norm(curNode.getLocation() - curNode.getChildren()[0].getLocation())
                curNode = curNode.getChildren()[0]
            return length
        # If inlet branch:
        else:
            # Walk from distal to proximal node
            curNode = branch.getDistal()
            length = 0
            while not(curNode == branch.getProximal()):
                length += np.linalg.norm(curNode.getLocation() - curNode.getParents()[0].getLocation())
                curNode = curNode.getParents()[0]
            return length

    def setResistance(self, mu):
        # Convert mu from centipoise to dyn-s/mm^2
        mu = mu * 0.0001
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

    def queryQin(self, RSL, Pin, Pout):
        self.setRSL(RSL)
        start = time.time()
        self.solveResistanceNetwork(Pin, Pout)
        end = time.time()
        for branch in self.branchList:
            if branch.getType() == "Inlet" and branch.isRoot():
                return branch.getFlowRate(), end - start

    def checkBuffers(self):
        vol = 0
        for branch in self.branchList:
            vol += branch.getBuffer()
        return vol
    
    def emptyBuffers(self):
        excess = 0.
        while self.checkBuffers() > 0:
            for branch in self.branchList:
                if branch.getBuffer() > 0 and not branch.isTerminal():
                    Q = branch.getFlowRate()
                    # For each child give its corresponding amount of excess fluid based on proportion of total flow rate
                    if len(branch.getChildren()) > 1:
                        for child in branch.getChildren():
                            ratio = child.getFlowRate() / Q
                            child.addFluid(branch.getBuffer() * ratio)
                    else:
                        # For branches with only one child (SL branches and Outlet branches), just pass all fluid directly to child
                        branch.getChildren()[0].addFluid(branch.getBuffer())
                    # Zero out the buffer from the parent branch
                    branch.setBuffer(0)
                else:
                    # Otherwise can just zero the buffer on the terminal nodes
                    excess += branch.getBuffer()
                    branch.setBuffer(0)
        return excess

    def getRNetworkDims(self):
        nBranch = len(self.branchList)
        nNode = 0
        nSL = len(self.slList)
        for node in self.nodeList:
            if node.isSL():
                nNode += 1
            else:
                if node.getType() == "Root":
                    nNode += 1
                elif node.getType() == "Terminal":
                    nNode += 1
                elif node.getType() == "Inlet" or node.getType() == "Outlet": # Furcation node
                    nNode += 1

        return nBranch, nNode, nSL

    def getRadiiLims(self):
        rmax = -np.inf
        rmin = np.inf
        for branch in self.branchList:
            rad = branch.getRadius()
            if rad > rmax:
                rmax = rad
            if rad < rmin:
                rmin = rad
        return rmin, rmax

    def setReynoldsNumber(self, rho, mu):
        # rho: g/mm^3
        mu = mu / 100 #convert cP to g/(mm-s)
        for branch in self.branchList:
            Q = branch.getFlowRate() / 1000 # mL/s
            Re = 2 * rho * Q / (mu * np.pi * branch.getRadius())
            branch.setReynolds(Re)

    def getReynoldsLims(self):
        remax = -np.inf
        remin = np.inf

        for branch in self.branchList:
            re = branch.getReynolds()
            if re > remax:
                remax = re
            if re < remin:
                remin = re
        return remin, remax

    def flowCheck(self):
        # Get the direction of the inlet branch
        for branch in self.branchList:
            if branch.getType() == "Inlet" and branch.isRoot():
                dir = branch.getFlowRate()
        for branch in self.branchList:
            if branch.getFlowRate() * dir < 0:
                raise RuntimeError('Backflow: computed flows have inconsistent direction')

    def setGeneration(self,omega,alpha):
        # Call on trees before merge
        # Walk through trees and assign generation numbers to branches according 
        # to proximity to root and adherence to the following heuristics:

        # Branch considered the same iff two conditions are met
        # Radius condition -> r_new/r_old > ω
        # Angle condition -> angle between parent & child branch > α

        # Maintain a list counting the number of merged branches at each generation level
        # This will be used to correct the n_branches when computing generational statistics
        mergedCounter = list([0])
        # Start with root branch
        self.branchList[0].setGeneration(1)
        for i in range(1, len(self.branchList)):
            ratio = self.branchList[i].getRadius() / self.branchList[i].getParents()[0].getRadius()
            angle = self.computeBranchAngle(self.branchList[i])
            if ratio >= omega and angle >= alpha:
                self.branchList[i].setGeneration(self.branchList[i].getParents()[0].getGeneration())
                # Add to merged branch counter
                mergedCounter[self.branchList[i].getParents()[0].getGeneration()-1] += 1
            else:
                gen = self.branchList[i].getParents()[0].getGeneration() + 1
                self.branchList[i].setGeneration(gen)
                if len(mergedCounter) < gen:
                    mergedCounter.append(0)

        return mergedCounter

    def generationStatistics(self, omega, alpha):
        mergedCounter = self.setGeneration(omega,alpha*np.pi/180)
        nominalCounter = [0] * len(mergedCounter)
        nominalLength = [0] * len(mergedCounter)
        nominalRadius = [0] * len(mergedCounter)
        for branch in self.branchList:
            genIdx = branch.getGeneration() - 1
            nominalCounter[genIdx] += 1
            nominalLength[genIdx] += self.branchLength(branch)
            nominalRadius[genIdx] += branch.getRadius()

        genLength = np.array(nominalLength) / (np.array(nominalCounter) - np.array(mergedCounter))
        genRadius = np.array(nominalRadius) / np.array(nominalCounter)

        for i in range(len(mergedCounter)):
            print("Gen #%d, Nominal = %d, Actual = %d, L = %.2f, R = %.2f" %(i+1, nominalCounter[i], nominalCounter[i]-mergedCounter[i] ,genLength[i], genRadius[i]))

    def generationBoxPlots(self, omega, alpha):
        mergedCounter = self.setGeneration(omega,alpha*np.pi/180)
        lengths = [[] for i in range(len(mergedCounter))]
        radii = [[] for i in range(len(mergedCounter))]
        for branch in self.branchList:
            genIdx = branch.getGeneration() - 1
            lengths[genIdx].append(self.branchLength(branch))
            radii[genIdx].append(branch.getRadius())

        plt.boxplot(lengths)
        plt.xlabel("Generation")
        plt.ylabel("Length")
        plt.show()

        plt.boxplot(radii)
        plt.xlabel("Generation")
        plt.ylabel("Radius")
        plt.show()

    def computeBranchAngle(self, branch):
        # Compute angle between branch and parent branch, mod pi
        v_branch = branch.getDistal().getLocation() - branch.getProximal().getLocation()
        v_parent = branch.getParents()[0].getProximal().getLocation() - branch.getParents()[0].getDistal().getLocation()
        return np.arccos(np.dot(v_branch, v_parent)/ (np.linalg.norm(v_branch) * np.linalg.norm(v_parent)))

    def branchAngleHist(self):
        alpha = list()
        for i in range(1, len(self.branchList)):
            alpha.append(self.computeBranchAngle(self.branchList[i])*180/np.pi)
        plt.hist(alpha,20)
        plt.ylabel("Frequency")
        plt.xlabel(r"$\alpha (deg)$")
        plt.show()

    def branchRadiusRatioHist(self):
        omega = list()
        for i in range(1, len(self.branchList)):
            omega.append(self.branchList[i].getRadius() / self.branchList[i].getParents()[0].getRadius())
        plt.hist(omega,20)
        plt.ylabel("Frequency")
        plt.xlabel(r"$\omega = \frac{r_{child}}{r_{parent}}$")
        plt.show()
        
    def branchAlphaOmegaHist(self):
        alpha = list()
        omega = list()
        for i in range(1, len(self.branchList)):
            alpha.append(self.computeBranchAngle(self.branchList[i])*180/np.pi)
            omega.append(self.branchList[i].getRadius() / self.branchList[i].getParents()[0].getRadius())
        h = plt.hist2d(np.array(alpha),np.array(omega),20)
        plt.colorbar(h[3])
        plt.xlabel(r"$\alpha (deg)$")
        plt.ylabel(r"$\omega = \frac{r_{child}}{r_{parent}}$")
        plt.show()

    def setSLVolume(self, bloodVolume):
        totalSLVolume = bloodVolume - self.totalVolume()
        SLVolume = totalSLVolume / len(self.slList)
        print(SLVolume)
        for sl in self.slList:
            sl.setVolume(SLVolume)

    def createResistanceMatrix(self):
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
                    VNodeList.append(node)
                    nNode += 1
                    nBC += 1
                elif node.getType() == "Terminal":
                    VNodeList.append(node)
                    nNode += 1
                    nBC += 1
                elif node.getType() == "Inlet" or node.getType() == "Outlet": # Furcation node
                    VNodeList.append(node)
                    nNode += 1

        A = np.zeros((nBranch+nNode+nSL, nBranch+nNode+nSL))

        # (1.1) Branch eqns: vprox - vdist - IR = 0 (Ohms law)
        for i, branch in enumerate(self.branchList):
            A[i,i] = -branch.getResistance()
            A[i, VNodeList.index(branch.getProximal()) + nBranch + nSL] = 1
            A[i, VNodeList.index(branch.getDistal()) + nBranch + nSL] = -1

        # (1.2) SL eqns: vprox - vdist - IR = 0 (Ohms law)
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
                    k += 1

        return A

    def setFillTime(self, t):
        self.fillTime = t

    def getFillTime(self):
        return self.fillTime

    def setOverfillTime(self, t):
        self.overfillTime = t

    def getOverfillTime(self):
        return self.overfillTime

    def getBranchNodes(self, branch):
        nodes = list()
        if branch.getType() == "Inlet":
            # If inlet tree go dist -> prox
            prox = branch.getProximal()
            nodeA = branch.getDistal()
            while nodeA != prox:
                nodes.append(nodeA)
                nodeA = nodeA.getParents()[0]
            nodes.append(nodeA)
            nodes.reverse()
        else:
            # If outlet tree go prox -> dist
            dist = branch.getDistal()
            nodeA = branch.getProximal()
            while nodeA != dist:
                nodes.append(nodeA)
                nodeA = nodeA.getChildren()[0]
            nodes.append(nodeA)
        return nodes

    def createSegments(self, lmax):
        rootflag = False
        for branch in self.branchList:
            nodes = self.getBranchNodes(branch)
            for i in range(len(nodes)-1):
                length = np.linalg.norm(nodes[i].getLocation() - nodes[i+1].getLocation())
                # Compute the number of segments for the following nodes
                segNum = np.ceil(length/lmax)
                prox = nodes[i].getLocation()
                dir = (nodes[i+1].getLocation() - nodes[i].getLocation())/np.linalg.norm(nodes[i+1].getLocation() - nodes[i].getLocation())
                step = dir * length/segNum
                for _ in range(int(segNum)):
                    dist = prox + step
                    self.segList.append(segment.Segment(prox, dist, branch))
                    # Check to see if we marked the root segment yet
                    if rootflag == False and branch.isRoot():
                        self.segList[-1].setRoot()
                        rootflag = True
                    prox = dist

    def connectSegments(self):
        for segment1 in self.segList:
            for segment2 in self.segList:
                if segment1 != segment2:
                    # Find segments sharing dist and proximal nodes
                    if np.linalg.norm(segment1.getDistal() - segment2.getProximal()) < 1e-6:
                        segment1.setChild(segment2)
                        segment2.setParent(segment1)

    def getTotalConcentration(self):
        vol = 0.
        conc = 0.
        for seg in self.segList:
            vol += seg.getVolume()
            conc += seg.getVolume() * seg.getConcentration()
        return conc/vol

    def findRootSegment(self):
        for seg in self.segList:
            if seg.isRoot():
                return seg
    
    def fillBySegment(self, dt):
        # We will pass all fluids between segments before
        time = 0.
        totalConc = self.getTotalConcentration()
        with imageio.get_writer('filling.gif', mode  = 'I') as writer:
            while totalConc < 1-1e-3:
                print("Total Concentration: ", totalConc)

                # Plot and Save Image
                path = os.getcwd() + ("/fill_%.3f.png" %(time))
                plotter.plotSegmentFilling(self,time,totalConc,path)
                image = imageio.imread('fill_%.3f.png' %time)
                writer.append_data(image)
                os.remove('fill_%.3f.png' %time)

                # Add fluid with concentration 1 to root segment
                rootSeg = self.findRootSegment()
                rootSeg.updateConcentration(Cin = 1, Vin = rootSeg.getFlowRate() * dt)
                # Build updated segment list
                newSegList = list()
                newSegList.append(rootSeg)
                for seg in self.segList:
                    # Now update concentrations in all segments of the i+1 segList except for the root segment (since we already updated it)
                    if seg != rootSeg:
                        newSegList.append(seg)
                        if len(seg.getParents()) > 1:
                            # If there are multiple parents weight their contribution to Cin by their proportional flow rates
                            Cin = 0.
                            flowRate = 0.
                            for parent in seg.getParents():
                                Cin += parent.getConcentration() * parent.getFlowRate()
                                flowRate += parent.getFlowRate()
                            Cin = Cin / flowRate
                        else:
                            Cin = seg.getParents()[0].getConcentration()
                            flowRate = seg.getParents()[0].getFlowRate()
                        newSegList[-1].updateConcentration(Cin = Cin, Vin = flowRate * dt)
                # Once complete, replace the old segList with the new segList
                self.segList = newSegList
                time += dt
                totalConc = self.getTotalConcentration()
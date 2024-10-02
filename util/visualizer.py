import numpy as np
import open3d as o3d
from util import bezier

def normalize(vec):
    norm = np.linalg.norm(vec)
    if norm == 0:
        return vec
    else:
        return vec / norm
    
def axisAngle(axis, theta):
    return np.array([
        [np.cos(theta)+axis[0]**2*(1-np.cos(theta)), axis[0]*axis[1]*(1-np.cos(theta)) - axis[2]*np.sin(theta), axis[0]*axis[2]*(1-np.cos(theta)) + axis[1]*np.sin(theta), 0],
        [axis[1]*axis[0]*(1-np.cos(theta)) + axis[2]*np.sin(theta), np.cos(theta) + axis[1]**2*(1-np.cos(theta)), axis[1]*axis[2]*(1-np.cos(theta)) - axis[0]*np.sin(theta), 0],
        [axis[2]*axis[0]*(1-np.cos(theta)) - axis[1]*np.sin(theta), axis[2]*axis[1]*(1-np.cos(theta)) + axis[0]*np.sin(theta), np.cos(theta) + axis[2]**2*(1-np.cos(theta)), 0],
        [0, 0, 0, 1]
    ])

def translation(point):
    t = np.eye(4)
    t[:3,3] = point
    return t

def createCircle(center, axis, radius, sliceCount):
    # Create a circle with specified radius about axis, located at center with number of vertices specified by slice count
    # Instantiate vertex list
    vertices = list()
    # Calculate theta increment
    dTheta = 2 * np.pi / sliceCount
    # Get the transformation matrix (move to origin, rotate, move back)
    M = translation(center) @ axisAngle(axis, dTheta) @ translation(-center)
    # Use a hard-coded random vector crossed with the axis to get a vector on the plane of the circle
    randomVec = normalize(np.array([0.2341, -0.32451, 0.3]))
    # NOTE: Code something to check that this random vector isnt collinear with the axis anywhere
    if np.linalg.norm(np.cross(axis, randomVec)) < 1e-6:
        raise ValueError("randomVec is collinear with axis of the cone")
    planeVec = normalize(np.cross(axis, randomVec))
    radiusVec = planeVec * radius
    # Generate the first point on the circle
    pt = center + radiusVec
    vertices.append(pt)
    # Add homogenous
    pt = np.concatenate((pt, np.array([1]))).reshape(-1,1)
    # Rotate around the circle
    for _ in range(sliceCount-1):
        pt = M @ pt
        vertices.append(pt[:3,0])

    return vertices
    
def rotateCircle(newCircle, vertex):
    # Align the zero index vertex of newCircle with the zero index vertex of the previous circle
    # Convert list of vertices to np array
    newCircle = np.asarray(newCircle)
    # Compute array of distances between vertex and newCircle vertices
    dists = np.linalg.norm(vertex - newCircle, axis = 1)
    # Find index of vertex on newCircle closest to vertex
    idx = np.argmin(dists)
    # Move the indices around, setting idx as the first vertex
    newCircle = np.vstack((newCircle[idx:], newCircle[:idx]))
    return [newCircle[i] for i in range(len(newCircle))]

def getBranchCoords(branch):
    coordList = list()
    # Check the type of branch. 
    # For Inlet branches go from distal to proximal using getParents()[0]
    # For Outlet branches go from proximal to distal using getChildren()[0]
    if branch.getType() == "Inlet":
        curNode = branch.getDistal()
        coordList.append(curNode.getLocation())
        while curNode != branch.getProximal():
            curNode = curNode.getParents()[0]
            coordList.append(curNode.getLocation())
        # For the inlet branches reverse the list so we go from root -> terminal
        coordList.reverse()
    elif branch.getType() == "Outlet":
        curNode = branch.getProximal()
        coordList.append(curNode.getLocation())
        while curNode != branch.getDistal():
            curNode = curNode.getChildren()[0]
            coordList.append(curNode.getLocation())
    else:
        raise TypeError("Branch has no type")
    return coordList

def skinPoints(points, stackCount, sliceCount, rad1, rad2, baseCap = False, topCap = False, color = [0.5, 0.5, 0.5]):
    # NOTE: would be smart if I had it create the not as a bunch of disjoint cones, but as a long tube (reduces nVertices)
    # Each stack corresponds to a circle along the branch, located at traj
    
    # Generate bezier segments for the branch
    bezierList = bezier.G1CubicSpline(points)
    # Discretize each bezier based on stack count
    t = np.linspace(0,1,stackCount)
    # First just compute the total length of the branch
    totalLength = 0
    for seg in bezierList:
        traj = seg.eval(t)
        for i in range(len(traj)-1):
            totalLength += np.linalg.norm(traj[i+1] - traj[i])
    # Create list to store all circles generated along the points
    circleList = list()
    # Track arc length of currently plotted points
    curLength = 0
    # Loop through each bezier to generate vertices
    for i, seg in enumerate(bezierList):
        traj, tang = seg.eval(t, tangent = True)
        for j in range(len(traj)):
            rad = rad1 + (rad2 - rad1) * curLength/totalLength
            circleList.append(createCircle(traj[j], normalize(tang[j]), rad, sliceCount))
            # Rotate the circle vertices to be as closely aligned with the previous one
            if i == 0 and j == 0:
                pass
            else:
                circleList[-1] = rotateCircle(circleList[-1], circleList[-2][0])
            # Update the current length
            if j < len(traj)-1:
                curLength += np.linalg.norm(traj[j+1] - traj[j])
    # Convert list of circles to list of vertices
    V = list()
    for circle in circleList:
        for vertex in circle:
            V.append(vertex)
    # Create faces
    totalStackCount = len(circleList)
    F = list()
    for i in range(totalStackCount-1):
        for j in range(sliceCount):
            F.append([sliceCount*i+j, sliceCount*i+(j+1)%sliceCount, sliceCount*(i+1)+(j+1)%sliceCount])
            F.append([sliceCount*i+j, sliceCount*(i+1)+(j+1)%sliceCount, sliceCount*(i+1)+j])
    if baseCap == True:
        V.append(points[0])
        baseIdx = len(V)-1
        for i in range(sliceCount):
            F.append([(i+1)%sliceCount, i, baseIdx])
    if topCap == True:
        V.append(points[-1])
        topIdx = len(V)-1
        for i in range(sliceCount):
            F.append([(totalStackCount-1)*sliceCount + i, (totalStackCount-1)*sliceCount + (i+1)%(sliceCount), topIdx])

    # Make the mesh
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(V)
    mesh.triangles = o3d.utility.Vector3iVector(F)
    mesh.compute_triangle_normals()
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color(color)

    return mesh

def computeVisRadii(branch):
    if branch.getType() == "Inlet":
        # Base radius
        if len(branch.getParents()) > 0:
            baseRad = (branch.getParents()[0].getRadius() + branch.getRadius())/2
        else:
            baseRad = branch.getRadius()
        # Top radius
        if len(branch.getChildren()) > 0 and not branch.isSL():
            childRadList = [child.getRadius() for child in branch.getChildren()]
            topRad = (max(childRadList) + branch.getRadius())/2
        else:
            topRad = branch.getRadius()
    elif branch.getType() == "Outlet":
        # Base radius
        if len(branch.getParents()) > 0 and not branch.isSL():
            parentRadList = [parent.getRadius() for parent in branch.getParents()]        
            baseRad = (max(parentRadList) + branch.getRadius())/2
        else:
            baseRad = branch.getRadius()
        # Top radius
        if len(branch.getChildren()) > 0:
            topRad = (branch.getChildren()[0].getRadius() + branch.getRadius())/2
        else:
            topRad = branch.getRadius()
    return baseRad, topRad

class visualizer:
    def __init__(self, meshList, fov_step):
        self.meshList = meshList
        self.fov_step = fov_step

    def skinTree(self, col, sliceCount, stackCount):
        meshList = list()
        for branch in col.branchList:
            if branch.getType() == "Inlet":
                # color = [1,0,0]
                color = [0,0,1]
            elif branch.getType() == "Outlet":
                color = [1,0,0]
            else:
                raise TypeError("Invalid branch type")
            coords = getBranchCoords(branch)
            baseRad, topRad = computeVisRadii(branch)
            meshList.append(skinPoints(coords, stackCount, sliceCount, baseRad, topRad, True, True, color))
        return meshList

    def customVisualizer(self):
        vis = o3d.visualization.Visualizer()
        vis.create_window()
        for mesh in self.meshList:
            vis.add_geometry(mesh)
        ctr = vis.get_view_control()
        print("FOV Initial: %.2f" %ctr.get_field_of_view())
        ctr.change_field_of_view(step = self.fov_step)
        print("New FOV: %.2f" %ctr.get_field_of_view())
        vis.run()
        vis.destroy_window()
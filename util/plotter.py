import open3d as o3d
import numpy as np
import copy

def plotByBranch(col):
    meshList = list()
    for branch in col.branchList:
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            if branch.getType() == "Inlet":
                meshList.append(createCylinder(coords[j],coords[j+1],branch.getRadius(), np.array([1,0,0])))
            if branch.getType() == "Outlet":
                meshList.append(createCylinder(coords[j],coords[j+1],branch.getRadius(), np.array([0,0,1])))

    print("Hold")
    o3d.visualization.draw_geometries(meshList)
    print("Done")

def plotBySegment(col):
    meshList = list()
    for seg in col.segList:
        if seg.isSL():
            meshList.append(createSphere(seg.getProximal(), seg.getRadius(), color = np.array([1,0,1])))
        elif seg.getAncestor().getType() == "Inlet":
            meshList.append(createCylinder(seg.getProximal(),seg.getDistal(),seg.getRadius(), np.array([1,0,0])))
        else: # seg.getAncestor().getType() == "Outlet"
            meshList.append(createCylinder(seg.getProximal(),seg.getDistal(),seg.getRadius(), np.array([0,0,1])))
    
    o3d.visualization.draw_geometries(meshList)

def plot2DSegment(col):
    pass

##### Helper Functions #####

def createConeAdv(baseCenter, topCenter, baseRadius, topRadius, baseAxis, topAxis, stackCount, sliceCount):
    # Functionally the same as createCone() but also accounts for the normal direction 
    # of each cone face and interpolates intermediate stack normals

    # Check slice and stack count is 2 or more
    if sliceCount < 2:
        raise ValueError("Slice count must be 2 or higher")
    if stackCount < 2:
        raise ValueError("Stack count must be 2 or higher")

    V = list()
    F = list()

    axis = topCenter - baseCenter
    coneLength = np.linalg.norm(axis)
    axis = axis/coneLength
    
    dTheta = 2 * np.pi / sliceCount
    # Normalize base and top axes
    baseAxis = baseAxis/np.linalg.norm(baseAxis)
    topAxis = topAxis/np.linalg.norm(topAxis)
    # Make sure random ray is orthogonal to the two axes
    randomRay = np.cross(baseAxis,topAxis)

    # Make the vertices
    for i in range(stackCount):
        # Compute center of the stack
        center = baseCenter + i/stackCount * coneLength * axis
        # Interpolate the radius based on stack
        stackRadius = baseRadius + (topRadius - baseRadius) * i/(stackCount-1)
        print(stackRadius)
        # Determine stack axis by interpolating between base and stack axes
        stackAxis = (1 - i/(stackCount-1)) * baseAxis + i/(stackCount-1) * topAxis
        stackAxis = stackAxis/np.linalg.norm(stackAxis)
        print(stackAxis)
        # Rotate around stack generating vertex for each slice
        for j in range(sliceCount):
            # TODO inconsistent rotation of face because of np.cross(np.cross())
            V.append(center + stackRadius * (np.cos(dTheta*j)*np.cross(stackAxis,randomRay) + np.sin(dTheta*j)*np.cross(np.cross(stackAxis,randomRay),axis)))
            # Assign the indices
            if i < stackCount-1 and j < sliceCount-1:
                # F.append([sliceCount*i+j, sliceCount*i+j+1, sliceCount*(i+1)+j+1])
                # F.append([sliceCount*i+j, sliceCount*(i+1)+j+1, sliceCount*(i+1)+j])
                F.append([sliceCount*i+j, sliceCount*(i+1)+j+1, sliceCount*i+j+1])
                F.append([sliceCount*i+j, sliceCount*(i+1)+j, sliceCount*(i+1)+j+1])
        if i < stackCount - 1:
            # F.append([sliceCount*(i+1)-1, sliceCount*i, sliceCount*(i+1)])
            # F.append([sliceCount*(i+1)-1, sliceCount*(i+1), sliceCount*(i+2)-1])
            F.append([sliceCount*(i+1)-1, sliceCount*(i+1), sliceCount*i])
            F.append([sliceCount*(i+1)-1, sliceCount*(i+2)-1, sliceCount*(i+1)])

    # Make the mesh
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(V)
    mesh.triangles = o3d.utility.Vector3iVector(F)
    mesh.compute_triangle_normals()
    mesh.compute_vertex_normals()

    return mesh

def createSphere(loc, radius, color):
    # Make template sphere
    sph = o3d.geometry.TriangleMesh.create_sphere(radius)
    # Translarte to correct position
    sph = copy.deepcopy(sph).translate(loc)
    # Make it look nice
    sph.paint_uniform_color(color)
    sph.compute_triangle_normals()
    sph.compute_vertex_normals()
    return sph

def createCylinder(a, b, radius, color):
    # Make template cylinder aligned with z-axis
    height = np.linalg.norm(b-a)
    cyl = o3d.geometry.TriangleMesh.create_cylinder(radius, height)
    # Need to rotate and translate the cylinder to its correct location
    # Compute axis and angle between vectors
    vec1 = np.array([0,0,1])
    vec2 = (b - a) / np.linalg.norm(b-a)
    axis = np.cross(vec1, vec2) / np.linalg.norm(np.cross(vec1, vec2))
    angle = np.arccos((np.dot(vec1, vec2)))
    # Apply rotation and translation
    R = cyl.get_rotation_matrix_from_axis_angle(axis * angle)
    cyl = copy.deepcopy(cyl).rotate(R, center = (0,0,0))
    cyl = copy.deepcopy(cyl).translate((a+b)/2)
    # Make it look nice
    cyl.paint_uniform_color(color)
    cyl.compute_triangle_normals()
    cyl.compute_vertex_normals()
    return cyl

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
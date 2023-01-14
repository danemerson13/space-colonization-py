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

##### Helper Functions #####

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
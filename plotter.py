import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from stl import mesh as np_mesh
from mpl_toolkits import mplot3d

def plotByBranch(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            if branch.getType() == "Inlet":
                ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'b-',lw = 4)
                if branch.isSL():
                    ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'o', color = 'magenta')
            if branch.getType() == "Outlet":
                ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'b-')
                if branch.isSL():
                    ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'o', color = 'magenta')

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    plt.show()

def plotPressure(col):
    # Convert from dyn/mm^2 to mmHg with P/13.3322
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min, max = getPressureLims(col)
    # Put pressure back in mmHg
    min = min/13.3322; max = max/13.3322
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'k-')

    for node in col.nodeList:
        if node.getType() != "Interior":
            nodeLoc = node.getLocation()
            ax.plot3D(nodeLoc[0], nodeLoc[1], nodeLoc[2], 'o', color = mymap.to_rgba(node.getPressure()/13.3322)[:-1])

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Pressure (mmHg)', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.savefig('pressure.png')
    plt.show()

def plotFlowRate(col):
    # Convert from mm3/s to mL/min by taking Q * 3/50 ... (same as Q/16.66)
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min, max = getFlowRateLims(col)
    min = min * 3/50; max = max * 3/50
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        color = mymap.to_rgba(branch.getFlowRate() * 3/50)[:-1]
        rad = branch.getRadius()
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Flow Rate (mL/min)', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.savefig('flowrate.png')
    plt.show()

def plotReynolds(col):
    # Convert from mm3/s to mL/min by taking Q * 3/50 ... (same as Q/16.66)
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min, max = getReynoldsLims(col)
    min = min; max = max
    cmap = matplotlib.cm.get_cmap('viridis')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        color = mymap.to_rgba(branch.getReynolds())[:-1]
        rad = branch.getRadius()
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Reynolds Number', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.savefig('reynolds.png')
    plt.show()

def plotFilling(col,time):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col.branchList:
        length = col.branchLength(branch)
        cumLength = 0
        coords = getBranchCoords(branch)
        filledLen = branch.percentFull() * length

        for j in range(len(coords) - 1):
            # Update the length to see how much we should color as filled
            cumLength += np.linalg.norm(coords[j+1] - coords[j])
            # If the cumLen is shorter than what is filled, color red
            if cumLength < filledLen or (np.abs(cumLength-filledLen) < 1e-6):
                ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', c = np.array([1,0,0]), lw = branch.getRadius())
            elif ((cumLength - filledLen) < np.linalg.norm(coords[j+1] - coords[j])):
                diff = cumLength - filledLen
                vec = coords[j+1] - coords[j]
                vecLen = np.linalg.norm(vec)
                unitVec = vec/vecLen
                step = coords[j] + unitVec * (diff/vecLen)
                ax.plot3D([coords[j][0],step[0]],[coords[j][1],step[1]],[coords[j][2],step[2]],'-', c = np.array([1,0,0]), lw = branch.getRadius())
                ax.plot3D([step[0],coords[j+1][0]],[step[1],coords[j+1][1]],[step[2],coords[j+1][2]],'-', c = np.array([0,0,1]), lw = branch.getRadius())
            else:
                ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', c = np.array([0,0,1]), lw = branch.getRadius())

        # for j in range(len(coords) - 1):
        #     # Update the length to see how much we should color as filled
        #     cumLength += np.linalg.norm(coords[j+1] - coords[j])
        #     # If the cumLen is shorter than what is filled, color red
        #     if cumLength < filledLen:
        #         ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', c = np.array([1,0,0]), lw = branch.getRadius())
        #     elif cumLength >= filledLen:
        #         diff = cumLength - filledLen
        #         vec = coords[j+1] - coords[j]
        #         vecLen = np.linalg.norm(vec)
        #         unitVec = vec/vecLen
        #         step = coords[j] + unitVec * (diff/vecLen)
        #         ax.plot3D([coords[j][0],step[0]],[coords[j][1],step[1]],[coords[j][2],step[2]],'-', c = np.array([1,0,0]), lw = branch.getRadius())
        #         ax.plot3D([step[0],coords[j+1][0]],[step[1],coords[j+1][1]],[step[2],coords[j+1][2]],'-', c = np.array([0,0,1]), lw = branch.getRadius())
            # else:
            #     ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', c = np.array([0,0,1]), lw = branch.getRadius())
            
    bound_mesh(ax, mesh)
    ax.set_title("T = %.2f" %time)
    ax.view_init(elev = -180., azim = -90.)
    ax.text(25,0,175,'Percent Full: %.0f' %(col.percentFull() * 100))
    ax.set_axis_off()
    plt.savefig("fill_%.2f.png" %time, bbox_inches='tight')
    plt.close(fig)

def plotRadiiHist(col):
    radii = list() # Branch diameter in um (rad * 2 * 1000)
    for branch in col.branchList:
        radii.append(branch.getRadius() * 2 * 1000)
    plt.figure()
    plt.hist(radii, bins = 10)
    plt.show()

def plotGenerations(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    colorList = list(["red","green","blue","orange","violet","turquoise","darkgrey"])
    # colorList = list()
    for branch in col.branchList:
        gen = branch.getGeneration()
        if gen > len(colorList):
            colorList.append(np.random.rand(3))
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-',color = colorList[gen-1],lw = branch.getRadius())

    bound_mesh(ax, mesh)
    ax.view_init(elev = -180., azim = -90.)
    ax.set_axis_off()
    plt.show()



##### HELPER FUNCTIONS #####

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

def getLiverSTL():
    # Load the STL files and add the vectors to the plot
    mesh = np_mesh.Mesh.from_file('Models/liver159.stl')
    liver = mplot3d.art3d.Poly3DCollection(mesh.vectors)
    liver.set_alpha(0.05)
    liver.set_facecolor([0.5, 0.5, 0.5, 0.25])
    liver.set_edgecolor([0., 0., 0., 0.25])

    return liver, mesh

def bound_mesh(ax, mesh):
    pts = mesh.points.reshape(len(mesh.points)*3,3)

    x_range = abs(max(pts[:,0]) - min(pts[:,0]))
    x_middle = np.mean([max(pts[:,0]), min(pts[:,0])])
    y_range = abs(max(pts[:,1]) - min(pts[:,1]))
    y_middle = np.mean([max(pts[:,1]), min(pts[:,1])])
    z_range = abs(max(pts[:,2]) - min(pts[:,2]))
    z_middle = np.mean([max(pts[:,2]), min(pts[:,2])])

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def getPressureLims(col):
    max = -np.inf
    min = np.inf
    for node in col.nodeList:
        if node.getType() != "Interior":
            PNode = node.getPressure()
            if PNode > max:
                max = PNode
            if PNode < min:
                min = PNode
    return min, max

def getFlowRateLims(col):
    max = -np.inf
    min = np.inf
    for branch in col.branchList:
        QBranch = branch.getFlowRate()
        if QBranch > max:
            max = QBranch
        if QBranch < min:
            min = QBranch
    return min, max

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
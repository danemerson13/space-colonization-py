import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from stl import mesh as np_mesh
from mpl_toolkits import mplot3d
import colony

def getLiverSTL():
    # Load the STL files and add the vectors to the plot
    mesh = np_mesh.Mesh.from_file('Models/liver.stl')
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

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
    ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def plotNodeList(col):
    print("Hello")
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for i in range(1, len(col.nodeList)):
        prox = col.nodeList[i-1].getLocation()
        dist = col.nodeList[i].getLocation()
        print("Norm: ", np.linalg.norm(prox - dist))
        
        ax.plot3D([prox[0],dist[0]],[prox[1],dist[1]],[prox[2],dist[2]],'b-')

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col.D, col.dk, col.di))
    ax.view_init(elev=25., azim=-55.)
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col.D, col.dk, col.di))
    plt.show()

def plotByBranch(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        color = np.random.rand(3)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', c = color)

        if branch.isRoot():
            ax.plot3D(coords[0][0],coords[0][1],coords[0][2],'go')
        elif branch.isTerminal():
            ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'ro')

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col.D, col.dk, col.di))
    ax.view_init(elev=25., azim=-55.)
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col.D, col.dk, col.di))
    plt.show()

def plotRadius(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        color = np.random.rand(3)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'b-', lw = branch.getRadius())

        if branch.isRoot():
            ax.plot3D(coords[0][0],coords[0][1],coords[0][2],'go')
        elif branch.isTerminal():
            ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'ro')

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col.D, col.dk, col.di))
    ax.view_init(elev=25., azim=-55.)
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col.D, col.dk, col.di))
    plt.show()

def plotPressure(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min, max = getPressureLims(col)
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'k-')

    for node in col.nodeList:
        if node.isFurcation() or node.isTerminal() or node.isRoot():
            nodeLoc = node.getLocation()
            ax.plot3D(nodeLoc[0], nodeLoc[1], nodeLoc[2], 'o', color = mymap.to_rgba(node.getPressure())[:-1])

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col.D, col.dk, col.di))
    ax.view_init(elev=25., azim=-55.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Pressure', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col.D, col.dk, col.di))
    plt.show()

def plotFlowRate(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min, max = getFlowRateLims(col)
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col.branchList:
        coords = getBranchCoords(branch)
        color = mymap.to_rgba(branch.getFlowRate())[:-1]
        rad = branch.getRadius()
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col.D, col.dk, col.di))
    ax.view_init(elev=0., azim=-90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Flow Rate', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col.D, col.dk, col.di))
    plt.show()

def plotGeneration(col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Get generations
    gens = list()
    for branch in col.branchList:
        gens.append(branch.getGeneration())
    gens = list(set(gens)) # Make the list unique

    for gen in gens:
        color = np.random.rand(3)
        for branch in col.branchList:
            if branch.getGeneration() == gen:
                coords = getBranchCoords(branch)
                rad = branch.getRadius()
                for j in range(len(coords) - 1):
                    ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col.D, col.dk, col.di))
    ax.view_init(elev=0., azim=-90.)
    ax.set_axis_off()
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col.D, col.dk, col.di))
    plt.show()

def plotFlowRate2x(col1, col2):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min1, max1 = getFlowRateLims(col1); min2, max2 = getFlowRateLims(col2)
    min = np.min([min1, min2]); max = np.max([max1, max2])
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col1.branchList:
        coords = getBranchCoords(branch)
        color = mymap.to_rgba(branch.getFlowRate())[:-1]
        rad = branch.getRadius()
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)
    for branch in col2.branchList:
        coords = getBranchCoords(branch)
        color = mymap.to_rgba(branch.getFlowRate())[:-1]
        rad = branch.getRadius()
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', lw = rad, color = color)

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col1.D, col1.dk, col1.di))
    ax.view_init(elev=0., azim=-90.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Flow Rate', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col1.D, col1.dk, col1.di))
    plt.show()

def plotPressure2x(col1, col2):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    # Set up cmap
    min1, max1 = getPressureLims2x(col1); min2, max2 = getPressureLims2x(col2)
    min = np.min([min1, min2]); max = np.max([max1, max2])
    cmap = matplotlib.cm.get_cmap('cool')
    norm = matplotlib.colors.Normalize(vmin = min, vmax = max)
    mymap = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)

    for branch in col1.branchList:
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'k-')
    for node in col1.nodeList:
        if node.isFurcation() or node.isTerminal() or node.isRoot():
            nodeLoc = node.getLocation()
            ax.plot3D(nodeLoc[0], nodeLoc[1], nodeLoc[2], 'o', color = mymap.to_rgba(node.getPressure())[:-1])

    for branch in col2.branchList:
        coords = getBranchCoords(branch)
        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'k-')
    for node in col2.nodeList:
        if node.isFurcation() or node.isRoot():
            nodeLoc = node.getLocation()
            ax.plot3D(nodeLoc[0], nodeLoc[1], nodeLoc[2], 'o', color = mymap.to_rgba(node.getPressure())[:-1])

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col1.D, col1.dk, col1.di))
    ax.view_init(elev=25., azim=-55.)
    ax.set_axis_off()
    fig.colorbar(mymap,  label = 'Pressure', orientation = 'vertical', pad = -0.1, shrink = 0.5)
    plt.savefig("spaceColonization_%d_%d_%d.png" % (col1.D, col1.dk, col1.di))
    plt.show()

def getFlowRateLims(col):
    max = -np.inf
    min = np.inf
    for branch in col.branchList:
        if branch.getFlowRate() > max:
            max = branch.getFlowRate()
        if branch.getFlowRate() < min:
            min = branch.getFlowRate()
    return min, max

def getPressureLims2x(col):
    max = -np.inf
    min = np.inf
    for node in col.nodeList:
        if node.isFurcation() or node.isRoot():
            if node.getPressure() > max:
                max = node.getPressure()
            if node.getPressure() < min:
                min = node.getPressure()
    return min, max

def getPressureLims(col):
    max = -np.inf
    min = np.inf
    for node in col.nodeList:
        if node.isFurcation() or node.isTerminal() or node.isRoot():
            if node.getPressure() > max:
                max = node.getPressure()
            if node.getPressure() < min:
                min = node.getPressure()
    return min, max

def plotByBranchMulti(init, col):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'b-')

        if branch.isRoot():
            ax.plot3D(coords[0][0],coords[0][1],coords[0][2],'go')
        elif branch.isTerminal():
            ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'ro')

    for branch in init.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'c-')

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col.D, col.dk, col.di))
    # Views
    # ax.view_init(elev=25., azim=-55.) # iso
    # ax.view_init(elev = -90., azim = -90.) # top
    ax.view_init(elev = 0., azim = -90.) # front
    # ax.view_init(elev = -180., azim = -180.) # right

    plt.savefig("spaceColonization_%d_%d_%d.png" % (col.D, col.dk, col.di))
    plt.show()

def plotByBranchDouble(init1, col1, init2, col2):
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    liver, mesh = getLiverSTL()
    ax.add_collection3d(liver)

    for branch in col1.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'b-')

        if branch.isRoot():
            ax.plot3D(coords[0][0],coords[0][1],coords[0][2],'go')
        elif branch.isTerminal():
            ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'ro')

    for branch in init1.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'c-')

    for branch in col2.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'r-')

        if branch.isRoot():
            ax.plot3D(coords[0][0],coords[0][1],coords[0][2],'go')
        elif branch.isTerminal():
            ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'mo')

    for branch in init2.branchList:
        coords = getBranchCoords(branch)

        for j in range(len(coords) - 1):
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'c-')

    bound_mesh(ax, mesh)
    ax.set_title("D = %d, dk = %d, di = %d" % (col1.D, col1.dk, col1.di))
    # Views
    # ax.view_init(elev=25., azim=-55.) # iso
    # ax.view_init(elev = -90., azim = -90.) # top
    ax.view_init(elev = 0., azim = -90.) # front
    # ax.view_init(elev = -180., azim = -180.) # right

    plt.savefig("spaceColonization_%d_%d_%d.png" % (col1.D, col1.dk, col1.di))
    plt.show()

def getBranchCoords(branch):
    coordList = list()
    curNode = branch.getDistal()
    coordList.append(curNode.getLocation())
    while curNode != branch.getProximal():
        curNode = curNode.getParent()
        coordList.append(curNode.getLocation())
    coordList.reverse()
    return coordList

def plotPressureDistribution(col):
    pressureDist = list()
    for node in col.nodeList:
        if node.isTerminal():
            pressureDist.append(node.getPressure())
            print('Pressure = %.2f' %pressureDist[-1])
    plt.hist(pressureDist, 10, (0,10))
    plt.axvline(np.asarray(pressureDist).mean(), color='k', linestyle='dashed', linewidth=1)
    plt.xlabel('Pressure')
    plt.ylabel('Frequency')
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

        for j in range(len(coords) - 1):
            # Update the length to see how much we should color as filled
            cumLength += np.linalg.norm(coords[j+1] - coords[j])
            if cumLength/length <= branch.percentFull():
                color = np.array([1,0,0])
            else:
                color = np.array([0,0,1])
            # Now plot
            ax.plot3D([coords[j][0],coords[j+1][0]],[coords[j][1],coords[j+1][1]],[coords[j][2],coords[j+1][2]],'-', c = color, lw = branch.getRadius())

        if branch.isTerminal():
            ax.plot3D(coords[-1][0],coords[-1][1],coords[-1][2],'ko')

    bound_mesh(ax, mesh)
    ax.set_title("T = %.2f" %time)
    ax.view_init(elev=25., azim=-55.)
    plt.savefig("fill_%.2f.png" %time, bbox_inches='tight')
    plt.close(fig)
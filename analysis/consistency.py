import pickle
import sys
import os

# Adding parent directory to path
sys.path.append(os.path.abspath(os.getcwd()))

from util import visualizer

def main():
    # Generate plots for 3 models
    nSL = 10
    # nSL = 100
    # nSL = 1000
    filename = '../Space Colonization/results/nSLConvergence/' + str(nSL) + 'SL' + '/model.pkl'
    with open(filename, 'rb') as handle:
        col = pickle.load(handle)

    stackCount = 20 # Set to 10 for 1000SL model
    sliceCount = 40 # Set to 20 for 1000SL model

    meshList = visualizer.skinTree(col, stackCount, sliceCount)

    visualizer.CustomVisualizer(meshList, -30)

if __name__ == "__main__":
    main()
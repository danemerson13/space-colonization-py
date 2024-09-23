import pickle
import time
import os
import sys

from src import colony

'''
In this script we perform a convergence study on the time step, dt, used during the filling of the model
We use a fixed model of 100SLs

We will look at the convergence of the concentration of the total model concentration, as well as the time for the last SL to fill
'''

def main():
    sys.setrecursionlimit(10**9)
    dts = [1e1, 5e0, 1e0, 5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
    # filename = 'results/100SL/sample0/model.pkl'
    filename = 'results/nSLConvergence/1000SL/model.pkl'
    for dt in dts:
        print('Simulating filling with dt = %.4f' %dt)
        # Load the original model
        with open(filename, 'rb') as handle:
            col = pickle.load(handle)
        # Clear the vessels
        col.clearConcentrations()
        start = time.time()
        col.fillTree(dt = dt, tStop = 200, gif = False)
        end = time.time()
        col.tFillSolve = end - start
        print('Filling simulation took %.3f secs' %(end-start))
        # Save the model
        # folder = 'results/tConvergence/dt' + str(dt)
        folder = 'results/tConvergenceFine/dt' + str(dt)
        if not os.path.exists(folder):
            os.mkdir(folder)
        col.saveModel(folder)
        print('Model Saved')

if __name__ == "__main__":
    main()
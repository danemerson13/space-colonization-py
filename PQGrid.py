import pickle
import time
import os
import sys

from src import colony

'''
In this script we perform a sweep over varying P and Q values on a fixed model
'''

def main():
    sys.setrecursionlimit(10**9)
    pressures = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    # flows = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    flows = [100]

    filename = 'results/100SL/sample0/model.pkl'
    for pressure in pressures:
        for flow in flows:
            print('Simulating P = %d mmHg, Q = %d mL/min' %(pressure, flow))
            # Load the base model
            with open(filename, 'rb') as handle:
                col = pickle.load(handle)
            # Clear the vessels concentrations
            col.clearConcentrations()
            start = time.time()
            # Run model for specific pressure and flow
            col.solveRSL(Qactual = flow, Pin = pressure, Pout = 0, n_iter=1000, tol=1e-6, verbosity=0, a=0, b=1e-4)
            col.setSLVolume(375000) # Setting SL volume with bloodVol = 375 mL
            end = time.time()
            print('Solving pressures and flows took %.3f secs' %(end - start))
            # Simulate filling
            start = time.time()
            col.fillTree(dt = 0.01, tStop = 1000, gif = False)
            end = time.time()
            print("Filling model took %.3f secs" %(end - start))
            # Save the model
            folder = 'results/PQGrid/P%dQ%d/' %(pressure, flow)
            if not os.path.exists(folder):
                os.mkdir(folder)
            col.saveModel(folder)
            print('Model Saved')

if __name__ == "__main__":
    main()


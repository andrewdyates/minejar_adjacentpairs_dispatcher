#!/usr/bin/python
import numpy as np
import sys
from scipy.spatial import distance

def main(npy_fname):
  M = np.load(npy_fname)
  Q = distance.squareform(M)
  Q += np.eye(np.size(Q,0))
  np.savetxt(npy_fname+".txt", Q, fmt="%.5f")


if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))

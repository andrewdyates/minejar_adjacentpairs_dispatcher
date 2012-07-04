#!/usr/bin/python
"""Save a copy of a matrix with all squared values."""
import sys
import numpy as np

def main():
  m_fname = sys.argv[1]
  print "Squaring %s..." % m_fname
  M = np.load(m_fname)
  M = M**2
  msq_fname = m_fname+".squared.npy"
  np.save(msq_fname, M)
  print "Saved squared at %s" % msq_fname

if __name__ == "__main__":
  main()

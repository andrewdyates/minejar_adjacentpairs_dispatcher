#!/usr/bin/python
import numpy as np
import sys
m_fname = sys.argv[1]
print "Loading %s..." % m_fname
M = np.load(sys.argv[1])
print "Argsorting %s..." % m_fname
Q = np.argsort(M)
m_fname_argsort = m_fname + ".argsorted.npy"
print "Saving argsorted as %s..." % m_fname_argsort
np.save(m_fname_argsort, Q)

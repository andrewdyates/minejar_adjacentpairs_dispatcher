#!/usr/bin/python
"""Merge numpy matrices. Use any 'True' value as merged value.

EXAMPLE USE:
python merge.py GSE2034.GPL96.eQTL.tab.compiled.values.npy GSE2034.GPL96.eQTL.tab.missing.values.npy GSE2034.GPL96.eQTL.tab.low25.values.npy
"""

import sys
import numpy as np
from py_symmetric_matrix import *

def main(m_files):
    M = np.load(m_files[0])
    for m_file in m_files[1:]:
        Q = np.load(m_file)
        for i, v in enumerate(Q):
            if v:
                M[i] = v
    np.save(m_files[0]+".merged.%d.npy" % len(m_files), M)

if __name__ == "__main__":
    main(sys.argv[1:])

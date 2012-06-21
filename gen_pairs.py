#!/usr/bin/python
"""From a boolean matrix of missing values, print pairs of missing variables.
varlist_fname order corresponds with bool_m_fname

USE:
python gen_pairs.py bool_m_fname, varlist_fname, tab_fname > missing_vars.tab

EXAMPLE:
python gen_pairs.py GSE2034.GPL96.eQTL.tab.compiled.bools.npy GSE2034.GPL96.eQTL.tab.varlist.txt GSE2034.GPL96.eQTL.tab
"""
import sys
import numpy as np
from py_symmetric_matrix import *

def main(bool_m_fname, varlist_fname, tab_fname):
  # load bool matrix
  B = np.load(bool_m_fname)
  # read tabfile as {varname => str}
  varlist = [s.strip() for s in open(varlist_fname)]
  M = {}
  fp = open(tab_fname)
  for line in fp:
    if line[0] == '#' or not line.strip():
      continue
    name,c,row = line.partition('\t')
    M[name] = line.strip('\n')
  fp.close()
  # assert that all variables have been accounted
  assert not (set(varlist) - set(M.keys()))
    
  # for each 0 in B, print pair of corresponding vars
  for i, v in enumerate(B):
    if v == 0:
      x, y = inv_sym_idx(i, len(varlist))
      print M[varlist[x]]
      print M[varlist[y]]

if __name__ == "__main__":
  main(*sys.argv[1:])

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
  M = {}
  fp = open(tab_fname)
  for line in fp:
    if line[0] in ('#' or '\n'):
      continue
    name,c,row = line.partition('\t')
    M[name] = line.strip('\n')
  fp.close()

  # Read varlist
  varlist = []
  fp = open(varlist_fname)
  for line in fp:
    if line[0] in ('#' or '\n'):
      continue
    name,c,row = line.partition('\t')
    varlist.append(name.strip('\n'))
  fp.close()

  n = len(varlist)
  # assert that all variables have been accounted
  assert not (set(varlist) - set(M.keys()))
  print (n*n-1)/2, len(B), (n*n-1)/2 - len(B)
  assert (n*n-1)/2 == len(B)
    
  # for each 0 in B, print pair of corresponding vars
  varset = set()
  for i in np.where(B == 0)[0]:
    x, y = inv_sym_idx(i, n)
    print M[varlist[x]]
    print M[varlist[y]]
    varset.add(varlist[x])
    varset.add(varlist[y])

  # Write list of variable names with at least one missing row.
  sys.stderr.write("%d variables with at least one missing pair." % len(varset))
  sys.stderr.write('\n'.join(varset))
  sys.stderr.write('\n')
  

if __name__ == "__main__":
  main(*sys.argv[1:])

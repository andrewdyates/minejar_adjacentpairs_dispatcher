#!/usr/bin/python
"""Compile MINE.jar results into a squareform numpy dependency matrix.

Revised June 20th, 2012
USG:

python script.py tab_fname.tab minefile1.out minefile2.out ...
"""
import sys
import numpy as np
from py_symmetric_matrix import *


# MIC value
V_COL = 2

def tab_to_varlist(tab_fname):
  """Transform a .tab gene expression value matrix into a list of variables.

  Args:
    tab_fname: str of path to .tab MINE.jar ready matrix file
  Returns:
    [str] of cleaned variable names IN FILE ORDER
  """ 
  fp = open(tab_fname)
  varlist = []
  for line in fp:
    if line[0] == "#": continue
    var,c,cc = line.partition('\t')
    varlist.append(var)
  fp.close()
  return varlist

def populate_matrix(varlist, fnames):
  # values
  V = NamedSymmetricMatrix(varlist, store_diagonal=False)
  # existence
  B = NamedSymmetricMatrix(varlist, store_diagonal=False, dtype=np.bool)
  n_set = 0
  for fname in fnames:
    fp = open(fname, 'r')
    for line in fp:
      row = line.strip().split(',')
      if not row or row[0] == "X var":
        continue
      x, y, v = row[0], row[1], np.float(row[V_COL])
      if not B.get(x,y):
        V.set(x,y,v)
        B.set(x,y,1)
        n_set += 1
    fp.close()
  return V, B, n_set

    

def main(tab_fname, minefiles):
  varlist = tab_to_varlist(tab_fname)
  V, B, n_set = populate_matrix(varlist, minefiles)
  # save value matrix to disk
  np.save("%s.compiled.values.npy" % tab_fname, V._m)
  np.save("%s.compiled.bools.npy" % tab_fname, B._m)
  fp = open("%s.varlist.txt" % tab_fname, 'w')
  for name in varlist:
    fp.write("%s\n" % name)
  print "Set %d of %d variable pair values (%.2f%%)." %\
    (n_set, len(V), float(n_set)/len(V)*100)


if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2:])

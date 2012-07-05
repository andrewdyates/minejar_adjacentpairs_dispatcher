#!/usr/bin/python
"""Batchfile for computing a list of pairwise MINE.jar rows.

Assume that tab_fname variable order corresponds to Boolean matrix.

SAMPLE USE:
python batch_mine_pairwise.py tab_fname=/fs/lustre/osu6683/GSE2034.GPL96.eQTL.normed.tab bool_fname=/fs/lustre/osu6683/GSE2034.GPL96.eQTL.normed.tab.varlist.txt.MINE.isset.npy offset=0 k=10000 
"""
import sys
import subprocess
import numpy as np
import os
from py_symmetric_matrix import *


def main(tab_fname=None, bool_fname=None, offset=0, k=200000, work_dir='/fs/lustre/osu6683', minejar_path="/fs/lustre/osu6683/MINE.jar"):
  offset = int(offset)
  k = int(k)
  assert tab_fname
  assert bool_fname
  assert k
  assert k > 0 and offset >= 0

  # read tabfile as {varname => str}
  M = {}
  varlist = []
  fp = open(tab_fname)
  for line in fp:
    if line[0] in ('#' or '\n'):
      continue
    name,c,row = line.partition('\t')
    M[name] = line.strip('\n')
    varlist.append(name)
  fp.close()
  n = len(M)
  B = np.load(bool_fname)
  assert n*(n-1)/2 == len(B)
  
  Q = np.where(B == 0)[0]
  if offset >= len(Q):
    print "Offset %d is greater than number of missing pairs %d. Exiting." % \
      (offset, len(Q))
    return 
  if k+offset > len(Q):
    k = len(Q)-offset
    print "k+offset greater than number of pairs. Setting k to %d." % k

  fname = os.path.basename(tab_fname).rpartition('.')[0] + "_offset_%d_%d.tab"%\
     (offset, k)
  tmp_filename = os.path.join(work_dir, fname)
  print "Generating pairwise .tab file '%s' for pairs %d through %d..." % \
  (tmp_filename, offset+1, offset+k)

  fp = open(tmp_filename, 'w')
  for i in Q[offset:offset+k]:
    x, y = inv_sym_idx(i, n)
    fp.write(M[varlist[x]] + "\n")
    fp.write(M[varlist[y]] + "\n")
  fp.close()

  cmd = 'java -jar %s "%s" -adjacentPairs' % (minejar_path, tmp_filename)
  print "Running MINE.jar job on %s..." % tmp_filename
  print cmd
  subprocess.call(cmd, shell=True)

  print "Deleting temporary input file %s..." % tmp_filename
  os.remove(tmp_filename)
  print "Complete."
  

if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))

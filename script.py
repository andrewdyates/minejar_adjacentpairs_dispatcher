#!/usr/bin/python
"""Compile a list of MINE.jar results into a squareform numpy dependency matrix.

USE:
  python script.py path/to/outdir path/to/tab_fname.tab minefile1.out minefile2.out ...
EXAMPLE:
  python script.py ~/gse2034 GSE2034.GPL96.eQTL.normed.tab.varlist.txt ~/gse2034/gse2034_mine/*Results.csv

MINE.jar headers in form like:
X var,Y var,MIC (strength),MIC-p^2 (nonlinearity),MAS (non-monotonicity),MEV (functionality),MCN (complexity),Linear regression (p)
"""
import sys
import numpy as np
import errno, os
from py_symmetric_matrix import *


# names of matrices in column order
M_NAMES = ['MIC', 'NONLIN', 'MAS', 'MEV', 'MCN', 'P']


def tab_to_varlist(tab_fname):
  """Transform a .tab gene expression value matrix into a list of variables.
  Also works with list of variable names, one per line
  
  Args:
    tab_fname: str of path to .tab MINE.jar ready matrix file
  Returns:
    [str] of cleaned variable names IN FILE ORDER
  """ 
  fp = open(tab_fname)
  varlist = []
  for line in fp:
    if line[0] == "#" or line[0] == '\n': continue
    var,c,cc = line.partition('\t')
    varlist.append(var.strip('\n'))
  fp.close()
  return varlist

def populate_matrix(fp, M_Dict, B):
  """Populate set of MINE matrices given an open file of results.

  Args:
    fp: [*str] of new, open file pointer to MINE.jar results file
    M_Dict: {str: `NamedSymmetricMatrix`} of matrices; str is name in M_NAMES
    B: `NamedSymmetricMatrix` of bool matrix indicating if a value has been set
  Returns:
    (int, int) of num pairs set and skipped respectively
  """
  n_set = 0
  n_dupe = 0
  for line in fp:
    row = line.strip().split(',')
    # skip blanks or header lines
    if not row or row[0] == "X var":
      continue
    x, y = row[0], row[1]
    v = map(np.float, row[2:])
    if not B.get(x,y):
      n_set += 1
      B.set(x,y,1)
      for i, m_name in enumerate(M_NAMES):
        M_Dict[m_name].set(x,y,v[i])
    else:
      n_dupe += 1

  return n_set, n_dupe
    

def main(outdir, tab_fname, minefiles):
  """Compile list of minefiles into a matrix and save it to disk."""
  assert outdir
  assert tab_fname
  assert minefiles
  assert len(minefiles) > 0 and type(minefiles) != str
  outdir = os.path.expanduser(outdir)
  
  varlist = tab_to_varlist(tab_fname)
  m = len(varlist)
  n_expected = m*(m-1)/2
  print "m=%d; expected num pairs=%d" % (m, n_expected)

  # for each matrix name, assign it a value matrix
  M_Dict = {}
  for name in M_NAMES:
    M_Dict[name] = NamedSymmetricMatrix(varlist, store_diagonal=False, dtype=np.float)
  # existence (boolean) matrix
  B = NamedSymmetricMatrix(varlist, store_diagonal=False, dtype=np.bool)

  # For each MINE.jar results file, load into the matrices
  n_set_total = 0
  n_dupe_total = 0
  for fname in minefiles:
    fp = open(fname, 'r')
    n_set, n_dupe = populate_matrix(fp, M_Dict, B)
    fp.close()
    print "Set %d pairs from MINE.jar results file %s. Skipped %d dupes." % \
      (n_set, fname, n_dupe)
    n_set_total += n_set
    n_dupe_total += n_dupe
  print "===================="
  print "Completed reading pairs from MINE.jar results"
  print "%d total pairs set (%d dupes skipped), expected %d." % \
    (n_set_total, n_dupe_total, n_expected)
  
    
  # Save matrices to disk.
  # Save only the underlying numpy matrix of the NamedMatrix. (._m)
  print "Saving matrices to disk..."
  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise
    
  for m_name in M_NAMES:
    m_fname = "%s.%s.npy" % (tab_fname, m_name)
    np.save(os.path.join(outdir, m_fname), M_Dict[m_name]._m)
  m_fname = "%s.MINE.isset.npy" % tab_fname
  np.save(os.path.join(outdir, m_fname), B._m)

  # Write variable list to disk.
  fp = open(os.path.join(outdir, "%s.MINE.varlist.txt" % tab_fname), 'w')
  for name in varlist:
    fp.write("%s\n" % name)


if __name__ == "__main__":
  main(outdir=sys.argv[1], tab_fname=sys.argv[2], minefiles=sys.argv[3:])

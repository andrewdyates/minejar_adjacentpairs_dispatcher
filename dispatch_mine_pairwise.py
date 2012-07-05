#!/usr/bin/python
"""Dispatch batches of pairwise MINE.jar jobs from a list of row pairs.


USE:
python dispatch_mine_pairwise.py tab_fname=/fs/lustre/osu6683/GSE2034.GPL96.eQTL.normed.tab bool_fname=/fs/lustre/osu6683/GSE2034.GPL96.eQTL.normed.tab.varlist.txt.MINE.isset.npy

BATCH SCRIPT USE:
  python batch_mine_pairwise.py tab_fname=/fs/lustre/osu6683/GSE2034.GPL96.eQTL.normed.tab bool_fname=/fs/lustre/osu6683/GSE2034.GPL96.eQTL.normed.tab.varlist.txt.MINE.isset.npy offset=0 k=10000 

Assume that tab file and bool matrix files are already in the temporary computing directory.
"""
import sys
import numpy as np
import subprocess
import os

TEMPLATE = \
"""#PBS -N %(jobname)s
#PBS -l nodes=%(n_nodes)d:ppn=%(n_ppn)d
#PBS -j oe
#PBS -S /bin/bash
#PBS -l walltime=%(walltime)s
#tdate=$(date +%%T)

set -x
cd /nfs/01/osu6683/
source .bash_profile
mpiexec parallel-command-processor %(dispatch_script)s
"""

CMD = "python %(script_path)s/batch_mine_pairwise.py tab_fname=%(tab_fname)s bool_fname=%(bool_fname)s offset=%(offset)d k=%(k)d minejar_path=%(minejar_path)s"

def main(k=100000, jobname='dispatch_mine', n_nodes=13, n_ppn=4, walltime='6:00:00', work_dir='/fs/lustre/osu6683', minejar_path="/fs/lustre/osu6683/MINE.jar", tab_fname=None, bool_fname=None, dry=False):
  assert tab_fname, bool_fname

  # read bool matrix to get number of missing pairs
  B = np.load(bool_fname)
  n = np.size(np.where(B == 0)[0])
  # generate job script 
  script_path = os.path.dirname(os.path.realpath(__file__))
  offset = 0
  dispatch_script_fname= os.path.join(work_dir, "tmp_script_%s.sh" % os.path.basename(tab_fname))
  print "Creating batch script '%s'..." % dispatch_script_fname
  fp = open(dispatch_script_fname, 'w')
  while offset < n:
    cmd = CMD % {
      'script_path': script_path, 
      'tab_fname': tab_fname, 
      'bool_fname': bool_fname, 
      'offset': offset, 
      'k': k,
      'minejar_path': minejar_path
      }
    fp.write(cmd + '\n')
  fp.close()
  
  # submit job script
  qsub_script = TEMPLATE % {'jobname': jobname, 'n_nodes': n_nodes, 'n_ppn': n_ppn, 'walltime': walltime, 'dispatch_script': dispatch_script_fname}

  print qsub_script
  if not dry:
    p = subprocess.Popen("qsub", stdin=subprocess.PIPE)
    p.communicate(input=qsub_script)
    p.stdin.close()
    print "Batch job submitted."
  else:
    print "Dry run."

  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))

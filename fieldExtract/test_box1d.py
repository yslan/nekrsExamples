"""Validate the box1d sampler (box mode, boxDims {1,1,41} -> line along z)."""
import glob
import sys

import fe_read

pattern = "box1d[0-9]*.vts"
fname = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob(pattern))[-1]
sys.exit(0 if fe_read.check(fname) else 1)

"""Validate the box2d sampler (box mode, boxDims {11,1,41} -> XZ plane)."""
import glob
import sys

import fe_read

pattern = "box2d[0-9]*.vts"
fname = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob(pattern))[-1]
sys.exit(0 if fe_read.check(fname) else 1)

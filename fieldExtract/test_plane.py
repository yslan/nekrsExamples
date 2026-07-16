"""Validate the plane sampler (box mode, boxDims {11,11} -> XY plane, "Plane" tag)."""
import glob
import sys

import fe_read

pattern = "planePlane*.vts"
fname = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob(pattern))[-1]
sys.exit(0 if fe_read.check(fname) else 1)

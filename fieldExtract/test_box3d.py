"""Validate the box3d sampler (box mode, boxDims {11,7,21} -> full box, "Box" tag)."""
import glob
import sys

import fe_read

pattern = "box3dBox*.vts"
fname = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob(pattern))[-1]
sys.exit(0 if fe_read.check(fname) else 1)

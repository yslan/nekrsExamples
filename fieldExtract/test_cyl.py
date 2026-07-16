"""Validate the cyl sampler (points mode, r-theta-z grid {6,9,11} -> "Box" tag)."""
import glob
import sys

import fe_read

pattern = "cylBox*.vts"
fname = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob(pattern))[-1]
sys.exit(0 if fe_read.check(fname) else 1)

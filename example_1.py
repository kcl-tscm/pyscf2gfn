import sys
import numpy as np
from pyscf2gfn import pyscfxtb
from pyscf import gto, scf, geomopt

#Geometry optimization using GFn-XTB with minimum parameters
pyscfxtb.opt("PATHTOTHEEXECUTABLE/xtb", "ferrocene.xyz")

#Geometry optimization using GFn-XTB with different parameters
pyscfxtb.opt("PATHTOTHEEXECUTABLE/xtb", "ferrocene.xyz", threshold="tight")

#Hessian calculation using GFN-XTB with minimum ammount of information
pyscfxtb.hess("PATHTOTHEEXECUTABLE/xtb", "xtbopt.xyz", etemp=0)

#Optimization + Hessian calculation using GFN-XTB
pyscfxtb.opthess("PATHTOTHEEXECUTABLE/xtb", "xtbopt.xyz",version=1)



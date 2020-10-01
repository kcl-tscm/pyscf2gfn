import sys
import numpy as np
from pyscf2gfn import pyscfxtb
from pyscf import gto, scf, geomopt

########################################
# Example 2                            #
# Workflow combining GFn-XTB and PySCF #
########################################

# Creating a PySCF object 
mf = gto.Mole()

# Defining the properties of this object (optimized geometry from Gfn-XTB) using
# pyscf2gfn own geometry module.

mf.atom = pyscfxtb.xyz_to_pyscf('xtbopt.xyz')

# PySCF calculation properties   
mf.basis = {'H': 'sto3g', 'C': 'sto3g', 'Fe': 'sto3g'}

# Relaxing using PySCF pyberny with a different level of theory included in PySCF
relax = mf.apply(scf.RHF)
mol_eq = geomopt.optimize(mf)

# Printing the obtained structure using pyscf2gfn instead of the normal PySCF xyz module
pyscfxtb.pyscf_to_xyz(mol_eq,"final_mol.xyz")
 

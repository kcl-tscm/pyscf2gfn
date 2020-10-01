import sys
import numpy as np
from pyscf import gto, scf, geomopt
import os
import subprocess as sp
import shutil


class pyscfxtb:
         
    @staticmethod     
    def xyz_to_pyscf(filename):
        """
           Returns a list with the coordinates of a molecule in .xyz
           format to be directly read by pyscf.

           Args:
           -----------
           filename : : (str) a .xyz file from any molecular system
     
           Returns: 
           --------
           :class: 'list'
    	   Coordinates of the molecular in pyscf format
        """

        with open(filename, "r") as f:
         num_atoms=int(f.readline())
         blank=f.readline()
         coord= list([i.strip('\n') for i in f])
         coord = list(filter(None, coord))
         total = []
         for i in coord:
            total.append(i.split())
            final = []
            for i in range(len(total)):
              final.append((total[i][0],total[i][1:4]))
         return final


    @staticmethod
    def pyscf_to_xyz(mol_eq, outfile):
        """
          Returns a file with the coordinates of an optimized molecule in .xyz
          format.

          Args:
          -----------
          mol : : 'np.darray' an array containing the coordinates reported by PySCF 
                which is provided by using the keyword NAMEOFYOURSYSTEM.build()
     
          Returns: 
          --------
          :'file':: final_geometry.xyz
                Optimized coordinates of the molecule in .xyz format
        """

        mol_eq.build()
        bohr2angs=0.529177
        coords = bohr2angs*mol_eq.atom_coords()
        with open(outfile, "w") as myfile:
         myfile.write('{}'.format(int(len(coords))))
         myfile.write('\n')
         myfile.write('\n')
         for i in range(len(mol_eq.atom)):
            myfile.write('{} {:.8f} {:.8f} {:.8f}'.format(str(mol_eq._atom[i][0]),float(coords[i][0]),float(coords[i][1]),float(coords[i][2])))
            myfile.write("\n")
        myfile.close()

    @staticmethod
    def opt(xtb_path, xyz, num_cores=None, threshold=None, cycles=None, charge=None,
            num_unpaired_electrons=None, etemp=None, version=None):
        """
         Returns a file with the coordinates of an optimized molecule in .xyz
         format.

         Args:
         -----------
         xtb_path : : 'string' real path to the folder where the binary for xtb is stored
         xyz : : 'np.darray' an array containing the coordinates reported by PySCF 
                  which is provided by using the keyword NAMEOFYOURSYSTEM.build()
         num_cores : : 'np.int' Number of cores used in the xtb calculation
         threshold : : 'string' Level of accuracy used to optimise the system as defined in xtb manual
         cycles    : : 'np.int' Number of SCF iterations used in the xtb calculation
         charge    : : 'np.int'   Charge as defined in the xtb manual
         num_unpaired_electrons : : 'np.int' Number of unpaired electrons in the calculation as defined in xtb manual
         etemp : : 'np.int'  Electronic temperature as defined in xtb manual
         version : : 'np.int' version of Hamiltonian used for the xtb calculation
     
         Returns: 
         --------
         :'file':: final_geometry.xyz

         Optimized coordinates of the molecule in .xyz format

        """

        num_cores = int(1) if num_cores is None else num_cores
        threshold = str("tight") if threshold is None else threshold
        cycles = int(500) if cycles is None else cycles
        charge = int(0) if charge is None else charge
        num_unpaired_electrons = int(0) if num_unpaired_electrons is None else num_unpaired_electrons
        etemp = int(0) if etemp is None else etemp
        version = int(2) if version is None else version
        
        cmd = (
            f'{xtb_path} --etemp {etemp} {xyz} '            
            f'--parallel {num_cores} '
            f'--opt {threshold} '
            f'--cycles {cycles} '
            f'--chrg {charge} '
            f'--uhf {num_unpaired_electrons} '
            f'--gfn {version} '
          )
        # Note that sp.call will hold the program until completion
        # of the calculation.
        sp.call(
           cmd,
           stdin=sp.PIPE,
           stderr=sp.PIPE,
           # Shell is required to run complex arguments.
           shell=True
          )
      
    @staticmethod
    def hess(xtb_path, xyz, num_cores=None, charge=None, num_unpaired_electrons=None, etemp=None, version=None):
      """
      Returns a file with the coordinates of an optimized molecule in .xyz
      format.

      Args:
      -----------
      xtb_path : : 'string' real path to the folder where the binary for xtb is stored
      xyz : : 'np.darray' an array containing the coordinates reported by PySCF 
           which is provided by using the keyword NAMEOFYOURSYSTEM.build()
      num_cores : : 'np.int' Number of cores used in the xtb calculation
      charge    : : 'np.int'   Charge as defined in the xtb manual
      num_unpaired_electrons : : 'np.int' Number of unpaired electrons in the calculation as defined in xtb manual
      etemp : : 'np.int'  Electronic temperature as defined in xtb manual
      version : : 'np.int' version of Hamiltonian used for the xtb calculation
     
      Returns: 
      --------
      :'file':: final_geometry.xyz
        Optimized coordinates of the molecule in .xyz format
      """
      
      num_cores = int(1) if num_cores is None else num_cores
      charge = int(0) if charge is None else charge
      num_unpaired_electrons = int(0) if num_unpaired_electrons is None else num_unpaired_electrons
      etemp = int(0) if etemp is None else etemp
      version = int(2) if version is None else version
      
      cmd = (
            f'{xtb_path} --etemp {etemp} {xyz} '
            f'--parallel {num_cores} '
            f'--hess '
            f'--chrg {charge} '
            f'--uhf {num_unpaired_electrons} '
            f'--gfn {version} '
        )
        # Note that sp.call will hold the program until completion
        # of the calculation.
      sp.call(
           cmd,
           stdin=sp.PIPE,
           stderr=sp.PIPE,
           # Shell is required to run complex arguments.
           shell=True
        )

    @staticmethod
    def opthess(xtb_path, xyz, num_cores=None, charge=None, num_unpaired_electrons=None, etemp=None, version=None):
      """
         Returns a file with the coordinates of an optimized molecule in .xyz
         format.

         Args:
         -----------
         xtb_path : : 'string' real path to the folder where the binary for xtb is stored
         xyz : : 'np.darray' an array containing the coordinates reported by PySCF 
           which is provided by using the keyword NAMEOFYOURSYSTEM.build()
         num_cores : : 'np.int' Number of cores used in the xtb calculation
         charge    : : 'np.int'   Charge as defined in the xtb manual
         num_unpaired_electrons : : 'np.int' Number of unpaired electrons in the calculation as defined in xtb manual
         etemp : : 'np.int'  Electronic temperature as defined in xtb manual
         version : : 'np.int' version of Hamiltonian used for the xtb calculation
     
        Returns: 
        --------
        :'file':: final_geometry.xyz
        Optimized coordinates of the molecule in .xyz format

      """
      
      num_cores = int(1) if num_cores is None else num_cores
      charge = int(0) if charge is None else charge
      num_unpaired_electrons = int(0) if num_unpaired_electrons is None else num_unpaired_electrons
      etemp = int(0) if etemp is None else etemp
      version = int(2) if version is None else version

      cmd = (
            f'{xtb_path} --etemp {etemp} {xyz} '
            f'--parallel {num_cores} '
            f'--ohess '
            f'--chrg {charge} '
            f'--uhf {num_unpaired_electrons} '
            f'--gfn {version} ' 
        )
        # Note that sp.call will hold the program until completion
        # of the calculation.
      sp.call(
           cmd,
           stdin=sp.PIPE,
           stderr=sp.PIPE,
           # Shell is required to run complex arguments.
           shell=True
        )



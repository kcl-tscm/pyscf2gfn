# pyscf2gfn
Parser code written in Python for using GFn-XTB as a relaxation engine for providing relaxed structures to PySCF directly within the same script.
One example has been added to display the use of the library. The three most important functionalities, namely, gfnxtb_opt, gfnxtb_hess and gfnxtb_opthess
have been added to be used as a single python code line. 

The functions must be provided with a valid path to the corresponding xtb executable (which can be downloaded here: https://github.com/grimme-lab/xtb)
and valid options, such as a valid .xyz file, number of desired cores for OMP parallelelization schemes used in xtb, formal charge number, number of 
unpaired electrons and electronic temperature (difficult cases of convergence like metallic systems). 

More functionalities are going to be added in the near future. Contributions are welcome!

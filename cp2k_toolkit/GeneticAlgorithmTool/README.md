DESCRIPTION in breif:

This folder contains files for running Genetic Algorithm implemented on ASE interfaced with several calculators.

To run GA smoothly, 

for cp2k user, if you want to run on middle or large system containing more than 100 atoms (may be fewer) and if you are using intel-mpi, you should re-compile cp2k because you need to substitute cp2k_shell.F in source code src/start/ by the one in ./code_modification/.  
It is because atomic positions will send from python ASE side to cp2k side by PIPE. However for intel-mpi, larger amount of information sending via PIPE is discouraged, information will be truncated in a rude manner, furtherly this will cause dead lock of cp2k.  
Also there are some changes in ASE, you should find where ASE stores its py code, enter the folder calculator and save a new file cp2k_sendFile.py. Then whenever you want to use GA interfaced with cp2k, import cp2k_sendFile, instead of cp2k. cp2k_sendFile.py is also provided in folder ./code_modification/  
  
If you want to run GA locally in series, it is very easy and almost no change should be made to adjust your computer,  
                          ... in parallel, I write a main program ga_p_local.py, you can edit input script ga.inp, and run GA with command "python ga_p_local.py ga.inp"  
                  ... on queue system, I write a mian program ga_slurm.py for slurm user. You need to adjust to your queue system in several aspects, such as function jtg for creating job submitting script, parameter lists relevant with "squeue" and "sbatch", should be adjust as "queue" or something else and "qsub", for PBS users.  
                      

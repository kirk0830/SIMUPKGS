from ase.optimize import BFGS
from ase.io import read, write
from ase.calculators.cp2k_sendFile import CP2K
from ase.ga.relax_attaches import VariansBreak
import sys

import os

fname = sys.argv[1]

print('Now relaxing {0}'.format(fname))
a = read(fname)
project_name = fname[11:-5]

os.system('mkdir {}'.format(project_name))
os.chdir(project_name)

cutoff_Ry = 400
a.calc = CP2K(
label=project_name,
debug=True,
        basis_set = 'DZVP-MOLOPT-SR-GTH',
        charge = 0.,
        cutoff = 13.605662285137*cutoff_Ry,
        force_eval_method = 'quickstep',
        pseudo_potential = 'GTH-PBE',
        stress_tensor = True,
        uks = False,
        xc = 'PBE',
        max_scf = 50,
        print_level = 'MEDIUM',
        inp = '''
&FORCE_EVAL
    &DFT
	&QS
            METHOD xTB
            &xTB
                DO_EWALD .TRUE.
                CHECK_ATOMIC_CHARGES .FALSE.
                &PARAMETER
                    DISPERSION_PARAMETER_FILE dftd3.dat
                &END PARAMETER
            &END xTB
	&END QS
	&POISSON
            &EWALD
                EWALD_TYPE SPME
            &END EWALD
	&END POISSON
        &MGRID
            REL_CUTOFF 30
            NGRIDS 5
        &END MGRID
        &SCF
            SCF_GUESS ATOMIC
            EPS_SCF 1.0E-5
            &OUTER_SCF
                EPS_SCF 1.0E-5
                MAX_SCF 10
            &END OUTER_SCF
            &OT
                ALGORITHM IRAC
                ENERGY_GAP 0.5
                PRECONDITIONER FULL_SINGLE_INVERSE
                MINIMIZER DIIS
            &END OT
        &END SCF
    &END DFT
&END FORCE_EVAL
        '''
        )

dyn = BFGS(a, trajectory=None, logfile=None)
vb = VariansBreak(a, dyn)
dyn.attach(vb.write)
dyn.run(fmax=0.05)

os.chdir('..')
a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()

write(fname[:-5] + '_done.traj', a)

print('Done relaxing {0}'.format(fname))


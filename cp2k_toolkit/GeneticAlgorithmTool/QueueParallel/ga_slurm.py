from random import random

from ase.io import write

from ase.ga.data import DataConnection
from ase.ga.population import Population

from ase.ga.standard_comparators import InteratomicDistanceComparator


from ase.ga.utilities import get_all_atom_types
from ase.ga.utilities import closest_distances_generator

from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation

from ase.ga.pbs_queue_run import PBSQueueRun
from time import sleep
import os

std_redirect_file = 'ga.out'
calc_config='ga.calc_cp2k_xtb.py'
database_name = 'gadb.db'
tmp_folder = 'tmp_folder'

if os.path.isfile(std_redirect_file):
    os.remove(std_redirect_file)

f_std_redirect_f = open(std_redirect_file, 'a+', encoding = 'utf-8')

def jtg(job_name, traj_file):

    s = '#!/bin/bash\n'
    s += '#SBATCH -p amd_256\n'
    s += '#SBATCH -N 1\n'
    s += '#SBATCH -n 64\n'
    s += '#SBATCH -J {0}\n'.format(job_name)
    s += '#SBATCH -x v1805\n'

    s += 'source /public4/soft/modules/module.sh\n'

    s += 'module load *******/17.0.5-cjj-public4-public4\n' 
    s += 'module load libxsmm/1.15-icc17-lcc-public4\n'
    s += 'module load libint/2.6.0-cp2k-icc17-lcc-public4\n'
    s += 'module load libxc/4.3.4-icc17-ls-public4\n'
    s += 'module load fftw/3.3.8-mpi-public4\n'

    s += 'export PATH=/*******/*******/*******/miniconda3/bin:$PATH\n'
    s += 'export PATH=/*******/*******/*******/cp2k_8.1_modified/cp2k-8.1.0/exe/Linux-x86-64-intel-minimal:$PATH\n'

    s += 'export ASE_CP2K_COMMAND=\"mpirun -n 64 cp2k_shell.psmp\"\n'
    s += 'export CP2K_DATA_DIR=\"/*******/*******/*******/POTENTIAL_AND_BASIS\"\n'
    s += 'python {} {}\n'.format(calc_config, traj_file)
    return s

def redirect_print(strinfo):
    '''
    redirect print function to file
    '''
    f_std_redirect_f.writelines(strinfo+'\n')

population_size = 20
mutation_probability = 0.3
n_to_test = 20

da = DataConnection(database_name)

job_prefix = database_name[5:-3]
pbs_run = PBSQueueRun(da,
                      tmp_folder=tmp_folder,
                      job_prefix=job_prefix,
                      n_simul=10,
                      job_template_generator=jtg,
                      qsub_command='sbatch',
                      qstat_command='squeue')
time_to_wait = 180

atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)

blmin = closest_distances_generator(all_atom_types,
                                    ratio_of_covalent_radii=0.7)

comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                     pair_cor_cum_diff=0.015,
                                     pair_cor_max=0.7,
                                     dE=0.02,
                                     mic=False)

pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
mutations = OperationSelector([1., 1., 1.],
                              [MirrorMutation(blmin, n_to_optimize),
                               RattleMutation(blmin, n_to_optimize),
                               PermutationMutation(n_to_optimize)])
redirect_print('\n')
redirect_print('='*50)
redirect_print('PREPARATION| structure relaxation procedure')
while da.get_number_of_unrelaxed_candidates() > 0:
    
    redirect_print('-(sub-motion)-> There are still unrelaxed structures, check if can submit relaxation jobs now...')
    if not pbs_run.enough_jobs_running():
        redirect_print('-(sub-motion)-> Yes')
        redirect_print('-(sub-motion)-> Presently there are {} jobs running'.format(pbs_run.number_of_jobs_running()))
        a = da.get_an_unrelaxed_candidate()
        pbs_run.relax(a)
        redirect_print('-(sub-motion)-> New job submitted.')
        redirect_print('-(sub-motion)-> Now the number of unrelaxed structures: {}'.format(da.get_number_of_unrelaxed_candidates()))
    else:
        redirect_print('-(sub-motion)-> No')
        redirect_print('-(sub-motion)-> wait till it is available to submit new jobs... I will sleep for {}s'.format(time_to_wait))
        sleep(time_to_wait)

if pbs_run.number_of_jobs_running() != 0:

    redirect_print('PREPARATION| Wait! There are still jobs unfinished.')
while pbs_run.number_of_jobs_running() != 0:

    redirect_print('PREPARATION| Wait for {} seconds and check again...'.format(time_to_wait))
    sleep(time_to_wait)

redirect_print('PREPARATION| All starting structures are relaxed!')
redirect_print('='*50)

population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

redirect_print('MUTATION| Now I come to the second job submitting section, to relax newly generated structures!')

for i in range(n_to_test):

    redirect_print('-(sub-motion)-> Now it is No.{} the new individual newly generated. Process ({}/{})'.format(i+1, i+1, n_to_test))
    a1, a2 = population.get_two_candidates()
    a3, desc = pairing.get_new_individual([a1, a2])
    if a3 is None:
        continue
    da.add_unrelaxed_candidate(a3, description=desc)

    if random() < mutation_probability:
        a3_mut, desc = mutations.get_new_individual([a3])
        if a3_mut is not None:
            da.add_unrelaxed_step(a3_mut, desc)
            a3 = a3_mut
    while pbs_run.enough_jobs_running():

        redirect_print('-(sub-motion)-> There are enough number of jobs running... I will wait for {} seconds and then check if I can submit new job(s)'.format(time_to_wait))
        sleep(time_to_wait)

    pbs_run.relax(a3)
    redirect_print('-(sub-motion)-> New job submitted.')

if pbs_run.number_of_jobs_running() != 0:

    redirect_print('MUTATION| Wait! There are still jobs unfinished.')
while pbs_run.number_of_jobs_running() != 0:

    redirect_print('MUTATION| Wait for {} seconds and check again...'.format(time_to_wait))
    sleep(time_to_wait)

coords_to_write = da.get_all_relaxed_candidates()
write('all_candidates.traj', coords_to_write)
if len(coords_to_write) == 0:

    redirect_print('ERROR| Something still goes wrong, no structures are generated, exit.')
    redirect_print('ERROR| execuating command scancel -p amd_256...')
    os.system('scancel -p amd_256')
    exit()
redirect_print('='*50)
redirect_print('POST-PROCESSING| Now all individuals are saved in the same database {}'.format(database_name))
idx_coord = 0
for icoord in coords_to_write:
    idx_coord += 1
    #write('candidate-{}.xyz'.format(idx_coord), icoord)
    s='-'*50
    s+='\n'
    s+='    Report:\n'
    s+='    Candidate number:     {}\n'.format(idx_coord)
    s+='    Candidate Energy:     {}\n'.format(icoord.get_potential_energy())
    print(s)
redirect_print('-'*50)
redirect_print('POST-PROCESSING| Only print the most stable structure as *.xyz file')
redirect_print('='*50)
write('ga-most-stable.xyz', coords_to_write[0])

f_std_redirect_f.close()

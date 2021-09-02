from random import random
from ase.io import write
from time import sleep
from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.parallellocalrun import ParallelLocalRun
from os import remove
import sys

inp_fname = sys.argv[1]
inp_f_tag = open(file=inp_fname, mode='r', encoding='utf-8')
line='START'
kws_list = {}

while line:
    line=inp_f_tag.readline()
    if line.startswith('#'):
        continue
    else:
        line = line[:-1]
        words = line.split(' ')
        words = [word for word in words if word != '']
        if len(words) == 3:
            key = words[0]
            if key == 'write_xyz' or key == 'population_size' or key == 'n_to_test' or key == 'n_parallel' or key == 't_sleep_then_check':
                words[-1] = int(words[-1])
            elif key == 'tmp_save' or key == 'parameterization':
                if words[-1] == 'False':
                    words[-1] = False
                else:
                    words[-1] = True
            elif key == 'mutation_probability':
                words[-1] = float(words[-1])
            else:
                words[-1] = words[-1][1:-1]
            kws_list[key] = words[-1]
inp_f_tag.close()
verbosity = kws_list['verbosity']
std_redirect = kws_list['std_redirect']
prefix = kws_list['prefix']
outdir = kws_list['outdir']
tmp_save = kws_list['tmp_save']
write_xyz = kws_list['write_xyz']
database_filename = kws_list['database_filename']
population_size = kws_list['population_size']
mutation_probability = kws_list['mutation_probability']
n_to_test = kws_list['n_to_test']
parameterization = kws_list['parameterization']
n_parallel = kws_list['n_parallel']
t_sleep_then_check = kws_list['t_sleep_then_check']
external_calc = kws_list['external_calc']
verbo_flag = 0
if verbosity == 'SILENT':
    verbo_flag = 0
elif verbosity == 'LOW':
    verbo_flag = 1
elif verbosity == 'MEDIUM':
    verbo_flag = 2
elif verbosity == 'HIGH':
    verbo_flag = 3
elif verbosity == 'DEBUG':
    verbo_flag = 999
else:
    verbo_flag = 0

if write_xyz < 0:
    write_xyz = population_size + n_to_test
f_name_out = prefix[:6]+'.out'
f_out = open(f_name_out, 'a+', encoding = 'utf-8')
def ga_print(vbsti_thr, str_out, mode = std_redirect):
    if verbo_flag >= vbsti_thr:
        if mode == 'file':
            f_out.writelines(str_out+'\n')
        elif mode == 'stdout':
            print(str_out)
def countdown(t):
    
    while t:
        mins, secs = divmod(t, 60)
        timer = '{:02d}:{:02d}'.format(mins, secs)
        print(timer, end="\r")
        sleep(1)
        t -= 1

ga_print(3, 'MOTION| output file is created: {}'.format(f_name_out))

da = DataConnection(database_filename)
if not da:
    ga_print(3, 'MOTION| database {} successfully opened.'.format(database_filename))

parallel_local_run = ParallelLocalRun(data_connection=da,
                                      tmp_folder=outdir,
                                      n_simul=n_parallel,
                                      calc_script=external_calc)
if not parallel_local_run:
    ga_print(3, 'MOTION| local parallel job running manager successfully loaded.')

atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
n_to_optimize = len(atom_numbers_to_optimize)
slab = da.get_slab()
all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
ga_print(2,
'''MOTION| Input parameters summary:
        information amount:             {}
        print information to:           {}
        project name:                   {}
        temporary folder:               {}
        save all temporary files:       {}
        number of xyz files to write:   {}
        starting population database:   {}
        starting population size:       {}
        mutation probability:           {}
        generating population size:     {}
        number of jobs run in parallel: {}
        time to wait between 2 checks:  {}
        external calculator config:     {}
        -----------------------------------------
        Information collection before run:
        atom types shown with index:    {}
        atom types need to optimize:    {}
        numbers of atom to optimize:    {}
        cell information: 
        {}
        -----------------------------------------
        *NOTE: if you want to run your external calculator in parallel,
         remember to adjust your core numbers used by every jobs.
         For example, if you in total have 32 cores, you want to run 4
         jobs in parallel, then you should manually change configuration
         of your external calculator.
         For cp2k and qe, you should change your environment variable(s),
         $ASE_CP2K_COMMAND or $ASE_ESPRESSO_COMMAND in this way:
         mpirun -n 32 /path/to/cp2k_shell/cp2k_shell.psmp
         ->
         mpirun -n 8 /path/to/cp2k_shell/cp2k_shell.psmp
         .
         If you forget to adjust this, terrible thing will happen...

        *NOTE: This program cannot run on Windows system, because process
         control command it use is designed for Linux. If you want to run
         this on Windows, you should modify code of function:
         parallellocalrun.py -> __clean_up__(), change contents:
         p = Popen(['ps -x -U `whoami`']
         to counterpart of Windows.

        *Hint: if you find there is something wrong, stop this program with
         Ctrl + C

        -> Program will wait for {} seconds now...'''
.format(
    verbosity, std_redirect, prefix, outdir, tmp_save, write_xyz, database_filename, population_size,
    mutation_probability, n_to_test, n_parallel, t_sleep_then_check, external_calc,
    all_atom_types, atom_numbers_to_optimize, n_to_optimize,
    slab.get_cell(), t_sleep_then_check
))
countdown(t_sleep_then_check)

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

while da.get_number_of_unrelaxed_candidates() > 0:
    a = da.get_an_unrelaxed_candidate()
    parallel_local_run.relax(a)
    ga_print(2, 'PREPARATION| start a new optimization of individual in starting population.')
while parallel_local_run.get_number_of_jobs_running() > 0:
    ga_print(3, 'PREPARATION| There are still some optimization jobs running, wait for {} seconds to see if they will finish.'.format(t_sleep_then_check))
    sleep(t_sleep_then_check)
ga_print(2, 'PREPARATION| all candidates in starting population have been labelled as relaxed. Now start to generate new individuals.')
population = Population(data_connection=da,
                        population_size=population_size,
                        comparator=comp)

for i in range(n_to_test):
    ga_print(2, 'MUTATION| Now starting configuration number {0}'.format(i))
    a1, a2 = population.get_two_candidates()
    a3, desc = pairing.get_new_individual([a1, a2])
    if a3 is None:
        continue
    da.add_unrelaxed_candidate(a3, description=desc)
    ga_print(3, 'MUTATION| new candidate generated.')

    if random() < mutation_probability:
        a3_mut, desc = mutations.get_new_individual([a3])
        if a3_mut is not None:
            da.add_unrelaxed_step(a3_mut, desc)
            a3 = a3_mut
            ga_print(3, 'MUTATION| mutation occurs on present newly generated candidate.')
    # Relax the new candidate
    parallel_local_run.relax(a3)
    ga_print(2, 'MUTATION| configuration number {0} is labelled as relaxed.'.format(i))
    population.update()
    ga_print(3, 'MUTATION| population update.')

# Wait until the last candidates are relaxed
while parallel_local_run.get_number_of_jobs_running() > 0:
    ga_print(3, 'There are still some optimization jobs running, wait for {} seconds to see if they will finish'.format(t_sleep_then_check))
    sleep(t_sleep_then_check)

write('all_candidates.traj', da.get_all_relaxed_candidates())
ga_print(2, 'SUMMARY| all candidates are save in ASE-specific format trajectory file \'all_candidates.traj\'')
if write_xyz > 0:
    ga_print(3, 'SUMMARY| write xyz files, candidates are sorted by energies from the lowest to the highest.')
    candidates_list = da.get_all_relaxed_candidates()
    idx_candi = 1
    for icandi in candidates_list:
        write('candidate-{}.xyz'.format(idx_candi), icandi)
f_out.close()
if std_redirect == 'stdout':
    remove(f_name_out)

from doi101038nphys2133 import q

q_val = q(
    lmp_trjfile = 'PdO_1000K_0.trj',
    elements = ['O', 'Pd'],
    axis = 'x',
    rslist = [0.2],
    idx_t_max = 3,
    nx = 100,
    ny = 100,
    nz = 100,
    mode = 'save_memory'
)
print(q_val)

from gro_conne_gen import topol_gen as topol_gen

[list1, list2, list3] = topol_gen(gro_file = 'template.gro', param_file = 'gro2topol_conf.dat', d_tol = 0.2)
print('list1')
print(list1)
print('list2')
print(list2)
print('list3')
print(list3)
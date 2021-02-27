# QE turbo_davidson.x drho-eign file post-processing script
echo "Welcome>>>---QE turbo_davidson.x drho-of-eign file correction for ibrav = 0 ---<<<Welcome"
echo "Script instruction:"
echo "This script is to solve problem emerges in procedure scf -> tddfpt -> pp, when charge response is needed to be visualize."
echo "For more detailed information and problem description, see: https://www.mail-archive.com/users@lists.quantum-espresso.org/msg40158.html"

outdir=a101
num_eign=5

cell_a=10.3867
cell_b=11.5382
cell_c=29.5215

e_a=1.0
e_b=$(echo "scale=4;$cell_b/$cell_a"|bc)
e_c=$(echo "scale=4;$cell_c/$cell_a"|bc)

echo "------------------------------------------------------------------------"
echo "Post-processing infomation collection:"
echo "outdir = ./$outdir"
echo "number of drho-of-eign files: $num_eign"
echo "normalized cell parameters (currently only orthogonal box is supported):"
echo "e_a = $e_a"
echo "e_b = $e_b"
echo "e_c = $e_c"
echo "------------------------------------------------------------------------"

cd $outdir
ieign=1
while [ $ieign -le $num_eign ]
do
    drho=Pd-1-H2-tddfpt.drho-of-eign-$ieign
    echo "Read and change file $drho ..."
    sed -i "4i 0.0 0.0 $e_c" $drho
    sed -i "4i 0.0 $e_b 0.0" $drho
    sed -i "4i $e_a 0.0 0.0" $drho
    let "ieign++"
    echo "complete"
done

echo "Quit>>>---QE turbo_davidson.x drho-of-eign file correction for ibrav = 0 ---<<<Quit"
echo "This is a trivial toolkit that will be collected into SIMUPKGS, a collection of toolkits which will be of use when perform coordinate-related simulations, mainly works on xyz-formated file and LAMMPS dump file. Its own simulation module will be published in future."
echo "See more on GitHub: https://github.com/kirk0830/SIMUPKGS"

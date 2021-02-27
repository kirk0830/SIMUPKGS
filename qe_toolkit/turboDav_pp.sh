# QE turbo_davidson.x drho-eign file post-processing script
echo "Welcome>>>---QE turbo_davidson.x drho-of-eign file correction for ibrav = 0 ---<<<Welcome"
echo "Script instruction:"
echo "This script is to solve problem emerges in procedure scf -> tddfpt -> pp, when charge response is needed to be visualize."
echo "For more detailed information and problem description, see: https://www.mail-archive.com/users@lists.quantum-espresso.org/msg40158.html"

outdir=c2h3cho
num_eign=5

echo "Post-processing infomation collection:"
echo "outdir = ./$outdir"
echo "number of drho-of-eign files: $num_eign"

cd outdir
ieign=1
while [ $ieign -le $num_eign ]
do
    drho=drho-of-eign-$ieign
    echo "Read and change file $drho ..."
    sed -i '4i 0.0 0.0 1.0' $drho
    sed -i '4i 0.0 1.0 0.0' $drho
    sed -i '4i 1.0 0.0 0.0' $drho
    let "ieign++"
    echo "complete"
done

echo "Quit>>>---QE turbo_davidson.x drho-of-eign file correction for ibrav = 0 ---<<<Quit"
echo "This is a trivial toolkit that will be collected into SIMUPKGS, a collection of toolkits which will be of use when perform coordinate-related simulations, mainly works on xyz-formated file and LAMMPS dump file. Its own simulation module will be published in future."
echo "See more on GitHub: https://github.com/kirk0830/SIMUPKGS"

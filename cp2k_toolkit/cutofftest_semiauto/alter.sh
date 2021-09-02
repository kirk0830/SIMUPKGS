echo "Settings>>>----------CP2K parameter convergence test shell----------<<<Settings"
echo "[BATTLE CONTROL ONLINE] welcome to convergence test parameter alteration script"
# user control >>>
inpfile=cutofftest.inp
# mode options: pw, gau, vdw, kpt
# pw: CUTOFF parameter, cutoff of planewaves, varying from 100 to 700, by step 100 as default (not implemented yet)
# gau: REL_CUTOFF parameter, cutoff of Gaussian function, varying from 10 to 50, by step 10 as default (not implemented yet)
# vdw: R_CUTOFF parameter, cutoff of vdw correction, varying from 10 to 50, by step 10 as default (not implemented yet)
# kpt: SCHEME parameter, kpoints level, varying from 1 to 8, by step 1 as default (not implemented yet)
echo "[BATTLE CONTROL ONLINE] select parameter to be updated:"
echo "                        [1] CUTOFF [2] REL_CUTOFF [3] R_CUTOFF [4] KPOINTS"
read option
echo "[BATTLE CONTROL ONLINE] write update value of parameter..."
read param_inp
echo "[BATTLE CONTROL ONLINE] parameter_input: $param_inp"
echo "                        cp2k input script: $inpfile"
if [ $option = 1 ]
then
header=pw
sed -i '/^CUTOFF/d' $inpfile
sed -i "/&MGRID/a CUTOFF $param_inp" $inpfile
elif [ $option = 2 ]
then
header=gau
sed -i '/^REL_CUTOFF/d' $inpfile
sed -i "/&MGRID/a REL_CUTOFF $param_inp" $inpfile
elif [ $option = 3 ]
then
header=vdw
sed -i '/^R_CUTOFF/d' $inpfile
sed -i "/&PAIR_POTENTIAL/a R_CUTOFF [angstrom] $param_inp" $inpfile
elif [ $option = 4 ]
then
header=kpt
sed -i '/^SCHEME/d' $inpfile
sed -i "/&KPOINTS/a SCHEME MONKHORST-PACK $param_inp $param_inp $param_inp" $inpfile
fi
foldername=$header-$param_inp
mkdir $foldername
echo "[BATTLE CONTROL TERMINATED] A folder named $foldername has been created. Remember to save produced output files before start the next run."
echo "Return>>>----------CP2K parameter convergence test shell----------<<<Return"
echo "This is a trivial toolkit that will be collected into SIMUPKGS, a collection of toolkits which will be of use when perform coordinate-related simulations, mainly works on xyz-formated file and LAMMPS dump file. Its own simulation module will be published in future."
echo "See more on GitHub: https://github.com/kirk0830/SIMUPKGS"
echo "Quit>>>----------CP2K parameter convergence test shell----------<<<Quit"

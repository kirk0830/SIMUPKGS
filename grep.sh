echo "Settings>>>----------CP2K parameter convergence test shell----------<<<Settings"
echo "[BATTLE CONTROL ONLINE] welcome to convergence test results collection script"
# user control >>>
inpfile=cutofftest.inp
outfile=cutofftest.out
# program processing and string parsering part
projectLine=$(grep "PROJECT" $inpfile)
rawNames=${projectLine#*PROJECT }
projectName=${rawNames% *}
echo "[BATTLE CONTROL ONLINE] PROJECT_NAME line in CP2K input: $projectName"
echo "[BATTLE CONTROL ONLINE] select parameter whose results to be collected:"
echo "                        [1] CUTOFF [2] REL_CUTOFF [3] R_CUTOFF [4] KPOINTS"
read option
echo "[BATTLE CONTROL ONLINE] input series of parameters you have entered when use alter.sh, or input 'linspace' if you want use linspace to generate parameter list:"
read parameters

if [ $parameters = linspace ]
then
echo "[BATTLE CONTROL ONLINE] linspace function activated! 3 parameters required: start[included], end[included], stepsize:"
read lininp
linset=($lininp)
linStart=${linset[0]}
linEnd=${linset[1]}
linStep=${linset[2]}
echo "[BATTLE CONTROL ONLINE] Linspace will start from $linStart to $linEnd by stepsize $linStep"
numberWrite=$linStart
paraReiter=$numberWrite
i=1
while(( $numberWrite<$linEnd ))
do
numberWrite=$(expr $numberWrite + $linStep)
paraReiter="$paraReiter $numberWrite"
if [ $i = 100 ]
then
break
fi
let "i++"
done
parameters=$paraReiter
echo "[BATTLE CONTROL ONLINE] linspace interpolation check: $parameters"
fi

echo "[BATTLE CONTROL ONLINE] input filename where to store results..."
read filename

if [ $option = 1 ]
then
header=pw
elif [ $option = 2 ]
then
header=gau
elif [ $option = 3 ]
then
header=vdw
elif [ $option = 4 ]
then
header=kpt
fi

for ipara in $parameters
do
folderName=$header-$ipara
echo "enter folder named $folderName..."
cd $folderName
grepLine=$(grep "ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):" $outfile)
ener=${grepLine:50}
linePrint=$ipara$ener
echo "data read: $linePrint"
cd ..
echo $linePrint >> $filename
done
echo "[BATTLE CONTROL TERMINATED] Energy information has been recored in $filename."
echo "Return>>>----------CP2K parameter convergence test shell----------<<<Return"
echo "This is a trivial toolkit that will be collected into SIMUPKGS, a collection of toolkits which will be of use when perform coordinate-related simulations, mainly works on xyz-formated file and LAMMPS dump file. Its own simulation module will be published in future."
echo "See more on GitHub: https://github.com/kirk0830/SIMUPKGS"
echo "Quit>>>----------CP2K parameter convergence test shell----------<<<Quit"

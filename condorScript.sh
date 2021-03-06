#!/bin/zsh
workdir=`pwd`

echo JOBSUB::START starting job in directory $workdir

cd $MINICALODIR
source env.sh
cp build/exampleB4a $workdir/runGeant
cp batchrun.mac $workdir/
cd $workdir
./runGeant -m batchrun.mac -f $1_out
exitstatus=$?
if [ $exitstatus != 0 ]
then
     echo JOBSUB::FAIL Geant failed with status $exitstatus
     
     rm -f $1_out.root
     exit $exitstatus
fi


cp $1_out.root /eos/cms/store/cmst3/group/dehep/miniCalo2/prod/
exitstatus=$?
if [ $exitstatus != 0 ]
then
     echo JOBSUB::FAIL eos copy failed with status $exitstatus
else
     echo JOBSUB::SUCC job ended sucessfully
fi
rm -f $1_out.root
exit $exitstatus

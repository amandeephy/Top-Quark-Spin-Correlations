#PBS -l nodes=1:ppn=20
#PBS -l walltime=00:30:00
#PBS -N Parafly_submit_ttbar_1
#PBS -o /depot/cms/users/bakshi3/CMSSW_9_1_0_pre3/delphes/Ntuplization_jobs/Parafly_submit_ttbar_1.log
#PBS -e /depot/cms/users/bakshi3/CMSSW_9_1_0_pre3/delphes/Ntuplization_jobs/Parafly_submit_ttbar_1.err
#PBS -q debug
#!/bin/sh

#Neccesary for parallizing the jobs
module load utilities
module load parafly

#set up the CMSSW environment
echo `hostname`
echo `date`
source /cvmfs/cms.cern.ch/cmsset_default.sh
cp /depot/cms/users/bakshi3/CMSSW_9_1_0_pre3/CMSSW_new.tgz .
echo 'Un-tarring the CMSSW environment'
echo 'Currently in ::'
echo $PWD
tar -zxf CMSSW_new.tgz
rm CMSSW_new.tgz
echo 'Done with that'
cd CMSSW_9_1_0_pre3/src
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

cd $PWD/../delphes/
echo 'Begin processing files in Parafly_inputs ::'
echo $PWD
ls -lhart

ParaFly -c /depot/cms/users/bakshi3/CMSSW_9_1_0_pre3/delphes/Ntuplization_jobs/Parafly_inputs_ttbar_1  -CPU 20 -failed_cmds rerun.txt

echo 'Done with Parafly inputs'
cd ../
rm -rf CMSSW_9_1_0_pre3
echo `date`

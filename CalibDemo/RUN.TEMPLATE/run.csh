#!/bin/csh

#PBS -N matthew_full_model_run
#PBS -A NRAL0017
#PBS -l walltime=06:00:00
#PBS -q regular
#PBS -j oe
#PBS -m abe
#PBS -M arezoo@ucar.edu
#PBS -l select=1:ncpus=36:mpiprocs=36


### Run the executable
mpiexec_mpt ./wrf_hydro.exe #> stdErrorOut

echo $?

set writeout=1

if ($writeout == 1) then
## clean up since wrf hydro dosent allow spec of an output folder.
set outFolder=OUTPUT
mkdir $outFolder
mkdir $outFolder/00_RUN
cp namelist.hrldas $outFolder/00_RUN/.
cp hydro.namelist $outFolder/00_RUN/.
cp CHANPARM.TBL $outFolder/00_RUN/.
cp HYDRO.TBL $outFolder/00_RUN/.
cp URBPARM.TBL $outFolder/00_RUN/.
cp DISTR_HYDRO_CAL_PARMS.TBL $outFolder/00_RUN/.
cp LAKEPARM.TBL $outFolder/00_RUN/.
cp VEGPARM.TBL $outFolder/00_RUN/.
cp GENPARM.TBL $outFolder/00_RUN/.
cp MPTABLE.TBL $outFolder/00_RUN/.
cp GWBUCKPARM.TBL $outFolder/00_RUN/.
cp SOILPARM.TBL $outFolder/00_RUN/.
## LSM
mv RESTART* $outFolder/.
mv *LDASOUT_DOMAIN* $outFolder/.
mv *LSMOUT_DOMAIN* $outFolder/.
## HYDRO
mv HYDRO*DOMAIN* $outFolder/.
mv diag_hydro* $outFolder/.
mv frxst_pts_out.txt* $outFolder/.
mv qstrmvolrt_accum.txt $outFolder/.
mv *CHANOBS_DOMAIN* $outFolder/.
mv *CHRTOUT_DOMAIN* $outFolder/.
mv *CHRTOUT_GRID* $outFolder/.
mv *RTOUT_DOMAIN* $outFolder/.
mv GW_*.txt $outFolder/.
mv qlink*.txt $outFolder/.
mv *LAKEOUT* $outFolder/.
endif

exit 0

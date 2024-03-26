#!/bin/bash

#EXPERIMENT=NetTMS.01 #environmental variable
source /etc/biac_sge.sh
EXPERIMENT=`findexp NetTMS.01`
EXPERIMENT=${EXPERIMENT:?"Returned NULL Experiment"}

# to execute these scripts, ssh -X mas51@cluster.biac.duke.edu, qinteract, 
# need to be in the folder that has these scripts too. 
# /mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/
# qstat -j JOBNUM


for subj in 5025 5026; do echo $subj;
	cd /mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/
	for day in 1 2 3 4; do echo $day;
	#for day in 1; do echo $day;
		for run in 1 2 3; do echo $run;
		#for run in 2; do echo $run;
			for trial in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do echo $trial;
			#for trial in 1 2; do echo $trial;
				qsub -v EXPERIMENT=$EXPERIMENT /mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/qsub_ST.sh ${subj} ${day} ${run} ${trial}
				echo "Job submitted for Subj${subj} Day${day} Run${run} Trial${trial}---"
			done
		done
	done
done
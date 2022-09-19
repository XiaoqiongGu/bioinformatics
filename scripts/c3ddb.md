https://github.com/mghpcc-projects/c3ddb/wiki


# c3ddb manual
Create three script files [c3ddb.sh, transfer_c3ddb2mac.sh, transfer_mac2c3ddb.sh] in the tools/script folder

c3ddb.sh:
	
	ssh -i /Users/xiaoqiong/c3ddb-cluster/linux/c3ddb-key -l xiaoqiong c3ddb01.mit.edu  #key: 921020

transfer_c3ddb2mac.sh:

	#!/bin/bash
	scp -i ~/c3ddb-cluster/linux/c3ddb-key c3ddb01.mit.edu:$1 $2

transfer_mac2c3ddb.sh:

	#!/bin/bash
	scp -i ~/c3ddb-cluster/linux/c3ddb-key $1 c3ddb01.mit.edu:$2


    rsync -av --progress xiaoqiong@c3ddb-globus.mit.edu:/home/xiaoqiong/xiaoqiong-ost22-files.txt ./c3ddb_demo -e 'ssh -i /data1/scripts/c3ddb/linux/c3ddb-key'


## set up in mac bash_profile

set tools/script folder in the bash\_profile and type below codes in the ~/.bash_profile
	nano ~/.bash_profile
	PATH="/Users/xiaoqiong/tools:${PATH}"
	
after that type in `ctrl + x` to exit the nano mode and save the changes

refresh the bash_profile in the console

	source ~/.bash_profile
	
## how to use c3ddb

everytime to log onto the system is through the command line 

	c3ddb.sh
	
after that you are in the c3ddb.sh server ssession
type in `scratch` to enter into my data directory (not home directory as the space is small)

## how to submit new jobs in c3ddb

put all the command line in the bash script with your prefered name such as 16s.sh

	#!/bin/bash
	all the scirpts here
	
>note	please add #!/bin/bash as the first headline in the bash script as the task schedular system can recognise and process only with this first headline
	
after that type in `sbatch.sh 16s.sh` or `sbatch_defq.sh 16s.sh`

the output will be:

for `sbatch.sh 16s.sh` :
	
	sbatch -p sched_mem1TB -c 40 -t 5-00:00:00 --mem=1000000 -J 16s -o 16s.out -e 16s.err 16s.sh

or for `sbatch_defq.sh 16s.sh` :

	sbatch -p defq -c 40 -t 5-00:00:00 --mem=250000 -J 16s -o 16s.out -e 16s.err 16s.sh
	
after that, copy paste either of the command line in the command prompt and type in enter.

> note1: can choose either of the commands line, as they use different processors, usually what I do is to use the first command line. 
> 
> note2: please put the current environment into the conda environment list you preferred and then type in the commands so that the processors can recognize the current software settings.

## transfer files between c3ddb and mac
    tansfer_mac2c3ddb.sh absolute_path_files_in_mac file_directory_in_c3ddb
	tansfer_c3ddb2mac.sh absolute_path_files_in_c3ddb file_directory_in_mac


## anni's tips in using c3ddb
Please remember to change “YOURUSERNAME” to your user name on c3ddb and /path/to/private-key to the path to your key to c3ddb

My trick for convenient login and transferring jobs:

#Add the following to the ~/.bashrc or ~/.bash_profile (mac) on your local computer

    c3ddb_key=/path/to/private-key
    function c3ddb() {
        ssh -i $c3ddb_key -l YOURUSERNAME c3ddb01.mit.edu
    }

    function c3ddb2() {
        ssh -i $c3ddb_key -l YOURUSERNAME c3ddb-globus.mit.edu
    }

    function c3push(){
        files="$1"
        dir="$2"
        rsync -av --progress $files YOURUSERNAME@c3ddb-globus.mit.edu:$dir -e 'ssh -i /path/to/private-key'
    }

    function c3pull ()
    {
        files="$1";
        dir="$2";
        rsync -av --progress YOURUSERNAME@c3ddb-globus.mit.edu:$files $dir -e 'ssh -i /path/to/private-key'
    }

Login

    c3ddb2 # only for submitting jobs and transferring files
    c3ddb # for small jobs

Upload files

    c3push your_file directory_on_c3ddb

Download files

    c3pull your_file_on_c3ddb directory_on_local_computer

Directories

    /scratch/users/YOURUSERNAME # main working space
    /home/YOURUSERNAME # very limited space

My trick for checking and submitting job:

Add the following to the ~/.bashrc on c3ddb
You can change the settings of -c number of threads (default 40), -t time estimation of the job (default 2 days, max 5 days). You can check out their capacities here: https://github-wiki-see.page/m/abiwaters/Getting-Started-on-c3ddb/wiki/How-to-Submit-a-Batch-Job

    function jobmit() {
        bashdir="$1"
        jobname="$2"
        jobtype="$3"
        if [ "$jobtype" == big ]; then
            sbatch -p sched_mem1TB -c 40 -t 2-00:00:00 -J $jobname -e $bashdir.err -o $bashdir.out $bashdir
        else
            if [ "$jobtype" == small ]; then
                sbatch -p defq  -c 40 -t 2-00:00:00 -J $jobname  -e $bashdir.err -o $bashdir.out $bashdir
            else
                sbatch -p defq,sched_mem1TB,sched_mem4TB -c 40 -t 0-12:00:00 -x node325,node310 -J $jobname  -e $bashdir.err -o $bashdir.out $bashdir
            fi
        fi
    }

    function checkjob() {
        jobtype="$1"
        if [ "$jobtype" == big ]; then
            squeue -a -p sched_mem1TB
        else
            if [ "$jobtype" == small ]; then
                squeue -a -p defq
            else
                squeue -u YOURUSERNAME
            fi
        fi
    }

#Then source ~/.bashrc  

    source ~/.bashrc

1.	Submit job requiring big memory or urgent jobs
jobmit yourbash.sh yourjobname big

2.	Submit normal jobs 
jobmit yourbash.sh yourjobname small
Or:
jobmit yourbash.sh yourjobname

3.	Check your jobs
Checkjob

4.	Check all big jobs
Checkjob big

5.	Check all small jobs
Checkjob small

How to load pre-installed programs:

module avail # gives you the list of all installed programs
module add /path/program # link the program to your environment – you can now access the program.

How to check the quota

        srun lfs quota -u YOURUSERNAME /scratch

Others

    1.  launch an interactive session e.g. on one node with 16 cores
        and exclusive use of the node

        salloc -N 1 -n 16 -p defq --time=1:00:00 –exclusive
        salloc -N 1 -n 16 -p sched_mem1TB_centos7 --time=1:00:00 --exclusive


Tutorials 
http://www.tchpc.tcd.ie/node/74




#!/bin/bash
# #SBATCH -D /home/my_user_my_dir/
# Job name,  will be displayed on the showq command
#SBATCH -J MMPY-CUDA
# Filename for standard output
# At end of job, it is in the directory from which sbatch was invoked
#SBATCH -o MMPY-CUDA.o%j
#SBATCH -e MMPY-CUDA.e%j
#SBATCH --get-user-env
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1

#  The requested wall clock job time limit in HH:MM:SS
#  Your job will end when it exceeds this time limit
#SBATCH --time=00:10:00

#In case there are problems with a node echo the hostname and GPU id to the log file
hostname
echo $CUDA_VISIBLE_DEVICES

#COMMANDS GO HERE



# Print out the environment
printenv


date

mkdir -p peak

for n in 256 512 1024 2048
do
for bx in 4 8 16 32
do
for by in 4 8 16 32
do
r=$[(10+200000/($n/120+1)/($n/120+1)/($n/120+1))*$bx*$by/16/16]
if [[ $[$bx*$by] -le $[16*16] ]]
then
make clean >/dev/null 2>&1
make "bx=$bx" "by=$by" >/dev/null 2>&1
echo -n "n=$n, r=$r, bx=$bx, by=$by, " | tee -a "peak/summary.txt"
./mmpy -n $n -r $r -R | tee "peak/result_${n}_${bx}_${by}.txt" | grep 'Device computation time' | tee -a "peak/summary.txt"
fi
done
done
done

echo ">>> Job Ends"

date

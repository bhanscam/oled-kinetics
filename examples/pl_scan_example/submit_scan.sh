#!/bin/bash

# Define the number of jobs
total_jobs=8
cores=8
nodes=$((total_jobs / cores))
remainder=$((total_jobs % cores))

echo $nodes
echo $remainder

# Loop through each job and submit it
for node in $(seq 1 $nodes)
do
	# Submit job using sbatch
	sbatch -J scan_$node << EOF
#!/bin/bash
#SBATCH --job-name=job_$node
#SBATCH --time=01:00:00
#SBATCH -o sbatch_$node.out
#SBATCH -e sbatch_$node.err
#SBATCH -c 16
#SBATCH -N 1

for core in {1..$cores}
do
python input_pl_scan.py $node \$core $cores > run-$node-\$core.log &
done

wait

EOF
done

if [ $remainder != 0 ]
then
	final_node=$((nodes+1))
	sbatch -J scan_$final_node << EOF
#!/bin/bash
#SBATCH --job-name=job_$final_node
#SBATCH --time=01:00:00
#SBATCH -o sbatch_$final_node.out
#SBATCH -e sbatch_$final_node.err
#SBATCH -c 16
#SBATCH -N 1

for core in {1..$remainder}
do
python input_pl_scan.py $final_node \$core $cores > run-$final_node-\$core.log &
done

wait

EOF
fi


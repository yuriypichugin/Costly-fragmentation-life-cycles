function SubmitJob(){
cat <<EOS | sbatch 
#!/bin/sh

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-12:00:00

# Set name of job shown in squeue
#SBATCH --job-name TrueRandDelay

# Request CPU resources
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -o ./ConsoleOut/$i.output.txt
#SBATCH -e ./Errors/$i.error.txt

python ./Screen_risk_costs.py $i 

EOS
}

# Read the arguments
# my_path = $1

# Purge and load OpenMPI module
module purge
module load python/3.5.0
module list

for i in $(seq 0 100)
do
	# Run the simulation
	SubmitJob
done
	
# show jobs status
squeue -l | grep -i "RFL_risk"

# Ended
echo "Finished."
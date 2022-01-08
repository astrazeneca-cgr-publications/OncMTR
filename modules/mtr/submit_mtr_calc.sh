#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=1


helpFunction() 
{
   echo ""
   echo "Usage: $0 -d dataset -w win -s out_suffix -f ab_filtered"
   echo -e "\t-d dataset {gnomad} [Required]"
   echo -e "\t-w win: number of codons in each side of the centre of the window [Required]"
   echo -e "\t-s output dir suffix"
   echo -e "\t-f ab_filtered (boolean) - {yes, true, 1} or {no, false, 0}"
   exit 1 # Exit script after printing help
}

while getopts "d:w:s:f:" opt; do
   case "$opt" in
      d ) dataset="$OPTARG" ;;
      w ) win="$OPTARG" ;;
      s ) out_suffix="$OPTARG" ;;
      f ) ab_filtered="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# Print helpFunction in case required parameters are empty
if [ -z "$dataset" ] || [ -z "$win" ]; then 
   echo "Some or all of the parameters are empty";
   helpFunction
fi

if [ -z "$out_suffix" ]; then
	out_suffix="None"
fi
if [ -z "$ab_filtered" ]; then
	ab_filtered="no"
fi


echo "dataset: $dataset"
echo "win: $win"
echo "out_suffix: $out_suffix"
echo "ab_filtered: $ab_filtered"



python -u calculate_mtr_scores.py -d $dataset -w $win -s $out_suffix -f $ab_filtered



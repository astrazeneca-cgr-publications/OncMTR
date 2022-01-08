#!/bin/bash
##SBATCH -o out.annotate_oncmtr.gnomad_win31.all_genes.oncology
#SBATCH -o out.annotate_oncmtr.gnomad_win31.all_genes
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 48:00:0



# Run for specific gene sets
#: '
top_percentile=10	#default: 10
declare -a genesets=(
		     'All_prioritised_oncogenes'
		     'Shortlist_prioritised_oncogenes'
          	     'Manually_inspected_prioritised_oncogenes'
		     'AM_Leukemia-genes'
		     'CL_Leukemia-genes'
		     'Rasopathy-genes'
		     'Hematological-gene-examples'
		    )

declare -a variant_annotations=("clinvar" "oncology")


mkdir -p logs

cnt=0
for geneset in "${genesets[@]}"; do
	for var_annot in "${variant_annotations[@]}"; do

		echo "> ${geneset} - ${var_annot} ..."
		python -u annotate_oncmtr_signals.py -d gnomad -w 31 -a $var_annot -g "$geneset" -t ${top_percentile} > "logs/${geneset}-${var_annot}.out" & 

		cnt=$((cnt+1))
		if [ "$cnt" == 4 ]; then
			cnt=0
			wait
		fi

	done
done

#wait
#exit
#'


# Run for All-genes (only with sbatch)
#python -u annotate_oncmtr_signals.py -d gnomad -w 31 -a clinvar #&
python -u annotate_oncmtr_signals.py -d gnomad -w 31 -a oncology #&
#wait

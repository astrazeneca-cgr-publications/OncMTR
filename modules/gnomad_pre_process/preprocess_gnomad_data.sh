#!/bin/bash
#SBATCH -t 1:0:0
#SBATCH -o out.preprocess_gnomad_data


gnomad_dir="../../data/gnomad_GRCh38/"
ensembl_dir="../../data/ensembl/"

#<< 'TEMP-MULTILINE-COMMENT'
echo -e "Keep PASS variants only from gnomAD (already liftovered to GRCh38)\n"
cat $gnomad_dir/gnomad.exomes.r2.0.2.sites.processed.snp.uniq.vcf | awk '{if($7 == "PASS" || $1 ~ /^#/) print $0}' > $gnomad_dir/gnomad.exomes.r2.0.2.sites.processed.snp.uniq.pass.vcf


echo -e "Subset gnomAD PASS variants that are part of exons only (as defined by Ensembl GRCh38.86 gtf annotation)\n"
intersectBed -a $gnomad_dir/gnomad.exomes.r2.0.2.sites.processed.snp.uniq.pass.vcf -b $ensembl_dir/hsa_exons.GRCh38.86.merged.bed -header > $gnomad_dir/gnomad_exons_mutations.pass.vcf


echo -e "Prepare snpEff VCF input with gnomad exon PASS variants\n"
cat $gnomad_dir/gnomad_exons_mutations.pass.vcf | awk '{print $1"\t"$2"\t.\t"$4"\t"$5"\t.\t.\t."}' | sort -k1,1 -k2,2n  > $gnomad_dir/snpeff_input.gnomad_exons_mutations.pass.vcf

echo -e "\n---------------------------\n"
#TEMP-MULTILINE-COMMENT


# ==================== AB_MEDIAN FILTERED ==========================
# 
# > Keep gnomAD entries with AB_MEDIAN > 0.30 (PASS variants only):
# 
# -- Input: gnomad.exomes.r2.0.2.sites.processed.snp.uniq.pass.vcf
# -- Output: gnomad.exomes.r2.0.2.sites.processed.snp.uniq.pass.ab_filtered.vcf
# __________________________________________________________________
echo -e "Select gnomAD variants with AB_MEDIAN > 0.30\n"
python filter_by_ab_median.py $gnomad_dir/gnomad.exomes.r2.0.2.sites.processed.snp.uniq.pass.vcf
	
echo -e "Subset gnomAD PASS variants that are part of exons only (& AB_MEDIAN>0.30)\n"
intersectBed -a $gnomad_dir/gnomad.exomes.r2.0.2.sites.processed.snp.uniq.pass.ab_filtered.vcf -b $ensembl_dir/hsa_exons.GRCh38.86.merged.bed -header > $gnomad_dir/gnomad_exons_mutations.pass.ab_filtered.vcf

echo -e "Prepare snpEff VCF input with gnomad exon PASS variants (& AB_MEDIAN>0.30)\n"
cat $gnomad_dir/gnomad_exons_mutations.pass.ab_filtered.vcf | awk '{print $1"\t"$2"\t.\t"$4"\t"$5"\t.\t.\t."}' | sort -k1,1 -k2,2n  > $gnomad_dir/snpeff_input.gnomad_exons_mutations.pass.ab_filtered.vcf


echo -e 'Done."

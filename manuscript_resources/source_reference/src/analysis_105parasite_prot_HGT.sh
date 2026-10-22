#!/bin/bash

Align_path="PROJECT_ROOT/align/protein_align"
Analysis_path="PROJECT_ROOT/analysis"
name2genome_file="PROJECT_ROOT/105_name2accfil_withTaxid.txt"


cd $Align_path

#cut -f2  all_accession2taxid.txt \
#              | sort \
#              | uniq \
#              | taxonkit lineage \
#              | taxonkit reformat \
#              > all_accession2taxid.lineage.txt

#while read line 
#do
#	name=`echo "$line" | awk -F ";" '{print $NF}'`
#	taxid=`echo "$line" | awk '{print $1}'`
#	genome=`grep "$name" ~/HGT/parasite/105_parasite/105_name2accfil_withTaxid.txt | cut -f2`
#	if [[ -n "$genome" ]]; then
#		echo $taxid"	"$genome"	"$name
#	fi
#done < all_accession2taxid.lineage.txt > Taxid2genome.txt

## combine the accession and genome, Then divide the alignments into different genomes.
#python combine_accWithgenome.py Taxid2genome.txt all_accession2taxid.txt | sort | uniq > all_accession_genome_prot.txt

Acc2genome="PROJECT_ROOT/align/protein_align/all_accession_genome_prot.txt"
Distribute_python="PROJECT_ROOT/src/ditribute_blast.py"
mkdir -p genome_divide
cd genome_divide

#python $Distribute_python $Acc2genome ../mmseqs_157parasite_results.out
#python $Distribute_python $Acc2genome ../157_parasite.taxonomy.blast.iden30.out

## extract all target taxonomy info
#cut -f4  ../157_parasite.taxonomy.blast.iden30.out \
#               | sort \
#               | uniq \
#               | taxonkit lineage \
#               | taxonkit reformat \
#               > target_nt_taxid_lineage_tmp.txt

#cat target_nt_taxid_lineage_tmp.txt ~/HGT/parasite/105_parasite/align/nucl_align/HGT_retest/target_nt_taxid_lineage.txt | sort | uniq  > target_nt_taxid_lineage.txt

#Protein HGT verify
HGT_VERIFY_SCRIPT="PROJECT_ROOT/src/HGT_verify_parasite.prot.py"
HGT_FILTER_ENTROPY="PROJECT_ROOT/src/HGT_filter_entropy.prot.py"
for tmp in `awk -F "\t" '{print $2"ABC"$3}' $name2genome_file`
do
	genome=`echo -e $tmp | awk -F "ABC" '{print $1}'`
	taxid=`echo -e $tmp | awk -F "ABC" '{print $2}'`
	echo $genome"       "$taxid 
#	python $HGT_VERIFY_SCRIPT target_nt_taxid_lineage.txt $genome.acc2taxid.out "$taxid" > $genome.HGT_verify_new.out
#	python $HGT_FILTER_ENTROPY $genome.HGT_verify_new.out -o $genome.HGT_verify_filtered.entropy.out
	
#	sed '1d' $genome.HGT_verify_filtered.entropy.out | awk '{print $1}' > tmp.acc
#	seqtk subseq EXTERNAL_DATA_ROOT/HGT/parasite/new_157_parasite/protein/sequence.fasta tmp.acc > $genome.HGT_verify_filtered.entropy.fa 
#	cd-hit -c 0.9 -i $genome.HGT_verify_filtered.entropy.fa -o $genome.HGT_verify_filtered.entropy.clustered.fa

done
# clu_the results with per50 and per90.
#cat *.entropy.clustered.fa > 105_parasite_prot_HGT.entropy.fa

mmseqs easy-linclust 105_parasite_prot_HGT.entropy.fa 105_parasite_prot_HGT.entropy_clu50 tmp  --cov-mode 0 -c 0.8 --min-seq-id 0.5 --threads 32





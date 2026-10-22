#!/bin/bash
Align_path="PROJECT_ROOT/align/nucl_align"
Analysis_path="PROJECT_ROOT/analysis"
nucl2taxid_cpp="PROJECT_ROOT/src/nucl2prot_acc2taxid"
name2genome_file="PROJECT_ROOT/105_name2accfil_withTaxid.txt"
Analysis_path="PROJECT_ROOT/analysis"
Genome_path="EXTERNAL_DATA_ROOT/Data/Genome/parasite/188_parasites"
##classify the align results into its belonging genome.
Acc2genome="PROJECT_ROOT/all_accession_genome.txt"
Distribute_python="PROJECT_ROOT/src/ditribute_blast.py"
#for file in `ls $Align_path/*.out`
#do
#	cd $Align_path
#	python $Distribute_python $Acc2genome $file
#done

#HGT analysis
mkdir -p $Analysis_path
mkdir -p $Analysis_path/nucl
cd $Align_path/genome_divided

#for genome in `awk -F "\t" '{print $2}' $name2genome_file`
#do
#	$nucl2taxid_cpp 12305 $genome.out > $genome.acc2taxid.out	
#done


## extract all target taxonomy info
#cut -f4  *.acc2taxid.out \
#               | sort \
#               | uniq \
#               | taxonkit lineage \
#               > target_nt_taxid_lineage_withmetazoa.txt



#cut -f1 target_nt_taxid_lineage_withmetazoa.txt \
#               | sort \
#               | uniq \
#               | taxonkit lineage \
#               | taxonkit reformat \
#               > target_nt_taxid_lineage.txt

#HGT_verify
#HGT_VERIFY_SCRIPT="PROJECT_ROOT/src/HGT_verify_new.py"
HGT_VERIFY_SCRIPT="PROJECT_ROOT/src/HGT_verify_parasite.py"
HGT_FILTER_ENTROPY="PROJECT_ROOT/src/HGT_filter_entropy.prot.py"
MERGE_SCRIPT="PROJECT_ROOT/src/merge_bed.py"
EXTRACT_SCRIPT="PROJECT_ROOT/src/extract_HGT_nucl_position.py"
for tmp in `awk -F "\t" '{print $2"ABC"$3}' $name2genome_file`
do
    genome=`echo -e $tmp | awk -F "ABC" '{print $1}'`
    taxid=`echo -e $tmp | awk -F "ABC" '{print $2}'`
#    echo $genome"	"$taxid 
#    python $HGT_VERIFY_SCRIPT target_nt_taxid_lineage.txt $genome.acc2taxid.out "$taxid" > $genome.HGT_verify.out
#    python $HGT_FILTER_ENTROPY $genome.HGT_verify.out -o $genome.HGT_verify_fitered.out
#    sed -i '1d' $genome.HGT_verify_fitered.out
#    python $EXTRACT_SCRIPT $genome.out $genome.HGT_verify_fitered.out > $genome.HGT_verify_fitered.bed
#    python $MERGE_SCRIPT $genome.HGT_verify_fitered.bed  $genome.HGT_verify_fitered_merged.bed
#    seqtk subseq $Genome_path/$genome $genome.HGT_verify_fitered_merged.bed > $genome.HGT_verify_fitered_merged.fa
#    cp $genome.HGT_verify_fitered_merged.fa $Analysis_path/nucl
#    cp $genome.HGT_verify_fitered.out $Analysis_path/nucl
#    cp $genome.HGT_verify_fitered.bed $Analysis_path/nucl
#    cp $genome.HGT_verify_fitered_merged.bed $Analysis_path/nucl


##    awk -F "\t" '{if($2>0.9 && $5 > 1 && ($4 < 5 || $9 >= 0)) print $0}' $genome.HGT_verify.out > $genome.HGT_verify.out_pct80.out
##    awk -F "\t" '{if($2 > 0.8) print $0}' $genome.HGT_verify.out > $genome.HGT_verify.out_pct80.out
##    sed -i '1d' $genome.HGT_verify.out_pct80.out

##    python $EXTRACT_SCRIPT $genome.out $genome.HGT_verify.out_pct80.out > $genome.HGT_verify.out_pct80.bed
##    mv $genome.HGT_verify.out_pct80_merged.bed $genome.HGT_verify.out_pct80.bed
##    python $MERGE_SCRIPT $genome.HGT_verify.out_pct80.bed $genome.HGT_verify.out_pct80_merged.bed
##    seqtk subseq $Genome_path/$genome $genome.HGT_verify.out_pct80_merged.bed > $genome.HGT_verify.out_pct80_merged.fa
#   # get the host  
##    awk -F "\t" '{if($10 >= $11 && $9 > 0) print $0}' $genome.HGT_verify.out_pct80.out > $genome.HGT_verify.out_pct80.host_origin.out
##    python $EXTRACT_SCRIPT $genome.out $genome.HGT_verify.out_pct80.host_origin.out > $genome.HGT_verify.out_pct80.host_origin.bed
##    python $MERGE_SCRIPT $genome.HGT_verify.out_pct80.host_origin.bed $genome.HGT_verify.out_pct80.host_origin_merged.bed
##    seqtk subseq $Genome_path/$genome $genome.HGT_verify.out_pct80.host_origin_merged.bed > $genome.HGT_verify.out_pct80.host_origin_merged.fa
done

#Using HPC to calculate the first filtered HGTs.

#processing the second alignment results.
diamond_second_align_file="PROJECT_ROOT/align/105_parasite.first_HGT_verify.out80_ai0.core_nt"


#$nucl2taxid_cpp 12305 $diamond_second_align_file.out >  $diamond_second_align_file.acc2taxid.out


cd $Align_path/..
mkdir -p second_align

cd second_align
##divide the second align result and acc2taxid file into 105 parasites.
#python $Distribute_python $Acc2genome $diamond_second_align_file.out
#mkdir -p acc2taxid
#cd acc2taxid
#python $Distribute_python $Acc2genome $diamond_second_align_file.acc2taxid.out

#extract the taxonomy information
#cut -f4  $diamond_second_align_file.acc2taxid.out \
#               | sort \
#               | uniq \
#               | taxonkit lineage \
#	       | taxonkit reformat \
#               > ../target_nt_taxid_lineage.txt

for tmp in `awk -F "\t" '{print $2"ABC"$3}' $name2genome_file`
do
    genome=`echo -e $tmp | awk -F "ABC" '{print $1}'`
    taxid=`echo -e $tmp | awk -F "ABC" '{print $2}'`
    echo $genome"	"$taxid 
#    mv acc2taxid/$genome.out $genome.acc2taxid.out
#    python $HGT_VERIFY_SCRIPT ../target_nt_taxid_lineage.txt $genome.acc2taxid.out "$taxid" > $genome.HGT_verify_new.out
#    python $HGT_FILTER_ENTROPY $genome.HGT_verify_new.out -o $genome.HGT_verify_filtered.out
#    sed '1d' $genome.HGT_verify_filtered.out | cut -f1 | awk -F ":" '{print $1"\t"$2}' | awk -F "-" '{print $1"\t"$2}' > $genome.HGT_verify_filtered.bed
    
#    seqtk subseq $Genome_path/$genome $genome.HGT_verify_filtered.bed > $genome.HGT_verify_filtered.fa

done
cat *HGT_verify_filtered.fa > all_genome_HGT.nucl.fa


if false; then
## Analysis the result of parasite nucl HGT

cd $Analysis_path/nucl
mkdir -p From_host

for tmp in `awk -F "\t" '{print $2"ABC"$3}' $name2genome_file`
do
    genome=`echo -e $tmp | awk -F "ABC" '{print $1}'`
    taxid=`echo -e $tmp | awk -F "ABC" '{print $2}'`
    echo $genome"   "$taxid

#    cd-hit-est -i $genome.HGT_verify_fitered_merged.fa -c 0.9 -o $genome.HGT_verify_fitered_merged.clustered.fa
    # get the host origin HGT seq
    awk -F "\t" '{if($10 >= $11 && $9 > 0) print $0}' $genome.HGT_verify_fitered.out > From_host/$genome.HGT_verify_fitered.from_host.out 

    python $EXTRACT_SCRIPT $Align_path/genome_divided/$genome.out From_host/$genome.HGT_verify_fitered.from_host.out > From_host/$genome.HGT_verify_fitered.from_host.bed
    python $MERGE_SCRIPT From_host/$genome.HGT_verify_fitered.from_host.bed  From_host/$genome.HGT_verify_fitered.from_host.merged.bed
    seqtk subseq $Genome_path/$genome From_host/$genome.HGT_verify_fitered.from_host.merged.bed > From_host/$genome.HGT_verify_fitered.from_host.merged.fa
    num=`grep ">" From_host/$genome.HGT_verify_fitered.from_host.merged.fa | wc -l`
    echo "共找到潜在来自宿主纲的HGT序列共"$num"条。"
done 
fi

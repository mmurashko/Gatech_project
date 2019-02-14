#!/bin/bash
#directory="/Users/Matvey/Desktop/GeorgiaTech_project/data/6"

#out_fold="/Users/Matvey/Desktop/GeorgiaTech_project/output_data/6hrs"

#path_to_trimmomatic="/Users/Matvey/Desktop/GeorgiaTech_project/Trimmomatic-0.38"


directory=$1
out_fold=$2
path_to_trimmomatic=$3
path_to_adapters=$4
path_to_transcripts=$5
path_to_transcript_index=$6


mkdir $out_fold

echo "Running FastQC for all raw data"

for dir in $(ls $directory)
do
	echo $dir
	mkdir $out_fold/$dir
	mkdir $out_fold/$dir/Raw_data
	mkdir $out_fold/$dir/raw_data/FastQC_raw_data

	fastqc $directory/$dir/* -o $out_fold/$dir/raw_data/FastQC_raw_data
done 

echo "Merging all files"

for dir in $(ls $directory)
do
	mkdir $out_fold/$dir/Merged_data
	
	cat $directory/$dir/*R1* > $out_fold/$dir/Merged_data/R1.fastq.gz
	cat $directory/$dir/*R2* > $out_fold/$dir/Merged_data/R2.fastq.gz
done 

echo "Running FastQC for all merged data"

for dir in $(ls $directory)
do
	mkdir $out_fold/$dir/Merged_data/FastQC_merged_data
	
	fastqc $out_fold/$dir/Merged_data/R1.fastq.gz -o $out_fold/$dir/Merged_data/FastQC_merged_data
	fastqc $out_fold/$dir/Merged_data/R2.fastq.gz -o $out_fold/$dir/Merged_data/FastQC_merged_data
done 

echo "Trimming all files"

path_to_trimmomatic="/Users/Matvey/Desktop/GeorgiaTech_project/Trimmomatic-0.38"

for dir in $(ls $directory)
do
	mkdir $out_fold/$dir/Trimmed_data

	java -jar $path_to_trimmomatic/trimmomatic-0.38.jar PE $out_fold/$dir/Merged_data/R1.fastq.gz $out_fold/$dir/Merged_data/R2.fastq.gz\
	$out_fold/$dir/Trimmed_data/output_forward_paired.fq.gz $out_fold/$dir/Trimmed_data/output_forward_unpaired.fq.gz\
	$out_fold/$dir/Trimmed_data/output_reverse_paired.fq.gz $out_fold/$dir/Trimmed_data/output_reverse_unpaired.fq.gz\
	ILLUMINACLIP:$path_to_adapters/adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36

done

echo "Running FastQC for Trimmed_data"

for dir in $(ls $directory)
do
	mkdir $out_fold/$dir/Trimmed_data/FastQC_trimmed_data
	
	fastqc $out_fold/$dir/Trimmed_data/output_forward_paired.fq.gz -o $out_fold/$dir/Trimmed_data/FastQC_trimmed_data
	fastqc $out_fold/$dir/Trimmed_data/output_reverse_paired.fq.gz -o $out_fold/$dir/Trimmed_data/FastQC_trimmed_data
done

echo $path_to_transcripts
echo "Running Salmon for Trimmed_data"

salmon index -t $path_to_transcripts -i transcripts_index --type quasi -k 31

for dir in $(ls $directory)
do
	mkdir $out_fold/$dir/Salmon 

	salmon quant -i $path_to_transcript_index/transcripts_index\
	 -l IU -1 $out_fold/$dir/Trimmed_data/output_forward_paired.fq.gz -2 $out_fold/$dir/Trimmed_data/output_reverse_paired.fq.gz --gcBias -o $out_fold/$dir/Salmon 
done

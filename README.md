# RNA-Seq pipeline

This pipeline was created to visualize and analyze bulk RNA-Seq data. It uses recently created up-to-date alignment-free tools. 

## Getting Started

These instructions will help you to install all necessary software and run our pipeline.

### Prerequisites
This pipeline requires installation of additional software. In this section there is an explanation how to download everything. All pathways to executables must be added to global PATH variable.

#### FastQC

You can download it from the website: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>. We used version v0.11.7.
You should modify .bash_profile (for Mac) or .bashrc (for Linux) from your home directory

```
#Example for Mac:
cd ~
vim .bash_profile
#add the string to the end of bash_profile
export PATH="/Applications/FastQC.app/Contents/MacOS:$PATH"
#save changes and quit(:wq)
source ./bash_profile
```

#### Trimmomatic

You can download Trimmomatic from the website (binary file): <http://www.usadellab.org/cms/?page=trimmomatic>. We used version Trimmomatic-0.38

#### Salmon

Instructions on how to install Salmon(Miniconda is needed) can be found on their website: <https://combine-lab.github.io/salmon/getting_started/>. We used version v0.10.2.
To run salmon you also need to have a file with transcript anotation and to give its' directory as a command line argument.


## Preparing data
Run a bash pipeline. The following command line arguments are needed for this: 
* Directory with raw data
* Output directory
* Directory of Trimmomatic folder
* Directory of a adapters.fa file
* File with transcripts
* Directory for transcripts index 

In the directory with raw data folders with .fastq or .fastq.gz files should be stored
```
#Here is an example how to run it on Mac:
./pipeline_for_preparation.sh /Users/Matvey/Desktop/GeorgiaTech_project/data/6hrs /Users/Matvey/Desktop/GeorgiaTech_project/output_data/6hrs /Users/Matvey/Desktop/GeorgiaTech_project/Trimmomatic-0.38 /Users/Matvey/Desktop/GeorgiaTech_project /Users/Matvey/Desktop/GeorgiaTech_project/Homo_sapiens.GRCh38.cdna.all.fa /Users/Matvey/Desktop/GeorgiaTech_project
```

How it works:
1. FastQC  for raw data
2. Merging raw data files
3. FastQC for merged data
4. Trimming
* Quality trimming
* Length trimming         
* Adapter trimming 
5. FastQC for trimmed data
6. Salmon estimates transcript-level abundances using quasi-mapping



## Running the main pipeline
This pipeline contains two modes: "pathway" and "set".

How it works:
1. Estimates abundance (TPM) for every gene using tximport to match transcripts with genes
2. Calculates Log2FoldChange for CCM1/CCM2/CCM3
2.1. Finds and saves genes that were not expressed in NT or CCM1/CCM2/CCM3 
3. Plots LFC distribution for all expressed genes 
4. Extraction of genes associated with a specified pathway, visualization of LFC for differentially expressed genes 
5. Visualization of expression levels of arbitrary set of genes  

### "Pathway" mode
You can run a pipeline from the terminal. Start it using Rscript RNA_Seq_pipeline.R. As a command line arguments you need to specify the following:
1. mode(pathway/set)
2. output directory of preparation pipeline
3. .csv file with names of all samples
4. label added to every new file.
5. threshold
6. GO id


```
#Here is an example how to run it on Mac:
Rscript user_friendly_pipeline.R pathway /Users/Matvey/Desktop/GeorgiaTech_project/output_data/6hrs samples_6hrs.csv _3D_6hrs 1 GO:0007219
```


### "Set" mode
You can run a pipeline from the terminal. Start it using Rscript RNA_Seq_pipeline.R. As a command line arguments you need to specify the following:
1. mode(pathway/set)
2. output directory of preparation pipeline
3. .csv file with names of all samples
4. label added to every new file. 
5. threshold
6. directory of a text file with all wanted genes is needed

```
#Here is an example how to run it on Mac:
Rscript user_friendly_pipeline.R set /Users/Matvey/Desktop/GeorgiaTech_project/output_data/6hrs samples_6hrs.csv _3D_6hrs 1 /Users/Matvey/Desktop/GeorgiaTech_project/Set_of_wanted_genes.txt
```

## Output of pipeline
Files:
* abundabce*.csv - information about all genes, their abundance and LogFoldChange
* differentially_expressed_genes_label.csv - genes with |LFC| > threshold
* differentially_expressed_genes_CCM*_label.csv - genes with |LFC_CCM*| > threshold
* differential_expression_in_pathway_CCM*_label.csv - genes from pathway with |LFC_CCM*| > threshold
* enrichment_label.csv - enrichment of Notch pathway (p-values)
* graphs*.png - images with barplots for differentially expressed genes
* no_expression_label.csv - genes that is not expressed in NT and expressed in CCM1/CCM2/CCM3 or vice versa
* own_graph_*.png - distibution of LFC for all genes
* pathway_graph_CCM*_*.png - distibution of LFC for genes from pathway
* results_of_differential_expression_in_pathway_*.csv - all genes from pathway that are differentially expressed at least in one of the cases

# EvaHLA
Nextflow pipeline for EvaHLA

## System requirements ##
The pipeline requires >= 256 GB memory. It is recommended to have at least 8 cores. We typically run the pipeline on a server with dual Epyc 7542 and 512 GB memory.

## Software dependencies ##
Dependencies  | Version
------------- | -------------
nextflow | 25.10.2
openjdk | 17.0.8
Perl | 5.32.1
python | 3.11.0
fastp | 1.3.6
bowtie2 | 2.5.5
t1k | 1.0.9
HLA-HD | 1.7.1

## Installation ##
Installation should finish within 10 minutes (Epyc 7542).  
(1) Setup the environment using conda or mamba:  
```
mamba create -n EvaHLA -c bioconda -c conda-forge python=3.11.0 openjdk=17.0.8 nextflow=25.10.2 fastp=1.3.6 bowtie2=2.5.5 t1k=1.0.9
```    
Then activate the conda environment:  
```
mamba activate EvaHLA
```
(2) Clone the repository with `git clone`, and execute
```
cd EvaHLA
```
(3) Install HLA-HD. Refer to their manuals to install and prepare HLA index files.
(4) Prepare HLA index files for t1k. Refer to their github pages.

## Overview ##
HLA-HD and t1k are used to perform HLA genotyping on EvaHLA libraries. We use a hybrid comparison method to decide the final HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1 and HLA-DRB3/4/5.

## Usage ##
(1) Prepare the `samplesheet.csv`. The csv file __must__ contain 3 columns with defined column names:  
`sample`: Name of the sequenced library. For example, `demo-1`. It will be the prefix of the output. Note: Different fastqs with same sample name will be merged before processing.  
`fastq_1`: Path to read 1.  
`fastq_2`: Path to read 2.  

(2) Change directory to EvaHLA with `cd EvaHLA`, and execute:
```
nextflow run main.nf \
  -output-dir your_output_dir \
  --input_csv samplesheet.csv \
  --run_hlahd true \
  --hlahd_linenum 400000 \
  --hlahd_refdir path_to_hlahd_basedir \
  --run_t1k true \
  --t1k_preset hla-wgs \
  --t1k_reffile hlaidx_dna_seq.fa \
  -bg -resume
```
`-output-dir`: Path to the output directory.  
`--input_csv`: Path to samplesheet.csv as described in **step (1)**.  
`--hostGenome`: the `--assembly` parameter for `pairtools parse`.  
`--run_hlahd`: Perform genotyping with HLA-HD. Default: true. Valid options: true, false.   
`--hlahd`: Path to hlahd.sh. If the PATH is already exported, path to hlahd.sh will be automatically detected using `which hlahd.sh`. Otherwise, specify the path to hlahd.sh explicitly.  
`--hlahd_linenum`: Number of lines to extract from fastq files as input to HLA-HD. Since HLA-HD is time-consuming with large targeted HLA dataset, we subset fastq to perform genotyping. Default: 400000 lines. Typical runtime: 6 hours. Set to 0 to disable subsetting.  
`--hlahd_refdir`: Directory that stores __HLA_gene.split.txt__, __freq_data/__ and __dictionary/__. Default: directory two levels up from the path of __--hlahd__.   
`--run_t1k`: Perform genotyping with t1k. Default: true. Valid options: true, false.  
`--t1k_preset`: Preset of t1k.  
`--t1k_reffile`: Index file of t1k.  

## Expected output ##
Go to `-output-dir`. The pipeline will create folders named according to the `sample` column in you csv file. Each folder contain 3 subfolders:  
`fastqs`: Merged and adapter-trimmed fastqs.  
`HLA_HD`: Genotyping results of HLA-HD.  
`t1k`: Genotyping results of t1k.  
if both HLA-HD and t1k are called, a merged putative result will be established under `-output-dir`.  













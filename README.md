# pipeline_salmonquant

## The pipeline performs the following:
   * Uses gff2fasta to convert a .gtf formatted geneset into a FASTA format,
     building a transcriptome.
   * Builds a salmon index from the FASTA formatted transcriptome generated
   * Quantifies .fastq formatted RNA-seq reads and
   * Generates gene expression estimates (TPM and counts) at the transcript and
     gene level, using  Salmon as an alignment-free expression estimation.

## Inputs needed
1. A geneset in a .gtf format. Unzip .gtf.gz files before running the pipeline_salmonquant.
2. RNA-seq reads in .fastq.gz format

## For accessing the pipeline through cgatflow command
Clone/download the repository
Copy the pipeline_salmonquant.py and pipeline_salmonquant folder to your cgat/cgat-flow/CGATPipelines/

## Configuration
The pipeline requires a configured :file: ###`pipeline.yml` file.

Make a directory with your project name, for Salmon quantification
Configure the pipeline with ###`cgatflow salmonquant config`
A pipeline.log and pipeline.yml file(s) will be added to your new directory
Modify the pipeline.yml according to your project (specify genome, genome directory, annotation database and directory, database for uploading the outputs; specify options for Salmon quantification)

## Pipeline use
Run the pipeline with ###`cgatflow salmonquant make full`

For running the pipeline on a large set of samples, submit the pipeline onto the cluster (sharc), using a submit_pipeline_cgtaflow custom script.




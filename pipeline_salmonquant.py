"""=====================================
RNA-Seq Differential expression pipeline
========================================


The Salmon Quantification pipeline is a tool for quantifying
the expression of transcripts using RNA-seq data.

It requires three inputs:

   1. A geneset in :term:`gtf` formatted file - for building a reference transcriptome
   2. A transcriptome in :term:`fasta` formatted file - for building a salmon index
   3. RNA-seq reads in :term:`fastq.1.gz` and `fastq.2.gz` formatted files - for quantification

This pipeline works on a single genome.

Overview
========

The pipeline performs the following:

   * Uses gff2fasta to convert a .gtf formatted geneset into a FASTA format,
     building a transcriptome.

   * Builds a salmon index from the FASTA formatted transcriptome generated

   * Quantifies .fastq formatted RNA-seq reads and

   * Generates gene expression estimates (TPM and counts) at the transcript and
     gene level, using  Salmon as an alignment-free expression estimation.

Background
============

Quantification:

Transcripts are the natural choice to measure expression. However
other quantities might be of interest. Some quantities are biological
meaningful, for example differential expression from a promotor shared
by several trancripts. Other quantities might not be biologically
meaningful but are necessary as a technical comprise.  For example,
the overlapping transcripts might be hard to resolve and thus might
need to be aggregated per gene. Furthermore, functional annotation is
primarily associated with genes and not individual transcripts.

This pipeline estimates transcript and gene-level expression and
performs differential expression analysis on both levels.

The quantification tools fall into two categories:
   * Alignment-free
      Quantification is performed directly from the raw reads against
      a reference transcriptome using "pseduo-alignment". In essence,
      the tool attempts to identify the compatible transcripts for
      each read without identifying the exact alignment position of
      the read on the transcript or genome. Following this the set of
      expression estimates which best explain the observed reads are
      obtained using an Expectation Maximisation approach

      Some of the available tools are:
      * kallisto_
      * salmon_
      * sailfish_

== This pipeline will only use Salmon for quantification. ==

   * Alignment-based
      Quantification is performed on the aligned reads using the
      position of features described in the reference geneset
      gtf. Reads are discretly assigned to one feature (may be
      performed at the transcript or gene-level).  It should be noted
      that transcript-level quantification with tag counting methods
      is inherrently inaccurate since a read which aligns to an exon
      present in multiple isoforms of the same gene can not be naively
      assigned to a single transcript.

      Some of the available tools are:
      * featurecounts_
      * gtf2table (in-house script)

The alignment-free methods should be preffered over featureCounts and
gtf2table in almost all instances. However, many analyses still use
tag counting so it may be neccessary to repeat other groups
analyses. In addition gtf2table provides more detailed tag counting
which may be useful when exploring problematic RNA-Seq
data. Alignment-free methods also provide estimated counts per
transcript which can be rounded to integer counts.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use cgat pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.yml` file (see :ref:`PipelineReporting`). To start
with, use the files supplied with the Example_ data.

Input
-----

Geneset
++++++++

1. The Geneset is specified by the "geneset" parameter
2. The generated .fa transcriptome
3. RNA-seq read in .fastq.gz format

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default cgat setup, the pipeline requires the following
software to be in the path:

+--------------+----------+------------------------------------+
|*Program*     |*Version* |*Purpose*                           |
+--------------+----------+------------------------------------+
|gff2fasta     |          |bulding transcriptome in .fa format |
+--------------+----------+------------------------------------+
|salmon_       |>=0.7.2   |building an index                   |
+--------------+----------+------------------------------------+
|salmon_       |>=0.7.2   |alignment-free quantification       |
+--------------+----------+------------------------------------+


Pipeline output
===============

Quantification
--------------

The quantification estimates from each method are outputted to:
[method].dir/[sample]/[level].tsv.gz,
where [method] is the quantification method, [sample] is the sample
name and [level] is the feature level (transcript or gene)

Each tool also generates specific log files etc which are outputted,
along with the raw quantification outfile in the directory:
[method].dir/[sample]

For each method, the merged counts are outputted to:
[method].dir/[level].tsv.gz

Glossary
========

.. glossary::

salmon
      salmon_ - alignment-free quantification

.. _salmon: https://combine-lab.github.io/salmon/

###########################################################################

Code
====

"""
# load modules
from ruffus import *

import sys
import os
import sqlite3

from cgatcore import pipeline as P
import cgat.GTF as GTF
import cgatcore.iotools as IOTools

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file

P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.

PARAMS = P.PARAMS
if os.path.exists(PARAMS["annotations_dir"]):
    PARAMS.update(P.peek_parameters(
        PARAMS["annotations_dir"],
        'genesets',
        prefix="annotations_",
        update_interface=True,
        restrict_interface=True))

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import cgatpipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except KeyError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise.

###############################################################################
# Utility function
###############################################################################

def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.
    This method also attaches to helper databases.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

###############################################################################
# build indexes
###############################################################################

@mkdir('geneset.dir')
@transform(PARAMS['geneset'],
           regex("(\S+).gtf.gz"),
           r"geneset.dir/\1.fa")
def buildReferenceTranscriptome(infile, outfile):
    '''
    Builds a reference transcriptome from the provided GTF geneset - generates
    a fasta file containing the sequence of each feature labelled as
    "exon" in the GTF.
    --fold-at specifies the line length in the output fasta file

    Parameters
    ----------
    infile: str
        path to the GTF file containing transcript and gene level annotations
    genome_dir: str
        :term: `PARAMS` the directory of the reference genome
    genome: str
        :term: `PARAMS` the filename of the reference genome (without .fa)
    outfile: str
        path to output file
    '''

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    statement = '''
    zcat %(infile)s |
    awk '$3=="exon"'|
    cgat gff2fasta
    --is-gtf --genome-file=%(genome_file)s --fold-at=60 -v 0
    --log=%(outfile)s.log > %(outfile)s;
    samtools faidx %(outfile)s
    '''

    P.run(statement)


@active_if(PARAMS["salmon_use_genome_decoy"])
@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           add_inputs(
               os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"),
               PARAMS.get("addition_transcripts", "")),
           r"_with_decoy.fa")
def add_decoy_to_transcriptome(infiles, outfile):
    '''Modern versions of salmon allow the addition of decay sequences to
    the index.  most comprehensive option is to add the whole genome
    as a decoy (after the transcriptome). This function does this if
    the :PARAM:use_genome_decoy is True

    '''

    transcriptome, decoy_sequence, additional_transcripts = infiles

    statement = ''' cat %(transcriptome)s 
                        %(decoy_sequence)s 
                        %(additional_transcripts)s > %(outfile)s'''
    P.run(statement)


@active_if(PARAMS["salmon_use_genome_decoy"])
@transform(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"),
           formatter(),
           "decoys.txt")
def build_decoy_list(infile, outfile):
    '''The salmon indexer also requres a list of the transcripts to use as a
    decoy. This can be parsed from the fasta'''

    statement = '''grep "^>" %(infile)s
                     | cut -d " " -f 1
                     | sed 's/>//g' > %(outfile)s '''

    P.run(statement)


if PARAMS["salmon_use_genome_decoy"]:
    transcriptome = add_decoy_to_transcriptome
else:
    transcriptome = buildReferenceTranscriptome


@transform(transcriptome,
           suffix(".fa"),
           add_inputs(build_decoy_list
                      if PARAMS["salmon_use_genome_decoy"]
                      else None),
           ".salmon.index")
def buildSalmonIndex(infiles, outfile):
    '''
    Builds a salmon index for the reference transriptome
    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    salmon_kmer: int
       :term: `PARAMS` kmer size for sailfish.  Default is 31.
       Salmon will ignores transcripts shorter than this.
    salmon_index_options: str
       :term: `PARAMS` string to append to the salmon index command to
       provide specific options e.g. --force --threads N
    outfile: str
       path to output file
    '''

    job_memory = PARAMS["salmon_index_memory"]
    job_threads = PARAMS["salmon_index_threads"]

    transcriptome, decoys = infiles

    salmon_index_options = ""

    if PARAMS["salmon_index_options"] is not None:
        salmon_index_options += " " + PARAMS["salmon_index_options"]

    if decoys is not None:
        salmon_index_options += " -d " + decoys

    # need to remove the index directory (if it exists) as ruffus uses
    # the directory timestamp which wont change even when re-creating
    # the index files
    statement = '''
    rm -rf %(outfile)s &&
    salmon index %(salmon_index_options)s
                 -t %(transcriptome)s
                 -i %(outfile)s
                 -k %(salmon_kmer)s
                 -p %(salmon_index_threads)s
    '''

    P.run(statement)


@originate("transcript2geneMap.tsv")
def getTranscript2GeneMap(outfile):
    ''' Extract a 1:1 map of transcript_id to gene_id from the geneset '''

    iterator = GTF.iterator(IOTools.open_file(PARAMS['geneset']))
    transcript2gene_dict = {}

    for entry in iterator:

        try:
            gene_id = entry[PARAMS["gene_id_field"]]
        except KeyError:
            gene_id = entry.gene_id

        try:
            transcript_id = entry[PARAMS["transcript_id_field"]]
        except KeyError:
            transcript_id = entry.transcript_id

        # Check the same transcript_id is not mapped to multiple gene_ids!
        if transcript_id in transcript2gene_dict:
            if not gene_id == transcript2gene_dict[transcript_id]:
                raise ValueError('''multipe gene_ids associated with
                the same transcript_id %s %s''' % (
                    gene_id,
                    transcript2gene_dict[transcript_id]))
        else:
            transcript2gene_dict[transcript_id] = gene_id

    with IOTools.open_file(outfile, "w") as outf:
        outf.write("transcript_id\tgene_id\n")
        for key, value in sorted(transcript2gene_dict.items()):
            outf.write("%s\t%s\n" % (key, value))


###################################################
# alignment-free quantifier - salmon
###################################################

@follows(mkdir("quantification.dir"), buildSalmonIndex)
@transform(["*.fastq.gz",
            "*.fastq.1.gz"],
           formatter(".+/(?P<TRACK>.+).fastq.+"),
           add_inputs(buildSalmonIndex),
           r"quantification.dir/{TRACK[0]}/quant.sf")
def quantifyWithSalmon(infiles, outfile):
    '''Quantify existing samples against genesets

    Computes read counts across transcripts and genes based on a fastq
    file and an indexed transcriptome using Salmon.

    Runs the salmon "quant" function across transcripts with the specified
    options.  Read counts across genes are counted as the total in all
    transcripts of that gene.

    Parameters
    ----------
    infiles: list
        list with three components
        0 - list of strings - paths to fastq files to merge then quantify
        across using sailfish
        1 - string - path to salmon index file

    salmon_threads: int
       :term: `PARAMS` the number of threads for salmon
    salmon_memory: str
       :term: `PARAMS` the job memory for salmon
    salmon_options: str
       :term: `PARAMS` string to append to the salmon quant command to
       provide specific
       options, see http://sailfish.readthedocs.io/en/master/salmon.html
    salmon_bootstrap: int
       :term: `PARAMS` number of bootstrap samples to run.
       Note, you need to bootstrap for differential expression with sleuth
       if there are no technical replicates. If you only need point estimates,
       set to 1.
    salmon_libtype: str
       :term: `PARAMS` salmon library type
       as for sailfish - use
       http://sailfish.readthedocs.io/en/master/library_type.html#fraglibtype
    outfiles: list
       paths to output files for transcripts and genes

    '''

    job_threads = PARAMS["salmon_threads"]
    job_memory = PARAMS["salmon_memory"]

    infile, salmonIndex = infiles
    basefile = os.path.basename(infile)
    outfile = os.path.dirname(outfile)

    if infile.endswith(".fastq.1.gz"):
        fastq1 = infile
        fastq2 = P.snip(infile, ".1.gz")+".2.gz"
        fastq_inputs = ''' -1 %(fastq1)s
                           -2 %(fastq2)s ''' % locals()
    else:
        fastq_inputs = "-r %(infile)s" % locals()

    salmon_options = ""

    if PARAMS["salmon_options"] is not None:
        salmon_options += " " + PARAMS["salmon_options"]

    if PARAMS["salmon_bootstraps"] > 0:
        if PARAMS["salmon_bootstrap_type"] in ("gibbs", "Gibbs"):
            salmon_options += " --numGibbsSamples " + \
                                        str(PARAMS["salmon_bootstraps"])
        elif PARAMS["salmon_bootstrap_type"] in ("bootstrap", "Bootstrap"):
            salmon_options += " --numBootstraps " + \
                                        str(PARAMS["salmon_bootstraps"])

    statement = '''
    salmon quant -i %(salmonIndex)s
        --libType %(salmon_libtype)s
        %(fastq_inputs)s
        -o %(outfile)s
        --threads %(job_threads)s
        %(salmon_options)s &&

    cp %(outfile)s/quant.sf %(outfile)s.sf;
    '''

    P.run(statement)

@collate(quantifyWithSalmon,
         formatter(),
         add_inputs(getTranscript2GeneMap),
         ["transcripts_tximport.RData",
          "genes_tximport.RData"])
def summarize_with_tximport(infiles, outfiles):
    '''Use tximport to build R objects containing the summarized
    results of the salmon quantification. Should probably recreate
    with tximeta'''

    pipeline_dir = os.path.dirname(__file__)

    job_memory="16G"
    statement = '''Rscript %(pipeline_dir)s/tximport.R
                   > tximport.log'''

    P.run(statement)

###################################################
# Loading into database
###################################################

@follows(quantifyWithSalmon)
@merge("quantification.dir/*.sf", "salmon_quant.load")
def mergeAllQuants(infiles, outfile):
    job_memory = "6G"
    P.concatenate_and_load(infiles, outfile,
                           regex_filename="quantification.dir/(.+).sf",
                           options="-i transcript_id"
                           " -i Name -i Length -i EffectiveLength"
                           " -i TPM -i NumReads -i track"
                           " -i source")

###################################################
# Generic pipeline tasks

@follows(mergeAllQuants,
         getTranscript2GeneMap,
         summarize_with_tximport)
def full():
    ''' Builds transcriptome and index, then quantifies'''

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


###################################################################
# Create tasks
###################################################################

#TARGETS = []
#Targets = {'connect': (connect,),
#                       'buildReferenceTranscriptome': (buildReferenceTranscriptome,),
#                       'buildSalmonIndex': (buildSalmonIndex,),
#                       'getTranscript2GeneMap': (getTranscript2GeneMap,),
#                       'quantifyWithSalmon': (quantifyWithSalmon,),
#		       'full': (full,)
#                       }

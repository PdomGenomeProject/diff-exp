#!/usr/bin/make -f
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2015   Indiana University
#
# This file is part of the Polistes dominula genome project and is licensed
# under the Creative Commons Attribution 4.0 International License.
# -----------------------------------------------------------------------------

SHELL := bash
SRAURL=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR

.SECONDARY:

all: pdom-rnaseq-q1-1.fq pdom-rnaseq-q1-2.fq \
     pdom-rnaseq-q2-1.fq pdom-rnaseq-q2-2.fq \
     pdom-rnaseq-q3-1.fq pdom-rnaseq-q3-2.fq \
     pdom-rnaseq-q4-1.fq pdom-rnaseq-q4-2.fq \
     pdom-rnaseq-q5-1.fq pdom-rnaseq-q5-2.fq \
     pdom-rnaseq-q6-1.fq pdom-rnaseq-q6-2.fq \
     pdom-rnaseq-w1-1.fq pdom-rnaseq-w1-2.fq \
     pdom-rnaseq-w2-1.fq pdom-rnaseq-w2-2.fq \
     pdom-rnaseq-w3-1.fq pdom-rnaseq-w3-2.fq \
     pdom-rnaseq-w4-1.fq pdom-rnaseq-w4-2.fq \
     pdom-rnaseq-w5-1.fq pdom-rnaseq-w5-2.fq \
     pdom-rnaseq-w6-1.fq pdom-rnaseq-w6-2.fq

pdom-rnaseq-q1-1.fq: SRR2133905.fastq
	@ ln -s $< $@
pdom-rnaseq-q1-2.fq: SRR2133906.fastq
	@ ln -s $< $@
pdom-rnaseq-q2-1.fq: SRR2133907.fastq
	@ ln -s $< $@
pdom-rnaseq-q2-2.fq: SRR2133909.fastq
	@ ln -s $< $@
pdom-rnaseq-q3-1.fq: SRR2133910.fastq
	@ ln -s $< $@
pdom-rnaseq-q3-2.fq: SRR2133911.fastq
	@ ln -s $< $@
pdom-rnaseq-q4-1.fq: SRR2133913.fastq
	@ ln -s $< $@
pdom-rnaseq-q4-2.fq: SRR2133914.fastq
	@ ln -s $< $@
pdom-rnaseq-q5-1.fq: SRR2133915.fastq
	@ ln -s $< $@
pdom-rnaseq-q5-2.fq: SRR2133916.fastq
	@ ln -s $< $@
pdom-rnaseq-q6-1.fq: SRR2133917.fastq
	@ ln -s $< $@
pdom-rnaseq-q6-2.fq: SRR2133918.fastq
	@ ln -s $< $@
pdom-rnaseq-w1-1.fq: SRR2133892.fastq
	@ ln -s $< $@
pdom-rnaseq-w1-2.fq: SRR2133893.fastq
	@ ln -s $< $@
pdom-rnaseq-w2-1.fq: SRR2133895.fastq
	@ ln -s $< $@
pdom-rnaseq-w2-2.fq: SRR2133896.fastq
	@ ln -s $< $@
pdom-rnaseq-w3-1.fq: SRR2133897.fastq
	@ ln -s $< $@
pdom-rnaseq-w3-2.fq: SRR2133898.fastq
	@ ln -s $< $@
pdom-rnaseq-w4-1.fq: SRR2133899.fastq
	@ ln -s $< $@
pdom-rnaseq-w4-2.fq: SRR2133900.fastq
	@ ln -s $< $@
pdom-rnaseq-w5-1.fq: SRR2133901.fastq
	@ ln -s $< $@
pdom-rnaseq-w5-2.fq: SRR2133902.fastq
	@ ln -s $< $@
pdom-rnaseq-w6-1.fq: SRR2133903.fastq
	@ ln -s $< $@
pdom-rnaseq-w6-2.fq: SRR2133904.fastq
	@ ln -s $< $@

%.fastq: %.sra
	@ echo [Convert accession $* from .sra to .fastq]
	@ which fastq-dump > /dev/null
	@ fastq-dump $<

%.sra:
	@ echo [Download accession $* from SRA]
	@ PREFIX=$$(echo $* | cut -c 1-6); wget ${SRAURL}/$$PREFIX/$*/$*.sra

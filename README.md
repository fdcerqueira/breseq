# SNVs from *E. coli* populations

The bash pipeline was made for experiments designed to study the tempo and mode of adaptation of a new strain of *E. coli* to the gut of obese mice, leptin deficient mice were colonized with a strain of *E. coli*. In addition to the invader strain, the mice were also colonized with another *E. coli* strain which typically resides in their normal gut microbiota.

Bash pipeline to process fastq files from whole-genome sequencing of *E. coli* populations, bbsplit.sh to remove host contamination, two plasmids from one of the *E. coli* strains  and run breseq. I may be run with merged pair-end reads.
(br/)
(br/)
(br/)
**breseq2.sh:**


1)Check and correct if necessary the names of the fastq files

1)fastp for reads QC

2)filter reads from the genomes used as reference

3)flash to merge pair-end reads (optional)

4)filter reads with bbsplit.sh (it is possible to add other genomes that are possible contaminants in the script options)

5)Run breseq

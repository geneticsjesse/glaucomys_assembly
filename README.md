# A <em> de novo </em> genome assembly and annotation of the southern flying squirrel (<em> Glaucomys volans</em>) (Wolf et al. 2021)
### This repository a bash script used to generate a <em> de novo </em> genome assembly for southern flying squirrels that was published in G3 Genes|Genomes|Genetics. Detailed below are the steps this script walks you through, starting with raw sequence data.

If you are interested in the manuscript, see https://doi.org/10.1093/g3journal/jkab373

1. Run FastQC on raw long linked reads
1. Run the Supernova genome assembler on the long linked reads (https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/welcome)
1. Create the fasta file from the assembly file
#### Note: I chose the raw output style but there are multiple options)
3. Run FastQC on raw short reads
4. Concatenate short reads across lanes and re-run fastqc
#### Note: The short reads were used as an extension of the paper mentioned above as I re-sequenced the short-read pairs
using the <em> de novo </em> assembly and generated some comparative analyses)
1. Trim the concatenated short reads using trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic) 
#### Note: There are various customizable paramaters here and will require specific input based on individual data)
3. Create a kraken database from which to screen reads for contaminants (http://ccb.jhu.edu/software/kraken/)
4. Check for contaminant reads using Kraken database
5. Count number of reads after trimming and contamination screening
6. Run BUSCO to assess genome completedness (https://busco.ezlab.org/)
7. Generate a repeat-masked version of the genome using RepeatMasker (https://www.repeatmasker.org/)
8. Annotate the genome using AUGUSTUS (https://bioinf.uni-greifswald.de/augustus/) both with and without hints. I generated hints using a transcriptome that was previously generated for flying squirrels

#1) Filtering raw reads, aligning to reference genome, creating and filtering bam files, and producing VCF files
#First need to trim raw reads, removing any adaptors and low quality bases
#Make sure you use the correct adaptor file - here we are using the 'TruSeq3-PE-2.fa'.
#Any read smaller than 36bp after trimming will be discarded

nohup java -jar /path/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 ./Name_R1.fastq.gz ./Name_R2.fastq.gz Name_R1_Trimmed.fq.gz Name_R1_Unpaired.fq.gz \
Name_R2_Trimmed.fq.gz Name_R2_Unpaired.fq.gz ILLUMINACLIP:/path/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36 &
 
#Align to the reference genome using bowtie2 

nohup bowtie2 --local -p 8 --un-conc ./FailedAlign_Name -x fasta_pseudohap.fasta \
-q -1 ./Name_R1_Trimmed.fq.gz -2 ./Name_R2_Trimmed.fq.gz -S NameMapped.sam &

#Convert Sam to Bam file
nohup samtools view -bS NameMapped.sam > NameMapped.bam &

#Sort the bam file
nohup samtools sort NameMapped.bam -o Name_Sorted.bam &

#Filter for mapping quality of reads to remove poorly mapped ones
nohup samtools view -b -q 20 Name_Sorted.bam > Name_Sorted_filtered.bam &

#Add Read Group information (Sones on a bigger cluster)
#To find out RGLB info, use:
zcat SFSCC1_1.fastq.gz | head
474


#Example of first line for NFS_50254
@E00389:474:HCF7VCCX2:1:1101:7649:993 1:N:0
#Here, the number we need for RGLB = 474
nohup java -jar picard.jar AddOrReplaceReadGroups I=SFSCC1_Sorted_filtered.bam O=SFSCC1_SortedRG_filtered.bam RGLB=474 RGPL=illumina RGPU=none RGSM=SFSCC1 &
nohup java -jar picard.jar AddOrReplaceReadGroups I=SFS_25428_Sorted_filtered.bam O=SFS_25428_SortedRG_filtered.bam RGLB=474 RGPL=illumina RGPU=none RGSM=SFS_25428 &
nohup java -jar picard.jar AddOrReplaceReadGroups I=NFS_6525_Sorted_filtered.bam O=NFS_6525_SortedRG_filtered.bam RGLB=474 RGPL=illumina RGPU=none RGSM=NFS_6525 &
nohup java -jar picard.jar AddOrReplaceReadGroups I=NFS_50254_Sorted_filtered.bam O=NFS_50254_SortedRG_filtered.bam RGLB=474 RGPL=illumina RGPU=none RGSM=NFS_50254 &

#Remove duplicates
nohup java -jar picard.jar MarkDuplicates I=SFSCC1_SortedRG_filtered.bam REMOVE_DUPLICATES=true O=SFSCC1_NoDups.bam M=marked_dup_metricsSFSCC1.txt &
nohup java -jar picard.jar MarkDuplicates I=SFS_25428_SortedRG_filtered.bam REMOVE_DUPLICATES=true O=SFS_25428_NoDups.bam M=marked_dup_metricsSFS_25428.txt &
nohup java -jar picard.jar MarkDuplicates I=NFS_6525_SortedRG_filtered.bam REMOVE_DUPLICATES=true O=NFS_6525_NoDups.bam M=marked_dup_metricsNFS_6525.txt &
nohup java -jar picard.jar MarkDuplicates I=NFS_50254_SortedRG_filtered.bam REMOVE_DUPLICATES=true O=NFS_50254_NoDups.bam M=marked_dup_metricsNFS_50254.txt &

#Clip overlapping regions 
nohup ./bam clipOverlap --in /home/ubuntu/vol2/6_mark_dup/SFSCC1_NoDups.bam --out /home/ubuntu/vol2/7_clipoverlap/SFSCC1_NoDups_clip.bam --stats &
nohup ./bam clipOverlap --in /home/ubuntu/vol2/6_mark_dup/SFS_25428_NoDups.bam --out /home/ubuntu/vol2/7_clipoverlap/SFS_25428_NoDups_clip.bam --stats &
nohup ./bam clipOverlap --in /home/ubuntu/vol2/6_mark_dup/NFS_6525_NoDups.bam --out NFS_6525_NoDups_clip.bam --stats > NFS_6525_NoDups_clip.out &
nohup ./bam clipOverlap --in /home/ubuntu/vol2/6_mark_dup/NFS_50254_NoDups.bam --out NFS_50254_NoDups_clip.bam --stats > NFS_50254_NoDups_clip.out &

#Indel realignment step 1 - requires the bam file to be indexed:
nohup samtools index SFSCC1_NoDups_clip.bam &
nohup samtools index SFS_25428_NoDups_clip.bam &
nohup samtools index NFS_6525_NoDups_clip.bam &
nohup samtools index NFS_50254_NoDups_clip.bam &

#Then run step 1 

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R /home/ubuntu/vol1/fasta_files/sfs_reference.fasta \
-I ./home/ubuntu/vol2/7_clipoverlap/SFSCC1_NoDups_clip.bam \
-o SFSCC1_indel_realigner.intervals \
-drf BadMate

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R /home/ubuntu/vol2/sfs_reference.fasta \
-I ./home/ubuntu/vol2/7_clipoverlap/SFS_25428_NoDups_clip.bam \
-o SFS_25428_indel_realigner.intervals \
-drf BadMate

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R /home/ubuntu/vol2/sfs_reference.fasta \
-I ./home/ubuntu/vol2/7_clipoverlap/NFS_6525_NoDups_clip.bam \
-o NFS_6525_indel_realigner.intervals \
-drf BadMate

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R /home/ubuntu/vol2/sfs_reference.fasta \
-I ./home/ubuntu/vol2/7_clipoverlap/NFS_50254_NoDups_clip.bam \
-o NFS_50254_indel_realigner.intervals \
-drf BadMate

#And then re-align step 2:

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T IndelRealigner \
-R /home/ubuntu/vol2/sfs_reference.fasta \
-I ./home/ubuntu/vol2/7_clipoverlap/SFSCC1_NoDups_clip.bam \
-targetIntervals SFSCC1_indel_realigner.intervals \
--consensusDeterminationModel USE_READS  \
--nWayOut _realigned.bam

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T IndelRealigner \
-R /home/ubuntu/vol2/sfs_reference.fasta \
-I ./home/ubuntu/vol2/7_clipoverlap/SFS_25428_NoDups_clip.bam \
-targetIntervals SFS_25428_indel_realigner.intervals \
--consensusDeterminationModel USE_READS  \
--nWayOut _realigned.bam

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T IndelRealigner \
-R /home/ubuntu/vol2/sfs_reference.fasta\
-I ./home/ubuntu/vol2/7_clipoverlap/NFS_6525_NoDups_clip.bam \
-targetIntervals NFS_6525_indel_realigner.intervals \
--consensusDeterminationModel USE_READS  \
--nWayOut _realigned.bam

nohup java -jar ./gatk/GenomeAnalysisTK.jar -T IndelRealigner \
-R /home/ubuntu/vol2/sfs_reference.fasta \
-I ./home/ubuntu/vol2/7_clipoverlap/NFS_50254_NoDups_clip.bam \
-targetIntervals NFS_50254_indel_realigner.intervals \
--consensusDeterminationModel USE_READS  \
--nWayOut _realigned.bam

#Re-sort the bam file 

nohup samtools sort SFSCC1_NoDups_clip_realigned.bam -o Name_NoDups_clip_realigned_sorted.bam &
nohup samtools sort SFS_25428_NoDups_clip_realigned.bam -o Name_NoDups_clip_realigned_sorted.bam &
nohup samtools sort NFS_6525_NoDups_clip_realigned.bam -o Name_NoDups_clip_realigned_sorted.bam &
nohup samtools sort NFS_50254_NoDups_clip_realigned.bam -o Name_NoDups_clip_realigned_sorted.bam &

#Calculate Average Depth of genome wide coverage for the BAM file
nohup samtools depth -a SFSCC1_re_sorted.bam > SFSCC1.coverage & 
nohup awk '{sum += $3} END {print "Average = ", sum/NR}' SFSCC1.coverage > SFSCC1.Average &
#Once you have got the SFSCC1.Average file, you can delete the (huge) SFSCC1.coverage file. 

#Calculate Average Depth of genome wide coverage for the BAM file
nohup samtools depth -a SFS_25428_re_sorted.bam > SFS_25428.coverage & 
nohup awk '{sum += $3} END {print "Average = ", sum/NR}' SFS_25428.coverage > SFS_25428.Average &
#Once you have got the SFS_25428.Average file, you can delete the (huge) SFS_25428.coverage file. 

#Calculate Average Depth of genome wide coverage for the BAM file
nohup samtools depth -a NFS_6525_re_sorted.bam > NFS_6525.coverage & 
nohup awk '{sum += $3} END {print "Average = ", sum/NR}' NFS_6525.coverage > NFS_6525.Average &
#Once you have got the NFS_6525.Average file, you can delete the (huge) NFS_6525.coverage file. 

#Calculate Average Depth of genome wide coverage for the BAM file
nohup samtools depth -a NFS_50254_re_sorted.bam > NFS_50254.coverage & 
nohup awk '{sum += $3} END {print "Average = ", sum/NR}' NFS_50254.coverage > NFS_50254.Average &
#Once you have got the NFS_50254.Average file, you can delete the (huge) NFS_50254.coverage file. 

#Now you have your BAM files which can be used for some analyses. Others need a VCF file. 
#Make bam index file  

nohup java -jar $EBROOTPICARD/path/Picard/picard/build/libs/picard.jar BuildBamIndex I=./SFSCC1_NoDups_clip_realigned_sorted.bam &
nohup java -jar $EBROOTPICARD/path/Picard/picard/build/libs/picard.jar BuildBamIndex I=./SFS_25428_NoDups_clip_realigned_sorted.bam &
nohup java -jar $EBROOTPICARD/path/Picard/picard/build/libs/picard.jar BuildBamIndex I=./NFS_6525_NoDups_clip_realigned_sorted.bam &
nohup java -jar $EBROOTPICARD/path/Picard/picard/build/libs/picard.jar BuildBamIndex I=./NFS_50254_NoDups_clip_realigned_sorted.bam &

#Make individual VCF files
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/ubuntu/vol1/fasta_files/sfs_genome.fasta -I Name_NoDups_clip_realigned_sorted.bam -o Name_NoDups_clip_realigned_sorted.g.vcf.gz -ERC GVCF

#Combine VCF files
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CombineGVCFs -R /home/ubuntu/vol1/fasta_files/sfs_genome.fasta --variant Name_NoDups_clip_realigned_sorted.g.vcf.gz --variant Name_NoDups_clip_realigned_sorted.g.vcf.gz --variant Name_NoDups_clip_realigned_sorted.g.vcf.gz --variant Name_NoDups_clip_realigned_sorted.g.vcf.gz -o CombinedVCFName.g.vcf.gz

#Perform joint genotyping
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /home/ubuntu/vol1/fasta_files/sfs_genome.fasta -V CombinedVCFName.g.vcf.gz -o CombinedVCFName_Genotyped.g.vcf.gz

#Filter VCF files - these may need to be adjusted for your needs. These steps you can do on the smaller cluster as VCFtools is quick to run.
#Step 1 - remove indels and calls with low scores (which are changed to missing). Also filter for minimum and max depth (adjusted so max is double the average depth, here set to 40 for genomes where the average depth is 20)
nohup vcftools --gzvcf ./CombinedVCF_FS_Genotyped.g.vcf.gz --remove-indels --minGQ 20 --min-meanDP 5 --max-meanDP 33 --minQ 20 --recode --recode-INFO-all --out CombinedVCF_FS_Genotyped_FilterStep1 &

#Step 2 - filter for missing data. How much missing you allow is up to you - here I am allowing 10% missing data (0.9)
nohup vcftools --gzvcf ./CombinedVCF_FS_Genotyped_FilterStep1.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out CombinedVCF_FS_Genotyped_FilterStep2 &

# countnumber of snps / find line number where to cut off for PCA
bcftools stats CombinedVCF_FS_Genotyped_FilterStep2.recode.vcf > finalvcf_stats
grep -n SFSCC1 CombinedVCF_FS_Genotyped_FilterStep2.recode.vcf > line_number.txt (just finds the line # where sample ID is)

#PSMC
generation time = 1.5 years
mutation rate/generation = 6.8  10 ^-5

#PCA prep
#remove 65,813 lines from vcf for PCA
sed '1,65,813d' CombinedVCF_FS_Genotyped_FilterStep2.recode.vcf > CombinedVCF_FS_Genotyped_FilterStep2_PCA.recode.vcf
#Then add the VCF header so the program knows what it is:
sed '1 i ##fileformat=VCFv4.2' -i nameoffile.vcf

# using angsd for heterozygosity 
./angsd -i my.bam -anc ref.fa -dosaf 1 -fold 1

nohup angsd/angsd -i NFS_50254_re_sorted.bam -gl 1 -anc sfs_reference.fasta -dosaf 1 -C 50 -ref sfs_reference.fasta -minQ 20 -minmapq 30 -out NFS50254_het &
nohup ./angsd -i /home/ubuntu/vol1/SFSCC1_re_sorted.bam -gl 1 -anc /home/ubuntu/vol1/sfs_reference.fasta -dosaf 1 -C 50 -ref /home/ubuntu/vol1/sfs_reference.fasta -minQ 20 -minmapq 30 -out SFSCC1_het &

And then:

nohup angsd/misc/realSFS NFS50254_het.saf.idx > NFS50254_het_est.ml

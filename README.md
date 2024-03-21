# BMMF_PBMC_breast_Project
BMMF_PBMC_breast_Projec

I downloaded the following .fastq files from the internet:
SRR14530610.fastq
SRR14530611.fastq
SRR14530612.fastq
SRR14530613.fastq
SRR14530614.fastq
SRR14530615.fastq
SRR14530616.fastq
SRR14530617.fastq

1) Quality check of.fastq files with fastqc in terminal

Terminal command:
fastqc *.fastq

2) Trim low quality reads with Trimmomatic-0.39 tool.

Terminal Command:
java -jar trimmomatic-0.39.jar SE -threads 4 inputFile.fastq outputFile_trimmed.fastq TRAILING:10 -phred33

3) Quality check of trimmed data.

Terminal Command:
fastqc outputFile_trimmed.fastq

4) Alignment of Reads with Reference Genome (grch38) using 
	HISAT2 tool. 
	
Terminal Command:
hisat2 -q --rna-strandness R -x  grch38/genome -U 	inputFile_trimmed.fastq | samtools sort -o 	outputFile_trimmed.bam

5) Use .bam file and Gene Transfer Format (GTF) file to build the 	Feature Count Matrix  (Subread)

Terminal Command:
featureCounts -S 1 -a Homo_sapiens.GRCh38.106.gtf -o 	outputFile_featurecounts.txt inputFile_trimmed.bam

I have 8 Count Matrix .txt files as output.
And combined them to one .csv file.




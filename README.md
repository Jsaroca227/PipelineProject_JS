## STEP 1: Preparing the Paired-End Reads from SRA

In order to process the SRA files into paired-end FASTQ reads, please follow t>

1. Retrieve the SRA files from NCBI:
- Go to the NCBI Data Access, and on the left of the search bar select 'SRA'
- Copy and paste the given SRA accession number
- Click on the SRR number under 'Runs'
- On the next screen, select 'Data access' tab
- Copy the 'SRA Normalized' downloadable link for each file

2. Set up the project directory:
- In your terminal -
 mkdir **directory name of your own liking**
cd **name of the directory made**
git init; this will initliaze Git repo

3. Download the SRA files using wget:
wget <SRA_normalized_link_1>
wget <SRA_normalized_link_2>

4. Convert the SRA files into a FASTQ using fasterq-dump:
fasterq-dump <SRA_file_1>
fasterq-dump <SRA_file_2>

5. Confirm the paried-end reads using ls:
ls
- This will print **hopefully** the paired end reads as so;
SRA_1.fastq SRA_2.fastq# PiplineProject

#STEP 1: Preparing the Paired-End Reads from SRA

1. Set up the project directory:
- git clone <GitHub_Repo_link>
- cd **name of the repository made**
- git init

2. Retrieve the SRA files from NCBI:
- Go to the NCBI Data Access, and on the left of the search bar select 'SRA'
- Copy and paste the given SRA accession number in the search bar and enter
- Click on the SRR number under 'Runs'
- On the next screen, select 'Data access' tab
- Copy the 'SRA Normalized' downloadable link for each file

3. Download the SRA files using wget:
wget <SRA_normalized_link_1>
wget <SRA_normalized_link_2>

4. Convert the SRA files into a FASTQ using fasterq-dump:
fasterq-dump <SRA_file_1>
fasterq-dump <SRA_file_2>

5. Confirm the paried-end reads using ls command:
- This will list the content in the directory
- If done correctly, there will be 2 output files per SRA

<SRA>_1.fastq <SRA>_2.fastq

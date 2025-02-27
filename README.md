Background:
- The purpose of this project is to automate through bioinformatic workflows by the integration of multiple software tools into a creation of a Python Wrapper script.

# STEP 1: Preparing the Paired-End Reads from SRA

1. Clone the repository:
- git clone <GitHub_Repo_link>
- cd **name of the repository made**

2. Retrieve the SRA files from NCBI:
- In the NCBI Data Access, and on the left of the search bar select 'SRA'
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

<SRA_1>_1.fastq <SRA_1>_2.fastq <SRA_2>_1.fastq <SRA_2>_2.fastq

# STEP 2: Install all necessary dependencies
- os
- SeqIO
- SPAdes
- bowtie2
- BLAST+
- NCBI Datasets

# STEP 3: Run the Pipeline_Wrapper
- To run the wrapper, ensure that you are in the correct directory, clone "PipelineProject_Junelle_Saroca"; https://github.com/Jsaroca227/PipelineProject_Junelle_Saroca.git
- Ensure that the right, correct inputs are provided into the terminal for each command in order to run the script properly
- When it asks to override the files in the dataset, please select "y" for everything but the README.md file select "n".


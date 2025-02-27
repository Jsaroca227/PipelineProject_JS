#######question 2########

import os
import argparse

#parse arguments, it will ask for these to be answered in the terminal prior to running the script
parser = argparse.ArgumentParser(description="Wrapper script for HCMV pipeline.")
parser.add_argument("--first_name", required=True)
parser.add_argument("--last_name", required=True)
parser.add_argument("--sra1_r1", required=True)
parser.add_argument("--sra1_r2", required=True)
parser.add_argument("--sra2_r1", required=True)
parser.add_argument("--sra2_r2", required=True)
args = parser.parse_args()

#sets up the project directory
project_dir = f"PipelineProject_{args.first_name}_{args.last_name}"
os.makedirs(project_dir, exist_ok=True)
os.chdir(project_dir)
log_file = "PipelineProject_Junelle_Saroca/PipelineProject.log"

#downloads the dataset
download_db = "datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report"
os.system(download_db)
#unzips the dataset
unzip = "unzip ncbi_dataset.zip"
os.system(unzip)

#builds the bowtie2 index
genome_fasta = "ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna"
bowtie2_command = f"bowtie2-build {genome_fasta} HCMV_index"

os.system(bowtie2_command)

#calculates the paired end reads before filtering
sra1_read_pairs = int(os.popen(f'wc -l < {args.sra1_r1}').read().strip())
before_sra1 = int(sra1_read_pairs) // 4
sra2_read_pairs = int(os.popen(f'wc -l < {args.sra2_r1}').read().strip())
before_sra2 = int(sra2_read_pairs) // 4


#runs bowtie2 mapping for each SRA with its respective paired end reads
bowtie2_1 = f"bowtie2 --quiet -x HCMV_index -1 {args.sra1_r1} -2 {args.sra1_r2} -S mapped_r1.sam"
os.system(bowtie2_1)

bowtie2_2 = f"bowtie2 --quiet -x HCMV_index -1 {args.sra2_r1} -2 {args.sra2_r2} -S mapped_r2.sam"
os.system(bowtie2_2)


#counts the mapped reads after bowtie2 filtering
sra1_mapped_reads = int(os.popen(f"samtools view -c -F 4 mapped_r1.sam").read().strip())
sra2_mapped_reads = int(os.popen(f"samtools view -c -F 4 mapped_r2.sam").read().strip())

#opens the log file and writes, appends the results to the file
with open (log_file, 'a') as log:
    log.write(f"Donor 1 (2dpi) had {before_sra1} read pairs before Bowtie2 filtering and {sra1_mapped_reads} read pairs after.\n")
    log.write(f"Donor 1 (6dpi) had {before_sra2} read pairs before Bowtie2 filtering and {sra2_mapped_reads} read pairs after.\n")


#########question 3##########
import os
import argparse

#parse arguments to run the file
parser = argparse.ArgumentParser(description="Run SPAdes assembly with Bowtie2 output.")
parser.add_argument("--sra1_r1", required=True, help="Path to first SRA file (R1).")
parser.add_argument("--sra1_r2", required=True, help="Path to first SRA file (R2).")
parser.add_argument("--sra2_r1", required=True, help="Path to second SRA file (R1).")
parser.add_argument("--sra2_r2", required=True, help="Path to second SRA file (R2).")
parser.add_argument("-k", type=int, default=99, help="K-mer size for SPAdes assembly (default: 99).")
parser.add_argument("-t", type=int, default=2, help="Number of threads for SPAdes (default: 2).")
parser.add_argument("-o", default="combined_assembly/", help="Output directory for SPAdes assembly (default: combined_assembly/).")

args = parser.parse_args()

#ensures the output directory exits
os.makedirs(args.o, exist_ok=True)

#the log file that will hold the answer
log_file = "PipelineProject.log"

#create the spades string
spades_command = (
    f"spades.py -k {args.k} -t {args.t} --only-assembler "
    f"-1 {args.sra1_r1} -2 {args.sra1_r2} "
    f"-1 {args.sra2_r1} -2 {args.sra2_r2} "
    f"-o {args.o}"
)

#run the spades assembly 
os.system(spades_command)

#open the log file and write and append the SPAdes command ran prior
with open(log_file, "a") as log:
    log.write(f"SPAdes assembly command: {spades_command}\n")

##########question 4#############
#parse arguments
parser = argparse.ArgumentParser(description="Count contigs >1000 bp and total bp in assembly.")
parser.add_argument("--contigs", required=True, help="Path to the contigs FASTA file")
parser.add_argument("--log", required=True, help="Path to the log file")

args = parser.parse_args()

#initalize varaibles that will hold the count for contigs and total base pairs
num_contigs = 0
total_bp = 0
sequence = ""

#opens the contig file and processes each line
with open(args.contigs, "r") as file:
    for line in file:
        line = line.strip()
        #if the line is a header, it will start with ('>')
        if line.startswith(">"):  
            #if the previous sequence is longer than 1000 bp, it will count it
            if len(sequence) > 1000:  
                num_contigs += 1
                total_bp += len(sequence)
            sequence = "" #resets the sequence for the next contig
        else:
            sequence += line  #appends the sequence line

#this checks the last sequence in the file after the loop and repeats the previous steps
#checks if it is longer than 1000 bp, it will count it
if len(sequence) > 1000:
    num_contigs += 1
    total_bp += len(sequence)

#use os.system to append the results of amount of contigs and bp in the assembly to the log file
os.system(f'echo "There are {num_contigs} contigs > 1000 bp in the assembly." >> {args.log}')
os.system(f'echo "There are {total_bp} bp in the assembly.\n" >> {args.log}')


###########question 5###############
from Bio import SeqIO
import os
import argparse

#parse arguments
parser = argparse.ArgumentParser(description="Find longest contig and run BLAST")
parser.add_argument("--contigs", required=True, help="Path to the contigs FASTA file")
parser.add_argument("--output_dir", required=True, help="Directory to save output files")
parser.add_argument("--log", default="PipelineProject.log", help="Log file path")

args = parser.parse_args()
#ensures that the directory exists
os.makedirs(args.output_dir, exist_ok=True)

#defines the output path files
longest_contig_file = os.path.join(args.output_dir, "longest_contig.fasta")
blast_output = os.path.join(args.output_dir, "blast.tsv")

#intialize new variable to hold the  longest contig
longest_contig = None
#intialize new variable for the maximum length
max_length = 0
#checks if there are existing contigs in the file
if os.path.exists(args.contigs):
    #loops throught the contig records in the input file
    for record in SeqIO.parse(args.contigs, "fasta"):
        #if the contig is longer than the previous longest contig, it will be updated
        if len(record.seq) > max_length:
            max_length = len(record.seq)
            longest_contig = record
    #if the longest contig was found, it will right it into the output file
    if longest_contig:
        with open(longest_contig_file, "w") as output:
            SeqIO.write(longest_contig, output, "fasta")
#downloading the refseq database
download_db = "datasets download virus genome taxon Betaherpesvirinae --include genome"
os.system(download_db)

#unzip the ncbi dataset
unzip = "unzip ncbi_dataset.zip"
os.system(unzip)

#rename the file to be more specific
os.system("mv ncbi_dataset/data/genomic.fna Betaherpesvirinae.fna")

#create a new BLAST database
db = "makeblastdb -in Betaherpesvirinae.fna -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl"
os.system(db)

#prepares the BLAST command to run with the longest contig as the query
blast_command = 'blastn -query problem5_contig_file.fasta -db Betaherpesvirinae -out myresults.tsv -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_target_seqs 10 -max_hsps 1'
#executes the BLAST command
os.system(blast_command)

with open('blast.tsv', 'r') as f:
    results =f.readlines()
    f.close()
#open and append results into the log file
with open(args.log, "a") as log: 
    log.write('sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n') #adds the column headers to the log file
    for i in results:
        log.write(i)

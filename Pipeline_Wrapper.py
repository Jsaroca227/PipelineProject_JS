#######question 2########
import os
from Bio import SeqIO

log_file = "PipelineProject.log"

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
sra1_read_pairs = int(os.popen(f'wc -l < sampleSRR5660030_1.fastq').read().strip())
before_sra1 = int(sra1_read_pairs) // 4
sra2_read_pairs = int(os.popen(f'wc -l < sampleSRR5660033_1.fastq').read().strip())
before_sra2 = int(sra2_read_pairs) // 4


#runs bowtie2 mapping for each SRA with its respective paired end reads
bowtie2_1 = f"bowtie2 --quiet -x HCMV_index -1 sampleSRR5660030_1.fastq -2 sampleSRR5660030_2.fastq -S mapped_r1.sam"
os.system(bowtie2_1)

bowtie2_2 = f"bowtie2 --quiet -x HCMV_index -1 sampleSRR5660033_1.fastq -2 sampleSRR5660033_2.fastq -S mapped_r2.sam"
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

#the log file that will hold the answer
log_file = "PipelineProject.log"

#create the spades string
spades_command = ("spades.py -k 99 -t 2 --only-assembler -1 sampleSRR5660030_1.fastq -2 sampleSRR5660030_2.fastq -1 sampleSRR5660033_1.fastq -2 sampleSRR5660033_2.fastq -o combined_assembly/"
)

#run the spades assembly 
os.system(spades_command)

#open the log file and write and append the SPAdes command ran prior
with open(log_file, "a") as log:
    log.write(f"SPAdes assembly command: {spades_command}\n")

##########question 4#############
import os

#gets the directory to where the file is located
script_dir = os.path.dirname(os.path.abspath(__file__))

#this constructs the relative path to the contigs.fasta file
contig_file = os.path.join(script_dir, "combined_assembly", "contigs.fasta")

#initalize variables that will hold the count for contigs and total base pairs
log_file = "PipelineProject.log"
num_contigs = 0
total_bp = 0
sequence = ""

#opens the contig file and processes each line
with open(contig_file, "r") as file:
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
with open(log_file, "a") as log:
    log.write(f"There are {num_contigs} contigs > 1000 bp in the assembly.\n")
    log.write(f"There are {total_bp} bp in the assembly.\n")

###########question 5###############
from Bio import SeqIO
import os

log_file = 'PipelineProject.log'
#defines the output path files
longest_contig_file = "longest_contig.fasta"
blast_output = "blast.tsv"

#gets the directory to where the file is located
script_dir = os.path.dirname(os.path.abspath(__file__))

#this constructs the relative path to the contigs.fasta file
contig_file = os.path.join(script_dir, "combined_assembly", "contigs.fasta")

#intialize new variable to hold the  longest contig
longest_contig = None
#intialize new variable for the maximum length
max_length = 0
#checks if there are existing contigs in the file
if os.path.exists(contig_file):
    #loops throught the contig records in the input file
    for record in SeqIO.parse(contig_file, "fasta"):
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
blast_command = 'blastn -query longest_contig.fasta -db Betaherpesvirinae -out blast.tsv -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_target_seqs 10 -max_hsps 1'
#executes the BLAST command
os.system(blast_command)

with open('blast.tsv', 'r') as f:
    results =f.readlines()
    f.close()
#open and append results into the log file
with open(log_file, "a") as log: 
    log.write('sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n') #adds the column headers to the log file
    for i in results:
        log.write(i)

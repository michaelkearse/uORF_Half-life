
# How to run the scripts  

In order to run the scripts, you will need to download the RefSeq Transcript and RefSeq Reference Genome Annotation files for human genome build 38 which can be found here (https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml). Additionally, a spreadsheet file will need to be downloaded that contains information about transcript half-lives. The file can be found here https://genome.cshlp.org/content/27/3/407.full by clicking on Supplemental Material then by downloading Supplemental_Table_S3.xlsx. Note that this file must be converted to a tab delimited text file in order to be compatible with the scripts. To run the script, make sure Python and Biopython are installed on your system. Also make sure to change the 2-3 lines of code that are commented as “needing to be changed” near the top of both scripts. generate_db.py must be run before calculate_uorf.py.  

  

# How the scripts work  

Two scripts were created to compute the uORF length and determine the half-life for the given transcripts. The first script utilized the RefSeq Transcript and RefSeq Reference Genome Annotation files for human genome build 38 to store relevant information about all known transcripts in a local database. The first script iterated through the transcriptome file and stored each sequence, along with its transcript id, into a local SQLite database. The transcript id was used as a primary key for the database. The script then iterated though the annotation file to calculate total exon length and CDS start position for all transcripts.  This information was stored in the same database. CDS start position was calculated by counting the number of nucleotides from the start of the first exon until it reached the CDS start position.  

A second script was created to iterate through the tab delimited text file containing the half-lives and transcript ids for the different genes. For each transcript id listed in a row, a query was sent to the local database created by the first script to get the information for the transcript. If multiple transcripts were listed the one with the longest length that did not contain an exon mismatch was chosen. Additionally, if a row only contained transcripts with exon mismatches or it could not find a transcript match in the database that gene was not included in the calculation. An exon mismatch is determined by comparing the length of all exons that make up a transcript in the annotation file to the length of the transcript in the transcriptome file. If these lengths are not equal, then a discrepancy between the sequence in the genome and transcriptome file exists. Thus, the CDS start position cannot be reliably determined as the calculation to determine CDS start position was based off the genomic coordinates in the annotation file. uORF length was calculated by cycling through all three reading frames for the given transcript using a while loop. For every AUG start codon in the 5ʹ UTR, a flag was set for that reading frame and the length in codons was calculated until reaching one of the three stop codons.  Results for total uORF length and half-life for each mRNA were then outputted to two different tab delineated text files for the two different sources for half-lives. There was also an additional column in this text file that listed the transcript id used from the tab delimited text file, and the transcript id used from the database for the given data point, separated by a comma. 

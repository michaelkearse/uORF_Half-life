import sqlite3
from Bio import SeqIO
import os

# CHANGE THE FOLLOWING 2 LINES OF CODE:
DATADIRECTORY = '/....' # path to folder containing annotation and transcript file
nih_custom_transcript_db = '/..../hg38_nih_custom.db' # path to db to be created

GFFFILENAME = 'GRCh38_latest_genomic.gff' # annotation file name
RNAFILENAME = 'GRCh38_latest_rna.fna' # transcript file name

os.chdir(DATADIRECTORY)


conn = sqlite3.connect(nih_custom_transcript_db)
cur = conn.cursor()

cur.execute("""CREATE TABLE TRANSCRIPTS (
            TRANSCRIPT_ID text not null,
            CDS_START_INDEX integer,
            EXON_LENS integer,
            CDS_LEN integer,
            SEQ text,
            primary key(TRANSCRIPT_ID)
)""")


for seq_record in SeqIO.parse(RNAFILENAME, "fasta"):
    cur.execute("insert into TRANSCRIPTS(TRANSCRIPT_ID,SEQ,EXON_LENS) values('{}','{}',0)".format(seq_record.id,
                                                                                                str(seq_record.seq)))
print("done adding info from rna file")


typeidx = 2
infoidx = 8
start_idx = 3
end_idx = 4
strand_idx = 6

id_info_idx = 0
parent_rna_idx = 1

with open(GFFFILENAME, 'r') as f:
    for whole_line in f:
        if not whole_line.startswith('#'):
            line = whole_line.strip().split('\t')
            info = line[infoidx].split(';')

            # if current line in file is of type transcript
            if "RNA" in line[typeidx].upper() or 'transcript' in line[typeidx]:
                rna_id = info[id_info_idx].replace("ID=rna-", "")
                exons = []
                found_cds_start = False

            # if current line in file is of type exon
            if line[typeidx] == 'exon' and rna_id == info[parent_rna_idx].replace("Parent=rna-", ""):
                exons.append([int(line[start_idx]), int(line[end_idx])])
                cur.execute("UPDATE TRANSCRIPTS SET EXON_LENS = EXON_LENS + {} WHERE TRANSCRIPT_ID = '{}'".format(((abs(int(line[end_idx]) - int(line[start_idx]))) + 1), rna_id))

            # if current line in file is of type CDS
            elif line[typeidx] == 'CDS' and rna_id == info[parent_rna_idx].replace("Parent=rna-", ""):
                if not found_cds_start:
                    cds_start = 0
                    if line[strand_idx] == '+':
                        # calculate cds start position for the positive strand using the transcripts exons to determine the start of the 5' end of the transcript
                        for index in range(0, len(exons)):
                            if exons[index][0] <= int(line[start_idx]) <= exons[index][1]:
                                cds_start += abs(int(line[start_idx]) - exons[index][0])
                                break
                            else:
                                cds_start += abs(exons[index][1] - exons[index][0]) + 1
                    else:
                        # calculate cds start position for the minus strand using the transcripts exons to determine the start of the 5' end of the transcript
                        for index in range(0, len(exons)):
                            if exons[index][1] >= int(line[end_idx]) >= exons[index][0]:
                                cds_start += abs(exons[index][1] - int(line[end_idx]))
                                break
                            else:
                                cds_start += abs(exons[index][0] - exons[index][1]) + 1


                    cur.execute("UPDATE TRANSCRIPTS SET CDS_START_INDEX = {} WHERE TRANSCRIPT_ID = '{}'".format(cds_start, rna_id))

                    exons = []
                    found_cds_start = True


                cur.execute("UPDATE TRANSCRIPTS SET CDS_LEN = CDS_LEN + {} WHERE TRANSCRIPT_ID = '{}'".format(abs(int(line[end_idx]) - int(line[start_idx])) + 1, rna_id))


f.close()
print("done with annotation file")


conn.commit()
conn.close()

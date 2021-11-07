import os
import sqlite3

# CHANGE THE FOLLOWING 3 LINES OF CODE:
half_lives_file = '/..../Supplemental_Table_S3_half lifves BRIC seq.csv'  # Input the path to the csv data for transcript half-lives
nih_custom_transcript_db = '/.../hg38_nih_custom.db'  # Input the path to the db generated from other script
RESULTSDIR = '/...'  # Input path to folder where the results will be saved

START_ROW = 4  # row in csv that data entries start, Don't
BIN_SIZE = 5  # size of uorf bins

conn = sqlite3.connect(nih_custom_transcript_db)
cur = conn.cursor()


def parse_file(filename):
    # indexes for data in csv's
    id_idx = 1
    tani_idx = 12
    maekawa_idx = 17

    # lists for storing results
    uorf_len_bins_tani = []
    uorf_len_bins_maekawa = []
    tani_values = []
    maekawa_values = []

    no_results = 0
    count = 1

    with open(filename, 'r') as f:
        for line in f:
            if count >= START_ROW:
                line = line.split('\t')
                ids = line[id_idx].replace(' ', '').split(',')  # get all transcript ids for the given line and convert them into a string array

                max_len = 0
                no_results = 0
                for index in range(0, len(ids)):
                    cur.execute("SELECT SEQ, CDS_START_INDEX,EXON_LENS, TRANSCRIPT_ID From TRANSCRIPTS WHERE TRANSCRIPT_ID  LIKE ?", ['%' + str(ids[index]) + '%'])
                    transcript_results = cur.fetchall()
                    # if there are multiple transcripts that come back from the sql query
                    if len(transcript_results) > 1:
                        for result in transcript_results:
                            # ensure the given result is an isophorm of the transcript we are looking for and that it has data for cds start position
                            if result[3].split('.')[0] == ids[index] and result[1] is not None:
                                seq_len = len(result[0])
                                # choose the isophorm with the longest length and make sure that exon length is equal to transcript length
                                if seq_len > max_len and seq_len == result[2]:
                                    max_len = seq_len
                                    max_seq = result[0]
                                    max_cds_start_pos = result[1]
                                    max_id = ids[index]
                                    max_trans_id = result[3]
                    elif len(transcript_results) == 1:
                        # make sure the result from the sql query has data for cds start postion
                        if transcript_results[0][1] is not None:
                            seq_len = len(transcript_results[0][0])
                            # choose the isophorm with the longest length and make sure that exon length is equal to transcript length
                            if seq_len > max_len and seq_len == transcript_results[0][2]:
                                max_len = seq_len
                                max_seq = transcript_results[0][0]
                                max_cds_start_pos = transcript_results[0][1]
                                max_id = ids[index]
                                max_trans_id = transcript_results[0][3]

                if max_len > 0:
                    # calulate uorf length
                    uorf_len = calculate_uorfs(max_seq, max_cds_start_pos)

                    # calulate bin
                    if uorf_len == 0:
                        uorf_bin = 0
                    else:
                        uorf_bin = int((uorf_len - 1) / BIN_SIZE) + 1

                    # get half life and save the results to an array
                    if line[tani_idx] != 'NA':
                        tani_half_life = float(line[tani_idx])
                        uorf_len_bins_tani = add_to_bin(uorf_len_bins_tani, uorf_bin, tani_half_life)
                        tani_values.append([uorf_len, tani_half_life, max_id + ', ' + max_trans_id])

                    if line[maekawa_idx] != 'NA':
                        maekawa_half_life = float(line[maekawa_idx])
                        uorf_len_bins_maekawa = add_to_bin(uorf_len_bins_maekawa, uorf_bin, maekawa_half_life)
                        maekawa_values.append([uorf_len, maekawa_half_life, max_id + ', ' + max_trans_id])
                else:
                    no_results += 1

                if count % 1000 == 0:
                    print(str(count) + ' rows have been processed so far...')

            count += 1

    write_bins_to_output_file(uorf_len_bins_tani, 'Tani Bins.csv', BIN_SIZE)
    write_bins_to_output_file(uorf_len_bins_maekawa, 'Maekawa Bins.csv', BIN_SIZE)

    write_values_to_output_file(tani_values, 'Tani Scatter Plot.csv')
    write_values_to_output_file(maekawa_values, 'Maekawa Scatter Plot.csv')

    print('Could not determine results for ' + str(no_results) + ' transcripts')
    print("done")


# add item to list bin at given index
def add_to_bin(bin, index, item):
    while index >= len(bin):
        bin.append([])
    bin[index].append(item)

    return bin


# calculate uorf length given a sequence and a cds start position
def calculate_uorfs(seq, cds_start_pos):
    is_uorf = [False, False, False] # each index represents one of the 3 diffrent reading frames
    uorf_codons = 0

    in_cds = False
    index = 0
    # while the current postion is not in the cds region or one of the current reading frames is still in the middle of a uorf
    while index < len(seq) - 2 and (not in_cds or is_uorf[0] or is_uorf[1] or is_uorf[2]):
        frame = index % 3 # calculate current reading frame

        # check to see if the current codon in the current reading is start codon and that the current codon is in the cds
        if seq[index] == 'A' and seq[index + 1] == 'T' and seq[index + 2] == 'G' and not in_cds:
            is_uorf[frame] = True
            uorf_codons += 1
        elif is_uorf[frame]:
            # the codon in the current fram is part of a uorf
            uorf_codons += 1
            # check to see if the current codon is a stop codon
            if ((seq[index] == 'T' and seq[index + 1] == 'G' and seq[index + 2] == 'A') or (seq[index] == 'T' and seq[index + 1] == 'A' and seq[index + 2] == 'G')
                    or (seq[index] == 'T' and seq[index + 1] == 'A' and seq[index + 2] == 'A')):
                is_uorf[frame] = False

        index += 1

        #check if the next codon will be in the cds
        if index + 2 >= cds_start_pos:
            in_cds = True

    return uorf_codons

# write results from bins to a tab delineated file
def write_bins_to_output_file(bins, filename, bin_size):
    os.chdir(RESULTSDIR)

    output_file = open(filename, "w")
    output_file.write('uORF Len\tHalf Lives\tNumber\n')

    for index in range(0, len(bins)):
        if index == 0:
            line = '0\t'
        else:
            line = str((index - 1) * bin_size + 1) + '-' + str((index) * bin_size) + '\t'
        if len(bins[index]) == 0:
            line = line + '0\t'
        else:
            line = line + str(sum(bins[index]) / len(bins[index])) + '\t'
        line = line + str(len(bins[index])) + '\n'
        output_file.write(line)

# write results from values to a tab delineated file.
def write_values_to_output_file(values, filename):
    os.chdir(RESULTSDIR)

    output_file = open(filename, "w")
    output_file.write('uORF Len\tHalf Lives\tTranscript Id\n')

    for pair in values:
        line = str(pair[0]) + '\t' + str(pair[1]) + '\t' + str(pair[2]) + '\n'
        output_file.write(line)

parse_file(half_lives_file)
conn.close()

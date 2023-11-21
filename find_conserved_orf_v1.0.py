########################################################################################################################
########################################################################################################################
#Ruither Arthur Loch Gomes
#ruitherarthur@gmail.com
#01/2023
########################################################################################################################
#Using a group of similar genomes, the algorithm searches for open reading frames that are conserved among them, using
# muscle and pairwise alignments.
########################################################################################################################
#lIRABRIES
########################################################################################################################
########################################################################################################################
import re
import subprocess
import time
import os
import concurrent.futures
from datetime import datetime
import numpy as np
import platform
import argparse


########################################################################################################################
########################################################################################################################
#AUXILIAR AND MODYFYNG DATA FUNCTIONS
########################################################################################################################
########################################################################################################################
# Function to format time and print
def format_time(time_in_sec):
    minutes = int(time_in_sec // 60)
    seconds = int(time_in_sec % 60)
    if minutes == 0:
        return f"{seconds} seconds"
    elif seconds == 0:
        return f"{minutes} minutes"
    else:
        return f"{minutes} minutes and {seconds} seconds"
########################################################################################################################
########################################################################################################################
# write to an output file
def write_output(path, output):
    with open(path, "a+") as outfile:
        outfile.write(output)

########################################################################################################################
########################################################################################################################
# Function to open a FASTA file and read its contents
def open_fasta(path):
    fasta_dic = {}
    fasta_seq = ""
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if len(fasta_dic) == 0:
                    seq_id = line[1:].replace("\n", "")
                    fasta_dic[seq_id] = ""
                    continue
                else:
                    fasta_dic[seq_id] = fasta_seq
                    fasta_seq = ""
                seq_id = line[1:].replace("\n", "")
                fasta_dic[seq_id] = fasta_seq
            else:
                fasta_seq += line.replace("\n", "")
        else:
            fasta_dic[seq_id] = fasta_seq
    return fasta_dic


########################################################################################################################
########################################################################################################################
# Function to generate reverse complement of a DNA sequence
def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp_seq = ''.join([complement[base] for base in dna_sequence[::-1]])
    return reverse_comp_seq


########################################################################################################################
########################################################################################################################
# Function to translate DNA sequence to protein sequence
def translate_dna_to_protein(dna_sequence):
    global genetic_code
    genetic_code = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    protein_sequence = ""
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i + 3]
        if codon in genetic_code:
            protein_sequence += genetic_code[codon]
        else:
            protein_sequence += "X"

    return protein_sequence


########################################################################################################################
########################################################################################################################
# Function to replace forbidden characters in a identifier
def change_chars(input_text):
    searched_chars = r'[\\/:"*?<>|,\n\t^.;|@_# ]'
    # new_text = re.sub(searched_chars, '-', input_text)
    return re.sub(searched_chars, '-', input_text)


########################################################################################################################
########################################################################################################################
# Function to find open reading frames with translations
def find_orfs_with_trans(seq):
    if not CHECK_ALT_CODON:
        start_codons = ['ATG']
        start_aas = ['M']
    else:
        start_codons = ['ATG', 'CTG', 'GTG', 'TTG', 'ATT', 'ATC']
        start_aas = ['M', 'L', 'V', 'I']

    protein_orfs_resp = []
    seq = seq.replace("\n", "")
    seq_len = len(seq)
    for strand, cur_nt_seq in [(+1, seq), (-1, reverse_complement(seq))]:
        for reading_frame in range(3):
            trans_seq = translate_dna_to_protein(cur_nt_seq[reading_frame:])
            trans_seq_len = len(trans_seq)
            aa_orf_start_pos = 0
            while aa_orf_start_pos < trans_seq_len:
                aa_orf_end_pos = trans_seq.find("*", aa_orf_start_pos)
                if aa_orf_end_pos == -1:
                    aa_orf_end_pos = trans_seq_len

                aa_seq = trans_seq[aa_orf_start_pos:aa_orf_end_pos]

                nt_start_codon_pos = -1
                aaseq_start_pos = -1
                for aa_let_pos, letter in enumerate(aa_seq):
                    if letter in start_aas:
                        nt_sc_pos = reading_frame + (aa_orf_start_pos + aa_let_pos) * 3

                        if cur_nt_seq[nt_sc_pos:nt_sc_pos + 3] in start_codons:
                            nt_start_codon_pos = nt_sc_pos
                            aaseq_start_pos = aa_let_pos

                            break  # only search until the first start codon

                # save when it has the minimum size
                if nt_start_codon_pos != -1 and (len(aa_seq) - aaseq_start_pos) >= MIN_PEP_SIZE:
                    if strand == 1:
                        sense_nt_start_pos = nt_start_codon_pos
                        sense_nt_end_pos = min(seq_len, reading_frame + aa_orf_end_pos * 3 + 3)

                        antisense_nt_end_pos = seq_len - sense_nt_start_pos
                        antisense_nt_start_pos = seq_len - sense_nt_end_pos

                        sense_nt_seq = cur_nt_seq[sense_nt_start_pos:sense_nt_end_pos]
                        antisense_nt_seq = reverse_complement(cur_nt_seq)[antisense_nt_start_pos:antisense_nt_end_pos]

                        new_aa_seq = aa_seq[aaseq_start_pos:]

                    else:
                        antisense_nt_start_pos = nt_start_codon_pos
                        antisense_nt_end_pos = min(seq_len, reading_frame + aa_orf_end_pos * 3 + 3)

                        sense_nt_end_pos = seq_len - antisense_nt_start_pos
                        sense_nt_start_pos = seq_len - antisense_nt_end_pos

                        sense_nt_seq = seq[sense_nt_start_pos:sense_nt_end_pos]
                        antisense_nt_seq = cur_nt_seq[antisense_nt_start_pos:antisense_nt_end_pos]

                        new_aa_seq = aa_seq[aaseq_start_pos:]

                    protein_orfs_resp.append((strand, reading_frame,
                                              (sense_nt_start_pos, sense_nt_end_pos, sense_nt_seq),
                                              (antisense_nt_start_pos, antisense_nt_end_pos, antisense_nt_seq),
                                              new_aa_seq))

                aa_orf_start_pos = aa_orf_end_pos + 1

    return protein_orfs_resp


########################################################################################################################
########################################################################################################################
# Function to check sequence length
def check_len(tested_seq_len, fixed_len):
    lower_bound = fixed_len * ((1 - MAX_PROT_LEN_PERC))
    upper_bound = fixed_len * ((1 + MAX_PROT_LEN_PERC))
    return lower_bound <= tested_seq_len <= upper_bound


########################################################################################################################
########################################################################################################################
#LOADING AA SEQUENCES
########################################################################################################################
########################################################################################################################
# Parallel processing (PP) sequence parsing function to increase the analysis spped
def do_parallel_ident_parsing(record_batch_list):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        result_info = executor.map(parse_identifier, record_batch_list)
    for result in result_info:
        ORFS_LIST.extend(result)


########################################################################################################################
########################################################################################################################
# Function called in the PP, where analyses each sequence getting its ORFs and saving the information
def parse_identifier(cur_identifier):
    cur_sequence = INPUT_GENOMES_DIC[cur_identifier]
    orfslist = find_orfs_with_trans(cur_sequence)
    aux_orfs_list = []
    orf_counter = 0
    for strand, frame, sense_info, antisense_info, aa_seq in orfslist:
        orf_counter+=1
        aux_orfs_list.append((cur_identifier+"__"+str(orf_counter), strand, frame, sense_info, antisense_info, aa_seq))
    return aux_orfs_list


########################################################################################################################
########################################################################################################################
# Parallel processing all ORFs gathered before to identify conservation among them for the given sequences
def do_parallel_seq_parsing(tested_batch_set):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        result_info = executor.map(parse_sequences, tested_batch_set)

    for results_index, result in enumerate(result_info):
        if result[0]:
            if result[1][0][-1] in ORFS_DIC.keys():
                ORFS_DIC[result[1][0][-1]].append(result[1][1])
                ORF_IN_ID_LIST.append(result[1][1][0])
            else:
                ORFS_DIC[result[1][0][-1]] = [result[1][0]]
                ORFS_DIC[result[1][0][-1]].append(result[1][1])
                ORF_IN_ID_LIST.append(result[1][0][0])
                ORF_IN_ID_LIST.append(result[1][1][0])

########################################################################################################################
########################################################################################################################
#Function called in the PP, compare ORFs by pairs using muscle and alignment without gaps
def parse_sequences(tested_set):

    name1 = tested_set[1][0]
    name2 = tested_set[0][0]
    name11 = name1.split("-")[0]+"--"+name1.split("__")[-1]
    name22 = name2.split("-")[0]+"--"+name2.split("__")[-1]

    write_output(OUTPUT_ORFS_PATH + name11 + "--" + name22 + ".temp", ">" + tested_set[0][0] + "\n" + tested_set[0][-1] + "\n")
    write_output(OUTPUT_ORFS_PATH + name11 + "--" + name22 + ".temp", ">" + tested_set[1][0] + "\n" + tested_set[1][-1] + "\n")

    try:
        subprocess.run([MUSCLE_PATH, '-align', OUTPUT_ORFS_PATH + name11 + "--" + name22 + ".temp", '-output',
                        OUTPUT_ORFS_PATH + name11 + "--" + name22 + ".aln", "--quiet"])

    except FileNotFoundError:
        print(f"Muscle not find at: {MUSCLE_PATH}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred at Muscle: {e}")

    align_result_dic = open_fasta(OUTPUT_ORFS_PATH + name11 + "--" + name22 + ".aln")

    os.remove(OUTPUT_ORFS_PATH + name11 + "--" + name22 + ".temp")
    os.remove(OUTPUT_ORFS_PATH + name11 + "--" + name22 + ".aln")

    hit_counter = 0
    gap_counter = 0

    seq1,seq2 = list(align_result_dic.values())[0],list(align_result_dic.values())[1]
    for idx,letter in enumerate(seq1):
        if seq1[idx]==seq2[idx]:
            hit_counter+=1
        if letter=="-":
            gap_counter+=1

    if not GAPCHECK:
        identity = hit_counter/(len(seq1)-gap_counter)
    else:
        identity = hit_counter/len(seq1)

    if identity > IDENTITY_THRESHOLD:
        return (True,tested_set)
    else:
        return (False,())


########################################################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Main
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
########################################################################################################################
DEFAULT_PATH = 'examples/EYMA_12.fasta'

cpus_avail = os.cpu_count()

parser = argparse.ArgumentParser(description='This tool receives a fasta file containing multiple genomes, find all '
                                             'open read frames for both standard or alternative start codons, and '
                                             'group out the obtained ORFs, using the identity with or without gaps '
                                             'between them for that, after a muscle alignment, and realising this groups'
                                             ' as all ORFs found for each input sequence.')
parser.add_argument('-input', default=DEFAULT_PATH, type=str, help='The path to the input file.'
                                                                       'Default: the example file of 12 "Geminiviridae"'
                                                                       ' alphasatellites.')
parser.add_argument('-out', default='output/', type=str, help='A output folder name, that will be created '
                                                               'if dont exists. Default: a folder named "output" will be'
                                                               ' created for the results.')
parser.add_argument('-muscle_path', default="Standard", type=str, help='A path to an alternative muscle '
                                                                        'version. Default: 5.1 at the bin folder.')
parser.add_argument('-pep_min_size', default=6, type=int, help='The number of amino acids a peptide must '
                                                                'have to be considered valid. Default: 6 aa.')
parser.add_argument('-alt_codon', default=False, type=bool, help='Use alternative codons when searching '
                                                                  'for the start codon (CTG, GTG, TTG, ATT, ATC). '
                                                                  'Default: False, using ATG only, change to True,'
                                                                 ' but it may take longer to run.')
parser.add_argument('-use_gaps', action='store_true', default=False, help='Include to use gaps for the '
                                                                           'identity between '
                                                            'proteins. Default: "False", not considering gaps.')
parser.add_argument('-prot_len_perc', default=0.6, type=float, help="The percentage (in decimals) of the "
                                                                     "ORF's length that will be accepted to vary (for "
                                                                     "more or less than the representative) between "
                                                                     "the first found sequence and the others. Default:"
                                                                     " 0.6, which means the a sequence must have more "
                                                                     "or less then 60 percent of the length of the first found "
                                                                     "sequence to be compared to it.")
parser.add_argument('-identity_threshold', default=0.7, type=float, help='The percentage of similarity, '
                                                                          'identity, between to ORFs to consider to be '
                                                                          'in the same group. Default: 0.7, which is 70'
                                                                          ' percent of similarity for grouping.')
parser.add_argument('-nproc', default=cpus_avail, type=int, help='The number of processing units that will'
                                                                  ' be used for parallel processing. Default: the number'
                                                                  ' of available processors.')

parser.add_argument('-comp_seq', default=False, type=bool, help='Define if the complementary sequence, both'
                                                                 ' for sense and anti-sense frames, will be shown together'
                                                                 'in the complete output. Default: False, use true to show.')

args = parser.parse_args()

PATH_TO_INPUT_FILE = args.input
OUTPUT_PATH = args.out
NEW_MUSCLE_PATH = args.muscle_path
MIN_PEP_SIZE = args.pep_min_size
CHECK_ALT_CODON = args.alt_codon
GAPCHECK = args.use_gaps
MAX_PROT_LEN_PERC = args.prot_len_perc
IDENTITY_THRESHOLD = args.identity_threshold
BATCH_SIZE = args.nproc
COMP_SEQ = args.comp_seq

os_dir_sep = os.path.sep
if not OUTPUT_PATH.endswith(os_dir_sep):
    OUTPUT_PATH+=os_dir_sep
OUTPUT_ORFS_PATH = OUTPUT_PATH + "ORFs/"
OUTPUT_NAME = PATH_TO_INPUT_FILE.split(".f")[0].split("/")[-1]

#Get the correct muscle path for the user OS
MUSCLE_PATH = ""
if NEW_MUSCLE_PATH == "Standard":
    #get the OS running the system
    opSys = platform.system()[:3].lower()
    #change for macOS specs
    if opSys == "dar":
        running_arch = platform.processor()
        if 'arm' in running_arch.lower():
            opSys = 'macos_arm64'
        else:
            opSys = 'macos_intel64'
    #get the path
    for file in os.listdir("bin/muscle"):
        if opSys in file:
            MUSCLE_PATH = "bin/muscle/"+file
else:
    MUSCLE_PATH = NEW_MUSCLE_PATH
    print("Changed muscle path!")

print("Muscle path: ", MUSCLE_PATH)
if not GAPCHECK:
    print("Not considering gaps to the identity between proteins.")
else:
    print("Considering gaps to the identity between proteins.")
########################################################################################################################
#creating output directory and remove file in case of previous crash

os.makedirs(OUTPUT_PATH, exist_ok=True)
os.makedirs(OUTPUT_ORFS_PATH, exist_ok=True)
for file1 in os.listdir(OUTPUT_ORFS_PATH):
    if ".aln" in file1 or ".temp" in file1:
        os.remove(OUTPUT_ORFS_PATH + file1)
        continue
########################################################################################################################
########################################################################################################################
#Starting
print("Starting at: "+(datetime.now()).strftime("%d/%m/%Y %H:%M")+".")
starttime = time.time()

INIT_INPUT_GENOMES_DIC = open_fasta(PATH_TO_INPUT_FILE)
#Change the identifier to be used as a file name in the OS, removing forbidden characters and check for unacceptable
# characters, changing them to N.
INPUT_GENOMES_DIC = {}
for k,v in INIT_INPUT_GENOMES_DIC.items():
    k1 = change_chars(k)
    if not all(base in 'ACTG' for base in v):
        print("Genome with non-canonical base, removed: ",k)
        continue
    INPUT_GENOMES_DIC[k1] = v
identifiers_total = len(INPUT_GENOMES_DIC)


########################################################################################################################
########################################################################################################################
#Loading data
print("First loading the data information.")
print("Total of sequences in file: " + str(identifiers_total) + ".")
ORFS_LIST = []
total_ids = 0
len_seq = 0
batch_counter = 0
batch_id_list = []
id_list = list(INPUT_GENOMES_DIC.keys())
end_bool = False
print("\rLoading data...",end="")
for idx_ident, ident1 in enumerate(id_list):
    batch_counter+=1
    batch_id_list.append(ident1)
    total_ids += 1
    #Remove output file if exists
    if os.path.exists(OUTPUT_ORFS_PATH +OUTPUT_NAME+ "_" + ident1 + "_nt.fasta"):
        os.remove(OUTPUT_ORFS_PATH +OUTPUT_NAME+ "_" + ident1 + "_nt.fasta")
    if os.path.exists(OUTPUT_ORFS_PATH +OUTPUT_NAME+ "_" + ident1 + "_AA.fasta"):
        os.remove(OUTPUT_ORFS_PATH +OUTPUT_NAME+ "_" + ident1 + "_AA.fasta")

    if id_list[idx_ident] == id_list[-1]:
        end_bool = True

    if batch_counter == BATCH_SIZE or end_bool:
        do_parallel_ident_parsing(batch_id_list)
        batch_id_list = []
        batch_counter = 0
        print('\rLoading data... %d%%' % (total_ids * (100 / len(id_list))), end="")

########################################################################################################################
########################################################################################################################
#Parsing loaded data
print("\n\rParsing data...",end="")
ORFS_DIC = {}
ORF_IN_ID_LIST = []

batch_counter = 0
batch_id_list = []
end_bool = False
parse_time = time.time()
for info_idx, seq_info in enumerate(ORFS_LIST):
    if seq_info[0] in ORF_IN_ID_LIST:
        continue
    for info_idx2, seq_info2 in enumerate(ORFS_LIST):
        id_info_list1 = seq_info2[0].split("__")[0]

        if info_idx2 <= info_idx or seq_info2[0] in ORF_IN_ID_LIST:

            continue

        if not check_len(len(seq_info[-1]),len(seq_info2[-1])):
            continue

        batch_counter += 1
        batch_id_list.append((seq_info,seq_info2))

        if ORFS_LIST[info_idx] == ORFS_LIST[-1]:
            end_bool = True

        if batch_counter == BATCH_SIZE or end_bool:
            do_parallel_seq_parsing(batch_id_list)
            batch_id_list = []
            batch_counter = 0
    print('\rParsing data... %d%% (Time elapsed: %.2f secs)' % (len(ORF_IN_ID_LIST) * (100 / len(ORFS_LIST)),time.time()-parse_time), end="")
print('\rParsing data... %d%% (Time elapsed: %.2f secs)' % (len(ORFS_LIST) * (100 / len(ORFS_LIST)),time.time()-parse_time), end="")

####################################################################################################################
####################################################################################################################
#saving results

print("\nTotal of possible unique sequences after parsing: "+str(len(ORFS_DIC))+".")

full_outfile = open(OUTPUT_PATH + OUTPUT_NAME + "_All_info.txt", "w")
only_conserved_orfs_outfile = open(OUTPUT_PATH + OUTPUT_NAME + "_ConservedOnly.txt", "w")

full_outfile.write("Starting at: "+(datetime.now()).strftime("%d/%m/%Y %H:%M")+".")
full_outfile.write("*Identifiers information: Length, Strand, Frame and Nucleotide (start and end) positions. *\n")
consev_counter = 0

orfs_by_len_list = sorted((list(ORFS_DIC.keys())), key=len, reverse=True)

for unique_aa_seq in orfs_by_len_list:
    orf_info_list = ORFS_DIC[unique_aa_seq]

    conserved_bool = False
    if identifiers_total <= len(orf_info_list):
        consev_counter+=1
        full_outfile.write("\n### Conserved ORF among the genomes. ###\n")
        if identifiers_total < len(orf_info_list):
            full_outfile.write("\n*  The ORF is  or has an early stopping in at least one genome!  *\n")
        conserved_bool = True

    len_values_list = [len(orf_set[-1]) for orf_set in orf_info_list]

    median_len = round(np.median(len_values_list),1)
    mean_len = round(np.mean(len_values_list),1)
    max_len = max(len_values_list)
    min_len = min(len_values_list)

    full_outfile.write("\n*** ORF length metrics ***\n")
    full_outfile.write("Mean length: " + str(mean_len) + "\n")
    full_outfile.write("Median length: " + str(median_len) + "\n")
    full_outfile.write("Max length: " + str(max_len) + "\n")
    full_outfile.write("Min length: " + str(min_len) + "\n")
    full_outfile.write("Total found: " + str(len(len_values_list)) + "\n")
    full_outfile.write("***  List of sequences ***" + "\n")

    if conserved_bool:
        only_conserved_orfs_outfile.write("\n*** ORF length metrics ***\n")
        only_conserved_orfs_outfile.write("Mean length: " + str(mean_len) + "\n")
        only_conserved_orfs_outfile.write("Median length: " + str(median_len) + "\n")
        only_conserved_orfs_outfile.write("Max length: " + str(max_len) + "\n")
        only_conserved_orfs_outfile.write("Min length: " + str(min_len) + "\n")
        only_conserved_orfs_outfile.write("Total found: " + str(len(len_values_list)) + "\n")
        only_conserved_orfs_outfile.write("***  List of sequences ***" + "\n")

    in_ids_aux_list = []
    for cur_genome_orf_info in orf_info_list:

        #Get the identifier and theorf position
        id_info_list = cur_genome_orf_info[0].split("__")
        first_id1, orf_pos = id_info_list[0], id_info_list[1]

        if cur_genome_orf_info[1] == 1:
            full_info_id = (first_id1 + "__L" + str(len(cur_genome_orf_info[-1])) +
                        "_S" + str(cur_genome_orf_info[1]) + "_F" + str(cur_genome_orf_info[2]) +
                        "_NtPos" + str(cur_genome_orf_info[3][0]) +"-" + str(cur_genome_orf_info[3][1]))
            sense_seq = cur_genome_orf_info[3][2]
            comp_seq = cur_genome_orf_info[4][2]
        else:
            full_info_id = (first_id1 + "__L" + str(len(cur_genome_orf_info[-1])) +
                        "_S" + str(cur_genome_orf_info[1]) + "_F" + str(cur_genome_orf_info[2]) +
                        "_NtPos" + str(cur_genome_orf_info[4][0]) + "-" + str(cur_genome_orf_info[4][1]))
            sense_seq = cur_genome_orf_info[4][2]
            comp_seq = cur_genome_orf_info[3][2]

            #Add orfs to fasta file using the identifier as the name
        write_output(OUTPUT_ORFS_PATH +OUTPUT_NAME+ "_" + first_id1 + "_AA.fasta",
                     ">" + full_info_id +"_ORFpos"+orf_pos+"\n" + cur_genome_orf_info[-1] +"\n")

        write_output(OUTPUT_ORFS_PATH +OUTPUT_NAME+ "_" + first_id1 + "_nt.fasta",
                     ">" + full_info_id + "\n" + sense_seq + "\n")
        if COMP_SEQ:
            write_output(OUTPUT_ORFS_PATH +OUTPUT_NAME+ "_" + first_id1 + "_nt.fasta",
                         ">" + full_info_id + "_Comp\n" + comp_seq + "\n")

        in_ids_aux_list.append(first_id1)
        full_outfile.write(">" + full_info_id + "\n"+ cur_genome_orf_info[-1] + "\n")

        if conserved_bool:
            only_conserved_orfs_outfile.write(">" + full_info_id + "\n"+str(cur_genome_orf_info[-1]) + "\n")

    full_outfile.write("***\n")
    #if there is duplicated orfs
    orf_set = set(in_ids_aux_list)
    if len(orf_set) > len(id_list):
        full_outfile.write("Genomes with multiple instances of the ORF:\n")
        for i in in_ids_aux_list:
            count = in_ids_aux_list.count(i)
            if count > 1:
                full_outfile.write(i+"\n")
                break
    elif len(orf_set) < len(id_list):
        #if there is some genome dont have the ORF
        full_outfile.write("Genomes without the ORF:\n")
        for i in id_list:
            if i not in in_ids_aux_list:
                full_outfile.write(i + "\n")
    full_outfile.write("--------------------------------\n")

full_outfile.write("*** Used parameters ***\n")
full_outfile.write(f'PATH_TO_INPUT_FILE: {PATH_TO_INPUT_FILE}\n')
full_outfile.write(f'OUTPUT_PATH: {OUTPUT_PATH}\n')
full_outfile.write(f'NEW_MUSCLE_PATH: {NEW_MUSCLE_PATH}\n')
full_outfile.write(f'MIN_PEP_SIZE: {MIN_PEP_SIZE}\n')
full_outfile.write(f'CHECK_ALT_CODON: {CHECK_ALT_CODON}\n')
full_outfile.write(f'MAX_PROT_LEN_PERC: {MAX_PROT_LEN_PERC}\n')
full_outfile.write(f'IDENTITY_THRESHOLD: {IDENTITY_THRESHOLD}\n')
full_outfile.write(f'BATCH_SIZE: {BATCH_SIZE}\n')
full_outfile.write(f'COMP_SEQ: {COMP_SEQ}\n')
full_outfile.write("***")
full_outfile.write("Finished at: "+(datetime.now()).strftime("%d/%m/%Y %H:%M")+".")

time_elapsed = format_time((time.time()-starttime))
full_outfile.write("Total time taken : "+time_elapsed+".")
print("Total Time taken: "+ time_elapsed +".")
only_conserved_orfs_outfile.write("Total conserved orfs: " + str(consev_counter) + "\n--------------------------------")
full_outfile.close()
only_conserved_orfs_outfile.close()
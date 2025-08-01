########################################################################################################################
#Ruither Arthur Loch Gomes
#ruitherarthur@gmail.com
#07/2025
########################################################################################################################
#Using a group of similar genomes, the algorithm searches for open reading frames that are conserved among them, using
# muscle and pairwise alignments.
########################################################################################################################
#lIRABRIES
########################################################################################################################
import time
import os

from datetime import datetime
import platform
import argparse

from src.bio_utils import open_fasta, format_sequences
from src.parse_utils import get_orf_list, compare_all_orfs, save_results
from src.utils import format_time, create_dirs

########################################################################################################################

cpus_avail = os.cpu_count()

parser = argparse.ArgumentParser(description='This tool receives a fasta file containing multiple genomes, find all '
                                             'open read frames for both standard or alternative start codons, and '
                                             'group out the obtained ORFs, using the identity with or without gaps '
                                             'between them for that, after a muscle alignment, and realising this groups'
                                             ' as all ORFs found for each input sequence.')
parser.add_argument('--input_file', default='examples/EYMA_12.fasta', type=str, help='The path to the input file.'
                                                                       'Default: the example file of 12 "Geminiviridae"'
                                                                       ' alphasatellites.')
parser.add_argument('--out', default='output/', type=str, help='A output folder name, that will be created '
                                                               'if dont exists. Default: a folder named "output" will be'
                                                               ' created for the results.')
parser.add_argument('--muscle_path', default="Standard", type=str, help='A path to an alternative muscle '
                                                                        'version. Default: 5.1 at the bin folder.')
parser.add_argument('--pep_min_size', default=6, type=int, help='The number of amino acids a peptide must '
                                                                'have to be considered valid. Default: 6 aa.')
parser.add_argument('--alt_codon', default=False, type=bool, help='Use alternative codons when searching '
                                                                  'for the start codon (CTG, GTG, TTG, ATT, ATC). '
                                                                  'Default: False, using ATG only, change to True,'
                                                                 ' but it may take longer to run.')
parser.add_argument('--use_gaps', action='store_true', default=False, help='Include to use gaps for the '
                                                                           'identity between '
                                                            'proteins. Default: "False", not considering gaps.')
parser.add_argument('--prot_len_perc', default=0.6, type=float, help="The percentage (in decimals) of the "
                                                                     "ORF's length that will be accepted to vary (for "
                                                                     "more or less than the representative) between "
                                                                     "the first found sequence and the others. Default:"
                                                                     " 0.6, which means the a sequence must have more "
                                                                     "or less then 60 percent of the length of the first found "
                                                                     "sequence to be compared to it.")
parser.add_argument('--identity_threshold', default=0.7, type=float, help='The percentage of similarity, '
                                                                          'identity, between to ORFs to consider to be '
                                                                          'in the same group. Default: 0.7, which is 70'
                                                                          ' percent of similarity for grouping.')
parser.add_argument('--nproc', default=cpus_avail, type=int, help='The number of processing units that will'
                                                                  ' be used for parallel processing. Default: the number'
                                                                  ' of available processors.')

parser.add_argument('--comp_seq', default=False, type=bool, help='Define if the complementary sequence, both'
                                                                 ' for sense and anti-sense frames, will be shown together'
                                                                 'in the complete output. Default: False, use true to show.')

args = parser.parse_args()

path_to_input_file = args.input_file
output_path = args.out
new_muscle_path = args.muscle_path
min_pep_size = args.pep_min_size
check_alt_codon = args.alt_codon
gapcheck = args.use_gaps
max_prot_len_perc = args.prot_len_perc
identity_threshold = args.identity_threshold
batch_size = args.nproc
comp_seq_bool = args.comp_seq

os_dir_sep = os.path.sep
if not output_path.endswith(os_dir_sep):
    output_path+=os_dir_sep
output_orfs_path = output_path + "ORFs/"
output_name = path_to_input_file.split(".f")[0].split("/")[-1]

#Get the correct muscle path for the user OS
muscle_path = ""
if new_muscle_path == "Standard":
    #get the OS running the system
    opSys = platform.system()[:3].lower()

    if opSys == "dar":
        running_arch = platform.processor()
        if 'arm' in running_arch.lower():
            opSys = 'macos_arm64'
        else:
            opSys = 'macos_intel64'

    for file in os.listdir("bin/muscle"):
        if opSys in file:
            muscle_path = "bin/muscle/" + file
else:
    muscle_path = new_muscle_path

#creating output directory and remove file in case of previous crash
create_dirs(output_path,output_orfs_path)

print("Starting at: "+(datetime.now()).strftime("%d/%m/%Y %H:%M")+".")
starttime = time.time()
init_input_genomes_dic = open_fasta(path_to_input_file)

#remove non canonical base characters
input_genomes_dic = format_sequences(init_input_genomes_dic)
identifiers_total = len(input_genomes_dic)

#Load data
orfs_set,id_list = get_orf_list(
    input_genomes_dic,
    batch_size,
    check_alt_codon,
    min_pep_size
)

#Parse loaded data
orfs_dic = compare_all_orfs(
    list(orfs_set),
    batch_size,
    output_orfs_path,
    muscle_path,
    gapcheck,
    identity_threshold,
    max_prot_len_perc
)

#save results
save_results(
    orfs_dic,
    output_path,
    output_name,
    identifiers_total,
    output_orfs_path,
    comp_seq_bool,
    id_list
)

#write run parameters
full_outfile = open(output_path + output_name + "_All_info.txt", "a")
full_outfile.write("*** Run parameters ***\n")
full_outfile.write(f'PATH_TO_INPUT_FILE: {path_to_input_file}\n')
full_outfile.write(f'OUTPUT_PATH: {output_path}\n')
full_outfile.write(f'NEW_MUSCLE_PATH: {new_muscle_path}\n')
full_outfile.write(f'MIN_PEP_SIZE: {min_pep_size}\n')
full_outfile.write(f'CHECK_ALT_CODON: {check_alt_codon}\n')
full_outfile.write(f'MAX_PROT_LEN_PERC: {max_prot_len_perc}\n')
full_outfile.write(f'IDENTITY_THRESHOLD: {identity_threshold}\n')
full_outfile.write(f'BATCH_SIZE: {batch_size}\n')
full_outfile.write(f'COMP_SEQ: {comp_seq_bool}\n')
full_outfile.write("***")
full_outfile.write("Finished at: " + (datetime.now()).strftime("%d/%m/%Y %H:%M") + ".")
full_outfile.write("Total time taken : "+format_time((time.time()-starttime))+".")
full_outfile.close()


print("Total Time taken: "+ format_time((time.time()-starttime)) +".")

import os
import time
import subprocess
import concurrent.futures
from datetime import datetime
import numpy as np
import tempfile
from src.bio_utils import find_orfs_with_trans,open_fasta
from src.utils import write_output, check_len, ensure_executable


def get_orf_list(input_genomes_dic,batch_size,check_alt_codon,min_pep_size):
    """
        Extracts ORFs (Open Reading Frames) from each genome in the input dictionary and aggregates them into a list.

        Parameters:
        - input_genomes_dic (dict): Dictionary of genome sequences, with genome identifiers as keys.
        - batch_size (int): Number of genomes to process in parallel batches.
        - check_alt_codon (bool): Whether to allow alternative start codons.
        - min_pep_size (int): Minimum length of peptides to consider.

        Returns:
        - orfs_set (set): List of extracted ORFs with metadata.
        - id_list (list): List of all genome identifiers processed.
        """
    orfs_set = set()
    batch_counter = 0
    batch_id_list = []
    id_list = list(input_genomes_dic.keys())

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()-1)
    for idx_ident, ident1 in enumerate(id_list):
        batch_counter += 1
        batch_id_list.append(ident1)

        if batch_counter == batch_size or idx_ident == len(id_list) - 1:
            resp_list = do_parallel_ident_parsing(
                batch_id_list,
                orfs_set,
                input_genomes_dic,
                check_alt_codon,
                min_pep_size,
                executor
            )
            orfs_set.update(resp_list)
            batch_id_list = []
            batch_counter = 0
    executor.shutdown()
    return orfs_set,id_list


# Parallel processing (PP) sequence parsing function to increase the analysis sped
def do_parallel_ident_parsing(
        record_batch_list,
        orfs_set,
        input_genomes_dic,
        check_alt_codon,
        min_pep_size,
        executor = None):
    """
        Processes a batch of genome identifiers in parallel to extract ORFs from each.

        Parameters:
        - record_batch_list (list): List of genome identifiers in the current batch.
        - orfs_set (set): Set to be updated with new ORFs found.
        - input_genomes_dic (dict): Dictionary containing genome sequences.
        - check_alt_codon (bool): Whether to allow alternative start codons.
        - min_pep_size (int): Minimum peptide size to consider.

        Returns:
        - orfs_set (set): Updated set containing extracted ORFs.
        """

    data1 = [(ident,input_genomes_dic,check_alt_codon,min_pep_size) for ident in record_batch_list]

    result_info = executor.map(parse_identifier, data1)
    for result in result_info:
        orfs_set.update(result)
    return orfs_set

# Function called in the PP, where analyses each sequence getting its ORFs and saving the information
def parse_identifier(cur_data):
    """
        Parses a single genome entry to identify ORFs and format them for downstream analysis.

        Parameters:
        - cur_data (tuple): Tuple containing (identifier, input_genomes_dic, check_alt_codon, min_pep_size)

        Returns:
        - aux_orfs_list (list): List of extracted ORFs for the given genome, each with metadata.
        """
    cur_identifier, input_genomes_dic, check_alt_codon, min_pep_size = cur_data
    cur_sequence = input_genomes_dic[cur_identifier]
    orfslist = find_orfs_with_trans(cur_sequence,check_alt_codon,min_pep_size)
    aux_orfs_set = set()
    orf_counter = 0
    for strand, frame, sense_info, antisense_info, aa_seq in orfslist:
        orf_counter+=1
        aux_orfs_set.add((cur_identifier+"_orf_"+str(orf_counter), strand, frame, sense_info, antisense_info, aa_seq))
    return aux_orfs_set

def compare_all_orfs(orfs_list,
                     batch_size,
                     output_orfs_path,
                     muscle_path,
                     gapcheck,
                     identity_threshold,
                     max_prot_len_perc
                     ):
    """
       Compares all ORFs pairwise to identify conserved ORFs across genomes using MUSCLE alignment.

       Parameters:
       - orfs_list (list): List of ORF metadata from multiple genomes.
       - batch_size (int): Number of ORF pairs to process per batch.
       - output_orfs_path (str): Directory where alignment files will be stored.
       - muscle_path (str): Path to MUSCLE executable.
       - gapcheck (bool): Whether to include gaps in identity calculation.
       - identity_threshold (float): Minimum identity value (0-1) to consider ORFs as conserved.
       - max_prot_len_perc (float): Maximum length difference allowed between proteins for comparison (as a fraction).

       Returns:
       - orfs_dic (dict): Dictionary mapping unique protein sequences to their ORF information.
       - orf_in_id_list (list): List of ORF identifiers already processed.
       """
    ensure_executable(muscle_path)

    orfs_dic = {}
    orf_in_id_set = set()

    batch_counter = 0
    batch_id_list = []
    parse_time = time.time()
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()-1)
    for info_idx, seq_info in enumerate(orfs_list):
        if seq_info[0] in orf_in_id_set:
            continue
        for info_idx2, seq_info2 in enumerate(orfs_list):

            if info_idx2 <= info_idx or seq_info2[0] in orf_in_id_set:
                continue

            if not check_len(len(seq_info[-1]), len(seq_info2[-1]),max_prot_len_perc):
                continue

            batch_counter += 1
            batch_id_list.append((seq_info, seq_info2))

            if batch_counter == batch_size or info_idx == len(orfs_list) - 1:
                orfs_dic, orf_in_id_set = do_parallel_seq_parsing(
                    batch_id_list,
                    orfs_dic,
                    orf_in_id_set,
                    output_orfs_path,
                    muscle_path,
                    gapcheck,
                    identity_threshold,
                    executor=executor
                )
                batch_id_list = []
                batch_counter = 0
        print('\rParsing data... %d%% (Time elapsed: %.2f secs)' % (len(orf_in_id_set) * (100 / len(orfs_list)),
                                                                    time.time() - parse_time), end="")
    executor.shutdown()
    return orfs_dic


# Parallel processing all ORFs gathered before to identify conservation among them for the given sequences
def do_parallel_seq_parsing(
    tested_batch_set,
    orfs_dic,
    orf_in_id_list,
    output_orfs_path,
    muscle_path,
    gapcheck,
    identity_threshold,
    executor=None
    ):
    """
    Compares ORF pairs in parallel to identify conserved sequences based on sequence identity.

    Parameters:
    - tested_batch_set (list): List of ORF pairs to be aligned and compared.
    - orfs_dic (dict): Dictionary of conserved ORFs.
    - orf_in_id_list (list): List of already processed ORF identifiers.
    - output_orfs_path (str): Path to save temporary alignment files.
    - muscle_path (str): Path to MUSCLE executable.
    - gapcheck (bool): Whether to include gaps in identity calculation.
    - identity_threshold (float): Minimum identity value (0-1) to consider sequences as conserved.
    Returns:
    - orfs_dic (dict): Updated dictionary with conserved ORFs.
    - orf_in_id_list (list): Updated list of processed ORF identifiers.
    """

    data1 = [(tset,output_orfs_path,muscle_path,gapcheck,identity_threshold ) for tset in tested_batch_set]
    result_info = executor.map(parse_sequences,data1)

    for results_index, result in enumerate(result_info):
        if result[0]:
            if result[1][0][-1] in orfs_dic.keys():
                orfs_dic[result[1][0][-1]].append(result[1][1])
                orf_in_id_list.add(result[1][1][0])
            else:
                orfs_dic[result[1][0][-1]] = [result[1][0]]
                orfs_dic[result[1][0][-1]].append(result[1][1])
                orf_in_id_list.add(result[1][0][0])
                orf_in_id_list.add(result[1][1][0])
    return orfs_dic,orf_in_id_list


#Function called in the PP, compare ORFs by pairs using muscle and alignment without gaps
def parse_sequences(cur_data1):
    """
        Aligns a pair of ORFs using MUSCLE and evaluates their conservation based on identity threshold.

        Parameters:
        - cur_data1 (tuple): Tuple containing (tested_set, output_orfs_path, muscle_path, gapcheck, identity_threshold)

        Returns:
        - result (tuple): (True, tested_set) if ORFs are conserved, (False, ()) otherwise.
        """
    tested_set, output_orfs_path, muscle_path, gapcheck, identity_threshold = cur_data1
    try:
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_in, \
                tempfile.NamedTemporaryFile(mode="r+", delete=False) as temp_out:

                    temp_in.write(f">{tested_set[0][0]}\n{tested_set[0][-1]}\n")
                    temp_in.write(f">{tested_set[1][0]}\n{tested_set[1][-1]}\n")
                    temp_in.flush()

                    subprocess.run([muscle_path,
                                    '-align', temp_in.name,
                                    '-output', temp_out.name,
                                    '--quiet'],
                                   check=True)

                    temp_out.seek(0)
                    aln_content = temp_out.read()

    except subprocess.CalledProcessError as e:
        print(f"[MUSCLE ERROR] {e}")
        return False, ()

    finally:
        os.remove(temp_in.name)
        os.remove(temp_out.name)

    align_result_dic = {}
    current_id = None
    current_seq = []

    for line in aln_content.strip().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                align_result_dic[current_id] = "".join(current_seq)
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        align_result_dic[current_id] = "".join(current_seq)

    seq1,seq2 = list(align_result_dic.values())[0],list(align_result_dic.values())[1]

    arr1 = np.frombuffer(seq1.encode(), dtype='S1')
    arr2 = np.frombuffer(seq2.encode(), dtype='S1')

    hit_counter = (arr1 == arr2).sum()
    gap_counter = (arr1 == b'-').sum()
    if not gapcheck:
        identity = hit_counter/(len(seq1)-gap_counter)
    else:
        identity = hit_counter/len(seq1)

    if identity > identity_threshold:
        return True,tested_set
    else:
        return False,()


def save_results(orfs_dic,output_path,output_name,identifiers_total,output_orfs_path,comp_seq_bool,id_list):
    """
        Saves the results of the ORF comparison, including detailed information on conserved ORFs and sequences.

        Parameters:
        - orfs_dic (dict): Dictionary of conserved ORFs with metadata.
        - output_path (str): Path to save final result files.
        - output_name (str): Base name to use in output files.
        - identifiers_total (int): Total number of genomes analyzed.
        - output_orfs_path (str): Directory where individual ORF FASTA files will be written.
        - comp_seq_bool (bool): Whether to save the complementary nucleotide sequence.
        - id_list (list): List of all genome identifiers.

        Returns:
        - None: Writes output files to disk.
        """
    full_outfile = open(output_path + output_name + "_All_info.txt", "w")
    only_conserved_orfs_outfile = open(output_path + output_name + "_ConservedOnly.txt", "w")

    full_outfile.write("Run done at: " + (datetime.now()).strftime("%d/%m/%Y %H:%M") + ".")
    full_outfile.write("*Identifiers information: Length, Strand, Frame and Nucleotide (start and end) positions. *\n")
    consev_counter = 0

    orfs_by_len_list = sorted((list(orfs_dic.keys())), key=len, reverse=True)

    for unique_aa_seq in orfs_by_len_list:
        orf_info_list = orfs_dic[unique_aa_seq]

        conserved_bool = False
        if identifiers_total <= len(orf_info_list):
            consev_counter += 1
            full_outfile.write("\n### Conserved ORF among the genomes. ###\n")
            if identifiers_total < len(orf_info_list):
                full_outfile.write("\n*  The ORF is  or has an early stopping in at least one genome!  *\n")
            conserved_bool = True

        len_values_list = [len(orf_set[-1]) for orf_set in orf_info_list]

        median_len = round(np.median(len_values_list), 1)
        mean_len = round(np.mean(len_values_list), 1)
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

            # Get the identifier and theorf position
            id_info_list = cur_genome_orf_info[0].split("_orf_")
            first_id1, orf_pos = id_info_list[0], id_info_list[1]

            if cur_genome_orf_info[1] == 1:
                full_info_id = (first_id1 + "__L" + str(len(cur_genome_orf_info[-1])) +
                                "_S" + str(cur_genome_orf_info[1]) + "_F" + str(cur_genome_orf_info[2]) +
                                "_NtPos" + str(cur_genome_orf_info[3][0]) + "-" + str(cur_genome_orf_info[3][1]))
                sense_seq = cur_genome_orf_info[3][2]
                comp_seq = cur_genome_orf_info[4][2]
            else:
                full_info_id = (first_id1 + "__L" + str(len(cur_genome_orf_info[-1])) +
                                "_S" + str(cur_genome_orf_info[1]) + "_F" + str(cur_genome_orf_info[2]) +
                                "_NtPos" + str(cur_genome_orf_info[4][0]) + "-" + str(cur_genome_orf_info[4][1]))
                sense_seq = cur_genome_orf_info[4][2]
                comp_seq = cur_genome_orf_info[3][2]

                # Add orfs to fasta file using the identifier as the name
            write_output(output_orfs_path + output_name + "_" + first_id1 + "_AA.fasta",
                         ">" + full_info_id + "_ORFpos" + orf_pos + "\n" + cur_genome_orf_info[-1] + "\n")

            write_output(output_orfs_path + output_name + "_" + first_id1 + "_nt.fasta",
                         ">" + full_info_id + "\n" + sense_seq + "\n")
            if comp_seq_bool:
                write_output(output_orfs_path + output_name + "_" + first_id1 + "_nt.fasta",
                             ">" + full_info_id + "_Comp\n" + comp_seq + "\n")

            in_ids_aux_list.append(first_id1)
            full_outfile.write(">" + full_info_id + "\n" + cur_genome_orf_info[-1] + "\n")

            if conserved_bool:
                only_conserved_orfs_outfile.write(">" + full_info_id + "\n" + str(cur_genome_orf_info[-1]) + "\n")

        full_outfile.write("***\n")
        # if there is duplicated orfs
        orf_set = set(in_ids_aux_list)
        if len(orf_set) > len(id_list):
            full_outfile.write("Genomes with multiple instances of the ORF:\n")
            for i in in_ids_aux_list:
                count = in_ids_aux_list.count(i)
                if count > 1:
                    full_outfile.write(i + "\n")
                    break
        elif len(orf_set) < len(id_list):
            # if there is some genome don't have the ORF
            full_outfile.write("Genomes without the ORF:\n")
            for i in id_list:
                if i not in in_ids_aux_list:
                    full_outfile.write(i + "\n")
        full_outfile.write("--------------------------------\n")
    full_outfile.close()
    only_conserved_orfs_outfile.close()

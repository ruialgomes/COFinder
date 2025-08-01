# Function to open a FASTA file and read its contents
from src.utils import change_chars

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


# Function to generate reverse complement of a DNA sequence
def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp_seq = ''.join([complement[base] for base in dna_sequence[::-1]])
    return reverse_comp_seq


# Function to translate DNA sequence to protein sequence
def translate_dna_to_protein(dna_sequence):
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

def format_sequences(init_input_genomes_dic):
    input_genomes_dic = {}
    for k, v in init_input_genomes_dic.items():
        k1 = change_chars(k)
        if not all(base in 'ACTG' for base in v):
            print("Genome with non-canonical base, removed: ", k)
            continue
        input_genomes_dic[k1] = v
    return input_genomes_dic

# Function to find open reading frames with translations
def find_orfs_with_trans(seq, check_alt_codon, min_pep_size):
    if not check_alt_codon:
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
                if nt_start_codon_pos != -1 and (len(aa_seq) - aaseq_start_pos) >= min_pep_size:
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
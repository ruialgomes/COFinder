# COFinder
Conserved ORFs Finder is a versatile tool designed for grouping conserved Open Reading Frames (ORFs) within a sequence set. It facilitates the identification of conserved regions based on sequence similarity, offering customizable parameters such as minimum protein size adjustment for detection, or the option to include non-canonical start codons. 

An open reading frame (ORF) is defined as an in-frame DNA sequence of codons (three consecutive nucleotides) between a start codon and a stop codon. This structural conservation made ORF detection the purpose of some of the first computational tools developed for DNA sequence analysis. In general, these tools' default parameters search for ORFs encoding at least 100 amino acids (aa), in consequence of their higher statistical confidence.
Recently, discoveries related to small ORFs are revealing them to be more present than thought before on the organism's proteome. They are, as a class, also enriched in sequence properties previously assumed to be unusual, including non-AUG start codons. 

Even with the recently developed tools using ML approaches to find this sORFs, there is still a challenge to find them. Trying to facilitate this research from another 'side', here is a simple tool for finding these sORFs using their existence in all available genomes within a group of exemplars, as suggested in some papers as an interesting feature to be evaluated. This is done by extracting all ORFs from each given sequence, using both aug or non-aug start codons and the first found stop codon as the ORF 'borders', checking for the minimum protein length defined, and saving them. Then all found ORF are sorted and compared, using comparative size as an excluding factor, by alignment with muscle followed by obtaining the percentage of identity/similarity between a pair of proteins. At the end, all similar proteins are saved (with some information about it on the identifier) in a file of groups of proteins in the output directory, and another file with only ORFs conserved among all sequences given as input. A directory of all found orfs for each input sequence is also created in a directory named "ORFs".

For any suggestions or doubts, please contact me at ruitherarthur@gmail.com.

# Installing and usage
*  Download the folder “COFinder”, and extract it from the zip file at the desired directory. To run the script open the terminal, and navigate to the extraction directory. With Python3 already installed, inside the folder, you will create a new environment for TensorFlow installation. You can do that by using the command:

*python -m venv ’environment_name’* for **Windows**. or *python3 -m venv ’environment_name’* for **Unix/Mac**.

*  Where you will change “environment_name” to a name of your preference, like “concret_venv”. Then you need activate the environment using:

*’environment_name’/Scripts/activate* for **Windows** or *source ’environment_name’/bin/activate* for **Unix/Mac**

*  Then you can use the command line options, with more options below:

*python COFinder_v1.py --in “yourFasta.fasta”* for **Windows** OR *python3 COFinder_v1.py --in “yourFasta.fasta”* for **Unix/Mac**

*  When you finish, deactivate the venv with:

*deactivate*

# Command line options

--input: The path to the input file. Defaults to a file containing 12 "Geminiviridae" alphasatellites if not specified.

--out: Specifies the output folder name. If the folder doesn't exist, it creates a folder named "output" to store results.

--muscle_path: Path to an alternative muscle version. The default is version 5.1 in the bin folder.

--pep_min_size: Minimum number of amino acids required for a valid peptide. Defaults to 6 amino acids.

--alt_codon: Set to True to use alternative start codons (CTG, GTG, TTG, ATT, ATC). The default is to use only ATG.

--use_gaps: Includes gaps for identity calculation between proteins if used. The default is removing the gaps from.

--prot_len_perc: Percentage (in decimals) of the ORF's length accepted for variation compared to the representative sequence. Default: 0.6.

--identity_threshold: Percentage of similarity required between ORFs to be considered in the same group. Default: 0.7 (70% similarity).

--nproc: Number of processing units used for parallel processing. Defaults to the number of available processors.

--comp_seq: Determines if the complementary sequence for both sense and anti-sense frames will be shown in the complete output. Default: False. Set to True to show.

Command line example

*python3 COFinder_v1.py --in ’path/to/your/file.fasta’*

import numpy as np
import os, sys
import os.path

def read_fasta_file(filename):
    """
    Reads the given FASTA file f and returns a dictionary of sequences.

    Lines starting with ';' in the FASTA file are ignored.
    """
    sequences_lines = {}
    current_sequence_lines = None
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
                current_sequence_lines = []
                sequences_lines[sequence_name] = current_sequence_lines
            else:
                if current_sequence_lines is not None:
                    current_sequence_lines.append(line)
    sequences = {}
    for name, lines in sequences_lines.items():
        sequences[name] = ''.join(lines)
    return sequences

def read_score_matrix_and_alphabet(filename):
    score_matrix = {}
    symbols = []
    with open(filename) as fp:
        first_line = True
        lines = []
        for line in fp:
            if first_line:
                first_line = False
                continue
            lines.append(line)
            symbols.append(line[0])

        for i in range(0,len(lines)):
            char_set = lines[i].split("  ")
            for j in range(1, len(char_set)):
                score_matrix[(symbols[i], symbols[j-1])] = int(char_set[j])
    return score_matrix,symbols

def parse_arguments(args):
    score_matrix_file = args[1]
    gap_cost = int(args[2])
    file1 = args[3]
    file2 = args[4]
    should_output_allignment = int(args[5])
    score_matrix,alphabet = read_score_matrix_and_alphabet(score_matrix_file)
    fastaSeq1 = ""
    fastaSeq2 = ""

    if os.path.isfile(file1):
        fastaSeq1 = read_fasta_file(file1)
    else:
        fastaSeq1 = file1

    if os.path.isfile(file2):
        fastaSeq2 = read_fasta_file(file2)
    else:
        fastaSeq2 = file2

    return score_matrix, gap_cost, should_output_allignment, alphabet, fastaSeq1, fastaSeq2
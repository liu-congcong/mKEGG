import os
from math import inf
from multiprocessing import Process

from .fasta_utility import split_fasta
from .hmmsearch_utility import run_hmmsearch
from .make_file import make_file
from .prodigal_utility import combine_files, predict_protein, train_model


def read_hit_file(input_file, score_hash):
    open_file = open(input_file, 'r')
    for line in open_file:
        if not line.startswith('#'):
            lines = line.rstrip('\n').split(maxsplit = 9)
            score_index, score_threshold = score_hash[lines[2]]
            if float(lines[score_index]) >= score_threshold:
                yield (lines[0], lines[2])
    open_file.close()
    return None


def read_ko_list(input_file):
    k_number2score = dict()
    open4r = open(input_file, 'r')
    open4r.readline()
    # knum threshold score_type profile_type F-measure nseq nseq_used alen mlen eff_nseq re/pos definition #
    for line in open4r:
        k_number, score, score_type, _ = line.rstrip('\n').split('\t', maxsplit = 3)
        if score_type == 'full':
            k_number2score[k_number] = (5, float(score))
        elif score_type == 'domain':
            k_number2score[k_number] = (8, float(score))
        else:
            k_number2score[k_number] = (5, -inf)
    open4r.close()
    return k_number2score


def output_kmap(hit_files, k_number2score, sequence_id2position, output_file):
    open_file = open(output_file, 'w')
    open_file.write('\t'.join(['Sequence ID', 'Start', 'End', 'K number']) + '\n')
    for hit_file in hit_files:
        for sequence_id, k_number in read_hit_file(hit_file, k_number2score):
            start, end = sequence_id2position[sequence_id]
            open_file.write('\t'.join([sequence_id.rsplit('_', maxsplit = 1)[0], start, end, k_number]) + '\n')
        os.remove(hit_file)
    open_file.close()
    return None


def mkegg_mapping(parameters):
    # Train model using Prodigal. #
    model_file = train_model(parameters.prodigal, parameters.assembly)

    processes = list()
    protein_files = list()
    for fasta_file in split_fasta(parameters.assembly, parameters.threads):
        protein_files.append(make_file())
        # Predict protein sequences using Prodigal. #
        processes.append(
            Process(
                target = predict_protein,
                args = (parameters.prodigal, model_file, fasta_file, protein_files[-1])
            )
        )
        processes[-1].start()
    for process in processes:
        process.join()
    os.remove(model_file)
    sequence_id2position = combine_files(protein_files, parameters.output + '.proteins')

    hit_files = list()
    for fasta_file in split_fasta(parameters.output + '.proteins', parameters.threads):
        hit_files.append(make_file())
        # Map proteins to KEGG hmms using hmmsearch. #
        processes.append(
            Process(
                target = run_hmmsearch,
                args = (parameters.hmmsearch, parameters.k_hmm, fasta_file, hit_files[-1])
            )
        )
        processes[-1].start()
    for process in processes:
        process.join()
    k_number2score = read_ko_list(parameters.ko_list)
    output_kmap(hit_files, k_number2score, sequence_id2position, parameters.output)
    return None
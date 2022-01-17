import os
from subprocess import DEVNULL, run

from .make_file import make_file


def train_model(prodigal, input_fasta):
    output_model = make_file()
    run_process = run([prodigal, '-g', '11', '-i', input_fasta, '-m', '-p', 'single', '-t', output_model], stdout = DEVNULL, stderr = DEVNULL)
    assert not run_process.returncode, 'Run Prodigal field.'
    return output_model


def predict_protein(prodigal, model_file, input_fasta, output_fasta):
    run_process = run([prodigal, '-a', output_fasta, '-g', '11', '-i', input_fasta, '-m', '-p', 'single', '-t', model_file], stdout = DEVNULL, stderr = DEVNULL)
    assert not run_process.returncode, 'Run Prodigal field.'
    os.remove(input_fasta)
    return None


def combine_files(input_files, output_file):
    sequence_id2position = dict()
    open4w = open(output_file, 'w')
    for input_file in input_files:
        open4r = open(input_file, 'r')
        for line in open4r:
            line = line.rstrip('\n')
            if line.startswith('>'):
                sequence_id, start, end, _ = line[1 : ].split(' # ', maxsplit = 3)
                sequence_id2position[sequence_id] = (start, end)
            open4w.write(line + '\n')
        open4r.close()
        os.remove(input_file)
    open4w.close()
    return sequence_id2position

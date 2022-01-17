import os
from math import ceil
from multiprocessing import Process
from subprocess import DEVNULL, run

from .make_file import make_file
from .samtools_utility import run_samtools


def split_kmap(input_kmap, output_files):
    open4r = open(input_kmap, 'r')
    open4r.readline()
    total_lines = 0
    for line in open4r:
        total_lines += 1
    open4r.seek(0, 0)
    open4r.readline()
    lines_per_file = ceil(total_lines / output_files)
    for output_file_index in range(0, total_lines, lines_per_file):
        output_file = make_file()
        open4w = open(output_file, 'w')
        for read_line in range(lines_per_file):
            line = open4r.readline().rstrip('\n')
            if not line:
                break
            sequence_id, start, end, k_number  = line.split('\t')
            open4w.write('\t'.join([sequence_id, str(int(start) - 1), end, k_number]) + '\n')
        open4w.close()
        yield output_file
    open4r.close()
    return None


def output_kmap_with_abundance(input_files, output_file):
    open4w = open(output_file, 'w')
    open4w.write('\t'.join(['Sequence ID', 'Start', 'End', 'K number', 'Abundance']) + '\n')
    for input_file in input_files:
        open4r = open(input_file, 'r')
        for line in open4r:
            lines = line.rstrip('\n').split('\t')
            start = int(lines[1])
            end = int(lines[2])
            open4w.write('\t'.join([lines[0], str(start + 1), lines[2], lines[3], str(float(lines[4]) / (end - start))]) + '\n')
        open4r.close()
        os.remove(input_file)
    open4w.close()
    return None


def mkegg_abundance(parameters):
    processes = list()
    kmap_files_with_abundances = list()
    for bed_file in split_kmap(parameters.kmap, parameters.threads):
        kmap_files_with_abundances.append(make_file())
        # Calculate abundance using samtools. #
        processes.append(
            Process(
                target = run_samtools,
                args = (parameters.samtools, parameters.bam, bed_file, kmap_files_with_abundances[-1])
            )
        )
        processes[-1].start()
    for process in processes:
        process.join()
    output_kmap_with_abundance(kmap_files_with_abundances, parameters.output)
    return None

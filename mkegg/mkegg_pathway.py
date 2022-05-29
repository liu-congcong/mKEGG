from json import load

import numpy

from .fasta_utility import read_fastas


def read_kmap(input_file):
    k_index2sequence_id2abundance = list()
    k_number2k_index = dict()
    k_index = 0
    open_file = open(input_file, 'r')
    open_file.readline()
    for line in open_file:
        lines = line.rstrip('\n').split('\t')
        sequence_id = lines[0]
        k_number = lines[3]
        if k_number not in k_number2k_index: # a new k_number #
            k_number2k_index[k_number] = k_index
            k_index2sequence_id2abundance.append(dict())
            k_index += 1
        sequence_id2abundance = k_index2sequence_id2abundance[k_number2k_index[k_number]]
        sequence_id2abundance[sequence_id] = sequence_id2abundance.get(sequence_id, 0.0) + float(lines[4])
    open_file.close()
    return (k_number2k_index, k_index2sequence_id2abundance)


def read_ko_json(input_file, k_number2k_index):
    levela2levelb2levelc2k_indices = dict()
    open_file = open(input_file, 'r')
    root_node = load(open_file)
    open_file.close()
    for node_a in root_node.get('children', list()):
        level_a = node_a['name'].split(maxsplit = 1)[1]
        levela2levelb2levelc2k_indices[level_a] = dict()
        levelb2levelc2k_indices = levela2levelb2levelc2k_indices[level_a]
        for node_b in node_a.get('children', list()):
            level_b = node_b['name'].split(maxsplit = 1)[1]
            levelb2levelc2k_indices[level_b] = dict()
            levelc2k_indices = levelb2levelc2k_indices[level_b]
            for node_c in node_b.get('children', list()):
                level_c = node_c['name'].split(maxsplit = 1)[1]
                levelc2k_indices[level_c] = list()
                k_indices = levelc2k_indices[level_c]
                for node_d in node_c.get('children', list()):
                    k_number = node_d['name'].split(maxsplit = 1)[0]
                    if k_number in k_number2k_index:
                        k_indices.append(k_number2k_index[k_number])
    return levela2levelb2levelc2k_indices


def calculate_relative_abundance(abundance, k_index2sequence_id2abundance, sequence_id2cluster_index):
    for k_index, sequence_id2abundance in enumerate(k_index2sequence_id2abundance):
        for sequence_id, abundance_ in sequence_id2abundance.items():
            if sequence_id in sequence_id2cluster_index:
                cluster_index = sequence_id2cluster_index[sequence_id]
                abundance[k_index][cluster_index] += abundance_
    abundance /= numpy.sum(abundance, axis = 0) + numpy.finfo(numpy.float64).eps
    abundance *= 100
    return abundance


def output_pathway(abundance, levela2levelb2levelc2k_indices, cluster_names, output_file, output_list):
    levela_k_indices = list()
    levelb_k_indices = list()
    open_file = open(output_file, 'w')
    if output_list:
        open_file.write('\t'.join(['Cluster', 'KEGG level', 'KEGG level 1', 'KEGG level 2', 'KEGG level 3', 'Abundance']) + '\n')
    else:
        open_file.write('\t'.join(['KEGG level', 'KEGG level 1', 'KEGG level 2', 'KEGG level 3'] + cluster_names) + '\n')
    for levela, levelb2levelc2k_indices in levela2levelb2levelc2k_indices.items():
        for levelb, levelc2k_indices in levelb2levelc2k_indices.items():
            for levelc, k_indices in levelc2k_indices.items():
                levelb_k_indices.extend(k_indices)
                abundance_c = numpy.sum(abundance[k_indices], axis = 0)
                if output_list:
                    for index, cluster_name in enumerate(cluster_names):
                        open_file.write('\t'.join([cluster_name, '3', levela, levelb, levelc, str(abundance_c[index])]) + '\n')
                else:
                    open_file.write('\t'.join(['3', levela, levelb, levelc] + abundance_c.astype(numpy.str_).tolist())+ '\n')
            abundance_b = numpy.sum(abundance[list(set(levelb_k_indices))], axis = 0)
            levela_k_indices.extend(levelb_k_indices)
            levelb_k_indices.clear()
            if output_list:
                for index, cluster_name in enumerate(cluster_names):
                    open_file.write('\t'.join([cluster_name, '2', levela, levelb, '-', str(abundance_b[index])]) + '\n')
            else:
                open_file.write('\t'.join(['2', levela, levelb, '-'] + abundance_b.astype(numpy.str_).tolist())+ '\n')
        abundance_a = numpy.sum(abundance[list(set(levela_k_indices))], axis = 0)
        levela_k_indices.clear()
        if output_list:
            for index, cluster_name in enumerate(cluster_names):
                open_file.write('\t'.join([cluster_name, '1', levela, '-', '-', str(abundance_a[index])]) + '\n')
        else:
            open_file.write('\t'.join(['1', levela, '-', '-'] + abundance_a.astype(numpy.str_).tolist())+ '\n')
    open_file.close()
    return None


def mkegg_pathway(parameters):
    k_number2k_index, k_index2sequence_id2abundance = read_kmap(parameters.kmap)
    levela2levelb2levelc2k_indices = read_ko_json(parameters.ko_json, k_number2k_index)
    sequence_id2cluster_index = read_fastas(parameters.clusters)
    abundance = numpy.zeros((len(k_number2k_index), len(parameters.clusters)), dtype = numpy.float64)
    calculate_relative_abundance(abundance, k_index2sequence_id2abundance, sequence_id2cluster_index)
    output_pathway(abundance, levela2levelb2levelc2k_indices, parameters.clusters, parameters.output, parameters.output_list)
    return None

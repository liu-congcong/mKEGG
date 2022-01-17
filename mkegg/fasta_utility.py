import os
from math import ceil

from .make_file import make_file


def read_fastas(input_files):
    sequence_id2file_index = dict()
    for file_index, input_file in enumerate(input_files):
        open_file = open(input_file, 'r')
        for line in open_file:
            if line.startswith('>'):
                sequence_id2file_index[line.rstrip('\n').split(maxsplit = 1)[0][1 : ]] = file_index
        open_file.close()
    return sequence_id2file_index


def split_fasta(input_file, output_files):
    total_size = os.path.getsize(input_file)
    block_size = ceil(total_size / output_files) # block_size <= total_size #
    file_position = 0
    file_position_ = 0
    open4r = open(input_file, 'rb')
    while file_position < total_size:
        line = open4r.readline()
        file_position += len(line)
        if line.startswith(b'>'):
            file_position -= len(line)
            if file_position > 0:
                open4r.seek(file_position_, os.SEEK_SET)
                output_file = make_file()
                open4w = open(output_file, 'wb')
                while file_position_ < file_position:
                    file_position_ += open4w.write(open4r.read(min(10485760, file_position - file_position_)))
                open4w.close()
                yield output_file
                # file_position_ will be equal to file_position, open4r.tell() will be equal to file_position #
            file_position = open4r.seek(min(file_position + block_size, total_size), os.SEEK_SET)
    open4r.seek(file_position_, os.SEEK_SET)
    output_file = make_file()
    open4w = open(output_file, 'wb')
    while file_position_ < file_position:
        file_position_ += open4w.write(open4r.read(min(10485760, file_position - file_position_)))
    open4w.close()
    yield output_file
    open4r.close()
    return None

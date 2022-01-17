import os
from tempfile import mkstemp


def make_file():
    file_descriptor, file_name = mkstemp(dir = os.getcwd())
    os.close(file_descriptor)
    os.remove(file_name)
    return file_name
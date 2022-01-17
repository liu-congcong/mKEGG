import os
from subprocess import run


def run_samtools(samtools, input_bam, input_bed, output_file):
    open_file = open(output_file, 'wb')
    run_process = run([samtools, 'bedcov', input_bed, input_bam], stdout = open_file)
    assert not run_process.returncode, 'Run samtools field.'
    os.remove(input_bed)
    return None

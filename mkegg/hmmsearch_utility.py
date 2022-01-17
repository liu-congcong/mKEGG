import os
from subprocess import DEVNULL, run


def run_hmmsearch(hmmersearch, hmm_file, protein_file, hit_file):
    run_process = run([hmmersearch, '--tblout', hit_file, '--noali', '--notextw', '-T', '0', '--cpu', '1', hmm_file, protein_file], stdout = DEVNULL)
    assert not run_process.returncode, 'Run hmmsearch field.'
    os.remove(protein_file)
    return None
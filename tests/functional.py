import os
import sys
import tempfile
import subprocess
import tempfile
import shutil
import multiprocessing
import pytest

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
EXAMPLES_DIR = os.path.join(CURRENT_DIR, '..', 'examples')
DB_DIR = os.path.join(CURRENT_DIR, '..', 'db')

sys.path.append(SCRIPTS_DIR)

from binary_utils import Binary

def get_matam_stats(f):
    res_handler = open(f, 'r')
    lines = res_handler.readlines()
    i = None
    for i, line in enumerate(lines):
        if 'One-line stats' in line: break
    stats = lines[i+1].strip()
    stats = stats.split('\t')
    stats = stats[-8:]
    #stats = stats[:4] + stats[5:]
    formated_stats = []
    for s in stats:
        s = s.strip()
        if '%' in s:
            s = s.replace('%','')
        if '.' in s:
            s = float(s)
        else: s = int(s)
        formated_stats.append(s)

    res_handler.close()
    return formated_stats

@pytest.mark.skipif(not os.path.isdir(DB_DIR),
                    reason= 'DB_DIR is missing:%s' % DB_DIR)
@pytest.mark.skipif(Binary.which('exonerate') is None,
                    reason= 'Exonerate have to be in your PATH:%s' % sys.path)
def test_matam():

    out = tempfile.mkdtemp(dir='/tmp/', prefix='matam_functionnal_test_')
    p = {
        'bin' : os.path.join(SCRIPTS_DIR, 'matam_assembly.py'),
        'reads' : os.path.join(EXAMPLES_DIR, '16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq'),
        'db' : os.path.join(DB_DIR, 'SILVA_128_SSURef_NR95'),
        'true_ref' : os.path.join(EXAMPLES_DIR, '16sp_simulated_dataset/16sp.fasta'),
        'true_ref_taxo' : os.path.join(EXAMPLES_DIR, '16sp_simulated_dataset/16sp.taxo.tab'),
        'out' : out,
        'log' : os.path.join(out, 'matam.log'),
        'cpu' : multiprocessing.cpu_count()
    }
    cmd = '{bin} -i {reads} -d {db} -o {out} --true_references {true_ref} --true_ref_taxo {true_ref_taxo} --cpu {cpu} --max_memory 10000 --debug --perform_taxonomic_assignment > {log} 2>&1'.format(**p)
    print(cmd)
    try:
        completed_process = subprocess.run(cmd, shell=True)
        assert completed_process.returncode == 0
        stats = get_matam_stats(p['log'])
        error_rate, error_rate2 = stats[5:7]
        assert error_rate <= 0.16
        assert error_rate2 <= 0.13

    finally:
        if os.path.isdir(out):
            shutil.rmtree(out)

import os
import sys
import tempfile
import subprocess
import shutil
import multiprocessing
import pytest
import collections

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
EXAMPLES_DIR = os.path.join(CURRENT_DIR, '..', 'examples')
DB_DIR = os.path.join(CURRENT_DIR, '..', 'db')

sys.path.append(SCRIPTS_DIR)

from binary_utils import Binary

# skip all module tests if needed
pytestmark = pytest.mark.skipif(
    not os.path.isdir(DB_DIR),
    reason='DB_DIR is missing:%s' % DB_DIR
)


@pytest.fixture(scope='module')
def matam_results():
    out = tempfile.mkdtemp(dir='/tmp/', prefix='matam_functionnal_test_')
    p = {
        'bin': os.path.join(SCRIPTS_DIR, 'matam_assembly.py'),
        'reads': os.path.join(
            EXAMPLES_DIR,
            '16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq'
        ),
        'db': os.path.join(DB_DIR, 'SILVA_128_SSURef_NR95'),
        'out': out,
        'cpu': multiprocessing.cpu_count()
    }
    cmd = '{bin} -i {reads} -d {db} -o {out} --cpu {cpu} --max_memory 3000 \
        --debug --coverage_threshold 2000 \
        --perform_taxonomic_assignment'.format(**p)

    completed_process = subprocess.run(cmd, shell=True)
    return_code = completed_process.returncode
    fasta = os.path.join(out, 'final_assembly.fa')
    krona_html = os.path.join(out, 'krona.html')
    krona_tab = os.path.join(out, 'krona.tab')
    rdp_tab = os.path.join(out, 'rdp.tab')

    MatamResults = collections.namedtuple(
        "MatamResults",
        "return_code fasta krona_html krona_tab rdp_tab"
    )
    results = MatamResults(
        return_code=return_code,
        fasta=fasta,
        krona_html=krona_html,
        krona_tab=krona_tab,
        rdp_tab=rdp_tab
    )

    yield results

    if os.path.isdir(out):
        shutil.rmtree(out)


def exists_and_not_empty(fpath):
    return os.path.isfile(fpath) and os.stat(fpath).st_size != 0


def test_return_code(matam_results):
    assert matam_results.return_code == 0


def test_final_fasta_file(matam_results):
    assert exists_and_not_empty(matam_results.fasta)


def test_krona_html(matam_results):
    assert exists_and_not_empty(matam_results.krona_html)


def test_krona_tab(matam_results):
    assert exists_and_not_empty(matam_results.krona_tab)


def test_rdp_tab(matam_results):
    assert exists_and_not_empty(matam_results.rdp_tab)


def extract_metaquast_val(tsv):
    with open(tsv, 'r') as tsv_handler:
        lines = tsv_handler.readlines()
        return float(lines[1].split('\t')[1].strip())


@pytest.mark.skipif(
    not Binary.which('metaquast.py'),
    reason="requires metaquast.py to be in PATH"
)
def test_metaquast(matam_results):
    data_directory = tempfile.mkdtemp(dir='/tmp/', prefix='metaquast_')
    fasta = matam_results.fasta
    true_ref = os.path.join(EXAMPLES_DIR, '16sp_simulated_dataset/16sp.fasta')
    cmd = "metaquast.py -a all --ambiguity-score 1 --min-identity 97 -x 500 \
        --unaligned-part-size 200 -R %s %s" % (true_ref, fasta)
    subprocess.run(cmd, shell=True, cwd=data_directory)

    genome_fraction_file = os.path.join(
        data_directory,
        'quast_results/latest/summary/TSV/Genome_fraction.tsv'
    )
    mismatches_file = os.path.join(
        data_directory,
        'quast_results/latest/summary/TSV/num_mismatches_per_100_kbp.tsv'
    )
    indels_file = os.path.join(
        data_directory,
        'quast_results/latest/summary/TSV/num_Ns_per_100_kbp.tsv'
    )
    ns_file = os.path.join(
        data_directory,
        'quast_results/latest/summary/TSV/num_Ns_per_100_kbp.tsv'
    )

    genome_fraction = extract_metaquast_val(genome_fraction_file)
    mismatches = extract_metaquast_val(mismatches_file)
    indels = extract_metaquast_val(indels_file)
    ns = extract_metaquast_val(ns_file)
    error_rate = (mismatches + indels + ns) / 1000  # 100000bp * 100

    assert genome_fraction > 86.4
    assert error_rate < 0.15

    if os.path.isdir(data_directory):
        shutil.rmtree(data_directory)

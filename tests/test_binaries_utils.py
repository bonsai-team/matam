import os
import sys

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
sys.path.append(SCRIPTS_DIR)


import pytest
from binary_utils import Binary

@pytest.mark.parametrize("valid_binary", ['componentsearch',
                                            'ovgraphbuild',
                                            'sga',
                                            'indexdb_rna',
                                            'sortmerna',
                                            'vsearch',
                                            'java', #rdp
                                            'ktImportText', #krona
                                            ])

def test_binary_which_ok(valid_binary):
    assert Binary.which(valid_binary) is not None

def test_binary_which_ko():
    assert Binary.which('this_is_not_a_valid_bin') is None

def test_binary_assert_which_ko():
    with pytest.raises(SystemExit) as se:
        Binary.assert_which('this_is_not_a_valid_bin')
    assert 'No valid' in str(se.value)

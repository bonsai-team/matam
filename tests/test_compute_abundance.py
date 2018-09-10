import os
import sys
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
sys.path.append(SCRIPTS_DIR)


SAMPLE_DIR = os.path.join(CURRENT_DIR, 'sample')

from compute_abundance import abundance_calculation
import pytest

@pytest.mark.parametrize('blast,expected_abundance',

    [    # Basic test
        ['scaffolds.blast',
         {'scaff1': 1.75,
          'scaff2': 3.75,
          'scaff3': 1.25,
          'scaff4': 1.25
         }
        ],
        # Test the case where a read can be found several times on the same scaffolds.
        # this behavior is tolerated but is not intented to occur often
        ['scaffolds_multiple_reads.blast',
         { '159': 0.5,
           '161': 1,
           '175': 0.5,
           '240': 3
         },
        ]
    ]
)

def test_abundance_calculation(blast, expected_abundance):
    blast_file = os.path.join(SAMPLE_DIR, blast)
    abundance = abundance_calculation(blast_file)

    assert set(abundance.keys()) == set(expected_abundance.keys())
    assert sorted(abundance.values()) == pytest.approx(sorted(expected_abundance.values()))

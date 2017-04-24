#!/usr/bin/env python3

import os
import sys


class Binary:
    """
    This class aims to retrieve binaries by name.
    First, look into matam directories (cf below).
    Then look into the PATH.
    If no match occured, return None

    Matam organisation:

    matam_root_dir/
        /bin/
        /scripts/
            binary_utils.py
        /componentsearch/
        /ovgraphbuild/bin/
        /sga/src/
        /sortmerna/
        /vsearch/bin/
        /RDPTools/
        /Krona/KronaTools/scripts/
    """

    _matam_root, _ = os.path.split(os.path.dirname(os.path.realpath(__file__)))
    _matam_bin_dirs = [ 'bin/',
                        'scripts/',
                        'componentsearch/',
                        'ovgraphbuild/bin/',
                        'sga/src/SGA',
                        'sortmerna/',
                        'vsearch/bin/',
                        'RDPTools/',
                        'Krona/KronaTools/scripts/'
    ]


    def _is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


    def which(program):
        path_to_look_for = [ os.path.join(Binary._matam_root, d) for d in Binary._matam_bin_dirs ]
        path_to_look_for.extend(os.environ["PATH"].split(os.pathsep))

        for path in path_to_look_for:
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if Binary._is_exe(exe_file):
                return exe_file

    def assert_which(program):
        p = Binary.which(program)
        if p:
            return p
        sys.exit('No valid binary found for %s' % program)

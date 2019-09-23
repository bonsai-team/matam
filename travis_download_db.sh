#!/bin/bash

wget --quiet http://bioinfo.lifl.fr/matam/SILVA_128_SSURef_NR95_indexed_max_mem_3G.tar.xz{.md5,}
md5sum -c SILVA_128_SSURef_NR95_indexed_max_mem_3G.tar.xz.md5
tar Jxvf SILVA_128_SSURef_NR95_indexed_max_mem_3G.tar.xz

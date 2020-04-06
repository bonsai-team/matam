#!/usr/bin/env python3

def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    # Variables initialization
    count = 0
    header = ''
    seq = ''
    qual = ''
    # Reading input file
    for line in (l.strip() for l in fastq_file_handle if l.strip()):
        count += 1
        if count % 4 == 1:
            if header:
                yield header, seq, qual
            header = line[1:].split()[0]
        elif count % 4 == 2:
            seq = line
        elif count % 4 == 0:
            qual = line
    # yield last fastq sequence
    yield header, seq, qual
    # Close input file
    fastq_file_handle.close()


def is_phred33(fastq, number_of_reads_to_test=-1):
    """
    Determine if the fastq use an offset of +33.
    See for more details: https://en.wikipedia.org/wiki/FASTQ_format
    By default, test until end of file.
    """
    min_ascii_code = 100
    max_ascii_code = 0
    with open(fastq, 'r') as handler:
        for i, (header, seq, qual) in enumerate(read_fastq_file_handle(handler)):
            if number_of_reads_to_test > 0 and i >= number_of_reads_to_test:
                break

            ascii_codes = set([ord(qual_char) for qual_char in qual])
            current_min_ascii_code = min(ascii_codes)
            current_max_ascii_code = max(ascii_codes)
            if current_max_ascii_code >= 75:
                return False

            min_ascii_code = min(current_min_ascii_code, min_ascii_code)
            max_ascii_code = max(current_max_ascii_code, max_ascii_code)

    return min_ascii_code < 59 and max_ascii_code <=74

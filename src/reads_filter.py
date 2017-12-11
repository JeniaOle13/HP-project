import argparse
import gzip
import bz2
import itertools
from contextlib import contextmanager


def skip_comments(it):
    for l in it:
        if l.startswith('#'):
            continue
        yield l.strip()


def parse_qual(file_obj):
    file_obj = enumerate(skip_comments(file_obj))
    while True:
        try:
            n_tag, tag_line = next(file_obj)
            n_qual, qual_line = next(file_obj)
            try:
                quals = [int(x) for x in qual_line.split()]
            except ValueError as e:
                msg = "Error on line #{lineno}".format(lineno=n_tag)
                print(msg)
                log(msg)
                raise e
            yield tag_line, quals
        except StopIteration:
            break


def parse_csfasta(file_obj):
    file_obj = skip_comments(file_obj)
    while True:
        try:
            tag_line = next(file_obj)
            colors = next(file_obj)
            yield tag_line, colors
        except StopIteration:
            break


@contextmanager
def open_file(filename, mode):
    if filename.endswith('.gz'):
        openf = gzip.open
    elif filename.endswith('.bz2'):
        openf = bz2.open
    else:
        openf = open

    f = openf(filename, mode)
    yield f
    f.close()


def log(message):
    with open(snakemake.log[0], 'a') as f:
        f.write(message)
        f.write('\n')


def filter_reads(quality_threshold, csfasta_file_name, qual_file_name, out_csfasta_file, out_qual_file):
    reads_total = 0
    reads_filtered = 0
    missing_qualities = 0
    extra_qualities = 0

    with open_file(csfasta_file_name, 'rt') as fasta_in, \
         open_file(qual_file_name, 'rt') as qual_in, \
         open(out_csfasta_file, 'wt') as fasta_out, \
         open(out_qual_file, 'wt') as qual_out:

        csfasta_records = parse_csfasta(fasta_in)
        qual_records = parse_qual(qual_in)

        for csfasta_tag, colors in csfasta_records:
            reads_total += 1
            while True:
                # TODO: rewrite into "for in" expression
                try:
                    qual_tag, quals = next(qual_records)
                    if not qual_tag or qual_tag < csfasta_tag:
                        extra_qualities += 1
                        continue
                    break
                except StopIteration:
                    break

            if csfasta_tag != qual_tag:
                missing_qualities += 1
                reads_filtered += 1
                continue

            neg_qualities = any(x < 0 for x in quals)
            avg_qual = float(sum(quals)) / len(quals)
            if not neg_qualities and \
                    (quality_threshold == 0 or avg_qual >= quality_threshold):
                fasta_out.write('{}\n'.format(csfasta_tag))
                fasta_out.write('{}\n'.format(colors))
                qual_out.write('{}\n'.format(qual_tag))
                qual_out.write('{}\n'.format(' '.join(map(str, quals))))
            else:
                reads_filtered += 1


    if missing_qualities > 0:
        log("WARNING: {missing_qualities} reads filtered due to missing quality scores.".format(
            missing_qualities=missing_qualities
        ))

    if extra_qualities > 0:
        log("WARNING: {extra_qualities} extra quality tags ignored.".format(
            extra_qualities=extra_qualities
        ))

    log("Filtered {reads_filtered} of {reads_total} reads.".format(
        reads_filtered=reads_filtered,
        reads_total=reads_total,
    ))


def main():
    filter_reads(
        quality_threshold=snakemake.params.quality_threshold, 
        csfasta_file_name=snakemake.input.csfasta, 
        qual_file_name=snakemake.input.qual, 
        out_csfasta_file=snakemake.output.csfasta, 
        out_qual_file=snakemake.output.qual,
    )

try:
    main()
except Exception as e:
    import traceback
    log(traceback.format_exc())
    exit(1)

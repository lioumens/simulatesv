#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
simulatesv simulates Structural Variations and SNPs in artificial genomes.

Copyright © 2017 Michael Liou

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
from operator import itemgetter
from random import choice, betavariate, seed
import argparse as ap
import math
import multiprocessing as mp
import re
import numpy as np

# Parse arguments into global namespace
PARSER = ap.ArgumentParser(
    prog="simulate.py",
    description="Simulate DNA template sequences with known SNPs and structural variations.")

PARSER.add_argument(
    "-n",
    "--number",
    default=3,
    metavar="N",
    type=int,
    help="Number of simulated genome sequences to generate. (default=3)")
PARSER.add_argument(
    "-s",
    "--seed",
    metavar="SEED",
    type=int,
    help="Pseudorandom number generator seed. Set number to reproduce genomes and changes.")

GENOME_GROUP = PARSER.add_mutually_exclusive_group()
GENOME_GROUP.add_argument(
    "-l",
    "--length",
    default=50000,
    metavar="N",
    type=int,
    help="Total length of reference genome sequence. (default=50000)")
GENOME_GROUP.add_argument(
    "-t",
    "--template",
    default=None,
    metavar="FILE",
    help="Template to use as a base for generating SNP/SV in Fasta format.")

OUTPUT_GROUP = PARSER.add_argument_group("Output File Options")
OUTPUT_GROUP.add_argument(
    "-b",
    "--basename",
    default="sim_",
    metavar="BASE",
    help="Basename of the simulated fake sequences. (default=sim_)")
OUTPUT_GROUP.add_argument(
    "-o",
    "--output",
    default="reference.fna",
    metavar="FILE",
    help="Reference file for the mutated template (default=reference.fna)")
OUTPUT_GROUP.add_argument(
    "-c",
    "--changes",
    default="changes_",
    metavar="BASE",
    help="Basename of the changes files (default=changes_)")

SNP_GROUP = PARSER.add_argument_group("SNP Options")
SNP_GROUP.add_argument(
    "-se",
    "--snp-error-rate",
    default=.001,
    metavar="ERR",
    type=float,
    help="Error rate for SNPs (default=.001)")

INDEL_GROUP = PARSER.add_argument_group("Indel Options")
INDEL_GROUP.add_argument(
    "-ie",
    "--insertion-error-rate",
    default=.0001,
    metavar="ERR",
    type=float,
    help="Error rate for insertions (default=.0001)")
INDEL_GROUP.add_argument(
    "-de",
    "--deletion-error-rate",
    default=.0001,
    metavar="ERR",
    type=float,
    help="Error rate for deletions (default=.0001)")
INDEL_GROUP.add_argument(
    "-lis",
    "--largest-insertion-size",
    default=500,
    metavar="N",
    type=int,
    help="Largest insertion size (default=500)")
INDEL_GROUP.add_argument(
    "-lds",
    "--largest-deletion-size",
    default=500,
    metavar="N",
    type=int,
    help="Largest deletion size (default=500)")
INDEL_GROUP.add_argument(
    "-sis",
    "--smallest-insertion-size",
    default=1,
    metavar="N",
    type=int,
    help="Smallest insertion size (default=1)")
INDEL_GROUP.add_argument(
    "-sds",
    "--smallest-deletion-size",
    default=1,
    metavar="N",
    type=int,
    help="Smallest deletion size(default=1)")

TRANS_GROUP = PARSER.add_argument_group("Trans Options")
TRANS_GROUP.add_argument(
    "-te",
    "--trans-error-rate",
    default=.0001,
    metavar="ERR",
    type=float,
    help="Error rate for translocations (default=.0001)")
TRANS_GROUP.add_argument(
    "-lts",
    "--largest-trans-size",
    default=500,
    metavar="N",
    type=int,
    help="Largest translocation size (default=500)")
TRANS_GROUP.add_argument(
    "-sts",
    "--smallest-trans-size",
    default=10,
    metavar="N",
    type=int,
    help="Smallest translocation size (default=10)")

CNV_GROUP = PARSER.add_argument_group("CNV Options")
CNV_GROUP.add_argument(
    "-ce",
    "--cnv-error-rate",
    default=.0001,
    metavar="ERR",
    type=float,
    help="Error rate for CNVs (default=.0001)")
CNV_GROUP.add_argument(
    "-lcn",
    "--largest-cnv-number",
    default=100,
    metavar="N",
    type=int,
    help="Maximum number of times the CNV can repeat")
CNV_GROUP.add_argument(
    "-scn",
    "--smallest-cnv-number",
    default=2,
    metavar="N",
    type=int,
    help="Minimum number of times the CNV will repeat")
CNV_GROUP.add_argument(
    "-lcs",
    "--largest-cnv-size",
    default=50,
    metavar="N",
    type=int,
    help="Largest size for CNVs (default=300)")
CNV_GROUP.add_argument(
    "-scs",
    "--smallest-cnv-size",
    default=2,
    metavar="N",
    type=int,
    help="Smallest size for CNVs (default=2)")

ARGS = PARSER.parse_args()


def generate_dna(num_bp):
    """Generate a DNA template of length num_bp"""
    bases = "GATC"
    template_array = map(
        lambda x: bases[x], np.floor(np.random.random(num_bp) * 4).astype(int))
    return list(template_array)


def write_fasta(filename, sequence, sequence_name):
    """Write a sequence to FASTA file format"""
    header_line = ">" + sequence_name + "\n"
    # Convert sequence to str
    if isinstance(sequence, list):
        sequence = "".join(sequence)
    sequence = [sequence[n:n+80] for n in range(0, len(sequence), 80)]
    with open(filename, "w") as outfile:
        outfile.write(header_line)
        for line in sequence:
            outfile.write(line + '\n')


def validate_template(template):
    """Check that the provided template only consists of DNA nucleotides"""
    for base in template:
        if base not in "GATC":
            raise ValueError("template must consist of only GATC")


def output_template(seedval=None):
    """Either generate a template or use the one provided"""
    if seedval:
        np.random.seed(seedval)

    if ARGS.template:
        template, header = get_first_fasta(ARGS.template)
    else:
        template = generate_dna(ARGS.length)
        header = "Simulated Reference Genome"
    write_fasta(ARGS.output, template, header)
    return template


def mutate_template(template):
    """
    Mutate a given template with SNPs and SVs. The mutations work in three
    stages: all deletions (trans, del), snps, then all insertions (cnvs, ins,
    trans). The reason for this logical split is so that none of the SVs or
    SNPs occur within an inserted element. This would make tracking the changes
    very difficult. The list of changes are lists with the format

    [SV, ref_idx, alt_idx, size, ref, alt]

    The alt_idx is only for translocations to track where the translocated
    element moved in the simulated genome.
    """
    changes = []
    template, trans_queue, changes = _mutate_deletions(template, changes)
    template, changes = _mutate_snps(template, trans_queue, changes)
    template, changes = _mutate_insertions(template, trans_queue, changes)
    return template, changes


def _mutate_deletions(template, changes):
    """Returns template with all deletions. (Translocations and Deletions)"""
    def get_deletions_list():
        """Helper function to provide all deletions in sorted list"""
        num_deletions = np.random.binomial(len(template), ARGS.deletion_error_rate)
        deletions_pos = np.random.rand(num_deletions) * len(template)
        deletions_pos = np.unique(np.floor(deletions_pos)).astype(int)
        num_trans = np.random.binomial(len(template), ARGS.trans_error_rate)
        trans_pos = np.random.rand(num_trans) * len(template)
        trans_pos = np.unique(np.floor(trans_pos)).astype(int)

        # Combine and sort the locations of deletions
        del_array = [("DEL", x) for x in deletions_pos]
        trans_array = [("TRANS", x) for x in trans_pos]

        deletion_list = _merge_sorted(del_array, trans_array)
        return deletion_list

    deletion_list = get_deletions_list()

    trans_queue = []
    # Reference position changes as deletions occur
    deletion_offset = 0
    for sv_type, pos in deletion_list:
        pos = int(pos)
        if pos >= len(template):
            # deletions are already past end of genome
            break
        elif sv_type == "DEL":
            deletion_size = translate_random_number(
                betavariate(.5, 2), ARGS.smallest_deletion_size, ARGS.largest_deletion_size)
        else:
            deletion_size = translate_random_number(
                betavariate(.5, 2), ARGS.smallest_trans_size, ARGS.largest_trans_size)
            # 3rd entry "pos" is saved solely for comparison for getting deletion offset
            trans_queue.append(
                [pos + deletion_offset + 1, template[pos:pos + deletion_size], pos + 1])

        # If index extends past the end, the actual deletion will be shorter
        dna_fragment = "".join(template[pos:pos + deletion_size])
        real_deletion_size = len(dna_fragment)
        template = template[:pos] + template[pos + deletion_size:]

        if sv_type == "DEL":
            ref_idx = pos + deletion_offset + 1
            # Reference index starts from 1
            entry = [sv_type, ref_idx, pos + 1, real_deletion_size, dna_fragment, "."]
            changes.append(entry)

        deletion_offset += real_deletion_size

    return template, trans_queue, changes


def _mutate_snps(template, trans_queue, changes):
    """Returns template with snps"""
    new_template = template
    num_snps = np.random.binomial(len(template), ARGS.snp_error_rate)
    snps_pos = np.random.rand(num_snps) * len(template)
    snps_pos = np.unique(np.floor(snps_pos)).astype(int)

    for pos in snps_pos:
        original_base = new_template[pos]
        new_template[pos] = mutate_single_base(new_template[pos])
        # Reference index starts from 1
        deletion_offset = get_deletion_offset(pos, trans_queue, changes)
        entry = ["SNP", pos + deletion_offset + 1, pos + 1, "1", original_base, new_template[pos]]
        changes.append(entry)
    return new_template, changes


def _merge_sorted(*arrays):
    """Helper function to merge k sorted lists of structural variations by position number"""
    if not arrays:
        return []
    elif len(arrays) == 1:
        return arrays[0]

    array_a = arrays[0]
    array_b = arrays[1]

    result = []
    i = 0
    j = 0
    while i < len(array_a) or j < len(array_b):
        if i >= len(array_a):
            result.append(array_b[j])
            j += 1
        elif j >= len(array_b):
            result.append(array_a[i])
            i += 1
        elif array_a[i][1] < array_b[j][1]:
            result.append(array_a[i])
            i += 1
        elif array_a[i][1] > array_b[j][1]:
            result.append(array_b[j])
            j += 1
        else:
            result.append(array_a[i])
            i += 1
            j += 1
    return _merge_sorted(result, *arrays[2:])


def get_deletion_offset(alt_pos, trans_queue, changes):
    """Helper function to calculate original position in reference genome"""
    offset = 0
    for struct_var in changes:
        if struct_var[2] > alt_pos and struct_var[0] != "TRANS":
            break
        elif struct_var[0] == "DEL":
            offset += struct_var[3]
    for trans in trans_queue:
        if trans[2] > alt_pos:
            return offset
        offset += len(trans[1])
    return offset


def _mutate_insertions(template, trans_queue, changes):
    """Returns temmplate with insertions (Translocations, CNVs, and Insertions)"""
    def find_ref_index(pos):
        """Helper function to get insert changes in sorted order by reference position"""
        for i, struct_var in enumerate(changes):
            if struct_var[1] > pos:
                return i
        return len(changes)

    def get_insert_list():
        """Helper function to generate one sorted list of all insertions"""
        num_insertions = np.random.binomial(
            len(template), ARGS.insertion_error_rate)
        insertion_pos = np.random.rand(num_insertions) * len(template)
        insertion_pos = np.unique(np.floor(insertion_pos)).astype(int)
        num_cnv = np.random.binomial(len(template), ARGS.cnv_error_rate)
        cnv_pos = np.random.rand(num_cnv) * len(template)
        cnv_pos = np.unique(np.floor(cnv_pos)).astype(int)
        trans_pos = np.random.rand(len(trans_queue)) * len(template)
        trans_pos = np.unique(np.floor(trans_pos)).astype(int)
        trans_data = zip(trans_pos, trans_queue)

        insert_array = [("INS", x) for x in insertion_pos]
        # trans[0] is the original position, trans[1] is the fragment deleted
        trans_array = [("TRANS", x, trans[0], trans[1]) for x, trans in trans_data]
        cnv_array = [("CNV", x) for x in cnv_pos]

        insert_list = _merge_sorted(insert_array, cnv_array, trans_array)
        return insert_list

    def update_snp_changes(pos, insertion_length):
        """Update snp alt index as insertions are being made"""
        for i, struct_var in enumerate(changes):
            if pos < struct_var[2] and (struct_var[0] == "SNP" or struct_var[0] == "DEL"):
                changes[i][2] += insertion_length
        for i, trans in enumerate(trans_queue):
            if pos < trans[2]:
                trans_queue[i][2] += insertion_length

    changes = sorted(changes, key=itemgetter(1))
    insert_list = get_insert_list()

    insertion_offset = 0
    for struct_var in insert_list:
        pos = struct_var[1]
        alt_pos = pos + insertion_offset + 1
        ref_pos = pos + get_deletion_offset(alt_pos, trans_queue, changes) + 1
        if struct_var[0] == "INS":
            dna_fragment = "".join(generate_insertion())
            entry = ["INS", ref_pos, alt_pos, len(dna_fragment), ".", dna_fragment]
        elif struct_var[0] == "CNV":
            dna_fragment = "".join(generate_cnv())
            entry = ["CNV", ref_pos, alt_pos, len(dna_fragment), ".", dna_fragment]
        else:
            dna_fragment = "".join(struct_var[3])
            ref_pos = struct_var[2]
            alt_pos = struct_var[1] + insertion_offset + 1
            entry = ["TRANS", ref_pos, alt_pos, len(dna_fragment), dna_fragment, dna_fragment]
        changes.insert(find_ref_index(ref_pos), entry)

        template[pos + insertion_offset:pos + insertion_offset] = dna_fragment
        update_snp_changes(alt_pos, len(dna_fragment))
        insertion_offset += len(dna_fragment)
    return template, changes


def generate_insertion():
    """Generate a random insertion sequence within the script arguments"""
    sv_size = translate_random_number(
        betavariate(.5, 2), ARGS.smallest_insertion_size, ARGS.largest_insertion_size)
    return generate_dna(sv_size)


def generate_cnv():
    """Generate a random CNV DNA sequence within the script arguments"""
    sv_size = translate_random_number(
        betavariate(.5, 2), ARGS.smallest_cnv_size, ARGS.largest_cnv_size)
    num_cnv = translate_random_number(
        betavariate(.5, 2), ARGS.smallest_cnv_number, ARGS.largest_cnv_number)
    return generate_dna(sv_size) * num_cnv


def write_changes(changes_file, changes):
    """Write SVs and SNPs to a file"""
    with open(changes_file, "w") as outfile:
        header = "\t".join(["type", "ref_idx", "alt_idx",
                            "size", "ref", "alt", "\n"])
        outfile.write(header)
        for change in changes:
            change = "\t".join(map(str, change)) + "\n"
            outfile.write(change)


def translate_random_number(number, start, stop):
    """
    Helper function to translate a random number 0 <= x <= 1 to integer
    between start <= x <= stop
    """
    return int(math.floor(number * (stop - start) + start))


def mutate_single_base(base):
    """Takes single nucleic acid and changes it to a different nucleic acid."""
    bases = "GATC"
    idx = bases.index(base)
    return choice(bases[:idx] + bases[idx + 1:])


def get_first_fasta(filename):
    """Validates the given file as a proper FASTA file. Returns the first sequence as string."""
    with open(filename) as file_in:
        # Recognize FASTA files with > or ; as first letter
        if file_in.read(1) not in ">;":
            raise TypeError(
                "{0} is invalid FASTA file, must start with '>' or ';'".format(filename))
        header = file_in.readline()
        seq = []
        for num, line in enumerate(file_in, 2):
            line = line.strip()
            # regular expression checking each line is made of AGCT's
            if re.fullmatch("^[AGCT]+", line):
                seq.append(line.strip())
            # empty line indicating end of first record or end of file
            elif not line:
                break
            else:
                raise ValueError(
                    "{0}:{1} sequence in FASTA file has invalid character".format(
                        filename, num))
    return list("".join(seq)), header.strip()


def generate_mutated_template(template, outfile, changes_file, header, seedval=None):
    """Generate a mutated template and writes fasta and changes file"""
    if seedval:
        np.random.seed(seedval)
        seed(seedval)
    else:
        np.random.seed()
        seed()
    mutated_template, changes = mutate_template(template)
    write_fasta(outfile, mutated_template, header)
    write_changes(changes_file, changes)
    return


def main():
    """Main function for script"""
    template = output_template(ARGS.seed)

    # Setup multiprocessing jobs, one for each genome to be generated
    procs = []
    for i in range(ARGS.number):
        outfile = "{0}{1}.fna".format(ARGS.basename, str(i))
        changes_file = "{0}{1}.txt".format(ARGS.changes, str(i))
        header = "Simulated_{0}".format(str(i))
        seedval = ARGS.seed + i if ARGS.seed else ARGS.seed
        procs.append(
            mp.Process(
                target=generate_mutated_template,
                args=(
                    template,
                    outfile,
                    changes_file,
                    header,
                    seedval)))
    for process in procs:
        process.start()
    for process in procs:
        process.join()


if __name__ == "__main__":
    main()

#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

# LIBRARY IMPORTS
import argparse
import sys
import gzip
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, List
from loguru import logger
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw


# METADATAS
__author__ = "Essmay Touami"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Essmay Touami"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Essmay Touami"
__email__ = "essmay.touami@etu.u-paris.fr"
__status__ = "Developpement"


# FUNCTIONS
def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    seq = ""
    for line in gzip.open(amplicon_file, 'rt'):
        if line.startswith(">"): # it's a header
            if len(seq) >= minseqlen:
                yield seq
            # We start a new sequence
            seq = ""
        else:
            seq += line.strip()
    if len(seq) >= minseqlen:  # Yield the last sequence
        yield seq


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequences and return unique sequences sorted by count.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum count for the sequence to be considered.
    :return: A generator object that provides a [sequence, count] list with sequences
             having a count >= mincount, sorted by occurrence in descending order.
    """
    # Count the number of sequences in the file
    seq_counts = Counter(read_fasta(amplicon_file, minseqlen))
    # Sort the sequences by count in descending order
    sorted_sequences = sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)
    # Yield only the sequences with counts >= mincount
    for seq, count in sorted_sequences:
        if count >= mincount:
            yield [seq, count]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences
                            in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The percentage of identity between the two sequences.
    """
    # Get sequences aligned
    seq1, seq2 = alignment_list
    # Compute matches between the two sequences by ignoring gaps
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')
    # Compute the identity rate
    identity = matches / len(seq1) * 100

    return identity


def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int,
                                chunk_size: int = 0, kmer_size: int = 0) -> List:
    """Compute an abundance greedy clustering
        regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    logger.info("Running abundance greedy clustering...")
    otu_list = []
    # Iterate over the dereplicated sequences
    for sequence, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        is_new_otu = True
        for otu in otu_list:
            # Align the sequence with the OTU
            alignment = nw.global_align(sequence, otu[0], gap_open=-1, gap_extend=-1,
                                        matrix=str(Path(__file__).parent / "MATCH"))
            # Compute the identity rate
            identity = get_identity(alignment)
            # If the identity rate is greater than 97%, the sequence is not a new OTU
            if identity >= 97:
                is_new_otu = False
                break
        if is_new_otu:
            otu_list.append([sequence, count])

    logger.success(f"{len(otu_list)} OTUs have been identified successfully ! \n")
    return otu_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    logger.info("Saving the OTU sequences in fasta format...")
    # create results folder if it does not exist
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w", encoding= "utf-8") as file_out:
        for i, (seq, count) in enumerate(OTU_list):
            file_out.write(f">OTU_{i+1} occurrence:{count}\n")

            # Write the sequence in fasta format
            wrapped_seq = textwrap.fill(seq, width=80)
            file_out.write(f"{wrapped_seq}\n")

    logger.success(f"Saving the OTU sequences in {output_file} sucessfully ! \n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Compute the OTU sequences
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount)

    # Write the OTU sequences in a fasta file
    write_OTU(OTU_list, args.output_file)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Simple demultiplexer: read R1 FASTQ and index FASTQ in lockstep and route reads
to per-index gzipped FASTQ files.

Creates output file per index sequence and an unmatched file for reads that do not
match any index within the allowed Hamming distance.
"""

import gzip
import sys
import os
from Bio import SeqIO
from pathlib import Path
import logging
import argparse

_log = logging.getLogger("demultiplex")

def set_logging(log_file, log_level):
    _log.setLevel(log_level)
    # create file handler that logs debug and higher level messages
    fh = logging.FileHandler(log_file)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add the handlers to logger
    _log.addHandler(ch)
    _log.addHandler(fh)

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two sequences."""
    if len(s1) != len(s2):
        raise ValueError("Sequences must be of equal length to calculate Hamming distance.")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def demultiplex(r1_path, index_path, out_dir, index_seqs, max_hamming):

    index_seqs = [seq.upper() for seq in index_seqs.split(",")]
    _log.info(f"Demultiplexing with index sequences: {index_seqs}")
    _log.info(f"Maximum Hamming distance allowed: {max_hamming}")

    # derive a base name for the r1 file without any compressed/fastq suffixes
    r1_name = Path(r1_path).name
    if r1_name.endswith('.fastq.gz'):
        r1_base = r1_name[:-9]  # strip '.fastq.gz'
    elif r1_name.endswith('.fq.gz'):
        r1_base = r1_name[:-6]  # strip '.fq.gz'
    else:
        # fallback: use stem (may still include '.fastq' if double-suffixed)
        r1_base = Path(r1_path).stem

    out_paths_map = {
        seq.upper(): os.path.join(
            os.path.dirname(out_dir),
            Path(out_dir).stem + "_" + seq.upper(),
            f"{r1_base}_{seq.upper()}.fastq.gz",
        )
        for seq in index_seqs
    }

    # make sure output directories exist
    for path in out_paths_map.values():
        Path(os.path.dirname(path)).mkdir(parents=True, exist_ok=True)  

    unmatched_out_path = os.path.join(out_dir, f"{r1_base}_Unmatched.fastq.gz")

    # create output directory if it doesn't exist
    Path(os.path.dirname(unmatched_out_path)).mkdir(parents=True, exist_ok=True)

    # iterate through R1 and index FASTQ in lockstep
    with gzip.open(r1_path, "rt") as r1_handle, \
        gzip.open(index_path, "rt") as index_handle:

        try:
            # prepare output handles
            unmatched_handle = gzip.open(unmatched_out_path, "wt")
            out_handles = {seq: gzip.open(path, "wt") for seq, path in out_paths_map.items()}

            r1_iter = SeqIO.parse(r1_handle, "fastq")
            index_iter = SeqIO.parse(index_handle, "fastq")

            total_reads = 0
            matched_reads = {seq.upper(): 0 for seq in index_seqs}
            unmatched_reads = 0
            for r1_record, index_record in zip(r1_iter, index_iter, strict=True):
                total_reads += 1
                if total_reads % 100_000 == 0:
                    _log.info(f"Processed {total_reads} reads")
                index_seq = str(index_record.seq)

                best_match = None
                best_distance = max_hamming + 1  # initialize to max_hamming + 1
                for target_seq in index_seqs:
                    dist = hamming_distance(index_seq.upper(), target_seq)
                    if dist < best_distance:
                        best_distance = dist
                        best_match = target_seq

                if best_match is not None:
                    SeqIO.write(r1_record, out_handles[best_match], "fastq")
                    matched_reads[best_match] += 1
                else:
                    SeqIO.write(r1_record, unmatched_handle, "fastq")
                    unmatched_reads += 1
        finally:
            # close all output handles
            for handle in out_handles.values():
                handle.close()
            unmatched_handle.close()

    _log.info(f"Demultiplexing complete. Total reads: {total_reads}")
    for seq in index_seqs:
        _log.info(f"Reads matched to index {seq}: {matched_reads[seq]}")
        _log.info(f"Reads saved to file: {out_paths_map[seq]}")
    _log.info(f"Unmatched reads: {unmatched_reads}")
    _log.info(f"Unmatched reads saved to file: {unmatched_out_path}")

def main():
    parser = argparse.ArgumentParser(description="Demultiplex FASTQ files")
    parser.add_argument("--r1", required=True, help="Path to R1 FASTQ file")
    parser.add_argument("--index", required=True, help="Path to index FASTQ file")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--index-seqs", required=True, help="Comma-separated list of index sequences")
    parser.add_argument("--max-hamming", type=int, default=0, help="Maximum Hamming distance for index matching")
    parser.add_argument("--log", help="Path to log file", default="demultiplex.log")
    parser.add_argument("--log-level", default="INFO", help="Log level")
    args = parser.parse_args()

    set_logging(args.log, args.log_level)
    demultiplex(args.r1, args.index, args.out_dir, args.index_seqs, args.max_hamming)

if __name__ == "__main__":
    # handle both script and snakemake execution
    if 'snakemake' in globals():
        # snakemake execution
        args = [
            "--r1", snakemake.input.r1,
            "--index", snakemake.input.index,
            "--out-dir", snakemake.output[0],
            "--index-seqs", ",".join(snakemake.params.index_seqs),
        ]
        if hasattr(snakemake.params, 'max_hamming'):
            args += ["--max-hamming", str(snakemake.params.max_hamming)]
        if hasattr(snakemake, 'log') and snakemake.log:
            args += ["--log", snakemake.log[0]]
        if hasattr(snakemake.params, 'log_level'):
            args += ["--log-level", snakemake.params.log_level]
        sys.argv[1:] = args
        main()
    else:
        # script execution
        main()
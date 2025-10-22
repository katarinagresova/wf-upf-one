#!/usr/bin/env python3
"""
Simple demultiplexer: read R1 FASTQ and index FASTQ in lockstep and route reads
to per-index gzipped FASTQ files defined in snakemake.params.index_seqs.

Expects to be run as a Snakemake script. Uses snakemake.input, snakemake.output,
snakemake.params and snakemake.log.
"""
import gzip
import sys
import os
from Bio import SeqIO
from pathlib import Path

# Snakemake context
r1_path = snakemake.input.r1
index_path = snakemake.input.index
out_dir = snakemake.output[0]
index_seqs = snakemake.params.index_seqs
log_path = snakemake.log[0] if snakemake.log else None

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
     gzip.open(index_path, "rt") as index_handle, \
     gzip.open(unmatched_out_path, "wt") as unmatched_handle:
    # prepare output handles
    out_handles = {seq: gzip.open(path, "wt") for seq, path in out_paths_map.items()}

    r1_iter = SeqIO.parse(r1_handle, "fastq")
    index_iter = SeqIO.parse(index_handle, "fastq")

    total_reads = 0
    matched_reads = 0

    for r1_record, index_record in zip(r1_iter, index_iter, strict=True):
        total_reads += 1
        index_seq = str(index_record.seq)

        if index_seq.upper() in out_paths_map.keys():
            SeqIO.write(r1_record, out_handles[index_seq.upper()], "fastq")
            matched_reads += 1
        else:
            SeqIO.write(r1_record, unmatched_handle, "fastq")

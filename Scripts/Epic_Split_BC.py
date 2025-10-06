#!/usr/bin/env python3
"""
Faster FASTQ splitter (gz input) that produces two FASTQ outputs:
 - R1: sequence from seq_start to end
 - R2: barcode made from BC2 + BC1 (as in your original code)

Usage:
    python split_fastq_fast.py -i reads_R2.fastq.gz -o1 out_R1.fastq.gz -o2 out_R2.fastq.gz
"""
from gzip import open as gzopen
import argparse
import os
import sys

def parse_args():
    ap = argparse.ArgumentParser(description="Fast FASTQ splitter (gz input).")
    ap.add_argument("-i", "--input", required=True, help="input gzipped FASTQ file")
    ap.add_argument("-o1", "--output_R1", required=True, help="output file R1 (can be .gz)")
    ap.add_argument("-o2", "--output_R2", required=True, help="output file R2 (can be .gz)")
    ap.add_argument("--seq_start", type=int, default=117, help="start index for R1 seq/qual (0-based)")
    ap.add_argument("--bc2_start", type=int, default=22, help="BC2 start (0-based, inclusive)")
    ap.add_argument("--bc2_end", type=int, default=30, help="BC2 end (0-based, exclusive)")
    ap.add_argument("--bc1_start", type=int, default=60, help="BC1 start (0-based, inclusive)")
    ap.add_argument("--bc1_end", type=int, default=68, help="BC1 end (0-based, exclusive)")
    return ap.parse_args()

def open_out(path):
    # If filename ends with .gz, open gzip in binary write mode, else plain binary file
    if path.endswith(".gz"):
        return gzopen(path, "wb")
    else:
        # Use buffering to reduce syscalls: buffering=1<<20 (~1MiB)
        return open(path, "wb", buffering=1 << 20)

def main():
    args = parse_args()

    # local copies for speed
    seq_start = args.seq_start
    bc2_start = args.bc2_start
    bc2_end = args.bc2_end
    bc1_start = args.bc1_start
    bc1_end = args.bc1_end

    # Sanity: warn if start > end (will produce empty slices)
    if (bc2_start > bc2_end) or (bc1_start > bc1_end):
        sys.stderr.write(
            "Warning: one of bc start > end; slices will be empty. "
            "Check bc2_start/bc2_end and bc1_start/bc1_end.\n"
        )

    # open files (binary)
    with gzopen(args.input, "rb") as fh_in, open_out(args.output_R1) as fh_out1, open_out(args.output_R2) as fh_out2:
        write1 = fh_out1.write
        write2 = fh_out2.write

        # process 4 lines at a time
        rec = 0
        while True:
            header = fh_in.readline()
            if not header:
                break  # EOF
            seq = fh_in.readline()
            plus = fh_in.readline()
            qual = fh_in.readline()

            # Basic FASTQ validation (if file truncated, break)
            if not qual:
                # bad/incomplete record; stop
                break

            # strip trailing newline bytes (\n or \r\n)
            header = header.rstrip(b'\r\n')
            seq = seq.rstrip(b'\r\n')
            plus = plus.rstrip(b'\r\n')
            qual = qual.rstrip(b'\r\n')

            # Extract title without leading '@' to mirror FastqGeneralIterator behaviour
            # (but keep the rest of header text)
            if header.startswith(b'@'):
                title = header[1:]
            else:
                title = header

            # slices (bytes)
            new_seq_R1 = seq[seq_start:]
            new_qual_R1 = qual[seq_start:]
            barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end]  # BC2 + BC1
            new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]

            # write outputs as standard FASTQ (binary)
            # minimize number of writes by concatenating the record bytes
            rec1 = b'@' + title + b'\n' + new_seq_R1 + b'\n+\n' + new_qual_R1 + b'\n'
            rec2 = b'@' + title + b'\n' + barcode + b'\n+\n' + new_qual_R2 + b'\n'
            write1(rec1)
            write2(rec2)

            rec += 1
            # optional: progress print every N records (uncomment to enable)
            # if rec % 100000 == 0:
            #     sys.stderr.write(f"Processed {rec} records\n")

    # done
    sys.stderr.write(f"Finished. Processed {rec} records.\n")

if __name__ == "__main__":
    main()
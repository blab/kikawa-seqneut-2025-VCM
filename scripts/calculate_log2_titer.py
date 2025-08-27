#!/usr/bin/env python3
import argparse
import csv
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--titers", required=True, help="TSV file of titers to calculate log2 value for")
    parser.add_argument("--titer-column", default="titer", help="name of titer column in the input titers table")
    parser.add_argument("--log2-titer-column", default="log2_titer", help="name of the log2 titer column in the output table")
    parser.add_argument("--output", required=True, help="TSV file of titers with log2 value added")

    args = parser.parse_args()

    with open(args.titers, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames + [args.log2_titer_column]

        with open(args.output, "w", encoding="utf-8") as oh:
            writer = csv.DictWriter(oh, fieldnames=header, delimiter="\t", lineterminator="\n")
            writer.writeheader()

            for row in reader:
                row[args.log2_titer_column] = np.log2(float(row[args.titer_column]))
                writer.writerow(row)

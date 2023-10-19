#!/usr/bin/env python

"""Script to merge Bracken output."""

from os import path
from sys import argv, exit
import pandas as pd

if __name__ == "__main__":
    files = argv[1:]

    data = []
    for fi in files:
        try:
            ab = pd.read_csv(fi, sep="\t")
        except Exception:
            print(f"Could not parse {fi}. Is that a valid bracken file?")
        ab["sample_id"] = path.basename(fi).split("_kraken2")[0]
        data.append(ab)
    merged = pd.concat(data)
    merged.to_csv("bracken_merged.csv", index=False)
    print(f"Merged {len(data)} files into `bracken_merged.csv`.")

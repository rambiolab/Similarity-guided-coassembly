#!/usr/bin/env python3
"""Convert a sourmash compare similarity matrix (.npy + .npy.labels.txt) into a
distance matrix TSV: distance = 1 - similarity (Jaccard, or angular similarity
for abundance sketches), diagonal set to 0.

Usage: sourmash_npy_to_tsv.py matrix.npy distance_matrix.tsv
"""

import os
import sys

import numpy as np
import pandas as pd

npy_file, out_file = sys.argv[1:3]

sim = np.load(npy_file)
dist = 1 - sim
np.fill_diagonal(dist, 0)

with open(npy_file + ".labels.txt") as f:
    samples = [os.path.basename(line.strip()).replace(".sig", "") for line in f if line.strip()]

if len(samples) != dist.shape[0]:
    sys.exit(f"{len(samples)} labels for a {dist.shape[0]} x {dist.shape[1]} matrix")

df = pd.DataFrame(dist, index=samples, columns=samples)
df.index.name = "Sample"
df.to_csv(out_file, sep="\t")

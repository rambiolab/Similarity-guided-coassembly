#!/usr/bin/env python3
"""Convert `mash dist` output (query, reference, distance, p-value, shared hashes)
into a square distance matrix TSV. Sample names are the file names without
directory and .fq.gz; samples are sorted by name; diagonal set to 0.

Usage: mash_dist_to_matrix.py mash_distances.tab distance_matrix.tsv
"""

import os
import re
import sys

import numpy as np
import pandas as pd

tab_file, out_file = sys.argv[1:3]


def clean(name):
    return re.sub(r"\.fq\.gz$", "", os.path.basename(name))


tab = pd.read_csv(tab_file, sep="\t", header=None, usecols=[0, 1, 2],
                  names=["a", "b", "dist"], dtype={"a": str, "b": str, "dist": float})
tab["a"] = tab["a"].map(clean)
tab["b"] = tab["b"].map(clean)

samples = sorted(set(tab["a"]) | set(tab["b"]))
mat = tab.pivot_table(index="a", columns="b", values="dist", aggfunc="first")
mat = mat.reindex(index=samples, columns=samples)
# fill any one-sided pair from its mirror
mat = mat.combine_first(mat.T).reindex(index=samples, columns=samples)
values = mat.to_numpy(copy=True)
np.fill_diagonal(values, 0)
mat = pd.DataFrame(values, index=samples, columns=samples)

n_missing = int(mat.isna().sum().sum())
if n_missing:
    print(f"warning: {n_missing} missing distances written as NA", file=sys.stderr)

mat.index.name = "Sample"
mat.columns.name = None
mat.to_csv(out_file, sep="\t", na_rep="NA")
print(f"{len(samples)} samples written to {out_file}", file=sys.stderr)

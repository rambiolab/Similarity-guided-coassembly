#!/usr/bin/env python3
"""Collect all MetaBAT2 bins of a partition for dRep.

Each bin <bins>/<cluster>/bin.N.fa is symlinked as <linkdir>/<cluster>_bin.N.fa,
so genome names are unique across clusters. Writes
  --bin-list  one path per line (dRep -g)
  --info      dRep --genomeInfo CSV: genome, completeness, contamination
  --quality   all bins with their CheckM2 values (TSV)
"""

import argparse
import os
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--bins", required=True, help="folder with one subfolder of bins per cluster")
ap.add_argument("--checkm", required=True, help="folder with <cluster>/quality_report.tsv")
ap.add_argument("--linkdir", required=True)
ap.add_argument("--bin-list", required=True)
ap.add_argument("--info", required=True)
ap.add_argument("--quality", required=True)
ap.add_argument("--clusters", nargs="+", required=True)
args = ap.parse_args()

os.makedirs(args.linkdir, exist_ok=True)
rows = []
for cluster in args.clusters:
    report = os.path.join(args.checkm, cluster, "quality_report.tsv")
    q = pd.read_csv(report, sep="\t")
    if q.empty:
        print(f"{cluster}: no bins", file=sys.stderr)
        continue
    for _, r in q.iterrows():
        src = os.path.abspath(os.path.join(args.bins, cluster, f"{r['Name']}.fa"))
        if not os.path.exists(src):
            sys.exit(f"Bin listed by CheckM2 not found: {src}")
        genome = f"{cluster}_{r['Name']}.fa"
        link = os.path.join(args.linkdir, genome)
        if os.path.lexists(link):
            os.remove(link)
        os.symlink(src, link)
        row = {"genome": genome, "cluster": cluster, "bin": r["Name"],
               "completeness": r["Completeness"], "contamination": r["Contamination"]}
        for col in ("Contig_N50", "Genome_Size", "Total_Contigs"):
            if col in q.columns:
                row[col.lower()] = r[col]
        rows.append(row)

if not rows:
    sys.exit("No bins in any cluster: nothing to dereplicate.")

df = pd.DataFrame(rows)
with open(args.bin_list, "w") as fh:
    for g in df["genome"]:
        fh.write(os.path.abspath(os.path.join(args.linkdir, g)) + "\n")
df[["genome", "completeness", "contamination"]].to_csv(args.info, index=False)
df.to_csv(args.quality, sep="\t", index=False)
print(f"{len(df)} bins from {df['cluster'].nunique()} clusters", file=sys.stderr)

#!/usr/bin/env python3
"""Summary tables for one partition.

--genomes-out: one row per dereplicated genome (= species at the dRep secondary
  ANI threshold): cluster of origin, CheckM2 quality, GTDB-Tk classification,
  number of bins in its species cluster.
--summary-out: one row for the partition: clusters, bins, dereplicated
  genomes, HQ (>= 90% completeness, < 5% contamination),
  MQ (>= 50%, < 10%, not HQ), and novel species (no GTDB species name).
"""

import argparse
import glob
import os

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--quality", required=True)
ap.add_argument("--cdb", required=True)
ap.add_argument("--wdb", required=True)
ap.add_argument("--gtdbtk-dir", required=True)
ap.add_argument("--partition", required=True)
ap.add_argument("--genomes-out", required=True)
ap.add_argument("--summary-out", required=True)
args = ap.parse_args()


def quality(comp, cont):
    if comp >= 90 and cont < 5:
        return "HQ"
    if comp >= 50 and cont < 10:
        return "MQ"
    return "LQ"


q = pd.read_csv(args.quality, sep="\t")
q["quality"] = [quality(a, b) for a, b in zip(q["completeness"], q["contamination"])]

cdb = pd.read_csv(args.cdb)
wdb = pd.read_csv(args.wdb)
cluster_size = cdb.groupby("secondary_cluster").size().rename("bins_in_species")

g = wdb[["genome", "cluster"]].rename(columns={"cluster": "species_cluster"})
g = g.merge(cluster_size, left_on="species_cluster", right_index=True, how="left")
g = g.merge(q.rename(columns={"cluster": "assembly_cluster"}), on="genome", how="left")

tax = []
for f in sorted(glob.glob(os.path.join(args.gtdbtk_dir, "gtdbtk.*.summary.tsv"))):
    t = pd.read_csv(f, sep="\t", usecols=["user_genome", "classification"])
    tax.append(t)
if tax:
    tax = pd.concat(tax, ignore_index=True)
    tax["genome"] = tax["user_genome"].astype(str) + ".fa"
    g = g.merge(tax[["genome", "classification"]], on="genome", how="left")
else:
    g["classification"] = pd.NA
species = g["classification"].fillna("").str.extract(r"s__([^;]*)$")[0].fillna("")
g["novel_species"] = species.str.strip() == ""

g = g.sort_values(["assembly_cluster", "genome"])
g.to_csv(args.genomes_out, sep="\t", index=False)

part = pd.read_csv(args.partition, sep="\t", dtype=str)
summary = pd.DataFrame([{
    "partition": args.partition,
    "samples": len(part),
    "clusters": part["cluster"].nunique(),
    "clusters_with_bins": q["cluster"].nunique(),
    "bins": len(q),
    "bins_MQ_or_better": int((q["quality"] != "LQ").sum()),
    "bins_HQ": int((q["quality"] == "HQ").sum()),
    "species": len(g),
    "species_HQ": int((g["quality"] == "HQ").sum()),
    "species_MQ": int((g["quality"] == "MQ").sum()),
    "species_novel": int(g["novel_species"].sum()),
}])
summary.to_csv(args.summary_out, sep="\t", index=False)
print(summary.T.to_string(header=False))

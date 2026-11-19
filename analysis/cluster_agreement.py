#!/usr/bin/env python3
"""Agreement between sample partitions and metadata labels (adjusted mutual information).

For every partition and every label column (e.g. site, city, Season), compute the AMI
(Vinh, Epps & Bailey 2010; scikit-learn, arithmetic-mean normalisation) and test it with
permutation nulls that keep the partition fixed and shuffle the labels:

  global        labels shuffled across all samples: does the partition follow the label at all?
  within <COL>  labels shuffled only within each value of COL (--within LABEL=COL, e.g.
                Season=site): does it follow the label beyond what COL explains?

AMI is chance-corrected by construction (expected 0 for random partitions); the
permutations give z-scores and empirical p-values. Also writes the AMI between all
pairs of partitions.

Example (VEO data, all runs):
  cluster_agreement.py --cluster_summary cluster_summary.tsv \
      --metadata VEO_sites.tsv --derive city=site:_:0 \
      --labels site city Season --within Season=site --outdir agreement/all_runs

Copenhagen only, one run per physical sample:
  cluster_agreement.py ... --subset city=Copenhagen --collapse_by ENA_ALIAS \
      --labels site Season --within Season=site --outdir agreement/Copenhagen_samples
"""
import argparse
import collections
import os
import random
import sys

import numpy as np
from sklearn.metrics import adjusted_mutual_info_score as AMI

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import add_metadata_args, add_partition_args, read_metadata, read_partitions, select_samples  # noqa: E402


def perm_test(clust, labels, n_perm, rng, strata=None):
    """Observed AMI, mean AMI under the permutation null, z and empirical p (one-sided)."""
    obs = AMI(labels, clust)
    groups = collections.defaultdict(list)
    for i, s in enumerate(strata or [0] * len(labels)):
        groups[s].append(i)
    null = []
    for _ in range(n_perm):
        lab = list(labels)
        for idx in groups.values():
            vals = [labels[i] for i in idx]
            rng.shuffle(vals)
            for i, v in zip(idx, vals):
                lab[i] = v
        null.append(AMI(lab, clust))
    null = np.array(null)
    sd = null.std()
    return {"AMI": obs, "null_mean": null.mean(),
            "z": (obs - null.mean()) / sd if sd > 0 else float("nan"),
            "p": (1 + (null >= obs).sum()) / (n_perm + 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_partition_args(ap)
    add_metadata_args(ap)
    ap.add_argument("--labels", nargs="+", required=True, help="metadata columns to compare with")
    ap.add_argument("--within", action="append", default=[], metavar="LABEL=COL",
                    help="also test LABEL with a null shuffled within each value of COL (repeatable)")
    ap.add_argument("--n_perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    part = read_partitions(a)
    meta = read_metadata(a)
    samples = select_samples(a, meta, part)
    within = dict(spec.split("=", 1) for spec in a.within)
    rng = random.Random(a.seed)

    rows = []
    for name, p in part.items():
        clust = [p[s] for s in samples]
        for lab in a.labels:
            labels = [meta[s][lab] for s in samples]
            rows.append(dict(partition=name, label=lab, null="global",
                             **perm_test(clust, labels, a.n_perm, rng)))
            if lab in within:
                strata = [meta[s][within[lab]] for s in samples]
                rows.append(dict(partition=name, label=lab, null=f"within_{within[lab]}",
                                 **perm_test(clust, labels, a.n_perm, rng, strata)))

    cols = ["partition", "label", "null", "AMI", "null_mean", "z", "p"]
    with open(os.path.join(a.outdir, "agreement.tsv"), "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

    names = list(part)
    with open(os.path.join(a.outdir, "partition_vs_partition.tsv"), "w") as o:
        o.write("AMI\t" + "\t".join(names) + "\n")
        for n1 in names:
            c1 = [part[n1][s] for s in samples]
            o.write(n1 + "\t" + "\t".join(
                f"{AMI(c1, [part[n2][s] for s in samples]):.3f}" for n2 in names) + "\n")

    with open(os.path.join(a.outdir, "samples_used.txt"), "w") as o:
        o.write("\n".join(samples) + "\n")
    print(f"{len(samples)} samples, {len(names)} partitions, {a.n_perm} permutations -> {a.outdir}")


if __name__ == "__main__":
    main()

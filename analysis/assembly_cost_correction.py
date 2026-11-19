#!/usr/bin/env python3
"""Correct per-cluster MEGAHIT CPU hours for resumed assemblies, and sum them per partition.

Input: a per-cluster resource table (TSV) with the columns
  cluster, input_gb, cpu_h, peak_ram_gb, resumed (True/False)
and optionally `methods` (comma-separated partition names that contain the cluster).

An assembly that was interrupted and resumed with `megahit --continue` can have a time
log that only covers the resumed part. For those clusters the CPU hours are replaced by
a prediction from a log-log fit of CPU hours on input size, made on all non-resumed
clusters with input_gb >= --min_gb:

    log(cpu_h) = a + b * log(input_gb)

The prediction's uncertainty is given as x exp(+-1 residual SD).

Note: Snakefile_assembly appends the time log on a restart (`/usr/bin/time -a`), so for
runs made with it both attempts are logged; the correction is then only needed if the
tool that builds the resource table reads just the last block of the log.

Outputs (in --outdir):
  resources_corrected.tsv   input table + cpu_h_fit, cpu_h_corrected, corrected
  fit.tsv                   intercept, slope, residual SD, n clusters in the fit
  totals.tsv                per partition: measured and corrected CPU hours, peak RAM
                            of the largest job (needs the `methods` column)
"""
import argparse
import collections
import csv
import math
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resources", required=True)
    ap.add_argument("--min_gb", type=float, default=1.0, help="smallest input size used in the fit")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows = list(csv.DictReader(open(a.resources), delimiter="\t"))
    resumed = lambda r: str(r["resumed"]).strip().lower() in ("true", "1", "yes")  # noqa: E731
    fit = [(math.log(float(r["input_gb"])), math.log(float(r["cpu_h"]))) for r in rows
           if not resumed(r) and float(r["input_gb"]) >= a.min_gb and float(r["cpu_h"]) > 0]
    x, y = np.array(fit).T
    slope, icpt = np.polyfit(x, y, 1)
    sd = float(np.std(y - (icpt + slope * x)))
    with open(os.path.join(a.outdir, "fit.tsv"), "w") as o:
        o.write("intercept\tslope\tresidual_sd\tn_clusters\tmin_gb\n")
        o.write(f"{icpt:.4f}\t{slope:.4f}\t{sd:.4f}\t{len(fit)}\t{a.min_gb}\n")

    cols = list(rows[0].keys()) + ["cpu_h_fit", "cpu_h_fit_low", "cpu_h_fit_high", "cpu_h_corrected", "corrected"]
    totals = collections.defaultdict(lambda: {"measured": 0.0, "corrected": 0.0, "ram": 0.0, "n": 0, "n_corr": 0})
    with open(os.path.join(a.outdir, "resources_corrected.tsv"), "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            gb = float(r["input_gb"])
            pred = math.exp(icpt + slope * math.log(gb)) if gb > 0 else 0.0
            corr = resumed(r)
            r.update(cpu_h_fit=f"{pred:.1f}", cpu_h_fit_low=f"{pred * math.exp(-sd):.1f}",
                     cpu_h_fit_high=f"{pred * math.exp(sd):.1f}",
                     cpu_h_corrected=f"{pred:.1f}" if corr else r["cpu_h"], corrected=str(corr))
            o.write("\t".join(str(r[c]) for c in cols) + "\n")
            for m in (r.get("methods") or "").split(","):
                if m:
                    t = totals[m]
                    t["measured"] += float(r["cpu_h"]); t["corrected"] += float(r["cpu_h_corrected"])
                    t["ram"] = max(t["ram"], float(r["peak_ram_gb"])); t["n"] += 1; t["n_corr"] += corr

    if totals:
        with open(os.path.join(a.outdir, "totals.tsv"), "w") as o:
            o.write("partition\tassemblies\tresumed_corrected\tcpu_h_measured\tcpu_h_corrected\tlargest_job_peak_ram_gb\n")
            for m, t in totals.items():
                o.write(f"{m}\t{t['n']}\t{t['n_corr']}\t{t['measured']:.0f}\t{t['corrected']:.0f}\t{t['ram']:.1f}\n")
    print(f"fit: log(cpu_h) = {icpt:.2f} + {slope:.2f} log(input_gb), residual SD {sd:.2f}, "
          f"{len(fit)} clusters; {sum(resumed(r) for r in rows)} resumed clusters corrected -> {a.outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""CPU hours per recovered species for one partition, optionally for a subset of samples
and against a reference approach (e.g. per-site co-assembly with predicted cost).

Species
  Without --combined, species are the partition's own dRep secondary clusters (Cdb.csv).
  With --combined (a dRep Cdb.csv over several approaches together), species are the
  combined clusters, so the partition and the reference are counted on one species
  definition. Each own species is mapped to the combined table through its winner genome
  (Wdb.csv). A species counts for a group of assemblies if ANY bin from those assemblies
  falls in it (all bins from the own Cdb.csv, not only the winners).

Genome names are expected as <prefix><cluster>_bin.N.fa (as written by
scripts/prepare_drep_input.py, plus an optional prefix such as `s_fkm_`), so the assembly
of origin is the part before `_bin.`. --genome_prefix / --ref_prefix are stripped to get
the names used in the combined table.

Subset
  --subset COL=VALUE keeps the clusters whose samples all have that value (e.g.
  city=Copenhagen). --mixed decides what happens with clusters that contain the subset
  only partly: fail (default), exclude, include (count in full), or majority (count in
  full if most samples are in the subset).

CPU hours
  From --resources (column cpu_h_corrected if present, e.g. from
  assembly_cost_correction.py, else cpu_h). The uncertainty of CPU-h per species is a
  bootstrap over the partition's clusters (resampled with replacement; CPU hours summed,
  species united). For the reference, --ref_cost gives cpu_h, cpu_h_low, cpu_h_high and
  peak_ram_gb per reference group (e.g. predicted per-site cost), and its range is used.

Example (VEO, Copenhagen, s_fkm vs per-site):
  cpu_per_species.py --cluster_summary cluster_summary.tsv --method s_fkm \
      --metadata VEO_sites.tsv --derive city=site:_:0 --subset city=Copenhagen \
      --resources cost/resources_corrected.tsv \
      --cdb s_fkm/drep/data_tables/Cdb.csv --wdb s_fkm/drep/data_tables/Wdb.csv \
      --genome_prefix s_fkm_ --combined combined/drep/data_tables/Cdb.csv \
      --ref_name per-site --ref_cdb paper_site/drep/data_tables/Cdb.csv --ref_prefix paper_site_ \
      --ref_groups site_Copenhagen_RA site_Copenhagen_RD site_Copenhagen_RL \
      --ref_cost site_assembly_prediction.tsv --ref_cost_id site --ref_cost_strip site_ \
      --outdir cpu_per_species/Copenhagen_s_fkm
"""
import argparse
import collections
import csv
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import add_metadata_args, add_partition_args, read_metadata, read_partitions  # noqa: E402


def strip(name, prefix):
    return name[len(prefix):] if prefix and name.startswith(prefix) else name


def origin(name):
    return name.split("_bin.")[0]


def species_by_group(cdb, wdb, prefix, combined):
    """{assembly group: set of species}. Species = combined cluster if `combined`, else own cluster."""
    rows = list(csv.DictReader(open(cdb)))
    if combined is None:
        to_sp = lambda k: k  # noqa: E731
    else:
        if wdb:
            win = {r["cluster"]: strip(r["genome"], prefix) for r in csv.DictReader(open(wdb))}
        else:  # winner = the member of each own cluster that is present in the combined table
            win = {}
            for r in rows:
                g = strip(r["genome"], prefix)
                if g in combined:
                    if r["secondary_cluster"] in win:
                        sys.exit(f"{cdb}: several members of {r['secondary_cluster']} in the combined table; give the Wdb.csv")
                    win[r["secondary_cluster"]] = g
        missing = [k for k in {r["secondary_cluster"] for r in rows} if k not in win or win[k] not in combined]
        if missing:
            sys.exit(f"{cdb}: {len(missing)} species have no winner in the combined table (check prefixes)")
        to_sp = lambda k: combined[win[k]]  # noqa: E731
    out = collections.defaultdict(set)
    for r in rows:
        out[origin(strip(r["genome"], prefix))].add(to_sp(r["secondary_cluster"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_partition_args(ap)
    add_metadata_args(ap)
    ap.add_argument("--method", required=True, help="name of the partition to evaluate")
    ap.add_argument("--mixed", choices=["fail", "exclude", "include", "majority"], default="fail")
    ap.add_argument("--resources", required=True, help="per-cluster table: cluster, input_gb, cpu_h[_corrected], peak_ram_gb")
    ap.add_argument("--cdb", required=True); ap.add_argument("--wdb", required=True)
    ap.add_argument("--genome_prefix", default="")
    ap.add_argument("--cluster_prefix", default="",
                    help="prefix that turns partition cluster ids into assembly names, e.g. cluster_ "
                         "for partitions from Snakefile_clustering (ids 1, 2, ...)")
    ap.add_argument("--combined", help="combined dRep Cdb.csv (needed for --ref_*)")
    ap.add_argument("--ref_name", default="reference")
    ap.add_argument("--ref_cdb"); ap.add_argument("--ref_wdb")
    ap.add_argument("--ref_prefix", default="")
    ap.add_argument("--ref_groups", nargs="+", help="reference assembly groups to count (default: all)")
    ap.add_argument("--ref_cost", help="TSV with cpu_h, cpu_h_low, cpu_h_high, peak_ram_gb per reference group")
    ap.add_argument("--ref_cost_id", default="group", help="id column of --ref_cost")
    ap.add_argument("--ref_cost_strip", default="", help="prefix removed from group names to match --ref_cost ids")
    ap.add_argument("--n_boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    a.collapse_by = None
    os.makedirs(a.outdir, exist_ok=True)

    part = read_partitions(a)
    if a.method not in part:
        sys.exit(f"--method {a.method} not among partitions {list(part)}")
    members = collections.defaultdict(list)
    for s, c in part[a.method].items():
        members[a.cluster_prefix + c].append(s)

    # which clusters belong to the subset
    meta = read_metadata(a)
    def in_subset(s):
        return all(meta[s][col] == val for col, val in (x.split("=", 1) for x in a.subset))
    used = sorted(c for c, ss in members.items() if all(in_subset(s) for s in ss))
    mixed = sorted(c for c, ss in members.items() if c not in used and any(in_subset(s) for s in ss))
    if mixed and a.mixed == "fail":
        sys.exit(f"{a.method} clusters only partly in the subset: {mixed}; choose --mixed")
    if a.mixed == "include":
        used = sorted(used + mixed)
    elif a.mixed == "majority":
        used = sorted(used + [c for c in mixed if 2 * sum(in_subset(s) for s in members[c]) > len(members[c])])

    res = {r["cluster"]: r for r in csv.DictReader(open(a.resources), delimiter="\t")}
    cpu_col = "cpu_h_corrected" if "cpu_h_corrected" in next(iter(res.values())) else "cpu_h"
    absent = [c for c in used if c not in res]
    if absent:
        sys.exit(f"no resource rows for clusters {absent[:5]}... (check --cluster_prefix)")
    cpu = {c: float(res[c][cpu_col]) for c in used}

    combined = {r["genome"]: r["secondary_cluster"] for r in csv.DictReader(open(a.combined))} if a.combined else None
    sp = species_by_group(a.cdb, a.wdb, a.genome_prefix, combined)
    s_all = set().union(*sp.values())
    s_used = set().union(*(sp.get(c, set()) for c in used))
    s_other = set().union(*(sp[c] for c in sp if c not in used))

    with open(os.path.join(a.outdir, "clusters.tsv"), "w") as o:
        o.write(f"cluster\tn_samples\tinput_gb\tcpu_h ({cpu_col})\tpeak_ram_gb\tspecies\tspecies_only_in_this_cluster\n")
        for c in sorted(used, key=lambda c: -len(members[c])):
            rest = set().union(*(sp[d] for d in sp if d != c))
            o.write(f"{c}\t{len(members[c])}\t{res[c]['input_gb']}\t{cpu[c]:.0f}\t{res[c]['peak_ram_gb']}\t"
                    f"{len(sp.get(c, ()))}\t{len(sp.get(c, set()) - rest)}\n")

    rng = random.Random(a.seed)
    boot = []
    for _ in range(a.n_boot):
        pick = [rng.choice(used) for _ in used]
        n = len(set().union(*(sp.get(c, set()) for c in pick)))
        boot.append(sum(cpu[c] for c in pick) / max(n, 1))
    lo, hi = np.quantile(boot, [0.025, 0.975])
    total = sum(cpu.values())
    ram = max(float(res[c]["peak_ram_gb"]) for c in used)

    summary = [f"{a.method} (mixed={a.mixed})\t{len(used)}\t{total:.0f}\t\t{ram:.0f}\t{len(s_used)}\t"
               f"{total / len(s_used):.2f}\t{lo:.2f}-{hi:.2f} (bootstrap over clusters)"]
    overlap = [(f"{a.method}, clusters used", s_used), (f"{a.method}, all clusters", s_all)]
    if a.subset:
        overlap.append((f"{a.method}, clusters used, not found in its other clusters", s_used - s_other))

    if a.ref_cdb:
        if combined is None:
            sys.exit("--ref_cdb needs --combined")
        ref = species_by_group(a.ref_cdb, a.ref_wdb, a.ref_prefix, combined)
        groups = a.ref_groups or sorted(ref)
        r_sp = set().union(*(ref.get(g, set()) for g in groups))
        line = f"{a.ref_name}\t{len(groups)}"
        if a.ref_cost:
            cost = {r[a.ref_cost_id]: r for r in csv.DictReader(open(a.ref_cost), delimiter="\t")}
            ids = [strip(g, a.ref_cost_strip) for g in groups]
            c, clo, chi = (sum(float(cost[i][k]) for i in ids) for k in ("cpu_h", "cpu_h_low", "cpu_h_high"))
            rram = max(float(cost[i]["peak_ram_gb"]) for i in ids)
            line += (f"\t{c:.0f}\t{clo:.0f}-{chi:.0f}\t{rram:.0f}\t{len(r_sp)}\t{c / len(r_sp):.2f}\t"
                     f"{clo / len(r_sp):.2f}-{chi / len(r_sp):.2f} (cost range)")
        else:
            line += f"\t\t\t\t{len(r_sp)}\t\t"
        summary.append(line)
        overlap += [(a.ref_name, r_sp), ("both", s_used & r_sp),
                    (f"only {a.method}", s_used - r_sp), (f"only {a.ref_name}", r_sp - s_used)]
        overlap += [(f"{a.ref_name} {g}", ref.get(g, set())) for g in groups]

    with open(os.path.join(a.outdir, "summary.tsv"), "w") as o:
        o.write("approach\tassemblies\tcpu_h\tcpu_h_range\tlargest_job_peak_ram_gb\tspecies\t"
                "cpu_h_per_species\tcpu_h_per_species_range\n")
        o.write("\n".join(summary) + "\n")
    with open(os.path.join(a.outdir, "overlap.tsv"), "w") as o:
        o.write("group\tspecies\n")
        for k, v in overlap:
            o.write(f"{k}\t{len(v)}\n")
    print(f"{a.method}: {len(used)} clusters used, mixed clusters {mixed} ({a.mixed}) -> {a.outdir}")


if __name__ == "__main__":
    main()

"""Shared input handling for the analysis scripts.

Partitions can be given in two ways (both may be combined):
  --partition NAME=PATH     a TSV with columns `sample` and `cluster`, as written by
                            Snakefile_clustering (results/{method}/{algorithm}_clusters.tsv)
  --cluster_summary PATH    a TSV with columns cluster_name, samples (comma-separated)
                            and found_in (labels like `sourmash:FKM`, several separated
                            by `;`); each label becomes one partition, named via
                            LABEL_NAMES (e.g. sourmash:FKM -> s_fkm)

Metadata is a TSV with one row per sample (sequencing run) and any label columns.
"""
import collections
import csv
import sys

LABEL_NAMES = {
    "mash:KM": "mash_km", "mash:FKM": "mash_fkm", "mash:DBSCAN": "mash_dbscan",
    "sourmash:KM": "s_km", "sourmash:FKM": "s_fkm", "sourmash:DBSCAN": "s_dbscan",
    "sourmash_abundance:KM": "sa_km", "sourmash_abundance:FKM": "sa_fkm",
    "sourmash_abundance:DBSCAN": "sa_dbscan",
}


def add_partition_args(ap):
    ap.add_argument("--partition", action="append", default=[], metavar="NAME=PATH",
                    help="partition TSV with columns sample, cluster (repeatable)")
    ap.add_argument("--cluster_summary", metavar="PATH",
                    help="cluster summary TSV (cluster_name, samples, found_in)")


def add_metadata_args(ap):
    ap.add_argument("--metadata", required=True, help="sample metadata TSV")
    ap.add_argument("--id_col", default="run_accession",
                    help="metadata column with the sample names used in the partitions")
    ap.add_argument("--derive", action="append", default=[], metavar="NEW=COL:SEP:INDEX",
                    help="add a column from part of another one, e.g. city=site:_:0 "
                         "(Copenhagen_RA -> Copenhagen); repeatable")
    ap.add_argument("--subset", action="append", default=[], metavar="COL=VALUE",
                    help="keep only samples with this value (repeatable)")
    ap.add_argument("--collapse_by", metavar="COL",
                    help="keep one sample per value of COL (e.g. a physical-sample id), "
                         "the one with the lowest --replica_col")
    ap.add_argument("--replica_col", default="REPLICA")


def read_partitions(a):
    part = collections.OrderedDict()
    for spec in a.partition:
        name, path = spec.split("=", 1)
        p = {}
        for r in csv.DictReader(open(path), delimiter="\t"):
            p[r["sample"]] = str(r["cluster"])
        part[name] = p
    if a.cluster_summary:
        found = collections.OrderedDict()
        for r in csv.DictReader(open(a.cluster_summary), delimiter="\t"):
            for lab in r["found_in"].strip('"').split(";"):
                name = LABEL_NAMES.get(lab, lab.replace(":", "_"))
                for s in r["samples"].split(","):
                    found.setdefault(name, {})[s] = r["cluster_name"]
        order = [n for n in LABEL_NAMES.values() if n in found] + \
                [n for n in found if n not in LABEL_NAMES.values()]
        for n in order:
            part[n] = found[n]
    if not part:
        sys.exit("No partitions given (use --partition and/or --cluster_summary).")
    return part


def read_metadata(a):
    meta = {}
    for r in csv.DictReader(open(a.metadata), delimiter="\t"):
        for spec in a.derive:
            new, rest = spec.split("=", 1)
            col, sep, idx = rest.split(":")
            r[new] = r[col].split(sep)[int(idx)]
        meta[r[a.id_col]] = r
    return meta


def select_samples(a, meta, part):
    """Samples present in the metadata and in every partition, after --subset and --collapse_by."""
    samples = sorted(s for s in meta if all(s in p for p in part.values()))
    missing = [n for n, p in part.items() if set(meta) - set(p)]
    if missing:
        print(f"note: metadata samples missing from partitions {missing} are left out",
              file=sys.stderr)
    for spec in a.subset:
        col, val = spec.split("=", 1)
        samples = [s for s in samples if meta[s][col] == val]
    if a.collapse_by:
        best = {}
        for s in samples:
            k = meta[s][a.collapse_by]
            if k not in best or int(meta[s][a.replica_col]) < int(meta[best[k]][a.replica_col]):
                best[k] = s
        samples = sorted(best.values())
    if not samples:
        sys.exit("No samples left after --subset/--collapse_by.")
    return samples

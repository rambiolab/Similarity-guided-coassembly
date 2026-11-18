import pandas as pd


# -----------------------------------------------------------------------------
# Sample sheet: tab-separated, columns  sample  read  path
#   one row per read file, read = r1 or r2 (paths to trimmed reads)
# -----------------------------------------------------------------------------

def read_sample_sheet(path):
    sheet = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    missing = {"sample", "read", "path"} - set(sheet.columns)
    if missing:
        raise ValueError(f"Sample sheet {path} lacks column(s): {', '.join(sorted(missing))}")
    sheet["read"] = sheet["read"].str.lower()
    bad = sheet.loc[~sheet["read"].isin(["r1", "r2"]), "read"].unique()
    if len(bad):
        raise ValueError(f"Sample sheet {path}: 'read' must be r1 or r2, found {list(bad)}")
    dup = sheet[sheet.duplicated(["sample", "read"], keep=False)]
    if len(dup):
        raise ValueError(f"Sample sheet {path}: duplicated sample/read rows:\n{dup}")
    wide = sheet.pivot(index="sample", columns="read", values="path")
    incomplete = wide[wide.isna().any(axis=1)].index.tolist()
    if incomplete:
        raise ValueError(f"Sample sheet {path}: samples without both r1 and r2: {incomplete}")
    # keep the order of first appearance in the sheet
    order = list(dict.fromkeys(sheet["sample"]))
    return {s: {"r1": wide.at[s, "r1"], "r2": wide.at[s, "r2"]} for s in order}


# -----------------------------------------------------------------------------
# Partition: tab-separated, columns  sample  cluster
#   one row per sample (output of the clustering workflow)
# -----------------------------------------------------------------------------

def read_partition(path, samples, prefix="cluster_"):
    """Return {cluster_name: [samples]} in order of first appearance.
    Cluster names are <prefix><cluster id>; samples labelled -1 (noise) are
    skipped. Every sample must be in the sample sheet."""
    part = pd.read_csv(path, sep="\t", dtype=str)
    missing = {"sample", "cluster"} - set(part.columns)
    if missing:
        raise ValueError(f"Partition {path} lacks column(s): {', '.join(sorted(missing))}")
    dup = part[part["sample"].duplicated(keep=False)]
    if len(dup):
        raise ValueError(f"Partition {path}: samples listed more than once:\n{dup}")
    unknown = sorted(set(part["sample"]) - set(samples))
    if unknown:
        raise ValueError(f"Partition {path}: samples not in the sample sheet: {unknown}")
    part = part[part["cluster"].str.strip() != "-1"]
    clusters = {}
    for sample, cid in zip(part["sample"], part["cluster"].str.strip()):
        clusters.setdefault(f"{prefix}{cid}", []).append(sample)
    return clusters

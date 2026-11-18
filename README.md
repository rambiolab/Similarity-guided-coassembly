# Similarity-guided co-assembly

Two Snakemake workflows:

1. **`Snakefile_clustering`** groups metagenomic samples into co-assembly groups (9 partitions).
2. **`Snakefile_assembly`** co-assembles, bins, dereplicates and classifies one chosen partition.

Both read the same sample sheet (see *Input*).

# Part 1: sample clustering

The workflow groups metagenomic samples for co-assembly. Reads are
sketched with three k-mer distances, and each distance matrix is clustered with
three algorithms (3 × 3 = 9 partitions).

| method | distance |
|---|---|
| `mash` | Mash distance (k = 31, sketch size 10,000) |
| `sourmash` | 1 − Jaccard similarity (k = 31, scaled = 1,000) |
| `sourmash_abundance` | 1 − angular similarity of abundance-tracking sketches (k = 31, scaled = 1,000) |

| algorithm | how the number of clusters is chosen |
|---|---|
| `KM` | k-medoids (PAM, BUILD initialisation, scikit-learn-extra). Elbow rule on the sum of distances to the nearest medoid. |
| `FKM` | Fuzzy k-medoids (FCMdd, m = 2, BUILD initialisation). Elbow rule on the fuzzy objective Σᵢ Σⱼ uᵢⱼ² d(i, medoid j); each sample is assigned to the cluster with its highest membership. |
| `DBSCAN` | DBSCAN with `min_samples = 1` (equivalent to a single-linkage cut at eps). eps is the value with the highest silhouette score (precomputed distances) on a grid of 1,000 values from 10⁻⁶ to 1. If no eps gives a valid partition, eps = 0.9 when the standard deviation of the distances is < 0.3, otherwise 0.1. |

**Elbow rule:** k is scanned from 2 to n − 1 (n = number of samples); the chosen
k is the point of the objective curve farthest from the straight line between the
first (k = 2) and last (k = n − 1) points.

All three algorithms are deterministic (no random initialisation). UMAP is used
for the cluster plots only (`random_state = 0`).

## Requirements

- Snakemake ≥ 8 and conda/mamba (`--use-conda`)
- GNU `time` at `/usr/bin/time` (used for per-step resource logs)
- `fuzzy_k_medoids_algorithm.py` from
  [rambiolab/Fuzzy-K-Medoids](https://github.com/rambiolab/Fuzzy-K-Medoids),
  included as a git submodule in `external/fuzzy_k_medoids/`. Clone with

  ```bash
  git clone --recurse-submodules https://github.com/rambiolab/Similarity-guided-coassembly
  ```

  or, in an existing clone, run `git submodule update --init`.

## Input

A tab-separated sample sheet with one row per read file (trimmed paired-end reads):

```
sample	read	path
sampleA	r1	/path/to/sampleA_1.fastq.gz
sampleA	r2	/path/to/sampleA_2.fastq.gz
```

Set its path as `samples:` in `config/config_clustering.yaml`.

## Run

```bash
snakemake -s Snakefile_clustering --configfile config/config_clustering.yaml \
    --use-conda --cores 16
```

Use `-n` for a dry run, or `--config methods=[sourmash] algorithms=[FKM]` to run
a subset.

## Output (`results/`)

| file | content |
|---|---|
| `{method}/sketches/` | per-sample sketches (`.msh` / `.sig`) |
| `{method}/distance_matrix.tsv` | square distance matrix; first column and header are sample names |
| `{method}/{algorithm}_clusters.tsv` | partition: columns `sample`, `cluster` (ids start at 1) |
| `{method}/{algorithm}_tuning_plot.png` | elbow curve (KM, FKM) or silhouette curve (DBSCAN) |
| `{method}/{algorithm}_clustering_plot.png` | UMAP of the samples coloured by cluster |

Each `{algorithm}_clusters.tsv` is a complete partition of all samples and is the
input for the co-assembly workflow.

Logs (`logs/`) include `/usr/bin/time -v` output per step; benchmarks are in
`benchmarks/`.

# Part 2: co-assembly of one partition

Choose one of the 9 partitions (e.g. `results/sourmash/FKM_clusters.tsv`) and set
it as `partition:` in `config/config_assembly.yaml`, with a separate `outdir:` per
partition. Each cluster `<id>` becomes `cluster_<id>`.

| step | tool | settings |
|---|---|---|
| co-assembly | MEGAHIT | default preset; all member samples' reads of a cluster |
| mapping | BWA-MEM, samtools | every member sample against its cluster's co-assembly |
| binning | MetaBAT2 | defaults; coverage from all member samples |
| bin quality | CheckM2 | `predict` |
| dereplication | dRep | `-comp 50 -con 10 -pa 0.9 -sa 0.95` on all bins of the partition |
| taxonomy | GTDB-Tk | `classify_wf` on the dereplicated genomes |

```bash
snakemake -s Snakefile_assembly --configfile config/config_assembly.yaml \
    --use-conda --cores 64
```

MEGAHIT needs a lot of memory for large clusters; set `mem_mb: megahit:` to what
your largest cluster needs and pass `--resources mem_mb=<memory of the node>` so
Snakemake does not start more assemblies than fit. An interrupted MEGAHIT run is
resumed with `--continue` when the workflow is restarted.

MetaBAT2 is run without a fixed seed (its default), so bins can differ slightly
between runs. Set `metabat2_extra: "--seed 1"` for reproducible bins.

Output (in `outdir`):

| file | content |
|---|---|
| `assembly/<cluster>/final.contigs.fa` | co-assembly |
| `bins/<cluster>/bin.N.fa` | MetaBAT2 bins |
| `checkm2/<cluster>/quality_report.tsv` | CheckM2 per bin |
| `drep_input/all_bins_quality.tsv` | all bins with cluster, completeness, contamination |
| `drep/` | dRep output; `dereplicated_genomes/` holds one genome per species |
| `gtdbtk/` | GTDB-Tk output |
| `final/genomes.tsv` | one row per dereplicated genome: cluster of origin, quality, taxonomy, novel species |
| `final/summary.tsv` | counts for the partition: clusters, bins, species (HQ, MQ, novel) |
| `logs/`, `benchmarks/` | per-step logs, `/usr/bin/time -v` output, Snakemake benchmarks |

HQ: ≥ 90% completeness and < 5% contamination; MQ: ≥ 50% and < 10%. A species
is novel when GTDB-Tk gives it no species name.

BAM files are deleted once the cluster is binned (`temp`).

## Data

The manuscript applies the workflow to, among others, the VEO sewage dataset
(ENA PRJEB68319; Becsei et al. 2024, *Nat Commun* 15:7551).

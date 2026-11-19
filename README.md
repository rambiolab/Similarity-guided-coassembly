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
| `DBSCAN` | DBSCAN with `min_samples = 1` (equivalent to a single-linkage cut at eps). eps is the value with the highest silhouette score (precomputed distances) on a grid of 1,000 values from 10⁻⁶ to 1.|

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


# Analysis scripts

Post-hoc analyses of the partitions and co-assemblies, as used for the manuscript.
They run outside Snakemake and need only numpy and scikit-learn; the clustering
environment has both:

```bash
conda env create -f envs/clustering.yaml -n coassembly-clustering
conda activate coassembly-clustering
```

| script | what it does |
|---|---|
| `cluster_agreement.py` | AMI of each partition vs metadata labels (site, city, season, …) with permutation nulls, and AMI between partitions |
| `season_split_search.py` | best split of the calendar year into 2–4 periods vs the partitions, with a null that repeats the search on shuffled dates |
| `assembly_cost_correction.py` | replaces the CPU hours of resumed MEGAHIT runs by a log-log fit on input size; totals per partition |
| `cpu_per_species.py` | CPU hours per recovered species for one partition (optionally a subset of samples, e.g. one city), with a bootstrap over clusters, and against a reference approach |
| `common.py` | shared input handling (partitions, metadata, subsets, replicas) |

## Inputs

**Partitions**, in either form (both can be combined):

- `--partition NAME=PATH`: the workflow's `results/{method}/{algorithm}_clusters.tsv`
  (columns `sample`, `cluster`). Repeat for several partitions, e.g.
  `--partition s_fkm=results/sourmash/FKM_clusters.tsv`.
- `--cluster_summary PATH`: one table of all clusters (columns `cluster_name`,
  `samples` comma-separated, `found_in` labels such as `sourmash:FKM`, several
  separated by `;`). Labels are renamed to short names (`sourmash:FKM` → `s_fkm`,
  `sourmash_abundance:KM` → `sa_km`, …; see `common.py`).

**Metadata**: a TSV with one row per sample (sequencing run) and the label columns.
Useful options:

- `--id_col` (default `run_accession`): column matching the partition sample names
- `--derive city=site:_:0`: build a column from part of another (`Copenhagen_RA` → `Copenhagen`)
- `--subset city=Copenhagen`: keep only those samples
- `--collapse_by ENA_ALIAS --replica_col REPLICA`: one run per physical sample (lowest replica)

## `cluster_agreement.py`

AMI (adjusted mutual information; Vinh, Epps & Bailey 2010, *JMLR* 11:2837–2854) is
chance-corrected: 0 for random agreement, 1 for identical partitions. It suits
unbalanced references and many small clusters better than ARI (Romano et al. 2016,
*JMLR* 17(134)).

Significance comes from permutation nulls that keep the partition fixed and shuffle
the labels (`--n_perm`, default 1,000, seed 1):

- `global`: labels shuffled across all samples. Does the partition follow the label at all?
- `within_<COL>` (from `--within LABEL=COL`): labels shuffled within each value of COL.
  For season within site: does the partition follow season beyond what site explains?

```bash
# all runs
python analysis/cluster_agreement.py --cluster_summary cluster_summary.tsv \
    --metadata VEO_sites.tsv --derive city=site:_:0 \
    --labels site city Season --within Season=site --outdir agreement/all_runs

# Copenhagen only
python analysis/cluster_agreement.py --cluster_summary cluster_summary.tsv \
    --metadata VEO_sites.tsv --derive city=site:_:0 --subset city=Copenhagen \
    --labels site Season --within Season=site --outdir agreement/Copenhagen_runs
```

Add `--collapse_by ENA_ALIAS` for one run per physical sample.

Output: `agreement.tsv` (partition, label, null, AMI, null mean, z, p),
`partition_vs_partition.tsv` (AMI between all pairs of partitions), and
`samples_used.txt`.

## `assembly_cost_correction.py`

Input: a per-cluster table with `cluster`, `input_gb`, `cpu_h`, `peak_ram_gb`,
`resumed` (True/False), and optionally `methods` (comma-separated partitions containing
the cluster). An assembly resumed with `megahit --continue` may have a time log that
covers only the resumed part. Those clusters get the CPU hours predicted by a fit on all
non-resumed clusters with at least `--min_gb` input:

    log(cpu_h) = a + b · log(input_gb)

`Snakefile_assembly` appends the time log on a restart (`/usr/bin/time -a`), so newer
runs keep both attempts; check that the table is built from all blocks of the log.

```bash
python analysis/assembly_cost_correction.py --resources assembly_resources.per_cluster.tsv \
    --outdir cost
```

Output: `resources_corrected.tsv` (adds `cpu_h_fit`, its ±1 residual-SD range,
`cpu_h_corrected`, `corrected`), `fit.tsv`, and `totals.tsv` (measured and corrected
CPU hours and largest-job peak RAM per partition).

## `cpu_per_species.py`

CPU hours per recovered species for one partition. Species are dRep secondary clusters
(95% ANI). A species counts for a set of assemblies if any of their bins falls in it,
using all bins from the partition's own `Cdb.csv` (not only the dereplicated winners,
which would undercount).

- Without `--combined`: species are the partition's own dRep clusters.
- With `--combined` (a dRep run over several approaches together): species are counted
  on the combined clusters, so the partition and a reference approach (`--ref_cdb`,
  e.g. per-site co-assembly) use one species definition. Own species are mapped to the
  combined table through their winner genome (`--wdb`; for the reference, the winner is
  taken to be the member present in the combined table if no `--ref_wdb` is given).
- `--subset COL=VALUE` restricts to clusters whose samples all match; `--mixed`
  (`fail`, `exclude`, `include`, `majority`) handles clusters that match only partly.
- Uncertainty: bootstrap over the partition's clusters (`--n_boot`, default 2,000);
  for the reference, the cost range from `--ref_cost` (columns `cpu_h`, `cpu_h_low`,
  `cpu_h_high`, `peak_ram_gb`).

Genome names are `<prefix><cluster>_bin.N.fa` (`scripts/prepare_drep_input.py` writes
`<cluster>_bin.N.fa`); give the prefix with `--genome_prefix` / `--ref_prefix`. For
partitions from `Snakefile_clustering` (cluster ids 1, 2, …), add `--cluster_prefix cluster_`.

```bash
# whole partition, own species
python analysis/cpu_per_species.py --partition s_fkm=results/sourmash/FKM_clusters.tsv \
    --cluster_prefix cluster_ --method s_fkm --metadata samples_meta.tsv \
    --resources cost/resources_corrected.tsv \
    --cdb out_s_fkm/drep/data_tables/Cdb.csv --wdb out_s_fkm/drep/data_tables/Wdb.csv \
    --outdir cpu_per_species/s_fkm

# Copenhagen, s_fkm vs per-site co-assembly on the combined species
python analysis/cpu_per_species.py --cluster_summary cluster_summary.tsv --method s_fkm \
    --metadata VEO_sites.tsv --derive city=site:_:0 --subset city=Copenhagen \
    --resources cost/resources_corrected.tsv \
    --cdb s_fkm/Cdb.csv --wdb s_fkm/Wdb.csv --genome_prefix s_fkm_ \
    --combined combined/Cdb.csv \
    --ref_name per-site --ref_cdb paper_site/Cdb.csv --ref_prefix paper_site_ \
    --ref_groups site_Copenhagen_RA site_Copenhagen_RD site_Copenhagen_RL \
    --ref_cost site_assembly_prediction.tsv --ref_cost_id site --ref_cost_strip site_ \
    --outdir cpu_per_species/Copenhagen_s_fkm
```

Output: `summary.tsv` (assemblies, CPU hours, largest-job peak RAM, species, CPU hours
per species with range), `overlap.tsv` (species counts and overlaps), and
`clusters.tsv` (per cluster: samples, input, CPU hours, peak RAM, species, species found
in no other cluster of the partition).


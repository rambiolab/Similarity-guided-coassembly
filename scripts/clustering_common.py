"""Shared helpers for the clustering scripts: loading, elbow rule, plotting, output."""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # non-GUI backend for saving images
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# --- Input ---------------------------------------------------------------------

def load_distance_matrix(path):
    """Read a square distance matrix TSV (first column and header = sample names)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.replace("NA", np.nan).astype(float)
    if list(df.index) != list(df.columns):
        raise ValueError(f"{path}: row and column names differ or are in a different order")
    d = df.values
    if np.isnan(d).any():
        rows, cols = np.where(np.isnan(d))
        pairs = [(df.index[r], df.columns[c]) for r, c in zip(rows[:5], cols[:5])]
        raise ValueError(f"{path}: {len(rows)} missing distances, e.g. {pairs}")
    return d, df.index.astype(str).tolist()


# --- Choosing k ------------------------------------------------------------------

def find_elbow(k_range, objective):
    """Elbow rule: the k whose point lies farthest from the straight line (chord)
    between the first and the last point of the objective curve."""
    ks = np.array(list(k_range), dtype=float)
    ys = np.array(objective, dtype=float)
    p1 = np.array([ks[0], ys[0]])
    p2 = np.array([ks[-1], ys[-1]])
    v = p2 - p1
    norm = np.linalg.norm(v)
    distances = []
    for k, y in zip(ks, ys):
        w = p1 - np.array([k, y])
        # magnitude of the 2-D cross product / chord length
        distances.append(abs(v[0] * w[1] - v[1] * w[0]) / norm)
    idx = int(np.argmax(distances))
    return int(ks[idx]), ys[idx], idx, p1, p2


def plot_elbow(k_range, objective, best_k, p1, p2, ylabel, title, path):
    plt.figure(figsize=(8, 5))
    plt.plot(list(k_range), objective, marker="o")
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color="gray", linestyle="--", label="chord")
    plt.axvline(best_k, color="red", linestyle="--", label=f"elbow k={best_k}")
    plt.xlabel("k (number of clusters)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


# --- Cluster plot -----------------------------------------------------------------

def make_palette():
    colors = []
    for name in ("tab20", "tab20b", "tab20c"):
        colors.extend(matplotlib.colormaps[name].colors)
    return colors


def umap_coords(distance_matrix):
    import umap  # imported here: only needed for the plot
    return umap.UMAP(metric="precomputed", random_state=0).fit_transform(distance_matrix)


def plot_clusters(distance_matrix, samples, labels, title, path):
    """UMAP of the distance matrix, coloured by cluster (plot only, no effect on clusters)."""
    coords = umap_coords(distance_matrix)
    palette = make_palette()
    unique_labels = sorted(set(labels))
    non_noise = [l for l in unique_labels if l != -1]
    cluster_colors = {lbl: palette[i % len(palette)] for i, lbl in enumerate(non_noise)}
    if -1 in unique_labels:
        cluster_colors[-1] = "gray"

    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    for (x, y), label, name in zip(coords, labels, samples):
        ax.scatter(x, y, s=1000, color=cluster_colors[label], edgecolors="black", alpha=0.8)
        ax.text(x, y, name, ha="center", va="center", fontsize=10)

    legend_elements = [
        Line2D([0], [0], marker="o", color="w",
               label=f"Cluster {label}" if label != -1 else "Noise",
               markerfacecolor=color, markeredgecolor="black", markersize=10)
        for label, color in cluster_colors.items()
    ]
    plt.legend(handles=legend_elements, loc="best")
    plt.title(title)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


# --- Output -------------------------------------------------------------------------

def write_clusters(samples, labels, path):
    """Write sample -> cluster. Cluster ids must already start at 1 (-1 = noise)."""
    pd.DataFrame({"sample": samples, "cluster": labels}).to_csv(path, sep="\t", index=False)

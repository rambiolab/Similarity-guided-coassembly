#!/usr/bin/env python3
"""DBSCAN (min_samples = 1) on a precomputed distance matrix; eps chosen by the
highest silhouette score over a grid of 1,000 values from 1e-6 to 1.

With min_samples = 1 every sample is a core point, so DBSCAN equals a
single-linkage cut at distance eps and no sample is labelled as noise.

Usage: cluster_DBSCAN.py distance_matrix.tsv clusters.tsv silhouette_plot.png cluster_plot.png
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score

import clustering_common as cc

EPS_GRID = np.linspace(1e-6, 1, 1000)

input_file, output_file, silhouette_plot, clust_plot = sys.argv[1:5]

distance_matrix, samples = cc.load_distance_matrix(input_file)


def fit(eps):
    return DBSCAN(eps=eps, min_samples=1, metric="precomputed").fit_predict(distance_matrix)


# --- Silhouette over the eps grid ---
scores = []
for eps in EPS_GRID:
    labels = fit(eps)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    _, counts = np.unique(labels[labels != -1], return_counts=True)
    # needs at least 2 clusters and not all singletons
    if n_clusters >= 2 and np.any(counts > 1):
        scores.append(silhouette_score(distance_matrix, labels, metric="precomputed"))
    else:
        scores.append(np.nan)
scores = np.array(scores)

if np.all(np.isnan(scores)):
    # no eps gives a valid partition: fall back to a fixed eps
    best_eps = 0.9 if np.nanstd(distance_matrix) < 0.3 else 0.1
    print(f"DBSCAN: no valid silhouette; fallback eps = {best_eps}")
else:
    best_eps = float(EPS_GRID[np.nanargmax(scores)])
    print(f"DBSCAN: best eps = {best_eps:.6f}, silhouette = {np.nanmax(scores):.4f}")

plt.figure(figsize=(8, 5))
plt.plot(EPS_GRID, scores, marker="o")
plt.axvline(best_eps, color="red", linestyle="--", label=f"eps={best_eps:.3f}")
plt.xlabel("eps")
plt.ylabel("Silhouette score")
plt.title("DBSCAN: eps chosen by maximum silhouette")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(silhouette_plot, dpi=300, bbox_inches="tight")
plt.close()

# --- Final clustering ---
labels = fit(best_eps)
labels = np.where(labels >= 0, labels + 1, -1)  # cluster ids start at 1, noise stays -1

cc.plot_clusters(distance_matrix, samples, labels,
                 title=f"DBSCAN clustering (eps={best_eps:.3f})", path=clust_plot)
cc.write_clusters(samples, labels, output_file)

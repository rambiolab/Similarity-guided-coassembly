#!/usr/bin/env python3
"""k-medoids (PAM, BUILD init) on a precomputed distance matrix; k by the elbow rule.

Usage: cluster_KM.py distance_matrix.tsv clusters.tsv elbow_plot.png cluster_plot.png
"""

import sys

from sklearn_extra.cluster import KMedoids

import clustering_common as cc

input_file, output_file, elbow_plot, clust_plot = sys.argv[1:5]

distance_matrix, samples = cc.load_distance_matrix(input_file)
n = distance_matrix.shape[0]


def fit(k):
    return KMedoids(n_clusters=k, metric="precomputed", method="pam",
                    init="build").fit(distance_matrix)


# --- Elbow: k = 2 .. n-1, objective = sum of distances to the nearest medoid ---
k_range = range(2, n)
objective = [fit(k).inertia_ for k in k_range]
best_k, _, _, p1, p2 = cc.find_elbow(k_range, objective)
print(f"KM: elbow at k = {best_k}")

cc.plot_elbow(k_range, objective, best_k, p1, p2,
              ylabel="Sum of distances to nearest medoid",
              title="k-medoids: elbow (largest distance to chord)",
              path=elbow_plot)

# --- Final clustering ---
labels = fit(best_k).labels_ + 1  # cluster ids start at 1

cc.plot_clusters(distance_matrix, samples, labels,
                 title=f"k-medoids clustering (k={best_k})", path=clust_plot)
cc.write_clusters(samples, labels, output_file)

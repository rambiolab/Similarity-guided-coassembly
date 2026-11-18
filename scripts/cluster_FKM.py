#!/usr/bin/env python3
"""Fuzzy k-medoids (FCMdd, m = 2, BUILD init) on a precomputed distance matrix;
k by the elbow rule on the fuzzy objective. Each sample is assigned to the
cluster with its highest membership.

Needs fuzzy_k_medoids_algorithm.py on the PYTHONPATH (separate repository).

Usage: cluster_FKM.py distance_matrix.tsv clusters.tsv elbow_plot.png cluster_plot.png
"""

import sys

import fuzzy_k_medoids_algorithm as FKM

import clustering_common as cc

M = 2.0
TOL = 1e-30

input_file, output_file, elbow_plot, clust_plot = sys.argv[1:5]

distance_matrix, samples = cc.load_distance_matrix(input_file)
n = distance_matrix.shape[0]

# --- Elbow: k = 2 .. n-1, objective J = sum_i sum_j u_ij^m * d(i, medoid_j) ---
k_range = range(2, n)
objective = []
for k in k_range:
    U, medoids = FKM.fuzzy_k_medoids(distance_matrix, n_clusters=k, m=M, tol=TOL)
    objective.append(FKM.fuzzy_objective(distance_matrix, U, medoids, m=M))
best_k, _, _, p1, p2 = cc.find_elbow(k_range, objective)
print(f"FKM: elbow at k = {best_k}")

cc.plot_elbow(k_range, objective, best_k, p1, p2,
              ylabel=r"Fuzzy objective $\sum_i \sum_j u_{ij}^2 \, d(i, m_j)$",
              title="Fuzzy k-medoids: elbow (largest distance to chord)",
              path=elbow_plot)

# --- Final clustering: hard assignment by highest membership ---
U, medoids = FKM.fuzzy_k_medoids(distance_matrix, n_clusters=best_k, m=M, tol=TOL)
labels = U.argmax(axis=1) + 1  # cluster ids start at 1

cc.plot_clusters(distance_matrix, samples, labels,
                 title=f"Fuzzy k-medoids clustering (k={best_k})", path=clust_plot)
cc.write_clusters(samples, labels, output_file)

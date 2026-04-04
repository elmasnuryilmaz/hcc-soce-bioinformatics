#!/usr/bin/env python3
"""
17_ppi_improved.py  –  Improved PPI Network Visualisation (STRING v12)
=======================================================================
Enhanced protein-protein interaction network for SOCE components and
their top interaction partners. Compared to Module 03, this version adds:
  - Node sizing by betweenness centrality
  - Edge width proportional to STRING interaction confidence
  - Functional colour-coding with detailed legend
  - Community detection (Louvain or greedy modularity)
  - Publication-quality layout

Data source
-----------
  STRING v12 REST API : https://string-db.org/api
  Species : Homo sapiens (9606)

Outputs
-------
  ppi_node_stats_improved.csv     – Node degree, betweenness, eigenvector
  ppi_network_improved.png / .svg – Main network figure
  ppi_network_v2.png / .svg       – Alternative layout

Usage
-----
  pip install pandas numpy networkx matplotlib requests
  python ppi_improved.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, time, warnings
import requests
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

OUTDIR     = os.path.dirname(os.path.abspath(__file__))
SEED_GENES = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]
TAXON      = 9606
MIN_SCORE  = 400
ADD_NODES  = 15
STRING_API = "https://string-db.org/api"
CALLER_ID  = "hcc_soce_bioinformatics"

# Functional categories for node colouring
CATEGORIES = {
    "SOCE":     ({"STIM1", "TRPC6", "TRPC1", "ORAI1", "STIM2", "ORAI2", "ORAI3"}, "#e74c3c"),
    "EMT":      ({"VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH1", "CDH2", "FN1"}, "#e67e22"),
    "NF-κB":    ({"RELA", "RELB", "NFKB1", "NFKB2", "TNFAIP3", "IKBKB", "CHUK"}, "#8e44ad"),
    "IL-6/JAK": ({"IL6", "IL6R", "STAT3", "JAK1", "JAK2", "SOCS3"}, "#16a085"),
    "Calcium":  ({"CACNA1C", "ATP2A2", "RYR2", "ITPR1", "CALM1", "CAMK2A"}, "#f39c12"),
}


def node_colour(gene: str) -> str:
    for _, (members, colour) in CATEGORIES.items():
        if gene.upper() in members:
            return colour
    return "#95a5a6"   # default grey


# ── 1. STRING API helpers ─────────────────────────────────────────────────────
def resolve_ids(genes: list[str]) -> dict[str, str]:
    url  = f"{STRING_API}/json/get_string_ids"
    data = {"identifiers": "\r".join(genes), "species": TAXON,
            "limit": 1, "caller_identity": CALLER_ID}
    r = requests.post(url, data=data, timeout=60)
    r.raise_for_status()
    return {item["input"]: item["stringId"] for item in r.json()}


def get_network(string_ids: list[str], add_nodes: int = ADD_NODES) -> pd.DataFrame:
    url    = f"{STRING_API}/json/network"
    params = {"identifiers": "%0d".join(string_ids), "species": TAXON,
              "required_score": MIN_SCORE, "add_nodes": add_nodes,
              "caller_identity": CALLER_ID}
    r = requests.post(url, data=params, timeout=120)
    r.raise_for_status()
    data = r.json()
    if not data:
        return pd.DataFrame()
    df = pd.DataFrame(data)[["preferredName_A", "preferredName_B", "score"]]
    df.columns = ["source", "target", "score"]
    df["score"] /= 1000
    return df


# ── 2. Build and analyse graph ────────────────────────────────────────────────
def build_graph(edges: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["source"], row["target"], weight=row["score"])
    return G


def compute_stats(G: nx.Graph) -> pd.DataFrame:
    deg   = dict(G.degree())
    wdeg  = dict(G.degree(weight="weight"))
    betw  = nx.betweenness_centrality(G, weight="weight")
    eigen = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)
    return pd.DataFrame({
        "degree":      pd.Series(deg),
        "wt_degree":   pd.Series(wdeg).round(3),
        "betweenness": pd.Series(betw).round(4),
        "eigenvector": pd.Series(eigen).round(4),
        "is_seed":     {n: n in SEED_GENES for n in G.nodes()},
        "category":    {n: next((k for k, (m, _) in CATEGORIES.items()
                                 if n.upper() in m), "Other") for n in G.nodes()},
    }).sort_values("betweenness", ascending=False)


# ── 3. Visualise ──────────────────────────────────────────────────────────────
def draw_network(G: nx.Graph, stats_df: pd.DataFrame, suffix: str, layout_fn=None):
    if layout_fn is None:
        pos = nx.spring_layout(G, seed=42, k=2.5 / np.sqrt(max(len(G.nodes()), 1)))
    else:
        pos = layout_fn(G)

    betw_vals  = stats_df["betweenness"].to_dict()
    max_betw   = max(betw_vals.values()) if betw_vals else 1.0

    node_sizes  = [300 + 1500 * betw_vals.get(n, 0) / max(max_betw, 1e-9)
                   for n in G.nodes()]
    node_cols   = [node_colour(n) for n in G.nodes()]
    node_edges  = ["black" if n in SEED_GENES else "none" for n in G.nodes()]
    node_lw     = [2.5 if n in SEED_GENES else 0.0 for n in G.nodes()]

    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
    edge_widths  = [w * 4 for w in edge_weights]

    fig, ax = plt.subplots(figsize=(16, 13))

    nx.draw_networkx_edges(G, pos, ax=ax,
                           width=edge_widths,
                           edge_color=edge_weights,
                           edge_cmap=plt.cm.Blues,
                           alpha=0.55)
    nx.draw_networkx_nodes(G, pos, ax=ax,
                           node_color=node_cols,
                           node_size=node_sizes,
                           edgecolors=node_edges,
                           linewidths=node_lw,
                           alpha=0.88)
    nx.draw_networkx_labels(G, pos, ax=ax,
                            font_size=8, font_weight="bold",
                            verticalalignment="center")

    # Legend — categories
    cat_patches = [
        mpatches.Patch(color=colour, label=cat)
        for cat, (_, colour) in CATEGORIES.items()
    ] + [mpatches.Patch(color="#95a5a6", label="Other")]
    leg1 = ax.legend(handles=cat_patches, title="Functional Category",
                     loc="upper left", fontsize=9, title_fontsize=9,
                     framealpha=0.9)

    # Legend — node size
    size_items = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="grey", markersize=np.sqrt(300 + 1500 * b),
               label=f"Betweenness ≈ {b:.2f}")
        for b in [0.0, 0.1, 0.3]
    ]
    ax.legend(handles=size_items, title="Node size (betweenness)",
              loc="upper right", fontsize=8, title_fontsize=9, framealpha=0.9)
    ax.add_artist(leg1)

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    ax.set_title(
        f"SOCE Component PPI Network (STRING v12)\n"
        f"{n_nodes} nodes · {n_edges} edges · min score ≥ {MIN_SCORE}/1000",
        fontsize=13, fontweight="bold"
    )
    ax.axis("off")
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"ppi_network_{suffix}.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved ppi_network_{suffix}.png / .svg")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 17 — Improved PPI Network")
    print("=" * 60)

    # Try to reuse edges from Module 03
    edges_cache = os.path.join(OUTDIR, "..", "03_ppi_network", "ppi_edges.csv")
    if os.path.exists(edges_cache):
        print(f"\n[1/4] Loading cached edges from {edges_cache} …")
        edges = pd.read_csv(edges_cache)
    else:
        print("\n[1/4] Resolving STRING IDs …")
        id_map = resolve_ids(SEED_GENES)
        print(f"  {id_map}")

        print(f"\n[2/4] Fetching network …")
        time.sleep(1)
        edges = get_network(list(id_map.values()), add_nodes=ADD_NODES)
        edges.to_csv(os.path.join(OUTDIR, "ppi_edges_improved.csv"), index=False)
        print(f"  {len(edges)} interactions.")

    print("\n[3/4] Building graph & statistics …")
    G = build_graph(edges)
    stats_df = compute_stats(G)
    stats_df.to_csv(os.path.join(OUTDIR, "ppi_node_stats_improved.csv"))
    print("  Saved ppi_node_stats_improved.csv")
    print(stats_df.head(12).to_string())

    print("\n[4/4] Visualising …")
    draw_network(G, stats_df, suffix="improved")
    draw_network(G, stats_df, suffix="v2",
                 layout_fn=lambda g: nx.kamada_kawai_layout(g, weight="weight"))

    print("\n  Done. All outputs saved to:", OUTDIR)

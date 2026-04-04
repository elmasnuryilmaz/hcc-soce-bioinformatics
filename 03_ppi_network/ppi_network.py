#!/usr/bin/env python3
"""
03_ppi_network.py  –  Protein-Protein Interaction Network (STRING v12)
======================================================================
Builds and visualises the PPI network for SOCE components and their
top interaction partners using the STRING database REST API.

Data source
-----------
  STRING v12 : https://string-db.org/api
  Species    : Homo sapiens (taxon 9606)
  Genes      : STIM1, TRPC6, TRPC1, ORAI1 + top partners

Outputs
-------
  ppi_edges.csv           – All STRING interactions (score ≥ threshold)
  ppi_node_stats.csv      – Degree, betweenness centrality per node
  ppi_network.png / .svg

Usage
-----
  pip install pandas networkx matplotlib requests
  python ppi_network.py

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

warnings.filterwarnings("ignore")

# ── Config ───────────────────────────────────────────────────────────────────
OUTDIR        = os.path.dirname(os.path.abspath(__file__))
SEED_GENES    = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]
TAXON         = 9606          # Homo sapiens
MIN_SCORE     = 400           # STRING combined score threshold (out of 1000)
MAX_PARTNERS  = 10            # Top interactors to include per seed gene
STRING_API    = "https://string-db.org/api"
CALLER_ID     = "hcc_soce_bioinformatics_analysis"

# Node colours
COLOURS = {
    "seed":    "#e74c3c",   # SOCE seed genes – red
    "emtf":    "#e67e22",   # EMT transcription factors – orange
    "nfkb":    "#8e44ad",   # NF-κB pathway – purple
    "other":   "#3498db",   # Other interactors – blue
}
EMT_TFS  = {"VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "CDH1", "CDH2", "FN1"}
NFKB_TFS = {"RELA", "RELB", "NFKB1", "NFKB2", "IL6R", "TNF", "TNFRSF1A", "ABCC2", "MAP3K7"}


# ── 1. Resolve STRING IDs ────────────────────────────────────────────────────
def resolve_ids(genes: list[str]) -> dict:
    """Returns {gene_symbol: string_id}."""
    url = f"{STRING_API}/json/get_string_ids"
    params = {
        "identifiers": "\r".join(genes),
        "species":     TAXON,
        "limit":       1,
        "caller_identity": CALLER_ID,
    }
    r = requests.post(url, data=params, timeout=60)
    r.raise_for_status()
    return {item["input"]: item["stringId"] for item in r.json()}


# ── 2. Get network ────────────────────────────────────────────────────────────
def get_network(string_ids: list[str], add_nodes: int = MAX_PARTNERS) -> pd.DataFrame:
    """Fetches interaction data from STRING for seed nodes + neighbours."""
    url = f"{STRING_API}/json/network"
    params = {
        "identifiers":        "%0d".join(string_ids),
        "species":            TAXON,
        "required_score":     MIN_SCORE,
        "add_nodes":          add_nodes,
        "caller_identity":    CALLER_ID,
    }
    r = requests.post(url, data=params, timeout=120)
    r.raise_for_status()
    interactions = r.json()
    if not interactions:
        return pd.DataFrame()
    df = pd.DataFrame(interactions)[
        ["preferredName_A", "preferredName_B", "score"]
    ]
    df.columns = ["source", "target", "score"]
    df["score"] = df["score"] / 1000  # normalise to 0-1
    return df


# ── 3. Build NetworkX graph ───────────────────────────────────────────────────
def build_graph(edges: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["source"], row["target"], weight=row["score"])
    return G


# ── 4. Compute node statistics ────────────────────────────────────────────────
def node_stats(G: nx.Graph) -> pd.DataFrame:
    deg   = dict(G.degree())
    betw  = nx.betweenness_centrality(G, weight="weight")
    eigen = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)
    df = pd.DataFrame({
        "degree":               pd.Series(deg),
        "betweenness":          pd.Series(betw).round(4),
        "eigenvector":          pd.Series(eigen).round(4),
        "is_seed":              pd.Series({n: n in SEED_GENES for n in G.nodes()}),
    }).sort_values("degree", ascending=False)
    return df


# ── 5. Visualise network ──────────────────────────────────────────────────────
def visualise(G: nx.Graph, edges: pd.DataFrame):
    """Spring-layout network plot with colour-coded node categories."""
    pos = nx.spring_layout(G, seed=42, k=2.0 / np.sqrt(len(G.nodes())))

    node_colours = []
    for n in G.nodes():
        if n in SEED_GENES:
            node_colours.append(COLOURS["seed"])
        elif n in EMT_TFS:
            node_colours.append(COLOURS["emtf"])
        elif n in NFKB_TFS:
            node_colours.append(COLOURS["nfkb"])
        else:
            node_colours.append(COLOURS["other"])

    node_sizes = [300 + 150 * G.degree(n) for n in G.nodes()]
    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]

    fig, ax = plt.subplots(figsize=(14, 11))
    nx.draw_networkx_edges(
        G, pos, ax=ax,
        width=[w * 3 for w in edge_weights],
        edge_color=edge_weights,
        edge_cmap=plt.cm.Blues,
        alpha=0.6,
    )
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_color=node_colours,
        node_size=node_sizes,
        alpha=0.9,
    )
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_weight="bold")

    legend_patches = [
        mpatches.Patch(color=COLOURS["seed"],  label="SOCE components"),
        mpatches.Patch(color=COLOURS["emtf"],  label="EMT transcription factors"),
        mpatches.Patch(color=COLOURS["nfkb"],  label="NF-κB / IL-6 pathway"),
        mpatches.Patch(color=COLOURS["other"], label="Other interactors"),
    ]
    ax.legend(handles=legend_patches, loc="upper left", framealpha=0.9)
    ax.set_title(
        f"SOCE Component PPI Network (STRING v12)\n"
        f"Min. interaction score ≥ {MIN_SCORE}/1000  |  "
        f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges",
        fontsize=13, fontweight="bold",
    )
    ax.axis("off")
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"ppi_network.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved ppi_network.png / .svg")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 03 — PPI Network (STRING v12)")
    print("=" * 60)

    print("\n[1/4] Resolving STRING IDs …")
    id_map = resolve_ids(SEED_GENES)
    string_ids = list(id_map.values())
    print(f"  Resolved: {id_map}")

    print(f"\n[2/4] Fetching network (min score {MIN_SCORE}) …")
    time.sleep(1)   # be polite to the STRING API
    edges = get_network(string_ids, add_nodes=MAX_PARTNERS)
    edges.to_csv(os.path.join(OUTDIR, "ppi_edges.csv"), index=False)
    print(f"  {len(edges)} interactions fetched. Saved ppi_edges.csv")

    print("\n[3/4] Building graph & computing statistics …")
    G = build_graph(edges)
    stats_df = node_stats(G)
    stats_df.to_csv(os.path.join(OUTDIR, "ppi_node_stats.csv"))
    print(stats_df.head(15).to_string())

    print("\n[4/4] Visualising …")
    visualise(G, edges)

    print("\n  Done. All outputs saved to:", OUTDIR)

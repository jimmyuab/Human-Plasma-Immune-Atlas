#!/usr/bin/env python
"""
build_data.py — rebuild the app's bundled data from the parent analysis project.

Reads   : I:\\Plasma immune atalas\\06_genetic_causality  (+ 09_tables)
Writes  : ./data/*.parquet  (compact, app-ready)

Run this whenever the upstream pipeline produces new results, then run
publish.sh / publish.bat to push the refreshed app to GitHub + Hugging Face.
"""
from __future__ import annotations
import os, sys, json, datetime
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)                      # parent analysis project
GC   = os.path.join(PROJ, "06_genetic_causality")
TB   = os.path.join(PROJ, "09_tables")
OUT  = os.path.join(HERE, "data")
os.makedirs(OUT, exist_ok=True)

manifest = {"built": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), "files": {}}


def log(msg):
    print(f"[build_data] {msg}", flush=True)


def save(df: pd.DataFrame, name: str, note: str = ""):
    p = os.path.join(OUT, name)
    df.to_parquet(p, compression="zstd", index=False)
    mb = round(os.path.getsize(p) / 1e6, 2)
    manifest["files"][name] = {"rows": int(len(df)), "MB": mb, "note": note}
    log(f"{name:34s} {len(df):>9,} rows  {mb:>7.2f} MB")


def read_tsv(path, **kw):
    if not os.path.exists(path):
        log(f"MISSING (skipped): {os.path.basename(path)}")
        return None
    return pd.read_csv(path, sep="\t", low_memory=False, **kw)


# ------------------------------------------------------------------ 1. pan-phenome cis-eQTL MR
src = os.path.join(GC, "cis_MR_ALL_finngen_results.tsv")
if os.path.exists(src):
    keep = ["gene_symbol", "SNP", "chr", "eqtl_Z", "eaf", "MR_beta", "MR_se", "MR_p",
            "OR", "OR_l95", "OR_u95", "phenocode", "phenotype", "num_cases",
            "num_controls", "category", "FDR", "immune_class"]
    log("streaming cis_MR_ALL_finngen_results.tsv (615 MB) ...")
    d = pd.concat(
        list(pd.read_csv(src, sep="\t", chunksize=500_000, low_memory=False, usecols=keep)),
        ignore_index=True)

    inst = (d.groupby("gene_symbol", as_index=False)
              .agg(SNP=("SNP", "first"), chr=("chr", "first"),
                   eqtl_Z=("eqtl_Z", "first"), eaf=("eaf", "first"),
                   immune_class=("immune_class", "first")))
    endp = (d.groupby("phenocode", as_index=False)
              .agg(phenotype=("phenotype", "first"), category=("category", "first"),
                   num_cases=("num_cases", "first"), num_controls=("num_controls", "first")))

    mr = d[["gene_symbol", "phenocode", "MR_beta", "MR_se", "MR_p",
            "OR", "OR_l95", "OR_u95", "FDR"]].copy()
    for c in ("gene_symbol", "phenocode"):
        mr[c] = mr[c].astype("category")
    for c in ("MR_beta", "MR_se", "OR", "OR_l95", "OR_u95"):
        mr[c] = mr[c].astype("float32")

    save(mr,   "mr_phenome_all.parquet", "cis-eQTL MR, every gene x every FinnGen R12 endpoint")
    save(inst, "instruments.parquet",    "one cis-eQTL instrument per gene (eQTLGen)")
    save(endp, "endpoints.parquet",      "FinnGen R12 endpoint dictionary")
    del d, mr

# ------------------------------------------------------------------ 2. colocalization
for fn, out, note in [
    ("coloc_ALL_finngen_results.tsv", "coloc_all.parquet", "coloc.abf on every FDR-significant locus"),
    ("coloc_phenome_results.tsv",     "coloc_phenome.parquet", "coloc, 28-endpoint core panel"),
    ("pqtl_MR_results.tsv",           "pqtl_mr.parquet", "protein-level MR (INTERVAL plasma pQTL)"),
    ("pqtl_coloc_results.tsv",        "pqtl_coloc.parquet", "protein-level coloc (INTERVAL)"),
    ("opengwas_replication.tsv",      "replication.parquet", "independent non-FinnGen GWAS replication"),
    ("novelty_engine_ranked.tsv",     "novelty_ranked.parquet", "novelty-priority engine score"),
    ("FINAL_evidence_tiers_repl.tsv", "evidence_tiers.parquet", "4-layer evidence tiering"),
    ("uk_cis_MR_results.tsv",         "uk_replication.parquet", "UK Biobank cross-population MR"),
    ("immune_cell_count_MR_results.tsv", "cell_count_mr.parquet", "gene -> blood immune cell count MR"),
    ("extended_cell_crp_MR_results.tsv", "cell_crp_mr.parquet", "gene -> cell/CRP/cytokine trait MR"),
    ("phenome_pleiotropy_axes.tsv",   "pleiotropy.parquet", "cross-category pleiotropic genes"),
    ("novelty_drug_direction.tsv",    "drug_direction.parquet", "genetic direction vs approved-drug mechanism"),
    ("novelty_class_enrichment.tsv",  "class_enrichment.parquet", "immune-class enrichment of causal targets"),
]:
    df = read_tsv(os.path.join(GC, fn))
    if df is not None:
        save(df, out, note)

# ------------------------------------------------------------------ 3. curated tables
for fn, out, note in [
    ("T6_disease_intelligence_final_table.tsv", "intelligence.parquet", "disease-intelligence final table"),
    ("T5_novelty_priority_targets.tsv", "t5_novelty.parquet", "top novelty-priority targets"),
    ("T4_druggable_immune_targets.tsv", "t4_druggable.parquet", "druggable plasma immune proteins"),
    ("T1_plasma_immune_universe.tsv",   "t1_universe.parquet", "plasma immune protein universe (HPA/Olink)"),
]:
    df = read_tsv(os.path.join(TB, fn))
    if df is not None:
        save(df, out, note)

# ------------------------------------------------------------------ manifest
with open(os.path.join(OUT, "manifest.json"), "w") as fh:
    json.dump(manifest, fh, indent=2)

total = sum(v["MB"] for v in manifest["files"].values())
log(f"DONE — {len(manifest['files'])} files, {total:.1f} MB total -> {OUT}")

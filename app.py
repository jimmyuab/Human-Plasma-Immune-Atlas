#!/usr/bin/env python
"""
Human Plasma Immune Atlas — interactive causal browser.

Open, genetics-anchored atlas of plasma immune proteins -> human disease phenome.
Every panel is computed live from the bundled result tables: pick, click Run, get
figures + a downloadable CSV. No account, no login, no upload required.

Layers
  1  cis-eQTL Mendelian randomisation   eQTLGen -> FinnGen R12 (672 genes x 2,466 endpoints)
  2  Bayesian colocalisation            coloc.abf, PP.H4
  3  Protein-level pQTL MR + coloc      INTERVAL (Sun 2018) plasma pQTL
  4  Independent replication            non-FinnGen consortium GWAS + UK Biobank cross-population
  5  Immune-cell / inflammation layer   gene -> blood cell counts, CRP, cytokines
  6  Novelty & therapeutic-direction engine
"""
from __future__ import annotations

import os
import re
import tempfile
import functools

# Hugging Face's free tier only offers ZeroGPU for Gradio Spaces (cpu-basic is now
# PRO-only), and the ZeroGPU runtime aborts with "No @spaces.GPU function detected
# during startup" unless at least one decorated function exists. This atlas is pure
# pandas/matplotlib and never needs a GPU, so the function below is only that startup
# handshake — it is deliberately not wired to any control, so no GPU is ever requested
# and no GPU quota is ever consumed. The decorator is a no-op off Hugging Face.
try:
    import spaces

    @spaces.GPU(duration=1)
    def _zerogpu_startup_handshake():
        return "ok"
except ImportError:
    pass

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import gradio as gr

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
TMP = os.path.join(tempfile.gettempdir(), "hpia_downloads")
os.makedirs(TMP, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 110, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "figure.facecolor": "white",
})

RISK = "#c0392b"      # OR > 1
PROT = "#1f6fb4"      # OR < 1
GREY = "#b8c0c8"
ACCENT = "#7d3c98"


# --------------------------------------------------------------------------- data
def _load(name: str) -> pd.DataFrame:
    p = os.path.join(DATA, name)
    if not os.path.exists(p):
        return pd.DataFrame()
    return pd.read_parquet(p)


print("[app] loading tables ...", flush=True)
MR = _load("mr_phenome_all.parquet")
INST = _load("instruments.parquet")
ENDP = _load("endpoints.parquet")
COLOC = _load("coloc_all.parquet")
PQ_MR = _load("pqtl_mr.parquet")
PQ_CO = _load("pqtl_coloc.parquet")
REPL = _load("replication.parquet")
UKREP = _load("uk_replication.parquet")
TWOPOP = _load("uk_twopop.parquet")
CELL = _load("cell_crp_mr.parquet")
CELLC = _load("cell_count_mr.parquet")
NOVEL = _load("novelty_ranked.parquet")
TIERS = _load("evidence_tiers.parquet")
INTEL = _load("intelligence.parquet")
PLEIO = _load("pleiotropy.parquet")
DRUGD = _load("drug_direction.parquet")

if not MR.empty:
    MR["gene_symbol"] = MR["gene_symbol"].astype(str)
    MR["phenocode"] = MR["phenocode"].astype(str)
    MR = MR.merge(ENDP, on="phenocode", how="left").merge(
        INST[["gene_symbol", "chr", "SNP", "immune_class"]], on="gene_symbol", how="left")
    MR["nlp"] = -np.log10(MR["MR_p"].clip(lower=1e-300))

GENES = sorted(MR["gene_symbol"].unique().tolist()) if not MR.empty else []
ENDP_LABEL = {}
if not ENDP.empty:
    ENDP = ENDP.sort_values("phenotype")
    ENDP_LABEL = {f"{r.phenotype}  [{r.phenocode}]": r.phenocode for r in ENDP.itertuples()}
ENDP_CHOICES = list(ENDP_LABEL.keys())
CATEGORIES = sorted(ENDP["category"].dropna().unique().tolist()) if not ENDP.empty else []
ICLASSES = sorted(INST["immune_class"].dropna().unique().tolist()) if not INST.empty else []

# gene x disease -> coloc PP.H4
CO_KEY = {}
if not COLOC.empty:
    CO_KEY = {(g, d): p for g, d, p in zip(COLOC["gene"], COLOC["disease_code"], COLOC["PP_H4"])}

N_SIG = int((MR["FDR"] < 0.05).sum()) if not MR.empty else 0
N_COLOC = int((COLOC["PP_H4"] >= 0.8).sum()) if not COLOC.empty else 0

print(f"[app] ready: {len(MR):,} MR tests | {len(GENES)} genes | {len(ENDP)} endpoints", flush=True)


def _csv(df: pd.DataFrame, stem: str) -> str:
    path = os.path.join(TMP, f"{stem}.csv")
    df.to_csv(path, index=False)
    return path


def _round(df: pd.DataFrame, n: int = 4) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_float_dtype(out[c]):
            small = out[c].abs().dropna()
            if len(small) and small[small > 0].min() < 1e-3:
                out[c] = out[c].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
            else:
                out[c] = out[c].round(n)
    return out


def _empty_fig(msg: str):
    fig, ax = plt.subplots(figsize=(7, 2.2))
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=11, color="#555")
    ax.axis("off")
    return fig


# FinnGen chapter names run to 109 characters ("III Diseases of the blood and
# blood-forming organs and certain disorders involving the immune mechanism (D3_)"),
# and there are 42 of them — unusable as axis labels. Key off the chapter code, which
# is stable, and fall back to tidying the free-text category names.
_CHAPTER = {
    "AB1_": "Infectious & parasitic", "CD2_": "Neoplasms (hospital)",
    "ICD-O-3": "Neoplasms (cancer reg.)", "D3_": "Blood & immune mechanism",
    "E4_": "Endocrine & metabolic", "F5_": "Mental & behavioural",
    "G6_": "Nervous system", "H7_": "Eye & adnexa", "H8_": "Ear & mastoid",
    "I9_": "Circulatory", "J10_": "Respiratory", "K11_": "Digestive",
    "L12_": "Skin & subcutaneous", "M13_": "Musculoskeletal",
    "N14_": "Genitourinary", "O15_": "Pregnancy & childbirth",
    "P16_": "Perinatal", "Q17": "Congenital", "R18_": "Symptoms & signs",
    "ST19_": "Injury & poisoning", "Z21_": "Health-service contact",
    "U22_": "Special-purpose codes",
}


def _short_cat(name: str) -> str:
    s = str(name).strip()
    for code, short in _CHAPTER.items():
        if f"({code})" in s or s.endswith(code):
            return short
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s)              # drop a trailing code
    s = re.sub(r"^(?:[IVXL]+)\s+", "", s)               # drop the roman numeral
    s = re.sub(r"\s+endpoints?$", "", s, flags=re.I)
    s = re.sub(r"^Comorbidities of\s+", "Comorbid: ", s, flags=re.I)
    s = re.sub(r"\s+from .*$", "", s)                   # "... from Katri Räikkönen"
    # no truncation: the labels are printed vertically, so full names always fit
    return s or "Unclassified"


# =========================================================== TAB 1 — atlas overview
@functools.lru_cache(maxsize=1)
def _overview_numbers():
    sig = MR[MR["FDR"] < 0.05]
    return dict(tests=len(MR), genes=MR["gene_symbol"].nunique(), endpoints=MR["phenocode"].nunique(),
                sig=len(sig), sig_genes=sig["gene_symbol"].nunique(),
                sig_diseases=sig["phenocode"].nunique(), coloc=N_COLOC)


def run_overview():
    if MR.empty:
        return _empty_fig("data not found"), pd.DataFrame(), None
    n = _overview_numbers()
    sig = MR[MR["FDR"] < 0.05].copy()

    fig = plt.figure(figsize=(13, 8.6))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)

    # A — evidence funnel
    ax = fig.add_subplot(gs[0, 0])
    steps = ["MR tests\nperformed", "FDR < 5%\ncausal pairs", "Colocalised\nPP.H4 >= 0.8",
             "Protein-level\npQTL support", "Independently\nreplicated"]
    vals = [n["tests"], n["sig"], n["coloc"],
            int((PQ_MR["MR_p"] < 0.05).sum()) if not PQ_MR.empty else 0,
            int(REPL["rep_sig"].sum()) if not REPL.empty and "rep_sig" in REPL else 0]
    y = np.arange(len(steps))[::-1]
    ax.barh(y, np.log10(np.maximum(vals, 1)), color=[GREY, RISK, ACCENT, "#e08214", "#2e8b57"],
            edgecolor="black", linewidth=.6)
    for yy, v in zip(y, vals):
        ax.text(np.log10(max(v, 1)) + .08, yy, f"{v:,}", va="center", fontsize=9, fontweight="bold")
    ax.set_yticks(y); ax.set_yticklabels(steps, fontsize=8)
    ax.set_xlabel("log10 count"); ax.set_title("A  Evidence funnel")
    ax.set_xlim(0, np.log10(max(vals) if max(vals) else 10) * 1.3)

    # B — causal pairs by ICD chapter
    ax = fig.add_subplot(gs[0, 1])
    cc = sig["category"].fillna("Unclassified").value_counts().head(12)[::-1]
    lbl = [c[:44] + ("..." if len(c) > 44 else "") for c in cc.index]
    ax.barh(range(len(cc)), cc.values, color="#4f77aa", edgecolor="black", linewidth=.5)
    ax.set_yticks(range(len(cc))); ax.set_yticklabels(lbl, fontsize=7)
    ax.set_xlabel("FDR<5% causal gene-disease pairs")
    ax.set_title("B  Disease chapters reached by the plasma immunome")

    # C — immune class
    ax = fig.add_subplot(gs[1, 0])
    tot = INST["immune_class"].value_counts()
    hit = sig.drop_duplicates("gene_symbol")["immune_class"].value_counts()
    order = tot.index.tolist()
    frac = [(hit.get(k, 0) / tot[k]) * 100 for k in order]
    ax.barh(range(len(order)), frac, color=ACCENT, alpha=.85, edgecolor="black", linewidth=.5)
    for i, k in enumerate(order):
        ax.text(frac[i] + .4, i, f"{hit.get(k,0)}/{tot[k]}", va="center", fontsize=7)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=7.5)
    ax.set_xlabel("% of class with >=1 causal disease link")
    ax.set_title("C  Which immune protein classes are causal")

    # D — direction of effect
    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(sig["MR_beta"], sig["nlp"], s=14, alpha=.6,
               c=np.where(sig["OR"] > 1, RISK, PROT), edgecolors="none")
    ax.axvline(0, color="k", lw=.7, ls="--")
    ax.set_xlabel("MR effect (log OR per 1-SD higher plasma protein)")
    ax.set_ylabel("-log10 P")
    ax.set_title("D  Effect direction of all causal pairs")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RISK, label="risk-increasing"),
                       Line2D([], [], marker="o", ls="", color=PROT, label="protective")],
              fontsize=8, frameon=False)
    lo, hi = np.percentile(sig["MR_beta"], [1, 99])
    ax.set_xlim(lo * 1.4, hi * 1.4)

    fig.suptitle("Human Plasma Immune Atlas — genome-anchored causal overview",
                 fontsize=13, fontweight="bold", y=0.985)

    summary = pd.DataFrame({
        "layer": ["1 cis-eQTL MR (eQTLGen -> FinnGen R12)", "2 Colocalisation (coloc.abf)",
                  "3 Plasma pQTL MR + coloc (INTERVAL)", "4 Independent replication",
                  "5 Two-population replication (Finland -> England)",
                  "6 Immune-cell / CRP layer", "7 Novelty engine"],
        "content": [f"{n['genes']} immune genes x {n['endpoints']} disease endpoints = {n['tests']:,} tests",
                    f"{len(COLOC):,} significant loci fine-mapped",
                    f"{len(PQ_MR):,} protein-level tests across "
                    f"{PQ_MR['disease'].nunique() if not PQ_MR.empty else 0} diseases",
                    f"{len(REPL)} non-FinnGen + {len(UKREP):,} UK Biobank cross-population tests",
                    (f"{len(TWOPOP)} pairs, {TWOPOP['phenocode'].nunique()} diseases, "
                     f"UK Biobank-only cohorts (no FinnGen overlap)" if not TWOPOP.empty else "not bundled"),
                    f"{len(CELL):,} gene x blood-trait tests",
                    f"{len(NOVEL):,} ranked candidate targets"],
        "result": [f"{n['sig']:,} pairs FDR<5% ({n['sig_genes']} genes, {n['sig_diseases']} diseases)",
                   f"{N_COLOC} loci with PP.H4 >= 0.8 (shared causal variant)",
                   f"{int((PQ_MR['FDR']<0.05).sum()) if not PQ_MR.empty else 0} significant at FDR<5% (protein level)",
                   f"{int(REPL['rep_sig'].sum()) if not REPL.empty else 0} replicated at P<0.05, all direction-concordant",
                   (f"{int(TWOPOP['concordant'].astype(bool).sum())} same-direction "
                    f"({100*TWOPOP['concordant'].astype(bool).mean():.0f}%), "
                    f"{int(TWOPOP['two_population_validated'].astype(bool).sum())} validated at P<0.05"
                    if not TWOPOP.empty else "-"),
                   f"{int((CELL['trait_FDR']<0.05).sum()) if not CELL.empty else 0} FDR<5% immune-cell links",
                   f"{int((NOVEL['category_label'].astype(str).str.contains('NOVEL')).sum())} novel nominations, "
                   f"{int((NOVEL['category_label'].astype(str)=='NOVEL protein-confirmed').sum())} protein-confirmed"],
    })
    return fig, summary, _csv(sig, "atlas_all_causal_pairs")


# =========================================================== TAB 2 — gene -> phenome
def run_gene(gene, fdr_only, top_n):
    if not gene:
        return _empty_fig("choose a gene"), _empty_fig(""), pd.DataFrame(), None
    d = MR[MR["gene_symbol"] == gene].copy()
    if d.empty:
        return _empty_fig(f"{gene}: no results"), _empty_fig(""), pd.DataFrame(), None
    d["category"] = d["category"].fillna("Unclassified")
    d = d.sort_values(["category", "phenotype"]).reset_index(drop=True)
    inst = INST[INST["gene_symbol"] == gene].iloc[0]

    # --- phenome-wide plot
    # FinnGen's 42 disease chapters differ ~100-fold in endpoint count, so laying the
    # x-axis out by raw endpoint index crushes the small chapters into a few pixels and
    # their names collide into an unreadable smear. Instead give every chapter an equal
    # one-unit slot (points spread inside it), split the scan across two stacked rows,
    # and print the chapter names vertically so each one is fully legible and can never
    # overlap its neighbour.
    cats = d["category"].unique().tolist()
    slot = {c: j for j, c in enumerate(cats)}
    d["x"] = (d.groupby("category").cumcount() + .5) / d["category"].map(
        d["category"].value_counts()) + d["category"].map(slot)

    colours = plt.cm.tab20(np.linspace(0, 1, max(len(cats), 2)))
    half = int(np.ceil(len(cats) / 2))
    rows = [cats[:half], cats[half:]] if len(cats) > 8 else [cats]

    thr = -np.log10(0.05 / len(d))
    # the top hits are named vertically above their point, so leave the upper part of
    # the axis free for that text instead of letting the names pile on top of each other
    ymax = max(float(d["nlp"].max()), thr) * 1.85
    sig = (d[d["FDR"] < 0.05].sort_values("nlp", ascending=False)
           .drop_duplicates("category")            # one name per chapter -> never collide
           .nlargest(min(int(top_n), 14), "nlp"))

    width = max(11.0, 0.62 * max(len(r) for r in rows))
    fig1, axes = plt.subplots(len(rows), 1, figsize=(width, 3.6 * len(rows) + 1.6))
    axes = np.atleast_1d(axes)
    for ax, row in zip(axes, rows):
        if not row:
            ax.axis("off"); continue
        x0, x1 = slot[row[0]], slot[row[-1]] + 1
        ticks, labels = [], []
        for j, c in enumerate(row):
            i = slot[c]
            sub = d[d["category"] == c]
            if j % 2:                      # alternate shading = visible chapter blocks
                ax.axvspan(i, i + 1, color="#f2f2f5", zorder=0, lw=0)
            ax.scatter(sub["x"], sub["nlp"], s=11, color=colours[i % len(colours)],
                       alpha=.85, edgecolors="none", zorder=2)
            ticks.append(i + .5)
            labels.append(f"{_short_cat(c)}  ({len(sub)})")
        ax.axhline(thr, color="k", ls="--", lw=.8, zorder=1)
        ax.text(x1, thr, "phenome-wide  ", va="bottom", ha="right", fontsize=7)
        for r in sig.itertuples():
            if x0 <= r.x <= x1:
                ax.annotate(str(r.phenotype)[:34], (r.x, r.nlp), fontsize=7, fontweight="bold",
                            xytext=(0, 5), textcoords="offset points", ha="center",
                            va="bottom", rotation=90,
                            color=RISK if r.OR > 1 else PROT, zorder=4)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, rotation=90, ha="center", va="top", fontsize=7.5)
        ax.tick_params(axis="x", length=2, pad=2)
        ax.set_ylabel("-log10 P (cis-eQTL MR)", fontsize=8.5)
        ax.set_xlim(x0, x1)
        ax.set_ylim(0, ymax)
    axes[0].set_title(f"{gene} — phenome-wide causal scan across {len(d):,} FinnGen R12 endpoints\n"
                      f"instrument {inst.SNP} (chr{inst.chr}, EAF {inst.eaf:.2f}) · class: {inst.immune_class}\n"
                      f"each disease chapter gets an equal slot; (n) = endpoints tested in it",
                      fontsize=10.5)
    fig1.tight_layout(h_pad=2.0)

    # --- forest of top hits
    tab = d[d["FDR"] < 0.05] if fdr_only else d
    tab = tab.nsmallest(int(top_n), "MR_p")
    if tab.empty:
        fig2 = _empty_fig(f"{gene}: no FDR<5% disease link (try unticking 'significant only')")
    else:
        t = tab.iloc[::-1]
        fig2, ax = plt.subplots(figsize=(9.5, max(2.4, 0.42 * len(t) + 1.1)))
        y = np.arange(len(t))
        col = [RISK if o > 1 else PROT for o in t["OR"]]
        ax.errorbar(t["OR"], y, xerr=[t["OR"] - t["OR_l95"], t["OR_u95"] - t["OR"]],
                    fmt="none", ecolor=col, lw=1.4, capsize=2.5)
        ax.scatter(t["OR"], y, s=52, c=col, edgecolors="black", linewidth=.6, zorder=3)
        ax.axvline(1, color="k", lw=.8, ls="--")
        lab = [f"{p[:46]}" + (f"  [H4={CO_KEY[(gene,c)]:.2f}]" if (gene, c) in CO_KEY else "")
               for p, c in zip(t["phenotype"].astype(str), t["phenocode"])]
        ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("OR per 1-SD genetically higher plasma protein (95% CI)")
        ax.set_title(f"{gene} — top causal disease effects" + ("  [H4 = colocalisation PP.H4]" if any("H4" in x for x in lab) else ""))
        fig2.tight_layout()

    show = tab[["phenotype", "phenocode", "category", "OR", "OR_l95", "OR_u95", "MR_p", "FDR",
                "num_cases", "num_controls"]].copy()
    show.insert(0, "direction", np.where(tab["OR"] > 1, "risk", "protective"))
    show["coloc_PP_H4"] = [CO_KEY.get((gene, c), np.nan) for c in tab["phenocode"]]
    return fig1, fig2, _round(show), _csv(d, f"{gene}_phenome_scan")


# =========================================================== TAB 3 — disease -> immunome
def run_disease(label, fdr_only, top_n):
    if not label:
        return _empty_fig("choose a disease endpoint"), _empty_fig(""), pd.DataFrame(), None
    code = ENDP_LABEL.get(label, label)
    d = MR[MR["phenocode"] == code].copy()
    if d.empty:
        return _empty_fig("no results for this endpoint"), _empty_fig(""), pd.DataFrame(), None
    meta = ENDP[ENDP["phenocode"] == code].iloc[0]

    # --- volcano
    fig1, ax = plt.subplots(figsize=(9.5, 5.6))
    ns = d[d["FDR"] >= 0.05]
    ss = d[d["FDR"] < 0.05]
    ax.scatter(ns["MR_beta"], ns["nlp"], s=13, color=GREY, alpha=.55, edgecolors="none")
    ax.scatter(ss["MR_beta"], ss["nlp"], s=34, c=np.where(ss["OR"] > 1, RISK, PROT),
               edgecolors="black", linewidth=.4, zorder=3)
    for r in ss.nlargest(min(int(top_n), 15), "nlp").itertuples():
        ax.annotate(r.gene_symbol, (r.MR_beta, r.nlp), fontsize=8, fontweight="bold",
                    xytext=(4, 3), textcoords="offset points")
    ax.axvline(0, color="k", lw=.7, ls="--")
    if len(ss):
        ax.axhline(ss["nlp"].min(), color="grey", lw=.8, ls=":")
        ax.text(ax.get_xlim()[1], ss["nlp"].min(), " FDR 5%", fontsize=7, va="bottom", ha="right")
    lo, hi = np.percentile(d["MR_beta"], [.5, 99.5])
    ax.set_xlim(lo * 1.5, hi * 1.5)
    ax.set_xlabel("MR effect (log OR per 1-SD higher plasma protein)")
    ax.set_ylabel("-log10 P")
    ax.set_title(f"{meta.phenotype}\n{len(d):,} plasma immune proteins tested · "
                 f"{len(ss)} causal at FDR<5% · {int(meta.num_cases):,} cases / {int(meta.num_controls):,} controls",
                 fontsize=10.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RISK, label="risk-increasing"),
                       Line2D([], [], marker="o", ls="", color=PROT, label="protective"),
                       Line2D([], [], marker="o", ls="", color=GREY, label="not significant")],
              fontsize=8, frameon=False, loc="upper left")
    fig1.tight_layout()

    # --- forest
    tab = (d[d["FDR"] < 0.05] if fdr_only else d).nsmallest(int(top_n), "MR_p")
    if tab.empty:
        fig2 = _empty_fig("no protein reaches FDR<5% for this endpoint")
    else:
        t = tab.iloc[::-1]
        fig2, ax = plt.subplots(figsize=(8.6, max(2.4, 0.4 * len(t) + 1.1)))
        y = np.arange(len(t))
        col = [RISK if o > 1 else PROT for o in t["OR"]]
        ax.errorbar(t["OR"], y, xerr=[t["OR"] - t["OR_l95"], t["OR_u95"] - t["OR"]],
                    fmt="none", ecolor=col, lw=1.4, capsize=2.5)
        ax.scatter(t["OR"], y, s=52, c=col, edgecolors="black", linewidth=.6, zorder=3)
        ax.axvline(1, color="k", lw=.8, ls="--")
        lab = [f"{g}  ({c})" + (f"  [H4={CO_KEY[(g,code)]:.2f}]" if (g, code) in CO_KEY else "")
               for g, c in zip(t["gene_symbol"], t["immune_class"].astype(str))]
        ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("OR per 1-SD genetically higher plasma protein (95% CI)")
        ax.set_title(f"{meta.phenotype} — causal plasma immune proteins")
        fig2.tight_layout()

    show = tab[["gene_symbol", "immune_class", "OR", "OR_l95", "OR_u95", "MR_p", "FDR"]].copy()
    show.insert(2, "direction", np.where(tab["OR"] > 1, "risk", "protective"))
    show["coloc_PP_H4"] = [CO_KEY.get((g, code), np.nan) for g in tab["gene_symbol"]]
    return fig1, fig2, _round(show), _csv(d, f"{code}_immunome_scan")


# =========================================================== TAB 4 — target browser
def run_browser(cats, classes, max_fdr, min_h4, coloc_only, novel_only, direction, n_show):
    if MR.empty:
        return _empty_fig("data not found"), pd.DataFrame(), None
    d = MR[MR["FDR"] < max(float(max_fdr), 1e-12)].copy()
    if cats:
        d = d[d["category"].isin(cats)]
    if classes:
        d = d[d["immune_class"].isin(classes)]
    d["coloc_PP_H4"] = [CO_KEY.get((g, c), np.nan) for g, c in zip(d["gene_symbol"], d["phenocode"])]
    if coloc_only:
        d = d[d["coloc_PP_H4"] >= float(min_h4)]
    if direction == "risk-increasing (block the protein)":
        d = d[d["OR"] > 1]
    elif direction == "protective (agonise / replace)":
        d = d[d["OR"] < 1]
    known = set(NOVEL.loc[NOVEL["category_label"].astype(str).str.contains("known", case=False),
                          "gene_symbol"].unique()) if not NOVEL.empty else set()
    d["known_drug_axis"] = d["gene_symbol"].isin(known)
    if novel_only:
        d = d[~d["known_drug_axis"]]
    if d.empty:
        return _empty_fig("no target matches these filters — relax them"), pd.DataFrame(), None

    d["nlfdr"] = -np.log10(d["FDR"].clip(lower=1e-300))
    d = d.sort_values("nlfdr", ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
    ax = axes[0]
    h4 = d["coloc_PP_H4"].fillna(0)
    sc = ax.scatter(d["nlfdr"], h4, s=30, c=np.where(d["OR"] > 1, RISK, PROT),
                    alpha=.75, edgecolors="black", linewidth=.35)
    ax.axhline(0.8, color="k", ls="--", lw=.8)
    ax.text(ax.get_xlim()[1], .81, "PP.H4 0.8 ", ha="right", fontsize=7)
    for r in d.head(12).itertuples():
        ax.annotate(f"{r.gene_symbol}", (r.nlfdr, (r.coloc_PP_H4 if pd.notna(r.coloc_PP_H4) else 0)),
                    fontsize=7.5, fontweight="bold", xytext=(4, 2), textcoords="offset points")
    ax.set_xlabel("-log10 FDR (causal strength)"); ax.set_ylabel("colocalisation PP.H4")
    ax.set_title("A  Causal strength vs shared-variant evidence")

    ax = axes[1]
    top = d.head(min(int(n_show), 22)).iloc[::-1]
    y = np.arange(len(top))
    ax.barh(y, np.log2(top["OR"]), color=[RISK if o > 1 else PROT for o in top["OR"]],
            edgecolor="black", linewidth=.5)
    ax.axvline(0, color="k", lw=.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g} -> {str(p)[:30]}" for g, p in zip(top["gene_symbol"], top["phenotype"])],
                       fontsize=7.5)
    ax.set_xlabel("log2 OR per 1-SD higher plasma protein")
    ax.set_title("B  Prioritised gene -> disease effects")
    fig.tight_layout()

    show = d[["gene_symbol", "immune_class", "phenotype", "category", "OR", "OR_l95", "OR_u95",
              "FDR", "coloc_PP_H4", "known_drug_axis"]].head(int(n_show)).copy()
    show.insert(4, "therapeutic_direction",
                np.where(show["OR"] > 1, "block / antagonise", "agonise / replace"))
    return fig, _round(show), _csv(d, "prioritised_targets")


# =========================================================== TAB 5 — evidence card
def gene_diseases(gene):
    if not gene or MR.empty:
        return gr.update(choices=[], value=None)
    d = MR[(MR["gene_symbol"] == gene)].nsmallest(60, "MR_p")
    ch = [f"{p}  [{c}]" for p, c in zip(d["phenotype"].astype(str), d["phenocode"])]
    return gr.update(choices=ch, value=ch[0] if ch else None)


def run_card(gene, dis_label):
    if not gene or not dis_label:
        return _empty_fig("choose a gene, then a disease"), pd.DataFrame(), None
    code = dis_label.split("[")[-1].rstrip("]")
    row = MR[(MR["gene_symbol"] == gene) & (MR["phenocode"] == code)]
    if row.empty:
        return _empty_fig("pair not found"), pd.DataFrame(), None
    r = row.iloc[0]
    h4 = CO_KEY.get((gene, code), np.nan)
    co = COLOC[(COLOC["gene"] == gene) & (COLOC["disease_code"] == code)]
    pq = PQ_MR[PQ_MR["gene"] == gene]
    pqc = PQ_CO[PQ_CO["gene"] == gene]
    rp = REPL[REPL["gene"] == gene]
    uk = UKREP[UKREP["gene_symbol"] == gene] if not UKREP.empty else pd.DataFrame()
    # two-population arm: same gene AND same endpoint, so this is a like-for-like
    # Finland-vs-England test rather than a match on gene alone
    tp = (TWOPOP[(TWOPOP["gene_symbol"] == gene) & (TWOPOP["phenocode"] == code)]
          if not TWOPOP.empty else pd.DataFrame())
    cl = CELL[(CELL["gene_symbol"] == gene) & (CELL["trait_FDR"] < 0.05)] if not CELL.empty else pd.DataFrame()

    fig = plt.figure(figsize=(13, 8.8))
    gs = fig.add_gridspec(2, 3, hspace=.45, wspace=.32)

    # A transcript-level MR
    ax = fig.add_subplot(gs[0, 0])
    ax.errorbar([r.OR], [0], xerr=[[r.OR - r.OR_l95], [r.OR_u95 - r.OR]], fmt="o",
                color=RISK if r.OR > 1 else PROT, capsize=4, ms=9, mec="black")
    ax.axvline(1, ls="--", color="k", lw=.8); ax.set_xscale("log")
    ax.set_yticks([]); ax.set_xlabel("OR (95% CI)")
    ax.set_title("A  Layer 1 — cis-eQTL MR")
    ax.text(.5, .78, f"OR = {r.OR:.2f}\nP = {r.MR_p:.2e}\nFDR = {r.FDR:.2e}",
            transform=ax.transAxes, ha="center", fontsize=9)

    # B coloc
    ax = fig.add_subplot(gs[0, 1])
    if not co.empty:
        pp = co.iloc[0][["PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4"]].astype(float).values
        ax.bar(["H0", "H1", "H2", "H3", "H4"], pp,
               color=[GREY, GREY, GREY, "#e08214", ACCENT], edgecolor="black", linewidth=.5)
        ax.axhline(.8, ls="--", color="k", lw=.8)
        ax.set_ylim(0, 1.05); ax.set_ylabel("posterior probability")
        ax.set_title(f"B  Layer 2 — colocalisation (PP.H4 = {pp[4]:.2f})")
    else:
        ax.axis("off"); ax.set_title("B  Layer 2 — colocalisation")
        ax.text(.5, .5, "not run\n(coloc restricted to\nsignificant loci)", ha="center", va="center", fontsize=9)

    # C protein-level pQTL
    ax = fig.add_subplot(gs[0, 2])
    if not pq.empty:
        p = pq.iloc[0]
        ax.errorbar([p.OR], [0], xerr=[[max(p.OR - p.OR_l95, 0)], [max(p.OR_u95 - p.OR, 0)]],
                    fmt="s", color="#e08214", capsize=4, ms=9, mec="black")
        ax.axvline(1, ls="--", color="k", lw=.8); ax.set_xscale("log"); ax.set_yticks([])
        conc = "concordant" if (p.OR - 1) * (r.OR - 1) > 0 else "DISCORDANT"
        h4p = pqc["PP_H4"].max() if not pqc.empty else np.nan
        ax.set_xlabel("protein OR (95% CI)")
        ax.text(.5, .78, f"INTERVAL plasma pQTL\nP = {p.MR_p:.2e}\nPP.H4 = {h4p:.2f}\n{conc} with transcript",
                transform=ax.transAxes, ha="center", fontsize=8.5)
    else:
        ax.axis("off")
        ax.text(.5, .5, "no plasma aptamer\nfor this protein\nin INTERVAL", ha="center", va="center", fontsize=9)
    ax.set_title("C  Layer 3 — protein-level pQTL MR")

    # D replication
    ax = fig.add_subplot(gs[1, 0])
    bars, cols, labs = [], [], []
    bars.append(-np.log10(max(r.MR_p, 1e-300))); cols.append(RISK); labs.append("FinnGen\n(discovery)")
    if not rp.empty and pd.notna(rp.iloc[0].get("rep_p", np.nan)):
        bars.append(-np.log10(max(float(rp.iloc[0]["rep_p"]), 1e-300)))
        cols.append("#2e8b57"); labs.append(str(rp.iloc[0]["rep_gwas"])[:16] + "\n(independent)")
    if not uk.empty:
        u = uk.nsmallest(1, "uk_p").iloc[0]
        bars.append(-np.log10(max(float(u["uk_p"]), 1e-300)))
        cols.append("#4f77aa"); labs.append("UK Biobank\n(cross-population)")
    if not tp.empty:
        t0 = tp.nsmallest(1, "uk_p").iloc[0]
        bars.append(-np.log10(max(float(t0["uk_p"]), 1e-300)))
        cols.append("#c0392b" if bool(t0["two_population_validated"]) else "#9ec5e8")
        labs.append("UKB same endpoint\n(two-population)")
    ax.bar(range(len(bars)), bars, color=cols, edgecolor="black", linewidth=.5)
    ax.axhline(-np.log10(0.05), ls="--", color="k", lw=.8)
    ax.set_xticks(range(len(bars))); ax.set_xticklabels(labs, fontsize=7.5)
    ax.set_ylabel("-log10 P"); ax.set_title("D  Layer 4 — replication")

    # E immune-cell layer
    ax = fig.add_subplot(gs[1, 1])
    if not cl.empty:
        t = cl.nsmallest(8, "trait_p").iloc[::-1]
        ax.barh(range(len(t)), t["trait_z"], color=[RISK if z > 0 else PROT for z in t["trait_z"]],
                edgecolor="black", linewidth=.5)
        ax.axvline(0, color="k", lw=.8)
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels([str(x)[:26] for x in t["trait_label"]], fontsize=7)
        ax.set_xlabel("MR Z (protein -> blood trait)")
    else:
        ax.axis("off")
        ax.text(.5, .5, "no FDR<5% blood\nimmune-cell / CRP link", ha="center", va="center", fontsize=9)
    ax.set_title("E  Layer 5 — immune-cell / inflammation")

    # F verdict
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    why = []
    ok_mr = bool(r.FDR < 0.05)
    ok_co = bool(pd.notna(h4) and h4 >= 0.8)
    ok_pq = bool(not pq.empty and pq.iloc[0]["MR_p"] < 0.05 and (pq.iloc[0]["OR"] - 1) * (r.OR - 1) > 0)
    ok_rp = bool(not rp.empty and bool(rp.iloc[0].get("rep_sig", False)))
    if ok_mr:
        why.append("MR FDR<5%")
    if ok_co:
        why.append(f"coloc PP.H4={h4:.2f}")
    if ok_pq:
        why.append("protein-level pQTL concordant")
    if ok_rp:
        why.append("independently replicated")
    # hierarchical: each tier requires every tier below it
    tier = 1
    if ok_mr:
        tier = 2
        if ok_co:
            tier = 3
            if ok_pq:
                tier = 4
                if ok_rp:
                    tier = 5
    verdict = {1: "not causal at FDR 5%", 2: "MR-only nomination", 3: "colocalised causal target",
               4: "protein-level causal target", 5: "replicated protein-level target"}[tier]
    act = "block / antagonise the protein" if r.OR > 1 else "agonise, replace or protect the protein"
    ax.text(0, 1, f"{gene}  →  {r.phenotype}", fontsize=12, fontweight="bold", va="top")
    ax.text(0, .84, f"Evidence tier {tier}/5\n{verdict}", fontsize=10.5, va="top", color=ACCENT,
            fontweight="bold")
    ax.text(0, .60, "Supported by:\n  · " + "\n  · ".join(why or ["nothing above threshold"]),
            fontsize=8.8, va="top")
    ax.text(0, .28, f"Direction: {'risk-increasing' if r.OR>1 else 'protective'}\n"
                    f"Therapeutic implication: {act}\nInstrument: {r.SNP} (chr{r.chr})\n"
                    f"Immune class: {r.immune_class}", fontsize=8.8, va="top")
    ax.set_title("F  Integrated verdict")

    fig.suptitle(f"Evidence card — {gene} and {r.phenotype}", fontsize=13, fontweight="bold", y=.985)

    card = pd.DataFrame({
        "layer": ["cis-eQTL MR (FinnGen R12)", "Colocalisation", "Plasma pQTL MR (INTERVAL)",
                  "Independent GWAS replication", "UK Biobank cross-population",
                  "Two-population (Finland -> England, same endpoint)",
                  "Immune-cell / CRP layer", "Verdict"],
        "result": [
            f"OR {r.OR:.3f} ({r.OR_l95:.3f}-{r.OR_u95:.3f}), P={r.MR_p:.2e}, FDR={r.FDR:.2e}",
            f"PP.H4 = {h4:.3f}" if pd.notna(h4) else "not evaluated at this locus",
            (f"OR {pq.iloc[0].OR:.3f}, P={pq.iloc[0].MR_p:.2e}" if not pq.empty else "no plasma aptamer available"),
            (f"{rp.iloc[0]['rep_gwas']}: P={rp.iloc[0]['rep_p']:.2e}, "
             f"{'concordant' if rp.iloc[0].get('rep_concordant') else 'discordant'}" if not rp.empty else "no independent GWAS matched"),
            (f"best P={uk['uk_p'].min():.2e} across {len(uk)} UKB endpoints" if not uk.empty else "not tested"),
            (f"{tp.iloc[0]['uk_trait']} ({tp.iloc[0]['match_method']}, {int(tp.iloc[0]['uk_ncase']):,} cases): "
             f"UK OR {tp.iloc[0]['uk_OR']:.3f}, P={tp.iloc[0]['uk_p']:.2e}, "
             f"{'VALIDATED — same direction' if bool(tp.iloc[0]['two_population_validated']) else ('same direction, not significant' if bool(tp.iloc[0]['concordant']) else 'opposite direction')}"
             if not tp.empty else "no independent UK Biobank GWAS of this endpoint"),
            (f"{len(cl)} blood traits at FDR<5% (top: {cl.nsmallest(1,'trait_p').iloc[0]['trait_label']})" if not cl.empty else "none at FDR<5%"),
            f"Tier {tier}/5 — {verdict}; {act}"],
    })
    return fig, card, _csv(card, f"{gene}_{code}_evidence_card")


# =========================================================== TAB 6 — immune-cell layer
def run_celllayer(gene, fdr_only):
    if CELL.empty:
        return _empty_fig("immune-cell layer not bundled"), pd.DataFrame(), None
    d = CELL[CELL["gene_symbol"] == gene].copy()
    if d.empty:
        return _empty_fig(f"{gene}: no blood-trait test"), pd.DataFrame(), None
    show = d[d["trait_FDR"] < 0.05] if fdr_only else d
    if show.empty:
        return (_empty_fig(f"{gene}: nothing at FDR<5% — untick to see all"),
                _round(d.nsmallest(30, "trait_p")), _csv(d, f"{gene}_blood_traits"))
    t = show.nsmallest(25, "trait_p").iloc[::-1]
    fig, ax = plt.subplots(figsize=(9.5, max(2.6, .34 * len(t) + 1.2)))
    ax.barh(range(len(t)), t["trait_z"], color=[RISK if z > 0 else PROT for z in t["trait_z"]],
            edgecolor="black", linewidth=.5)
    ax.axvline(0, color="k", lw=.8)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([f"{str(l)[:44]}  ({g})" for l, g in zip(t["trait_label"], t["trait_group"])],
                       fontsize=7.5)
    ax.set_xlabel("MR Z score (1-SD higher plasma protein -> trait)")
    ax.set_title(f"{gene} — causal effect on circulating immune-cell / inflammation traits\n"
                 f"{len(show)} of {len(d)} traits at FDR<5%", fontsize=10.5)
    fig.tight_layout()
    cols = ["trait_label", "trait_group", "trait_beta", "trait_se", "trait_z", "trait_p", "trait_FDR"]
    return fig, _round(show.nsmallest(40, "trait_p")[cols]), _csv(d, f"{gene}_blood_traits")


# =========================================================== UI
CSS = """
.gradio-container {max-width: 1400px !important}
footer {display:none !important}
#hdr h1 {margin-bottom:0}
"""

INTRO = f"""
# Human Plasma Immune Atlas
### An open, genetics-anchored causal map from the plasma immunome to the human disease phenome

**{len(GENES)} plasma immune proteins × {len(ENDP):,} FinnGen R12 disease endpoints = {len(MR):,} Mendelian-randomisation tests**,
of which **{N_SIG:,} are causal at FDR < 5 %** and **{N_COLOC} colocalise** (PP.H4 ≥ 0.8) with the disease signal.

Pick something on any tab and press **Run** — figures and a downloadable CSV are computed live from the
underlying result tables. Free and open: **no account, no login, no data upload.**
"""

METHODS = f"""
## How the atlas is built

Every layer below runs across the **whole FinnGen R12 phenome** — none is restricted to a curated
disease list. The one exception is stated explicitly in Layer 5.

**Layer 1 — cis-eQTL Mendelian randomisation.** One conditionally independent *cis*-eQTL per immune gene
(eQTLGen, n ≈ 31,684) is used as the instrument. Effects are transported to disease using the Zhu/Wald ratio
(β_MR = β_outcome / β_exposure) against every FinnGen R12 endpoint. Multiple testing is controlled with
Benjamini–Hochberg FDR across all tests.

**Layer 2 — Bayesian colocalisation.** For every significant locus, `coloc.abf` (Wakefield approximate Bayes
factors, per-SNP variance) tests whether the eQTL and the disease signal share one causal variant.
PP.H4 ≥ 0.8 = shared causal variant; this separates real causal effects from LD contamination.

**Layer 3 — protein-level confirmation.** Plasma pQTL instruments from INTERVAL (Sun et al. 2018,
fully public via the EBI GWAS Catalog) repeat the MR and colocalisation at the level of the circulating
protein itself. Discordance between transcript and protein is reported, not hidden — it is informative
(e.g. soluble decoy receptors act opposite to membrane signalling). This layer now spans
{PQ_MR['disease'].nunique() if not PQ_MR.empty else 0} diseases
({len(PQ_MR)} protein-level MR tests, {len(PQ_CO)} protein colocalisation loci).

**Layer 4 — independent replication.** Every nomination is re-tested in non-FinnGen consortium GWAS
(OpenGWAS) and in UK Biobank for cross-population support.

**Layer 5 — two-population replication (Finland → England).** The same instrument and the same Wald
ratio are applied to an independent UK Biobank GWAS of the *same* endpoint:
{len(TWOPOP)} pairs over {TWOPOP['phenocode'].nunique() if not TWOPOP.empty else 0} diseases, of which
{int(TWOPOP['concordant'].astype(bool).sum()) if not TWOPOP.empty else 0} agree in causal direction and
{int(TWOPOP['two_population_validated'].astype(bool).sum()) if not TWOPOP.empty else 0} also reach
P < 0.05. Only **single-cohort UK Biobank** datasets are eligible: several large public meta-analyses of
these endpoints silently include FinnGen, so using them would mean replicating FinnGen in FinnGen. This
layer cannot span the whole phenome, because it requires an independent UK GWAS of the same endpoint to
exist; the matching route used for each pair is recorded in the downloadable table.

**Layer 6 — immune-cell and inflammation layer.** The same instruments are tested against circulating
blood-cell counts, CRP and cytokine traits, to show *how* the protein changes the immune system.

**Layer 7 — novelty and therapeutic direction.** Targets are scored on causal strength, colocalisation,
pleiotropy across disease chapters, druggability, cell source, protein-level confirmation and
two-population replication, with a penalty for already-approved drug axes and for the MHC region.
Effect direction is translated into a therapeutic action:
OR > 1 → block the protein; OR < 1 → agonise or replace it.

### Interpretation and limits
* MR estimates are **lifelong genetically-proxied** effects, not the effect of a short drug course.
* An HLA/MHC-region hit is held at *nomination* — long-range LD there defeats colocalisation.
* eQTL instruments proxy transcript, not always circulating protein; Layer 3 is the arbiter.
* About half of the causal proteins have no public SomaScan aptamer (largely intracellular or
  MHC-region) and therefore cannot reach a protein-level tier from login-free data at all.
* Nothing here is clinical advice or a validated diagnostic.

### Data sources
eQTLGen Consortium *cis*-eQTL · FinnGen release R12 · INTERVAL plasma pQTL (Sun 2018, EBI GWAS Catalog) ·
OpenGWAS/IEU · UK Biobank GWAS · Human Protein Atlas · UK Biobank Olink Explore protein universe (coding 143) ·
MSigDB C7/C8 immunologic signatures.

### Citation
> Zhao J. *Human Plasma Immune Atlas: an open, genetics-anchored causal map of the plasma immunome
> across the human disease phenome.* Hugging Face Space, 2026.

Code and data: https://github.com/jimmyuab/Human-Plasma-Immune-Atlas — released under the MIT licence.
The underlying GWAS/eQTL datasets remain under their own licences and conditions of use.
"""

with gr.Blocks(title="Human Plasma Immune Atlas", theme=gr.themes.Soft(), css=CSS) as demo:
    gr.Markdown(INTRO, elem_id="hdr")

    # ---------------- overview
    with gr.Tab("1 · Atlas overview"):
        gr.Markdown("Press **Run** for a live summary of the whole atlas: evidence funnel, disease "
                    "chapters reached, which immune protein classes are causal, and effect directions.")
        b0 = gr.Button("Run atlas overview", variant="primary")
        p0 = gr.Plot(label="Atlas overview")
        t0 = gr.Dataframe(label="Layer-by-layer summary", wrap=True)
        f0 = gr.File(label="Download every FDR<5% causal pair (CSV)")
        b0.click(run_overview, None, [p0, t0, f0])

    # ---------------- gene -> phenome
    with gr.Tab("2 · Gene → phenome"):
        gr.Markdown("Choose a plasma immune protein and scan its causal effect across **every** disease "
                    "endpoint in FinnGen R12.")
        with gr.Row():
            g2 = gr.Dropdown(GENES, value=("IL6ST" if "IL6ST" in GENES else (GENES[0] if GENES else None)),
                             label="Plasma immune protein (gene)", filterable=True, scale=3)
            s2 = gr.Checkbox(True, label="table: significant only (FDR<5%)", scale=1)
            n2 = gr.Slider(5, 40, 15, step=5, label="rows / labels", scale=1)
        b2 = gr.Button("Run phenome scan", variant="primary")
        p2a = gr.Plot(label="Phenome-wide scan")
        p2b = gr.Plot(label="Top causal effects")
        t2 = gr.Dataframe(label="Results", wrap=True)
        f2 = gr.File(label="Download full phenome scan for this gene (CSV)")
        b2.click(run_gene, [g2, s2, n2], [p2a, p2b, t2, f2])

    # ---------------- disease -> immunome
    with gr.Tab("3 · Disease → immunome"):
        gr.Markdown("Choose a disease and see which plasma immune proteins causally drive it. "
                    "Type to search the 2,400+ endpoints.")
        with gr.Row():
            d3 = gr.Dropdown(ENDP_CHOICES,
                             value=next((c for c in ENDP_CHOICES if "Rheumatoid arthritis" in c),
                                        ENDP_CHOICES[0] if ENDP_CHOICES else None),
                             label="Disease endpoint (FinnGen R12)", filterable=True, scale=3)
            s3 = gr.Checkbox(True, label="table: significant only (FDR<5%)", scale=1)
            n3 = gr.Slider(5, 40, 15, step=5, label="rows / labels", scale=1)
        b3 = gr.Button("Run immunome scan", variant="primary")
        p3a = gr.Plot(label="Volcano across the plasma immunome")
        p3b = gr.Plot(label="Causal proteins")
        t3 = gr.Dataframe(label="Results", wrap=True)
        f3 = gr.File(label="Download full immunome scan for this disease (CSV)")
        b3.click(run_disease, [d3, s3, n3], [p3a, p3b, t3, f3])

    # ---------------- browser
    with gr.Tab("4 · Target browser"):
        gr.Markdown("Filter the whole atlas down to a therapeutic short-list, then download it.")
        with gr.Row():
            c4 = gr.Dropdown(CATEGORIES, multiselect=True, label="Disease chapter (blank = all)", scale=2)
            k4 = gr.Dropdown(ICLASSES, multiselect=True, label="Immune protein class (blank = all)", scale=2)
        with gr.Row():
            f4a = gr.Slider(1e-6, 0.05, 0.05, label="max FDR", scale=1)
            h4a = gr.Slider(0.5, 0.99, 0.8, step=0.01, label="min colocalisation PP.H4", scale=1)
            co4 = gr.Checkbox(True, label="colocalised only", scale=1)
            nv4 = gr.Checkbox(False, label="exclude known drug axes", scale=1)
        with gr.Row():
            dir4 = gr.Radio(["any", "risk-increasing (block the protein)", "protective (agonise / replace)"],
                            value="any", label="Effect direction", scale=2)
            n4 = gr.Slider(10, 300, 60, step=10, label="rows to return", scale=1)
        b4 = gr.Button("Run target browser", variant="primary")
        p4 = gr.Plot(label="Prioritisation")
        t4 = gr.Dataframe(label="Prioritised targets", wrap=True)
        fl4 = gr.File(label="Download the filtered short-list (CSV)")
        b4.click(run_browser, [c4, k4, f4a, h4a, co4, nv4, dir4, n4], [p4, t4, fl4])

    # ---------------- evidence card
    with gr.Tab("5 · Evidence card"):
        gr.Markdown("One gene, one disease, all six layers of evidence on a single page — "
                    "with an integrated tier and the implied therapeutic direction.")
        with gr.Row():
            g5 = gr.Dropdown(GENES, value=("IL6ST" if "IL6ST" in GENES else (GENES[0] if GENES else None)),
                             label="Protein (gene)", filterable=True, scale=2)
            d5 = gr.Dropdown([], label="Disease (top 60 for that gene)", filterable=True, scale=3)
        b5 = gr.Button("Build evidence card", variant="primary")
        p5 = gr.Plot(label="Evidence card")
        t5 = gr.Dataframe(label="Layer-by-layer evidence", wrap=True)
        f5 = gr.File(label="Download this evidence card (CSV)")
        g5.change(gene_diseases, g5, d5)
        demo.load(gene_diseases, g5, d5)
        b5.click(run_card, [g5, d5], [p5, t5, f5])

    # ---------------- cell layer
    with gr.Tab("6 · Immune-cell layer"):
        gr.Markdown("How does the protein change the immune system itself? Causal effects on circulating "
                    "blood-cell counts, CRP and cytokine traits.")
        with gr.Row():
            g6 = gr.Dropdown(sorted(CELL["gene_symbol"].unique().tolist()) if not CELL.empty else [],
                             value=("IL6ST" if not CELL.empty and "IL6ST" in set(CELL["gene_symbol"]) else None),
                             label="Protein (gene)", filterable=True, scale=3)
            s6 = gr.Checkbox(True, label="significant only (FDR<5%)", scale=1)
        b6 = gr.Button("Run immune-cell scan", variant="primary")
        p6 = gr.Plot(label="Blood immune traits")
        t6 = gr.Dataframe(label="Results", wrap=True)
        f6 = gr.File(label="Download all blood-trait results for this gene (CSV)")
        b6.click(run_celllayer, [g6, s6], [p6, t6, f6])

    # ---------------- methods
    with gr.Tab("7 · Methods & citation"):
        gr.Markdown(METHODS)

if __name__ == "__main__":
    # ssr_mode=False: Gradio 5's experimental server-side rendering answers the
    # Space's readiness probe with 405s, so the Space never leaves APP_STARTING.
    demo.queue(max_size=32).launch(server_name="0.0.0.0",
                                   server_port=int(os.environ.get("PORT", 7860)),
                                   ssr_mode=False)

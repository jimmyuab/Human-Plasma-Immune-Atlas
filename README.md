---
title: Human Plasma Immune Atlas
emoji: 🩸
colorFrom: red
colorTo: indigo
sdk: gradio
sdk_version: 5.49.1
app_file: app.py
pinned: true
license: mit
short_description: Open causal map from the plasma immunome to the human disease phenome
tags:
  - biology
  - genomics
  - mendelian-randomization
  - immunology
  - drug-target
---

# Human Plasma Immune Atlas

**An open, genetics-anchored causal map from the plasma immunome to the human disease phenome.**

672 plasma immune proteins × 2,466 FinnGen R12 disease endpoints = **1,656,872 Mendelian-randomisation
tests**, of which **1,016 are causal at FDR < 5 %** and **417 colocalise** (PP.H4 ≥ 0.8) with the disease
signal.

👉 **Live app — click and run, no login, no account, no upload:**
https://huggingface.co/spaces/jianlizhao/Human-Plasma-Immune-Atlas

---

## What you can do

| Tab | What it does |
|---|---|
| **1 · Atlas overview** | One click → evidence funnel, disease chapters reached by the immunome, which immune protein classes are causal, effect directions. Download every causal pair. |
| **2 · Gene → phenome** | Pick a plasma immune protein → phenome-wide causal scan over *every* FinnGen endpoint + forest plot + CSV. |
| **3 · Disease → immunome** | Pick a disease → volcano over all 672 plasma immune proteins, ranked causal proteins + CSV. |
| **4 · Target browser** | Filter by disease chapter, immune class, FDR, colocalisation, direction, novelty → therapeutic short-list + CSV. |
| **5 · Evidence card** | One gene × one disease, all six evidence layers on a single page, with an integrated 1–5 tier and the implied therapeutic action. |
| **6 · Immune-cell layer** | How the protein changes the immune system itself — blood cell counts, CRP, cytokines. |
| **7 · Methods & citation** | Full methodology, interpretation limits, data sources. |

Every figure is computed **live** from the bundled result tables — nothing is a static image, and every
panel has a matching CSV download.

## Evidence layers

1. **cis-eQTL Mendelian randomisation** — eQTLGen (n ≈ 31,684) instrument → FinnGen R12, Zhu/Wald ratio, BH-FDR.
2. **Bayesian colocalisation** — `coloc.abf` (Wakefield ABF), PP.H4 ≥ 0.8 = shared causal variant.
3. **Protein-level pQTL MR + coloc** — INTERVAL plasma pQTL (Sun 2018, public via EBI GWAS Catalog).
4. **Independent replication** — non-FinnGen consortium GWAS (OpenGWAS) + UK Biobank cross-population.
5. **Immune-cell / inflammation layer** — the same instruments against blood cell counts, CRP, cytokines.
6. **Novelty & therapeutic-direction engine** — causal strength + coloc + pleiotropy + druggability,
   penalised for known drug axes and the MHC region; OR > 1 → block the protein, OR < 1 → agonise/replace.

Positive controls recovered by the pipeline include IL6ST → rheumatoid arthritis (tocilizumab),
CTLA4 → RA / autoimmune thyroid disease (abatacept) and TNFRSF1A → ankylosing spondylitis (etanercept).

## Data sources

eQTLGen Consortium *cis*-eQTL · FinnGen release R12 · INTERVAL plasma pQTL (Sun et al. 2018) ·
OpenGWAS / IEU · UK Biobank GWAS · Human Protein Atlas · UK Biobank Olink Explore protein universe ·
MSigDB C7 / C8.

## Run it locally

```bash
git clone https://huggingface.co/spaces/jianlizhao/Human-Plasma-Immune-Atlas
cd Human-Plasma-Immune-Atlas
pip install -r requirements.txt
python app.py          # -> http://localhost:7860
```

## Rebuild the bundled data

```bash
python build_data.py   # re-reads the upstream analysis project into ./data
```

## Limits

* MR estimates are **lifelong genetically-proxied** effects, not the effect of a drug course.
* MHC/HLA-region hits are held at *nomination* — long-range LD defeats colocalisation there.
* eQTL instruments proxy transcript, not always circulating protein; layer 3 is the arbiter.
* Research resource only — **not** clinical advice or a validated diagnostic.

## Citation

> Zhao J. *Human Plasma Immune Atlas: an open, genetics-anchored causal map of the plasma immunome
> across the human disease phenome.* 2026.

## License

MIT. Underlying GWAS/eQTL datasets remain under their original licences and use conditions.

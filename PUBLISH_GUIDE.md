# Publishing guide — Human Plasma Immune Atlas

Everything for the app lives in this one folder (`I:\Plasma immune atalas\Huggeface`).
Publishing is: **create two empty repos once → run one command forever after.**

```
Huggeface/
├── app.py                     the Gradio app (7 tabs, everything computed live)
├── build_data.py              rebuilds ./data from ../06_genetic_causality and ../09_tables
├── data/                      20 parquet tables, 54 MB — bundled with the Space
├── publish.py                 one command: rebuild data → commit → push GitHub + HF
├── publish.bat / publish.sh   double-click wrappers
├── README.md                  Space front page (contains the HF YAML header)
├── requirements.txt           runtime deps (gradio version comes from README sdk_version)
├── .gitattributes             git-lfs rules (the 53 MB parquet needs LFS)
├── .github/workflows/         auto-mirror GitHub → Hugging Face on every push
├── LICENSE                    MIT
└── CITATION.cff
```

---

## Step 0 — install the two tools (once)

| Tool | Why | Get it |
|---|---|---|
| **git-lfs** | `data/mr_phenome_all.parquet` is 53 MB; both GitHub and Hugging Face require LFS above 10 MB | https://git-lfs.com |
| **GitHub CLI** (optional) | lets `publish.py` create the repo for you | https://cli.github.com |

Check: `git lfs version`

---

## Step 1 — create the two empty repos (once)

**GitHub** → https://github.com/new
* Repository name: `Human-Plasma-Immune-Atlas`
* Public · **do not** add a README, .gitignore or licence (this folder already has them)

**Hugging Face Space** → https://huggingface.co/new-space
* Owner `jianlizhao`, Space name: `Human-Plasma-Immune-Atlas`
* SDK: **Gradio** · Hardware: **CPU basic (free)**
* Visibility: **Public** ← this is what makes it usable by anyone with no login

---

## Step 2 — put your tokens where the script can read them (once)

Tokens are read from `..\.secrets\` (already git-ignored) or from environment variables.
**Create the tokens yourself in the browser and paste them into these files — nothing else reads them.**

| File | Token | Where to make it |
|---|---|---|
| `I:\Plasma immune atalas\.secrets\hf_token.txt` | Hugging Face, role **write** | https://huggingface.co/settings/tokens |
| `I:\Plasma immune atalas\.secrets\github_token.txt` | GitHub PAT, scope **repo** *(skip if `gh auth login` or Git Credential Manager already works)* | https://github.com/settings/tokens |

---

## Step 3 — publish

Double-click **`publish.bat`**, or:

```bash
cd "I:\Plasma immune atalas\Huggeface"
python publish.py
```

What it does:
1. re-runs `build_data.py` → refreshes `data/` from the latest analysis results
2. `git add -A` + commit with a timestamped message
3. pushes to GitHub, then to the Hugging Face Space
4. strips the tokens back out of `.git/config` so they are never stored

The Space rebuilds and restarts itself automatically, ~2–4 minutes. Then it is live at

> **https://huggingface.co/spaces/jianlizhao/Human-Plasma-Immune-Atlas**

Anyone can open it, pick a gene or disease, press **Run**, see the figures and download the CSV.
No Hugging Face account, no login, no upload.

---

## Auto-update from here on

**Manual, one command** — after any new analysis result upstream:

```bash
python publish.py                      # rebuild data + push both
python publish.py --no-build           # code/README changes only
python publish.py --only hf            # push to the Space only
python publish.py -m "add pQTL layer"  # custom message
```

**Automatic mirroring** — `.github/workflows/sync-to-huggingface.yml` pushes GitHub → Space on
every commit to `main`. Enable it once:

> GitHub repo → **Settings → Secrets and variables → Actions → New repository secret**
> name `HF_TOKEN`, value = your Hugging Face **write** token.

After that you can edit files directly on github.com and the Space redeploys by itself.

**Scheduled auto-update** (optional) — to refresh the app every week from the local analysis
folder, register a Windows scheduled task:

```powershell
schtasks /create /tn "PlasmaImmuneAtlas-Publish" /tr "\"I:\Plasma immune atalas\Huggeface\publish.bat\"" /sc weekly /d MON /st 09:00
```

---

## Making a change

| You want to… | Do this |
|---|---|
| add a tab / change a figure | edit `app.py`, test with `python app.py` → http://localhost:7860, then `python publish.py --no-build` |
| ship new analysis results | put the new `.tsv` in `../06_genetic_causality`, add it to the list in `build_data.py`, then `python publish.py` |
| change the Space title, emoji, gradio version | edit the YAML header at the top of `README.md` |
| add a python package | add it to `requirements.txt` |

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `remote: error: File data/mr_phenome_all.parquet is 53.00 MB; this exceeds…` | git-lfs is not active — `git lfs install`, `git lfs track "*.parquet"`, re-commit |
| HF push says `401` / `403` | token is read-only or expired → make a new **write** token |
| HF push says repository not found | the Space does not exist yet → create it (Step 1) with exactly the same name |
| Space build fails on `gradio` | the version in `README.md` `sdk_version:` must be a real release (currently `5.49.1`) |
| Space runs but shows "data not found" | `data/` was not pushed — check `git lfs ls-files` lists the parquet files |
| Space is slow on first click | free CPU hardware loads 1.66 M rows at boot; it sleeps after 48 h idle and wakes on the next visit |

---

## Local test before publishing

```bash
cd "I:\Plasma immune atalas\Huggeface"
pip install -r requirements.txt gradio==5.49.1
python app.py            # http://localhost:7860
```

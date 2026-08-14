#!/usr/bin/env python
"""
publish.py — one command to update and publish the Human Plasma Immune Atlas.

    python publish.py                 # rebuild data + push to GitHub + Hugging Face
    python publish.py --no-build      # skip the data rebuild (code/README changes only)
    python publish.py --only hf       # push only to Hugging Face
    python publish.py --only github   # push only to GitHub
    python publish.py --message "..." # custom commit message

Credentials (never stored in the repo)
--------------------------------------
Hugging Face : env var  HF_TOKEN            or  ../.secrets/hf_token.txt
GitHub       : env var  GITHUB_TOKEN        or  ../.secrets/github_token.txt
               (or just have `gh auth login` / a git credential helper set up)

Create the Hugging Face token at  https://huggingface.co/settings/tokens  with
role **write**. Create the GitHub token at  https://github.com/settings/tokens
with the **repo** scope. Paste it into the file yourself — this script only reads it.
"""
from __future__ import annotations

import argparse
import datetime
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SECRETS = os.path.join(os.path.dirname(HERE), ".secrets")

GH_USER = "jimmyuab"
GH_REPO = "Human-Plasma-Immune-Atlas"
HF_USER = "jianlizhao"
HF_SPACE = "Human-Plasma-Immune-Atlas"

GH_URL_PLAIN = f"https://github.com/{GH_USER}/{GH_REPO}.git"
HF_URL_PLAIN = f"https://huggingface.co/spaces/{HF_USER}/{HF_SPACE}"


# The project lives on an external drive that does not record ownership, so git
# refuses to operate on it ("dubious ownership"). Passing the exception through the
# environment keeps it scoped to this script — the user's global git config is untouched.
GIT_ENV = dict(os.environ,
               GIT_CONFIG_COUNT="1",
               GIT_CONFIG_KEY_0="safe.directory",
               GIT_CONFIG_VALUE_0=HERE.replace("\\", "/"))


def sh(cmd, check=True, quiet=False, **kw):
    if not quiet:
        printable = cmd if isinstance(cmd, str) else " ".join(cmd)
        for tok in ("hf_", "ghp_", "github_pat_"):
            if tok in printable:
                printable = printable.split(tok)[0] + tok + "***"
        print(f"  $ {printable}")
    return subprocess.run(cmd, cwd=HERE, shell=isinstance(cmd, str), env=GIT_ENV,
                          check=check, text=True, capture_output=True, **kw)


def step(msg):
    print(f"\n=== {msg} ===")


def token(env_name: str, filename: str) -> str | None:
    v = os.environ.get(env_name)
    if v:
        return v.strip()
    p = os.path.join(SECRETS, filename)
    if os.path.exists(p):
        with open(p) as fh:
            v = fh.read().strip()
        return v or None
    return None


def ensure_repo():
    if not os.path.isdir(os.path.join(HERE, ".git")):
        step("initialising git repository")
        sh("git init -b main")
    # local-only identity so the machine's global config is untouched
    if not sh("git config user.email", check=False, quiet=True).stdout.strip():
        sh('git config user.email "atlas@plasma-immunome.local"')
        sh('git config user.name "Human Plasma Immune Atlas"')
    if not sh("git lfs version", check=False, quiet=True).stdout.strip():
        print("  ! git-lfs not found. Install it (https://git-lfs.com) — the 53 MB "
              "parquet file needs LFS on both GitHub and Hugging Face.")
    else:
        sh("git lfs install --local", check=False, quiet=True)
        sh("git lfs track", check=False, quiet=True)


def commit(message: str) -> bool:
    sh("git add -A")
    if not sh("git diff --cached --quiet", check=False, quiet=True).returncode:
        print("  nothing new to commit")
        return False
    sh(["git", "commit", "-m", message])
    return True


def set_remote(name: str, url: str):
    if sh(f"git remote get-url {name}", check=False, quiet=True).returncode:
        sh(["git", "remote", "add", name, url], quiet=True)
    else:
        sh(["git", "remote", "set-url", name, url], quiet=True)


def push(name: str, url_with_auth: str, plain_url: str):
    set_remote(name, url_with_auth)
    r = sh(f"git push -u {name} main", check=False)
    # restore the token-free URL so it is never persisted in .git/config
    set_remote(name, plain_url)
    if r.returncode:
        print(r.stdout or "")
        print(r.stderr or "")
        return False
    print(f"  pushed -> {plain_url}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-build", action="store_true", help="skip rebuilding ./data")
    ap.add_argument("--only", choices=["hf", "github"], help="push to one target only")
    ap.add_argument("--message", "-m", default=None)
    a = ap.parse_args()

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    msg = a.message or f"Update Human Plasma Immune Atlas ({stamp})"

    if not a.no_build:
        step("rebuilding ./data from the analysis project")
        r = subprocess.run([sys.executable, os.path.join(HERE, "build_data.py")], cwd=HERE)
        if r.returncode:
            sys.exit("data rebuild failed — fix that first")

    step("git")
    ensure_repo()
    commit(msg)

    ok = True
    if a.only != "hf":
        step("pushing to GitHub")
        gt = token("GITHUB_TOKEN", "github_token.txt")
        url = (f"https://{GH_USER}:{gt}@github.com/{GH_USER}/{GH_REPO}.git" if gt else GH_URL_PLAIN)
        if not gt:
            print("  no GITHUB_TOKEN — relying on your git credential helper / gh login")
        if not push("github", url, GH_URL_PLAIN):
            ok = False
            print(f"  ! GitHub push failed. Create the empty repo first:\n"
                  f"    https://github.com/new  ->  name it  {GH_REPO}")

    if a.only != "github":
        step("pushing to the Hugging Face Space")
        ht = token("HF_TOKEN", "hf_token.txt")
        if not ht:
            print("  ! no HF_TOKEN found.\n"
                  f"    1) create a WRITE token at https://huggingface.co/settings/tokens\n"
                  f"    2) save it to {os.path.join(SECRETS, 'hf_token.txt')}  (git-ignored)\n"
                  f"       or  export HF_TOKEN=hf_xxx\n"
                  f"    3) create the Space (SDK: Gradio, visibility: PUBLIC) at\n"
                  f"       https://huggingface.co/new-space  ->  name it  {HF_SPACE}")
            ok = False
        else:
            url = f"https://{HF_USER}:{ht}@huggingface.co/spaces/{HF_USER}/{HF_SPACE}"
            if not push("hf", url, HF_URL_PLAIN):
                ok = False
                print(f"  ! Hugging Face push failed. Make sure the Space exists and is PUBLIC:\n"
                      f"    https://huggingface.co/new-space  ->  name  {HF_SPACE}, SDK Gradio")

    print("\n" + "=" * 62)
    print("DONE" if ok else "FINISHED WITH WARNINGS — see the messages above")
    print(f"  GitHub : {GH_URL_PLAIN[:-4]}")
    print(f"  Space  : {HF_URL_PLAIN}")
    print("  The Space rebuilds itself automatically after each push (~2-4 min).")
    print("=" * 62)


if __name__ == "__main__":
    main()

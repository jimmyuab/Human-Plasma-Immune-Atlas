#!/usr/bin/env bash
# Human Plasma Immune Atlas — one-click update and publish.
# Rebuilds ./data from the analysis project, then pushes to GitHub + Hugging Face.
set -e
cd "$(dirname "$0")"
python publish.py "$@"

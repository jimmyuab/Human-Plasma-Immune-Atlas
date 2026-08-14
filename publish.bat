@echo off
REM ===========================================================
REM  Human Plasma Immune Atlas - one-click update and publish
REM  Double-click this file. It rebuilds ./data from the
REM  analysis project, then pushes to GitHub + Hugging Face.
REM ===========================================================
cd /d "%~dp0"
python publish.py %*
echo(
pause

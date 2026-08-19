#!/bin/bash
# Launch this R/Shiny app locally; it opens in your web browser. Double-click me.
cd "$(dirname "$0")"
RS=""
[ -n "$R_RSCRIPT" ] && [ -x "$R_RSCRIPT" ] && RS="$R_RSCRIPT"
[ -z "$RS" ] && command -v Rscript >/dev/null 2>&1 && RS="$(command -v Rscript)"
[ -z "$RS" ] && [ -x /Library/Frameworks/R.framework/Resources/bin/Rscript ] && RS=/Library/Frameworks/R.framework/Resources/bin/Rscript
for C in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/mambaforge" /opt/miniconda3 /opt/anaconda3 /opt/homebrew/Caskroom/miniforge/base; do
  [ -z "$RS" ] && [ -x "$C/bin/Rscript" ] && RS="$C/bin/Rscript"
  [ -z "$RS" ] && [ -x "$C/envs/r/bin/Rscript" ] && RS="$C/envs/r/bin/Rscript"
done
if [ -z "$RS" ]; then
  echo "Could not find R (Rscript). Install R from https://cran.r-project.org/,"
  echo "or set R_RSCRIPT to your Rscript (e.g. inside a conda env) and re-run."
  read -r -p "Press Enter to close..." _ ; exit 1
fi
echo "Starting — a browser window will open shortly. (First run installs R packages.)"
"$RS" run_app.R

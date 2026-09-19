#!/usr/bin/env bash
# Build the precipitate specimen, shear it, and render the comparison.
# The clear-path half of the comparison is reused from ../shear_demo/out/prod,
# so run that demo first (or point --clean-dir somewhere else).
set -euo pipefail
cd "$(dirname "$0")"

RATE=${RATE:-1e9}
GMAX=${GMAX:-0.22}
FRAMES=${FRAMES:-240}
RADIUS=${RADIUS:-12.0}
XC=${XC:-12.0}

mkdir -p out
python3 src/build_precip.py --radius "$RADIUS" --xc "$XC" --out out/edge_precip.data

python3 ../shear_demo/src/run_shear.py edge \
    --data out/edge_precip.data --freeze-type 2 --tag edge_precip \
    --temp 300 --rate "$RATE" --gamma-max "$GMAX" --frames "$FRAMES" \
    --outdir out > out/edge_precip.run.log 2>&1

python3 src/make_summary_pair.py --rate "$RATE" --diameter "$(python3 -c "print(2*$RADIUS)")"
python3 src/render_pair.py --rate "$RATE"

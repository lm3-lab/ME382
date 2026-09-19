#!/usr/bin/env bash
# Reproduce the whole demonstration from scratch.
#   ./run_all.sh            production settings (~80 min on 4 cores)
#   ./run_all.sh quick      coarse preview (~6 min)
set -euo pipefail
cd "$(dirname "$0")"

if [[ "${1:-}" == "quick" ]]; then
  RATE=5e9; GMAX=0.18; FRAMES=90; OUT=out/quick
else
  RATE=1e9; GMAX=0.22; FRAMES=240; OUT=out/prod
fi

mkdir -p "$OUT"
for case in perfect edge screw; do
  python3 src/run_shear.py "$case" --temp 300 --rate "$RATE" \
      --gamma-max "$GMAX" --frames "$FRAMES" --outdir "$OUT" \
      > "$OUT/$case.run.log" 2>&1 &
done
wait

python3 src/make_summary.py  --outdir "$OUT" --rate "$RATE" \
        --png "$OUT/stress_strain.png"
python3 src/render_movie.py  --outdir "$OUT" --rate "$RATE" \
        --movie "$OUT/shear_dislocation_demo.mp4"

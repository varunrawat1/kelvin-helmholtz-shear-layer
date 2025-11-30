#!/usr/bin/env bash
set -e

# Find latest run directory in output/ or build/output/
run_dir=$(ls -d output/run_* build/output/run_* 2>/dev/null | sort | tail -n 1)
if [ -z "$run_dir" ]; then
  echo "No run_* directory found in output/ or build/output/"
  exit 1
fi

echo "Using run directory: $run_dir"

cd "$run_dir"

if ! ls vorticity_*.png >/dev/null 2>&1; then
  echo "No PNGs found. Run: python3 ../../scripts/plot_latest_run.py"
  exit 1
fi

convert -delay 22 -loop 0 vorticity_*.png kh_vorticity.gif
echo "Created $run_dir/kh_vorticity.gif"
EOF

chmod +x scripts/make_gif_latest_run.sh
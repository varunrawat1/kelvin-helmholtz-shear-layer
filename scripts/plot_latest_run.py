import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm


def find_latest_run():
    candidates = sorted(glob.glob("output/run_*") +
                        glob.glob("build/output/run_*"))
    if not candidates:
        raise SystemExit("No run_* directory found in output/ or build/output/")
    return candidates[-1]

def main():
    run_dir = find_latest_run()
    print("Plotting vorticity snapshots from:", run_dir)

    files = sorted(glob.glob(os.path.join(run_dir, "vorticity_*.dat")))
    if not files:
        raise SystemExit("No vorticity_*.dat files found in " + run_dir)

    max_abs = 0.0
    for fname in files:
        w = np.loadtxt(fname)
        max_abs = max(max_abs, np.max(np.abs(w)))
    vmax = max_abs
    vmin = -max_abs
    print(f"Using symmetric colour limits vmin={vmin:.3g}, vmax={vmax:.3g}")

    for fname in files:
        omega = np.loadtxt(fname)
        fig, ax = plt.subplots(figsize=(5, 4), dpi=150)

        im = ax.imshow(
        omega,
        origin="lower",
        cmap="coolwarm",
        extent=[0.0, 1.0, 0.0, 1.0],
        interpolation="bilinear",
        norm=SymLogNorm(linthresh=0.1 * vmax, vmin=vmin, vmax=vmax),
        )


        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.tick_params(axis="both", labelsize=8)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"$\omega$", fontsize=9)

        fig.tight_layout()

        png_name = fname.replace(".dat", ".png")
        fig.savefig(png_name)
        plt.close(fig)
        print("Saved", png_name)

if __name__ == "__main__":
    main()

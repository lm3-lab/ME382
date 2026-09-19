"""Static stress-strain summary figure (for slides) + printed numbers."""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, PANEL = "#1a1a19", "#232321"
TEXT_1, TEXT_2, TEXT_3, GRID = "#ffffff", "#c3c2b7", "#8a897f", "#33332f"
SERIES = {"perfect": "#3987e5", "edge": "#d95926", "screw": "#199e70"}
LABEL = {"perfect": "perfect crystal (no dislocation)",
         "edge": "one edge dislocation",
         "screw": "one screw dislocation"}
CASES = ["perfect", "edge", "screw"]


def flow_stress(g, t, gmin):
    """Median stress once the specimen is flowing, i.e. past gmin."""
    m = g > gmin
    return float(np.median(t[m])) if m.any() else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="out/prod")
    ap.add_argument("--png", default="out/stress_strain.png")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--rate", type=float, default=1.0e9)
    args = ap.parse_args()

    d = {c: np.loadtxt(os.path.join(args.outdir, f"{c}.stress.txt"), skiprows=2)
         for c in CASES}
    g0, t0 = d["perfect"][:, 1], d["perfect"][:, 2]
    ipk = int(np.argmax(t0))
    tau_ideal, gam_ideal = t0[ipk], g0[ipk]
    el = g0 < 0.6 * gam_ideal
    G = np.polyfit(g0[el], t0[el], 1)[0]

    fig, ax = plt.subplots(figsize=(11, 6.6), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(PANEL)
    stats = {}
    for c in CASES:
        g, t = d[c][:, 1], d[c][:, 2]
        ax.plot(g, t, color=SERIES[c], lw=2.0, label=LABEL[c],
                solid_capstyle="round")
        stats[c] = flow_stress(g, t, max(gam_ideal, 0.05))

    ax.axhline(tau_ideal, color=SERIES["perfect"], lw=1.1, ls=(0, (5, 4)), alpha=0.7)
    ax.annotate(f"ideal shear strength  τₘₐₓ = {tau_ideal:.2f} GPa",
                xy=(g0.max() * 0.995, tau_ideal * 1.015), color=TEXT_2,
                fontsize=11.5, ha="right", va="bottom")
    # both flow stresses sit almost on the axis, so label them with leaders
    # reaching into the empty space under the elastic branch
    for c, ty in (("screw", 0.32), ("edge", 0.18)):
        ax.annotate(f"{c}: {stats[c]:.2f} GPa  ({tau_ideal/stats[c]:.0f}× lower)",
                    xy=(g0.max() * 0.66, stats[c]),
                    xytext=(g0.max() * 0.34, tau_ideal * ty),
                    color=SERIES[c], fontsize=11.5, fontweight="bold",
                    va="center", ha="left",
                    arrowprops=dict(arrowstyle="-", color=SERIES[c], lw=1.0,
                                    alpha=0.55,
                                    connectionstyle="angle,angleA=0,angleB=90,rad=6"))


    ax.set_xlabel("shear strain  γ", color=TEXT_2, fontsize=12.5)
    ax.set_ylabel("resolved shear stress  τ  (GPa)", color=TEXT_2, fontsize=12.5)
    ax.tick_params(colors=TEXT_3, labelsize=11)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.set_xlim(0, g0.max())
    ax.set_ylim(0, tau_ideal * 1.15)
    ax.legend(loc="upper left", bbox_to_anchor=(0.015, 0.88), frameon=False,
              fontsize=11.5, labelcolor=TEXT_2)
    ax.set_title("BCC iron sheared on {110}⟨111⟩ — "
                 f"{args.temp:.0f} K, {args.rate:.0e} s⁻¹",
                 color=TEXT_1, fontsize=15, fontweight="bold", pad=14, loc="left")
    fig.tight_layout()
    fig.savefig(args.png, facecolor=SURFACE)

    print(f"elastic shear modulus  G          = {G:6.1f} GPa")
    print(f"ideal shear strength   tau_max    = {tau_ideal:6.2f} GPa "
          f"at gamma = {gam_ideal:.3f}   (= G/{G/tau_ideal:.0f})")
    for c in ("edge", "screw"):
        print(f"flow stress, {c:<7s}  tau_flow   = {stats[c]:6.2f} GPa "
              f"-> {tau_ideal/stats[c]:5.1f}x weaker than the perfect lattice")
    print("wrote", args.png)


if __name__ == "__main__":
    main()

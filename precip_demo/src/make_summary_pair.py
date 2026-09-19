"""Static comparison figure + the numbers, for the precipitate study."""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, PANEL = "#1a1a19", "#232321"
TEXT_1, TEXT_2, TEXT_3, GRID = "#ffffff", "#c3c2b7", "#8a897f", "#33332f"
CLEAN, PINNED = "#d95926", "#3987e5"

G_MEASURED = 68.7      # GPa, from shear_demo/src/check_modulus.py
B = 2.4728             # A, |a/2[111]|


def orowan_bks(G, b, D, L):
    """Bacon-Kocks-Scattergood critical stress for a periodic row of strong,
    impenetrable particles.  L is the free channel between particle surfaces,
    D the particle diameter, Dbar their harmonic mean."""
    Dbar = D * L / (D + L)
    return (G * b) / (2 * np.pi * L) * (np.log(Dbar / b) + 0.7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean", default="../shear_demo/out/prod/edge.stress.txt")
    ap.add_argument("--precip", default="out/edge_precip.stress.txt")
    ap.add_argument("--png", default="out/precipitate_stress_strain.png")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--rate", type=float, default=1.0e9)
    ap.add_argument("--diameter", type=float, default=24.0)
    ap.add_argument("--lz", type=float, default=48.959)
    ap.add_argument("--gamma-flow", type=float, default=0.05)
    args = ap.parse_args()

    dc = np.loadtxt(args.clean, skiprows=2)
    dp = np.loadtxt(args.precip, skiprows=2)
    n = min(len(dc), len(dp))
    dc, dp = dc[:n], dp[:n]

    def flow(d):
        m = d[:, 1] > args.gamma_flow
        return float(np.median(d[m, 2])), float(np.max(d[m, 2]))

    tc, pc = flow(dc)
    tp, pp = flow(dp)
    channel = args.lz - args.diameter
    t_or = orowan_bks(G_MEASURED, B, args.diameter, channel)

    fig, ax = plt.subplots(figsize=(11, 6.6), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(PANEL)
    ax.plot(dc[:, 1], dc[:, 2], color=CLEAN, lw=1.8,
            label="edge dislocation, clear path", solid_capstyle="round")
    ax.plot(dp[:, 1], dp[:, 2], color=PINNED, lw=1.8,
            label="same dislocation, one 2.4 nm hard particle",
            solid_capstyle="round")

    ax.axhline(tc, color=CLEAN, lw=1.0, ls=(0, (5, 4)), alpha=0.65)
    ax.axhline(tp, color=PINNED, lw=1.0, ls=(0, (5, 4)), alpha=0.65)
    ax.axhline(t_or, color=TEXT_3, lw=1.0, ls=(0, (2, 3)), alpha=0.8)

    ymax = max(pp, t_or) * 1.18
    ax.annotate(f"flow stress with the particle:  {tp:.2f} GPa",
                xy=(dp[-1, 1] * 0.99, tp), xytext=(dp[-1, 1] * 0.99, tp + ymax * 0.045),
                color=PINNED, fontsize=11.5, fontweight="bold", ha="right")
    ax.annotate(f"clear path:  {tc:.2f} GPa",
                xy=(dc[-1, 1] * 0.99, tc), xytext=(dc[-1, 1] * 0.99, tc + ymax * 0.035),
                color=CLEAN, fontsize=11.5, fontweight="bold", ha="right")
    ax.annotate(f"Orowan estimate (Bacon–Kocks–Scattergood)  {t_or:.2f} GPa",
                xy=(dp[-1, 1] * 0.02, t_or), xytext=(dp[-1, 1] * 0.02, t_or + ymax * 0.035),
                color=TEXT_2, fontsize=11, ha="left")

    ax.set_xlabel("shear strain  γ", color=TEXT_2, fontsize=12.5)
    ax.set_ylabel("resolved shear stress  τ  (GPa)", color=TEXT_2, fontsize=12.5)
    ax.tick_params(colors=TEXT_3, labelsize=11)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.set_xlim(0, dp[-1, 1])
    ax.set_ylim(0, ymax)
    ax.legend(loc="upper left", bbox_to_anchor=(0.015, 0.99), frameon=False,
              fontsize=11.5, labelcolor=TEXT_2)
    ax.set_title("One hard particle in the glide path — BCC iron, "
                 f"{args.temp:.0f} K, {args.rate:.0e} s⁻¹",
                 color=TEXT_1, fontsize=15, fontweight="bold", pad=14, loc="left")
    fig.tight_layout()
    fig.savefig(args.png, facecolor=SURFACE)

    print(f"particle            d = {args.diameter:.1f} A, "
          f"free channel L = {channel:.1f} A")
    print(f"flow stress, clear path      = {tc:6.2f} GPa   (peak {pc:.2f})")
    print(f"flow stress, with particle   = {tp:6.2f} GPa   (peak {pp:.2f})")
    print(f"strengthening increment      = {tp-tc:6.2f} GPa "
          f"-> {tp/tc:.0f}x the unobstructed flow stress")
    print(f"Orowan estimate (BKS)        = {t_or:6.2f} GPa "
          f"({100*tp/t_or:.0f} % of it measured)")
    print("wrote", args.png)


if __name__ == "__main__":
    main()

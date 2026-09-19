"""
Render the side-by-side shear movie: three 3D atomistic views over one live
resolved-shear-stress / shear-strain chart.

The three specimens are the same block of BCC Fe under the same loading; only
the dislocation they contain differs.  The chart is what carries the message --
the perfect crystal climbs to the ideal shear strength, the two crystals that
already contain one dislocation yield an order of magnitude lower.
"""

import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dumpio import read_frames

# ---------------------------------------------------------------- palette ---
SURFACE   = "#ffffff"
PANEL     = "#ffffff"
TEXT_1    = "#0b0b0b"
TEXT_2    = "#52514e"
TEXT_3    = "#78776f"
GRID      = "#e3e2dd"
BULK      = "#c9c8c1"          # ghost lattice (off by default)
GRIP      = "#b0afa5"          # the loading platens
FRAME_REF = "#c9c8c1"          # undeformed supercell
FRAME_DEF = "#6f6e66"          # the supercell as the grips have sheared it
PARTICLE  = "#7d7c72"          # an inert obstacle: a neutral, not a series hue
# categorical slots 1-3, light-mode steps (validated all-pairs)
SERIES = {"perfect": "#2a78d6", "edge": "#eb6834", "screw": "#1baf7a"}

LABEL = {
    "perfect": "perfect crystal",
    "edge":    "edge dislocation",
    "screw":   "screw dislocation",
}
SUBLABEL = {
    "perfect": "no dislocation present",
    "edge":    "line ⊥ b, glides along x",
    "screw":   "line ∥ b, glides along z",
}
CASES = ["perfect", "edge", "screw"]


# ------------------------------------------------------------ projection ---
def make_view(theta_deg=30.0, phi_deg=16.0):
    """Orthographic camera: cell +y is screen-up, view looks down +z tilted."""
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    view = np.array([np.sin(t) * np.cos(p), np.sin(p), np.cos(t) * np.cos(p)])
    view /= np.linalg.norm(view)
    right = np.cross([0.0, 1.0, 0.0], view)
    right /= np.linalg.norm(right)
    up = np.cross(view, right)
    return np.vstack([right, up, view])          # rows: right, up, depth


def project(p, M, centre):
    q = (p - centre) @ M.T
    return q[:, 0], q[:, 1], q[:, 2]


BOX_EDGES = [(0, 1), (1, 3), (3, 2), (2, 0), (4, 5), (5, 7), (7, 6), (6, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]


def edge_points(corners, n=24):
    """Box edges as one polyline, NaN-separated, subdivided so that the
    piecewise shear below shows up as a visible bend."""
    segs = []
    for i, j in BOX_EDGES:
        a, b = corners[i], corners[j]
        t = np.linspace(0.0, 1.0, n)[:, None]
        segs.append(a + (b - a) * t)
        segs.append(np.full((1, 3), np.nan))
    return np.vstack(segs)


def shear_frame(pts, gamma, Ly, slab, h):
    """Apply the imposed grip motion to a set of points.

    The lower grip is held, the upper grip has translated by gamma*h, and the
    crystal between them is sheared uniformly -- so the drawn cell shows the
    applied strain directly, with the two rigid grips staying square.
    """
    f = np.clip((pts[:, 1] - slab) / max(Ly - 2.0 * slab, 1e-9), 0.0, 1.0)
    out = pts.copy()
    out[:, 0] = out[:, 0] + gamma * h * f
    return out


def box_corners(L):
    return np.array([[i * L[0], j * L[1], k * L[2]]
                     for i in (0, 1) for j in (0, 1) for k in (0, 1)])


# ------------------------------------------------------------------ data ---
def load_case(outdir, case, slab, stride, rng=None):
    traj = os.path.join(outdir, f"{case}.shear.dump")
    frames = read_frames(traj, fields={"x", "y", "z", "c_cna", "type"})
    stress = np.loadtxt(os.path.join(outdir, f"{case}.stress.txt"), skiprows=2)
    packed = []
    for fr in frames:
        x, y, z, cna = fr["x"], fr["y"], fr["z"], fr["c_cna"]
        # atom type 2, where present, is a rigid second-phase particle
        part = (fr["type"] == 2) if "type" in fr else np.zeros(len(x), bool)
        ylo, yhi = y.min(), y.max()
        grip = (y < ylo + slab) | (y > yhi - slab)
        defect = (~grip) & (~part) & (cna != 3)
        bulk = (~grip) & (~part) & (~defect)
        # a *random* subsample -- taking every n-th atom of an id-sorted list
        # samples lattice planes and produces heavy moire fringes
        bi = np.where(bulk)[0]
        bi = rng.choice(bi, size=max(1, len(bi) // stride), replace=False)
        gi = np.where(grip)[0]
        gi = rng.choice(gi, size=max(1, len(gi) // max(1, stride // 2)),
                        replace=False)
        packed.append(dict(
            step=fr["step"],
            bulk=np.stack([x[bi], y[bi], z[bi]], 1),
            grip=np.stack([x[gi], y[gi], z[gi]], 1),
            defect=np.stack([x[defect], y[defect], z[defect]], 1),
            precip=np.stack([x[part], y[part], z[part]], 1),
            ndef=int(defect.sum())))
    return packed, stress


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="out/prod")
    ap.add_argument("--movie", default="out/shear_dislocation_demo.mp4")
    ap.add_argument("--slab", type=float, default=8.0)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--dpi", type=int, default=100)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--rate", type=float, default=1.0e9)
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--bulk", action="store_true",
                    help="also draw the (ghost) bulk lattice; off by default so "
                         "the dislocation and the supercell frame stand alone")
    args = ap.parse_args()

    data, stress = {}, {}
    for c in CASES:
        data[c], stress[c] = load_case(args.outdir, c, args.slab, args.stride,
                                       np.random.default_rng(7))
        print(f"  {c}: {len(data[c])} frames, {len(stress[c])} stress points", flush=True)

    nf = min(len(data[c]) for c in CASES)
    if args.max_frames:
        nf = min(nf, args.max_frames)

    # match each dump frame to its stress row by timestep
    curves = {}
    for c in CASES:
        s = stress[c]
        step2row = {int(r[0]): (r[1], r[2]) for r in s}
        g, t = [], []
        for fr in data[c][:nf]:
            gv, tv = step2row.get(fr["step"],
                                  (0.0, 0.0) if fr["step"] == 0 else (np.nan, np.nan))
            g.append(gv); t.append(tv)
        curves[c] = (np.array(g), np.array(t))

    gmax = np.nanmax([curves[c][0][-1] for c in CASES])
    tmax = np.nanmax([np.nanmax(curves[c][1]) for c in CASES])

    # geometry, fixed for every frame so the three panels share one scale
    L = np.array([data["perfect"][0]["bulk"][:, k].max() for k in range(3)])
    M = make_view()
    centre = L / 2.0
    corners = box_corners(L)
    hgrip = L[1] - args.slab                  # grip centre-to-centre, as loaded
    sheared = shear_frame(corners, gmax, L[1], args.slab, hgrip)
    cx, cy, _ = project(np.vstack([corners, sheared]), M, centre)
    pad = 0.04 * max(np.ptp(cx), np.ptp(cy))
    xlim = [cx.min() - pad, cx.max() + pad]
    ylim = [cy.min() - pad, cy.max() + pad]

    # ------------------------------------------------------------- figure --
    fig = plt.figure(figsize=(16, 9), dpi=args.dpi, facecolor=SURFACE)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.3, 1.0],
                  left=0.045, right=0.985, top=0.845, bottom=0.085,
                  wspace=0.045, hspace=0.20)
    axes3d = [fig.add_subplot(gs[0, i]) for i in range(3)]
    axc = fig.add_subplot(gs[1, :])

    fig.text(0.045, 0.962,
             "What a single dislocation does to the strength of a crystal",
             color=TEXT_1, fontsize=23, fontweight="bold", ha="left", va="center")
    fig.text(0.045, 0.921,
             f"BCC iron sheared on the {{110}}\u27e8111\u27e9 slip system at {args.temp:.0f} K, "
             f"strain rate {args.rate:.0e} s⁻¹ — identical blocks, "
             "identical loading, one dislocation apart",
             color=TEXT_2, fontsize=12.5, ha="left", va="center")

    # expand the shorter data axis so the projected cell fills the panel with
    # no letterboxing, while keeping the projection undistorted
    fig.canvas.draw()
    bb = axes3d[0].get_window_extent()
    target = bb.width / bb.height
    dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]
    if dx / dy < target:
        extra = (dy * target - dx) / 2.0
        xlim = [xlim[0] - extra, xlim[1] + extra]
    else:
        extra = (dx / target - dy) / 2.0
        ylim = [ylim[0] - extra, ylim[1] + extra]

    for ax, c in zip(axes3d, CASES):
        ax.set_facecolor(PANEL)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
        for sp in ax.spines.values():
            sp.set_color(GRID)
        # the loading: top grip slides along +x = b
        ax.annotate("", xy=(0.80, 0.955), xytext=(0.52, 0.955),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=TEXT_3, lw=1.6,
                                    shrinkA=0, shrinkB=0))
        ax.text(0.81, 0.955, "grip", transform=ax.transAxes, color=TEXT_3,
                fontsize=9.5, va="center", ha="left")

    axc.set_facecolor(PANEL)
    axc.set_xlim(0, gmax * 1.14)
    axc.set_ylim(0, tmax * 1.16)
    axc.set_xlabel("shear strain  γ", color=TEXT_2, fontsize=12.5)
    axc.set_ylabel("resolved shear stress  τ  (GPa)", color=TEXT_2, fontsize=12.5)
    axc.tick_params(colors=TEXT_3, labelsize=11)
    axc.grid(True, color=GRID, lw=0.8, alpha=0.9)
    axc.set_axisbelow(True)
    for s in ("top", "right"):
        axc.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        axc.spines[s].set_color(GRID)

    lines, heads, htexts = {}, {}, {}
    for c in CASES:
        lines[c], = axc.plot([], [], color=SERIES[c], lw=2.0, solid_capstyle="round")
        heads[c], = axc.plot([], [], "o", color=SERIES[c], ms=8,
                             mec=PANEL, mew=1.6)
        htexts[c] = axc.text(0, 0, "", color=SERIES[c], fontsize=11.5,
                             fontweight="bold", va="center", ha="left")
    axc.legend(handles=[Line2D([], [], color=SERIES[c], lw=2.6, label=LABEL[c])
                        for c in CASES],
               loc="upper left", frameon=False, fontsize=11.5,
               labelcolor=TEXT_2, ncol=3, handlelength=1.8, columnspacing=1.8)

    peak_line = axc.axhline(np.nan, color=SERIES["perfect"], lw=1.2, ls=(0, (5, 4)),
                            alpha=0.0)
    peak_txt = axc.text(0, 0, "", color=TEXT_2, fontsize=11, va="bottom", ha="right")
    caption = fig.text(0.045, 0.022, "", color=TEXT_3, fontsize=11.5, ha="left")

    panel_ttl, panel_sub, panel_val, arts = [], [], [], []
    for ax, c in zip(axes3d, CASES):
        panel_ttl.append(ax.text(0.5, 1.075, LABEL[c], transform=ax.transAxes,
                                 color=SERIES[c], fontsize=15, fontweight="bold",
                                 ha="center", va="center"))
        panel_sub.append(ax.text(0.5, 1.022, SUBLABEL[c], transform=ax.transAxes,
                                 color=TEXT_3, fontsize=10.5, ha="center", va="center"))
        panel_val.append(ax.text(0.03, 0.955, "", transform=ax.transAxes,
                                 color=TEXT_2, fontsize=11.5, ha="left", va="top",
                                 family="monospace"))
        arts.append({})

    # the undeformed supercell, drawn once, plus the sheared cell which is
    # redrawn every frame -- together they show the applied strain
    ref_pts = edge_points(corners)
    rx, ry, _ = project(ref_pts, M, centre)
    cell_art = []
    for ax in axes3d:
        ax.plot(rx, ry, color=FRAME_REF, lw=1.0, ls=(0, (4, 3)), zorder=0)
        cell_art.append(ax.plot([], [], color=FRAME_DEF, lw=1.5, zorder=1)[0])

    import imageio_ffmpeg
    writer = imageio_ffmpeg.write_frames(
        args.movie, (int(16 * args.dpi), int(9 * args.dpi)), fps=args.fps,
        quality=8, macro_block_size=1,
        output_params=["-movflags", "+faststart"])
    writer.send(None)

    peak_found = {"done": False, "tau": 0.0, "gam": 0.0}

    for k in range(nf):
        for ax, c, ttl, val, store in zip(axes3d, CASES, panel_ttl, panel_val, arts):
            fr = data[c][k]
            empty = np.empty((0, 3))
            groups = [fr["grip"], fr["bulk"] if args.bulk else empty,
                      fr["precip"], fr["defect"]]
            pts = np.vstack(groups)
            kind = np.concatenate([np.full(len(g), i) for i, g in enumerate(groups)])
            sx, sy, sd = project(pts, M, centre)
            order = np.argsort(-sd)
            sx, sy, kind = sx[order], sy[order], kind[order]
            col = np.array([GRIP, BULK, PARTICLE, SERIES[c]], dtype=object)[kind]
            siz = np.array([3.4, 2.0, 9.0, 26.0])[kind]
            alp = np.array([0.70, 0.45, 0.95, 1.0])[kind]

            for a in store.values():
                a.remove()
            store.clear()
            store["cloud"] = ax.scatter(sx, sy, s=siz, c=list(col), alpha=None,
                                        linewidths=0, zorder=2,
                                        edgecolors="none")
            store["cloud"].set_alpha(None)
            rgba = store["cloud"].get_facecolor()
            rgba[:, 3] = alp
            store["cloud"].set_facecolor(rgba)

            dsel = kind == 3
            if dsel.any():
                store["def"] = ax.scatter(sx[dsel], sy[dsel], s=26,
                                          c=SERIES[c], linewidths=0.4,
                                          edgecolors=PANEL, zorder=3)
            g, t = curves[c]
            dpts = shear_frame(ref_pts, g[k], L[1], args.slab, hgrip)
            dx_, dy_, _ = project(dpts, M, centre)
            cell_art[CASES.index(c)].set_data(dx_, dy_)
            val.set_text(f"τ = {t[k]:5.2f} GPa")

        # ---- chart
        for c in CASES:
            g, t = curves[c]
            lines[c].set_data(g[:k + 1], t[:k + 1])
            heads[c].set_data([g[k]], [t[k]])
            htexts[c].set_position((g[k] + gmax * 0.012, t[k]))
            htexts[c].set_text(f"{t[k]:.2f}")

        gp, tp = curves["perfect"]
        if not peak_found["done"] and k > 6:
            i = int(np.nanargmax(tp[:k + 1]))
            if tp[k] < 0.8 * tp[i] and i < k - 1:
                peak_found.update(done=True, tau=float(tp[i]), gam=float(gp[i]))
                peak_line.set_ydata([peak_found["tau"]] * 2)
                peak_line.set_alpha(0.75)
                peak_txt.set_position((gmax * 1.13, peak_found["tau"] * 1.012))
                peak_txt.set_text(
                    f"ideal shear strength  τₘₐₓ = "
                    f"{peak_found['tau']:.2f} GPa")
        if peak_found["done"]:
            fe = np.nanmedian(curves["edge"][1][max(0, k - 20):k + 1])
            caption.set_text(
                f"the perfect lattice had to shear every bond at once — "
                f"{peak_found['tau']:.2f} GPa;  one dislocation lets the same "
                f"crystal flow at ~{fe:.2f} GPa, a factor of "
                f"{peak_found['tau']/max(fe,1e-6):.0f} lower")
        else:
            caption.set_text("all three blocks climb the same elastic line — "
                             "watch where each one leaves it")

        writer.send(_rgb(fig))
        if k % 20 == 0:
            print(f"  frame {k+1}/{nf}", flush=True)

    writer.close()
    print("wrote", args.movie)


def _rgb(fig):
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return np.ascontiguousarray(buf[:, :, :3])


if __name__ == "__main__":
    main()

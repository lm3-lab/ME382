"""
Two-panel movie: the same edge dislocation, with and without a hard particle
in its glide path, over one shared resolved-shear-stress / shear-strain chart.

The left panel reuses the dislocation-only run from ../shear_demo verbatim --
same cell, same loading, same relaxed starting structure.  The only difference
on the right is that 612 atoms are held rigid.
"""

import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

_SHEAR_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "shear_demo", "src")
sys.path.insert(0, os.path.abspath(_SHEAR_SRC))
from render_movie import (SURFACE, PANEL, TEXT_1, TEXT_2, TEXT_3, GRID,   # noqa
                          BULK, GRIP, PARTICLE, FRAME_REF, FRAME_DEF,
                          make_view, project, box_corners, BOX_EDGES,
                          edge_points, shear_frame, load_case, _rgb)

CLEAN, PINNED = "#eb6834", "#2a78d6"     # light-mode categorical slots 2 and 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean-dir", default="../shear_demo/out/prod")
    ap.add_argument("--clean-case", default="edge")
    ap.add_argument("--precip-dir", default="out")
    ap.add_argument("--precip-case", default="edge_precip")
    ap.add_argument("--movie", default="out/precipitate_demo.mp4")
    ap.add_argument("--slab", type=float, default=8.0)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--dpi", type=int, default=100)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--rate", type=float, default=1.0e9)
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--bulk", action="store_true",
                    help="also draw the (ghost) bulk lattice")
    args = ap.parse_args()

    spec = [("clean", args.clean_dir, args.clean_case, CLEAN,
             "edge dislocation, clear path",
             "nothing in the way — it glides freely"),
            ("pinned", args.precip_dir, args.precip_case, PINNED,
             "same dislocation, one hard particle",
             "2.4 nm impenetrable particle on the glide plane")]

    data, curves_raw, colors, titles, subs, keys = {}, {}, {}, {}, {}, []
    rng = np.random.default_rng(7)
    for key, d, case, col, ttl, sub in spec:
        data[key], curves_raw[key] = load_case(d, case, args.slab, args.stride,
                                               np.random.default_rng(7))
        colors[key], titles[key], subs[key] = col, ttl, sub
        keys.append(key)
        print(f"  {key}: {len(data[key])} frames", flush=True)

    nf = min(len(data[k]) for k in keys)
    if args.max_frames:
        nf = min(nf, args.max_frames)

    curves = {}
    for k in keys:
        s = curves_raw[k]
        step2row = {int(r[0]): (r[1], r[2]) for r in s}
        g, t = [], []
        for fr in data[k][:nf]:
            gv, tv = step2row.get(fr["step"],
                                  (0.0, 0.0) if fr["step"] == 0 else (np.nan, np.nan))
            g.append(gv); t.append(tv)
        curves[k] = (np.array(g), np.array(t))

    gmax = np.nanmax([curves[k][0][-1] for k in keys])
    tmax = np.nanmax([np.nanmax(curves[k][1]) for k in keys])

    L = np.array([data["clean"][0]["bulk"][:, i].max() for i in range(3)])
    M = make_view()
    centre = L / 2.0
    corners = box_corners(L)
    hgrip = L[1] - args.slab
    sheared = shear_frame(corners, gmax, L[1], args.slab, hgrip)
    cx, cy, _ = project(np.vstack([corners, sheared]), M, centre)
    pad = 0.04 * max(np.ptp(cx), np.ptp(cy))
    xlim = [cx.min() - pad, cx.max() + pad]
    ylim = [cy.min() - pad, cy.max() + pad]

    fig = plt.figure(figsize=(16, 9), dpi=args.dpi, facecolor=SURFACE)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.3, 1.0],
                  left=0.05, right=0.98, top=0.845, bottom=0.085,
                  wspace=0.05, hspace=0.20)
    axes3d = [fig.add_subplot(gs[0, i]) for i in range(2)]
    axc = fig.add_subplot(gs[1, :])

    fig.text(0.05, 0.962, "How a precipitate strengthens a metal",
             color=TEXT_1, fontsize=23, fontweight="bold", ha="left", va="center")
    fig.text(0.05, 0.921,
             f"The same BCC-iron edge dislocation under the same shear at "
             f"{args.temp:.0f} K, {args.rate:.0e} s⁻¹ — "
             "the only difference is one hard particle in its path",
             color=TEXT_2, fontsize=12.5, ha="left", va="center")

    fig.canvas.draw()
    bb = axes3d[0].get_window_extent()
    target = bb.width / bb.height
    dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]
    if dx / dy < target:
        e = (dy * target - dx) / 2.0
        xlim = [xlim[0] - e, xlim[1] + e]
    else:
        e = (dx / target - dy) / 2.0
        ylim = [ylim[0] - e, ylim[1] + e]

    for ax, k in zip(axes3d, keys):
        ax.set_facecolor(PANEL)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
        for sp in ax.spines.values():
            sp.set_color(GRID)
        ax.annotate("", xy=(0.84, 0.955), xytext=(0.60, 0.955),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=TEXT_3, lw=1.6,
                                    shrinkA=0, shrinkB=0))
        ax.text(0.85, 0.955, "grip", transform=ax.transAxes, color=TEXT_3,
                fontsize=9.5, va="center", ha="left")
    ref_pts = edge_points(corners)
    rx, ry, _ = project(ref_pts, M, centre)
    cell_art = []
    for ax in axes3d:
        ax.plot(rx, ry, color=FRAME_REF, lw=1.0, ls=(0, (4, 3)), zorder=0)
        cell_art.append(ax.plot([], [], color=FRAME_DEF, lw=1.5, zorder=1)[0])

    axc.set_facecolor(PANEL)
    axc.set_xlim(0, gmax * 1.10)
    axc.set_ylim(0, tmax * 1.18)
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
    for k in keys:
        lines[k], = axc.plot([], [], color=colors[k], lw=2.0,
                             solid_capstyle="round")
        heads[k], = axc.plot([], [], "o", color=colors[k], ms=8,
                             mec=PANEL, mew=1.6)
        htexts[k] = axc.text(0, 0, "", color=colors[k], fontsize=11.5,
                             fontweight="bold", va="center", ha="left")
    axc.legend(handles=[Line2D([], [], color=colors[k], lw=2.6, label=titles[k])
                        for k in keys],
               loc="upper left", frameon=False, fontsize=11.5,
               labelcolor=TEXT_2, ncol=2, handlelength=1.8, columnspacing=2.4)

    caption = fig.text(0.05, 0.022, "", color=TEXT_3, fontsize=11.5, ha="left")
    panel_val, arts = [], []
    for ax, k in zip(axes3d, keys):
        ax.text(0.5, 1.075, titles[k], transform=ax.transAxes, color=colors[k],
                fontsize=15, fontweight="bold", ha="center", va="center")
        ax.text(0.5, 1.022, subs[k], transform=ax.transAxes, color=TEXT_3,
                fontsize=10.5, ha="center", va="center")
        panel_val.append(ax.text(0.03, 0.955, "", transform=ax.transAxes,
                                 color=TEXT_2, fontsize=11.5, ha="left",
                                 va="top", family="monospace"))
        arts.append({})

    import imageio_ffmpeg
    writer = imageio_ffmpeg.write_frames(
        args.movie, (int(16 * args.dpi), int(9 * args.dpi)), fps=args.fps,
        quality=8, macro_block_size=1)
    writer.send(None)

    for n in range(nf):
        for ax, k, val, store in zip(axes3d, keys, panel_val, arts):
            fr = data[k][n]
            groups = [fr["grip"], fr["bulk"] if args.bulk else np.empty((0, 3)),
                      fr["precip"], fr["defect"]]
            pts = np.vstack(groups)
            kind = np.concatenate([np.full(len(g), i) for i, g in enumerate(groups)])
            sx, sy, sd = project(pts, M, centre)
            order = np.argsort(-sd)
            sx, sy, kind = sx[order], sy[order], kind[order]
            col = np.array([GRIP, BULK, PARTICLE, colors[k]], dtype=object)[kind]
            siz = np.array([3.4, 2.0, 9.0, 26.0])[kind]
            alp = np.array([0.70, 0.45, 0.95, 1.0])[kind]

            for a in store.values():
                a.remove()
            store.clear()
            sc = ax.scatter(sx, sy, s=siz, c=list(col), linewidths=0,
                            edgecolors="none", zorder=2)
            rgba = sc.get_facecolor()
            rgba[:, 3] = alp
            sc.set_facecolor(rgba)
            store["cloud"] = sc
            d = kind == 3
            if d.any():
                store["def"] = ax.scatter(sx[d], sy[d], s=26, c=colors[k],
                                          linewidths=0.4, edgecolors=PANEL,
                                          zorder=3)
            dpts = shear_frame(ref_pts, curves[k][0][n], L[1], args.slab, hgrip)
            dx_, dy_, _ = project(dpts, M, centre)
            cell_art[keys.index(k)].set_data(dx_, dy_)
            val.set_text(f"τ = {curves[k][1][n]:5.2f} GPa")

        for k in keys:
            g, t = curves[k]
            lines[k].set_data(g[:n + 1], t[:n + 1])
            heads[k].set_data([g[n]], [t[n]])
            htexts[k].set_position((g[n] + gmax * 0.010, t[n]))
            htexts[k].set_text(f"{t[n]:.2f}")

        w = slice(max(0, n - 25), n + 1)
        tc = np.nanmedian(curves["clean"][1][w])
        tp = np.nanmedian(curves["pinned"][1][w])
        if n > 25 and tc > 1e-3:
            caption.set_text(
                f"flow stress just now:  {tc:.2f} GPa with a clear path, "
                f"{tp:.2f} GPa against the particle "
                f"— {tp/tc:.0f}× harder to keep it moving")
        else:
            caption.set_text("both dislocations start from the identical "
                             "relaxed crystal — watch the right-hand one "
                             "reach the particle")

        writer.send(_rgb(fig))
        if n % 20 == 0:
            print(f"  frame {n+1}/{nf}", flush=True)

    writer.close()
    print("wrote", args.movie)


if __name__ == "__main__":
    main()

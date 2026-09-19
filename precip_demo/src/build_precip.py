"""
The edge-dislocation specimen of ../shear_demo, with a hard precipitate placed
on the glide plane directly in the dislocation's path.

Everything else is identical -- same cell frame, same box, same slip system,
same loading -- so the only difference between the two movies is the particle.

The particle is a sphere of atoms held rigid for the whole run: a hard,
non-shearable second-phase particle (an oxide in an ODS steel is the closest
real analogue). The dislocation cannot cut it, so it must bow between the
particle and its periodic images and pinch off an Orowan loop -- the textbook
strengthening mechanism.  Modelling it as rigid, rather than as a chemically
distinct phase, isolates the purely geometric Orowan mechanism; a coherent,
shearable precipitate would instead be cut, and needs an alloy potential.

The particle atoms are given atom type 2. They are still ordinary Fe as far as
the potential is concerned (both types map to Fe), so there is no spurious
misfit strain -- they simply do not move.
"""

import os, sys
import numpy as np

_SHEAR_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "shear_demo", "src")
sys.path.insert(0, os.path.abspath(_SHEAR_SRC))
from build_cells import build_edge, write_data, period_lengths   # noqa: E402

# Measured from the relaxed clean-edge run (shear_demo/out/prod/edge.shear.dump):
# the core sits on y = 31.25 A and glides in -x at about 1.4 A/ps.
SLIP_Y = 31.25
DEFAULT_XC = 12.0       # ~20 A downstream of the core's start at x = 31.9


def build_edge_with_precipitate(nx=26, ny=8, nz=7, radius=12.0,
                                xc=DEFAULT_XC, yc=None, zc=None):
    """Edge specimen plus a rigid spherical particle centred on the glide plane.

    Returns (positions, box, types) with type 2 marking the particle.
    """
    pos, L = build_edge(nx, ny, nz)
    yc = SLIP_Y if yc is None else yc
    zc = 0.5 * L[2] if zc is None else zc

    # nearest-image separation along the periodic x and z
    dx = pos[:, 0] - xc
    dx -= np.round(dx / L[0]) * L[0]
    dz = pos[:, 2] - zc
    dz -= np.round(dz / L[2]) * L[2]
    dy = pos[:, 1] - yc

    inside = (dx ** 2 + dy ** 2 + dz ** 2) < radius ** 2
    types = np.where(inside, 2, 1).astype(int)
    return pos, L, types


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=26)
    ap.add_argument("--ny", type=int, default=8)
    ap.add_argument("--nz", type=int, default=7)
    ap.add_argument("--radius", type=float, default=12.0)
    ap.add_argument("--xc", type=float, default=DEFAULT_XC)
    ap.add_argument("--out", default="out/edge_precip.data")
    a = ap.parse_args()

    pos, L, types = build_edge_with_precipitate(a.nx, a.ny, a.nz, a.radius, a.xc)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    write_data(a.out, pos, L, title="BCC Fe edge dislocation + rigid particle",
               types=types)
    n = int((types == 2).sum())
    vol = 4.0 / 3.0 * np.pi * a.radius ** 3
    print(f"box {np.round(L,3)} A,  {len(pos)} atoms")
    print(f"particle: d = {2*a.radius:.1f} A, {n} atoms "
          f"({100.0*n/len(pos):.1f} % of the cell), centred at "
          f"({a.xc:.1f}, {SLIP_Y:.2f}, {L[2]/2:.1f})")
    print(f"expected from volume: {vol*2/2.855325**3:.0f} atoms")
    print(f"free channel between periodic particles: {L[2]-2*a.radius:.1f} A")
    print("wrote", a.out)

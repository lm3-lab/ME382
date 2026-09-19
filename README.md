# ME 382

Teaching material for ME 382.

## `phase-topology/` — how phase diagrams change shape

An interactive page on the basic question of why a eutectic forms. One slider controls
how strongly A and B atoms resist mixing in the solid; as it increases, a plain freezing
range (a lens) turns into a eutectic. A second preset shows the case where one component
melts far higher than the other, which produces a peritectic instead.

Both panels — the energy curves with their convex hull, and the resulting phase diagram —
are computed live as you drag. Nothing is drawn by hand.

**Open `phase-topology/index.html` in any browser.** No build step, no server, no
dependencies.

| file | what it is |
|---|---|
| `phase-topology/index.html` | the interactive page |
| `phase-topology/README.md` | the model, verified parameters, and implementation notes |
| `phase-topology/verify_topology.py` | dependency-free Python reference for every table in that README |
| `phase-topology/make_web_version.py` | derives the publishable copy from `index.html` |

Run `python3 phase-topology/verify_topology.py` to regenerate the parameter sweeps. The
JavaScript in the page reproduces its output exactly.

## Companion material

A molecular dynamics demonstration of the shear-strength anomaly in BCC Fe is being
developed on a separate branch and is not merged here yet.

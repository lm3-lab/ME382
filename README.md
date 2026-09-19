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

## `dislocations/` — why a crystal is weaker than its bonds

Two molecular dynamics demonstrations in BCC iron, published as movies on the
[course page](https://lm3-lab.github.io/ME382/dislocations/).

The first shears three identical blocks — one perfect, one holding an edge dislocation, one
holding a screw — and shows the perfect lattice reaching the ideal shear strength while the
other two flow one to two orders of magnitude below it. The second puts a hard particle in
the edge dislocation's path and turns free glide into a stick–slip cycle.

> **These are artificial, qualitative demonstrations.** They are built to make a mechanism
> visible, not to predict a material. The crystals are a few nanometres across; the shear is
> applied at about 10⁹ s⁻¹, some twelve orders of magnitude faster than a laboratory test;
> and one dislocation in a cell this size is a density far beyond any real microstructure.
> What happens on screen is real physics — the stresses below are **not** material
> properties of iron and should not be quoted as such. Treat every number as illustrative of
> a trend, not as a measurement.

| folder | what it is |
|---|---|
| `dislocations/index.html` | the published page, with both movies |
| `shear_demo/` | perfect vs edge vs screw: build, run, render, and the validation |
| `precip_demo/` | the same edge dislocation with a hard particle in its way |

Each folder's README covers the construction of the dislocations, the loading, and how to
reproduce every number. `shear_demo/run_all.sh` and `precip_demo/run_all.sh` regenerate
everything from scratch on CPU.

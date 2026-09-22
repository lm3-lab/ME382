# Phase-diagram topology: lens, eutectic, peritectic

An interactive companion to the MD movie. Open `index.html` in any browser — no build step,
no server, no dependencies beyond a webfont that degrades gracefully if offline.

Where the MD movie shows melting kinetically, this shows the equilibrium end of the same
story: which topology a binary phase diagram takes, and why.

## The model

Two phases, a liquid and one solid solution, both regular. The solid is referenced to zero at
each pure end, so the melting terms carry all the temperature dependence.

```
G_l(x) = (1-x) dHA (1 - T/TmA) + x dHB (1 - T/TmB) + Om_l x(1-x) - R T s(x)
G_s(x) =                                             Om_s x(1-x) - R T s(x)
s(x)   = -x ln x - (1-x) ln(1-x)
```

`s(x) >= 0` is the ideal mixing entropy in units of R, so `-R T s(x)` is the `-TS` term of
`F = E - TS`.

`Om_l = 0` throughout, so the liquid stays ideal. **`Om_s` is the only knob.**

Equilibrium at each temperature is the lower convex hull of the two sampled curves, pooled.
A straight hull segment is a tie-line. An *invariant* is a temperature at which one segment
carries three contacts at once.

## Classification

At the invariant, with solid contacts `x_alpha < x_beta` and liquid contact `x_L`:

```
x_alpha < x_L < x_beta   ->  eutectic     L -> alpha + beta
x_L outside that range   ->  peritectic   L + beta -> alpha
```

Same convex-hull event, same code path. Only the contact ordering differs.

## Verified parameter sets

**Branch 1 — lens -> eutectic.** TmA = 700 K, TmB = 900 K, dHA = dHB = 12 kJ/mol.

| Om_s (kJ/mol) | type | T (K) | x_alpha | x_L | x_beta |
|---|---|---|---|---|---|
| <= 10 | lens | — | — | — | — |
| 11 | eutectic | 600.0 | 0.245 | 0.388 | 0.755 |
| 14 | eutectic | 579.6 | 0.080 | 0.388 | 0.920 |

**Branch 2 — lens -> peritectic -> eutectic.** TmA = 700 K, TmB = 1400 K, dHA = 12, dHB = 40 kJ/mol.

| Om_s | type | T (K) | x_alpha | x_L | x_beta |
|---|---|---|---|---|---|
| <= 11 | lens | — | — | — | — |
| 12 | peritectic | 721.4 | 0.485 | 0.020 | 0.515 |
| 13 | peritectic | 746.0 | 0.318 | 0.043 | 0.682 |
| 16 | peritectic | 715.7 | 0.107 | 0.035 | 0.892 |
| 20 | peritectic | 701.2 | 0.040 | 0.033 | 0.960 |
| 21–22 | degenerate | ~698 | ~0.030 | ~0.028 | ~0.970 |
| 23 | eutectic | 695.4 | 0.022 | 0.028 | 0.978 |

The peritectic-to-eutectic crossover is exactly where `x_L` crosses `x_alpha`, as it must be.
At N = 401 the two contacts fall within one grid spacing of each other over Om_s = 21–22, so
that window is reported as degenerate rather than typed.

## The negative result

**A symmetric solid gap with equal melting enthalpies never gives a peritectic.** Holding
dHA = dHB = 12 kJ/mol and TmA = 700 K, sweeping TmB over 900–2400 K against Om_s = 12–36
returns a eutectic in essentially every cell. Two isolated peritectic cells survive at
TmB >= 1800 K right at gap onset (Om_s = 12), and they are gone by Om_s = 14.

So the common claim that *very different melting points give a peritectic* is wrong as stated.
What matters is asymmetry in the melting **entropy**:

```
dS_m = dH_m / T_m      A: 12000/700  = 17.1 J/mol/K
                       B: 40000/1400 = 28.6 J/mol/K
```

A large `dS_m` on B makes the B-side liquidus steep, confining the liquid to the A-rich corner.
The solid miscibility gap, meanwhile, is centred at x = 0.5. When the widening gap reaches the
solidus, both solid contacts therefore lie to the B-rich side of `x_L`, and the reaction is
peritectic. Raising `Om_s` widens the gap until `x_alpha` slides past `x_L` and it converts to
a eutectic.

Also tested and **not** a route to a peritectic: a subregular skew
`Om_s(x) = Om0 + Om1 (1-2x)` with Om1 up to 24 kJ/mol. It shifts `x_L` and `x_alpha` together
and never flips their order.

## Implementation notes

- Composition on a uniform grid, N = 401. Both phases pooled as `(x, G, tag)` through a
  monotone-chain lower hull that carries the tag.
- Reduce to the pointwise lower envelope before hulling. A point sitting above another at the
  same composition can never be a hull vertex, so this leaves the geometry untouched while
  removing the duplicate-x degeneracy that otherwise makes monotone chain retain a non-hull
  point at x = 1.
- A hull edge counts as a tie-line when it spans more than 1.5 grid steps.
- Find the invariant by cooling until the first solid–solid tie-line appears on the hull.
- **Read `x_L` one step above that temperature, never at it.** At the invariant the liquid
  contact is on the point of vanishing and is numerically degenerate — at Branch 1 / Om_s = 14
  there are two L|S tie-lines at 581 K and none at 580 K, where the S|S tie-line appears.
- That same step is the lens discriminator: a solid miscibility gap that opens *below* the
  solidus produces a solid–solid tie-line with no liquid contact above it. That is a lens, not
  an invariant. Testing for liquid at the invariant temperature itself silently reports every
  eutectic as a lens.
- The phase-diagram raster is the per-x hull label stacked over 240 temperature rows. Each x is
  labelled by its bracketing hull edge: joined tags if that edge is a tie-line, otherwise the
  nearer endpoint's tag.

Sanity slice, Branch 2 at Om_s = 16, T = 716 K: liquid lies to the left of *both* solid
contacts. That is the peritectic signature.

## Files

| file | what it is |
|---|---|
| `index.html` | the interactive page — open it directly |
| `verify_topology.py` | pure-Python reference implementation (no dependencies) that produced every table above |

Run `python3 verify_topology.py` to regenerate the sweeps. The JavaScript in `index.html`
reproduces its output exactly; the page's "Recompute in browser" buttons re-derive the tables
live so you can confirm that yourself.

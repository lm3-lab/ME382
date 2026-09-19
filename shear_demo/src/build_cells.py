"""
Build three BCC-Fe specimens that are identical in every way except for the
dislocation they contain.

ONE cell frame is used for all three, because both directions we need are
axes of the same orthogonal frame:

    x = [111]    Burgers-vector direction, b = a/2[111]   (and the shear direction)
    y = [1-10]   normal of the (1-10) slip plane          (loading axis)
    z = [11-2]   the remaining in-plane direction

Slip system is (1-10)[111] throughout, and all three specimens are sheared the
same way -- the top grip is translated along +x -- so the driving resolved
shear stress is tau = sigma_xy in every case.

    perfect   no dislocation                     -> theoretical shear strength
    edge      line along z, glides along x       (line perpendicular to b)
    screw     line along x, glides along z       (line parallel to b)

Because the line directions [11-2] and [111] are both cell axes, the edge and
the screw live in the *same box, same orientation, same size*, and within ~2%
the same number of atoms.  The only difference between the movies is the
character of the dislocation.
"""

import numpy as np

A0 = 2.855325          # Mendelev Fe_mm.eam.fs equilibrium lattice constant (A)
MASS_FE = 55.845

# the single cell frame shared by all three specimens
ORIENT = [[1, 1, 1], [1, -1, 0], [1, 1, -2]]

# ----------------------------------------------------------------------------
# lattice generation
# ----------------------------------------------------------------------------

def _rotation(orient):
    """Rows of the rotation matrix taking cubic coords -> cell coords."""
    R = np.array(orient, dtype=float)
    R /= np.linalg.norm(R, axis=1)[:, None]
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-10), "orient vectors not orthogonal"
    return R


def period_lengths(orient=ORIENT, a=A0):
    """Smallest BCC lattice repeat along each oriented axis, in Angstrom.

    A BCC lattice translation is an integer triple, or an all-half-odd-integer
    triple.  Along [uvw] the repeat is a|[uvw]|/k with k = 2 when a/2[uvw] is
    itself a lattice translation.
    """
    out = []
    for v in np.asarray(orient, dtype=float):
        half = v / 2.0
        allint = np.allclose(half, np.round(half))
        allhalf = np.allclose(np.abs(half - np.round(half)), 0.5)
        k = 2.0 if (allint or allhalf) else 1.0
        out.append(a * np.linalg.norm(v) / k)
    return np.array(out)


BURGERS = period_lengths()[0]          # |a/2[111]| = 2.4728 A


def make_block(n_cells, a=A0, xscale=1.0, orient=ORIENT):
    """BCC atom positions filling n_cells periods of the oriented cell.

    xscale uniformly rescales x (used for the misfit block that carries the
    extra half plane of the edge dislocation).
    Returns (positions Nx3, box lengths 3).
    """
    R = _rotation(orient)
    p = period_lengths(orient, a)
    L = p * np.asarray(n_cells, dtype=float)

    basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) * a
    reach = int(np.ceil(np.linalg.norm(L) / a)) + 2
    rng = np.arange(-reach, reach + 1)
    gx, gy, gz = np.meshgrid(rng, rng, rng, indexing="ij")
    cells = np.stack([gx, gy, gz], axis=-1).reshape(-1, 3) * a
    pts = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 3)

    r = pts @ R.T
    tol = 1e-6
    r = r[np.all((r > -tol) & (r < L[None, :] - tol), axis=1)]

    expected = round(2.0 * np.prod(L) / a ** 3)
    assert len(r) == expected, (
        f"atom count {len(r)} != expected {expected}; box not commensurate")

    r[:, 0] *= xscale
    L = L.copy()
    L[0] *= xscale
    return r, L


# ----------------------------------------------------------------------------
# LAMMPS data file
# ----------------------------------------------------------------------------

def write_data(path, pos, L, mass=MASS_FE, title="BCC Fe", types=None):
    """Write a LAMMPS data file.  `types` is an optional per-atom type array."""
    if types is None:
        types = np.ones(len(pos), dtype=int)
    types = np.asarray(types, dtype=int)
    ntypes = int(types.max())
    with open(path, "w") as f:
        f.write(f"{title}\n\n{len(pos)} atoms\n{ntypes} atom types\n\n")
        f.write(f"0.0 {L[0]:.10f} xlo xhi\n")
        f.write(f"0.0 {L[1]:.10f} ylo yhi\n")
        f.write(f"0.0 {L[2]:.10f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t in range(1, ntypes + 1):
            f.write(f"{t} {mass}\n")
        f.write("\nAtoms # atomic\n\n")
        for i, ((x, y, z), t) in enumerate(zip(pos, types), 1):
            f.write(f"{i} {t} {x:.8f} {y:.8f} {z:.8f}\n")


# ----------------------------------------------------------------------------
# specimens
# ----------------------------------------------------------------------------

def build_perfect(nx, ny, nz, a=A0):
    """Defect-free reference crystal.  ny is the HALF height, in y periods."""
    return make_block((nx, 2 * ny, nz), a)


def build_edge(nx, ny, nz, a=A0):
    """One a/2[111] edge dislocation on (1-10); line along z=[11-2].

    Standard misfit construction: the lower half of the crystal holds nx
    periods of [111] along x, the upper half holds nx-1 periods stretched to
    the same length Lx.  The halves therefore differ by exactly one Burgers
    vector, and on relaxation that misfit localises into a single edge
    dislocation sitting on the joining (1-10) plane.
    """
    px = period_lengths(ORIENT, a)[0]
    lower, Ll = make_block((nx, ny, nz), a)
    upper, Lu = make_block((nx - 1, ny, nz), a,
                           xscale=(nx * px) / ((nx - 1) * px))
    assert abs(Lu[0] - Ll[0]) < 1e-8
    upper[:, 1] += Ll[1]                 # stack on top; the y stacking continues
    return np.vstack([lower, upper]), np.array([Ll[0], Ll[1] + Lu[1], Ll[2]])


def screw_displacement(p, q, pc, qc, b, Lp):
    """Screw displacement along the line, for a periodic row of dislocations.

    `p` is the glide-direction coordinate (periodic, length Lp), `q` the
    slip-plane-normal coordinate.  Complex potential of an infinite array of
    screws of spacing Lp:

        u = (b/2pi) * arg( sin(pi*(p' + i q')/Lp) )

    smooth, and reducing to the isolated (b/2pi)*theta field near the core.
    The bare arg() field leaves a +-b/2 offset between the two sides of the
    p-boundary, which is NOT a lattice translation and would seed a spurious
    fault; adding the uniform shear b*p'/(2*Lp) moves the whole discontinuity
    to one side, where it becomes exactly one Burgers vector -- a genuine
    lattice translation, and therefore invisible to the crystal.
    """
    pp = p - pc
    ap = np.pi * pp / Lp
    aq = np.pi * (q - qc) / Lp
    u = (b / (2 * np.pi)) * np.arctan2(np.cos(ap) * np.sinh(aq),
                                       np.sin(ap) * np.cosh(aq))
    return u + b * pp / (2.0 * Lp)


def easy_core_centre(cols, target):
    """Centre of the triangle of three <111> atomic columns nearest `target`.

    In BCC the a/2<111> screw relaxes into the "easy" core at the centroid of
    a triangle of <111> columns.  Putting the Volterra centre there also keeps
    the branch cut off every atomic row -- an atom sitting exactly on the cut
    would be displaced by b and land on top of its neighbour.
    """
    cols = np.unique(np.round(cols, 4), axis=0)
    c0 = cols[np.argmin(np.hypot(*(cols - target).T))]
    d = np.hypot(*(cols - c0).T)
    order = np.argsort(d)
    nn, dnn = cols[order[1:7]], d[order[1:7]]
    best, bestd = None, np.inf
    for i in range(len(nn)):
        for j in range(i + 1, len(nn)):
            s = np.hypot(*(nn[i] - nn[j])) + dnn[i] + dnn[j]
            if s < bestd:
                bestd, best = s, (nn[i], nn[j])
    return (c0 + best[0] + best[1]) / 3.0


def build_screw(nx, ny, nz, a=A0):
    """One a/2[111] screw dislocation; line along x=[111], glides along z.

    Same box and orientation as the edge specimen -- only the character of the
    dislocation differs.  The displacement is u_x(y, z), so the line is
    trivially periodic along x, and the field is made periodic along the glide
    direction z by the array potential above.
    """
    pos, L = make_block((nx, 2 * ny, nz), a)
    b = BURGERS
    zc, yc = easy_core_centre(pos[:, [2, 1]], np.array([0.5 * L[2], 0.5 * L[1]]))
    pos[:, 0] += screw_displacement(pos[:, 2], pos[:, 1], zc, yc, b, L[2])
    return pos, L

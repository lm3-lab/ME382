"""Fast reader for LAMMPS text dump files."""
import numpy as np


def read_frames(fn, fields=None, max_frames=None):
    """Read all frames of a LAMMPS dump.

    Returns a list of dicts with 'step', 'box' and one float32 array per
    requested column (default: every column in the file).
    """
    with open(fn) as f:
        txt = f.read()
    chunks = txt.split("ITEM: TIMESTEP\n")[1:]
    out = []
    for ch in chunks:
        head, _, body = ch.partition("ITEM: ATOMS ")
        if not body:
            continue                       # truncated trailing frame
        hl = head.split("\n")
        step = int(hl[0])
        n = int(hl[2])
        box = np.array([[float(v) for v in hl[4 + k].split()[:2]] for k in range(3)])
        cols = body.split("\n", 1)[0].split()
        rows = body.split("\n", 1)[1]
        vals = np.array(rows.split(), dtype=np.float32)
        if vals.size < n * len(cols):
            continue                       # frame still being written
        vals = vals[: n * len(cols)].reshape(n, len(cols))
        fr = dict(step=step, box=box)
        for k, c in enumerate(cols):
            if fields is None or c in fields:
                fr[c] = vals[:, k]
        out.append(fr)
        if max_frames and len(out) >= max_frames:
            break
    return out

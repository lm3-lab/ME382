"""
Strain-controlled simple shear of BCC Fe on the (1-10)[111] slip system.

Loading: the crystal is sandwiched between two rigid slabs normal to y.
The lower slab is held; the upper slab is translated at constant velocity
along the Burgers-vector direction, so the shear strain

    gamma(t) = v * t / h          (h = slab centre-to-centre separation)

increases linearly -- this is displacement (strain) controlled loading.
The resolved shear stress on the slip plane is measured two ways:

    tau_wall   = (total force on the upper rigid slab) / (Lx*Lz)
    tau_virial = virial stress of the mobile region

Both are reported; tau_virial is the smoother of the two and is what the
movie plots.
"""

import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_cells import (A0, ORIENT, BURGERS, build_perfect, build_edge,
                         build_screw, write_data, period_lengths)

POT = "/usr/local/lib/python3.11/dist-packages/lammps/share/lammps/potentials/Fe_mm.eam.fs"
BAR_TO_GPA = 1.0e-4
CNA_CUT = 0.5 * (1.0 + np.sqrt(2.0)) * A0      # 3.447 A, the standard BCC CNA cutoff


BUILDERS = {"perfect": build_perfect, "edge": build_edge, "screw": build_screw}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=["perfect", "edge", "screw"])
    ap.add_argument("--nx", type=int, default=26)
    ap.add_argument("--ny", type=int, default=8)      # half-height, in periods
    ap.add_argument("--nz", type=int, default=7)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--rate", type=float, default=1.0e9)   # 1/s
    ap.add_argument("--gamma-max", type=float, default=0.16)
    ap.add_argument("--dt", type=float, default=0.002)     # ps
    ap.add_argument("--slab", type=float, default=8.0)     # rigid slab thickness, A
    ap.add_argument("--frames", type=int, default=240)
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--relax-only", action="store_true")
    ap.add_argument("--data", default=None,
                    help="read this prebuilt data file instead of building one")
    ap.add_argument("--freeze-type", type=int, default=0,
                    help="hold every atom of this type rigid (a hard particle)")
    ap.add_argument("--tag", default=None, help="output basename")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or args.case
    if args.data:
        datafile = args.data
        L = np.array([float(l.split()[1]) for l in open(datafile)
                      if l.rstrip().endswith(("xhi", "yhi", "zhi"))])
        npos = int(next(l for l in open(datafile) if l.rstrip().endswith("atoms")).split()[0])
        print(f"[{tag}] {npos} atoms from {datafile}, box = {np.round(L,3)} A",
              flush=True)
    else:
        datafile = os.path.join(args.outdir, f"{tag}.data")
        pos, L = BUILDERS[args.case](args.nx, args.ny, args.nz)
        write_data(datafile, pos, L, title=f"BCC Fe {args.case} (1-10)[111]")
        print(f"[{tag}] {len(pos)} atoms, box = {np.round(L,3)} A", flush=True)

    from lammps import lammps
    lmp = lammps(cmdargs=["-log", os.path.join(args.outdir, f"{tag}.log"),
                          "-screen", "none", "-nocite"])
    c = lmp.commands_string

    ylo_hi = args.slab
    yhi_lo = L[1] - args.slab
    h = L[1] - args.slab                      # slab centre-to-centre separation
    area = L[0] * L[2]                        # slip-plane area
    vol_mobile = area * (yhi_lo - ylo_hi)

    if args.freeze_type:
        # The particle is NOT frozen during minimisation -- it relaxes with the
        # rest of the crystal, so the starting structure is identical to the
        # dislocation-only run.  It is excluded from `mobile`, so from the
        # thermal equilibration onwards it is never integrated, i.e. rigid.
        freeze_cmds = (f"group precip type {args.freeze_type}\n"
                       f"group mobile subtract interior precip")
        pair_map = " ".join(["Fe"] * args.freeze_type)
    else:
        freeze_cmds = "group mobile union interior"
        pair_map = "Fe"

    c(f"""
units metal
dimension 3
boundary p s p
atom_style atomic
atom_modify map array sort 0 0.0
read_data {datafile}
pair_style eam/fs
pair_coeff * * {POT} {pair_map}

region rlo block INF INF INF {ylo_hi} INF INF units box
region rhi block INF INF {yhi_lo} INF INF INF units box
group lo region rlo
group hi region rhi
group boundary union lo hi
group interior subtract all boundary
{freeze_cmds}

compute cna all cna/atom {CNA_CUT}
compute cs  all centro/atom bcc
compute peat all pe/atom
compute sa interior stress/atom NULL
compute svir interior reduce sum c_sa[1] c_sa[2] c_sa[3] c_sa[4] c_sa[5] c_sa[6]

timestep {args.dt}
velocity all set 0.0 0.0 0.0
fix frz boundary setforce 0.0 0.0 0.0

min_style cg
min_modify dmax 0.05
minimize 1e-12 1e-12 20000 40000
""")

    natoms = lmp.get_natoms()
    print(f"[{tag}] relaxed, PE = {lmp.get_thermo('pe'):.4f} eV "
          f"({lmp.get_thermo('pe')/natoms:.6f} eV/atom)", flush=True)

    # report the defect content of the relaxed cell
    c("run 0")
    cna = np.array(lmp.numpy.extract_compute("cna", 1, 1))   # per-atom, vector
    ndef = int(np.sum(cna != 3))                             # 3 = BCC
    print(f"[{tag}] non-BCC (CNA) atoms after relaxation: {ndef} / {natoms}", flush=True)

    dumpfile = os.path.join(args.outdir, f"{tag}.relaxed.dump")
    c(f"""
dump d0 all custom 1 {dumpfile} id type x y z c_cna c_cs c_peat
dump_modify d0 sort id
run 0
undump d0
""")
    if args.relax_only:
        lmp.close()
        return

    # ------------------------------------------------------------------
    # shear
    # ------------------------------------------------------------------
    # v [A/ps] = rate [1/s] * h [A] * 1e-12 [s/ps]
    v = args.rate * h * 1.0e-12
    total_time = args.gamma_max / (args.rate * 1.0e-12)      # ps
    nsteps = int(round(total_time / args.dt))
    # fix ave/time requires nfreq to be a multiple of nevery, so round the
    # output interval to a multiple of NEVERY and the run length to a whole
    # number of intervals.
    NEVERY = 10
    every = max(NEVERY, (nsteps // args.frames) // NEVERY * NEVERY)
    nsteps = (nsteps // every) * every
    # All three specimens share one cell frame and one loading direction:
    # the top grip slides along +x = [111] = b, so the resolved shear stress
    # on the (1-10) slip plane is tau = sigma_xy for every case.
    vx, vy, vz = v, 0.0, 0.0
    scomp, fcomp = 4, 1                  # c_sa[4] = xy ; force component x
    tauvar = f"(c_svir[{scomp}]/{vol_mobile})*{BAR_TO_GPA}"

    freeze_shear = ("fix frz_p precip setforce 0.0 0.0 0.0"
                    if args.freeze_type else "")
    thermofile = os.path.join(args.outdir, f"{tag}.stress.txt")
    trajfile = os.path.join(args.outdir, f"{tag}.shear.dump")

    c(f"""
velocity mobile create {args.temp} {args.seed} mom yes rot no dist gaussian
fix integ mobile nve
fix thermo_l mobile langevin {args.temp} {args.temp} 1.0 {args.seed + 7} zero yes
run 2000

reset_timestep 0
unfix frz
fix frz_lo lo setforce 0.0 0.0 0.0
velocity lo set 0.0 0.0 0.0
fix mv_hi hi move linear {vx} {vy} {vz} units box
fix frc_hi hi setforce 0.0 0.0 0.0
{freeze_shear}

variable gam  equal (step*{args.dt}*{v})/{h}
variable tauv equal {tauvar}
variable tauw equal (f_frc_hi[{fcomp}]/{area})*160.21766208
variable tmp  equal temp

fix out all ave/time {NEVERY} {every//NEVERY} {every} &
    v_gam v_tauv v_tauw v_tmp &
    file {thermofile} &
    title2 "# step gamma tau_virial_GPa tau_wall_GPa T_K (interval averages)"

dump dt all custom {every} {trajfile} id type x y z c_cna
dump_modify dt sort id format line "%d %d %.3f %.3f %.3f %.0f" 

thermo {every*10}
thermo_style custom step v_gam v_tauv v_tauw temp pe
run {nsteps}
""")
    print(f"[{tag}] done: {nsteps} steps, gamma_max={args.gamma_max}, "
          f"v={v:.5f} A/ps, frames every {every} steps", flush=True)
    lmp.close()


if __name__ == "__main__":
    main()

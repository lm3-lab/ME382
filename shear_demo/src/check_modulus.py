"""Reference shear modulus of the (1-10)[111] system for this potential.

Fully periodic perfect crystal, small xy tilt applied by hand, energy
minimised at fixed strain.  This is the G the elastic branch of the
stress-strain curves should reproduce.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_cells import make_block, write_data
from lammps import lammps

POT = "/usr/local/lib/python3.11/dist-packages/lammps/share/lammps/potentials/Fe_mm.eam.fs"

pos, L = make_block((8, 6, 4))
write_data("/tmp/gref.data", pos, L)
gammas, taus = [], []
for g in [0.0, 0.002, 0.004, 0.006, 0.008]:
    lmp = lammps(cmdargs=["-log", "none", "-screen", "none", "-nocite"])
    lmp.commands_string(f"""
units metal
boundary p p p
atom_style atomic
read_data /tmp/gref.data
pair_style eam/fs
pair_coeff * * {POT} Fe
change_box all triclinic
change_box all xy delta {g*L[1]} remap units box
min_style cg
minimize 1e-14 1e-14 10000 20000
variable tau equal -pxy*1.0e-4
run 0
""")
    taus.append(lmp.extract_variable("tau"))
    gammas.append(g)
    lmp.close()
g = np.array(gammas); t = np.array(taus)
G = np.polyfit(g, t, 1)[0]
print("gamma :", np.round(g, 4))
print("tau   :", np.round(t, 4), "GPa")
print(f"\nG(1-10)[111] = {G:.1f} GPa   (fully periodic, athermal reference)")

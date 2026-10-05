#!/bin/bash -eE
# Union-of-spheres centers = all solvent-exposed protein heavy atoms (SASA > cutoff) of a structure.
# usage: make_centers.sh <structure.gro> [sasa_cut_nm2 (default 0)]
cd "$(dirname "$0")"   # needs gmx (any version) and MDAnalysis on the PATH
gro=$(realpath "$1"); cut=${2:-0.0}
gmx sasa -f "$gro" -s "$gro" -o sasa_total.xvg -oa sasa_atom.xvg \
    -surface 'group "Protein-H"' -output 'group "Protein-H"' -probe 0.14 > sasa.log 2>&1
python3 - "$gro" "$cut" <<'PY'
import sys, numpy as np, MDAnalysis as mda
gro, cut = sys.argv[1], float(sys.argv[2])
u = mda.Universe(gro)
heavy = u.select_atoms('protein and not name H*')
area = np.array([float(l.split()[1]) for l in open('sasa_atom.xvg') if l[0] not in '#@'])
assert len(area) == heavy.n_atoms, (len(area), heavy.n_atoms)
exposed = heavy[area > cut]
np.savetxt('centers.pos', exposed.positions / 10.0, fmt='%.5f',
           header=f'union-of-spheres centers (nm): protein heavy atoms with SASA > {cut} nm^2 in {gro}')
np.savetxt('center_serials.txt', exposed.ids, fmt='%d', header='1-based serials of center atoms')
sol = u.select_atoms('resname SOL and name OW')
print(f"protein heavy atoms: {heavy.n_atoms}; exposed (SASA > {cut}): {exposed.n_atoms}; buried: {heavy.n_atoms - exposed.n_atoms}")
print(f"water oxygens: {sol.n_atoms}, first serial {sol.ids[0]}, stride {sol.ids[1]-sol.ids[0]}, last {sol.ids[-1]}")
open('ow_range.txt','w').write(f"{sol.ids[0]}-{sol.ids[-1]}:{sol.ids[1]-sol.ids[0]}\n")
PY

# Dewetting hydrophobin HFBII with a union of spheres

Inputs used for the example in the top-level README. The simulation system (HFBII chain A
of PDB 2B97, AMBER99SB-ILDN, SPC/E, 5.8 nm box, protein heavy atoms position-restrained) is
prepared with standard GROMACS tools; `test/sample_traj/hfbii/` holds the equilibrated
structure and the centers file.

| File | Purpose |
| --- | --- |
| `select_exposed_atoms.sh` | `gmx sasa` + MDAnalysis: writes the positions of all protein heavy atoms with nonzero SASA (the sphere centers) and the water-oxygen index range |
| `indus.input` | INDUS input (PLUMED mode): water oxygens as targets, `union_of_spheres` with r = 0.6 nm on those centers |
| `plumed.dat` | `INDUS` action plus a `MOVINGRESTRAINT` that ramps N* from 542 to 0 over 2 ns and holds it for 1 ns (kappa = 0.5 kJ/mol) |
| `biased.mdp` | NPT production settings used for the biased run (2 fs, C-rescale, position restraints on) |
| `plot_ntilde.py` | plots Ntilde(t) against the moving target from `plumed.out` |

Run, with a PLUMED-patched GROMACS whose kernel contains INDUS (see
`scripts/install/install_indus_stack.sh`):

```
gmx_mpi grompp -f biased.mdp -c equilibrated.gro -r equilibrated.gro -t equilibrated.cpt -p topol.top -o biased.tpr
mpirun -np 1 gmx_mpi mdrun -s biased.tpr -deffnm biased -plumed plumed.dat -ntomp 12
python3 plot_ntilde.py plumed.out ntilde_vs_time.png
```

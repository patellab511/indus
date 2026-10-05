#!/usr/bin/env python3
"""Check a PLUMED-driver probe-volume test against its standalone twin in test/indus.

usage: compare_with_standalone.py <plumed test dir> <standalone test dir> <first target serial> <stride>
Compares Ntilde per frame (plumed.out vs ref/time_samples_indus.out) and the bias force on
every target atom (forces.out vs ref/forces_debug.out). Both runs use kappa = 1, N* = 0.
"""
import sys
from pathlib import Path
import numpy as np

pdir, sdir = Path(sys.argv[1]), Path(sys.argv[2])
first, stride = int(sys.argv[3]), int(sys.argv[4])

nt_plumed = np.array([float(l.split()[2]) for l in (pdir / "plumed.out").read_text().splitlines() if l and not l.startswith("#")])
nt_ref = np.array([float(l.split()[2]) for l in (sdir / "ref/time_samples_indus.out").read_text().splitlines() if l and not l.startswith("#")])

# plumed --dump-forces: per frame "natoms", a box/virial line, then one "X fx fy fz" line per atom in order
lines = (pdir / "forces.out").read_text().splitlines()
frames, i = [], 0
while i < len(lines):
    n = int(lines[i].split()[0]); block = lines[i + 2:i + 2 + n]
    frames.append(np.array([[float(x) for x in l.split()[1:4]] for l in block])); i += 2 + n
serials = np.arange(first, frames[0].shape[0] + 1, stride)

ref_frames, cur = [], None
for l in (sdir / "ref/forces_debug.out").read_text().splitlines():
    if l.startswith("t="): cur = {}; ref_frames.append(cur)
    elif cur is not None and l[:1].isdigit():
        t = l.split(); cur[int(t[0])] = np.array(t[1:4], float)

assert len(nt_plumed) == len(nt_ref) == len(frames) == len(ref_frames), (len(nt_plumed), len(nt_ref), len(frames), len(ref_frames))
d_nt = np.abs(nt_plumed - nt_ref).max()
d_f = 0.0
for fp, fr in zip(frames, ref_frames):
    f_ref = np.array([fr.get(int(s), np.zeros(3)) for s in serials])
    d_f = max(d_f, float(np.abs(fp[serials - 1] - f_ref).max()))
    others = np.delete(fp, serials - 1, axis=0)
    assert np.abs(others).max() == 0.0, "force on a non-target atom"
print(f"{pdir.name}: frames {len(frames)}, max |dNtilde| = {d_nt:.2e}, max |dF| = {d_f:.2e} kJ/mol/nm")
ok = d_nt <= 1.5e-3 and d_f <= 2.0e-4
print("  " + ("PASS" if ok else "FAIL")); sys.exit(0 if ok else 1)

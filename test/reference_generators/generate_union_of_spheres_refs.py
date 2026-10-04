#!/usr/bin/env python3
"""Independent PyTorch reference for the static union_of_spheres probe volume.

Reads a test directory's indus.input (files, target range, radii, centers file,
bias) and recomputes, for every trajectory frame,

    htilde_k = INDUS switching function of |x_i - c_k| for each sphere center c_k
    htilde_v = 1 - prod_k (1 - htilde_k)               (union of spheres)
    Ntilde   = sum_i htilde_v(i),   N = # targets with h_k = 1 for some k
    U_bias   = 0.5 * kappa * (Ntilde - x_star)^2

with the bias forces on every target atom from autograd, plus the virial in the
driver's convention. It then compares against the driver's time_samples_indus.out
and forces_debug.out in the test directory.

This script is outside the default CTest path because it needs MDAnalysis and torch.

usage:
  generate_union_of_spheres_refs.py [--test-dir DIR] [--indus-exe EXE] [--update-ref]
    --indus-exe   run the driver in the test directory first
    --update-ref  copy the validated driver outputs into DIR/ref (only if they pass)
"""
from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import torch

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TEST_DIR = REPO_ROOT / "test/indus/union_of_spheres/hfbii"
DTYPE = torch.float64
REF_FILES = ["time_samples_indus.out", "forces_debug.out"]


# ----- INDUS switching function (same form as SwitchingFunction_Indus_x) -----
def switching_constants(sigma: float, alpha_c: float):
    two_sigma_sq = 2.0 * sigma * sigma
    sqrt_two_sigma = math.sqrt(2.0) * sigma
    k = math.sqrt(math.pi * two_sigma_sq) * math.erf(alpha_c / sqrt_two_sigma)
    k -= 2.0 * alpha_c * math.exp(-(alpha_c * alpha_c) / two_sigma_sq)
    return sqrt_two_sigma, math.sqrt(0.5 * math.pi * sigma * sigma) / k, math.exp(-(alpha_c * alpha_c) / two_sigma_sq) / k


def switch_x(x, x_min, x_max, sigma, alpha_c):
    sqrt_two_sigma, k1, k2 = switching_constants(sigma, alpha_c)
    h = ((x >= x_min) & (x <= x_max)).to(x.dtype)
    dx_min, dx_max = x - x_min, x - x_max
    htilde = torch.ones_like(x)
    htilde = torch.where(torch.abs(dx_min) <= alpha_c, 0.5 + k1 * torch.erf(dx_min / sqrt_two_sigma) - k2 * dx_min, htilde)
    htilde = torch.where(torch.abs(dx_max) <= alpha_c, 0.5 + k1 * torch.erf(-dx_max / sqrt_two_sigma) + k2 * dx_max, htilde)
    htilde = torch.where((x < x_min - alpha_c) | (x > x_max + alpha_c), torch.zeros_like(x), htilde)
    return h, htilde


def minimum_image(dx, box):
    return dx - box * torch.round(dx / box)


# ----- indus.input parsing (just the keys this test uses) -----
def parse_indus_input(path: Path) -> dict:
    text = path.read_text()
    def key(name, default=None, cast=str):
        m = re.search(rf"^\s*{name}\s*=\s*([^\s#]+)", text, re.M)
        return cast(m.group(1)) if m else default
    m = re.search(r"Target\s*=\s*\[\s*atom_index\s+(\d+)-(\d+):(\d+)\s*\]", text)
    if not m:
        raise SystemExit("expected 'Target = [ atom_index first-last:stride ]' in indus.input")
    return {
        "gro": path.parent / key("GroFile"),
        "xtc": path.parent / key("XtcFile"),
        "first": int(m.group(1)), "last": int(m.group(2)), "stride": int(m.group(3)),
        "centers": path.parent / key("reference_positions_file"),
        "r_max": key("r_max", cast=float), "r_min": key("r_min", -1.0, float),
        "sigma": key("sigma", 0.01, float), "alpha_c": key("alpha_c", 0.02, float),
        "kappa": key("kappa", 1.0, float), "x_star": key("x_star", 0.0, float),
    }


# ----- per-frame union calculation -----
def compute_frame(x_nm: np.ndarray, box_nm: np.ndarray, centers: torch.Tensor, p: dict):
    box = torch.tensor(box_nm, dtype=DTYPE)
    x = torch.tensor(x_nm, dtype=DTYPE, requires_grad=True)
    dx = minimum_image(x[:, None, :] - centers[None, :, :], box)
    r = torch.sqrt((dx * dx).sum(-1) + 1e-300)
    h_k, htilde_k = switch_x(r, p["r_min"], p["r_max"], p["sigma"], p["alpha_c"])
    htilde_v = 1.0 - torch.prod(1.0 - htilde_k, dim=1)
    n = int((h_k.max(dim=1).values == 1.0).sum())
    ntilde = htilde_v.sum()
    (0.5 * p["kappa"] * (ntilde - p["x_star"]) ** 2).backward()
    forces = -x.grad.detach()
    virial = 0.5 * torch.einsum("ia,ib->ab", torch.remainder(x.detach(), box), -forces)
    return n, float(ntilde), forces.numpy(), virial.numpy()


# ----- driver output parsing -----
def parse_time_samples(path: Path):
    rows = {}
    for line in path.read_text().splitlines():
        if line and not line.startswith("#"):
            t, n, nt = line.split()[:3]
            rows[round(float(t), 3)] = (int(n), float(nt))
    return rows


def parse_forces(path: Path):
    frames, cur = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(("#", "N_total")):
            continue
        if line.startswith("t="):
            cur = {"forces": {}, "virial": None}
            frames[round(float(line.split()[1]), 3)] = cur
        elif line.startswith("virial="):
            cur["virial"] = np.array(line.split()[1:], dtype=float).reshape(3, 3)
        elif not line.startswith("box="):
            tok = line.split()
            cur["forces"][int(tok[0])] = np.array(tok[1:4], dtype=float)
    return frames


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-dir", type=Path, default=DEFAULT_TEST_DIR)
    ap.add_argument("--indus-exe", type=Path, help="run this driver in the test dir first")
    ap.add_argument("--update-ref", action="store_true", help="copy validated outputs into ref/")
    ap.add_argument("--ntilde-tol", type=float, default=1.0e-3)
    ap.add_argument("--force-tol", type=float, default=2.0e-4)
    ap.add_argument("--virial-tol", type=float, default=2.0e-4)
    args = ap.parse_args()
    test_dir = args.test_dir.resolve()

    if args.indus_exe:
        subprocess.run([str(args.indus_exe.resolve()), "indus.input"], cwd=test_dir, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    p = parse_indus_input(test_dir / "indus.input")
    centers = torch.tensor(np.loadtxt(p["centers"], comments="#", ndmin=2), dtype=DTYPE)
    u = mda.Universe(str(p["gro"]), str(p["xtc"]), convert_units=False)  # nm
    idx0 = np.arange(p["first"] - 1, p["last"], p["stride"])
    serials = idx0 + 1
    print(f"test dir: {test_dir}\ncenters: {len(centers)}  r_max={p['r_max']} r_min={p['r_min']}  targets: {len(idx0)}")

    ts_ref = parse_time_samples(test_dir / "time_samples_indus.out")
    f_ref = parse_forces(test_dir / "forces_debug.out")

    worst = dict(ntilde=0.0, force=0.0, virial=0.0)
    n_mismatch = missing = compared = 0
    print("\n   t(ps)   N_torch  Ntilde_torch   N_cpp  Ntilde_cpp    max|dF|  max|F_cpp|")
    for ts in u.trajectory:
        t = round(float(ts.time), 3)
        n, ntilde, forces, virial = compute_frame(ts.positions[idx0].astype(float), ts.dimensions[:3].astype(float), centers, p)
        if t not in ts_ref or t not in f_ref:
            missing += 1
            print(f"{t:8.1f}  {n:7d}  {ntilde:12.3f}   (no driver output for this frame)")
            continue
        n_ref, nt_ref = ts_ref[t]
        fr = f_ref[t]
        f_cpp = np.array([fr["forces"].get(int(s), np.zeros(3)) for s in serials])
        d_f = float(np.abs(forces - f_cpp).max())
        worst["ntilde"] = max(worst["ntilde"], abs(ntilde - nt_ref))
        worst["force"] = max(worst["force"], d_f)
        worst["virial"] = max(worst["virial"], float(np.abs(virial - fr["virial"]).max()))
        n_mismatch += int(n != n_ref)
        compared += 1
        print(f"{t:8.1f}  {n:7d}  {ntilde:12.3f}  {n_ref:6d}  {nt_ref:10.3f}  {d_f:9.2e}  {np.abs(f_cpp).max():10.3f}")

    ok = (compared > 0 and missing == 0 and n_mismatch == 0 and worst["ntilde"] <= args.ntilde_tol
          and worst["force"] <= args.force_tol and worst["virial"] <= args.virial_tol)
    print(f"\nframes compared: {compared}  missing: {missing}  N mismatches: {n_mismatch}")
    print(f"max |dNtilde| = {worst['ntilde']:.3e} (tol {args.ntilde_tol})")
    print(f"max |dForce|  = {worst['force']:.3e} (tol {args.force_tol}) kJ/mol/nm")
    print(f"max |dVirial| = {worst['virial']:.3e} (tol {args.virial_tol})")
    print("RESULT:", "PASS" if ok else "FAIL")

    if args.update_ref:
        if not ok:
            sys.exit("refusing to update ref/: validation failed")
        ref = test_dir / "ref"
        ref.mkdir(exist_ok=True)
        for name in REF_FILES + ["indus.input", "output_files_to_check.input"]:
            shutil.copy2(test_dir / name, ref / name)
        print(f"updated {ref}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

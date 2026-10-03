"""Shared helpers for INDUS probe-volume reference generators."""

from __future__ import annotations

import math
from pathlib import Path

import MDAnalysis as mda
import torch


REPO_ROOT = Path(__file__).resolve().parents[2]
GRO = REPO_ROOT / "test/sample_traj/bulk_water_278K/initconf.gro"
XTC = REPO_ROOT / "test/sample_traj/bulk_water_278K/traj_water_278K_b1000_e1010.xtc"

SIGMA = 0.01
ALPHA_C = 0.02
X_STAR = 0.0
KAPPA = 1.0

TARGET_INDICES = list(range(0, 16500, 4))
TARGET_SERIALS = [index + 1 for index in TARGET_INDICES]

DTYPE = torch.float64


def put_in_box(x: torch.Tensor, box_lengths: torch.Tensor) -> torch.Tensor:
    return torch.remainder(x, box_lengths)


def minimum_image_vector(
    x1: torch.Tensor, x2: torch.Tensor, box_lengths: torch.Tensor
) -> torch.Tensor:
    dx = x2 - x1
    half = 0.5 * box_lengths
    dx = torch.where(dx > half, dx - box_lengths, dx)
    dx = torch.where(dx < -half, dx + box_lengths, dx)
    return dx


def shifted_to_reference(
    x: torch.Tensor, ref: torch.Tensor, box_lengths: torch.Tensor
) -> torch.Tensor:
    dx = x - ref
    half = 0.5 * box_lengths
    shift = torch.zeros_like(dx)
    shift = torch.where(dx > half, -box_lengths, shift)
    shift = torch.where(dx < -half, box_lengths, shift)
    return x + shift


def switching_constants() -> tuple[float, float, float, float]:
    two_sigma_sq = 2.0 * SIGMA * SIGMA
    sqrt_two_sigma = math.sqrt(2.0) * SIGMA
    k = math.sqrt(math.pi * two_sigma_sq) * math.erf(ALPHA_C / sqrt_two_sigma)
    k -= 2.0 * ALPHA_C * math.exp(-(ALPHA_C * ALPHA_C) / two_sigma_sq)
    k1 = math.sqrt(0.5 * math.pi * SIGMA * SIGMA) / k
    k2 = math.exp(-(ALPHA_C * ALPHA_C) / two_sigma_sq) / k
    return two_sigma_sq, sqrt_two_sigma, k1, k2


TWO_SIGMA_SQ, SQRT_TWO_SIGMA, K1, K2 = switching_constants()


def switch_r(r: torch.Tensor, r_max: float) -> tuple[torch.Tensor, torch.Tensor]:
    r_inner = r_max - ALPHA_C
    r_outer = r_max + ALPHA_C
    h = torch.where(r <= r_max, torch.ones_like(r), torch.zeros_like(r))
    delta = r_max - r
    htilde_buffer = 0.5 + K1 * torch.erf(delta / SQRT_TWO_SIGMA) - K2 * delta
    htilde = torch.where(r > r_outer, torch.zeros_like(r), htilde_buffer)
    htilde = torch.where(r <= r_inner, torch.ones_like(r), htilde)
    h = torch.where(r > r_outer, torch.zeros_like(r), h)
    h = torch.where(r <= r_inner, torch.ones_like(r), h)
    return h, htilde


def switch_x(
    x: torch.Tensor, x_min: float | torch.Tensor, x_max: float | torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    x_lowest = x_min - ALPHA_C
    x_highest = x_max + ALPHA_C
    h = torch.where((x >= x_min) & (x <= x_max), torch.ones_like(x), torch.zeros_like(x))
    dx_min = x - x_min
    dx_max = x - x_max
    left = torch.abs(dx_min) <= ALPHA_C
    right = torch.abs(dx_max) <= ALPHA_C
    htilde_left = 0.5 + K1 * torch.erf(dx_min / SQRT_TWO_SIGMA) - K2 * dx_min
    htilde_right = 0.5 + K1 * torch.erf(-dx_max / SQRT_TWO_SIGMA) + K2 * dx_max
    htilde = torch.ones_like(x)
    htilde = torch.where(left, htilde_left, htilde)
    htilde = torch.where(right, htilde_right, htilde)
    htilde = torch.where((x < x_lowest) | (x > x_highest), torch.zeros_like(x), htilde)
    return h, htilde


def build_frame(
    time: float,
    box_nm,
    positions: torch.Tensor,
    h: torch.Tensor,
    htilde: torch.Tensor,
) -> dict:
    ntilde = htilde.sum()
    bias = 0.5 * KAPPA * (ntilde - X_STAR) * (ntilde - X_STAR)
    bias.backward()

    target_forces = -positions.grad.detach()
    d_u_target = -target_forces
    target_positions_box = put_in_box(positions.detach(), torch.tensor(box_nm, dtype=DTYPE))
    virial = 0.5 * torch.einsum("ia,ib->ab", target_positions_box, d_u_target)
    return {
        "time": time,
        "box": [float(x) for x in box_nm],
        "n": int(torch.count_nonzero(h == 1.0).item()),
        "ntilde": float(ntilde.detach()),
        "forces": target_forces.numpy(),
        "virial": virial.numpy(),
    }


def expected_frames(compute_frame):
    offset_cache_pattern = f".{XTC.name}_offsets.*"
    existing_offset_caches = set(XTC.parent.glob(offset_cache_pattern))
    try:
        universe = mda.Universe(str(GRO), str(XTC), convert_units=False)
        frames = []
        for ts in universe.trajectory:
            frames.append(
                compute_frame(
                    float(ts.time),
                    universe.atoms.positions,
                    ts.dimensions[:3],
                )
            )
        return frames
    finally:
        for cache_path in XTC.parent.glob(offset_cache_pattern):
            if cache_path not in existing_offset_caches:
                cache_path.unlink(missing_ok=True)


def parse_time_samples(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) == 3:
            rows.append({"time": float(fields[0]), "n": int(fields[1]), "ntilde": float(fields[2])})
    return rows


def parse_forces(path: Path) -> list[dict]:
    frames = []
    current = None
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("N_total="):
            continue
        fields = line.split()
        if fields[0] == "t=":
            if current is not None:
                frames.append(current)
            current = {"time": float(fields[1]), "virial": None, "forces": []}
        elif fields[0] == "virial=":
            current["virial"] = [float(x) for x in fields[1:]]
        elif fields[0] == "box=":
            continue
        else:
            current["forces"].append((int(fields[0]), [float(x) for x in fields[1:4]]))
    if current is not None:
        frames.append(current)
    return frames


def write_generated_files(output_dir: Path, frames: list[dict], serials: list[int] = TARGET_SERIALS):
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "time_samples_indus.out").open("w") as handle:
        handle.write("# Generated INDUS reference values\n")
        handle.write("# Time[ps]\tN_v\tNtilde_v\n")
        for frame in frames:
            handle.write(f"{frame['time']:g}\t{frame['n']}\t{frame['ntilde']:.6g}\n")

    with (output_dir / "forces_debug.out").open("w") as handle:
        handle.write("# Generated biasing forces on OrderParameters atoms\n")
        handle.write(f"N_total= 16500  N_biased= {len(serials)}\n")
        for frame in frames:
            box = frame["box"]
            handle.write(f"t= {frame['time']:g}\n")
            handle.write(
                "box= "
                + "".join(f" {value:.5f}" for value in [box[0], 0.0, 0.0, 0.0, box[1], 0.0, 0.0, 0.0, box[2]])
                + "\n"
            )
            handle.write(
                "virial= "
                + "".join(f" {value:.5f}" for value in frame["virial"].reshape(9))
                + "\n"
            )
            for serial, force in zip(serials, frame["forces"]):
                handle.write(f"{serial}  {force[0]:.5f}  {force[1]:.5f}  {force[2]:.5f}\n")


def compare_generated_to_ref(
    generated_dir: Path,
    test_dir: Path,
    ntilde_tol: float = 1.0e-3,
    force_tol: float = 2.0e-4,
    virial_tol: float = 2.0e-4,
) -> dict:
    ref_dir = test_dir / "ref"
    generated_time = parse_time_samples(generated_dir / "time_samples_indus.out")
    ref_time = parse_time_samples(ref_dir / "time_samples_indus.out")
    generated_forces = parse_forces(generated_dir / "forces_debug.out")
    ref_forces = parse_forces(ref_dir / "forces_debug.out")

    result = {
        "frames": len(ref_time),
        "n_mismatches": 0,
        "serial_mismatches": 0,
        "max_ntilde_diff": 0.0,
        "max_force_diff": 0.0,
        "max_virial_diff": 0.0,
    }

    if len(generated_time) != len(ref_time) or len(generated_forces) != len(ref_forces):
        result["length_mismatch"] = True
        result["passed"] = False
        return result

    result["length_mismatch"] = False
    for generated, ref in zip(generated_time, ref_time):
        if generated["n"] != ref["n"]:
            result["n_mismatches"] += 1
        result["max_ntilde_diff"] = max(
            result["max_ntilde_diff"], abs(generated["ntilde"] - ref["ntilde"])
        )

    for generated, ref in zip(generated_forces, ref_forces):
        generated_virial = torch.tensor(generated["virial"], dtype=DTYPE)
        ref_virial = torch.tensor(ref["virial"], dtype=DTYPE)
        result["max_virial_diff"] = max(
            result["max_virial_diff"],
            float(torch.max(torch.abs(generated_virial - ref_virial)).item()),
        )
        if len(generated["forces"]) != len(ref["forces"]):
            result["serial_mismatches"] += abs(len(generated["forces"]) - len(ref["forces"]))
            continue
        for (generated_serial, generated_force), (ref_serial, ref_force) in zip(
            generated["forces"], ref["forces"]
        ):
            if generated_serial != ref_serial:
                result["serial_mismatches"] += 1
            generated_tensor = torch.tensor(generated_force, dtype=DTYPE)
            ref_tensor = torch.tensor(ref_force, dtype=DTYPE)
            result["max_force_diff"] = max(
                result["max_force_diff"],
                float(torch.max(torch.abs(generated_tensor - ref_tensor)).item()),
            )

    result["passed"] = (
        not result["length_mismatch"]
        and result["n_mismatches"] == 0
        and result["serial_mismatches"] == 0
        and result["max_ntilde_diff"] <= ntilde_tol
        and result["max_force_diff"] <= force_tol
        and result["max_virial_diff"] <= virial_tol
    )
    result["ntilde_tol"] = ntilde_tol
    result["force_tol"] = force_tol
    result["virial_tol"] = virial_tol
    return result


def print_comparison_report(name: str, test_dir: Path, output_dir: Path, result: dict):
    status = "PASS" if result["passed"] else "FAIL"
    print(f"{name}: {status}")
    print(f"  test_dir: {test_dir}")
    print(f"  generated_dir: {output_dir}")
    print(f"  frames: {result['frames']}")
    print(f"  n_mismatches: {result['n_mismatches']}")
    print(f"  serial_mismatches: {result['serial_mismatches']}")
    print(f"  max_ntilde_diff: {result['max_ntilde_diff']:.8g} (tol {result['ntilde_tol']})")
    print(f"  max_force_diff: {result['max_force_diff']:.8g} (tol {result['force_tol']})")
    print(f"  max_virial_diff: {result['max_virial_diff']:.8g} (tol {result['virial_tol']})")

#!/usr/bin/env python3
"""Generate cylinder probe-volume references and compare them to existing refs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import torch

from probe_volume_ref_common import (
    DTYPE,
    REPO_ROOT,
    TARGET_INDICES,
    build_frame,
    compare_generated_to_ref,
    expected_frames,
    print_comparison_report,
    shifted_to_reference,
    switch_r,
    switch_x,
    write_generated_files,
)


NAME = "cylinder"
DEFAULT_TEST_DIR = REPO_ROOT / "test/indus/cylinder/bulk_water"
DEFAULT_OUTPUT_DIR = Path("/tmp/indus-reference-check/cylinder")

AXIS_INDEX = 1
R_MAX = 0.9
HEIGHT = 0.8
BASE = torch.tensor([1.1, 1.2, 2.3], dtype=DTYPE)
ZETA_MIN = float(BASE[AXIS_INDEX])
ZETA_MAX = ZETA_MIN + HEIGHT
CENTER = BASE.clone()
CENTER[AXIS_INDEX] += 0.5 * HEIGHT


def compute_frame(time: float, positions_nm, box_nm) -> dict:
    box = torch.tensor(box_nm, dtype=DTYPE)
    positions = torch.tensor(positions_nm[TARGET_INDICES], dtype=DTYPE, requires_grad=True)
    positions_box = torch.remainder(positions, box)
    positions_shifted = shifted_to_reference(positions_box, CENTER, box)

    transverse = [d for d in range(3) if d != AXIS_INDEX]
    radial_dx = positions_shifted[:, transverse] - CENTER[transverse]
    r = torch.linalg.norm(radial_dx, dim=1)
    h_r, htilde_r = switch_r(r, R_MAX)
    h_axis, htilde_axis = switch_x(positions_shifted[:, AXIS_INDEX], ZETA_MIN, ZETA_MAX)

    h = h_r * h_axis
    htilde = htilde_r * htilde_axis
    htilde = torch.where(htilde > 1.0, torch.ones_like(htilde), htilde)
    return build_frame(time, box_nm, positions, h, htilde)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test-dir", type=Path, default=DEFAULT_TEST_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ntilde-tol", type=float, default=1.0e-3)
    parser.add_argument("--force-tol", type=float, default=2.0e-4)
    parser.add_argument("--virial-tol", type=float, default=2.0e-4)
    args = parser.parse_args()

    frames = expected_frames(compute_frame)
    write_generated_files(args.output_dir, frames)
    result = compare_generated_to_ref(
        args.output_dir, args.test_dir, args.ntilde_tol, args.force_tol, args.virial_tol
    )
    print_comparison_report(NAME, args.test_dir, args.output_dir, result)
    if not result["passed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()

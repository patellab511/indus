#!/usr/bin/env python3
"""Generate sphere probe-volume references and compare them to existing refs."""

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
    minimum_image_vector,
    print_comparison_report,
    put_in_box,
    switch_x,
    write_generated_files,
)


NAME = "sphere"
DEFAULT_TEST_DIR = REPO_ROOT / "test/indus/sphere_near_box_edge/bulk_water"
DEFAULT_OUTPUT_DIR = Path("/tmp/indus-reference-check/sphere")

CENTER = torch.tensor([0.2, 0.3, 0.4], dtype=DTYPE)
R_MIN = -1.0
R_MAX = 1.0


def compute_frame(time: float, positions_nm, box_nm) -> dict:
    box = torch.tensor(box_nm, dtype=DTYPE)
    positions = torch.tensor(positions_nm[TARGET_INDICES], dtype=DTYPE, requires_grad=True)
    positions_box = put_in_box(positions, box)
    dx_center = minimum_image_vector(CENTER, positions_box, box)
    r = torch.linalg.norm(dx_center, dim=1)
    h, htilde = switch_x(r, R_MIN, R_MAX)
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

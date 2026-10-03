#!/usr/bin/env python3
"""Generate box probe-volume references and compare them to existing refs."""

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


NAME = "box"
DEFAULT_TEST_DIR = REPO_ROOT / "test/indus/box/bulk_water"
DEFAULT_OUTPUT_DIR = Path("/tmp/indus-reference-check/box")

X_RANGE = (1.0, 2.0)
Y_RANGE = (1.8, 3.1)
Z_RANGE = (3.5, 4.0)
RANGES = (X_RANGE, Y_RANGE, Z_RANGE)


def compute_frame(time: float, positions_nm, box_nm) -> dict:
    box = torch.tensor(box_nm, dtype=DTYPE)
    positions = torch.tensor(positions_nm[TARGET_INDICES], dtype=DTYPE, requires_grad=True)
    positions_box = put_in_box(positions, box)

    center = torch.tensor([0.5 * (low + high) for low, high in RANGES], dtype=DTYPE)
    half_lengths = [0.5 * (high - low) for low, high in RANGES]
    dx_center = minimum_image_vector(center, positions_box, box)

    h_components = []
    htilde_components = []
    for d, half_length in enumerate(half_lengths):
        h_d, htilde_d = switch_x(dx_center[:, d], -half_length, half_length)
        h_components.append(h_d)
        htilde_components.append(htilde_d)

    h = h_components[0] * h_components[1] * h_components[2]
    htilde = htilde_components[0] * htilde_components[1] * htilde_components[2]
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

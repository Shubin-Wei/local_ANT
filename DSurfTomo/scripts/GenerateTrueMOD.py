#!/usr/bin/env python3

"""Generate a checkerboard MOD.true file from an existing DSurfTomo MOD."""

import os

import numpy as np

from dsurftomo_utils import (
    print_banner,
    print_kv,
    print_section,
    prompt,
    prompt_choice,
    prompt_float,
)


def load_mod(path):
    with open(path, "r") as fobj:
        lines = [line.strip() for line in fobj if line.strip()]

    if not lines:
        raise ValueError("MOD file is empty: %s" % path)

    dep1 = np.array([float(item) for item in lines[0].split()], dtype=float)
    if dep1.size == 0:
        raise ValueError("No depth values found in MOD header: %s" % path)

    nz = int(dep1.size)
    data_rows = []
    nx = None
    for raw_line in lines[1:]:
        row = [float(item) for item in raw_line.split()]
        if not row:
            continue
        if nx is None:
            nx = len(row)
        elif len(row) != nx:
            raise ValueError("Inconsistent MOD row length in %s" % path)
        data_rows.append(row)

    if nx is None or not data_rows:
        raise ValueError("No velocity rows found in MOD file: %s" % path)
    if len(data_rows) % nz != 0:
        raise ValueError(
            "Cannot infer ny from MOD rows: %d rows for %d depth layers" % (len(data_rows), nz)
        )

    ny = len(data_rows) // nz
    mod = np.array(data_rows, dtype=float).reshape(nz, ny, nx)
    return dep1, mod


def resolve_mod_path(path):
    normalized = os.path.expanduser(path)
    if os.path.isdir(normalized):
        return os.path.join(normalized, "MOD")
    if path.endswith(os.sep):
        return os.path.join(normalized, "MOD")
    return normalized


def build_checkerboard_anomaly(nx, ny, nz, anosize, amplitude):
    x = np.arange(1, nx + 1, dtype=float)
    y = np.arange(1, ny + 1, dtype=float)
    anomaly_xy = np.outer(np.sin(anosize * y), np.sin(anosize * x))
    anomaly = np.broadcast_to(anomaly_xy, (nz, ny, nx)).copy()
    return anomaly * amplitude


def generate_true_mod(background_mod, anosize, amplitude):
    nz, ny, nx = background_mod.shape
    anomaly = build_checkerboard_anomaly(nx, ny, nz, anosize, amplitude)
    return background_mod + anomaly


def write_true_mod(path, true_mod):
    nz, ny, _ = true_mod.shape
    parent_dir = os.path.dirname(os.path.abspath(path))
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)
    with open(path, "w") as fobj:
        for k in range(nz):
            for j in range(ny):
                for value in true_mod[k, j, :]:
                    fobj.write("%7.3f" % value)
                fobj.write("\n")


def resolve_output_path(path):
    normalized = os.path.expanduser(path)
    if os.path.isdir(normalized):
        return os.path.join(normalized, "MOD.true")
    if path.endswith(os.sep):
        base = normalized.rstrip(os.sep)
        if os.path.exists(base) and not os.path.isdir(base):
            raise ValueError(
                "Output path %s ends with '/' but %s is an existing file, not a directory."
                % (path, base)
            )
        return os.path.join(normalized, "MOD.true")
    return normalized


def main():
    print_banner(
        "GenerateTrueMOD.py",
        [
            "Build a checkerboard MOD.true using an existing MOD as background model.",
            "Output format matches the MOD.true file expected by DSurfTomo synthetic runs.",
        ],
    )

    print_section("Inputs")
    mod_path = resolve_mod_path(prompt("background MOD path", "MOD"))
    if not os.path.exists(mod_path):
        raise ValueError("File not found: %s" % mod_path)
    if os.path.isdir(mod_path):
        raise ValueError("Expected a MOD file but got a directory: %s" % mod_path)

    dep1, background_mod = load_mod(mod_path)
    nz, ny, nx = background_mod.shape
    print_kv("background MOD", mod_path)
    print_kv("nx ny nz", "%d %d %d" % (nx, ny, nz))
    print_kv("depth range", "%.3f-%.3f km" % (dep1[0], dep1[-1]))

    print_section("Checkerboard Parameters")
    parameter_mode = prompt_choice(
        "checkerboard parameter mode (default/manual)",
        ["default", "manual"],
        default="default",
    )
    anosize = 0.5
    amplitude = 0.4
    if parameter_mode == "manual":
        anosize = prompt_float("checkerboard angular size", "%.3f" % anosize)
        amplitude = prompt_float("checkerboard amplitude (km/s)", "%.3f" % amplitude)
    else:
        print_kv("checkerboard angular size", "%.3f" % anosize)
        print_kv("checkerboard amplitude", "%.3f km/s" % amplitude)

    true_mod = generate_true_mod(background_mod, anosize=anosize, amplitude=amplitude)
    output_path = resolve_output_path(prompt("output MOD.true path", "MOD.true"))
    write_true_mod(output_path, true_mod)

    print_section("Complete")
    print_kv("MOD.true", output_path)
    print_kv("velocity range", "%.3f-%.3f km/s" % (float(np.min(true_mod)), float(np.max(true_mod))))


if __name__ == "__main__":
    main()

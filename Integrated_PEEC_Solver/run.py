"""
Integrated PEEC Solver
======================
Combines Python L/P matrix extraction with C++ PEEC frequency domain solver.

Input:  input/Node.txt, input/Branch.txt, input/Port.txt
Output: workspace/ directory (S/Z/Y parameter result files)

Usage:
    python run.py
"""

import os
import sys
import shutil
import subprocess
import numpy as np

import config
from core.kernels import cal_p_self, cal_p_oth, cal_l_oth_approx
from core.solver import PEECSolver
from utils.geometry import get_branch_info, split_grid
from utils.io_handler import (load_node_data, load_branch_data,
                               save_matrix_cpp_format, save_b2n,
                               convert_port_file)


def main():
    print("=" * 60)
    print("  Integrated PEEC Solver")
    print("  Node/Branch -> L/P Extraction -> Circuit Simulation")
    print("=" * 60)

    # Create workspace directory
    os.makedirs(config.WORKSPACE_DIR, exist_ok=True)

    input_dir = config.INPUT_DIR
    workspace = config.WORKSPACE_DIR

    node_file = os.path.join(input_dir, '68/Node.txt')
    branch_file = os.path.join(input_dir, '68/Branch.txt')
    port_file = os.path.join(input_dir, '68/Port.txt')

    # Verify input files exist
    for fpath, name in [(node_file, '68/Node.txt'), (branch_file, '68/Branch.txt'),
                         (port_file, '68/Port.txt')]:
        if not os.path.exists(fpath):
            print("  Error: {} not found at {}".format(name, os.path.abspath(fpath)))
            return

    # ================================================================
    # Step 1: Compute L and P Matrices
    # ================================================================
    print("\n[Step 1/3] Computing L and P matrices...")
    scale = config.UNIT_SCALE
    points, node_sizes = load_node_data(node_file, scale)
    connects, branch_sizes = load_branch_data(branch_file, node_sizes)

    num_nodes = len(points)
    num_branches = len(connects)
    print("  Mesh: {} nodes, {} branches".format(num_nodes, num_branches))

    # P matrix (coefficient of potential)
    print("  Computing P matrix...")
    P = np.zeros((num_nodes, num_nodes))
    for i in range(num_nodes):
        P[i, i] = cal_p_self(node_sizes[i])
        for j in range(i + 1, num_nodes):
            P[i, j] = P[j, i] = cal_p_oth(points[i], points[j], node_sizes[j])

    # L matrix (partial inductance)
    print("  Computing L matrix...")
    L = np.zeros((num_branches, num_branches))
    peec_solver = PEECSolver(gauss_order=config.DEFAULT_GAUSS_ORDER)

    for i in range(num_branches):
        region, axis, dims = get_branch_info(i, connects, branch_sizes, points)
        subs = split_grid(region, config.DEFAULT_DIV_X, config.DEFAULT_DIV_Y)

        total_A1 = total_A2 = total_A3 = 0.0
        for s_i in subs:
            for s_j in subs:
                _, t_a1, t_a2, t_a3 = peec_solver.compute_pair_integral(
                    s_i, s_j, axis, axis, dims, dims)
                total_A1 += t_a1
                total_A2 += t_a2
                total_A3 += t_a3

        L[i, i] = total_A1 #+ total_A2 + total_A3

        for j in range(i + 1, num_branches):
            L[i, j] = L[j, i] = cal_l_oth_approx(points, connects[i], connects[j])

        sys.stdout.write("\r  Progress: {}/{}".format(i + 1, num_branches))
        sys.stdout.flush()
    print()

    # ================================================================
    # Step 2: Format Output for C++ Solver
    # ================================================================
    print("\n[Step 2/3] Writing solver input files...")
    save_matrix_cpp_format(L, os.path.join(workspace, 'L.txt'))
    save_matrix_cpp_format(P, os.path.join(workspace, 'P.txt'))
    save_b2n(connects, os.path.join(workspace, 'B2N.txt'))
    convert_port_file(port_file, os.path.join(workspace, 'PORT.txt'))

    # Copy solver config files to workspace
    for fname in ['set.txt', 'DIELECTRIC.txt']:
        src = os.path.join(input_dir, fname)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(workspace, fname))

    print("  Files written to {}".format(os.path.abspath(workspace)))

    # ================================================================
    # Step 3: Run C++ Frequency Solver
    # ================================================================
    print("\n[Step 3/3] Running C++ frequency domain solver...")
    solver_exe = os.path.abspath(config.CPP_SOLVER_EXE)
    workspace_abs = os.path.abspath(workspace)

    if not os.path.exists(solver_exe):
        print("  Warning: C++ solver executable not found at:")
        print("    " + solver_exe)
        print("  Please compile the Quasi_Static_Solver project first.")
        return

    # Pass workspace path as arg; the modified C++ exe reads L/P/B2N/PORT
    # directly from this directory (no mesh file needed).
    env = os.environ.copy()
    extra_paths = []
    if hasattr(config, 'MKL_DLL_DIR') and os.path.isdir(config.MKL_DLL_DIR):
        extra_paths.append(config.MKL_DLL_DIR)
    extra_paths.append(os.path.dirname(solver_exe))
    env["PATH"] = os.pathsep.join(extra_paths) + os.pathsep + env.get("PATH", "")

    cmd = [solver_exe, workspace_abs]
    print("  Executing: {} \"{}\"".format(os.path.basename(solver_exe), workspace_abs))
    print("-" * 60)

    result = subprocess.run(cmd, env=env)

    print("-" * 60)
    if result.returncode == 0:
        print("\n" + "=" * 60)
        print("  Simulation completed successfully!")
        print("  Results directory: " + workspace_abs)
        print("  Output files:")
        print("    - map_Final.txt   (S parameters)")
        print("    - Z_in_Final.txt  (Z parameters)")
        print("    - Y_in_Final.txt  (Y parameters)")
        print("=" * 60)
    else:
        print("\n  Solver exited with return code {}".format(result.returncode))


if __name__ == "__main__":
    main()

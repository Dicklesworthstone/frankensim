//! A clamped quadratic-tetrahedron bar under a nodal force pulse and release.
//! Run `cargo run -p fs-solid --example spectral_elastic_motion`.
//!
//! Geometry and material values are explicitly synthetic: 1 x 0.2 x 0.2 m,
//! E=100 MPa, nu=0.27, rho=1000 kg/m^3. Twelve affine TET10 cells retain the
//! consistent mass and all 108 free coordinates. A triangular force pulse has
//! total peak 0.1 N, equally distributed over end-face z DOFs; it is a nodal
//! loading specification, not an assertion of exact uniform surface traction.
//!
//! The existing dense native eigensolver supplies the COMPLETE basis during
//! this small, bounded setup; it has no internal cancellation checkpoints.
//! The spectral owner validates that basis against the original K and M, then
//! advances every represented mode with exact affine forces. This example
//! demonstrates temporal integration, not spatial or experimental validation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::{CompleteSpectralSolid, CompleteSpectralSolidControls};
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetMaterialField};
use fs_sparse::Csr;

fn dense(matrix: &Csr) -> Vec<f64> {
    let n = matrix.nrows();
    let mut values = vec![0.0; n * n];
    for row in 0..n {
        let (columns, entries) = matrix.row(row);
        for (&column, &value) in columns.iter().zip(entries) {
            values[row * n + column] = value;
        }
    }
    values
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cancel = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &cancel,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 36,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let mut corners = Vec::new();
        for i in 0..=2 {
            for j in 0..=1 {
                for k in 0..=1 {
                    corners.push([0.5 * i as f64, 0.2 * j as f64, 0.2 * k as f64]);
                }
            }
        }
        let mut cells = Vec::new();
        for base in [0, 4] {
            let v = [
                base,
                base + 4,
                base + 2,
                base + 6,
                base + 1,
                base + 5,
                base + 3,
                base + 7,
            ];
            cells.extend([
                [v[0], v[1], v[3], v[7]],
                [v[0], v[3], v[2], v[7]],
                [v[0], v[2], v[6], v[7]],
                [v[0], v[6], v[4], v[7]],
                [v[0], v[4], v[5], v[7]],
                [v[0], v[5], v[1], v[7]],
            ]);
        }
        let (nodes, elements) = promote_p1_to_p2(&corners, &cells, None, &cx)?;
        let fixed: Vec<_> = nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0] == 0.0)
            .flat_map(|(node, _)| (0..3).map(move |axis| 3 * node + axis))
            .collect();
        let material = TetElasticMaterial::try_new(1000.0, 1e8, 0.27, [36; 32])?;
        let assembly = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget {
                maximum_nodes: 64,
                maximum_tetrahedra: 12,
                maximum_free_dofs: 128,
                minimum_scaled_jacobian: 1e-10,
            },
            quadrature: TetQuadrature {
                stiffness_points_per_dir: 3,
                mass_points_per_dir: 4,
            },
        }
        .assemble(&cx)?;
        let n = assembly.free_dofs.len();
        let modes =
            fs_modal::eigh_gen_dense(&dense(&assembly.stiffness), &dense(&assembly.mass), n)?;
        let eigenvalues: Vec<_> = modes.iter().map(|mode| mode.lambda).collect();
        let basis: Vec<_> = (0..n)
            .flat_map(|row| modes.iter().map(move |mode| mode.phi[row]))
            .collect();
        let controls = CompleteSpectralSolidControls {
            maximum_free_dofs: 128,
            maximum_preparation_bytes: 64 * 1024 * 1024,
            maximum_preparation_multiply_terms: 100_000_000,
            maximum_step_multiply_terms: 1_000_000,
            maximum_relative_operator_error: 1e-9,
            maximum_absolute_eigen_residual: 0.0,
            maximum_relative_state_error: 1e-10,
            maximum_work_error_j: 1e-10,
            maximum_steps: 240,
        };
        let mut solid = CompleteSpectralSolid::prepare_with_basis(
            &assembly,
            &eigenvalues,
            &basis,
            1.0 / 48000.0,
            &vec![0.0; n],
            &vec![0.0; n],
            controls,
            &cx,
        )?;
        let end_rows: Vec<_> = assembly
            .free_dofs
            .iter()
            .enumerate()
            .filter(|(_, dof)| **dof % 3 == 2 && nodes[**dof / 3][0] == 1.0)
            .map(|(row, _)| row)
            .collect();
        let pulse = |step: usize| -> f64 {
            match step {
                0..=40 => step as f64 / 40.0,
                41..=80 => (80 - step) as f64 / 40.0,
                _ => 0.0,
            }
        };
        eprintln!(
            "synthetic TET10 bar: {} elements, {n} free DOFs, {} retained modes",
            elements.len(),
            modes.len()
        );
        println!("time_s,mean_end_z_m,stored_j,kinetic_j,cumulative_work_j,balance_error_j");
        let mut start = vec![0.0; n];
        let mut end = vec![0.0; n];
        for step in 0..240 {
            for &row in &end_rows {
                start[row] = 0.1 * pulse(step) / end_rows.len() as f64;
                end[row] = 0.1 * pulse(step + 1) / end_rows.len() as f64;
            }
            let result = solid.advance_affine(&start, &end, &cx)?;
            let mean_end = end_rows
                .iter()
                .map(|&row| solid.displacement_m()[row])
                .sum::<f64>()
                / end_rows.len() as f64;
            println!(
                "{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e}",
                result.time_s,
                mean_end,
                result.energy.stored_after_j,
                result.energy.kinetic_after_j,
                solid.cumulative_work_j(),
                result.energy.cumulative_balance_error_j
            );
        }
        Ok(())
    })
}

//! Actual quadratic-element assembly -> complete native spectrum -> motion.
//! This is an integration test of the represented spatial pencil, not a mesh
//! convergence or measured-material claim.

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

fn quadratic(matrix: &Csr, x: &[f64]) -> f64 {
    let mut y = vec![0.0; x.len()];
    matrix.spmv(x, &mut y);
    0.5 * x.iter().zip(y).map(|(a, b)| a * b).sum::<f64>()
}

#[test]
fn g1_actual_tet10_bar_preserves_all_modes_and_original_force_work() {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 62, kernel_id: 37, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        // Two axial cells, each split into six consistently oriented tets.
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
            let v = [base, base + 4, base + 2, base + 6, base + 1, base + 5, base + 3, base + 7];
            cells.extend([
                [v[0], v[1], v[3], v[7]], [v[0], v[3], v[2], v[7]],
                [v[0], v[2], v[6], v[7]], [v[0], v[6], v[4], v[7]],
                [v[0], v[4], v[5], v[7]], [v[0], v[5], v[1], v[7]],
            ]);
        }
        let (nodes, elements) = promote_p1_to_p2(&corners, &cells, None, &cx).unwrap();
        let fixed: Vec<_> = nodes.iter().enumerate().filter(|(_, p)| p[0] == 0.0)
            .flat_map(|(node, _)| (0..3).map(move |axis| 3 * node + axis)).collect();
        let material = TetElasticMaterial::try_new(1000.0, 1e8, 0.27, [37; 32]).unwrap();
        let assembly = QuadraticTetProblem {
            nodes_m: &nodes, elements: &elements,
            materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed,
            budget: TetAssemblyBudget { maximum_nodes: 64, maximum_tetrahedra: 12,
                maximum_free_dofs: 128, minimum_scaled_jacobian: 1e-10 },
            quadrature: TetQuadrature { stiffness_points_per_dir: 3, mass_points_per_dir: 4 },
        }.assemble(&cx).unwrap();
        let n = assembly.free_dofs.len();
        assert_eq!(n, 108);
        assert!((0..n).any(|row| assembly.mass.row(row).1.iter().any(|&entry| entry < 0.0)),
            "TET10 consistent signed mass entries must survive assembly");
        let modes = fs_modal::eigh_gen_dense(&dense(&assembly.stiffness), &dense(&assembly.mass), n).unwrap();
        assert_eq!(modes.len(), n, "no modal truncation is admissible");
        assert!(modes.iter().all(|mode| mode.lambda > 0.0));
        let eigenvalues: Vec<_> = modes.iter().map(|mode| mode.lambda).collect();
        let basis: Vec<_> = (0..n).flat_map(|row| modes.iter().map(move |mode| mode.phi[row])).collect();
        // Excite both ends of the spectrum: a low-mode-only implementation
        // cannot satisfy the complete physical trajectory comparison.
        let initial_u: Vec<_> = (0..n).map(|row| 1e-6 * modes[0].phi[row]
            + 2e-8 * modes[n - 1].phi[row]).collect();
        let initial_v: Vec<_> = (0..n).map(|row| -1e-3 * modes[n - 1].phi[row]).collect();
        let controls = CompleteSpectralSolidControls {
            maximum_free_dofs: 128, maximum_preparation_bytes: 64 * 1024 * 1024,
            maximum_preparation_multiply_terms: 100_000_000,
            maximum_step_multiply_terms: 1_000_000,
            maximum_relative_operator_error: 1e-9,
            maximum_absolute_eigen_residual: 0.0,
            maximum_relative_state_error: 1e-10,
            maximum_work_error_j: 1e-10, maximum_steps: 128,
        };
        let dt = 1.0 / 48000.0;
        let mut solid = CompleteSpectralSolid::prepare_with_basis(&assembly,
            &eigenvalues, &basis, dt, &initial_u, &initial_v, controls, &cx).unwrap();
        let energy0 = quadratic(&assembly.stiffness, &initial_u) + quadratic(&assembly.mass, &initial_v);
        let zero = vec![0.0; n];
        let omega0 = modes[0].lambda.sqrt();
        let omega_last = modes[n - 1].lambda.sqrt();
        let mut maximum_free_energy_error = 0.0_f64;
        for _ in 0..64 {
            let report = solid.advance(&zero, &cx).unwrap();
            let t = report.time_s;
            let (s0, c0) = (omega0 * t).sin_cos();
            let (s1, c1) = (omega_last * t).sin_cos();
            let q_last = 2e-8 * c1 - 1e-3 * s1 / omega_last;
            let v_last = -2e-8 * omega_last * s1 - 1e-3 * c1;
            for row in 0..n {
                let exact_u = 1e-6 * c0 * modes[0].phi[row] + q_last * modes[n - 1].phi[row];
                let exact_v = -1e-6 * omega0 * s0 * modes[0].phi[row] + v_last * modes[n - 1].phi[row];
                assert!((solid.displacement_m()[row] - exact_u).abs() < 1e-14);
                assert!((solid.velocity_mps()[row] - exact_v).abs() < 1e-10);
            }
            let energy = quadratic(&assembly.stiffness, solid.displacement_m())
                + quadratic(&assembly.mass, solid.velocity_mps());
            maximum_free_energy_error = maximum_free_energy_error.max((energy - energy0).abs());
            assert!((energy - energy0).abs() < 1e-12 * energy0 + 1e-16);
        }
        // Physical force work is checked directly against the sparse original
        // operators, without using the spectral owner's energy calculation.
        let end_rows: Vec<_> = assembly.free_dofs.iter().enumerate()
            .filter(|(_, dof)| **dof % 3 == 2 && nodes[**dof / 3][0] == 1.0)
            .map(|(row, _)| row).collect();
        let mut force = vec![0.0; n];
        for &row in &end_rows { force[row] = 0.1 / end_rows.len() as f64; }
        let mut independent_work = 0.0;
        for step in 0..64 {
            let load = if step < 32 { &force } else { &zero };
            let before = solid.displacement_m().to_vec();
            let report = solid.advance(load, &cx).unwrap();
            let work: f64 = load.iter().enumerate()
                .map(|(row, force)| force * (solid.displacement_m()[row] - before[row])).sum();
            independent_work += work;
            assert!((report.energy.supplied_work_j - work).abs() < 1e-15);
            let energy = quadratic(&assembly.stiffness, solid.displacement_m())
                + quadratic(&assembly.mass, solid.velocity_mps());
            assert!((energy - energy0 - independent_work).abs() < 1e-14);
        }
        assert!(independent_work.abs() > 1e-10, "forcing must exercise live motion");
        println!("actual TET10: elements={}, free_dofs={n}, retained_modes={}, lowest_hz={:e}, highest_hz={:e}, maximum_free_energy_error_j={maximum_free_energy_error:e}, force_work_j={independent_work:e}",
            elements.len(), modes.len(), omega0 / core::f64::consts::TAU, omega_last / core::f64::consts::TAU);
    });
}

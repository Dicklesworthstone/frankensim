//! Original pressure-only Newton blocks and independently known volume fields.

use super::*;
use crate::linear3::tet10::promote_p1_to_p2;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::hyper::{HyperelasticModel, OgdenTerm, Volumetric};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

fn run(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..Default::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: 12.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .position(|&c| c == column)
        .map_or(0.0, |i| values[i])
}

#[test]
fn prescribed_curved_maxwell_pressure_block_matches_the_full_original_tangent() {
    run(|cx| {
        let corners = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (mut nodes, elements) = promote_p1_to_p2(&corners, &[[0, 1, 2, 3]], None, cx).unwrap();
        for p in &mut nodes {
            p[0] += 0.125 * p[0] * p[0];
        }
        let displacement: Vec<_> = nodes
            .iter()
            .map(|p| [0.08 * p[1] * p[2], 0.03 * p[0], -0.02 * p[2]])
            .collect();
        let prescribed: Vec<_> = displacement
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .collect();
        for volumetric in [
            Volumetric::Quadratic { kappa: 120.0 },
            Volumetric::LogSquared { kappa: 120.0 },
        ] {
            let elastic = Hyperelastic::new(
                HyperelasticModel::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 12.0,
                        alpha: 2.0,
                    }],
                    volumetric,
                },
                4.0,
            )
            .unwrap();
            let law = ViscoHyperelastic::new(
                elastic,
                vec![FiniteStrainBranch {
                    energy: BranchEnergy::NeoHookean { mu_pa: 36.0 },
                    tau_s: 0.07,
                }],
            )
            .unwrap();
            let old_f = [1.0, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
            let history = law
                .step_remainder_response(&law.initial_remainder_state(), &old_f, 0.02)
                .unwrap()
                .state;
            let frozen = law.freeze_remainder(&history, 0.01).unwrap();
            let problem = MixedUpProblem {
                nodes_m: &nodes,
                elements: &elements,
                material: MixedLawField::Uniform(&frozen),
                volumetric: VolumetricTreatment::OriginalMaterial,
                prescribed_m: &prescribed,
                nodal_forces_n: &[],
                body_force_n_m3: [0.0; 3],
                pressure_gauge: None,
                budget: TetAssemblyBudget::standard(),
                quadrature_points_per_dir: 4,
                maximum_quadrature_points: 64,
                settings: settings(),
            };
            let prepared = problem.prepare(cx).unwrap();
            let pressure = [-1.0, 3.0, 0.5, -0.25];
            let only = prepared
                .pressure_tangent(&displacement, &pressure, cx)
                .unwrap();
            let full = prepared
                .evaluate(&displacement, &pressure, true, cx)
                .unwrap()
                .tangent
                .unwrap();
            let nu = 3 * nodes.len();
            for row in 0..nu {
                assert!(only.row(row).1.is_empty());
            }
            for a in 0..4 {
                for b in 0..4 {
                    assert_eq!(entry(&only, nu + a, nu + b), entry(&full, nu + a, nu + b));
                }
            }
            // Independent derivative of the original pressure residual. The
            // frozen Maxwell remainder and actual curved reference weights
            // remain present in both residual evaluations.
            let epsilon = 1e-4;
            for b in 0..4 {
                let mut plus = pressure;
                let mut minus = pressure;
                plus[b] += epsilon;
                minus[b] -= epsilon;
                let plus = prepared.evaluate(&displacement, &plus, false, cx).unwrap();
                let minus = prepared.evaluate(&displacement, &minus, false, cx).unwrap();
                for a in 0..4 {
                    let derivative = (plus.pressure_residual_m3[a] - minus.pressure_residual_m3[a])
                        / (2.0 * epsilon);
                    assert!((derivative - entry(&only, nu + a, nu + b)).abs() < 2e-11);
                }
            }
        }
    });
}

#[test]
fn prescribed_quadratic_volume_solves_its_original_linear_pressure_field() {
    run(|cx| {
        let corners = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (nodes, elements) = promote_p1_to_p2(&corners, &[[0, 1, 2, 3]], None, cx).unwrap();
        let law = Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6.0,
                c01: 0.0,
                kappa: 120.0,
            },
            4.0,
        )
        .unwrap();
        // X'=X+.1X², Y'=Y+.05XZ: det F=1+.2X, so the original
        // quadratic bulk law has the EXACT P1 field p=24X, not uniform p.
        let displacement: Vec<_> = nodes
            .iter()
            .map(|p| [0.1 * p[0] * p[0], 0.05 * p[0] * p[2], 0.0])
            .collect();
        let prescribed: Vec<_> = displacement
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .collect();
        let problem = MixedUpProblem {
            nodes_m: &nodes,
            elements: &elements,
            material: MixedLawField::Uniform(&law),
            volumetric: VolumetricTreatment::OriginalMaterial,
            prescribed_m: &prescribed,
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            pressure_gauge: None,
            budget: TetAssemblyBudget::standard(),
            quadrature_points_per_dir: 4,
            maximum_quadrature_points: 64,
            settings: settings(),
        };
        let prepared = problem.prepare(cx).unwrap();
        let solution = prepared.solve(cx).unwrap();
        assert_eq!(solution.displacement_m, displacement);
        for (&node, &pressure) in solution.pressure_nodes.iter().zip(&solution.pressure_pa) {
            assert!((pressure - 24.0 * nodes[node][0]).abs() < 2e-8);
        }
        let tangent = prepared
            .pressure_tangent(&displacement, &solution.pressure_pa, cx)
            .unwrap();
        let nu = 3 * nodes.len();
        for a in 0..4 {
            for b in 0..4 {
                let expected = -1.0 / (6.0 * 120.0 * if a == b { 10.0 } else { 20.0 });
                assert!((entry(&tangent, nu + a, nu + b) - expected).abs() < 2e-16);
            }
        }
    });
}

//! G1/G3 consistent reference inertia and signed prestressed truss vibration.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::{ModalError, MultiSliceReport, SliceOptions, WindowPlan};
use fs_sparse::Csr;
use fs_truss::elastic::modal::{
    ModalAssemblyControls, ModalInertia, TrussModalError, TrussModalPencil, assemble,
};
use fs_truss::elastic::{ElasticError, ElasticImposedState};
use fs_truss::{GroundLimits, GroundStructure, LayoutCase, TrussConstructionError};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs().max(1e-12),
        "{actual:e} != {expected:e}",
    );
}

fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (columns, values) = matrix.row(row);
    columns.binary_search(&column).map_or(0.0, |i| values[i])
}

fn action(matrix: &Csr, vector: &[f64]) -> Vec<f64> {
    (0..matrix.nrows())
        .map(|row| {
            let (columns, values) = matrix.row(row);
            columns
                .iter()
                .zip(values)
                .map(|(&j, &v)| v * vector[j])
                .sum()
        })
        .collect()
}

fn window_plan(windows: Vec<(f64, f64)>) -> (WindowPlan, SliceOptions) {
    let plan = WindowPlan {
        maximum_windows: windows.len(),
        windows,
        threads: 1,
        maximum_bytes: 64 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-7,
        coordinates: None,
    };
    let mut options = SliceOptions::default();
    options.ldlt.max_factor_bytes = Some(8 * 1024 * 1024);
    options.mass_diagonal_equilibration = true;
    (plan, options)
}

fn slice(pencil: &TrussModalPencil, windows: Vec<(f64, f64)>, cx: &Cx<'_>) -> MultiSliceReport {
    let (plan, options) = window_plan(windows);
    let result = pencil.slice_windows(&plan, &options, cx).unwrap();
    for report in &result.reports {
        assert_eq!(report.modes.len(), report.expected);
        for mode in &report.modes {
            let mass_phi = action(pencil.mass(), &mode.phi);
            let norm = mode.phi.iter().zip(mass_phi).map(|(p, mp)| p * mp).sum();
            near(norm, 1.0);
            assert!(mode.residual.is_finite());
        }
    }
    result
}

fn bar(cx: &Cx<'_>, transverse_free: bool) -> (GroundStructure, LayoutCase) {
    let ground = GroundStructure::try_from_parts(
        &[[0.0, 0.0], [2.0, 0.0]],
        &[(0, 1)],
        &[2.0],
        GroundLimits::default(),
        cx,
    )
    .unwrap();
    let case = LayoutCase::try_new(
        vec![[true; 2], [false, !transverse_free]],
        vec![[0.0; 2]; 2],
        2,
    )
    .unwrap();
    (ground, case)
}

#[test]
fn g1_bar_consistent_mass_frequency_and_explicit_joint_mass() {
    with_cx(|cx, _| {
        let (ground, case) = bar(cx, false);
        for area in [0.5, 2.0] {
            for density in [3.0, 6.0] {
                for joints in [None, Some([5.0, 2.0])] {
                    let pencil = assemble(
                        &ground,
                        &case,
                        &[area],
                        &[12.0],
                        ElasticImposedState::default(),
                        ModalInertia {
                            member_density_kg_m3: &[density],
                            joint_mass_kg: joints.as_ref().map(|m| &m[..]),
                        },
                        ModalAssemblyControls::default(),
                        cx,
                    )
                    .unwrap();
                    let member_mass = density * area * 2.0;
                    let added_free = joints.map_or(0.0, |m| m[1]);
                    let reduced_mass = member_mass / 3.0 + added_free;
                    near(entry(pencil.mass(), 0, 0), reduced_mass);
                    near(
                        pencil.total_reference_mass_kg(),
                        member_mass + joints.map_or(0.0, |m| m[0] + m[1]),
                    );
                    let expected = 12.0 * area / 2.0 / reduced_mass;
                    let result = slice(&pencil, vec![(0.5 * expected, 1.5 * expected)], cx);
                    near(result.reports[0].modes[0].lambda, expected);
                    assert_eq!(pencil.dof_map(), &[None, None, Some(0), None]);
                }
            }
        }
    });
}

#[test]
fn g1_serial_member_mass_coupling_matches_the_two_element_spectrum() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
            &[(0, 1), (1, 2)],
            &[1.0; 2],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        let case = LayoutCase::try_new(
            vec![[true; 2], [false, true], [false, true]],
            vec![[0.0; 2]; 3],
            3,
        )
        .unwrap();
        let pencil = assemble(
            &ground,
            &case,
            &[0.5; 2],
            &[12.0; 2],
            ElasticImposedState::default(),
            ModalInertia {
                member_density_kg_m3: &[3.0; 2],
                joint_mass_kg: None,
            },
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        // M=[[1,1/4],[1/4,1/2]], K=[[12,-6],[-6,6]].
        near(entry(pencil.mass(), 0, 0), 1.0);
        near(entry(pencil.mass(), 0, 1), 0.25);
        near(entry(pencil.mass(), 1, 1), 0.5);
        let result = slice(&pencil, vec![(0.1, 40.0)], cx);
        assert_eq!(result.reports[0].modes.len(), 2);
        for (mode, sign) in result.reports[0].modes.iter().zip([-1.0, 1.0]) {
            near(mode.lambda, (120.0 + sign * 72.0 * 2.0_f64.sqrt()) / 7.0);
        }
    });
}

#[test]
fn g3_rigid_transverse_inertia_and_rotated_member_kinetic_energy_are_retained() {
    with_cx(|cx, _| {
        for angle in [0.0_f64, 0.37, 1.1] {
            let (s, c) = angle.sin_cos();
            let ground = GroundStructure::try_from_parts(
                &[[0.0, 0.0], [2.0 * c, 2.0 * s]],
                &[(0, 1)],
                &[2.0],
                GroundLimits::default(),
                cx,
            )
            .unwrap();
            let case = LayoutCase::try_new(vec![[false; 2]; 2], vec![[0.0; 2]; 2], 2).unwrap();
            let pencil = assemble(
                &ground,
                &case,
                &[0.5],
                &[12.0],
                ElasticImposedState::default(),
                ModalInertia {
                    member_density_kg_m3: &[3.0],
                    joint_mass_kg: None,
                },
                ModalAssemblyControls::default(),
                cx,
            )
            .unwrap();
            let velocity = [1.0, 2.0, 1.0, 2.0];
            let mass_v = action(pencil.mass(), &velocity);
            let kinetic = 0.5
                * velocity
                    .iter()
                    .zip(mass_v)
                    .map(|(v, mv)| v * mv)
                    .sum::<f64>();
            near(kinetic, 0.5 * 3.0 * 5.0);
            let transverse = [-s, c, -s, c];
            let mass_v = action(pencil.mass(), &transverse);
            near(
                transverse.iter().zip(mass_v).map(|(v, mv)| v * mv).sum(),
                3.0,
            );
            assert_eq!(pencil.dof_map(), &[Some(0), Some(1), Some(2), Some(3)]);
        }
        let (ground, case) = bar(cx, true);
        let pencil = assemble(
            &ground,
            &case,
            &[0.5],
            &[12.0],
            ElasticImposedState::default(),
            ModalInertia {
                member_density_kg_m3: &[3.0],
                joint_mass_kg: None,
            },
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        assert_eq!(entry(pencil.tangent(), 1, 1), 0.0);
        near(entry(pencil.mass(), 1, 1), 1.0);
        // Query the axial mode. Inertia still counts the zero transverse mode
        // below this window; it was not silently discarded during assembly.
        let result = slice(&pencil, vec![(2.0, 4.0)], cx);
        assert_eq!(result.reports[0].below_low, 1);
        near(result.reports[0].modes[0].lambda, 3.0);
    });
}

#[test]
fn g3_thermal_settlement_and_dead_load_prestress_shift_signed_rotated_modes() {
    with_cx(|cx, _| {
        for angle in [0.0_f64, 0.37] {
            let (s, c) = angle.sin_cos();
            let rotate = |p: [f64; 2]| [c * p[0] - s * p[1], s * p[0] + c * p[1]];
            let ground = GroundStructure::try_from_parts(
                &[[-2.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, -1.0]].map(rotate),
                &[(0, 1), (1, 2), (1, 3)],
                &[2.0, 1.0, 1.0],
                GroundLimits::default(),
                cx,
            )
            .unwrap();
            for (strain, settlement, force) in
                [(0.0, 0.0, 0.0), (1e-5, 2e-5, 20.0), (1e-4, 0.0, 0.0)]
            {
                let case = LayoutCase::try_new(
                    vec![[true; 2], [false; 2], [true; 2], [true; 2]],
                    vec![[0.0; 2], rotate([force, 0.0]), [0.0; 2], [0.0; 2]],
                    4,
                )
                .unwrap();
                let pencil = assemble(
                    &ground,
                    &case,
                    &[1.0; 3],
                    &[2e6, 1e6, 100.0],
                    ElasticImposedState {
                        member_stress_free_strain: Some(&[strain, strain, 0.0]),
                        prescribed_displacements_m: Some(&[
                            [0.0; 2],
                            [0.0; 2],
                            rotate([settlement, 0.0]),
                            [0.0; 2],
                        ]),
                    },
                    ModalInertia {
                        member_density_kg_m3: &[3.0, 6.0, 9.0],
                        joint_mass_kg: None,
                    },
                    ModalAssemblyControls::default(),
                    cx,
                )
                .unwrap();
                let n0 = 0.5e6 * (settlement - 3.0 * strain);
                let expected = (100.0 + 1.5 * n0 - force / 4.0) / 7.0;
                let radius = expected.abs() / 3.0;
                let result = slice(&pencil, vec![(expected - radius, expected + radius)], cx);
                let mode = &result.reports[0].modes[0];
                near(mode.lambda, expected);
                assert_eq!(mode.lambda.is_sign_negative(), expected.is_sign_negative());
                near(
                    (-s * mode.phi[0] + c * mode.phi[1]).abs(),
                    1.0 / 7.0_f64.sqrt(),
                );
                near(pencil.total_reference_mass_kg(), 21.0);
                near(entry(pencil.mass(), 0, 0), 7.0);
                near(entry(pencil.mass(), 1, 1), 7.0);
                for axis in 0..2 {
                    near(
                        pencil.equilibrium().displacements_m[2][axis],
                        rotate([settlement, 0.0])[axis],
                    );
                }
            }
        }
    });
}

#[test]
fn g4_pruned_massless_coordinates_are_preserved_and_native_mass_admission_refuses() {
    with_cx(|cx, _| {
        let (ground, case) = bar(cx, true);
        let pencil = assemble(
            &ground,
            &case,
            &[0.0],
            &[12.0],
            ElasticImposedState::default(),
            ModalInertia {
                member_density_kg_m3: &[3.0],
                joint_mass_kg: None,
            },
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        assert_eq!(pencil.dof_map(), &[None, None, Some(0), Some(1)]);
        assert_eq!(pencil.mass().nrows(), 2);
        assert_eq!(pencil.mass().nnz(), 0);
        assert_eq!(pencil.tangent().nnz(), 0);
        assert_eq!(pencil.total_reference_mass_kg(), 0.0);
        let (plan, options) = window_plan(vec![(-1.0, 1.0)]);
        assert!(matches!(
            pencil.slice_windows(&plan, &options, cx),
            Err(TrussModalError::Modal(_))
        ));
    });
}

#[test]
fn g4_inertia_work_and_cancellation_refuse_without_a_partial_pencil_or_spectrum() {
    with_cx(|cx, gate| {
        let (ground, case) = bar(cx, false);
        let controls = ModalAssemblyControls::default();
        for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                assemble(
                    &ground,
                    &case,
                    &[0.5],
                    &[12.0],
                    ElasticImposedState::default(),
                    ModalInertia {
                        member_density_kg_m3: &[density],
                        joint_mass_kg: None
                    },
                    controls,
                    cx,
                )
                .is_err()
            );
        }
        for masses in [&[1.0][..], &[-1.0, 0.0][..], &[0.0, f64::NAN][..]] {
            assert!(
                assemble(
                    &ground,
                    &case,
                    &[0.5],
                    &[12.0],
                    ElasticImposedState::default(),
                    ModalInertia {
                        member_density_kg_m3: &[3.0],
                        joint_mass_kg: Some(masses)
                    },
                    controls,
                    cx,
                )
                .is_err()
            );
        }
        let inertia = ModalInertia {
            member_density_kg_m3: &[3.0],
            joint_mass_kg: None,
        };
        assert!(matches!(
            assemble(
                &ground,
                &case,
                &[0.5],
                &[12.0],
                ElasticImposedState::default(),
                inertia,
                ModalAssemblyControls {
                    maximum_staged_entries: 1,
                    ..controls
                },
                cx,
            ),
            Err(TrussModalError::Elastic(ElasticError::Construction(
                TrussConstructionError::WorkBudget { .. }
            )))
        ));
        let pencil = assemble(
            &ground,
            &case,
            &[0.5],
            &[12.0],
            ElasticImposedState::default(),
            inertia,
            controls,
            cx,
        )
        .unwrap();
        let (mut plan, options) = window_plan(vec![(2.0, 4.0)]);
        plan.maximum_bytes = 1;
        assert!(matches!(
            pencil.slice_windows(&plan, &options, cx),
            Err(TrussModalError::Modal(ModalError::MemoryBudget { .. }))
        ));
        plan.maximum_bytes = 64 * 1024 * 1024;
        gate.request();
        assert!(matches!(
            pencil.slice_windows(&plan, &options, cx),
            Err(TrussModalError::Modal(ModalError::Cancelled))
        ));
        assert!(matches!(
            assemble(
                &ground,
                &case,
                &[0.5],
                &[12.0],
                ElasticImposedState::default(),
                inertia,
                controls,
                cx,
            ),
            Err(TrussModalError::Elastic(ElasticError::Construction(
                TrussConstructionError::Cancelled { .. }
            )))
        ));
    });
}

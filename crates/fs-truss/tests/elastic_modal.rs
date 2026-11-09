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
// G1: sizing sensitivities retain actual mass and use one isolated native mode.
#[test]
fn g1_axial_mode_design_derivatives_include_member_and_joint_mass() {
    use fs_truss::elastic::modal::sensitivity::{ModalSensitivityControls, PreparedModalDesign};
    with_cx(|cx, _| {
        let (ground, case) = bar(cx, false);
        let areas = [0.5];
        let young = [12.0];
        let density = [3.0];
        let joints = [5.0, 2.0];
        let design = PreparedModalDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            ElasticImposedState::default(),
            ModalInertia {
                member_density_kg_m3: &density,
                joint_mass_kg: Some(&joints),
            },
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        let (plan, options) = window_plan(vec![(0.5, 1.5)]);
        let result = design
            .solve_window(&plan, &options, ModalSensitivityControls::default(), cx)
            .unwrap();
        // K=3, M=1+2=3, lambda=1. The member's area affects both K and M.
        near(result.eigenvalue_s2, 1.0);
        let member = result.members[0].unwrap();
        near(member.area_derivative, 4.0 / 3.0);
        near(member.young_modulus_derivative, 1.0 / 12.0);
        near(member.density_derivative, -1.0 / 9.0);
        assert_eq!(member.stress_free_strain_derivative, 0.0);
        near(result.joint_mass_derivative[1], -1.0 / 3.0);
        assert_eq!(result.joint_mass_derivative[0], 0.0);
        // lambda = (EA/L)/(rho*A*L/3 + m_joint): both length factors vary.
        near(result.node_coordinate_gradient[0][0], 2.0 / 3.0);
        near(result.node_coordinate_gradient[1][0], -2.0 / 3.0);
        assert_eq!(result.node_coordinate_gradient[0][1], 0.0);
        assert_eq!(result.node_coordinate_gradient[1][1], 0.0);
        assert_eq!(result.load_gradient, vec![[0.0; 2]; 2]);
        assert_eq!(result.adjoint_iterations, 0);
        // Reusing the same prepared design repeats no physical equilibrium.
        let again = design
            .solve_window(&plan, &options, ModalSensitivityControls::default(), cx)
            .unwrap();
        assert_eq!(result.members, again.members);
        assert_eq!(design.pencil().equilibrium().iterations, 0);
    });
}

#[allow(clippy::too_many_arguments)]
fn transverse_eigenvalue(
    ground: &GroundStructure,
    case: &LayoutCase,
    areas: &[f64],
    young: &[f64],
    strain: &[f64],
    motion: &[[f64; 2]],
    density: &[f64],
    joints: &[f64],
    angle: f64,
    _cx: &Cx<'_>,
) -> f64 {
    // Independent series-spring equilibrium, followed by exact transverse
    // geometric stiffness. Orthogonal member families retain these principal
    // directions under every tested physical-input perturbation.
    let axial = [angle.cos(), angle.sin()];
    let project = |v: [f64; 2]| v[0] * axial[0] + v[1] * axial[1];
    let lengths = ground.lengths();
    let kl = young[0] * areas[0] / lengths[0];
    let kr = young[1] * areas[1] / lengths[1];
    let left = project(motion[0]) + lengths[0] * strain[0];
    let right = project(motion[2]) - lengths[1] * strain[1];
    let center = (project(case.loads()[1]) + kl * left + kr * right) / (kl + kr);
    let nl = kl * (center - left);
    let nr = kr * (right - center);
    let stiffness = young[2] * areas[2] / lengths[2] + nl / lengths[0] + nr / lengths[1];
    let mass = density
        .iter()
        .zip(areas)
        .zip(lengths)
        .map(|((&rho, &area), &length)| rho * area * length / 3.0)
        .sum::<f64>()
        + joints[1];
    stiffness / mass
}

fn derivative_near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-5 * expected.abs().max(1e-5),
        "{actual:e} != {expected:e}",
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_prestress_design_adjoint_matches_load_thermal_support_and_inertia_differences() {
    use fs_truss::elastic::modal::sensitivity::{ModalSensitivityControls, PreparedModalDesign};
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
            let case = LayoutCase::try_new(
                vec![[true; 2], [false; 2], [true; 2], [true; 2]],
                vec![[0.0; 2], rotate([30.0, 0.005]), [0.0; 2], [0.0; 2]],
                4,
            )
            .unwrap();
            let areas = [1.1, 0.9, 1.3];
            let young = [2e6, 1e6, 100.0];
            let strain = [1e-5, 2e-5, 0.0];
            let motion = [[-1e-5, 0.0], [0.0; 2], [4e-5, 0.0], [0.0, 2e-6]].map(rotate);
            let density = [3.0, 6.0, 9.0];
            let joints = [0.2, 0.7, 0.3, 0.5];
            let design = PreparedModalDesign::try_new(
                &ground,
                &case,
                &areas,
                &young,
                ElasticImposedState {
                    member_stress_free_strain: Some(&strain),
                    prescribed_displacements_m: Some(&motion),
                },
                ModalInertia {
                    member_density_kg_m3: &density,
                    joint_mass_kg: Some(&joints),
                },
                ModalAssemblyControls::default(),
                cx,
            )
            .unwrap();
            let expected = transverse_eigenvalue(
                &ground, &case, &areas, &young, &strain, &motion, &density, &joints, angle, cx,
            );
            assert!(expected > 0.0);
            let (plan, options) = window_plan(vec![(0.5 * expected, 1.5 * expected)]);
            let result = design
                .solve_window(&plan, &options, ModalSensitivityControls::default(), cx)
                .unwrap();
            near(result.eigenvalue_s2, expected);
            let mode = &result.spectrum.reports[0].modes[0];
            let kphi = action(design.pencil().tangent(), &mode.phi);
            let mphi = action(design.pencil().mass(), &mode.phi);
            let residual = kphi
                .iter()
                .zip(mphi)
                .fold(0.0_f64, |r, (k, m)| r.hypot(k - mode.lambda * m));
            assert!(residual < 1e-7 * mode.lambda.abs());
            assert!(result.adjoint_relative_equilibrium_residual < 1e-9);
            assert!(result.load_gradient[1].iter().any(|v| v.abs() > 1e-3));
            assert!(
                result.members[0]
                    .unwrap()
                    .stress_free_strain_derivative
                    .abs()
                    > 1e4
            );
            for index in 0..3 {
                let derivative = result.members[index].unwrap();
                for field in 0..4 {
                    let mut ap = areas;
                    let mut am = areas;
                    let mut ep = young;
                    let mut em = young;
                    let mut sp = strain;
                    let mut sm = strain;
                    let mut rp = density;
                    let mut rm = density;
                    let (step, predicted) = match field {
                        0 => {
                            let h = 1e-5 * areas[index];
                            ap[index] += h;
                            am[index] -= h;
                            (h, derivative.area_derivative)
                        }
                        1 => {
                            let h = 1e-5 * young[index];
                            ep[index] += h;
                            em[index] -= h;
                            (h, derivative.young_modulus_derivative)
                        }
                        2 => {
                            sp[index] += 1e-9;
                            sm[index] -= 1e-9;
                            (1e-9, derivative.stress_free_strain_derivative)
                        }
                        _ => {
                            let h = 1e-5 * density[index];
                            rp[index] += h;
                            rm[index] -= h;
                            (h, derivative.density_derivative)
                        }
                    };
                    let difference = (transverse_eigenvalue(
                        &ground, &case, &ap, &ep, &sp, &motion, &rp, &joints, angle, cx,
                    ) - transverse_eigenvalue(
                        &ground, &case, &am, &em, &sm, &motion, &rm, &joints, angle, cx,
                    )) / (2.0 * step);
                    derivative_near(predicted, difference);
                }
            }
            for node in 0..4 {
                let mut plus = joints;
                let mut minus = joints;
                let step = 1e-5;
                plus[node] += step;
                minus[node] -= step;
                let difference = (transverse_eigenvalue(
                    &ground, &case, &areas, &young, &strain, &motion, &density, &plus, angle, cx,
                ) - transverse_eigenvalue(
                    &ground, &case, &areas, &young, &strain, &motion, &density, &minus, angle, cx,
                )) / (2.0 * step);
                derivative_near(result.joint_mass_derivative[node], difference);
                for axis in 0..2 {
                    let mut plus = if node == 1 {
                        case.loads().to_vec()
                    } else {
                        motion.to_vec()
                    };
                    let mut minus = plus.clone();
                    let step = if node == 1 { 1e-3 } else { 1e-9 };
                    plus[node][axis] += step;
                    minus[node][axis] -= step;
                    let (difference, predicted) = if node == 1 {
                        let cp = LayoutCase::try_new(case.supported().to_vec(), plus, 4).unwrap();
                        let cm = LayoutCase::try_new(case.supported().to_vec(), minus, 4).unwrap();
                        (
                            (transverse_eigenvalue(
                                &ground, &cp, &areas, &young, &strain, &motion, &density, &joints,
                                angle, cx,
                            ) - transverse_eigenvalue(
                                &ground, &cm, &areas, &young, &strain, &motion, &density, &joints,
                                angle, cx,
                            )) / (2.0 * step),
                            result.load_gradient[node][axis],
                        )
                    } else {
                        (
                            (transverse_eigenvalue(
                                &ground, &case, &areas, &young, &strain, &plus, &density, &joints,
                                angle, cx,
                            ) - transverse_eigenvalue(
                                &ground, &case, &areas, &young, &strain, &minus, &density, &joints,
                                angle, cx,
                            )) / (2.0 * step),
                            result.prescribed_displacement_gradient[node][axis],
                        )
                    };
                    derivative_near(predicted, difference);
                }
            }
            let no_iterations = ModalSensitivityControls {
                elastic: fs_truss::elastic::ElasticControls {
                    maximum_iterations: 0,
                    ..Default::default()
                },
                ..Default::default()
            };
            assert!(matches!(
                design.solve_window(&plan, &options, no_iterations, cx),
                Err(TrussModalError::Elastic(ElasticError::NotConverged { .. }))
            ));
        }
    });
}

#[test]
fn g1_consistent_off_diagonal_mass_sizing_derivative_matches_the_exact_pencil() {
    use fs_truss::elastic::modal::sensitivity::{ModalSensitivityControls, PreparedModalDesign};
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
        let areas = [0.5, 0.8];
        let young = [12.0, 19.0];
        let density = [3.0, 5.0];
        let exact = |a: &[f64], rho: &[f64]| {
            let k0 = young[0] * a[0];
            let k1 = young[1] * a[1];
            let m0 = rho[0] * a[0];
            let m1 = rho[1] * a[1];
            let aa = (m0 + m1) * m1 / 9.0 - m1 * m1 / 36.0;
            let bb = (k0 + k1) * m1 / 3.0 + k1 * (m0 + m1) / 3.0 + k1 * m1 / 3.0;
            let cc = k0 * k1;
            2.0 * cc / (bb + (bb * bb - 4.0 * aa * cc).sqrt())
        };
        let expected = exact(&areas, &density);
        let design = PreparedModalDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            ElasticImposedState::default(),
            ModalInertia {
                member_density_kg_m3: &density,
                joint_mass_kg: None,
            },
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        let (plan, options) = window_plan(vec![(0.5 * expected, 1.5 * expected)]);
        let result = design
            .solve_window(&plan, &options, ModalSensitivityControls::default(), cx)
            .unwrap();
        near(result.eigenvalue_s2, expected);
        for index in 0..2 {
            let mut plus = areas;
            let mut minus = areas;
            plus[index] += 1e-5;
            minus[index] -= 1e-5;
            derivative_near(
                result.members[index].unwrap().area_derivative,
                (exact(&plus, &density) - exact(&minus, &density)) / 2e-5,
            );
            let mut plus = density;
            let mut minus = density;
            plus[index] += 1e-5;
            minus[index] -= 1e-5;
            derivative_near(
                result.members[index].unwrap().density_derivative,
                (exact(&areas, &plus) - exact(&areas, &minus)) / 2e-5,
            );
        }
    });
}

#[test]
fn g4_signed_mode_derivatives_refuse_nonisolated_windows_budgets_and_cancellation() {
    use fs_truss::elastic::modal::sensitivity::{ModalSensitivityControls, PreparedModalDesign};
    with_cx(|cx, gate| {
        let (ground, case) = bar(cx, true);
        let areas = [0.5];
        let young = [12.0];
        let density = [3.0];
        // The unloaded transverse mechanism remains in the original pencil.
        let design = PreparedModalDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            ElasticImposedState::default(),
            ModalInertia {
                member_density_kg_m3: &density,
                joint_mass_kg: None,
            },
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        let (plan, options) = window_plan(vec![(2.0, 4.0)]);
        let no_bytes = ModalSensitivityControls {
            maximum_bytes: 1,
            ..Default::default()
        };
        assert!(matches!(
            design.solve_window(&plan, &options, no_bytes, cx),
            Err(TrussModalError::Elastic(ElasticError::Construction(
                TrussConstructionError::WorkBudget { .. }
            )))
        ));
        let (cluster, _) = window_plan(vec![(-1.0, 4.0)]);
        assert!(
            design
                .solve_window(&cluster, &options, ModalSensitivityControls::default(), cx)
                .is_err()
        );
        let (empty, _) = window_plan(vec![(4.0, 5.0)]);
        assert!(
            design
                .solve_window(&empty, &options, ModalSensitivityControls::default(), cx)
                .is_err()
        );
        let (multiple, _) = window_plan(vec![(1.0, 2.0), (2.0, 4.0)]);
        assert!(
            design
                .solve_window(&multiple, &options, ModalSensitivityControls::default(), cx)
                .is_err()
        );
        let compressed =
            LayoutCase::try_new(case.supported().to_vec(), vec![[0.0; 2], [-0.02, 0.0]], 2)
                .unwrap();
        let negative = PreparedModalDesign::try_new(
            &ground,
            &compressed,
            &areas,
            &young,
            ElasticImposedState::default(),
            ModalInertia {
                member_density_kg_m3: &density,
                joint_mass_kg: None,
            },
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        let (negative_window, _) = window_plan(vec![(-0.015, -0.005)]);
        let response = negative
            .solve_window(
                &negative_window,
                &options,
                ModalSensitivityControls::default(),
                cx,
            )
            .unwrap();
        near(response.eigenvalue_s2, -0.01);
        near(response.members[0].unwrap().area_derivative, 0.02);
        near(response.members[0].unwrap().density_derivative, 0.01 / 3.0);
        near(response.load_gradient[1][0], 0.5);
        // Fixed axial compression gives lambda proportional to L^-2,
        // preserving the derivative's sign on the negative spectral branch.
        near(response.node_coordinate_gradient[0][0], -0.01);
        near(response.node_coordinate_gradient[1][0], 0.01);
        gate.request();
        assert!(
            design
                .solve_window(&plan, &options, ModalSensitivityControls::default(), cx)
                .is_err()
        );
    });
}

// Lowest root of det(K-lambda M)=0 for the native two-coordinate pencil.
// This bypasses the spectral derivative path and retains off-diagonal mass.
fn lower_two_coordinate_eigenvalue(pencil: &TrussModalPencil) -> f64 {
    assert_eq!(pencil.mass().nrows(), 2);
    let k00 = entry(pencil.tangent(), 0, 0);
    let k01 = entry(pencil.tangent(), 0, 1);
    let k11 = entry(pencil.tangent(), 1, 1);
    let m00 = entry(pencil.mass(), 0, 0);
    let m01 = entry(pencil.mass(), 0, 1);
    let m11 = entry(pencil.mass(), 1, 1);
    let a = m00 * m11 - m01 * m01;
    let b = k00 * m11 + k11 * m00 - 2.0 * k01 * m01;
    let c = k00 * k11 - k01 * k01;
    2.0 * c / (b + (b * b - 4.0 * a * c).sqrt())
}

fn coordinate_ground(
    nodes: &[[f64; 2]],
    members: &[(usize, usize)],
    cx: &Cx<'_>,
) -> GroundStructure {
    let lengths: Vec<_> = members
        .iter()
        .map(|&(a, b)| (nodes[b][0] - nodes[a][0]).hypot(nodes[b][1] - nodes[a][1]))
        .collect();
    GroundStructure::try_from_parts(nodes, members, &lengths, GroundLimits::default(), cx).unwrap()
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_coordinate_adjoints_include_heated_settled_geometry_and_transform_covariantly() {
    use fs_truss::elastic::modal::sensitivity::{ModalSensitivityControls, PreparedModalDesign};
    with_cx(|cx, _| {
        let reference_nodes = [[-1.7, -0.3], [0.15, 0.22], [1.4, 0.4], [-0.35, -1.25]];
        let members = [(0, 1), (1, 2), (1, 3)];
        let areas = [1.1, 0.9, 1.3];
        let young = [2e6, 1e6, 0.7e6];
        let strain = [1e-4, 3e-5, -2e-5];
        let reference_motion = [[-2e-5, 1e-5], [0.0; 2], [3e-5, -1e-5], [-1e-5, 4e-5]];
        let density = [3.0, 6.0, 9.0];
        let joints = [0.2, 0.7, 0.3, 0.5];
        let mut reference_gradient: Option<Vec<[f64; 2]>> = None;
        let mut reference_eigenvalue = None;
        for angle in [0.0_f64, 0.43] {
            let (sn, cs) = angle.sin_cos();
            let rotate = |p: [f64; 2]| [cs * p[0] - sn * p[1], sn * p[0] + cs * p[1]];
            let nodes = reference_nodes.map(rotate);
            let motion = reference_motion.map(rotate);
            let ground = coordinate_ground(&nodes, &members, cx);
            let case = LayoutCase::try_new(
                vec![[true; 2], [false; 2], [true; 2], [true; 2]],
                vec![[0.0; 2], rotate([40.0, -25.0]), [0.0; 2], [0.0; 2]],
                nodes.len(),
            )
            .unwrap();
            let imposed = ElasticImposedState {
                member_stress_free_strain: Some(&strain),
                prescribed_displacements_m: Some(&motion),
            };
            let inertia = ModalInertia {
                member_density_kg_m3: &density,
                joint_mass_kg: Some(&joints),
            };
            let design = PreparedModalDesign::try_new(
                &ground,
                &case,
                &areas,
                &young,
                imposed,
                inertia,
                ModalAssemblyControls::default(),
                cx,
            )
            .unwrap();
            let expected = lower_two_coordinate_eigenvalue(design.pencil());
            let (plan, options) = window_plan(vec![(0.8 * expected, 1.2 * expected)]);
            let result = design
                .solve_window(&plan, &options, ModalSensitivityControls::default(), cx)
                .unwrap();
            near(result.eigenvalue_s2, expected);
            assert!(result.adjoint_iterations > 0);
            assert!(result.adjoint_relative_equilibrium_residual < 1e-9);
            assert!(
                design
                    .pencil()
                    .equilibrium()
                    .member_force_n
                    .iter()
                    .any(|n| n.abs() > 1.0)
            );
            let eigenvalue = |points: &[[f64; 2]]| {
                let geometry = coordinate_ground(points, &members, cx);
                let pencil = assemble(
                    &geometry,
                    &case,
                    &areas,
                    &young,
                    imposed,
                    inertia,
                    ModalAssemblyControls::default(),
                    cx,
                )
                .unwrap();
                lower_two_coordinate_eigenvalue(&pencil)
            };
            for node in 0..nodes.len() {
                for axis in 0..2 {
                    let h = 2e-5;
                    let mut plus = nodes;
                    let mut minus = nodes;
                    plus[node][axis] += h;
                    minus[node][axis] -= h;
                    derivative_near(
                        result.node_coordinate_gradient[node][axis],
                        (eigenvalue(&plus) - eigenvalue(&minus)) / (2.0 * h),
                    );
                }
            }
            // Reference position changes of strong supports are not masked.
            assert!(
                result.node_coordinate_gradient[0]
                    .iter()
                    .any(|v| v.abs() > 1.0)
            );
            // Rigid translation changes no member and gives a zero covector sum.
            for axis in 0..2 {
                let sum: f64 = result
                    .node_coordinate_gradient
                    .iter()
                    .map(|g| g[axis])
                    .sum();
                assert!(sum.abs() < 1e-10 * expected);
            }
            let shifted = nodes.map(|p| [p[0] + 0.73, p[1] - 1.31]);
            near(eigenvalue(&shifted), expected);
            // Scaling reference coordinates AND prescribed displacements while
            // scaling point masses like member mass gives lambda(s)=lambda/s².
            let mut homothety = 0.0;
            let mut rotation = 0.0;
            for node in 0..nodes.len() {
                for axis in 0..2 {
                    homothety += result.node_coordinate_gradient[node][axis] * nodes[node][axis]
                        + result.prescribed_displacement_gradient[node][axis] * motion[node][axis];
                }
                homothety += result.joint_mass_derivative[node] * joints[node];
                let turn = |v: [f64; 2]| [-v[1], v[0]];
                let dx = turn(nodes[node]);
                let df = turn(case.loads()[node]);
                let du = turn(motion[node]);
                for axis in 0..2 {
                    rotation += result.node_coordinate_gradient[node][axis] * dx[axis]
                        + result.load_gradient[node][axis] * df[axis]
                        + result.prescribed_displacement_gradient[node][axis] * du[axis];
                }
            }
            near(homothety, -2.0 * expected);
            assert!(rotation.abs() < 2e-8 * expected);
            if let Some(reference) = &reference_gradient {
                for (g, original) in result.node_coordinate_gradient.iter().zip(reference) {
                    let rotated = rotate(*original);
                    derivative_near(g[0], rotated[0]);
                    derivative_near(g[1], rotated[1]);
                }
                near(expected, reference_eigenvalue.unwrap());
            } else {
                reference_gradient = Some(result.node_coordinate_gradient.clone());
                reference_eigenvalue = Some(expected);
            }
        }
    });
}

#[test]
fn g1_coordinate_adjoints_retain_the_serial_consistent_mass_coupling() {
    use fs_truss::elastic::modal::sensitivity::{ModalSensitivityControls, PreparedModalDesign};
    with_cx(|cx, _| {
        let nodes = [[0.0, 0.0], [0.8, 0.0], [2.1, 0.0]];
        let members = [(0, 1), (1, 2)];
        let ground = coordinate_ground(&nodes, &members, cx);
        let case = LayoutCase::try_new(
            vec![[true; 2], [false, true], [false, true]],
            vec![[0.0; 2]; 3],
            3,
        )
        .unwrap();
        let areas = [0.5, 0.8];
        let young = [12.0, 19.0];
        let density = [3.0, 5.0];
        let inertia = ModalInertia {
            member_density_kg_m3: &density,
            joint_mass_kg: None,
        };
        let design = PreparedModalDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            ElasticImposedState::default(),
            inertia,
            ModalAssemblyControls::default(),
            cx,
        )
        .unwrap();
        assert!(entry(design.pencil().mass(), 0, 1) > 0.0);
        let expected = lower_two_coordinate_eigenvalue(design.pencil());
        let (plan, options) = window_plan(vec![(0.8 * expected, 1.2 * expected)]);
        let result = design
            .solve_window(&plan, &options, ModalSensitivityControls::default(), cx)
            .unwrap();
        near(result.eigenvalue_s2, expected);
        assert_eq!(result.adjoint_iterations, 0);
        for node in 0..nodes.len() {
            let mut values = [0.0; 2];
            let h = 1e-5;
            for (value, sign) in values.iter_mut().zip([1.0, -1.0]) {
                let mut changed = nodes;
                changed[node][0] += sign * h;
                let geometry = coordinate_ground(&changed, &members, cx);
                let pencil = assemble(
                    &geometry,
                    &case,
                    &areas,
                    &young,
                    ElasticImposedState::default(),
                    inertia,
                    ModalAssemblyControls::default(),
                    cx,
                )
                .unwrap();
                *value = lower_two_coordinate_eigenvalue(&pencil);
            }
            derivative_near(
                result.node_coordinate_gradient[node][0],
                (values[0] - values[1]) / (2.0 * h),
            );
            assert_eq!(result.node_coordinate_gradient[node][1], 0.0);
        }
        let derivative: f64 = result
            .node_coordinate_gradient
            .iter()
            .zip(nodes)
            .map(|(g, x)| g[0] * x[0])
            .sum();
        near(derivative, -2.0 * expected);
    });
}

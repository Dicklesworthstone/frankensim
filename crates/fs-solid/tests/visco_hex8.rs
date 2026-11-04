//! Native full-tensor HEX8 histories, B-bar work and consistent-mass motion.
//! All material numbers below are explicitly synthetic software fixtures.

use fs_exec::Cx;
use fs_material::visco::GeneralizedMaxwell;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseControls;
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticError, TetElasticMaterial, TetMaterialField, TetStaticSolveConfig,
};

const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const CELLS: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const FORMS: [Formulation; 2] = [Formulation::Standard, Formulation::BBar];

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 814,
                    kernel_id: 814,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    near(
        actual,
        expected,
        tolerance * actual.abs().max(expected.abs()).max(1e-20),
    );
}
fn diagonal(d: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { d[i] } else { 0.0 }))
}
fn outer(v: [f64; 6], scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| scale * v[i] * v[j]))
}
fn isotropic(k: f64, g: f64) -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i < 3 && j < 3 {
                k + if i == j {
                    4.0 * g / 3.0
                } else {
                    -2.0 * g / 3.0
                }
            } else if i == j {
                2.0 * g
            } else {
                0.0
            }
        })
    })
}
fn material(law: TensorProny) -> TetViscoMaterial {
    TetViscoMaterial {
        law,
        density_kg_m3: 27.0,
        material_state_identity: [81; 32],
    }
}
fn coupled() -> TetViscoMaterial {
    let mut c = diagonal([90.0, 70.0, 50.0, 60.0, 30.0, 40.0]);
    c[0][3] = -8.0;
    c[3][0] = -8.0;
    c[1][3] = 5.0;
    c[3][1] = 5.0;
    material(
        TensorProny::new_mandel(
            c,
            vec![
                PronyBranch::new(outer([2.0, 0.0, 1.0, 3.0, 0.0, -1.0], 1.0), 0.03).unwrap(),
                PronyBranch::new(outer([1.0, -1.0, 0.0, -1.0, 2.0, 3.0], 0.5), 0.3).unwrap(),
            ],
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(outer([1.0, 0.0, 0.0, -2.0, 0.0, 1.0], 0.2))
        .unwrap(),
    )
}
fn mesh(nodes: &[[f64; 3]], formulation: Formulation, order: usize) -> SolidMeshRef<'_> {
    SolidMeshRef::Hex8 {
        nodes_m: nodes,
        elements: &CELLS,
        formulation,
        stiffness_points_per_dir: order,
    }
}
fn system(
    nodes: &[[f64; 3]],
    material: &TetViscoMaterial,
    formulation: Formulation,
    fixed: &[usize],
    order: usize,
    cx: &Cx<'_>,
) -> ViscoElements {
    ViscoElements::new(
        mesh(nodes, formulation, order),
        TetViscoField::Uniform(material),
        fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}
fn config() -> ViscoDynamicConfig {
    ViscoDynamicConfig {
        maximum_steps: 200,
        maximum_free_dofs: 100,
        solve: TetStaticSolveConfig {
            maximum_iterations: 300,
            relative_residual_tolerance: 2e-11,
        },
        absolute_energy_tolerance_j: 1e-12,
        relative_energy_tolerance: 2e-9,
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn warped_history_and_mass_use_the_original_elastic_pencil() {
    with_cx(|cx, _| {
        let nodes = NODES.map(|[x, y, z]| [x + 0.1 * y * z, y + 0.07 * x * z, z + 0.05 * x * y]);
        let m = material(TensorProny::new_mandel(*coupled().law.c_eq(), vec![], 0.02).unwrap());
        let elastic = TetElasticMaterial::try_new_mandel(
            m.density_kg_m3,
            *m.law.c_eq(),
            m.material_state_identity,
        )
        .unwrap();
        for formulation in FORMS {
            for order in [2, 3] {
                let solid = system(&nodes, &m, formulation, &[], order, cx);
                let pencil = Hex8ElasticProblem {
                    nodes_m: &nodes,
                    hexahedra: &CELLS,
                    materials: TetMaterialField::Uniform(&elastic),
                    formulation,
                    fixed_dofs: &[],
                    budget: TetAssemblyBudget::standard(),
                    stiffness_points_per_dir: order,
                }
                .assemble(cx)
                .unwrap();
                for i in 0..24 {
                    for j in 0..24 {
                        assert_eq!(
                            solid.equilibrium().stiffness.get(i, j),
                            pencil.stiffness.get(i, j)
                        );
                        assert_eq!(solid.equilibrium().mass.get(i, j), pencil.mass.get(i, j));
                    }
                }
                let retained = solid.reference_material(0).unwrap();
                assert_eq!(retained.material_state_identity, m.material_state_identity);
                assert_eq!(retained.density_kg_m3, m.density_kg_m3);
                let u: Vec<_> = (0..24)
                    .map(|i| 0.0001 * ((7 * i % 11) as f64 - 5.0))
                    .collect();
                let mut expected = vec![0.0; 24];
                pencil.stiffness.spmv(&u, &mut expected);
                let response = solid
                    .state(cx)
                    .unwrap()
                    .step_internal_force(&u, 0.01, cx)
                    .unwrap();
                for (actual, expected) in response.full_force_n.iter().zip(&expected) {
                    near(*actual, *expected, 2e-13);
                }
                relative(response.energy.stored_j, 0.5 * dot(&u, &expected), 2e-12);
                assert_eq!(response.step_dissipated_j, 0.0);
            }
        }
    });
}

#[test]
fn zero_mean_dilatation_does_not_activate_bulk_relaxation_in_bbar() {
    with_cx(|cx, _| {
        let g = 23.0;
        let a = 0.001;
        // u_x=a*x*(2y-1): mean dilatation is zero. Its integrated deviatoric
        // strain norm squared is 8a^2/9 on the unit cube, independently of K.
        let u: Vec<_> = NODES
            .iter()
            .flat_map(|[x, y, _]| [a * x * (2.0 * y - 1.0), 0.0, 0.0])
            .collect();
        let expected = 8.0 * g * a * a / 9.0;
        for k in [100.0 * g, 1e6 * g] {
            let m = material(
                TensorProny::isotropic(
                    GeneralizedMaxwell::new(k, vec![(2.0 * k, 0.05)]).unwrap(),
                    GeneralizedMaxwell::new(g, vec![]).unwrap(),
                )
                .unwrap(),
            );
            let bbar = system(&NODES, &m, Formulation::BBar, &[], 2, cx);
            let response = bbar
                .state(cx)
                .unwrap()
                .step_internal_force(&u, 0.03, cx)
                .unwrap();
            relative(response.energy.stored_j, expected, 2e-10);
            relative(response.step_work_j, expected, 2e-10);
            assert!(response.step_dissipated_j < 1e-18 * expected);
            relative(dot(&u, &response.full_force_n), 2.0 * expected, 2e-9);
            let standard = system(&NODES, &m, Formulation::Standard, &[], 2, cx);
            let unprojected = standard
                .state(cx)
                .unwrap()
                .step_internal_force(&u, 0.03, cx)
                .unwrap();
            assert!(unprojected.energy.stored_j > expected + k * a * a / 6.0);
            assert!(unprojected.step_dissipated_j > 0.0);
        }
    });
}

// Independent exact strain for u_x=a*x*y on a unit cube; the source law is
// evaluated at the two-point product rule, with no native shape/B-bar calls.
fn xy_strain(formulation: Formulation, x: f64, y: f64, a: f64) -> [f64; 6] {
    let normals = if formulation == Formulation::BBar {
        [
            a * (1.0 / 6.0 + 2.0 * y / 3.0),
            a * (1.0 / 6.0 - y / 3.0),
            a * (1.0 / 6.0 - y / 3.0),
        ]
    } else {
        [a * y, 0.0, 0.0]
    };
    [
        normals[0],
        normals[1],
        normals[2],
        a * x * core::f64::consts::FRAC_1_SQRT_2,
        0.0,
        0.0,
    ]
}

#[test]
fn anisotropic_memory_and_kelvin_work_match_independent_point_histories() {
    with_cx(|cx, _| {
        let m = coupled();
        let q = [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()];
        let points: Vec<_> = q
            .iter()
            .flat_map(|&x| q.iter().map(move |&y| (x, y)))
            .collect();
        let mode: Vec<_> = NODES
            .iter()
            .flat_map(|[x, y, _]| [x * y, 0.0, 0.0])
            .collect();
        for formulation in FORMS {
            let solid = system(&NODES, &m, formulation, &[], 2, cx);
            assert_eq!(solid.internal_variables(), 8 * 6 * 2);
            let mut native = solid.state(cx).unwrap();
            let mut oracle: Vec<_> = points.iter().map(|_| m.law.state()).collect();
            let mut old_a = 0.0;
            for (a, dt) in [(0.001, 0.02), (0.0015, 0.1), (-0.0003, 0.05)] {
                let u: Vec<_> = mode.iter().map(|v| a * v).collect();
                let actual = native.step_internal_force(&u, dt, cx).unwrap();
                let (mut stored, mut work, mut heat, mut force, mut mean) =
                    (0.0, 0.0, 0.0, 0.0, 0.0);
                for (state, &(x, y)) in oracle.iter_mut().zip(&points) {
                    let response = m
                        .law
                        .step_with_mean(state, xy_strain(formulation, x, y, a), dt)
                        .unwrap();
                    stored += 0.25 * state.stored;
                    work += 0.25 * state.work;
                    heat += 0.25 * state.dissipated;
                    let b = xy_strain(formulation, x, y, 1.0);
                    force += 0.25 * dot(&b, &response.stress);
                    mean += 0.25 * dot(&b, &response.mean_stress);
                }
                relative(actual.energy.stored_j, stored, 3e-12);
                relative(actual.energy.work_j, work, 3e-12);
                relative(actual.energy.dissipated_j, heat, 3e-12);
                relative(dot(&mode, &actual.full_force_n), force, 3e-12);
                relative(dot(&mode, &actual.full_mean_force_n), mean, 3e-12);
                relative(actual.step_work_j, (a - old_a) * mean, 3e-12);
                relative(
                    actual.energy.work_j,
                    actual.energy.stored_j + actual.energy.dissipated_j,
                    3e-12,
                );
                old_a = a;
            }
        }
    });
}

#[test]
fn near_incompressible_kelvin_brick_dynamics_converges_to_analytic_ringdown() {
    with_cx(|cx, _| {
        let (k, g, bulk_viscosity, shear_viscosity) = (1e4, 1.0, 10.0, 0.3);
        let m = material(
            TensorProny::new_mandel(isotropic(k, g), vec![], 0.02)
                .unwrap()
                .with_kelvin_viscosity_mandel(isotropic(bulk_viscosity, shear_viscosity))
                .unwrap(),
        );
        let fixed: Vec<_> = (0..24).filter(|&d| d != 3).collect();
        let solid = system(&NODES, &m, Formulation::BBar, &fixed, 2, cx);
        // Only u_1x is free: integral N_1^2=1/27, mean N_1,x=1/4,
        // and integral ||dev(sym(e_x tensor grad N_1))||^2=5/27.
        let mass = m.density_kg_m3 / 27.0;
        let stiffness = k / 16.0 + 10.0 * g / 27.0;
        let viscosity = bulk_viscosity / 16.0 + 10.0 * shear_viscosity / 27.0;
        relative(solid.equilibrium().mass.get(0, 0), mass, 2e-13);
        relative(solid.equilibrium().stiffness.get(0, 0), stiffness, 2e-12);
        let (u0, v0, end) = (0.0001, -0.001, 0.2);
        let alpha = viscosity / (2.0 * mass);
        let omega = (stiffness / mass - alpha * alpha).sqrt();
        let b = (v0 + alpha * u0) / omega;
        let phase = omega * end;
        let shape = u0 * phase.cos() + b * phase.sin();
        let exact_u = (-alpha * end).exp() * shape;
        let exact_v = (-alpha * end).exp()
            * (-alpha * shape - u0 * omega * phase.sin() + b * omega * phase.cos());
        let mut errors = Vec::new();
        for count in [16, 32, 64] {
            let dt = end / f64::from(count);
            let prepared = solid.dynamic_operator(dt, config(), cx).unwrap();
            let mut initial = solid.state(cx).unwrap();
            let mut u = [0.0; 24];
            u[3] = u0;
            initial.step_internal_force(&u, 0.1, cx).unwrap();
            let previous_heat = initial.energy().dissipated_j;
            let mut v = [0.0; 24];
            v[3] = v0;
            let mut state = prepared.start(initial, &v, cx).unwrap();
            let energy = state.kinetic_energy_j() + state.material_history().energy().stored_j;
            let (mut heat, mut previous_u) = (0.0, u0);
            for _ in 0..count {
                let step = state.advance(&[0.0; 24], cx).unwrap();
                let next = state.material_history().displacement_m()[3];
                heat += viscosity * (next - previous_u).powi(2) / dt;
                relative(
                    step.constitutive.energy.dissipated_j - previous_heat,
                    heat,
                    2e-9,
                );
                relative(
                    step.kinetic_energy_j + step.constitutive.energy.stored_j + heat,
                    energy,
                    2e-9,
                );
                previous_u = next;
            }
            errors.push(
                (state.material_history().displacement_m()[3] - exact_u)
                    .hypot((state.velocity_mps()[3] - exact_v) / omega),
            );
        }
        assert!(errors[1] < 0.28 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.28 * errors[1], "{errors:?}");
        assert!(
            errors[2] < 0.003 * exact_u.hypot(exact_v / omega),
            "{errors:?}"
        );
    });
}

#[test]
fn anisotropic_quasistatic_preload_releases_through_the_same_dynamic_history() {
    with_cx(|cx, gate| {
        let m = coupled();
        let fixed: Vec<_> = NODES
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0] == 0.0)
            .flat_map(|(n, _)| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect();
        let solid = system(&NODES, &m, Formulation::BBar, &fixed, 2, cx);
        let u: Vec<_> = NODES
            .iter()
            .flat_map(|[x, y, z]| [0.001 * x * (1.0 + 0.2 * y), -0.0003 * x * z, 0.0002 * x * y])
            .collect();
        let forces = solid
            .state(cx)
            .unwrap()
            .step_internal_force(&u, 0.04, cx)
            .unwrap()
            .full_force_n;
        let load_operator = solid.quasistatic_operator(0.04, cx).unwrap();
        let mut preload = solid.state(cx).unwrap();
        let loaded = preload
            .step_quasistatic(&load_operator, &forces, &[0.0; 24], config().solve, cx)
            .unwrap();
        for (&actual, &expected) in loaded.displacement_m.iter().zip(&u) {
            near(actual, expected, 2e-12);
        }
        let previous_heat = preload.energy().dissipated_j;
        let energy = preload.energy().stored_j;
        let prepared = solid.dynamic_operator(0.01, config(), cx).unwrap();
        let mut state = prepared.start(preload, &[0.0; 24], cx).unwrap();
        for _ in 0..40 {
            let step = state.advance(&[0.0; 24], cx).unwrap();
            assert!(step.constitutive.step_dissipated_j >= 0.0);
            relative(
                step.kinetic_energy_j
                    + step.constitutive.energy.stored_j
                    + step.constitutive.energy.dissipated_j
                    - previous_heat,
                energy,
                2e-9,
            );
        }
        assert!(state.material_history().energy().dissipated_j > previous_heat);
        let before = state.material_history().energy();
        let displacement = state.material_history().displacement_m().to_vec();
        let velocity = state.velocity_mps().to_vec();
        let steps = state.steps();
        assert!(matches!(
            state.advance(&[1e10; 24], cx),
            Err(ViscoElementError::Material(_))
        ));
        assert_eq!(state.material_history().energy(), before);
        assert_eq!(state.material_history().displacement_m(), displacement);
        assert_eq!(state.velocity_mps(), velocity);
        assert_eq!(state.steps(), steps);
        gate.request();
        assert!(matches!(
            state.advance(&[0.0; 24], cx),
            Err(ViscoElementError::Elastic(TetElasticError::Cancelled))
        ));
        assert_eq!(state.material_history().energy(), before);
        assert_eq!(state.steps(), steps);
    });
}

#[test]
fn compatible_strain_is_checked_in_time_and_over_the_full_harmonic_cycle() {
    with_cx(|cx, _| {
        let amplitude = 0.01;
        let m = material(
            TensorProny::new_mandel(diagonal([50.0; 6]), vec![], 1.23 * amplitude)
                .unwrap()
                .with_kelvin_viscosity_mandel(diagonal([0.2; 6]))
                .unwrap(),
        );
        let solid = system(&NODES, &m, Formulation::BBar, &[], 2, cx);
        let bad: Vec<_> = NODES
            .iter()
            .flat_map(|[x, y, _]| [amplitude * x * (2.0 * y - 1.0), 0.0, 0.0])
            .collect();
        let q = [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()];
        // At the largest Gauss point the projected norm is 1.211a, while
        // the compatible norm is 1.256a. The original limit lies between them.
        for x in q {
            for y in q {
                let d = amplitude * (2.0 * y - 1.0);
                m.law
                    .validate_mandel_strain(&[
                        2.0 * d / 3.0,
                        -d / 3.0,
                        -d / 3.0,
                        core::f64::consts::SQRT_2 * amplitude * x,
                        0.0,
                        0.0,
                    ])
                    .unwrap();
            }
        }
        assert!(matches!(
            solid.validate_displacement(&bad, || false),
            Err(ViscoElementError::Material(_))
        ));
        let mut state = solid.state(cx).unwrap();
        let good: Vec<_> = bad.iter().map(|v| 0.1 * v).collect();
        state.step_internal_force(&good, 0.1, cx).unwrap();
        let before = state.energy();
        let mut untouched = state.clone();
        assert!(matches!(
            state.step_internal_force(&bad, 0.1, cx),
            Err(ViscoElementError::Material(_))
        ));
        assert_eq!(state.energy(), before);
        assert_eq!(state.displacement_m(), good);
        let actual = state.step_internal_force(&[0.0; 24], 0.1, cx).unwrap();
        let expected = untouched.step_internal_force(&[0.0; 24], 0.1, cx).unwrap();
        assert_eq!(actual.full_force_n, expected.full_force_n);
        assert_eq!(actual.energy, expected.energy);

        let omega = 0.37;
        let complex = solid.complex_stiffness(omega, cx).unwrap();
        let force: Vec<_> = (0..24)
            .map(|i| {
                C64::new(
                    (0..24)
                        .map(|j| {
                            (complex.real.get(i, j)
                                - omega * omega * solid.equilibrium().mass.get(i, j))
                                * bad[j]
                        })
                        .sum(),
                    -(0..24)
                        .map(|j| complex.imag.get(i, j) * bad[j])
                        .sum::<f64>(),
                )
            })
            .collect();
        let harmonic = solid
            .prepare_harmonic_response(
                omega,
                HarmonicResponseControls {
                    maximum_coordinates: 32,
                    maximum_payload_bytes: 1024 * 1024,
                    maximum_relative_residual: 2e-10,
                },
                cx,
            )
            .unwrap();
        assert!(matches!(
            harmonic.solve(&force, cx),
            Err(ViscoHarmonicError::Element(ViscoElementError::Material(_)))
        ));
    });
}

#[test]
fn hex8_admission_counts_every_memory_point_and_bounds_reduced_projection() {
    with_cx(|cx, _| {
        let m = coupled();
        let create = |order, cap| {
            ViscoElements::new(
                mesh(&NODES, Formulation::BBar, order),
                TetViscoField::Uniform(&m),
                &[],
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: cap,
                },
                cx,
            )
        };
        assert!(matches!(
            create(2, 95),
            Err(ViscoElementError::InternalStateBudget {
                needed: 96,
                cap: 95
            })
        ));
        assert!(matches!(
            create(3, 323),
            Err(ViscoElementError::InternalStateBudget {
                needed: 324,
                cap: 323
            })
        ));
        assert!(matches!(
            create(1, 1000),
            Err(ViscoElementError::Elastic(
                TetElasticError::UnderIntegrated { .. }
            ))
        ));
        assert!(matches!(
            create(usize::MAX, 1000),
            Err(ViscoElementError::Elastic(
                TetElasticError::QuadratureInadmissible { .. }
            ))
        ));
        let budget = SolidProjectionBudget {
            max_coordinates: 1,
            max_basis_scalars: 24,
            max_material_points: 8,
            max_scalar_terms: 100_000,
        };
        let bbar = create(2, 96).unwrap();
        let projected = bbar.project_basis(&[vec![1.0; 24]], budget, cx).unwrap();
        assert!(projected.material_points().iter()
            .all(|point| point.compatible_strain_columns().is_some()));
        let standard = system(&NODES, &m, Formulation::Standard, &[], 2, cx);
        let projected = standard.project_basis(&[vec![1.0; 24]], budget, cx).unwrap();
        assert!(projected.material_points().iter()
            .all(|point| point.compatible_strain_columns().is_none()));
    });
}

#[test]
fn hex8_bbar_projection_preserves_original_work_and_local_strain_admission() {
    with_cx(|cx, gate| {
        let m = coupled();
        let solid = system(&NODES, &m, Formulation::BBar, &[], 2, cx);
        let columns: Vec<Vec<f64>> = (0..2)
            .map(|axis| {
                NODES
                    .iter()
                    .flat_map(|p| {
                        let mut u = [0.0; 3];
                        u[axis] = (p[axis] - 0.5) * (p[axis + 1] - 0.5);
                        u
                    })
                    .collect()
            })
            .collect();
        let budget = SolidProjectionBudget {
            max_coordinates: 2,
            max_basis_scalars: 48,
            max_material_points: 8,
            max_scalar_terms: 100_000,
        };
        let projected = solid.project_basis(&columns, budget, cx).unwrap();
        let native = solid.equilibrium();
        let mut applied = vec![0.0; 24];
        for j in 0..2 {
            native.stiffness.spmv(&columns[j], &mut applied);
            for i in 0..2 {
                let expected: f64 = columns[i].iter().zip(&applied).map(|(a, b)| a * b).sum();
                let actual: f64 = projected
                    .material_points()
                    .iter()
                    .map(|point| point.equilibrium_stiffness()[2 * i + j])
                    .sum();
                near(actual, expected, 1e-12);
            }
            native.mass.spmv(&columns[j], &mut applied);
            for i in 0..2 {
                let expected: f64 = columns[i].iter().zip(&applied).map(|(a, b)| a * b).sum();
                near(projected.mass_matrix_kg()[2 * i + j], expected, 1e-12);
            }
        }

        let norm = |strain: &[f64; 6]| {
            (strain[..3].iter().map(|x| x * x).sum::<f64>()
                + 2.0 * strain[3..].iter().map(|x| x * x).sum::<f64>())
            .sqrt()
        };
        let assumed = projected
            .material_points()
            .iter()
            .map(|point| norm(&point.strain_columns()[0]))
            .fold(0.0, f64::max);
        let compatible = projected
            .material_points()
            .iter()
            .map(|point| norm(&point.compatible_strain_columns().unwrap()[0]))
            .fold(0.0, f64::max);
        assert!(compatible > assumed);
        let amplitude = m.law.maximum_strain_norm() / (assumed * compatible).sqrt();
        for point in projected.material_points() {
            let unguarded = m
                .law
                .project_strain(
                    point.strain_columns(),
                    point.volume_m3(),
                    fs_material::visco_tensor::TensorProjectionBudget {
                        max_coordinates: 2,
                        max_scalar_terms: 100_000,
                    },
                    || false,
                )
                .unwrap();
            assert!(unguarded.strain(&[amplitude, 0.0]).is_ok());
        }
        assert!(
            projected
                .material_points()
                .iter()
                .any(|point| point.strain(&[amplitude, 0.0]).is_err())
        );
        let rebuilt = projected
            .with_element_laws(&[m.law.clone()], budget, gate)
            .unwrap();
        for (old, new) in projected
            .material_points()
            .iter()
            .zip(rebuilt.material_points())
        {
            assert_eq!(
                old.compatible_strain_columns(),
                new.compatible_strain_columns()
            );
            assert_eq!(
                old.strain(&[amplitude, 0.0]).is_ok(),
                new.strain(&[amplitude, 0.0]).is_ok()
            );
        }
    });
}

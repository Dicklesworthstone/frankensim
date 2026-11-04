//! Independent P2 mass, conservative loading and constrained-mode identities.

use std::cell::Cell;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::HyperTetError;
use fs_solid::hyper3::gravity::{GravityField3, GravitySample3};
use fs_solid::hyper3::mixed::pressure::UniformBoundaryPressure;
use fs_solid::hyper3::mixed::vibration::MixedVibrationSettings;
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpError, MixedUpProblem, MixedUpSettings};
use fs_solid::hyper3::vibration::VibrationSettings;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

fn with_gate(gate: &CancelGate, f: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 61,
                kernel_id: 137,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    with_gate(&CancelGate::new(), f);
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn law(bulk: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 1.0,
            c01: 0.0,
            kappa: bulk,
        },
        5.0,
    )
    .unwrap()
}

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

fn clamp_face(nodes: &[[f64; 3]]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, x)| x[2] == 0.0)
        .flat_map(|(node, _)| (0..3).map(move |i| (3 * node + i, 0.0)))
        .collect()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a Hyperelastic,
    supports: &'a [(usize, f64)],
    forces: &'a [(usize, f64)],
    volumetric: VolumetricTreatment,
) -> MixedUpProblem<'a> {
    MixedUpProblem {
        nodes_m: nodes,
        elements,
        material: MixedLawField::Uniform(material),
        volumetric,
        prescribed_m: supports,
        nodal_forces_n: forces,
        body_force_n_m3: [0.0; 3],
        pressure_gauge: None,
        budget: TetAssemblyBudget::standard(),
        quadrature_points_per_dir: 4,
        maximum_quadrature_points: 20_000,
        settings: MixedUpSettings {
            load_steps: 1,
            maximum_newton_iterations: 30,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-11,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            length_scale_m: 1.0,
            stress_scale_pa: 2.0,
            maximum_assembly_bytes: 8 * 1024 * 1024,
            maximum_factor_bytes: 8 * 1024 * 1024,
        },
    }
}

fn modal_settings() -> MixedVibrationSettings {
    MixedVibrationSettings {
        spectrum: VibrationSettings {
            maximum_free_dofs: 128,
            relative_residual_tolerance: 1e-8,
        },
        maximum_mass_assembly_bytes: 8 * 1024 * 1024,
        maximum_condensation_bytes: 16 * 1024 * 1024,
        pressure_residual_tolerance: 1e-10,
    }
}

// Independent rational barycentric integrals on the unit tetrahedron:
// integral lambda^alpha = product(alpha!) / (3+sum(alpha))!.
// N_i=lambda_i(2lambda_i-1), N_ij=4lambda_i lambda_j. In particular,
// corner row sums are NEGATIVE, -rho V/20; edge row sums are rho V/5.
fn rational_mass(tet: &[usize; 10], density: f64) -> [[f64; 10]; 10] {
    let mut mass = [[0.0; 10]; 10];
    for a in 0..10 {
        for b in 0..10 {
            let coefficient = match (a < 4, b < 4) {
                (true, true) => {
                    if a == b {
                        1.0 / 70.0
                    } else {
                        1.0 / 420.0
                    }
                }
                (true, false) => {
                    if TET10_EDGES[b - 4].contains(&a) {
                        -1.0 / 105.0
                    } else {
                        -1.0 / 70.0
                    }
                }
                (false, true) => {
                    if TET10_EDGES[a - 4].contains(&b) {
                        -1.0 / 105.0
                    } else {
                        -1.0 / 70.0
                    }
                }
                (false, false) => {
                    if a == b {
                        8.0 / 105.0
                    } else if TET10_EDGES[a - 4]
                        .iter()
                        .any(|i| TET10_EDGES[b - 4].contains(i))
                    {
                        4.0 / 105.0
                    } else {
                        2.0 / 105.0
                    }
                }
            };
            mass[tet[a]][tet[b]] = density * coefficient / 6.0;
        }
    }
    mass
}

fn affine_field(x: [f64; 3]) -> Result<GravitySample3, HyperTetError> {
    let g0 = [0.3, -0.5, -1.1];
    let gradient = [[0.2, 0.03, -0.01], [0.03, -0.15, 0.02], [-0.01, 0.02, 0.1]];
    let acceleration =
        core::array::from_fn(|i| g0[i] + (0..3).map(|j| gradient[i][j] * x[j]).sum::<f64>());
    let potential = 4.0
        - (0..3).map(|i| g0[i] * x[i]).sum::<f64>()
        - 0.5
            * (0..3)
                .map(|i| (0..3).map(|j| x[i] * gradient[i][j] * x[j]).sum::<f64>())
                .sum::<f64>();
    Ok(GravitySample3 {
        potential_j_kg: potential,
        acceleration_m_s2: acceleration,
        acceleration_gradient_s2: gradient,
    })
}

#[test]
fn current_quadratic_gravity_matches_rational_mass_force_tangent_and_potential() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(20.0);
        let supports = clamp_face(&nodes);
        let p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        let density = [7.0];
        let mass = rational_mass(&elements[0], density[0]);
        let gravity = p
            .prepare_with_gravity(&density, &affine_field, None, cx)
            .unwrap();
        let u: Vec<_> = nodes
            .iter()
            .map(|x| {
                [
                    0.03 * x[0] + 0.015 * x[0] * x[2],
                    -0.01 * x[1] * x[1],
                    0.02 * x[0] * x[1],
                ]
            })
            .collect();
        let pressure = [0.2, -0.1, 0.3, -0.2];
        let plain = p
            .prepare(cx)
            .unwrap()
            .evaluate(&u, &pressure, true, cx)
            .unwrap();
        let loaded = gravity.evaluate(&u, &pressure, true, cx).unwrap();
        let positions: Vec<[f64; 3]> = nodes
            .iter()
            .zip(&u)
            .map(|(x, u)| core::array::from_fn(|i| x[i] + u[i]))
            .collect();
        let samples: Vec<_> = positions
            .iter()
            .map(|&x| affine_field(x).unwrap())
            .collect();
        let gradient = samples[0].acceleration_gradient_s2;
        let k0 = plain.tangent.as_ref().unwrap();
        let k = loaded.tangent.as_ref().unwrap();
        for a in 0..10 {
            for i in 0..3 {
                let expected_force: f64 = (0..10)
                    .map(|b| mass[a][b] * samples[b].acceleration_m_s2[i])
                    .sum();
                near(
                    plain.internal_forces_n[3 * a + i] - loaded.internal_forces_n[3 * a + i],
                    expected_force,
                    2e-13,
                );
                for b in 0..10 {
                    for j in 0..3 {
                        near(
                            k.get(3 * a + i, 3 * b + j) - k0.get(3 * a + i, 3 * b + j),
                            -mass[a][b] * gradient[i][j],
                            3e-12,
                        );
                    }
                }
                let mut plus = u.clone();
                let mut minus = u.clone();
                plus[a][i] += 1e-5;
                minus[a][i] -= 1e-5;
                let derivative = (gravity
                    .loading_potentials(&plus, cx)
                    .unwrap()
                    .gravitational_j
                    - gravity
                        .loading_potentials(&minus, cx)
                        .unwrap()
                        .gravitational_j)
                    / 2e-5;
                near(derivative, -expected_force, 2e-9);
            }
        }
        let g0 = affine_field([0.0; 3]).unwrap().acceleration_m_s2;
        let mut expected_potential = 4.0 * density[0] / 6.0;
        for a in 0..10 {
            for b in 0..10 {
                for i in 0..3 {
                    expected_potential -= mass[a][b] * g0[i] * positions[b][i];
                    for j in 0..3 {
                        expected_potential -=
                            0.5 * mass[a][b] * positions[a][i] * gradient[i][j] * positions[b][j];
                    }
                }
            }
        }
        near(
            gravity.loading_potentials(&u, cx).unwrap().gravitational_j,
            expected_potential,
            3e-13,
        );
        assert_eq!(loaded.potential_j, plain.potential_j);
        assert_eq!(loaded.pressure_residual_m3, plain.pressure_residual_m3);
        for row in 30..34 {
            assert_eq!(k.row(row), k0.row(row));
            for col in 0..30 {
                assert_eq!(k.get(col, row), k0.get(col, row));
            }
        }
        near(mass[0].iter().sum(), -density[0] / 120.0, 2e-17);
    });
}

fn harmonic_field(x: [f64; 3]) -> Result<GravitySample3, HyperTetError> {
    let gamma = 1.7;
    Ok(GravitySample3 {
        potential_j_kg: 0.5 * gamma * x.iter().map(|v| v * v).sum::<f64>(),
        acceleration_m_s2: x.map(|v| -gamma * v),
        acceleration_gradient_s2: [[-gamma, 0.0, 0.0], [0.0, -gamma, 0.0], [0.0, 0.0, -gamma]],
    })
}

fn harmonic_balance(nodes: &[[f64; 3]], tet: &[usize; 10], density: f64) -> Vec<(usize, f64)> {
    let mass = rational_mass(tet, density);
    (0..30)
        .map(|d| {
            (
                d,
                1.7 * (0..10)
                    .map(|b| mass[d / 3][b] * nodes[b][d % 3])
                    .sum::<f64>(),
            )
        })
        .collect()
}

#[test]
fn finite_compliance_modes_include_the_exact_reference_mass_gravity_shift() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let supports = clamp_face(&nodes);
        let material = law(20.0);
        let density = [7.0];
        let forces = harmonic_balance(&nodes, &elements[0], density[0]);
        let field: &GravityField3<'_> = &harmonic_field;
        let p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        let base = p.solve_vibration(&density, modal_settings(), cx).unwrap();
        let loaded_problem = MixedUpProblem {
            nodal_forces_n: &forces,
            ..p
        };
        let prepared = loaded_problem
            .prepare_with_gravity(&density, field, None, cx)
            .unwrap();
        let loaded = prepared.solve_vibration(modal_settings(), cx).unwrap();
        assert_eq!(loaded.mass_kg(), base.mass_kg());
        assert_eq!(loaded.reference_density_kg_m3(), &density);
        assert!(core::ptr::eq(
            loaded.gravity().unwrap().reference_density_kg_m3(),
            &density
        ));
        assert!(core::ptr::eq(loaded.gravity().unwrap().field(), field));
        assert_eq!(loaded.problem().nodal_forces_n, forces);
        for ((&k, &k0), &m) in loaded
            .stiffness_n_m()
            .iter()
            .zip(base.stiffness_n_m())
            .zip(base.mass_kg())
        {
            near(k - k0, 1.7 * m, 2e-12);
        }
        for (a, b) in loaded.modes().iter().zip(base.modes()) {
            near(
                a.displacement.lambda,
                b.displacement.lambda + 1.7,
                2e-8 * b.displacement.lambda,
            );
        }
        let replay = prepared
            .solve_vibration_from(
                &loaded.equilibrium().displacement_m,
                &loaded.equilibrium().pressure_pa,
                modal_settings(),
                cx,
            )
            .unwrap();
        assert_eq!(loaded.stiffness_n_m(), replay.stiffness_n_m());
    });
}

#[test]
fn exact_modes_keep_bulk_independence_gauge_and_gravity_frequency_shift() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let density = [7.0];
        let forces = harmonic_balance(&nodes, &elements[0], density[0]);
        for pin in [None, Some(0), Some(1)] {
            let supports = if pin.is_some() {
                (0..30)
                    .filter(|d| ![0, 3, 7, 11].contains(d))
                    .map(|d| (d, 0.0))
                    .collect()
            } else {
                clamp_face(&nodes)
            };
            let mut previous = None;
            for bulk in [20.0, 1e300] {
                let material = law(bulk);
                let mut p = problem(
                    &nodes,
                    &elements,
                    &material,
                    &supports,
                    &[],
                    VolumetricTreatment::ExactIncompressible,
                );
                p.pressure_gauge = pin.map(|node| (node, 0.0));
                let base = p
                    .solve_incompressible_vibration(&density, modal_settings(), cx)
                    .unwrap();
                let loaded_problem = MixedUpProblem {
                    nodal_forces_n: &forces,
                    ..p
                };
                let prepared = loaded_problem
                    .prepare_with_gravity(&density, &harmonic_field, None, cx)
                    .unwrap();
                let loaded = prepared
                    .solve_incompressible_vibration(modal_settings(), cx)
                    .unwrap();
                assert_eq!(loaded.mass_kg(), base.mass_kg());
                assert_eq!(loaded.problem().pressure_gauge, p.pressure_gauge);
                assert_eq!(
                    loaded.gravity().unwrap().reference_density_kg_m3(),
                    &density
                );
                for ((&k, &k0), &m) in loaded
                    .projected_stiffness_n_m()
                    .iter()
                    .zip(base.projected_stiffness_n_m())
                    .zip(base.projected_mass_kg())
                {
                    near(k - k0, 1.7 * m, 3e-12);
                }
                for (a, b) in loaded.modes().iter().zip(base.modes()) {
                    near(
                        a.coordinates.lambda,
                        b.coordinates.lambda + 1.7,
                        2e-8 * b.coordinates.lambda,
                    );
                    assert!(
                        a.constraint_relative_residual
                            <= modal_settings().pressure_residual_tolerance
                    );
                    assert!(
                        a.momentum_relative_residual
                            <= modal_settings().spectrum.relative_residual_tolerance
                    );
                    if let Some(node) = pin {
                        assert_eq!(a.pressure_shape[node], 0.0);
                    }
                }
                let pencil = loaded.projected_stiffness_n_m().to_vec();
                if let Some(ref previous) = previous {
                    assert_eq!(&pencil, previous);
                }
                previous = Some(pencil);
            }
        }
    });
}

#[test]
fn deformed_exact_equilibrium_combines_weight_pressure_and_full_support_reactions() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let supports = clamp_face(&nodes);
        let material = law(1e300);
        let density = [7.0];
        let field = |x: [f64; 3]| {
            Ok(GravitySample3 {
                potential_j_kg: 0.02 * x[2],
                acceleration_m_s2: [0.0, 0.0, -0.02],
                acceleration_gradient_s2: [[0.0; 3]; 3],
            })
        };
        let p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let prepared = p.prepare_with_gravity(&density, &field, None, cx).unwrap();
        let base = prepared.solve(cx).unwrap();
        let state = base.equilibrium();
        assert!(state.displacement_m.iter().any(|u| u[2] < -1e-5));
        let evaluation = prepared
            .evaluate(&state.displacement_m, &state.pressure_pa, false, cx)
            .unwrap();
        for r in evaluation.pressure_residual_m3 {
            near(r, 0.0, 2e-11);
        }
        for component in 0..3 {
            let reaction: f64 = state
                .reactions_n
                .iter()
                .filter(|(d, _)| d % 3 == component)
                .map(|(_, r)| r)
                .sum();
            near(
                reaction,
                if component == 2 {
                    density[0] * 0.02 / 6.0
                } else {
                    0.0
                },
                2e-10,
            );
        }
        let q = UniformBoundaryPressure::new(0.1).unwrap();
        let combined = p
            .prepare_with_gravity(&density, &field, Some(q), cx)
            .unwrap();
        let shifted: Vec<_> = state
            .pressure_pa
            .iter()
            .map(|v| v - q.pressure_pa())
            .collect();
        let result = combined
            .solve_from(&state.displacement_m, &shifted, cx)
            .unwrap();
        // In the exact formulation a uniform closed-domain q is balanced by
        // p -> p-q at the SAME deformation, without changing the weight.
        for (a, b) in result
            .equilibrium()
            .displacement_m
            .iter()
            .zip(&state.displacement_m)
        {
            for i in 0..3 {
                near(a[i], b[i], 2e-11);
            }
        }
        for (&a, &b) in result.equilibrium().pressure_pa.iter().zip(&shifted) {
            near(a, b, 2e-10);
        }
        for ((d, a), (d0, b)) in result
            .equilibrium()
            .reactions_n
            .iter()
            .zip(&state.reactions_n)
        {
            assert_eq!(d, d0);
            near(*a, *b, 2e-11);
        }
        assert_eq!(result.boundary_pressure(), Some(q));
        near(
            result.loading_potentials().gravitational_j,
            base.loading_potentials().gravitational_j,
            2e-12,
        );
        near(
            result.loading_potentials().boundary_pressure_j,
            q.pressure_pa() / 6.0,
            2e-11,
        );
        let modes = combined
            .solve_incompressible_vibration_from(
                &state.displacement_m,
                &shifted,
                modal_settings(),
                cx,
            )
            .unwrap();
        assert_eq!(modes.boundary_pressure(), Some(q));
        assert!(modes.gravity().is_some());
    });
}

#[test]
fn gravity_domain_density_gradient_and_cancellation_refusals_preserve_inputs() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(20.0);
        let supports = clamp_face(&nodes);
        let p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        for density in [
            vec![],
            vec![0.0],
            vec![-1.0],
            vec![f64::NAN],
            vec![f64::from_bits(1)],
            vec![1.0, 2.0],
        ] {
            assert!(matches!(
                p.prepare_with_gravity(&density, &affine_field, None, cx),
                Err(MixedUpError::InvalidInput { .. })
            ));
        }
        let density = [7.0];
        let refusing = Cell::new(false);
        let calls = Cell::new(0);
        let field = |x| {
            calls.set(calls.get() + 1);
            if refusing.get() {
                Err(HyperTetError::MaterialRefused {
                    element: 711,
                    what: "original source domain".into(),
                })
            } else {
                affine_field(x)
            }
        };
        let prepared = p.prepare_with_gravity(&density, &field, None, cx).unwrap();
        let u = vec![[0.0; 3]; nodes.len()];
        let pressure = vec![0.0; 4];
        let original_u = u.clone();
        let original_pressure = pressure.clone();
        let first = prepared.evaluate(&u, &pressure, true, cx).unwrap();
        refusing.set(true);
        assert!(matches!(
            prepared.solve_from(&u, &pressure, cx),
            Err(MixedUpError::Gravity {
                element: 0,
                point: 0,
                error: HyperTetError::MaterialRefused { element: 711, .. },
            })
        ));
        refusing.set(false);
        let retry = prepared.evaluate(&u, &pressure, true, cx).unwrap();
        assert_eq!(retry.internal_forces_n, first.internal_forces_n);
        assert_eq!(u, original_u);
        assert_eq!(pressure, original_pressure);
        let invalid_field = |x| {
            let mut sample = affine_field(x)?;
            sample.acceleration_gradient_s2[0][1] = 2.0;
            Ok(sample)
        };
        let invalid = p
            .prepare_with_gravity(&density, &invalid_field, None, cx)
            .unwrap();
        assert!(matches!(
            invalid.evaluate(&u, &pressure, true, cx),
            Err(MixedUpError::Gravity {
                error: HyperTetError::InvalidInput { .. },
                ..
            })
        ));
        let gate = CancelGate::new();
        let count = Cell::new(0);
        let cancelling = |x| {
            count.set(count.get() + 1);
            gate.request();
            affine_field(x)
        };
        let cancelled = p
            .prepare_with_gravity(&density, &cancelling, None, cx)
            .unwrap();
        with_gate(&gate, |cancel_cx| {
            assert!(matches!(
                cancelled.evaluate(&u, &pressure, true, cancel_cx),
                Err(MixedUpError::Cancelled)
            ));
        });
        assert_eq!(
            count.get(),
            1,
            "poll again before the next quadrature point"
        );
        // The callback can cancel without changing Cx. A later Newton trial
        // must not backtrack, call the source again and swallow cancellation.
        // One Tet10 at the declared 4-point rule has 64 points: 129 follows
        // its original residual and first full material/load linearization.
        for stop in [1, 129] {
            let callback_calls = Cell::new(0);
            let callback_cancel = |x: [f64; 3]| {
                callback_calls.set(callback_calls.get() + 1);
                if callback_calls.get() == stop {
                    Err(HyperTetError::Cancelled)
                } else {
                    Ok(GravitySample3 {
                        potential_j_kg: 0.02 * x[2],
                        acceleration_m_s2: [0.0, 0.0, -0.02],
                        acceleration_gradient_s2: [[0.0; 3]; 3],
                    })
                }
            };
            let direct = p
                .prepare_with_gravity(&density, &callback_cancel, None, cx)
                .unwrap();
            assert!(matches!(
                direct.solve_from(&u, &pressure, cx),
                Err(MixedUpError::Cancelled)
            ));
            assert_eq!(
                callback_calls.get(),
                stop,
                "a cancelled source is not retried"
            );
        }
        assert_eq!(u, original_u);
        assert_eq!(pressure, original_pressure);
        assert!(calls.get() > 1);
    });
}

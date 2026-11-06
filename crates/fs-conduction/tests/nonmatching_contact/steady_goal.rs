//! Independent two-slab goals: live residuals, full transpose, finite remainder.
use super::steady::{contacts, coupling, scalar_heat, solver_config};
use super::temperature::source_card;
use super::*;
use fs_conduction::adjoint::compare_discrete_goal;
use fs_conduction::interface::temperature::{TemperatureContact, TemperatureContacts};
use fs_conduction::radiation::{EnclosureRadiationGradient, OpenDiffuseEnclosure};

fn contact_mean(mesh: &ConductionMesh, split: usize) -> Vec<f64> {
    let (faces, _) = sides(mesh, split, 0.0);
    let mut weights = vec![0.0; mesh.vertex_count()];
    for index in faces {
        let face = &mesh.boundary()[index];
        for vertex in face.vertices {
            weights[vertex as usize] += face.area / 3.0;
        }
    }
    weights
}

fn slab_field(mesh: &ConductionMesh, split: usize, a: f64, b: f64, hot: f64) -> Vec<f64> {
    mesh.positions()
        .iter()
        .enumerate()
        .map(|(i, p)| {
            if i < split {
                300.0 + (a - 300.0) * p[0]
            } else {
                hot + (b - hot) * (2.0 - p[0])
            }
        })
        .collect()
}

fn goal_radiation_patch(name: &str, epsilon: f64, bath: f64) -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: epsilon,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 350.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic four-state radiation oracle, not measurements".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = fs_matdb::MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic surface".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    AmbientRadiationPatch::new(
        name,
        SurfaceEmissivity::from_card(name, &card, 320.0, SelectionPolicy::SingleClaimOnly).unwrap(),
        bath,
    )
    .unwrap()
}

fn card_with_value(
    value: PropertyValue,
    interpolation: InterpolationPolicy,
) -> InterfaceSystemCard {
    let dims = AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
            value,
            validity: ValidityDomain::unconstrained().with("T", 280.0, 380.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation,
            observations: vec![],
            provenance: Provenance {
                source: "analytic goal-comparison fixture, not measurements".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        texture_frame: "reference normal".into(),
    };
    InterfaceSystemCard::assemble(
        side("A"),
        side("B"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "test".into(),
            history: "unaged".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn mean_contacts<'m, 'c>(
    cx: &fs_exec::Cx<'_>,
    mesh: &'m ConductionMesh,
    bc: &ThermalBoundary,
    card: &'c InterfaceSystemCard,
) -> TemperatureContacts<'m, 'c> {
    let law = TemperatureContact::new(
        "bond",
        ThermalInterfaces::coincident_face_pairs(mesh).unwrap(),
        card,
        &QueryPoint::new().with("T", 320.0).unwrap(),
        "T",
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    TemperatureContacts::new(cx, mesh, bc, vec![], vec![law], 10_000).unwrap()
}

#[test]
fn g3_live_contact_goal_matches_independent_two_state_residual_and_complete_dual() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = boundary(&mesh);
            let card = source_card(1.0, false);
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &bc,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let heat = scalar_heat(false);
            let (a, b) = (300.0 + 0.1 * heat, 330.0 - 0.1 * heat);
            let reference = slab_field(&mesh, split, a, b, 330.0);
            let weights = contact_mean(&mesh, split);
            let resistance = |a, b| 0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0;
            let r = resistance(a, b);
            let qa = 1.0 / r - (a - b) / (1024.0 * r * r);
            let qb = -1.0 / r - 3.0 * (a - b) / (1024.0 * r * r);
            let determinant = 100.0 + 10.0 * (qa - qb);
            let lambda = [(10.0 - qb) / determinant, -qb / determinant];
            let nodal_residual = |a, b, r| {
                let flux = (a - b) / r;
                [10.0 * (a - 300.0) + flux, 10.0 * (b - 330.0) - flux]
            };
            let original = nodal_residual(a, b, r);
            let mut remainders = Vec::new();
            for fraction in [1.0, 0.5, 0.25] {
                let (aa, bb) = (a + 3.0 * fraction, b - 5.0 * fraction);
                let approximate = slab_field(&mesh, split, aa, bb, 330.0);
                let before = approximate.clone();
                let report = set
                    .compare_steady_goal(
                        cx,
                        problem,
                        &reference,
                        &approximate,
                        solver_config().linear,
                        1e-8,
                        1e-8,
                        &weights,
                        10_000,
                    )
                    .unwrap();
                let other = nodal_residual(aa, bb, resistance(aa, bb));
                let expected =
                    lambda[0] * (original[0] - other[0]) + lambda[1] * (original[1] - other[1]);
                near(report.signed_goal_change, -3.0 * fraction, 3e-12);
                near(report.signed_residual_change, expected, 2e-10);
                near(
                    report.signed_residual_change + report.linearization_remainder,
                    report.signed_goal_change,
                    3e-14,
                );
                near(report.nodal_contributions.iter().sum(), expected, 2e-10);
                assert!(report.uses_nonlinear_jacobian);
                assert!(report.primal_relative_residual < 1e-12);
                assert!(report.dual_relative_residual < solver_config().linear.tolerance);
                assert!(report.dual_iterations > 0);
                for &(vertex, _) in bc.dirichlet() {
                    assert_eq!(report.nodal_contributions[vertex], 0.0);
                }
                assert_eq!(approximate, before);
                remainders.push(report.linearization_remainder.abs());
                if fraction == 1.0 {
                    // Either freezing the residual law or freezing the dual
                    // loses a real first-order feedback contribution.
                    let frozen_other = nodal_residual(aa, bb, r);
                    let frozen_residual = lambda[0] * (original[0] - frozen_other[0])
                        + lambda[1] * (original[1] - frozen_other[1]);
                    assert!((expected - frozen_residual).abs() > 1e-3);
                    let physical = set.evaluate(cx, &bc, &reference).unwrap();
                    let frozen = compare_discrete_goal(
                        cx,
                        problem,
                        Some(physical.frozen_interfaces()),
                        solver_config().linear,
                        &reference,
                        &approximate,
                        &weights,
                    )
                    .unwrap();
                    assert!(
                        (frozen.signed_residual_change - report.signed_residual_change).abs()
                            > 1e-3
                    );
                    assert!(frozen.linearization_remainder.abs() < 1e-9);
                }
            }
            assert!(remainders[2] > 1e-5);
            for pair in remainders.windows(2) {
                assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{remainders:?}");
            }
            let zero = set
                .compare_steady_goal(
                    cx,
                    problem,
                    &reference,
                    &reference,
                    solver_config().linear,
                    1e-8,
                    1e-8,
                    &weights,
                    10_000,
                )
                .unwrap();
            assert_eq!(zero.signed_goal_change, 0.0);
            assert_eq!(zero.signed_residual_change, 0.0);
            assert_eq!(zero.linearization_remainder, 0.0);
        }
    });
}

#[test]
fn g3_live_contact_goal_keeps_nonlinear_bulk_feedback_in_the_same_transpose() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 1, 1], 0.0);
        let bc = boundary(&mesh);
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, true);
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(250.0, 5.0), (400.0, 20.0)]).unwrap(),
        );
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let solved = set
            .solve_steady(
                cx,
                problem,
                &vec![315.0; mesh.vertex_count()],
                solver_config(),
                coupling(),
            )
            .unwrap();
        let reference = solved.conduction.temperature;
        let dofs = DofMap::new(&bc, mesh.vertex_count()).unwrap();
        let weights = contact_mean(&mesh, split);
        let mut remainders = Vec::new();
        for fraction in [0.8, 0.4, 0.2] {
            let mut approximate = reference.clone();
            for &vertex in dofs.free() {
                approximate[vertex] += fraction * (1.0 + (vertex % 5) as f64);
            }
            let report = set
                .compare_steady_goal(
                    cx,
                    problem,
                    &reference,
                    &approximate,
                    solver_config().linear,
                    1e-8,
                    1e-8,
                    &weights,
                    10_000,
                )
                .unwrap();
            let observed: f64 = reference
                .iter()
                .zip(&approximate)
                .zip(&weights)
                .map(|((r, a), w)| w * (r - a))
                .sum();
            near(report.signed_goal_change, observed, 2e-12);
            near(
                report.signed_residual_change + report.linearization_remainder,
                observed,
                2e-12,
            );
            assert!(report.uses_nonlinear_jacobian);
            remainders.push(report.linearization_remainder.abs());
        }
        // Omitting K'(T), either contact partial, or the transpose leaves an
        // O(delta) defect. The correct complete tangent leaves O(delta^2).
        assert!(remainders[2] > 1e-6, "{remainders:?}");
        for pair in remainders.windows(2) {
            assert!((3.6..4.4).contains(&(pair[0] / pair[1])), "{remainders:?}");
        }
    });
}

#[test]
fn g1_constant_source_contact_goal_reduces_to_the_linear_identity() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 1, 1], 0.0);
        let bc = boundary(&mesh);
        let card = card_with_value(
            PropertyValue::Scalar {
                value: 0.5,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
        );
        let set = mean_contacts(cx, &mesh, &bc, &card);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let heat = 30.0 / 0.7;
        let reference = slab_field(
            &mesh,
            split,
            300.0 + heat / 10.0,
            330.0 - heat / 10.0,
            330.0,
        );
        let approximate = slab_field(&mesh, split, 309.0, 321.0, 330.0);
        let weights = contact_mean(&mesh, split);
        let report = set
            .compare_steady_goal(
                cx,
                problem,
                &reference,
                &approximate,
                solver_config().linear,
                1e-8,
                1e-8,
                &weights,
                10_000,
            )
            .unwrap();
        near(
            report.signed_residual_change,
            300.0 + heat / 10.0 - 309.0,
            2e-11,
        );
        near(report.linearization_remainder, 0.0, 2e-11);
        // Source tangents consistently use the complete FGMRES owner even
        // when this particular law has exactly zero temperature partials.
        assert!(report.uses_nonlinear_jacobian);
    });
}

#[test]
fn g0_live_contact_goal_rejects_dirty_reference_invalid_fields_and_budgets() {
    support::with_gate(|gate, cx| {
        let (mesh, split) = mesh([1, 2, 3], 0.0);
        let bc = boundary(&mesh);
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, false);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let heat = scalar_heat(false);
        let reference = slab_field(
            &mesh,
            split,
            300.0 + heat / 10.0,
            330.0 - heat / 10.0,
            330.0,
        );
        let weights = contact_mean(&mesh, split);
        let compare = |r: &[f64], a: &[f64], w: &[f64], linear, watts, cap| {
            set.compare_steady_goal(
                cx,
                problem,
                r,
                a,
                linear,
                watts,
                coupling().energy_tolerance_w,
                w,
                cap,
            )
        };
        let linear = solver_config().linear;
        let free = DofMap::new(&bc, mesh.vertex_count()).unwrap().free()[0];
        let mut dirty = reference.clone();
        dirty[free] += 1e-3;
        assert!(compare(&dirty, &reference, &weights, linear, 1e-8, 10_000).is_err());
        let mut invalid = reference.clone();
        let fixed = bc.dirichlet()[0].0;
        invalid[fixed] = invalid[fixed].next_up();
        assert!(compare(&reference, &invalid, &weights, linear, 1e-8, 10_000).is_err());
        assert!(compare(&invalid, &reference, &weights, linear, 1e-8, 10_000).is_err());
        for value in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
            invalid = reference.clone();
            invalid[free] = value;
            assert!(compare(&reference, &invalid, &weights, linear, 1e-8, 10_000).is_err());
        }
        assert!(
            compare(
                &reference,
                &reference[..reference.len() - 1],
                &weights,
                linear,
                1e-8,
                10_000
            )
            .is_err()
        );
        assert!(
            compare(
                &reference[..reference.len() - 1],
                &reference,
                &weights,
                linear,
                1e-8,
                10_000
            )
            .is_err()
        );
        assert!(
            compare(
                &reference,
                &reference,
                &weights[..weights.len() - 1],
                linear,
                1e-8,
                10_000
            )
            .is_err()
        );
        let mut bad_weights = weights.clone();
        bad_weights[free] = f64::NAN;
        assert!(compare(&reference, &reference, &bad_weights, linear, 1e-8, 10_000).is_err());
        for bad in [
            LinearConfig {
                restart: 0,
                ..linear
            },
            LinearConfig {
                max_iterations: 0,
                ..linear
            },
            LinearConfig {
                tolerance: f64::NAN,
                ..linear
            },
        ] {
            assert!(compare(&reference, &reference, &weights, bad, 1e-8, 10_000).is_err());
        }
        assert!(compare(&reference, &reference, &weights, linear, 0.0, 10_000).is_err());
        assert!(compare(&reference, &reference, &weights, linear, 1e-8, 0).is_err());
        // Source extrapolation at the approximate field must never silently
        // reuse the valid reference's resistance.
        invalid = slab_field(&mesh, split, 450.0, 450.0, 330.0);
        assert!(compare(&reference, &invalid, &weights, linear, 1e-8, 10_000).is_err());
        gate.request();
        assert!(matches!(
            compare(&reference, &reference, &weights, linear, 1e-8, 10_000),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

#[test]
fn g0_source_corner_refuses_at_reference_but_approximate_corner_is_evaluated() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 1, 1], 0.0);
        let card = card_with_value(
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: TEMPERATURE_DIMS,
                knots: vec![(280.0, 0.4), (315.0, 0.8), (380.0, 0.9)],
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            InterpolationPolicy::LinearInside,
        );
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let weights = contact_mean(&mesh, split);
        let bc = boundary(&mesh);
        let set = mean_contacts(cx, &mesh, &bc, &card);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let reference = slab_field(&mesh, split, 303.0, 327.0, 330.0);
        // R(315)=0.8 and H=30 satisfy the original heat equation exactly,
        // but the resistance slopes on the two sides of 315 K differ.
        assert!(set.evaluate(cx, &bc, &reference).is_ok());
        assert!(
            set.compare_steady_goal(
                cx,
                problem,
                &reference,
                &reference,
                solver_config().linear,
                1e-8,
                1e-8,
                &weights,
                10_000
            )
            .is_err()
        );
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "cold",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(300.0).unwrap(),
            )
            .unwrap()
            .region(
                "hot",
                |f| on_box_face(f.centroid[0], 2.0),
                ThermalBc::dirichlet(320.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let set = mean_contacts(cx, &mesh, &bc, &card);
        let problem = ConductionProblem {
            boundary: &bc,
            ..problem
        };
        let r = 0.4 + 30.0 * (0.8 - 0.4) / 35.0;
        let heat = 20.0 / (0.2 + r);
        let reference = slab_field(
            &mesh,
            split,
            300.0 + heat / 10.0,
            320.0 - heat / 10.0,
            320.0,
        );
        let approximate = slab_field(&mesh, split, 315.0, 315.0, 320.0);
        let report = set
            .compare_steady_goal(
                cx,
                problem,
                &reference,
                &approximate,
                solver_config().linear,
                1e-8,
                1e-8,
                &weights,
                10_000,
            )
            .unwrap();
        near(
            report.signed_goal_change,
            300.0 + heat / 10.0 - 315.0,
            2e-12,
        );
        near(
            report.signed_residual_change + report.linearization_remainder,
            report.signed_goal_change,
            2e-12,
        );
    });
}

// A four-temperature series circuit is independent of the FEM/contact/Robin
// assemblers. Its nonsymmetric transpose can be solved directly at this size.
fn solve_four_transpose(matrix: [[f64; 4]; 4], mut rhs: [f64; 4]) -> [f64; 4] {
    let mut a: [[f64; 4]; 4] = std::array::from_fn(|i| std::array::from_fn(|j| matrix[j][i]));
    for column in 0..4 {
        let pivot = (column..4)
            .max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))
            .unwrap();
        a.swap(column, pivot);
        rhs.swap(column, pivot);
        assert!(a[column][column].abs() > 1e-12);
        for row in column + 1..4 {
            let factor = a[row][column] / a[column][column];
            for j in column..4 {
                a[row][j] -= factor * a[column][j];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut solution = [0.0; 4];
    for row in (0..4).rev() {
        let rest: f64 = (row + 1..4).map(|j| a[row][j] * solution[j]).sum();
        solution[row] = (rhs[row] - rest) / a[row][row];
    }
    solution
}

fn check_four_state_goal(cx: &fs_exec::Cx<'_>, radiative: bool) {
    let resistance = |a, b| 0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0;
    let (mut lo, mut hi) = (0.0_f64, 75.0);
    for _ in 0..80 {
        let heat = 0.5 * (lo + hi);
        let r = resistance(300.0 + 0.2 * heat, 330.0 - 0.2 * heat);
        if heat * (0.4 + r) > 30.0 {
            hi = heat;
        } else {
            lo = heat;
        }
    }
    let heat = 0.5 * (lo + hi);
    let x = [
        300.0 + 0.1 * heat,
        300.0 + 0.2 * heat,
        330.0 - 0.2 * heat,
        330.0 - 0.1 * heat,
    ];
    let feedback = [0.07, 0.17, -0.09, -0.04];
    let epsilon = if radiative { [0.4, 0.6] } else { [0.0, 0.0] };
    let baths: [f64; 2] = [280.0, 350.0];
    // Direct fourth powers, not the production secant or its derivative.
    let radiation = |i: usize, wall: f64| {
        epsilon[i] * STEFAN_BOLTZMANN_W_M2_K4 * (wall.powi(4) - baths[i].powi(4))
    };
    let radiation_slope =
        |i: usize, wall: f64| 4.0 * epsilon[i] * STEFAN_BOLTZMANN_W_M2_K4 * wall.powi(3);
    // Manufacture exact convection reservoirs so the same physical
    // contact/bulk state balances with and without the radiative fluxes.
    let rc = [
        300.0 + radiation(0, x[0]) / 10.0,
        330.0 + radiation(1, x[3]) / 10.0,
    ];
    let points = [[10.0, rc[0], 0.25, 0.11], [10.0, rc[1], -0.35, 0.23]];
    let patches = if radiative {
        vec![
            goal_radiation_patch("cold reservoir", epsilon[0], baths[0]),
            goal_radiation_patch("hot reservoir", epsilon[1], baths[1]),
        ]
    } else {
        vec![]
    };
    let values = |state: [f64; 4]| {
        let u = state[0] - x[0];
        let v = state[3] - x[3];
        [
            [
                10.0 + 0.25 * u + 0.03 * u * u,
                rc[0] + (0.11 + feedback[0]) * u + feedback[1] * v,
            ],
            [
                10.0 - 0.35 * v + 0.02 * v * v,
                rc[1] + (0.23 + feedback[3]) * v + feedback[2] * u,
            ],
        ]
    };
    let residual = |state: [f64; 4], rows: [[f64; 2]; 2]| {
        let [c, a, b, d] = state;
        let q = (a - b) / resistance(a, b);
        [
            10.0 * (c - a) + rows[0][0] * (c - rows[0][1]) + radiation(0, c),
            10.0 * (a - c) + q,
            10.0 * (b - d) - q,
            10.0 * (d - b) + rows[1][0] * (d - rows[1][1]) + radiation(1, d),
        ]
    };
    let r = resistance(x[1], x[2]);
    let qa = 1.0 / r - (x[1] - x[2]) / (1024.0 * r * r);
    let qb = -1.0 / r - 3.0 * (x[1] - x[2]) / (1024.0 * r * r);
    let jacobian = [
        [
            20.0 + 0.25 * (x[0] - rc[0]) - 10.0 * (0.11 + feedback[0]) + radiation_slope(0, x[0]),
            -10.0,
            0.0,
            -10.0 * feedback[1],
        ],
        [-10.0, 10.0 + qa, qb, 0.0],
        [0.0, -qa, 10.0 - qb, -10.0],
        [
            -10.0 * feedback[2],
            0.0,
            -10.0,
            20.0 - 0.35 * (x[3] - rc[1]) - 10.0 * (0.23 + feedback[3]) + radiation_slope(1, x[3]),
        ],
    ];
    let lambda = solve_four_transpose(jacobian, [0.0, 1.0, 0.0, 0.0]);
    let reference_residual = residual(x, values(x));
    for matching in [true, false] {
        let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
        // Deliberately stale original values: both physical evaluations
        // must be bound from the supplied actual law, including the reference.
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "cold reservoir",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::robin(7.0, 285.0).unwrap(),
            )
            .unwrap()
            .region(
                "hot reservoir",
                |f| on_box_face(f.centroid[0], 2.0),
                ThermalBc::robin(14.0, 345.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let boundary_before = bc.clone();
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, matching);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let field = |state: [f64; 4]| {
            mesh.positions()
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    if i < split {
                        state[0] + (state[1] - state[0]) * p[0]
                    } else {
                        state[3] + (state[2] - state[3]) * (2.0 - p[0])
                    }
                })
                .collect::<Vec<_>>()
        };
        let reference = field(x);
        let weights = contact_mean(&mesh, split);
        let regions = ["cold reservoir", "hot reservoir"];
        let entries = 8 * mesh.vertex_count();
        let mut remainders = Vec::new();
        for fraction in [1.0, 0.5, 0.25] {
            let y = std::array::from_fn(|i| x[i] + fraction * [1.7, 3.0, -2.0, -1.1][i]);
            let approximate = field(y);
            let approximate_before = approximate.clone();
            let compare = |actual_values: &[[f64; 2]], map: &[f64]| {
                if radiative {
                    return set
                        .compare_steady_goal_with_ambient_radiation(
                            cx,
                            problem,
                            &reference,
                            &approximate,
                            solver_config().linear,
                            1e-8,
                            1e-8,
                            &regions,
                            &points,
                            actual_values,
                            map,
                            &patches,
                            &weights,
                            entries,
                            10_000,
                        )
                        .unwrap();
                }
                set.compare_steady_goal_with_robin_feedback(
                    cx,
                    problem,
                    &reference,
                    &approximate,
                    solver_config().linear,
                    1e-8,
                    1e-8,
                    &regions,
                    &points,
                    actual_values,
                    map,
                    &weights,
                    entries,
                    10_000,
                )
                .unwrap()
            };
            let report = compare(&values(y), &feedback);
            let other = residual(y, values(y));
            let expected: f64 = (0..4)
                .map(|i| lambda[i] * (reference_residual[i] - other[i]))
                .sum();
            near(report.signed_residual_change, expected, 3e-10);
            near(report.signed_goal_change, -3.0 * fraction, 3e-12);
            near(report.nodal_contributions.iter().sum(), expected, 3e-10);
            near(
                report.signed_residual_change + report.linearization_remainder,
                report.signed_goal_change,
                3e-14,
            );
            assert!(report.uses_nonlinear_jacobian);
            assert!(report.primal_relative_residual < 1e-12);
            assert!(report.dual_relative_residual < solver_config().linear.tolerance);
            remainders.push(report.linearization_remainder.abs());
            if fraction == 1.0 {
                let frozen_boundary = compare(&values(x), &feedback);
                assert!((frozen_boundary.signed_residual_change - expected).abs() > 0.1);
                let missing_map = compare(&values(y), &[0.0; 4]);
                assert!((missing_map.signed_residual_change - expected).abs() > 1e-3);
                if radiative {
                    // Retain the exact reference tangent/contact, but freeze
                    // only the approximate radiative secants. A plausible
                    // combined-Robin surrogate must fail this same oracle.
                    let mut combined_points = [[0.0; 4]; 2];
                    let mut frozen_values = [[0.0; 2]; 2];
                    let mut combined_map = [0.0; 4];
                    for i in 0..2 {
                        let wall = x[3 * i];
                        let scale = epsilon[i] * STEFAN_BOLTZMANN_W_M2_K4;
                        let hr = scale * (wall + baths[i]) * (wall * wall + baths[i] * baths[i]);
                        let dhr = scale
                            * (3.0 * wall * wall + 2.0 * wall * baths[i] + baths[i] * baths[i]);
                        let [hc, tc, dhc, dtc] = points[i];
                        let h = hc + hr;
                        let t = (hc / h) * tc + (hr / h) * baths[i];
                        combined_points[i] = [
                            h,
                            t,
                            dhc + dhr,
                            (dhc / h) * (tc - t) + (dhr / h) * (baths[i] - t) + (hc / h) * dtc,
                        ];
                        for j in 0..2 {
                            combined_map[2 * i + j] = (hc / h) * feedback[2 * i + j];
                        }
                        let [ha, ta] = values(y)[i];
                        frozen_values[i] =
                            [ha + hr, (ha / (ha + hr)) * ta + (hr / (ha + hr)) * baths[i]];
                    }
                    let wrong = set
                        .compare_steady_goal_with_robin_feedback(
                            cx,
                            problem,
                            &reference,
                            &approximate,
                            solver_config().linear,
                            1e-8,
                            1e-8,
                            &regions,
                            &combined_points,
                            &frozen_values,
                            &combined_map,
                            &weights,
                            entries,
                            10_000,
                        )
                        .unwrap();
                    assert!((wrong.signed_residual_change - expected).abs() > 1e-3);
                }
            }
            assert_eq!(approximate, approximate_before);
        }
        assert!(remainders[2] > 1e-4);
        for pair in remainders.windows(2) {
            assert!((3.6..4.5).contains(&(pair[0] / pair[1])), "{remainders:?}");
        }
        assert_eq!(bc, boundary_before);
    }
}

#[test]
fn g3_live_contact_robin_goal_matches_complete_four_state_transpose_and_both_law_values() {
    with_cx(|cx| check_four_state_goal(cx, false));
}

#[test]
fn g3_live_contact_ambient_goal_matches_four_state_transpose_and_live_radiative_residual() {
    with_cx(|cx| check_four_state_goal(cx, true));
}

#[test]
fn g3_live_contact_radiation_only_goal_accepts_physical_source_endpoint_but_retains_gates() {
    support::with_gate(|gate, cx| {
        let (mesh, split) = mesh([1, 1, 1], 0.0);
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "cold",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(300.0).unwrap(),
            )
            .unwrap()
            .region(
                "exterior",
                |f| on_box_face(f.centroid[0], 2.0),
                ThermalBc::neumann(0.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, true);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let patch = goal_radiation_patch("exterior", 0.6, 330.0);
        // Independent series balance: H crosses both bulk slabs and the
        // original two-sided source, then equals the direct fourth-power heat.
        let states = |heat: f64| {
            let a = 300.0 + 0.1 * heat;
            let b = (a * (1.0 + heat / 1024.0) - 0.75 * heat) / (1.0 - 3.0 * heat / 1024.0);
            [a, b, b + 0.1 * heat]
        };
        let (mut lo, mut hi) = (0.0_f64, 50.0);
        for _ in 0..80 {
            let heat = 0.5 * (lo + hi);
            let [_, _, wall] = states(heat);
            if heat > 0.6 * STEFAN_BOLTZMANN_W_M2_K4 * (330.0_f64.powi(4) - wall.powi(4)) {
                hi = heat;
            } else {
                lo = heat;
            }
        }
        let [a, b, wall] = states(0.5 * (lo + hi));
        let reference = slab_field(&mesh, split, a, b, wall);
        let approximate = slab_field(&mesh, split, a, b, 350.0);
        assert!(patch.secant_coefficient_w_m2_k(350.0).is_ok());
        assert!(patch.secant_partials_w_m2_k2(350.0).is_err());
        let weights = contact_mean(&mesh, split);
        let point = [0.0, 330.0, 0.0, 0.0];
        let value = [0.0, 330.0];
        let entries = 4 * mesh.vertex_count();
        let compare = |reference: &[f64],
                       approximate: &[f64],
                       points: &[[f64; 4]],
                       values: &[[f64; 2]],
                       map: &[f64],
                       patches: &[AmbientRadiationPatch],
                       weights: &[f64],
                       cap| {
            set.compare_steady_goal_with_ambient_radiation(
                cx,
                problem,
                reference,
                approximate,
                solver_config().linear,
                1e-8,
                1e-8,
                &["exterior"],
                points,
                values,
                map,
                patches,
                weights,
                cap,
                10_000,
            )
        };
        let patches = [patch.clone()];
        let report = compare(
            &reference,
            &approximate,
            &[point],
            &[value],
            &[0.0],
            &patches,
            &weights,
            entries,
        )
        .unwrap();
        near(report.signed_goal_change, 0.0, 1e-14);
        assert!(report.signed_residual_change.abs() > 0.1);
        near(
            report.signed_residual_change + report.linearization_remainder,
            0.0,
            1e-14,
        );
        assert!(report.uses_nonlinear_jacobian);
        assert!(report.primal_relative_residual < 1e-12);
        assert!(report.dual_relative_residual < solver_config().linear.tolerance);
        // This trace has no convection. Arbitrary finite positive placeholders
        // cannot overwrite the retained radiative reservoir or create air heat.
        let placeholder = compare(
            &reference,
            &approximate,
            &[[0.0, 290.0, 0.0, 0.0]],
            &[[0.0, 345.0]],
            &[0.0],
            &patches,
            &weights,
            entries,
        )
        .unwrap();
        assert_eq!(report, placeholder);
        for bad in [
            [1.0, 330.0, 0.0, 0.0],
            [0.0, 330.0, 1.0, 0.0],
            [0.0, 330.0, 0.0, 1.0],
            [0.0, f64::NAN, 0.0, 0.0],
        ] {
            assert!(
                compare(
                    &reference,
                    &approximate,
                    &[bad],
                    &[value],
                    &[0.0],
                    &patches,
                    &weights,
                    entries
                )
                .is_err()
            );
        }
        assert!(
            compare(
                &reference,
                &approximate,
                &[point],
                &[[1.0, 330.0]],
                &[0.0],
                &patches,
                &weights,
                entries
            )
            .is_err()
        );
        assert!(
            compare(
                &reference,
                &approximate,
                &[point],
                &[value],
                &[1.0],
                &patches,
                &weights,
                entries
            )
            .is_err()
        );
        assert!(
            compare(
                &reference,
                &approximate,
                &[point],
                &[value],
                &[0.0],
                &[],
                &weights,
                entries
            )
            .is_err()
        );
        assert!(
            compare(
                &reference,
                &approximate,
                &[point],
                &[value],
                &[0.0],
                &[patch.clone(), patch],
                &weights,
                entries
            )
            .is_err()
        );
        assert!(
            compare(
                &reference,
                &approximate,
                &[point],
                &[value],
                &[0.0],
                &patches,
                &weights,
                entries - 1
            )
            .is_err()
        );
        assert!(
            compare(
                &reference,
                &slab_field(&mesh, split, a, b, 351.0),
                &[point],
                &[value],
                &[0.0],
                &patches,
                &weights,
                entries
            )
            .is_err()
        );
        let zero = vec![0.0; mesh.vertex_count()];
        // Zero objectives do not launder an unbalanced reference or permit
        // an unsupported physical approximate source query.
        assert!(
            compare(
                &slab_field(&mesh, split, a, b, wall + 1.0),
                &reference,
                &[point],
                &[value],
                &[0.0],
                &patches,
                &zero,
                entries
            )
            .is_err()
        );
        assert!(
            compare(
                &reference,
                &slab_field(&mesh, split, a, b, 351.0),
                &[point],
                &[value],
                &[0.0],
                &patches,
                &zero,
                entries
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            compare(
                &reference,
                &approximate,
                &[point],
                &[value],
                &[0.0],
                &patches,
                &weights,
                entries
            ),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

#[test]
fn g0_live_contact_robin_goal_rejects_malformed_stale_and_cancelled_feedback() {
    support::with_gate(|gate, cx| {
        let (mesh, split) = mesh([1, 1, 1], 0.0);
        let bc = super::steady::reservoir_boundary(&mesh, true);
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, true);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let heat = scalar_heat(true);
        let reference = slab_field(
            &mesh,
            split,
            300.0 + 0.1 * heat,
            330.0 - 0.2 * heat,
            330.0 - 0.1 * heat,
        );
        let weights = contact_mean(&mesh, split);
        let point = [10.0, 330.0, 0.0, 0.0];
        let value = [10.0, 330.0];
        let entries = 4 * mesh.vertex_count();
        let compare = |names: &[&str],
                       points: &[[f64; 4]],
                       values: &[[f64; 2]],
                       map: &[f64],
                       weights: &[f64],
                       cap| {
            set.compare_steady_goal_with_robin_feedback(
                cx,
                problem,
                &reference,
                &reference,
                solver_config().linear,
                1e-8,
                1e-8,
                names,
                points,
                values,
                map,
                weights,
                cap,
                10_000,
            )
        };
        let region = ["hot reservoir"];
        let good = compare(&region, &[point], &[value], &[0.0], &weights, entries).unwrap();
        assert_eq!(good.signed_residual_change, 0.0);
        assert_eq!(good.linearization_remainder, 0.0);
        assert!(compare(&region, &[point], &[value], &[0.0], &weights, entries - 1).is_err());
        assert!(compare(&region, &[], &[value], &[0.0], &weights, entries).is_err());
        assert!(compare(&region, &[point], &[], &[0.0], &weights, entries).is_err());
        assert!(compare(&region, &[point], &[value], &[], &weights, entries).is_err());
        assert!(compare(&region, &[point], &[value], &[0.0, 0.0], &weights, entries).is_err());
        assert!(compare(&region, &[point], &[value], &[f64::NAN], &weights, entries).is_err());
        assert!(compare(&["absent"], &[point], &[value], &[0.0], &weights, entries).is_err());
        assert!(compare(&["cold"], &[point], &[value], &[0.0], &weights, entries).is_err());
        assert!(
            compare(
                &["hot reservoir"; 2],
                &[point; 2],
                &[value; 2],
                &[0.0; 4],
                &weights,
                2 * entries
            )
            .is_err()
        );
        assert!(
            compare(
                &["hot reservoir"; 65],
                &[point; 65],
                &[value; 65],
                &[0.0; 65 * 65],
                &weights,
                65 * entries
            )
            .is_err()
        );
        for index in 0..4 {
            let mut bad = point;
            bad[index] = f64::NAN;
            assert!(compare(&region, &[bad], &[value], &[0.0], &weights, entries).is_err());
        }
        for bad in [
            [0.0, 330.0],
            [-1.0, 330.0],
            [f64::INFINITY, 330.0],
            [10.0, f64::NAN],
        ] {
            assert!(compare(&region, &[point], &[bad], &[0.0], &weights, entries).is_err());
        }
        assert!(
            compare(
                &region,
                &[[0.0, 330.0, 0.0, 0.0]],
                &[value],
                &[0.0],
                &weights,
                entries
            )
            .is_err()
        );
        // The reference's original watt gates still run when the objective
        // vanishes; stale law values do not become a checked physical state.
        assert!(
            compare(
                &region,
                &[[10.0, 340.0, 0.0, 0.0]],
                &[value],
                &[0.0],
                &vec![0.0; mesh.vertex_count()],
                entries
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            compare(&region, &[point], &[value], &[0.0], &weights, entries),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

// Independent closed exchange and explicit two-by-two open radiosity. The
// open model includes both self reflection and exchange through a black bath.
fn enclosure_oracle(open: bool, walls: [f64; 2], bath: f64) -> ([f64; 2], [[f64; 2]; 2]) {
    let sigma = STEFAN_BOLTZMANN_W_M2_K4;
    let epsilon = [0.4, 0.6];
    if !open {
        let scale = sigma / (1.0 / epsilon[0] + 1.0 / epsilon[1] - 1.0);
        let heat = scale * (walls[0].powi(4) - walls[1].powi(4));
        let slope = [
            4.0 * scale * walls[0].powi(3),
            -4.0 * scale * walls[1].powi(3),
        ];
        return ([heat, -heat], [slope, [-slope[0], -slope[1]]]);
    }
    let factors = [[0.1, 0.25], [0.25, 0.15]];
    let escape = [0.65, 0.60];
    let matrix: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (if i == j { 1.0 } else { 0.0 }) - (1.0 - epsilon[i]) * factors[i][j]
        })
    });
    let determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    let inverse = [
        [matrix[1][1] / determinant, -matrix[0][1] / determinant],
        [-matrix[1][0] / determinant, matrix[0][0] / determinant],
    ];
    let rhs: [f64; 2] = std::array::from_fn(|i| {
        epsilon[i] * sigma * walls[i].powi(4)
            + (1.0 - epsilon[i]) * escape[i] * sigma * bath.powi(4)
    });
    let radiosity: [f64; 2] =
        std::array::from_fn(|i| inverse[i][0] * rhs[0] + inverse[i][1] * rhs[1]);
    let heat = std::array::from_fn(|i| {
        radiosity[i]
            - factors[i][0] * radiosity[0]
            - factors[i][1] * radiosity[1]
            - escape[i] * sigma * bath.powi(4)
    });
    let derivative = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (inverse[i][j] - factors[i][0] * inverse[0][j] - factors[i][1] * inverse[1][j])
                * epsilon[j]
                * 4.0
                * sigma
                * walls[j].powi(3)
        })
    });
    (heat, derivative)
}

fn check_four_state_enclosure_goal(cx: &fs_exec::Cx<'_>, open: bool, live: bool, matching: bool) {
    let resistance = |a, b| 0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0;
    let (mut lo, mut hi) = (0.0_f64, 75.0);
    for _ in 0..80 {
        let heat = 0.5 * (lo + hi);
        if heat * (0.4 + resistance(300.0 + 0.2 * heat, 330.0 - 0.2 * heat)) > 30.0 {
            hi = heat;
        } else {
            lo = heat;
        }
    }
    let heat = 0.5 * (lo + hi);
    let x = [
        300.0 + 0.1 * heat,
        300.0 + 0.2 * heat,
        330.0 - 0.2 * heat,
        330.0 - 0.1 * heat,
    ];
    let bath = 310.0;
    let (radiation, radiation_jacobian) = enclosure_oracle(open, [x[0], x[3]], bath);
    let rc = [300.0 + radiation[0] / 10.0, 330.0 + radiation[1] / 10.0];
    let feedback = [0.07, 0.17, -0.09, -0.04];
    let points = [[10.0, rc[0], 0.25, 0.11], [10.0, rc[1], -0.35, 0.23]];
    let values = |state: [f64; 4]| {
        let u = state[0] - x[0];
        let v = state[3] - x[3];
        [
            [
                10.0 + 0.25 * u + 0.03 * u * u,
                rc[0] + (0.11 + feedback[0]) * u + feedback[1] * v,
            ],
            [
                10.0 - 0.35 * v + 0.02 * v * v,
                rc[1] + (0.23 + feedback[3]) * v + feedback[2] * u,
            ],
        ]
    };
    let r = resistance(x[1], x[2]);
    let residual = |state: [f64; 4], rows: [[f64; 2]; 2]| {
        let [c, a, b, d] = state;
        let q = (a - b) / if live { resistance(a, b) } else { r };
        let (radiation, _) = enclosure_oracle(open, [c, d], bath);
        [
            10.0 * (c - a) + rows[0][0] * (c - rows[0][1]) + radiation[0],
            10.0 * (a - c) + q,
            10.0 * (b - d) - q,
            10.0 * (d - b) + rows[1][0] * (d - rows[1][1]) + radiation[1],
        ]
    };
    let qa = 1.0 / r
        - if live {
            (x[1] - x[2]) / (1024.0 * r * r)
        } else {
            0.0
        };
    let qb = -1.0 / r
        - if live {
            3.0 * (x[1] - x[2]) / (1024.0 * r * r)
        } else {
            0.0
        };
    let jacobian = [
        [
            20.0 + 0.25 * (x[0] - rc[0]) - 10.0 * (0.11 + feedback[0]) + radiation_jacobian[0][0],
            -10.0,
            0.0,
            -10.0 * feedback[1] + radiation_jacobian[0][1],
        ],
        [-10.0, 10.0 + qa, qb, 0.0],
        [0.0, -qa, 10.0 - qb, -10.0],
        [
            -10.0 * feedback[2] + radiation_jacobian[1][0],
            0.0,
            -10.0,
            20.0 - 0.35 * (x[3] - rc[1]) - 10.0 * (0.23 + feedback[3]) + radiation_jacobian[1][1],
        ],
    ];
    let lambda = solve_four_transpose(jacobian, [0.0, 1.0, 0.0, 0.0]);
    let reference_residual = residual(x, values(x));
    let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
    let bc = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "cold reservoir",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::robin(7.0, 285.0).unwrap(),
        )
        .unwrap()
        .region(
            "hot reservoir",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::robin(14.0, 345.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let boundary_before = bc.clone();
    let card = source_card(1.0, false);
    let set = contacts(cx, &mesh, &bc, &card, split, matching);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &bc,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let surfaces = [0.0, 2.0]
        .into_iter()
        .enumerate()
        .map(|(i, coordinate)| {
            let name = ["cold reservoir", "hot reservoir"][i];
            RadiationSurface::new(
                &mesh,
                name,
                |f| on_box_face(f.centroid[0], coordinate),
                goal_radiation_patch(name, [0.4, 0.6][i], bath)
                    .emissivity()
                    .clone(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let closed = GrayDiffuseEnclosure::new(
        surfaces.clone(),
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
    )
    .unwrap();
    let opened = OpenDiffuseEnclosure::new(
        cx,
        surfaces,
        vec![vec![0.1, 0.25], vec![0.25, 0.15]],
        ViewFactorEvidence::Analytic {
            geometry: "synthetic reciprocal two-patch goal oracle".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    let field = |state: [f64; 4]| {
        mesh.positions()
            .iter()
            .enumerate()
            .map(|(i, p)| {
                if i < split {
                    state[0] + (state[1] - state[0]) * p[0]
                } else {
                    state[3] + (state[2] - state[3]) * (2.0 - p[0])
                }
            })
            .collect::<Vec<_>>()
    };
    let reference = field(x);
    // This baseline deliberately represents a constant physical resistance.
    // The live cases must differ: they query the two-sided source at both fields.
    let frozen = set.evaluate(cx, &bc, &reference).unwrap();
    let weights = contact_mean(&mesh, split);
    let regions = ["cold reservoir", "hot reservoir"];
    let entries = 12 * mesh.vertex_count(); // 4*n*2 Robin + 2*n*2 radiation.
    let mut remainders = Vec::new();
    for fraction in [1.0, 0.5, 0.25] {
        let y = std::array::from_fn(|i| x[i] + fraction * [1.7, 3.0, -2.0, -1.1][i]);
        let approximate = field(y);
        let approximate_before = approximate.clone();
        let compare = |rows: &[[f64; 2]], map: &[f64], live_contacts| {
            match (open, live_contacts) {
                (false, true) => set.compare_steady_goal_with_gray_diffuse_enclosure(
                    cx,
                    problem,
                    &reference,
                    &approximate,
                    solver_config().linear,
                    1e-8,
                    1e-8,
                    &regions,
                    &points,
                    rows,
                    map,
                    &closed,
                    1e-9,
                    &weights,
                    entries,
                    10_000,
                ),
                (true, true) => set.compare_steady_goal_with_open_gray_diffuse_enclosure(
                    cx,
                    problem,
                    &reference,
                    &approximate,
                    solver_config().linear,
                    1e-8,
                    1e-8,
                    &regions,
                    &points,
                    rows,
                    map,
                    &opened,
                    bath,
                    1e-9,
                    &weights,
                    entries,
                    10_000,
                ),
                (false, false) => EnclosureRadiationGradient::compare_goal_at(
                    cx,
                    problem,
                    Some(frozen.frozen_interfaces()),
                    solver_config().linear,
                    &reference,
                    &approximate,
                    1e-8,
                    1e-8,
                    &regions,
                    &points,
                    rows,
                    map,
                    &closed,
                    1e-9,
                    &weights,
                    entries,
                ),
                (true, false) => EnclosureRadiationGradient::compare_open_goal_at(
                    cx,
                    problem,
                    Some(frozen.frozen_interfaces()),
                    solver_config().linear,
                    &reference,
                    &approximate,
                    1e-8,
                    1e-8,
                    &regions,
                    &points,
                    rows,
                    map,
                    &opened,
                    bath,
                    1e-9,
                    &weights,
                    entries,
                ),
            }
            .unwrap()
        };
        let report = compare(&values(y), &feedback, live);
        let other = residual(y, values(y));
        let difference: [f64; 4] = std::array::from_fn(|i| reference_residual[i] - other[i]);
        let expected: f64 = (0..4).map(|i| lambda[i] * difference[i]).sum();
        near(report.signed_residual_change, expected, 3e-10);
        near(report.nodal_contributions.iter().sum(), expected, 3e-10);
        near(report.signed_goal_change, -3.0 * fraction, 3e-12);
        near(
            report.signed_residual_change + report.linearization_remainder,
            report.signed_goal_change,
            3e-14,
        );
        assert!(report.primal_relative_residual < 1e-12);
        assert!(report.dual_relative_residual < solver_config().linear.tolerance);
        if live {
            assert!(report.uses_nonlinear_jacobian);
        }
        remainders.push(report.linearization_remainder.abs());
        if fraction == 1.0 {
            assert!(
                (compare(&values(x), &feedback, live).signed_residual_change - expected).abs()
                    > 0.05
            );
            assert!(
                (compare(&values(y), &[0.0; 4], live).signed_residual_change - expected).abs()
                    > 1e-3
            );
            if live {
                assert!(
                    (compare(&values(y), &feedback, false).signed_residual_change - expected).abs()
                        > 1e-4
                );
            }
            // A local-only radiation tangent misses reflected cross-patch heat.
            let mut local_only = jacobian;
            local_only[0][3] -= radiation_jacobian[0][1];
            local_only[3][0] -= radiation_jacobian[1][0];
            let wrong_dual = solve_four_transpose(local_only, [0.0, 1.0, 0.0, 0.0]);
            let wrong: f64 = (0..4).map(|i| wrong_dual[i] * difference[i]).sum();
            assert!((wrong - expected).abs() > 1e-3);
        }
        assert_eq!(approximate, approximate_before);
    }
    assert!(remainders[2] > 1e-5);
    for pair in remainders.windows(2) {
        assert!((3.5..4.6).contains(&(pair[0] / pair[1])), "{remainders:?}");
    }
    assert_eq!(reference, field(x));
    assert_eq!(bc, boundary_before);
    if open && live && matching {
        // A shifted Robin reservoir rounds to 300 K despite a real outward
        // load. Zero goal weights cannot bypass the physical watt admission.
        let isothermal = vec![300.0; mesh.vertex_count()];
        let radiation = enclosure_oracle(true, [300.0; 2], 290.0).0;
        for q in radiation {
            assert!(q > 1.0);
            assert_eq!(300.0 - q / 1e20, 300.0);
        }
        assert!(
            set.compare_steady_goal_with_open_gray_diffuse_enclosure(
                cx,
                problem,
                &isothermal,
                &isothermal,
                solver_config().linear,
                1e-8,
                1e-8,
                &regions,
                &[[1e20, 300.0, 0.0, 0.0]; 2],
                &[[1e20, 300.0]; 2],
                &[0.0; 4],
                &opened,
                290.0,
                1e-9,
                &vec![0.0; mesh.vertex_count()],
                entries,
                10_000
            )
            .is_err()
        );
    }
}

#[test]
fn g3_live_contact_enclosure_goals_match_reflected_four_state_oracles() {
    with_cx(|cx| {
        for open in [false, true] {
            for matching in [true, false] {
                check_four_state_enclosure_goal(cx, open, true, matching);
            }
        }
    });
}

#[test]
fn g3_fixed_contact_enclosure_goals_retain_both_physical_fields_and_reflection() {
    with_cx(|cx| {
        for open in [false, true] {
            check_four_state_enclosure_goal(cx, open, false, true);
        }
    });
}

#[test]
fn g3_open_enclosure_goal_keeps_prescribed_flux_without_robin_and_physical_source_endpoints() {
    support::with_gate(|gate, cx| {
        let (mesh, split) = mesh([1, 1, 1], 0.0);
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "cold",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(300.0).unwrap(),
            )
            .unwrap()
            .region(
                "exterior",
                |f| on_box_face(f.centroid[0], 2.0),
                ThermalBc::neumann(2.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let boundary_before = bc.clone();
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, true);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let patch = goal_radiation_patch("exterior", 0.6, 330.0);
        let enclosure = |mesh: &ConductionMesh| {
            let surface = RadiationSurface::new(
                mesh,
                "exterior",
                |f| on_box_face(f.centroid[0], 2.0),
                patch.emissivity().clone(),
            )
            .unwrap();
            OpenDiffuseEnclosure::new(
                cx,
                vec![surface],
                vec![vec![0.3]],
                ViewFactorEvidence::Analytic {
                    geometry: "synthetic self-reflecting aperture".into(),
                },
                ViewFactorTolerance::default(),
            )
            .unwrap()
        };
        let model = enclosure(&mesh);
        // Eliminate the one-patch radiosity including repeated self reflection.
        // The authored 2 W outward flux remains in the original heat balance.
        let effective = 0.6 * 0.7 / (1.0 - 0.4 * 0.3);
        let radiation =
            |wall: f64| effective * STEFAN_BOLTZMANN_W_M2_K4 * (wall.powi(4) - 330.0_f64.powi(4));
        let states = |heat: f64| {
            let a = 300.0 + 0.1 * heat;
            let b = (a * (1.0 + heat / 1024.0) - 0.75 * heat) / (1.0 - 3.0 * heat / 1024.0);
            [a, b, b + 0.1 * heat]
        };
        let (mut lo, mut hi) = (0.0_f64, 50.0);
        for _ in 0..80 {
            let heat = 0.5 * (lo + hi);
            if heat + 2.0 + radiation(states(heat)[2]) > 0.0 {
                hi = heat;
            } else {
                lo = heat;
            }
        }
        let heat = 0.5 * (lo + hi);
        let [a, b, wall] = states(heat);
        near(heat + 2.0 + radiation(wall), 0.0, 1e-11);
        let reference = slab_field(&mesh, split, a, b, wall);
        let approximate = slab_field(&mesh, split, a, b, 350.0);
        assert!(patch.secant_coefficient_w_m2_k(350.0).is_ok());
        assert!(patch.secant_partials_w_m2_k2(350.0).is_err());
        let weights = contact_mean(&mesh, split);
        let zero = vec![0.0; mesh.vertex_count()];
        let entries = 2 * mesh.vertex_count(); // No Robin controls, one radiation patch.
        let compare = |reference: &[f64],
                       approximate: &[f64],
                       model: &OpenDiffuseEnclosure,
                       weights: &[f64],
                       cap,
                       pairs| {
            set.compare_steady_goal_with_open_gray_diffuse_enclosure(
                cx,
                problem,
                reference,
                approximate,
                solver_config().linear,
                1e-8,
                1e-8,
                &[],
                &[],
                &[],
                &[],
                model,
                330.0,
                1e-9,
                weights,
                cap,
                pairs,
            )
        };
        let report = compare(&reference, &approximate, &model, &weights, entries, 10_000).unwrap();
        near(report.signed_goal_change, 0.0, 1e-14);
        assert!(report.signed_residual_change.abs() > 0.1);
        near(
            report.signed_residual_change + report.linearization_remainder,
            0.0,
            1e-14,
        );
        assert!(report.uses_nonlinear_jacobian);
        assert!(report.primal_relative_residual < 1e-12);
        assert!(report.dual_relative_residual < solver_config().linear.tolerance);
        let same = compare(&reference, &reference, &model, &zero, entries, 10_000).unwrap();
        assert_eq!(same.signed_goal_change, 0.0);
        assert_eq!(same.signed_residual_change, 0.0);
        assert!(
            compare(
                &reference,
                &approximate,
                &model,
                &weights,
                entries - 1,
                10_000
            )
            .is_err()
        );
        assert!(compare(&reference, &approximate, &model, &weights, entries, 0).is_err());
        let unsupported = slab_field(&mesh, split, a, b, 351.0);
        assert!(compare(&reference, &unsupported, &model, &zero, entries, 10_000).is_err());
        let unbalanced = slab_field(&mesh, split, a, b, wall + 1.0);
        assert!(compare(&unbalanced, &reference, &model, &zero, entries, 10_000).is_err());
        // Names only select original Robin convection. A radiating Neumann
        // trace must not be promoted to an invented zero-coefficient film.
        for name in ["exterior", "absent"] {
            assert!(
                set.compare_steady_goal_with_open_gray_diffuse_enclosure(
                    cx,
                    problem,
                    &reference,
                    &approximate,
                    solver_config().linear,
                    1e-8,
                    1e-8,
                    &[name],
                    &[[10.0, 330.0, 0.0, 0.0]],
                    &[[10.0, 330.0]],
                    &[0.0],
                    &model,
                    330.0,
                    1e-9,
                    &weights,
                    6 * mesh.vertex_count(),
                    10_000
                )
                .is_err()
            );
        }
        let (other_mesh, _) = super::mesh([1, 2, 3], 0.0);
        let stale = enclosure(&other_mesh);
        assert!(compare(&reference, &approximate, &stale, &zero, entries, 10_000).is_err());
        // A perfectly balanced reference at the emissivity support endpoint
        // still lacks the required tangent, even when the objective vanishes.
        let edge_bc = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "cold",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(350.0).unwrap(),
            )
            .unwrap()
            .region(
                "exterior",
                |f| on_box_face(f.centroid[0], 2.0),
                ThermalBc::neumann(0.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let edge_contacts = contacts(cx, &mesh, &edge_bc, &card, split, true);
        let edge = vec![350.0; mesh.vertex_count()];
        assert!(
            edge_contacts
                .compare_steady_goal_with_open_gray_diffuse_enclosure(
                    cx,
                    ConductionProblem {
                        boundary: &edge_bc,
                        ..problem
                    },
                    &edge,
                    &edge,
                    solver_config().linear,
                    1e-8,
                    1e-8,
                    &[],
                    &[],
                    &[],
                    &[],
                    &model,
                    350.0,
                    1e-9,
                    &zero,
                    entries,
                    10_000
                )
                .is_err()
        );
        assert_eq!(reference, slab_field(&mesh, split, a, b, wall));
        assert_eq!(approximate, slab_field(&mesh, split, a, b, 350.0));
        assert_eq!(bc, boundary_before);
        gate.request();
        assert!(matches!(
            compare(&reference, &approximate, &model, &weights, entries, 10_000),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

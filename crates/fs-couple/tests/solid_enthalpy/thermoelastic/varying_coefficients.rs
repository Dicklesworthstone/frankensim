//! Material-data changes reach the original spatial and inertial owners.
//! Authored numerical sources, not qualification of a measured specimen.
use super::*;
use fs_solid::loads3::transient::TransientConfig3;

#[derive(Clone, Copy)]
struct Coefficients {
    young: [f64; 2],
    poisson: [f64; 2],
    expansion: [f64; 2],
}

const FIRST: Coefficients = Coefficients {
    young: [1000.0, 2.0],
    poisson: [0.15, 0.0015],
    expansion: [0.0003, 0.000006],
};
const SECOND: Coefficients = Coefficients {
    young: [1600.0, -0.5],
    poisson: [0.25, -0.001],
    expansion: [0.0002, -0.000002],
};

fn varying_fixture(
    values: Coefficients,
    capacity: f64,
) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, [value, slope]) in [
        ("young_modulus", Pressure::DIMS, values.young),
        ("poisson_ratio", Dims([0; 6]), values.poisson),
        (
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
            values.expansion,
        ),
        ("density", Density::DIMS, [6.0, 0.0]),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [capacity, 0.0],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, value - 20.0 * slope), (360.0, value + 60.0 * slope)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored variable thermoelastic coefficients".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic variable solid".into(),
            phase: "solid".into(),
            process: "numerical source".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot {
                temperature_k: 280.0,
                specific_enthalpy_j_kg: 0.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 6.0,
            },
            EnthalpyPhaseKnot {
                temperature_k: 360.0,
                specific_enthalpy_j_kg: 80.0 * capacity,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 6.0,
            },
        ],
    )
    .unwrap();
    (card, curve)
}

// Independent Lamé-coordinate derivatives; production uses bulk/shear and
// deviatoric invariants. No source resolver or runtime potential enters here.
fn lame(t: f64, coefficients: Coefficients) -> [f64; 6] {
    let delta = t - 300.0;
    let [e0, et] = coefficients.young;
    let [nu0, nut] = coefficients.poisson;
    let [alpha0, alphat] = coefficients.expansion;
    let e = e0 + et * delta;
    let nu = nu0 + nut * delta;
    let denominator = (1.0 + nu) * (1.0 - 2.0 * nu);
    let lambda = e * nu / denominator;
    let lambda_t =
        (et * nu + e * nut) / denominator + e * nu * nut * (1.0 + 4.0 * nu) / denominator.powi(2);
    let mu = e / (2.0 * (1.0 + nu));
    let mu_t = (et * (1.0 + nu) - e * nut) / (2.0 * (1.0 + nu).powi(2));
    let free = alpha0 * delta + 0.5 * alphat * delta * delta;
    [lambda, mu, lambda_t, mu_t, free, alpha0 + alphat * delta]
}

fn response(t: f64, axial_strain: f64, coefficients: Coefficients) -> [f64; 2] {
    let [lambda, mu, lambda_t, mu_t, free, alpha] = lame(t, coefficients);
    let trace = axial_strain - 3.0 * free;
    let square = (axial_strain - free).powi(2) + 2.0 * free * free;
    let psi = VOLUME * (0.5 * lambda * trace * trace + mu * square);
    let entropy = VOLUME
        * (alpha * (3.0 * lambda + 2.0 * mu) * trace
            - 0.5 * lambda_t * trace * trace
            - mu_t * square);
    [
        (lambda + 2.0 * mu) * axial_strain - (3.0 * lambda + 2.0 * mu) * free,
        psi + t * entropy,
    ]
}

fn compatible_displacement(t: [f64; 2], force: f64, coefficients: [Coefficients; 2]) -> f64 {
    let mut stiffness = 0.0;
    let mut thermal_force = 0.0;
    for i in 0..2 {
        let [lambda, mu, _, _, free, _] = lame(t[i], coefficients[i]);
        stiffness += VOLUME * (lambda + 2.0 * mu);
        thermal_force += VOLUME * (3.0 * lambda + 2.0 * mu) * free;
    }
    (force + thermal_force) / stiffness
}

#[test]
fn g1_affine_expansion_and_poisson_sources_drive_free_and_constrained_finite_heating() {
    with_cx(|cx| {
        for clamped in [false, true] {
            let (card, curve) = varying_fixture(FIRST, 2.0);
            let heat = thermal(&[&curve], &[2.0], &[20.0], 0.0);
            let context = context(300.0);
            let source = SolidEnthalpy::try_new(
                &heat,
                &NODES,
                &TETS,
                vec![source(30, &card, &curve, &context)],
                TetAssemblyBudget::standard(),
                0.05,
                mechanical().mesh_frame,
                cx,
            )
            .unwrap();
            let loads = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let fixed = if clamped {
                (0..12).collect::<Vec<_>>()
            } else {
                FREE.to_vec()
            };
            let body = ThermoelasticSolid::try_new(
                &source,
                &heat,
                &loads,
                &fixed,
                &[[280.0, 360.0]],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let next = body
                .advance(&heat, &(), 0.1, &loads, thermal_controls(), solve(), cx)
                .unwrap();
            let expected_t = if clamped {
                root(280.0, 359.0, |t| {
                    2.0 * (t - 300.0) + response(t, 0.0, FIRST)[1] - 2.0
                })
            } else {
                301.0
            };
            near(
                next.thermal.temperature_k(IDS[0]).unwrap(),
                expected_t,
                2e-8,
            );
            let free = lame(expected_t, FIRST)[4];
            let displacement = next
                .solid
                .equilibrium()
                .solution
                .displacement
                .displacement_m();
            for (node, position) in NODES.iter().enumerate() {
                for axis in 0..3 {
                    near(
                        displacement[node][axis],
                        if clamped { 0.0 } else { free * position[axis] },
                        2e-10,
                    );
                }
            }
            let expected_stress = if clamped {
                response(expected_t, 0.0, FIRST)[0]
            } else {
                0.0
            };
            let actual = next.solid.equilibrium().solution.stress.elements()[0].stress_mandel_pa();
            for stress in &actual[..3] {
                near(*stress, expected_stress, 2e-8);
            }
            near(
                next.thermal.energy_j() - heat.energy_j() + next.report.internal_energy_change_j,
                2.0,
                1e-9,
            );
            near(next.report.combined_residual_j, 0.0, 1e-9);
            assert_eq!(
                next.solid.state().samples()[0].elastic.card_identity(),
                card.content_hash()
            );
            assert_eq!(heat.time_s(), 0.0);
            assert_eq!(body.time_s(), 0.0);
        }
    });
}

// Independent two-variable thermal root, instead of the runtime heat fixed
// point. Each residual runs the specified four explicit conduction substeps.
fn pair_reference(
    old_t: [f64; 2],
    old_u: f64,
    force: f64,
    coefficients: [Coefficients; 2],
) -> ([f64; 2], f64, [f64; 2]) {
    let capacity = [2.0, 3.0];
    let power = [20.0, 1.0];
    let evaluate = |t: [f64; 2]| {
        let u = compatible_displacement(t, force, coefficients);
        let heat = [0, 1].map(|i| {
            let before = response(old_t[i], old_u, coefficients[i]);
            let after = response(t[i], u, coefficients[i]);
            0.5 * VOLUME * (before[0] + after[0]) * (u - old_u) - (after[1] - before[1])
        });
        let mut predicted = old_t;
        for _ in 0..4 {
            let before = predicted;
            predicted = [0, 1].map(|i| {
                before[i]
                    + 0.025 / capacity[i]
                        * (power[i] + heat[i] / 0.1 + 0.2 * (before[1 - i] - before[i]))
            });
        }
        ([0, 1].map(|i| capacity[i] * (t[i] - predicted[i])), u, heat)
    };
    let mut t = old_t;
    for _ in 0..30 {
        let (residual, u, heat) = evaluate(t);
        if residual[0].abs().max(residual[1].abs()) < 1e-12 {
            return (t, u, heat);
        }
        let mut tangent = [[0.0; 2]; 2];
        for col in 0..2 {
            let mut probe = t;
            probe[col] += 1e-4;
            let next = evaluate(probe).0;
            for row in 0..2 {
                tangent[row][col] = (next[row] - residual[row]) / 1e-4;
            }
        }
        let determinant = tangent[0][0] * tangent[1][1] - tangent[0][1] * tangent[1][0];
        t[0] -= (tangent[1][1] * residual[0] - tangent[0][1] * residual[1]) / determinant;
        t[1] -= (tangent[0][0] * residual[1] - tangent[1][0] * residual[0]) / determinant;
        assert!(t.iter().all(|v| (280.0..360.0).contains(v)));
    }
    panic!("independent variable-coefficient thermal root did not converge");
}

#[test]
fn g1_nonuniform_material_coefficient_changes_reach_compatible_stress_heat_and_displacement() {
    with_cx(|cx| {
        let mut outcomes = Vec::new();
        for varying in [false, true] {
            let mut coefficients = [FIRST, SECOND];
            if !varying {
                for values in &mut coefficients {
                    values.poisson[1] = 0.0;
                    values.expansion[1] = 0.0;
                }
            }
            let first = varying_fixture(coefficients[0], 2.0);
            let second = varying_fixture(coefficients[1], 3.0);
            let cards = [first.0, second.0];
            let curves = [first.1, second.1];
            let context = context(300.0);
            let mut heat = thermal(&[&curves[0], &curves[1]], &[2.0, 3.0], &[20.0, 1.0], 0.2);
            let source = SolidEnthalpy::try_new(
                &heat,
                &PAIR_NODES,
                &PAIR_TETS,
                vec![
                    source(30, &cards[0], &curves[0], &context),
                    source(10, &cards[1], &curves[1], &context),
                ],
                TetAssemblyBudget::standard(),
                0.05,
                mechanical().mesh_frame,
                cx,
            )
            .unwrap();
            let fixed: Vec<_> = (0..15).filter(|d| *d != 3).collect();
            let initial_load = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
                .unwrap();
            let mut body = ThermoelasticSolid::try_new(
                &source,
                &heat,
                &initial_load,
                &fixed,
                &[[280.0, 360.0]; 2],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let mut old_t = [300.0; 2];
            let mut old_u = compatible_displacement(old_t, 0.8, coefficients);
            let initial_energy = heat.energy_j()
                + body
                    .potentials()
                    .iter()
                    .map(|p| p.internal_energy_j)
                    .sum::<f64>();
            let mut supplied = 0.0;
            for (index, force) in [1.0, 0.6, 0.8].into_iter().enumerate() {
                let expected = pair_reference(old_t, old_u, force, coefficients);
                let loads = source
                    .geometry()
                    .assemble(&[], VolumeLoad3::None, &[(3, force)], cx)
                    .unwrap();
                let next = body
                    .advance(
                        &heat,
                        &(),
                        (index + 1) as f64 * 0.1,
                        &loads,
                        thermal_controls(),
                        solve(),
                        cx,
                    )
                    .unwrap();
                for i in 0..2 {
                    near(
                        next.thermal.temperature_k(IDS[i]).unwrap(),
                        expected.0[i],
                        3e-8,
                    );
                    near(next.report.cells[i].reversible_heat_j, expected.2[i], 2e-9);
                    near(
                        next.solid.equilibrium().solution.stress.elements()[i].stress_mandel_pa()
                            [0],
                        response(expected.0[i], expected.1, coefficients[i])[0],
                        2e-7,
                    );
                    assert_eq!(
                        next.solid.state().samples()[i].elastic.card_identity(),
                        cards[i].content_hash()
                    );
                    near(
                        next.solid.state().samples()[i].reference_mass_kg,
                        1.0,
                        1e-14,
                    );
                }
                near(
                    next.solid
                        .equilibrium()
                        .solution
                        .displacement
                        .displacement_m()[1][0],
                    expected.1,
                    2e-10,
                );
                near(next.report.combined_residual_j, 0.0, 1e-9);
                supplied += next.report.mechanical_work_j + next.report.external_heat_j;
                near(
                    next.thermal.energy_j()
                        + next
                            .solid
                            .potentials()
                            .iter()
                            .map(|p| p.internal_energy_j)
                            .sum::<f64>(),
                    initial_energy + supplied,
                    4e-9,
                );
                let (next_heat, next_body, _, _) = next.into_parts();
                heat = next_heat;
                body = next_body;
                old_t = expected.0;
                old_u = expected.1;
            }
            assert!((old_t[0] - old_t[1]).abs() > 2.0);
            outcomes.push(old_u);
        }
        assert!((outcomes[0] - outcomes[1]).abs() > 5e-7);
    });
}

#[test]
fn g1_variable_coefficients_enter_actual_inertia_and_the_same_finite_heat_balance() {
    with_cx(|cx| {
        let (card, curve) = varying_fixture(FIRST, 2.0);
        let heat = thermal(&[&curve], &[2.0], &[20.0], 0.0);
        let context = context(300.0);
        let source = SolidEnthalpy::try_new(
            &heat,
            &NODES,
            &TETS,
            vec![source(30, &card, &curve, &context)],
            TetAssemblyBudget::standard(),
            0.05,
            mechanical().mesh_frame,
            cx,
        )
        .unwrap();
        let initial = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 0.8)], cx)
            .unwrap();
        let applied = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 1.2)], cx)
            .unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let equilibrium = ThermoelasticSolid::try_new(
            &source,
            &heat,
            &initial,
            &fixed,
            &[[280.0, 360.0]],
            mechanical(),
            solve(),
            cx,
        )
        .unwrap();
        let u0 = equilibrium
            .equilibrium()
            .solution
            .displacement
            .displacement_m()[1][0];
        let velocity = [[0.0; 3], [0.08, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let body = equilibrium
            .begin_transient(
                &initial,
                &velocity,
                TransientConfig3 {
                    maximum_steps: 1000,
                    maximum_displacement_m: 0.05,
                    absolute_force_tolerance_n: 1e-11,
                    absolute_energy_tolerance_j: 1e-10,
                    relative_energy_tolerance: 1e-8,
                    ..TransientConfig3::default()
                },
                cx,
            )
            .unwrap();
        let dt = 0.02_f64;
        let evaluate = |t: f64| {
            let [lambda0, mu0, _, _, _, _] = lame(300.0, FIRST);
            let [lambda, mu, _, _, free, _] = lame(t, FIRST);
            let k0 = VOLUME * (lambda0 + 2.0 * mu0);
            let k = VOLUME * (lambda + 2.0 * mu);
            let thermal_force = VOLUME * (3.0 * lambda + 2.0 * mu) * free;
            // One free consistent-mass tetrahedral DOF has mass 0.1 kg.
            let u =
                ((0.2 / dt.powi(2) - 0.5 * k0) * u0 + 0.2 / dt * 0.08 + 1.0 + 0.5 * thermal_force)
                    / (0.2 / dt.powi(2) + 0.5 * k);
            let v = 2.0 * (u - u0) / dt - 0.08;
            let before = response(300.0, u0, FIRST);
            let after = response(t, u, FIRST);
            let reversible =
                0.5 * VOLUME * (before[0] + after[0]) * (u - u0) - (after[1] - before[1]);
            (
                2.0 * (t - 300.0) - 20.0 * dt - reversible,
                [u, v, reversible],
            )
        };
        let expected_t = root(280.0, 359.0, |t| evaluate(t).0);
        let [u, v, reversible] = evaluate(expected_t).1;
        let next = body
            .advance(&heat, &(), dt, &applied, thermal_controls(), solve(), cx)
            .unwrap();
        let motion = next.solid.nodal_motion(cx).unwrap();
        near(
            next.thermal.temperature_k(IDS[0]).unwrap(),
            expected_t,
            3e-8,
        );
        near(motion[0][1][0], u, 1e-10);
        near(motion[1][1][0], v, 1e-9);
        near(next.report.reversible_heat_j, reversible, 2e-9);
        near(
            next.report.kinetic_energy_change_j,
            0.05 * (v * v - 0.08_f64.powi(2)),
            1e-10,
        );
        near(next.report.combined_residual_j, 0.0, 1e-9);
        assert_eq!(next.solid.step(), 1);
        assert_eq!(body.step(), 0);
        assert_eq!(heat.time_s(), 0.0);
    });
}

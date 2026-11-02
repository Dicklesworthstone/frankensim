//! Original-source derivatives, with an independent latent-plateau inverse.
use super::*;
use fs_conduction::radiation::pointwise::material::MaterialRadiationSmoothnessError;
use fs_conduction::radiation::pointwise::material::adjoint::MaterialEmissivityGradient;
use fs_conduction::transient::enthalpy::EnthalpyError;
use fs_conduction::{ConductionError, LinearConfig};

const SEED: [f64; 4] = [0.7, -0.3, 0.2, 0.1];
const DIRECTION: [f64; 4] = [0.3, -0.7, 0.9, -0.2];

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 64,
        restart: 4,
    }
}

fn source_card(axis: Option<&str>, values: &[f64]) -> MaterialCard {
    let (value, interpolation) = match axis {
        Some(axis) => {
            let coordinates: Vec<f64> = if axis == "T" {
                vec![300.0, 500.0]
            } else if values.len() == 3 {
                vec![0.0, 0.5, 1.0]
            } else {
                vec![0.0, 1.0]
            };
            assert_eq!(coordinates.len(), values.len());
            (
                PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: if axis == "T" {
                        Temperature::DIMS
                    } else {
                        Dims::NONE
                    },
                    knots: coordinates
                        .into_iter()
                        .zip(values.iter().copied())
                        .collect(),
                    dims: EMISSIVITY_DIMS,
                },
                InterpolationPolicy::LinearInside,
            )
        }
        None => {
            assert_eq!(values.len(), 1);
            (
                PropertyValue::Scalar {
                    value: values[0],
                    dims: EMISSIVITY_DIMS,
                },
                InterpolationPolicy::ConstantWithinValidity,
            )
        }
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value,
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, 500.0)
                .with("melt", 0.0, 1.0)
                .with("lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation,
            observations: Vec::new(),
            provenance: Provenance {
                source: "independent material radiation derivative oracle".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "chart".into(),
            process: "radiation reference".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

// Inverse of a I + b 11^T, derived directly from the three face mass moments.
fn rank_one_inverse(rhs: [f64; 3], a: f64, b: f64) -> [f64; 3] {
    let correction = b * rhs.iter().sum::<f64>() / (a * (a + 3.0 * b));
    rhs.map(|value| value / a - correction)
}

#[test]
fn complete_plateau_action_and_original_phase_controls_match_an_analytic_inverse() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let card = card(Some("melt"));
    let curve = curve(&card, 0.0);
    let context = context();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let old = [250.0, 500.0, 750.0, 600.0];
    let dt = 0.2;
    for bath in [300.0_f64, 450.0] {
        // dT/dh is exactly zero, but the full phase-emissivity Jacobian is not.
        // M_face = A/12 (I + 11^T), and all four reference masses are 1 kg.
        let c = STEFAN_BOLTZMANN_W_M2_K4 * (400.0_f64.powi(4) - bath.powi(4));
        let beta = dt * c * 0.6 / LATENT * AREA / 12.0;
        let expected = rank_one_inverse(
            std::array::from_fn(|i| old[i] - dt * c * 0.2 * AREA / 3.0),
            1.0 + beta,
            beta,
        );
        let lambda = rank_one_inverse([SEED[0], SEED[1], SEED[2]], 1.0 + beta, beta);
        with_cx(|cx| {
            let radiation = MaterialPointwiseAmbientRadiation::new(
                cx,
                &mesh,
                &boundary,
                &[patch(&card, &curve, &context, bath)],
            )
            .unwrap();
            let engine = HeterogeneousEnthalpyBackwardEuler::new(
                cx,
                &mesh,
                &[ReferenceEnthalpyMaterial {
                    curve: &curve,
                    reference_density_kg_m3: RHO,
                }],
                &[0],
                budget(),
            )
            .unwrap();
            let accepted = engine
                .advance_with_material_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &old,
                    dt,
                    config(),
                    &radiation,
                )
                .unwrap();
            for (actual, expected) in accepted.specific_enthalpy_j_kg[..3].iter().zip(expected) {
                close(*actual, expected, 3e-9);
            }
            close(accepted.specific_enthalpy_j_kg[3], old[3], 1e-12);
            let claim = accepted.pointwise_radiation.as_ref().unwrap().patches[0].material_states
                [0]
            .emissivity
            .receipt()
            .selected;
            let bound = engine
                .linearize_accepted_with_material_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &old,
                    dt,
                    config(),
                    &radiation,
                    accepted,
                    1e-8,
                    radiation.adjoint_workspace_entries().unwrap(),
                )
                .unwrap();
            let transport = bound.transport();
            assert_eq!(
                transport.temperature_pullback(cx, &SEED).unwrap(),
                vec![0.0; 4]
            );
            let action = transport.apply_jacobian(cx, &DIRECTION).unwrap();
            let transpose = transport.apply_jacobian_transpose(cx, &SEED).unwrap();
            for i in 0..3 {
                close(
                    action[i],
                    DIRECTION[i] + beta * (DIRECTION[i] + DIRECTION[..3].iter().sum::<f64>()),
                    3e-13,
                );
                close(
                    transpose[i],
                    SEED[i] + beta * (SEED[i] + SEED[..3].iter().sum::<f64>()),
                    3e-13,
                );
            }
            close(action[3], DIRECTION[3], 1e-15);
            close(transpose[3], SEED[3], 1e-15);
            close(dot(&SEED, &action), dot(&transpose, &DIRECTION), 3e-13);
            let gradient = bound
                .pullback_with_material_ambient_controls(cx, &SEED, linear())
                .unwrap();
            for (actual, expected) in gradient.thermal.transport.previous_specific_enthalpy[..3]
                .iter()
                .zip(lambda)
            {
                close(*actual, expected, 2e-12);
            }
            close(gradient.thermal.transport.adjoint[3], SEED[3], 2e-13);
            let patch = &gradient.radiation.patches[0];
            assert_eq!(patch.material_identity, card.content_hash());
            assert_eq!(patch.claim, claim);
            assert_eq!(patch.region, "surface");
            let MaterialEmissivityGradient::Curve {
                axis,
                knot_abscissae,
                knot_value_derivatives,
            } = &patch.emissivity
            else {
                panic!("original curve controls required");
            };
            assert_eq!(axis, "melt");
            assert_eq!(knot_abscissae, &[0.0, 1.0]);
            let mut knots = [0.0; 2];
            let mut reservoir = 0.0;
            for (i, &dual) in lambda.iter().enumerate() {
                for (j, &h) in expected.iter().enumerate() {
                    let mass = AREA / 12.0 * if i == j { 2.0 } else { 1.0 };
                    knots[0] -= dt * c * dual * mass * (1.0 - h / LATENT);
                    knots[1] -= dt * c * dual * mass * h / LATENT;
                    reservoir += dt
                        * 4.0
                        * STEFAN_BOLTZMANN_W_M2_K4
                        * bath.powi(3)
                        * dual
                        * mass
                        * (0.2 + 0.6 * h / LATENT);
                }
            }
            for (actual, expected) in knot_value_derivatives.iter().zip(knots) {
                close(*actual, expected, 3e-11);
            }
            close(patch.reservoir_temperature, reservoir, 3e-13);
            let density = -lambda
                .iter()
                .enumerate()
                .map(|(i, dual)| dual * (expected[i] - old[i]))
                .sum::<f64>()
                / RHO;
            close(
                gradient.thermal.material_reference_density[0],
                density,
                3e-11,
            );
            assert!(
                transport.pointwise_parameter_pullback(cx, &SEED).is_err(),
                "a material tangent must not expose frozen snapshot controls"
            );
        });
    }
    assert_eq!(old, [250.0, 500.0, 750.0, 600.0]);
}

// Rebuild the actual source and its identity-matched caloric owner for each
// perturbed control; this never differentiates a frozen nodal epsilon vector.
fn solve_controls(
    axis: Option<&str>,
    values: &[f64],
    bath: f64,
    history: f64,
    density: f64,
    reverse: bool,
) -> (f64, Vec<f64>) {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let card = source_card(axis, values);
    let curve = curve(&card, 0.0);
    let context = context();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let old: [f64; 4] =
        std::array::from_fn(|i| [-60.0, -45.0, -30.0, -20.0][i] + history * DIRECTION[i]);
    with_cx(|cx| {
        let radiation = MaterialPointwiseAmbientRadiation::new(
            cx,
            &mesh,
            &boundary,
            &[patch(&card, &curve, &context, bath)],
        )
        .unwrap();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[ReferenceEnthalpyMaterial {
                curve: &curve,
                reference_density_kg_m3: density,
            }],
            &[0],
            budget(),
        )
        .unwrap();
        let accepted = engine
            .advance_with_material_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                0.02,
                config(),
                &radiation,
            )
            .unwrap();
        let objective = dot(&SEED, &accepted.specific_enthalpy_j_kg);
        if !reverse {
            return (objective, Vec::new());
        }
        let bound = engine
            .linearize_accepted_with_material_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                0.02,
                config(),
                &radiation,
                accepted,
                1e-8,
                radiation.adjoint_workspace_entries().unwrap(),
            )
            .unwrap();
        let transport = bound.transport();
        let action = transport.apply_jacobian(cx, &DIRECTION).unwrap();
        let transpose = transport.apply_jacobian_transpose(cx, &SEED).unwrap();
        close(dot(&SEED, &action), dot(&transpose, &DIRECTION), 3e-13);
        let gradient = bound
            .pullback_with_material_ambient_controls(cx, &SEED, linear())
            .unwrap();
        let patch = &gradient.radiation.patches[0];
        let mut derivatives = match &patch.emissivity {
            MaterialEmissivityGradient::Scalar { derivative } => vec![*derivative],
            MaterialEmissivityGradient::Curve {
                axis,
                knot_abscissae,
                knot_value_derivatives,
            } => {
                assert_eq!(axis, "T");
                assert_eq!(knot_abscissae, &[300.0, 500.0]);
                knot_value_derivatives.clone()
            }
        };
        derivatives.extend([
            patch.reservoir_temperature,
            dot(
                &gradient.thermal.transport.previous_specific_enthalpy,
                &DIRECTION,
            ),
            gradient.thermal.material_reference_density[0],
        ]);
        (objective, derivatives)
    })
}

#[test]
fn sensible_scalar_and_temperature_curve_controls_match_resolved_source_finite_differences() {
    for (axis, values) in [(None, vec![0.5]), (Some("T"), vec![0.2, 0.8])] {
        let bath = 310.0;
        let (_, derivatives) = solve_controls(axis, &values, bath, 0.0, RHO, true);
        for (control, &expected) in derivatives.iter().enumerate() {
            let delta = if control < values.len() { 1e-4 } else { 0.02 };
            let objective = |offset: f64| {
                let mut changed = values.clone();
                let mut bath = bath;
                let mut history = 0.0;
                let mut density = RHO;
                if control < values.len() {
                    changed[control] += offset;
                } else if control == values.len() {
                    bath += offset;
                } else if control == values.len() + 1 {
                    history += offset;
                } else {
                    density += offset;
                }
                solve_controls(axis, &changed, bath, history, density, false).0
            };
            let difference = (objective(-2.0 * delta) - 8.0 * objective(-delta)
                + 8.0 * objective(delta)
                - objective(2.0 * delta))
                / (12.0 * delta);
            close(expected, difference, 2e-7 * expected.abs().max(1.0));
        }
    }
}

#[test]
fn moving_source_and_chart_corners_refuse_but_constant_coordinates_admit() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let context = context();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let card = source_card(Some("melt"), &[0.2, 0.6, 0.8]);
    let curve = curve(&card, 0.0);
    for (h, bath, expected) in [
        (500.0, 400.0, "source"),
        (0.0, 400.0, "chart"),
        (-50.0, 350.0, "constant"),
    ] {
        with_cx(|cx| {
            let radiation = MaterialPointwiseAmbientRadiation::new(
                cx,
                &mesh,
                &boundary,
                &[patch(&card, &curve, &context, bath)],
            )
            .unwrap();
            let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
            let old = [h; 4];
            let accepted = engine
                .advance_with_material_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &old,
                    0.1,
                    config(),
                    &radiation,
                )
                .unwrap();
            assert_eq!(accepted.specific_enthalpy_j_kg, old);
            let bound = engine
                .linearize_accepted_with_material_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &old,
                    0.1,
                    config(),
                    &radiation,
                    accepted,
                    1e-8,
                    radiation.adjoint_workspace_entries().unwrap(),
                );
            if expected == "constant" {
                // f=0 on an open solid h interval. The original fraction
                // endpoint is therefore a smooth constant composition.
                let gradient = bound.unwrap()
                    .pullback_with_material_ambient_controls(cx, &SEED, linear()).unwrap();
                let MaterialEmissivityGradient::Curve { knot_value_derivatives, .. } =
                    &gradient.radiation.patches[0].emissivity else { panic!("curve source") };
                assert_eq!(knot_value_derivatives, &[0.0; 3]);
                return;
            }
            let error = bound.unwrap_err();
            let EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
                ConductionError::RadiationMaterialSmoothness(error),
            )) = error
            else {
                panic!("expected original source/phase smoothness refusal, got {error:?}");
            };
            assert!(
                match expected {
                    "source" =>
                        matches!(error, MaterialRadiationSmoothnessError::SourceKink { .. }),
                    "chart" => matches!(error, MaterialRadiationSmoothnessError::ChartKink { .. }),
                    _ => matches!(
                        error,
                        MaterialRadiationSmoothnessError::SourceEndpoint { .. }
                    ),
                },
                "{error:?}"
            );
        });
    }
}

#[test]
fn bounded_binding_rechecks_receipts_and_never_relabels_a_cancelled_derivative() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let card = card(Some("melt"));
    let curve = curve(&card, 0.0);
    let context = context();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let old = [500.0; 4];
    with_cx(|cx| {
        let radiation = MaterialPointwiseAmbientRadiation::new(
            cx,
            &mesh,
            &boundary,
            &[patch(&card, &curve, &context, 300.0)],
        )
        .unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let accepted = engine
            .advance_with_material_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                0.1,
                config(),
                &radiation,
            )
            .unwrap();
        let entries = radiation.adjoint_workspace_entries().unwrap();
        assert!(matches!(
            engine.linearize_accepted_with_material_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                0.1,
                config(),
                &radiation,
                accepted.clone(),
                1e-8,
                entries - 1,
            ),
            Err(EnthalpyAdjointError::InvalidInput(_))
        ));
        let mut changed = accepted.clone();
        changed.pointwise_radiation.as_mut().unwrap().patches[0].material_states[0].vertex = 3;
        assert!(matches!(
            engine.linearize_accepted_with_material_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                0.1,
                config(),
                &radiation,
                changed,
                1e-8,
                entries,
            ),
            Err(EnthalpyAdjointError::InvalidInput(_))
        ));
        with_cancelled_cx(|cancelled| {
            assert!(
                engine
                    .linearize_accepted_with_material_ambient_radiation(
                        cancelled,
                        problem,
                        None,
                        &old,
                        0.1,
                        config(),
                        &radiation,
                        accepted.clone(),
                        1e-8,
                        entries,
                    )
                    .is_err()
            );
        });
        let bound = engine
            .linearize_accepted_with_material_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                0.1,
                config(),
                &radiation,
                accepted,
                1e-8,
                entries,
            )
            .unwrap();
        assert!(
            bound
                .pullback_with_material_ambient_controls(
                    cx,
                    &SEED,
                    LinearConfig {
                        max_iterations: 0,
                        ..linear()
                    }
                )
                .is_err()
        );
        let gradient = bound
            .pullback_with_material_ambient_controls(cx, &SEED, linear())
            .unwrap();
        assert!(gradient.thermal.relative_residual <= linear().tolerance);
        with_cancelled_cx(|cancelled| {
            assert!(
                bound
                    .pullback_with_material_ambient_controls(cancelled, &SEED, linear())
                    .is_err()
            );
        });
        assert_eq!(bound.primal().temperature, vec![400.0; 4]);
    });
    assert_eq!(old, [500.0; 4]);
}

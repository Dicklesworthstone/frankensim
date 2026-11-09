//! Actual SDIRK stage controls against independent latent equations and solves.
use super::*;
use fs_conduction::transient::enthalpy::adjoint::sdirk::{
    EnthalpySdirk2ParameterBudget, EnthalpySdirk2RadiationGradient,
};
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA as G};

const OLD_H: [f64; 4] = [250.0, 500.0, 750.0, 600.0];
const DT: f64 = 0.2;

fn difference(mut value: impl FnMut(f64) -> f64, delta: f64) -> f64 {
    (value(-2.0 * delta) - 8.0 * value(-delta) + 8.0 * value(delta) - value(2.0 * delta))
        / (12.0 * delta)
}

// Each face node has mass one. Exact triangle moments give a I + b 11^T
// on the latent plateau. This solves the two Butcher stages independently of
// the production Jacobian, source queries, quadrature and parameter pullback.
fn latent_stages(values: [[f64; 3]; 2]) -> [[f64; 4]; 2] {
    let load = |h: &[f64; 4], p: [f64; 3]| {
        let epsilon = [h[0], h[1], h[2]].map(|h| p[0] + (p[1] - p[0]) * h / LATENT);
        let coefficient =
            AREA / 12.0 * STEFAN_BOLTZMANN_W_M2_K4 * (400.0_f64.powi(4) - p[2].powi(4));
        std::array::from_fn::<_, 4, _>(|i| {
            if i == 3 {
                0.0
            } else {
                coefficient * (epsilon.iter().sum::<f64>() + epsilon[i])
            }
        })
    };
    let stage = |base: [f64; 4], p: [f64; 3]| {
        let coefficient =
            G * DT * AREA / 12.0 * STEFAN_BOLTZMANN_W_M2_K4 * (400.0_f64.powi(4) - p[2].powi(4));
        let beta = coefficient * (p[1] - p[0]) / LATENT;
        let face = rank_one_inverse(
            [base[0], base[1], base[2]].map(|h| h - 4.0 * coefficient * p[0]),
            1.0 + beta,
            beta,
        );
        [face[0], face[1], face[2], base[3]]
    };
    let first = stage(OLD_H, values[0]);
    let rate = load(&first, values[0]);
    let base = std::array::from_fn(|i| OLD_H[i] - (1.0 - G) * DT * rate[i]);
    [first, stage(base, values[1])]
}

#[test]
fn original_material_ordinates_and_baths_follow_each_latent_stage_without_extra_weights() {
    with_cx(|cx| {
        let mesh = mesh();
        let boundary = boundary(&mesh);
        let card = source_card(Some("melt"), &[0.2, 0.8]);
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
        let radiation = [300.0, 450.0].map(|bath| {
            MaterialPointwiseAmbientRadiation::new(
                cx,
                &mesh,
                &boundary,
                &[patch(&card, &curve, &context, bath)],
            )
            .unwrap()
        });
        let entries = radiation
            .iter()
            .map(|r| r.adjoint_workspace_entries().unwrap())
            .sum();
        let inputs = radiation.each_ref().map(|radiation| EnthalpySdirk2Input {
            problem,
            interfaces: None,
            contacts: None,
            radiation: EnthalpyRadiation::MaterialAmbient(radiation),
        });
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let bound = engine
            .linearize_sdirk2(cx, inputs, &OLD_H, DT, config(), entries)
            .unwrap();
        let caps = EnthalpySdirk2ParameterBudget {
            max_contact_face_pairs_per_stage: 0,
            max_radiation_entries: entries,
        };
        let ordinary = bound
            .pullback_stages(cx, [&SEED, &DIRECTION], linear())
            .unwrap();
        let gradient = bound
            .pullback_stages_with_parameters(cx, [&SEED, &DIRECTION], linear(), caps)
            .unwrap();
        assert_eq!(gradient.transport.stage_adjoint, ordinary.stage_adjoint);
        assert_eq!(gradient.transport.iterations, ordinary.iterations);
        assert_eq!(
            gradient.transport.previous_specific_enthalpy,
            ordinary.previous_specific_enthalpy
        );
        let controls = [[0.2, 0.8, 300.0], [0.2, 0.8, 450.0]];
        let expected = latent_stages(controls);
        for (actual, expected) in [&bound.primal().first_stage, &bound.primal().second_stage]
            .iter()
            .zip(expected)
        {
            for (&h, expected) in actual.specific_enthalpy_j_kg.iter().zip(expected) {
                close(h, expected, 3e-8);
            }
            assert_eq!(actual.temperature, [400.0; 4]);
        }
        for stage in 0..2 {
            assert!(gradient.stages[stage].resistance_multipliers.is_empty());
            let Some(EnthalpySdirk2RadiationGradient::MaterialAmbient(radiation)) =
                &gradient.stages[stage].radiation
            else {
                panic!("original material controls required")
            };
            let patch = &radiation.patches[0];
            assert_eq!(patch.material_identity, card.content_hash());
            assert_eq!(patch.region, "surface");
            let MaterialEmissivityGradient::Curve {
                axis,
                knot_abscissae,
                knot_value_derivatives,
            } = &patch.emissivity
            else {
                panic!("original phase curve ordinates required")
            };
            assert_eq!(axis, "melt");
            assert_eq!(knot_abscissae, &[0.0, 1.0]);
            let actual = [
                knot_value_derivatives[0],
                knot_value_derivatives[1],
                patch.reservoir_temperature,
            ];
            for parameter in 0..3 {
                let expected = difference(
                    |offset| {
                        let mut values = controls;
                        values[stage][parameter] += offset;
                        let h = latent_stages(values);
                        dot(&h[0], &SEED) + dot(&h[1], &DIRECTION)
                    },
                    if parameter == 2 { 0.01 } else { 1e-4 },
                );
                close(actual[parameter], expected, 2e-7 * expected.abs().max(1.0));
            }
        }
        let multipliers = [
            &gradient.transport.stage_adjoint[0][..],
            &gradient.transport.stage_adjoint[1][..],
        ];
        assert_eq!(
            bound.parameter_pullback(cx, multipliers, caps).unwrap(),
            gradient.stages
        );
        assert!(
            bound
                .parameter_pullback(
                    cx,
                    multipliers,
                    EnthalpySdirk2ParameterBudget {
                        max_radiation_entries: entries - 1,
                        ..caps
                    }
                )
                .is_err()
        );
        assert!(
            bound
                .parameter_pullback(cx, [&[0.0; 3], &[0.0; 4]], caps)
                .is_err()
        );
        assert!(
            bound
                .parameter_pullback(cx, [&[0.0; 4], &[f64::NAN; 4]], caps)
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(
                bound
                    .parameter_pullback(cancelled, multipliers, caps)
                    .is_err()
            )
        });
        assert_eq!(
            bound
                .pullback_stages(cx, [&SEED, &DIRECTION], linear())
                .unwrap()
                .stage_adjoint,
            ordinary.stage_adjoint
        );
    });
}

fn masked_boundary(mesh: &ConductionMesh, axis: usize, temperature: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "surface",
            |face| face.centroid[2] == 0.0,
            ThermalBc::adiabatic(),
        )
        .unwrap()
        .region(
            "bath",
            |face| face.centroid[axis] == 0.0,
            ThermalBc::dirichlet(temperature).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

#[test]
fn fixed_pointwise_controls_use_the_actual_changing_prescribed_rows() {
    with_cx(|cx| {
        let mesh = mesh();
        let boundaries = [
            masked_boundary(&mesh, 0, 340.0),
            masked_boundary(&mesh, 1, 360.0),
        ];
        let card = source_card(None, &[0.5]);
        let curve = curve(&card, 0.0);
        let context = context();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let old = [-50.0, -40.0, -30.0, -20.0];
        let dt = 0.005;
        let controls = [[0.4, 300.0], [0.7, 450.0]];
        let declarations = |values: [[f64; 2]; 2]| {
            values.map(|p| {
                let source = source_card(None, &[p[0]]);
                [AmbientRadiationPatch::new(
                    "surface",
                    SurfaceEmissivity::from_condition(
                        "surface",
                        binding(&source, &context),
                        MaterialPropertySelection::SingleClaimOnly,
                    )
                    .unwrap()
                    .0,
                    p[1],
                )
                .unwrap()]
            })
        };
        let patches = declarations(controls);
        let radiation = [0, 1].map(|i| {
            PointwiseAmbientRadiation::new(cx, &mesh, &boundaries[i], &patches[i]).unwrap()
        });
        let inputs = [0, 1].map(|i| EnthalpySdirk2Input {
            problem: ConductionProblem {
                mesh: &mesh,
                boundary: &boundaries[i],
                material: &material,
                source: &source,
                element_materials: None,
            },
            interfaces: None,
            contacts: None,
            radiation: EnthalpyRadiation::PointwiseAmbient(&radiation[i]),
        });
        let caps = EnthalpySdirk2ParameterBudget {
            max_contact_face_pairs_per_stage: 0,
            max_radiation_entries: 100_000,
        };
        let bound = engine
            .linearize_sdirk2(cx, inputs, &old, dt, config(), caps.max_radiation_entries)
            .unwrap();
        let result = bound
            .pullback_stages_with_parameters(cx, [&SEED, &DIRECTION], linear(), caps)
            .unwrap();
        let value = |values| {
            let patches = declarations(values);
            let radiation = [0, 1].map(|i| {
                PointwiseAmbientRadiation::new(cx, &mesh, &boundaries[i], &patches[i]).unwrap()
            });
            let inputs = [0, 1].map(|i| EnthalpySdirk2Input {
                radiation: EnthalpyRadiation::PointwiseAmbient(&radiation[i]),
                ..inputs[i]
            });
            let step = engine
                .advance_sdirk2(cx, inputs, &old, dt, config())
                .unwrap();
            dot(&step.first_stage.specific_enthalpy_j_kg, &SEED)
                + dot(&step.second_stage.specific_enthalpy_j_kg, &DIRECTION)
        };
        for stage in 0..2 {
            let Some(EnthalpySdirk2RadiationGradient::Pointwise(gradient)) =
                &result.stages[stage].radiation
            else {
                panic!("fixed pointwise controls required")
            };
            let actual = [gradient.emissivities[0], gradient.reservoir_temperatures[0]];
            for parameter in 0..2 {
                let expected = difference(
                    |offset| {
                        let mut values = controls;
                        values[stage][parameter] += offset;
                        value(values)
                    },
                    if parameter == 0 { 1e-4 } else { 0.01 },
                );
                close(actual[parameter], expected, 3e-7 * expected.abs().max(1.0));
            }
        }
        let mut changed = result.transport.stage_adjoint.clone();
        for stage in 0..2 {
            for &(vertex, _) in boundaries[stage].dirichlet() {
                changed[stage][vertex] += 1000.0;
            }
        }
        assert_eq!(
            bound
                .parameter_pullback(cx, [&changed[0], &changed[1]], caps)
                .unwrap(),
            result.stages,
            "prescribed-row multipliers must not enter radiation controls"
        );
    });
}

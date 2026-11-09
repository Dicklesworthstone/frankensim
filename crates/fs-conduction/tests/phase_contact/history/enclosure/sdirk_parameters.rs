//! Separate physical SDIRK controls retain coupled contact and reflected heat.
use super::*;
use fs_conduction::radiation::pointwise::material::adjoint::MaterialEmissivityGradient;
use fs_conduction::transient::enthalpy::adjoint::sdirk::{
    EnthalpySdirk2ParameterBudget, EnthalpySdirk2RadiationGradient,
};
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA as G};

const FIRST: [f64; 8] = [0.03, -0.01, 0.04, -0.02, 0.02, -0.03, 0.01, 0.05];
const LAST: [f64; 8] = [0.7, -0.4, 0.3, 0.2, -0.6, 0.9, -0.2, 0.5];
const DT: f64 = 0.0005;

fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 64,
        restart: 8,
    }
}

// Controls are contact scale, two original emissivity ordinates on each body,
// and open-bath kelvin. Direct P1 moments plus two scalar radiosity equations
// use the real phase inventory, with no production source/tangent calls.
fn outward(h: &[f64; 8], p: [f64; 6], mutual: f64) -> [f64; 8] {
    let fa = h[..3].iter().sum::<f64>() / 300.0;
    let fb = h[4..7].iter().sum::<f64>() / 300.0;
    let contact = AREA * (400.0 - 500.0) / (p[0] * (1.0 + 2.0 * fa + fb));
    let epsilon = [0, 1].map(|body| {
        RADIATING[body]
            .map(|node| p[1 + 2 * body] + (p[2 + 2 * body] - p[1 + 2 * body]) * h[node] / 100.0)
    });
    let means = epsilon.map(|e| e.iter().sum::<f64>() / 3.0);
    let black = [400.0_f64, 500.0].map(|t| STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4));
    let bath = STEFAN_BOLTZMANN_W_M2_K4 * p[5].powi(4);
    let rhs: [f64; 2] =
        std::array::from_fn(|i| means[i] * black[i] + (1.0 - means[i]) * (1.0 - mutual) * bath);
    let a = (1.0 - means[0]) * mutual;
    let b = (1.0 - means[1]) * mutual;
    let radiosity = [
        (rhs[0] + a * rhs[1]) / (1.0 - a * b),
        (rhs[1] + b * rhs[0]) / (1.0 - a * b),
    ];
    let mut load = [0.0; 8];
    for local in 0..3 {
        load[local] = contact / 3.0;
        load[4 + local] = -contact / 3.0;
    }
    for body in 0..2 {
        let irradiation = mutual * radiosity[1 - body] + (1.0 - mutual) * bath;
        for (local, node) in RADIATING[body].into_iter().enumerate() {
            load[node] += AREA / 12.0
                * (epsilon[body].iter().sum::<f64>() + epsilon[body][local])
                * (black[body] - irradiation);
        }
    }
    load
}

fn stages(parameters: [[f64; 6]; 2], mutual: f64) -> [[f64; 8]; 2] {
    let stage = |base: [f64; 8], parameters| {
        let mut h = base;
        for _ in 0..128 {
            let load = outward(&h, parameters, mutual);
            h = std::array::from_fn(|i| {
                base[i] - G * DT * load[i] * if i < 4 { 24.0 } else { 12.0 }
            });
        }
        assert!(h.iter().all(|h| *h > 0.0 && *h < 100.0));
        h
    };
    let first = stage(OLD, parameters[0]);
    let load = outward(&first, parameters[0], mutual);
    let base = std::array::from_fn(|i| {
        OLD[i] - (1.0 - G) * DT * load[i] * if i < 4 { 24.0 } else { 12.0 }
    });
    [first, stage(base, parameters[1])]
}

#[test]
fn native_contact_and_material_reflection_controls_match_independent_stage_equations() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let history = history(cx, &mesh);
        let boundary = history.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let contact = contacts(
            cx,
            &mesh,
            &boundary,
            &cards,
            &curves,
            &source,
            &context,
            sides(&mesh, 4, 2, 0.0),
            0,
            false,
        );
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let zero = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &conductivity,
            element_materials: None,
            source: &zero,
        };
        let engine = engine(cx, &mesh, &curves, 4);
        let closed = closed(&mesh, &cards, &context);
        let open = open(cx, &mesh, &cards, &context);
        let declarations = patches(&cards, &curves, &context);
        for mutual in [1.0, 0.4] {
            let bases = [350.0, 420.0].map(|bath| {
                if mutual == 1.0 {
                    PointwiseEnclosureRadiation::new(
                        cx,
                        &mesh,
                        &boundary,
                        &closed,
                        PointwiseEnclosureTolerance::default(),
                    )
                    .unwrap()
                } else {
                    PointwiseEnclosureRadiation::new_open(
                        cx,
                        &mesh,
                        &boundary,
                        &open,
                        bath,
                        PointwiseEnclosureTolerance::default(),
                    )
                    .unwrap()
                }
            });
            let radiation = bases.each_ref().map(|base| {
                MaterialPointwiseEnclosureRadiation::new(cx, base, &declarations).unwrap()
            });
            let entries = radiation
                .iter()
                .map(|r| r.adjoint_workspace_entries().unwrap())
                .sum();
            let inputs = radiation.each_ref().map(|r| EnthalpySdirk2Input {
                problem,
                interfaces: None,
                contacts: Some((&contact, config())),
                radiation: EnthalpyRadiation::MaterialEnclosure(r),
            });
            let bound = engine
                .linearize_sdirk2(cx, inputs, &OLD, DT, config().step, entries)
                .unwrap();
            let caps = EnthalpySdirk2ParameterBudget {
                max_contact_face_pairs_per_stage: 1,
                max_radiation_entries: entries,
            };
            let gradient = bound
                .pullback_stages_with_parameters(cx, [&FIRST, &LAST], linear(), caps)
                .unwrap();
            let ordinary = bound
                .pullback_stages(cx, [&FIRST, &LAST], linear())
                .unwrap();
            assert_eq!(gradient.transport.stage_adjoint, ordinary.stage_adjoint);
            assert_eq!(gradient.transport.iterations, ordinary.iterations);
            let parameters = [
                [1.0, 0.2, 0.8, 0.2, 0.8, 350.0],
                [1.0, 0.2, 0.8, 0.2, 0.8, 420.0],
            ];
            let expected = stages(parameters, mutual);
            for (stage, expected) in [&bound.primal().first_stage, &bound.primal().second_stage]
                .iter()
                .zip(expected)
            {
                for (&actual, expected) in stage.specific_enthalpy_j_kg.iter().zip(expected) {
                    near(actual, expected, 4e-8);
                }
            }
            for stage in 0..2 {
                let controls = &gradient.stages[stage];
                assert_eq!(controls.resistance_multipliers.len(), 1);
                let resistance = &controls.resistance_multipliers[0];
                assert_eq!(resistance.interface, "joint");
                assert_eq!(resistance.card_identity, source.content_hash());
                assert_eq!(resistance.face_pairs, 1);
                let mut actual = vec![resistance.derivative];
                let Some(EnthalpySdirk2RadiationGradient::MaterialEnclosure(radiation)) =
                    &controls.radiation
                else {
                    panic!("original enclosure controls required")
                };
                for (body, patch) in radiation.patches.iter().enumerate() {
                    assert_eq!(patch.region, REGIONS[body]);
                    assert_eq!(patch.material_identity, cards[body].content_hash());
                    let MaterialEmissivityGradient::Curve {
                        axis,
                        knot_abscissae,
                        knot_value_derivatives,
                    } = &patch.emissivity
                    else {
                        panic!("source ordinates must survive both stages")
                    };
                    assert_eq!(axis, "melt");
                    assert_eq!(knot_abscissae, &[0.0, 1.0]);
                    actual.extend(knot_value_derivatives);
                }
                assert_eq!(radiation.reservoir_temperature.is_none(), mutual == 1.0);
                actual.push(radiation.reservoir_temperature.unwrap_or(0.0));
                for (parameter, actual) in actual.into_iter().enumerate() {
                    let delta = if parameter == 5 { 0.01 } else { 0.0001 };
                    let value = |offset: f64| {
                        let mut p = parameters;
                        p[stage][parameter] += offset;
                        let h = stages(p, mutual);
                        dot(&h[0], &FIRST) + dot(&h[1], &LAST)
                    };
                    let expected = (value(-2.0 * delta) - 8.0 * value(-delta) + 8.0 * value(delta)
                        - value(2.0 * delta))
                        / (12.0 * delta);
                    near(actual, expected, 3e-7 * expected.abs().max(1.0));
                }
            }
            assert!(
                bound
                    .pullback_with_parameters(
                        cx,
                        &LAST,
                        linear(),
                        EnthalpySdirk2ParameterBudget {
                            max_contact_face_pairs_per_stage: 0,
                            ..caps
                        }
                    )
                    .is_err()
            );
            assert!(
                bound
                    .pullback_with_parameters(
                        cx,
                        &LAST,
                        linear(),
                        EnthalpySdirk2ParameterBudget {
                            max_radiation_entries: entries - 1,
                            ..caps
                        }
                    )
                    .is_err()
            );
            if mutual == 1.0 {
                let mass: [f64; 8] =
                    std::array::from_fn(|i| if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 });
                let conserved = bound
                    .pullback_with_parameters(cx, &mass, linear(), caps)
                    .unwrap();
                for stage in conserved.stages {
                    near(stage.resistance_multipliers[0].derivative, 0.0, 2e-10);
                    let Some(EnthalpySdirk2RadiationGradient::MaterialEnclosure(r)) =
                        stage.radiation
                    else {
                        unreachable!()
                    };
                    for patch in r.patches {
                        let MaterialEmissivityGradient::Curve {
                            knot_value_derivatives,
                            ..
                        } = patch.emissivity
                        else {
                            unreachable!()
                        };
                        for derivative in knot_value_derivatives {
                            near(derivative, 0.0, 2e-10);
                        }
                    }
                }
            }
        }
    });
}

#[test]
fn fixed_contact_parameters_survive_dropped_input_owners_and_keep_stage_weights() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let authority = InterfaceResistance::from_card(
            "joint",
            &source,
            &context,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let zero = ScalarField::Uniform(0.0);
        let engine = engine(cx, &mesh, &curves, 4);
        let resistances = [2.1, 3.2];
        let bound = {
            let interfaces = resistances.map(|r| {
                ThermalInterfaces::new(
                    &mesh,
                    &boundary,
                    vec![
                        InterfaceSurface::new(
                            "joint",
                            ThermalInterfaces::coincident_face_pairs(&mesh).unwrap(),
                            authority
                                .with_measured_value(
                                    r,
                                    UncertaintyModel::Unstated,
                                    "manufactured fixed-stage resistance control",
                                )
                                .unwrap(),
                        )
                        .unwrap(),
                    ],
                )
                .unwrap()
            });
            let inputs = interfaces.each_ref().map(|interfaces| EnthalpySdirk2Input {
                problem: ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &conductivity,
                    element_materials: None,
                    source: &zero,
                },
                interfaces: Some(interfaces),
                contacts: None,
                radiation: EnthalpyRadiation::None,
            });
            engine
                .linearize_sdirk2(cx, inputs, &OLD, DT, config().step, 0)
                .unwrap()
        };
        assert!(bound.primal().first_stage.physical_contacts.is_none());
        let caps = EnthalpySdirk2ParameterBudget {
            max_contact_face_pairs_per_stage: 1,
            max_radiation_entries: 0,
        };
        let gradient = bound
            .pullback_stages_with_parameters(cx, [&FIRST, &LAST], linear(), caps)
            .unwrap();
        for stage in 0..2 {
            let rate: [f64; 8] = std::array::from_fn(|i| {
                if i % 4 == 3 {
                    0.0
                } else {
                    AREA * 100.0 / (3.0 * resistances[stage]) * if i < 4 { 24.0 } else { -12.0 }
                }
            });
            let expected = if stage == 0 {
                -DT * (G * dot(&rate, &FIRST) + (1.0 - G) * dot(&rate, &LAST))
            } else {
                -G * DT * dot(&rate, &LAST)
            };
            assert!(gradient.stages[stage].radiation.is_none());
            let control = &gradient.stages[stage].resistance_multipliers[0];
            assert_eq!(control.card_identity, source.content_hash());
            near(control.derivative, expected, 2e-11);
        }
    });
}

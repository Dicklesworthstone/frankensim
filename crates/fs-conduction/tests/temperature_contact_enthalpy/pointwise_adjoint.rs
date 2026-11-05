//! G1/G3 complete contact/radiation derivatives against perturbed primal solves.
use super::pointwise::{ambient_boundary, emissivity, enclosure, enclosure_boundary};
use super::*;
use fs_conduction::transient::enthalpy::adjoint::pointwise::workspace_entries;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::radiation::{OpenDiffuseEnclosure, ViewFactorEvidence, ViewFactorTolerance};
use fs_conduction::{
    AmbientRadiationPatch, GrayDiffuseEnclosure, LinearConfig, MaterialAmbientRadiationPatch,
    MaterialEnclosureRadiationPatch, MaterialPointwiseAmbientRadiation,
    MaterialPointwiseEnclosureRadiation, PointwiseAmbientRadiation, PointwiseEnclosureRadiation,
    PointwiseEnclosureTolerance, RadiationSurface, SurfaceEmissivity, ViewFactorMatrix,
};
use fs_conduction::radiation::pointwise::material::adjoint::MaterialEmissivityGradient;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;

const OLD: [f64; 8] = [-100., -80., -40., -20., 320., 350., 400., 450.];
const HISTORY_DIRECTION: [f64; 8] = [0.3, -0.7, 0.1, 0.4, -0.2, 0.6, 0.5, -0.4];
const SOURCE_DIRECTION: [f64; 8] = [-0.2, 0.1, 0.7, 0.3, 0.6, -0.5, 0.4, 0.8];
const GOAL: [f64; 8] = [0.7, -0.4, 0.2, 1.1, -0.3, 0.6, 0.9, -0.2];
const DT: f64 = 0.025;

#[derive(Clone, Copy)]
struct Controls {
    resistance: f64,
    epsilon: f64,
    ambient: f64,
    history: f64,
    source: f64,
}

const BASE: Controls = Controls {
    resistance: 1.,
    epsilon: 0.4,
    ambient: 320.,
    history: 0.,
    source: 0.,
};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn five_point(step: f64, f: impl Fn(f64) -> f64) -> f64 {
    (-f(2. * step) + 8. * f(step) - 8. * f(-step) + f(-2. * step)) / (12. * step)
}

fn objective(
    accepted: &fs_conduction::interface::temperature::enthalpy::ContactEnthalpyStep,
) -> f64 {
    accepted
        .conduction
        .temperature
        .iter()
        .zip(GOAL)
        .map(|(t, w)| (t - 300.) * w)
        .sum()
}

// All initial and accepted states lie inside smooth sensible chart pieces;
// the separate primal tests exercise true latent plateaus. Both nonlinear
// conduction and R(T) remain active here, with unequal reference densities.
fn solve(cx: &Cx<'_>, reflected: bool, open: bool, controls: Controls, reverse: bool) -> (f64, [f64; 5]) {
    let mesh = mesh();
    let boundary = if reflected {
        enclosure_boundary(&mesh)
    } else {
        ambient_boundary(&mesh)
    };
    let card = card_with_curve(vec![
        (280., 0.1 * controls.resistance),
        (380., 0.7 * controls.resistance),
    ]);
    let contacts =
        TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1).unwrap();
    let cold = chart(false);
    let hot = chart(true);
    let engine = HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        &mesh,
        &[
            ReferenceEnthalpyMaterial {
                curve: &cold,
                reference_density_kg_m3: 1.,
            },
            ReferenceEnthalpyMaterial {
                curve: &hot,
                reference_density_kg_m3: 2.,
            },
        ],
        &[0, 1],
        EnthalpyBudget {
            max_vertices: 8,
            max_elements: 2,
        },
    )
    .unwrap();
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280., 1.6), (380., 3.6)]).unwrap(),
    );
    let source = ScalarField::Nodal(
        SOURCE_DIRECTION
            .iter()
            .enumerate()
            .map(|(i, d)| (if i < 4 { 0. } else { 60. }) + controls.source * d)
            .collect(),
    );
    let old = std::array::from_fn::<_, 8, _>(|i| OLD[i] + controls.history * HISTORY_DIRECTION[i]);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let config = coupled_newton();
    let (objective, linearization) = if reflected {
        let model = enclosure(&mesh, [controls.epsilon, 0.65]);
        let open_model = open.then(|| OpenDiffuseEnclosure::new(
            cx, model.surfaces().to_vec(), vec![vec![0.0, 0.4], vec![0.4, 0.0]],
            ViewFactorEvidence::Analytic { geometry: "manufactured reciprocal contact/bath exchange".into() },
            ViewFactorTolerance::default(),
        ).unwrap());
        let radiation = match &open_model {
            Some(model) => PointwiseEnclosureRadiation::new_open(
                cx, &mesh, &boundary, model, controls.ambient,
                PointwiseEnclosureTolerance::default(),
            ),
            None => PointwiseEnclosureRadiation::new(
                cx, &mesh, &boundary, &model, PointwiseEnclosureTolerance::default(),
            ),
        }.unwrap();
        let accepted = contacts
            .advance_enthalpy_newton_with_pointwise_enclosure_radiation(
                cx, &engine, problem, &old, DT, config, &radiation,
            )
            .unwrap();
        let objective = objective(&accepted);
        let linearization = reverse.then(|| {
            engine
                .linearize_temperature_contact_endpoint_with_radiation(
                    cx,
                    problem,
                    &contacts,
                    &old,
                    DT,
                    config.step,
                    EnthalpyRadiation::PointwiseEnclosure(&radiation),
                    accepted.conduction,
                    config.residual_tolerance_j,
                    workspace_entries(8, mesh.boundary().len(), 2).unwrap(),
                    1,
                )
                .unwrap()
        });
        (objective, linearization)
    } else {
        let patches = [AmbientRadiationPatch::new(
            "exterior",
            emissivity("exterior", controls.epsilon),
            controls.ambient,
        )
        .unwrap()];
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let accepted = contacts
            .advance_enthalpy_newton_with_pointwise_ambient_radiation(
                cx, &engine, problem, &old, DT, config, &radiation,
            )
            .unwrap();
        let objective = objective(&accepted);
        let linearization = reverse.then(|| {
            let entries = workspace_entries(8, mesh.boundary().len(), 1).unwrap();
            for (model, cap) in [
                (EnthalpyRadiation::None, entries),
                (EnthalpyRadiation::PointwiseAmbient(&radiation), entries - 1),
            ] {
                assert!(engine.linearize_temperature_contact_endpoint_with_radiation(
                    cx, problem, &contacts, &old, DT, config.step, model,
                    accepted.conduction.clone(), config.residual_tolerance_j, cap, 1,
                ).is_err());
            }
            engine
                .linearize_temperature_contact_endpoint_with_radiation(
                    cx,
                    problem,
                    &contacts,
                    &old,
                    DT,
                    config.step,
                    EnthalpyRadiation::PointwiseAmbient(&radiation),
                    accepted.conduction,
                    config.residual_tolerance_j,
                    entries,
                    1,
                )
                .unwrap()
        });
        (objective, linearization)
    };
    let mut derivatives = [0.; 5];
    if let Some(linearization) = linearization {
        let transport = linearization.transport();
        let u = [0.2, -0.4, 0.6, 0.3, -0.1, 0.7, -0.2, 0.5];
        let v = [-0.3, 0.8, 0.1, -0.5, 0.6, 0.2, 0.4, -0.7];
        let jv = transport.apply_jacobian(cx, &v).unwrap();
        let jtu = transport.apply_jacobian_transpose(cx, &u).unwrap();
        close(dot(&u, &jv), dot(&jtu, &v), 2e-12);
        let seed = transport.temperature_pullback(cx, &GOAL).unwrap();
        let gradient = linearization
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-11,
                    max_iterations: 128,
                    restart: 8,
                },
            )
            .unwrap();
        let load = transport.nodal_load_pullback(cx, &gradient.transport.adjoint).unwrap();
        let contact = contacts.evaluate_enthalpy(
            cx, &boundary, &engine, &transport.primal().specific_enthalpy_j_kg,
        ).unwrap();
        derivatives[0] = contact.frozen_interfaces().resistance_scale_pullback(
            cx, &transport.primal().temperature, &load, 1,
        ).unwrap()[0].derivative;
        let outward: Vec<f64> = load.iter().map(|value| -value).collect();
        let radiation = transport.pointwise_parameter_pullback(cx, &outward).unwrap();
        derivatives[1] = radiation.emissivities[0];
        if reflected && !open {
            assert!(radiation.reservoir_temperatures.is_empty());
        } else {
            assert_eq!(radiation.reservoir_temperatures.len(), 1);
            derivatives[2] = radiation.reservoir_temperatures[0];
        }
        derivatives[3] = dot(
            &gradient.transport.previous_specific_enthalpy,
            &HISTORY_DIRECTION,
        );
        derivatives[4] = dot(
            &gradient.transport.source_density,
            &SOURCE_DIRECTION,
        );
    }
    (objective, derivatives)
}

#[test]
fn g1_pointwise_contact_reverse_retains_resistance_radiation_history_and_source_feedback() {
    context(&CancelGate::new_clock_free(), |cx| {
        for (reflected, open) in [(false, false), (true, false), (true, true)] {
            let (_, derivatives) = solve(cx, reflected, open, BASE, true);
            for (control, &step) in [1e-3, 1e-4, 0.05, 0.02, 0.1].iter().enumerate() {
                if reflected && !open && control == 2 {
                    continue;
                }
                let difference = five_point(step, |delta| {
                    let mut changed = BASE;
                    match control {
                        0 => changed.resistance += delta,
                        1 => changed.epsilon += delta,
                        2 => changed.ambient += delta,
                        3 => changed.history += delta,
                        4 => changed.source += delta,
                        _ => unreachable!(),
                    }
                    solve(cx, reflected, open, changed, false).0
                });
                assert!(
                    (derivatives[control] - difference).abs() < 5e-6 * (1. + difference.abs()),
                    "reflected={reflected} open={open} control={control}: adjoint={} perturbation={difference}",
                    derivatives[control]
                );
            }
        }
    });
}

fn phase_emissivity_card(name: &str, ordinates: [f64; 2]) -> fs_matdb::MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::SURFACE_EMISSIVITY_PROPERTY, fs_conduction::EMISSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "melt".into(), abscissa_dims: fs_qty::Dims::NONE,
            knots: vec![(0.0, ordinates[0]), (1.0, ordinates[1])], dims: fs_conduction::EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 280., 380.).with("melt", 0., 1.),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "synthetic phase emissivity with live contact; not measurements".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    fs_matdb::MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "chart".into(),
        process: "contact radiation derivative fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}

// Actual latent h, phase emissivity, contact exchange and reflection are solved
// afresh for each source ordinate. Temperature stays flat, so a T-only or
// apparent-capacity derivative cannot supply this enthalpy objective.
fn material_solve(cx: &Cx<'_>, mode: u8, ordinate_delta: f64, reverse: bool) -> (f64, f64) {
    let mesh = mesh();
    let boundary = if mode == 0 { ambient_boundary(&mesh) } else { enclosure_boundary(&mesh) };
    let cards = [phase_emissivity_card("cold", [0.2 + ordinate_delta, 0.8]),
        phase_emissivity_card("hot", [0.25, 0.7])];
    let curves = [0, 1].map(|body| EquilibriumEnthalpyPhaseCurve::try_new(
        cards[body].content_hash(), chart(body == 1).knots().to_vec(),
    ).unwrap());
    let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &[
        ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: 1.0 },
        ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: 2.0 },
    ], &[0, 1], EnthalpyBudget { max_vertices: 8, max_elements: 2 }).unwrap();
    let point = QueryPoint::new().with("T", 320.).unwrap().with("melt", 0.5).unwrap();
    let binding = |body: usize| ThermalMaterialPoint { material: &cards[body], context: &point,
        temperature_axis: "T", liquid_fraction_axis: Some("melt") };
    let contact_card = card();
    let contacts = TemperatureContacts::new(cx, &mesh, &boundary, vec![],
        vec![law(&mesh, &contact_card)], 1).unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    let old = [145., 150., 155., 150., 540., 545., 550., 555.];
    let config = coupled_newton();
    let value = |accepted: &fs_conduction::transient::enthalpy::EnthalpyStepSolution| {
        assert!(accepted.specific_enthalpy_j_kg[..4].iter().all(|h| *h > 100. && *h < 200.));
        assert!(accepted.specific_enthalpy_j_kg[4..].iter().all(|h| *h > 500. && *h < 600.));
        dot(&accepted.specific_enthalpy_j_kg, &GOAL)
    };
    let linear = LinearConfig { tolerance: 1e-11, max_iterations: 128, restart: 8 };
    if mode == 0 {
        let radiation = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary, &[
            MaterialAmbientRadiationPatch { region: "exterior", ambient_temperature_k: 320.,
                curve: &curves[0], material: binding(0), selection: MaterialPropertySelection::SingleClaimOnly },
        ]).unwrap();
        let accepted = contacts.advance_enthalpy_newton_with_material_ambient_radiation(
            cx, &engine, problem, &old, DT, config, &radiation,
        ).unwrap().conduction;
        let objective = value(&accepted);
        if !reverse { return (objective, 0.0); }
        let entries = radiation.adjoint_workspace_entries().unwrap();
        let mut altered = accepted.clone();
        assert!(altered.pointwise_radiation.as_mut().unwrap().patches[0].material_states.pop().is_some());
        for (candidate, cap) in [(accepted.clone(), entries - 1), (altered, entries)] {
            assert!(engine.linearize_temperature_contact_endpoint_with_radiation(
                cx, problem, &contacts, &old, DT, config.step,
                EnthalpyRadiation::MaterialAmbient(&radiation), candidate,
                config.residual_tolerance_j, cap, 1,
            ).is_err());
        }
        let fixed_patches = [
            AmbientRadiationPatch::new("exterior", SurfaceEmissivity::from_condition(
                "exterior", binding(0), MaterialPropertySelection::SingleClaimOnly,
            ).unwrap().0, 320.).unwrap(),
        ];
        let fixed = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &fixed_patches).unwrap();
        assert!(engine.linearize_temperature_contact_endpoint_with_radiation(
            cx, problem, &contacts, &old, DT, config.step,
            EnthalpyRadiation::PointwiseAmbient(&fixed), accepted.clone(),
            config.residual_tolerance_j, workspace_entries(8, mesh.boundary().len(), 1).unwrap(), 1,
        ).is_err());
        let step = engine.linearize_temperature_contact_endpoint_with_radiation(
            cx, problem, &contacts, &old, DT, config.step, EnthalpyRadiation::MaterialAmbient(&radiation),
            accepted, config.residual_tolerance_j, entries, 1,
        ).unwrap();
        assert_eq!(step.transport().temperature_pullback(cx, &GOAL).unwrap(), vec![0.0; 8]);
        let gradient = step.pullback_with_material_ambient_controls(cx, &GOAL, linear).unwrap();
        let MaterialEmissivityGradient::Curve { knot_value_derivatives, .. } =
            &gradient.radiation.patches[0].emissivity else { panic!("original source ordinates") };
        return (objective, knot_value_derivatives[0]);
    }
    let surfaces = [0_usize, 1].map(|body| {
        let name = ["left", "right"][body];
        RadiationSurface::new(&mesh, name, |face| face.vertices.iter().all(|v|
            [4 * body, 4 * body + 1, 4 * body + 3].contains(&(*v as usize))),
            SurfaceEmissivity::from_condition(name, binding(body),
                MaterialPropertySelection::SingleClaimOnly).unwrap().0,
        ).unwrap()
    }).to_vec();
    let model = GrayDiffuseEnclosure::new(surfaces.clone(),
        ViewFactorMatrix::infinite_parallel_plates(0.5).unwrap()).unwrap();
    let open = OpenDiffuseEnclosure::new(cx, surfaces, vec![vec![0.0, 0.4], vec![0.4, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "manufactured phase/contact bath exchange".into() },
        ViewFactorTolerance::default()).unwrap();
    let base = if mode == 1 {
        PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &model, PointwiseEnclosureTolerance::default())
    } else {
        PointwiseEnclosureRadiation::new_open(cx, &mesh, &boundary, &open, 320., PointwiseEnclosureTolerance::default())
    }.unwrap();
    let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &[0_usize, 1].map(|body|
        MaterialEnclosureRadiationPatch { region: ["left", "right"][body], curve: &curves[body],
            material: binding(body), selection: MaterialPropertySelection::SingleClaimOnly },
    )).unwrap();
    let accepted = contacts.advance_enthalpy_newton_with_material_enclosure_radiation(
        cx, &engine, problem, &old, DT, config, &radiation,
    ).unwrap().conduction;
    let objective = value(&accepted);
    if !reverse { return (objective, 0.0); }
    let entries = radiation.adjoint_workspace_entries().unwrap();
    let mut altered = accepted.clone();
    assert!(altered.pointwise_enclosure_radiation.as_mut().unwrap().patches[0].material_states.pop().is_some());
    for (candidate, cap) in [(accepted.clone(), entries - 1), (altered, entries)] {
        assert!(engine.linearize_temperature_contact_endpoint_with_radiation(
            cx, problem, &contacts, &old, DT, config.step,
            EnthalpyRadiation::MaterialEnclosure(&radiation), candidate,
            config.residual_tolerance_j, cap, 1,
        ).is_err());
    }
    assert!(engine.linearize_temperature_contact_endpoint_with_radiation(
        cx, problem, &contacts, &old, DT, config.step,
        EnthalpyRadiation::PointwiseEnclosure(&base), accepted.clone(),
        config.residual_tolerance_j, workspace_entries(8, mesh.boundary().len(), 2).unwrap(), 1,
    ).is_err());
    let step = engine.linearize_temperature_contact_endpoint_with_radiation(
        cx, problem, &contacts, &old, DT, config.step, EnthalpyRadiation::MaterialEnclosure(&radiation),
        accepted, config.residual_tolerance_j, entries, 1,
    ).unwrap();
    assert_eq!(step.transport().temperature_pullback(cx, &GOAL).unwrap(), vec![0.0; 8]);
    let gradient = step.pullback_with_material_enclosure_controls(cx, &GOAL, linear).unwrap();
    let MaterialEmissivityGradient::Curve { knot_value_derivatives, .. } =
        &gradient.radiation.patches[0].emissivity else { panic!("original source ordinates") };
    (objective, knot_value_derivatives[0])
}

#[test]
fn g3_live_contact_pointwise_endpoint_retains_original_material_phase_emissivity_controls() {
    context(&CancelGate::new_clock_free(), |cx| {
        for mode in 0..3 {
            let (_, actual) = material_solve(cx, mode, 0.0, true);
            let expected = five_point(1e-4, |delta| material_solve(cx, mode, delta, false).0);
            assert!(actual.abs() > 1e-5);
            close(actual, expected, 5e-6 * (1.0 + expected.abs()));
        }
    });
}

//! G1/G3 complete contact/radiation derivatives against perturbed primal solves.
use super::pointwise::{ambient_boundary, emissivity, enclosure, enclosure_boundary};
use super::*;
use fs_conduction::transient::enthalpy::adjoint::pointwise::workspace_entries;
use fs_conduction::{
    AmbientRadiationPatch, LinearConfig, PointwiseAmbientRadiation, PointwiseEnclosureRadiation,
    PointwiseEnclosureTolerance,
};

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
fn solve(cx: &Cx<'_>, reflected: bool, controls: Controls, reverse: bool) -> (f64, [f64; 5]) {
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
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &mesh,
            &boundary,
            &model,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let accepted = contacts
            .advance_enthalpy_newton_with_pointwise_enclosure_radiation(
                cx, &engine, problem, &old, DT, config, &radiation,
            )
            .unwrap();
        let objective = objective(&accepted);
        let linearization = reverse.then(|| {
            contacts
                .linearize_accepted_enthalpy_with_pointwise_enclosure_radiation(
                    cx,
                    &engine,
                    problem,
                    &old,
                    DT,
                    config.step,
                    config.residual_tolerance_j,
                    &radiation,
                    accepted,
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
            contacts
                .linearize_accepted_enthalpy_with_pointwise_ambient_radiation(
                    cx,
                    &engine,
                    problem,
                    &old,
                    DT,
                    config.step,
                    config.residual_tolerance_j,
                    &radiation,
                    accepted,
                    workspace_entries(8, mesh.boundary().len(), 1).unwrap(),
                    1,
                )
                .unwrap()
        });
        (objective, linearization)
    };
    let mut derivatives = [0.; 5];
    if let Some(linearization) = linearization {
        let transport = linearization.transport().transport();
        let u = [0.2, -0.4, 0.6, 0.3, -0.1, 0.7, -0.2, 0.5];
        let v = [-0.3, 0.8, 0.1, -0.5, 0.6, 0.2, 0.4, -0.7];
        let jv = transport.apply_jacobian(cx, &v).unwrap();
        let jtu = transport.apply_jacobian_transpose(cx, &u).unwrap();
        close(dot(&u, &jv), dot(&jtu, &v), 2e-12);
        let seed = linearization.temperature_pullback(cx, &GOAL).unwrap();
        let gradient = linearization
            .pullback_with_pointwise_radiation_controls(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-11,
                    max_iterations: 128,
                    restart: 8,
                },
                1,
            )
            .unwrap();
        derivatives[0] = gradient.thermal.resistance_multipliers[0].derivative;
        derivatives[1] = gradient.radiation.emissivities[0];
        if reflected {
            assert!(gradient.radiation.reservoir_temperatures.is_empty());
        } else {
            derivatives[2] = gradient.radiation.reservoir_temperatures[0];
        }
        derivatives[3] = dot(
            &gradient
                .thermal
                .storage
                .transport
                .previous_specific_enthalpy,
            &HISTORY_DIRECTION,
        );
        derivatives[4] = dot(
            &gradient.thermal.storage.transport.source_density,
            &SOURCE_DIRECTION,
        );
    }
    (objective, derivatives)
}

#[test]
fn g1_pointwise_contact_reverse_retains_resistance_radiation_history_and_source_feedback() {
    context(&CancelGate::new_clock_free(), |cx| {
        for reflected in [false, true] {
            let (_, derivatives) = solve(cx, reflected, BASE, true);
            for (control, &step) in [1e-3, 1e-4, 0.05, 0.02, 0.1].iter().enumerate() {
                if reflected && control == 2 {
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
                    solve(cx, reflected, changed, false).0
                });
                assert!(
                    (derivatives[control] - difference).abs() < 5e-6 * (1. + difference.abs()),
                    "reflected={reflected} control={control}: adjoint={} perturbation={difference}",
                    derivatives[control]
                );
            }
        }
    });
}

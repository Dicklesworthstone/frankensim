//! Shared circuit storage, nonlinear field, electrical/thermal contacts and radiation.
use super::*;
use crate::contacts::{
    binding as electrical_binding, contact_mesh, interface as electrical_interface,
    interface_context,
};
use fs_couple::circuit_enthalpy::conduction::circuit::{
    ConductionCircuitTerminal, response::ConductionResponseCircuit,
};
use fs_em::conduction::ConductionContact;
use fs_material::electrical_interface::resolve_electrical_interface_state;
use fs_phs::circuit::{Branch, CircuitGraph};

const C: f64 = 0.5;
const L: f64 = 0.05;
const CONTACT_R: f64 = 0.5;
const FRACTION: f64 = 0.25;

fn circuit(field: PlanarConductionResponse) -> ConductionResponseCircuit {
    let contacted = !field.contacts().is_empty();
    ConductionResponseCircuit::new(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 2, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::Inductor { henries: L }),
                (1, 0, Branch::Capacitor { farads: C }),
            ],
            transformers: vec![],
        },
        ConductionCircuitTerminal {
            branch: 0,
            positive: if contacted { vec![5, 6] } else { vec![2, 5] },
            negative: vec![0, 3],
        },
        0.0,
        vec![0.0],
        vec![3.0],
        &[],
        limits(),
        field_controls(),
        controls(),
        || false,
    )
    .unwrap()
}

// Independent two-body energy elimination and midpoint RLC current equation.
// Radiation is the closed two-surface gray exchange, and the thermal contact
// retains the fixture's original first-surface temperature dependence.
fn reference(
    dt: f64,
    q: f64,
    flux: f64,
    old: [f64; 2],
    emissivity: f64,
    resistance: f64,
) -> (f64, [f64; 2], [f64; 2], f64) {
    moving_reference(dt, q, flux, old, emissivity, resistance, CONTACT_R, 0.0)
}

fn moving_reference(
    dt: f64,
    q: f64,
    flux: f64,
    old: [f64; 2],
    emissivity: f64,
    resistance: f64,
    contact_r: f64,
    velocity_y: f64,
) -> (f64, [f64; 2], [f64; 2], f64) {
    let endpoint = |j: f64| {
        let d = BETA.map(|b| 1.0 - b * j.abs());
        let a = [0, 1].map(|i| CP[i] - dt * j * j * ALPHA[i] / d[i]);
        let total = CP[0] * old[0]
            + CP[1] * old[1]
            + dt * j * j * (RHO[0] / d[0] + RHO[1] / d[1] + contact_r);
        let balance = |theta1: f64| {
            let theta2 = (total - a[0] * theta1) / a[1];
            let e1 = emissivity * (1.0 + 0.0001 * theta1);
            let e2 = 0.6 * (1.0 + 0.0001 * theta2);
            let radiated = 5.670_374_419e-8 * ((300.0 + theta1).powi(4) - (300.0 + theta2).powi(4))
                / (1.0 / e1 + 1.0 / e2 - 1.0);
            let transfer =
                (theta1 - theta2) / (1.0 + resistance * (1.0 + 0.001 * theta1)) + radiated;
            let power1 = j * j * (RHO[0] + ALPHA[0] * theta1) / d[0] + FRACTION * j * j * contact_r;
            (
                theta2,
                CP[0] * (theta1 - old[0]) - dt * (power1 - transfer),
                radiated,
            )
        };
        let (mut lo, mut hi) = (0.0, total / a[0]);
        assert!(balance(lo).1 <= 1e-12 && balance(hi).1 >= -1e-12);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if balance(mid).1 < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let theta1 = 0.5 * (lo + hi);
        let (theta2, _, radiated) = balance(theta1);
        let theta = [theta1, theta2];
        let fields = [0, 1].map(|i| j * (RHO[i] + ALPHA[i] * theta[i]) / d[i]);
        (theta, fields, radiated)
    };
    let rhs = q / C + 2.0 * flux / dt - 2.0 * velocity_y;
    let (mut lo, mut hi) = (0.0, rhs / (RHO[0] + RHO[1] + contact_r));
    for _ in 0..100 {
        let j = 0.5 * (lo + hi);
        let (_, fields, _) = endpoint(j);
        if fields.iter().sum::<f64>() + (contact_r + 2.0 * L / dt + dt / (2.0 * C)) * j < rhs {
            lo = j;
        } else {
            hi = j;
        }
    }
    let j = 0.5 * (lo + hi);
    let (theta, fields, radiated) = endpoint(j);
    (j, theta, fields, radiated)
}

#[test]
fn nonlinear_circuit_contact_and_radiating_materials_share_one_accepted_heat_history() {
    let mut temperatures = Vec::new();
    for (emissivity, resistance) in [(0.2, 0.2), (0.8, 0.2), (0.2, 1.5)] {
        with_heat(
            emissivity,
            resistance,
            BETA,
            Support::Ordinary,
            |cx, _, assembly, mut thermal, field| {
                let cards = [
                    field.sources()[0].card().clone(),
                    field.sources()[2].card().clone(),
                ];
                let interface = electrical_interface(&cards, CONTACT_R, 0.0, 1000.0);
                let context = interface_context(300.0);
                let contact = ConductionContact {
                    first_vertices: [1, 2],
                    second_vertices: [4, 7],
                    source: resolve_electrical_interface_state(
                        &interface,
                        &context,
                        MaterialPropertySelection::SingleClaimOnly,
                    )
                    .unwrap(),
                };
                let field = PlanarConductionResponse::new_with_phase_contacts(
                    contact_mesh(),
                    field.sources().to_vec(),
                    vec![contact],
                    vec![[None, None]],
                    limits(),
                    &mut || false,
                )
                .unwrap();
                let contacts = [electrical_binding(&cards, &interface, &context, FRACTION)];
                let mut actual = circuit(field);
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let (mut q, mut flux, mut old) = (3.0, 0.0, [0.0; 2]);
                for interval in 1..=3 {
                    let (j, theta, fields, radiation) =
                        reference(0.1, q, flux, old, emissivity, resistance);
                    let next = selected
                        .step_conduction_response_contact_circuit_candidate(
                            cx,
                            &thermal,
                            interval as f64 * 0.1,
                            &actual,
                            &[],
                            &[],
                            &contacts,
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 2],
                        )
                        .unwrap();
                    let report = next.report();
                    let electric = &report.electrical.circuit;
                    near(electric.spatial_current_a, j);
                    near(
                        electric.spatial_voltage_v,
                        fields[0] + fields[1] + CONTACT_R * j,
                    );
                    near(
                        electric.spatial.contacts[0].response.power_w,
                        j * j * CONTACT_R,
                    );
                    near(
                        report.deposited_joule_j,
                        0.1 * j * electric.spatial_voltage_v,
                    );
                    near(
                        report.source.thermal_contacts[0].conductance_w_k,
                        1.0 / (1.0 + resistance * (1.0 + 0.001 * theta[0])),
                    );
                    let emitted = report.source.radiation.as_ref().unwrap();
                    near(emitted.sources_w()[&CellId(30)], -radiation);
                    near(emitted.sources_w()[&CellId(10)], radiation);
                    assert!(radiation.abs() > 0.0);
                    near(
                        report.non_electrical_heat_j - 0.1 * report.source.radiation_net_power_w(),
                        0.0,
                    );
                    for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                        near(
                            next.thermal().temperature_k(cell).unwrap(),
                            300.0 + theta[i],
                        );
                        for triangle in 2 * i..2 * i + 2 {
                            near(electric.spatial.cell_power_w[triangle], 0.5 * j * fields[i]);
                        }
                    }
                    q -= 0.1 * j;
                    flux = 2.0 * L * j - flux;
                    near(next.circuit().charges_c()[0], q);
                    near(next.circuit().fluxes_wb()[0], flux);
                    near(
                        q * q / (2.0 * C)
                            + flux * flux / (2.0 * L)
                            + CP[0] * theta[0]
                            + CP[1] * theta[1],
                        9.0,
                    );
                    if interval == 3 {
                        temperatures.push(theta);
                    }
                    old = theta;
                    (thermal, actual, _) = next.into_parts();
                }
            },
        );
    }
    assert!((temperatures[0][0] - temperatures[1][0]).abs() > 1e-7);
    assert!((temperatures[0][0] - temperatures[2][0]).abs() > 1e-7);
}

#[test]
fn circuit_heat_rejects_unowned_field_boundaries_and_preserves_original_histories() {
    for support in [
        Support::EmissivityField,
        Support::ContactField,
        Support::ShortEmissivity,
    ] {
        with_heat(
            0.2,
            0.2,
            BETA,
            support,
            |cx, _, assembly, thermal, field| {
                let actual = circuit(field);
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let before = thermal.energy_j();
                assert!(
                    selected
                        .step_conduction_response_contact_circuit_candidate(
                            cx,
                            &thermal,
                            0.5,
                            &actual,
                            &[],
                            &[],
                            &[],
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 2],
                        )
                        .is_err()
                );
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(thermal.energy_j(), before);
                assert_eq!(actual.time_s(), 0.0);
                assert_eq!(actual.charges_c(), &[3.0]);
                assert_eq!(actual.fluxes_wb(), &[0.0]);
            },
        );
    }
}

#[test]
fn absent_surface_heat_paths_reuse_the_original_circuit_solution_exactly() {
    with_heat(
        0.2,
        0.2,
        BETA,
        Support::Ordinary,
        |cx, _, assembly, thermal, field| {
            let actual = circuit(field);
            let bare = MaterialHeatNetwork {
                electrical: assembly.electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let selected = bare.prepare(cx, &thermal).unwrap();
            let complete = selected
                .step_conduction_response_contact_circuit_candidate(
                    cx,
                    &thermal,
                    0.1,
                    &actual,
                    &[],
                    &[],
                    &[],
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            let original = network(assembly.electrical.bulk.cells)
                .step_response_circuit_candidate(
                    &thermal,
                    0.1,
                    &actual,
                    &[],
                    &[],
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            assert_eq!(complete.thermal().cells(), original.thermal().cells());
            assert_eq!(
                complete.circuit().charges_c(),
                original.circuit().charges_c()
            );
            assert_eq!(
                complete.circuit().fluxes_wb(),
                original.circuit().fluxes_wb()
            );
            assert_eq!(
                complete.report().electrical.circuit.spatial.cell_power_w,
                original.report().electrical.circuit.spatial.cell_power_w
            );
        },
    );
}

fn moving_circuit(
    field: PlanarConductionResponse,
    charge_c: f64,
    motion: fs_em::conduction::motion::ConductionMotion<'_>,
) -> ConductionResponseCircuit {
    ConductionResponseCircuit::new_with_motion(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 2, Branch::Resistor { ohms: 999.0 }),
                (2, 0, Branch::Inductor { henries: L }),
                (1, 0, Branch::Capacitor { farads: C }),
            ],
            transformers: vec![],
        },
        ConductionCircuitTerminal {
            branch: 0,
            positive: vec![2, 5],
            negative: vec![0, 3],
        },
        0.0,
        vec![0.0],
        vec![charge_c],
        &[],
        motion,
        limits(),
        field_controls(),
        controls(),
        1000,
        || false,
    )
    .unwrap()
}

#[test]
fn moving_nonlinear_rlc_and_complete_heat_share_motor_generator_storage() {
    use fs_couple::circuit_enthalpy::conduction::circuit::response::thermal::step_moving_response_circuit_candidate;
    use fs_em::conduction::motion::ConductionMotion;

    for (charge_c, velocity_y) in [(3.0, 0.5), (0.0, -2.0), (3.0, 0.0)] {
        with_heat(
            0.2,
            0.2,
            BETA,
            Support::Ordinary,
            |cx, _, assembly, mut thermal, field| {
                let velocity = [[0.0, velocity_y]; 4];
                let motion = ConductionMotion {
                    velocity_m_s: &velocity,
                    normal_flux_density_t: &[1.0; 4],
                };
                let mut actual = moving_circuit(field, charge_c, motion);
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let (mut q, mut flux, mut old) = (charge_c, 0.0, [0.0; 2]);
                let mut mechanical_total = 0.0;
                for interval in 1..=2 {
                    let (j, theta, fields, radiation) =
                        moving_reference(0.1, q, flux, old, 0.2, 0.2, 0.0, velocity_y);
                    let next = step_moving_response_circuit_candidate(
                        &selected,
                        cx,
                        &thermal,
                        interval as f64 * 0.1,
                        &actual,
                        &[],
                        motion,
                        &[],
                        &[],
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                    )
                    .unwrap();
                    let report = next.report();
                    let step = &report.electrical;
                    let electric = &step.circuit;
                    let moving = step.motion.as_ref().unwrap();
                    let terminal_v = fields[0] + fields[1] + 2.0 * velocity_y;
                    let heat_j = 0.1 * j * (fields[0] + fields[1]);
                    let mechanical_j = 0.1 * 2.0 * velocity_y * j;
                    near(electric.spatial_current_a, j);
                    near(electric.spatial_voltage_v, terminal_v);
                    near(moving.terminal_power_w, j * terminal_v);
                    near(step.mechanical_work_j, mechanical_j);
                    near(report.deposited_joule_j, heat_j);
                    near(electric.heat_j, heat_j);
                    near(electric.supplied_work_j, 0.0);
                    near(electric.stored_change_j + heat_j + mechanical_j, 0.0);
                    near(report.total_energy_residual_j, 0.0);
                    assert!(report.circuit_trials >= step.field_solves);
                    near(
                        report.source.thermal_contacts[0].conductance_w_k,
                        1.0 / (1.0 + 0.2 * (1.0 + 0.001 * theta[0])),
                    );
                    let emitted = report.source.radiation.as_ref().unwrap();
                    near(emitted.sources_w()[&CellId(30)], -radiation);
                    near(emitted.sources_w()[&CellId(10)], radiation);
                    for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                        near(
                            next.thermal().temperature_k(cell).unwrap(),
                            300.0 + theta[i],
                        );
                        for triangle in 2 * i..2 * i + 2 {
                            near(moving.conductor_field_v_m[triangle][0], -fields[i]);
                            near(electric.spatial.current_density_a_m2[triangle][0], -j);
                            near(electric.spatial.cell_power_w[triangle], 0.5 * j * fields[i]);
                            near(
                                step.evaluated[triangle].conductor.material().query().axes()["operating"],
                                fields[i],
                            );
                        }
                    }
                    if charge_c == 0.0 && interval == 1 {
                        assert!(moving.terminal_power_w < 0.0);
                        assert!(electric.stored_change_j > 0.0);
                        assert!(step.mechanical_work_j < 0.0);
                    }
                    q -= 0.1 * j;
                    flux = 2.0 * L * j - flux;
                    mechanical_total += mechanical_j;
                    near(next.circuit().charges_c()[0], q);
                    near(next.circuit().fluxes_wb()[0], flux);
                    near(
                        q * q / (2.0 * C)
                            + flux * flux / (2.0 * L)
                            + CP[0] * theta[0]
                            + CP[1] * theta[1]
                            + mechanical_total,
                        charge_c * charge_c / (2.0 * C),
                    );
                    old = theta;
                    (thermal, actual, _) = next.into_parts();
                }
            },
        );
    }
}

#[test]
fn moving_circuit_complete_heat_budget_cancel_and_source_refusal_are_atomic() {
    use fs_couple::circuit_enthalpy::conduction::{
        circuit::response::thermal::step_moving_response_circuit_candidate,
        heat::ConductionHeatError,
    };
    use fs_em::conduction::motion::ConductionMotion;

    for support in [
        Support::Ordinary,
        Support::EmissivityField,
        Support::ContactField,
    ] {
        with_heat(
            0.2,
            0.2,
            BETA,
            support,
            |cx, gate, assembly, thermal, field| {
                let motion = ConductionMotion {
                    velocity_m_s: &[[0.0, -2.0]; 4],
                    normal_flux_density_t: &[1.0; 4],
                };
                let actual = moving_circuit(field, 0.0, motion);
                let selected = assembly.prepare(cx, &thermal).unwrap();
                let before = thermal.energy_j();
                let call = |solve| {
                    step_moving_response_circuit_candidate(
                        &selected,
                        cx,
                        &thermal,
                        0.1,
                        &actual,
                        &[],
                        motion,
                        &[],
                        &[],
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        solve,
                        &[0.0; 2],
                    )
                };
                if support == Support::Ordinary {
                    let mut exhausted = config();
                    exhausted.maximum_circuit_trials = 2;
                    assert!(matches!(
                        call(exhausted),
                        Err(CircuitEnthalpyError::CircuitBudget)
                    ));
                    assert!(call(config()).is_ok());
                    gate.request();
                    assert!(matches!(
                        call(config()),
                        Err(CircuitEnthalpyError::Cancelled)
                    ));
                } else {
                    assert!(matches!(
                        call(config()),
                        Err(CircuitEnthalpyError::Source(ConductionHeatError::Heat(_)))
                    ));
                }
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(thermal.energy_j(), before);
                assert_eq!(actual.time_s(), 0.0);
                assert_eq!(actual.charges_c(), &[0.0]);
                assert_eq!(actual.fluxes_wb(), &[0.0]);
            },
        );
    }
}

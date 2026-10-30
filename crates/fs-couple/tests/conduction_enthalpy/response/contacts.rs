//! G1/G3/G4: nonlinear spatial material, contact sources and finite heat.
#[path = "contacts/field_storage.rs"]
mod field_storage;
use super::super::contacts::{binding, contact_mesh, interface, interface_context};
use super::*;
use fs_couple::circuit_enthalpy::conduction::circuit::{
    ConductionCircuitTerminal, response::ConductionResponseCircuit,
};
use fs_couple::circuit_enthalpy::{
    contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    material::MaterialNetwork,
};
use fs_em::conduction::ConductionContact;
use fs_exec::{CancelGate, Cx};
use fs_matdb::InterfaceSystemCard;
use fs_material::electrical_interface::resolve_electrical_interface_state;
use fs_phs::circuit::{Branch, CircuitGraph};

const FRACTION: f64 = 0.25;
const CONTACT_SLOPE: f64 = 0.002;
const BETA: [f64; 2] = [0.03, 0.05];

fn specimen(
    cards: &[MaterialCard; 2],
    card: &InterfaceSystemCard,
    coordinate: ConductorMagnitude,
) -> PlanarConductionResponse {
    let (_, sources) = problem(cards, coordinate).into_frozen_parts();
    let source = resolve_electrical_interface_state(
        card,
        &interface_context(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PlanarConductionResponse::new_with_phase_contacts(
        contact_mesh(),
        sources.to_vec(),
        vec![ConductionContact {
            first_vertices: [1, 2],
            second_vertices: [4, 7],
            source,
        }],
        vec![[None, None]],
        limits(),
        &mut || false,
    )
    .unwrap()
}

// Independent elimination of the two BE thermal equations at a specified J.
// Bulk rho_i=r_i+alpha_i*theta_i+beta_i*|E_i|, hence
// E_i=J*(r_i+alpha_i*theta_i)/(1-beta_i*J). Surface R=rc0*(1+k*(theta_a+theta_b)/2).
fn endpoint(j: f64, dt: f64, old: [f64; 2], r: [f64; 2], rc0: f64) -> ([f64; 2], [f64; 2], f64) {
    let inv = BETA.map(|b| 1.0 / (1.0 - b * j));
    let g = rc0 * CONTACT_SLOPE / 2.0;
    let p = dt * j * j;
    let a = CP[0] + dt - p * (ALPHA[0] * inv[0] + FRACTION * g);
    let b = -dt - p * FRACTION * g;
    let c = -dt - p * (1.0 - FRACTION) * g;
    let d = CP[1] + dt - p * (ALPHA[1] * inv[1] + (1.0 - FRACTION) * g);
    let x = CP[0] * old[0] + p * (r[0] * inv[0] + FRACTION * rc0);
    let y = CP[1] * old[1] + p * (r[1] * inv[1] + (1.0 - FRACTION) * rc0);
    let det = a * d - b * c;
    let t = [(x * d - b * y) / det, (a * y - c * x) / det];
    let e = std::array::from_fn(|i| j * (r[i] + ALPHA[i] * t[i]) * inv[i]);
    (t, e, rc0 * (1.0 + CONTACT_SLOPE * (t[0] + t[1]) / 2.0))
}

#[test]
fn nonlinear_bulk_and_contact_temperature_response_share_current_and_voltage_heat() {
    let coordinate = ConductorMagnitude::ElectricField;
    let r = [1.0, 2.0];
    for rc0 in [0.5, 1.5] {
        let (a, ca) = material(r[0], CP[0], ALPHA[0], BETA[0], coordinate, false);
        let (b, cb) = material(r[1], CP[1], ALPHA[1], BETA[1], coordinate, false);
        let (cards, curves) = ([a, b], [ca, cb]);
        let thermal = initial(&curves, 999.0, None);
        let card = interface(&cards, rc0, CONTACT_SLOPE, 1000.0);
        let original = specimen(&cards, &card, coordinate);
        let context = point(300.0, coordinate);
        let contact_context = interface_context(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let layout = network(&cells);
        let bindings = [binding(&cards, &card, &contact_context, FRACTION)];
        let (j, dt) = (2.0, 0.1);
        let (theta, e, rc) = endpoint(j, dt, [0.0; 2], r, rc0);
        let voltage = e[0] + e[1] + j * rc;
        let drive = ConductionDrive::Current {
            positive: &[5, 6],
            negative: &[0, 3],
            current_a: j,
        };
        let next = layout
            .step_response_contact_candidate(
                &thermal,
                dt,
                &original,
                drive,
                &bindings,
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let step = &next.report().electrical;
        for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
            near(
                next.thermal().temperature_k(cell).unwrap(),
                300.0 + theta[i],
            );
            for triangle in 2 * i..2 * i + 2 {
                near(step.field.electric_field_v_m[triangle][0], -e[i]);
                near(step.field.current_density_a_m2[triangle][0], -j);
                near(
                    step.field.evaluated[triangle]
                        .conductor
                        .material()
                        .query()
                        .axes()["operating"],
                    e[i],
                );
            }
        }
        near(step.field.voltage_v[5], voltage);
        near(step.contacts[0].response.power_w, j * j * rc);
        near(step.field.contact_response[0].power_w, j * j * rc);
        near(
            step.field.cell_power_w.iter().sum::<f64>(),
            j * (e[0] + e[1]),
        );
        near(step.heat_j, dt * j * voltage);
        near(next.report().deposited_joule_j, step.heat_j);
        near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
        near(next.problem().contacts()[0].source.resistance_ohm_m2(), rc);
        assert!(step.field.field_solves > 1);
        let boundary = [(0, 0.0), (3, 0.0), (5, voltage), (6, voltage)];
        let by_voltage = layout
            .step_response_contact_candidate(
                &thermal,
                dt,
                &original,
                ConductionDrive::Voltage(&boundary),
                &bindings,
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        near(
            by_voltage.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + theta[0],
        );
        near(
            by_voltage.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + theta[1],
        );
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &crate::LINKS,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 73,
                    kernel_id: 51,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let selected = heat.prepare(&cx, &thermal).unwrap();
            let complete = selected
                .step_conduction_response_contacts_candidate(
                    &cx,
                    &thermal,
                    dt,
                    &original,
                    drive,
                    &bindings,
                    &crate::REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            near(
                complete.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + theta[0],
            );
            near(
                complete.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + theta[1],
            );
            near(
                complete.report().electrical.contacts[0].response.power_w,
                j * j * rc,
            );
        });
        assert!(
            layout
                .step_response_candidate(
                    &thermal,
                    dt,
                    &original,
                    drive,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
    }
}

#[test]
fn continuing_rlc_current_and_contact_heat_follow_both_nonlinear_material_sources() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (r, dt, cap, inductance) = ([1.0, 2.0], 0.1, 0.5, 0.05);
    let mut first_current = Vec::new();
    for rc0 in [0.5, 1.5] {
        let (a, ca) = material(r[0], CP[0], ALPHA[0], BETA[0], coordinate, false);
        let (b, cb) = material(r[1], CP[1], ALPHA[1], BETA[1], coordinate, false);
        let (cards, curves) = ([a, b], [ca, cb]);
        let card = interface(&cards, rc0, CONTACT_SLOPE, 1000.0);
        let context = point(300.0, coordinate);
        let contact_context = interface_context(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let layout = network(&cells);
        let bindings = [binding(&cards, &card, &contact_context, FRACTION)];
        let mut thermal = initial(&curves, 999.0, None);
        let mut actual = ConductionResponseCircuit::new(
            specimen(&cards, &card, coordinate),
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (1, 2, Branch::Resistor { ohms: 999.0 }),
                    (
                        2,
                        0,
                        Branch::Inductor {
                            henries: inductance,
                        },
                    ),
                    (1, 0, Branch::Capacitor { farads: cap }),
                ],
                transformers: vec![],
            },
            ConductionCircuitTerminal {
                branch: 0,
                positive: vec![5, 6],
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
        .unwrap();
        let (mut q, mut flux, mut old) = (3.0, 0.0, [0.0; 2]);
        for step in 1..=3 {
            // Eliminate native midpoint q/flux to a scalar current balance,
            // independently of the spatial and constitutive solvers.
            let rhs = q / cap + 2.0 * flux / dt;
            let (mut lo, mut hi) = (0.0, rhs / (r[0] + r[1] + rc0));
            for _ in 0..100 {
                let j = 0.5 * (lo + hi);
                let (_, e, rc) = endpoint(j, dt, old, r, rc0);
                if e[0] + e[1] + j * rc + (2.0 * inductance / dt + dt / (2.0 * cap)) * j < rhs {
                    lo = j;
                } else {
                    hi = j;
                }
            }
            let j = 0.5 * (lo + hi);
            let (theta, e, rc) = endpoint(j, dt, old, r, rc0);
            let next = layout
                .step_response_contact_circuit_candidate(
                    &thermal,
                    step as f64 * dt,
                    &actual,
                    &[],
                    &[],
                    &bindings,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            q -= dt * j;
            flux = 2.0 * inductance * j - flux;
            near(next.circuit().charges_c()[0], q);
            near(next.circuit().fluxes_wb()[0], flux);
            let record = &next.report().electrical;
            near(record.circuit.spatial_current_a, j);
            near(record.circuit.spatial_voltage_v, e[0] + e[1] + j * rc);
            near(
                record.circuit.spatial.contacts[0].response.power_w,
                j * j * rc,
            );
            near(record.circuit.heat_j, dt * j * (e[0] + e[1] + j * rc));
            near(next.report().deposited_joule_j, record.circuit.heat_j);
            near(next.circuit().contacts()[0].source.resistance_ohm_m2(), rc);
            for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                near(
                    next.thermal().temperature_k(cell).unwrap(),
                    300.0 + theta[i],
                );
            }
            near(
                q * q / (2.0 * cap)
                    + flux * flux / (2.0 * inductance)
                    + CP[0] * theta[0]
                    + CP[1] * theta[1],
                9.0,
            );
            assert!(record.field_solves > 1);
            if step == 1 {
                first_current.push(j);
            }
            old = theta;
            (thermal, actual, _) = next.into_parts();
        }
    }
    assert!(first_current[0] > first_current[1]);
}

#[test]
fn nonlinear_contact_trials_keep_original_histories_on_budget_and_late_cancellation() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], BETA[0], coordinate, false);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], BETA[1], coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let thermal = initial(&curves, 999.0, None);
    let card = interface(&cards, 0.5, CONTACT_SLOPE, 1000.0);
    let original = specimen(&cards, &card, coordinate);
    let context = point(300.0, coordinate);
    let contact_context = interface_context(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let layout = network(&cells);
    let bindings = [binding(&cards, &card, &contact_context, FRACTION)];
    let drive = ConductionDrive::Current {
        positive: &[5, 6],
        negative: &[0, 3],
        current_a: 2.0,
    };
    let old = original.conductors().to_vec();
    let surface = original.contacts().to_vec();
    let energy = thermal.energy_j();
    let mut budget = config();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_response_contact_candidate(
            &thermal,
            0.1,
            &original,
            drive,
            &bindings,
            limits(),
            field_controls(),
            controls(),
            budget,
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut polls = 0;
    assert!(matches!(
        layout.step_response_contact_candidate(
            &thermal,
            0.1,
            &original,
            drive,
            &bindings,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                polls > 100
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(original.conductors(), old);
    assert_eq!(original.contacts(), surface);
    assert_eq!(thermal.energy_j().to_bits(), energy.to_bits());
    assert_eq!(thermal.time_s(), 0.0);
    layout
        .step_response_contact_candidate(
            &thermal,
            0.1,
            &original,
            drive,
            &bindings,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
}

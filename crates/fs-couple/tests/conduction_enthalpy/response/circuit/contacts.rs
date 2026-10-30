//! Independent RLC, original contact and finite-capacity heat balances.
use super::*;
use crate::contacts::{binding, contact_mesh, interface, interface_context};
use fs_em::conduction::ConductionContact;
use fs_matdb::InterfaceSystemCard;
use fs_material::electrical_interface::resolve_electrical_interface_state;

fn contacted_problem(
    cards: &[MaterialCard; 2],
    interface: &InterfaceSystemCard,
) -> PlanarConductionResponse {
    let bulk = problem(cards, ConductorMagnitude::CurrentDensity);
    let surface = resolve_electrical_interface_state(
        interface,
        &interface_context(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PlanarConductionResponse::new_with_phase_contacts(
        contact_mesh(),
        bulk.sources().to_vec(),
        vec![ConductionContact {
            first_vertices: [1, 2],
            second_vertices: [4, 7],
            source: surface,
        }],
        vec![[None, None]],
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn contacted_circuit(
    cards: &[MaterialCard; 2],
    interface: &InterfaceSystemCard,
    sign: f64,
) -> ConductionResponseCircuit {
    ConductionResponseCircuit::new(
        contacted_problem(cards, interface),
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
            positive: vec![5, 6],
            negative: vec![0, 3],
        },
        0.0,
        vec![0.0],
        vec![sign * 3.0],
        &[],
        limits(),
        field_controls(),
        controls(),
        || false,
    )
    .unwrap()
}

const BETA: [f64; 2] = [0.05, 0.08];
const CONTACT_SLOPE: f64 = 0.003;
const FRACTION: f64 = 0.25;

// For a given |I|, rho_i = r_i + alpha_i theta_i + beta_i |I| and
// R_contact = r0(1+gamma(theta1+theta2)/2). The two BE heat balances are
// linear. Solve them independently, then root the native midpoint RLC law.
fn reference(
    dt: f64,
    q: f64,
    flux: f64,
    old: [f64; 2],
    r: [f64; 2],
    rc: f64,
) -> (f64, [f64; 2], [f64; 3]) {
    let endpoint = |j: f64| {
        let k = dt * j * j;
        let slope = rc * CONTACT_SLOPE / 2.0;
        let a11 = CP[0] + dt - k * (ALPHA[0] + FRACTION * slope);
        let a12 = -dt - k * FRACTION * slope;
        let a21 = -dt - k * (1.0 - FRACTION) * slope;
        let a22 = CP[1] + dt - k * (ALPHA[1] + (1.0 - FRACTION) * slope);
        let rhs1 = CP[0] * old[0] + k * (r[0] + BETA[0] * j + FRACTION * rc);
        let rhs2 = CP[1] * old[1] + k * (r[1] + BETA[1] * j + (1.0 - FRACTION) * rc);
        let det = a11 * a22 - a12 * a21;
        let theta = [
            (rhs1 * a22 - a12 * rhs2) / det,
            (a11 * rhs2 - a21 * rhs1) / det,
        ];
        let drops = [
            j * (r[0] + ALPHA[0] * theta[0] + BETA[0] * j),
            j * (r[1] + ALPHA[1] * theta[1] + BETA[1] * j),
            j * rc * (1.0 + CONTACT_SLOPE * 0.5 * (theta[0] + theta[1])),
        ];
        (theta, drops)
    };
    let rhs = q / C + 2.0 * flux / dt;
    let (mut lo, mut hi) = (0.0, rhs / (r[0] + r[1] + rc));
    for _ in 0..100 {
        let j = 0.5 * (lo + hi);
        let (_, drops) = endpoint(j);
        if drops.iter().sum::<f64>() + (2.0 * L / dt + dt / (2.0 * C)) * j < rhs {
            lo = j;
        } else {
            hi = j;
        }
    }
    let j = 0.5 * (lo + hi);
    let (theta, drops) = endpoint(j);
    (j, theta, drops)
}

#[test]
fn original_bulk_and_contact_materials_drive_continuing_rlc_and_two_sided_finite_heat() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let mut currents = Vec::new();
    for (r, rc) in [([1.0, 2.0], 0.5), ([2.0, 2.0], 0.5), ([1.0, 2.0], 1.5)] {
        let (a, ca) = material(r[0], CP[0], ALPHA[0], BETA[0], coordinate, false);
        let (b, cb) = material(r[1], CP[1], ALPHA[1], BETA[1], coordinate, false);
        let (cards, curves) = ([a, b], [ca, cb]);
        let context = point(300.0, coordinate);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let layout = network(&cells);
        let surface = interface(&cards, rc, CONTACT_SLOPE, 1000.0);
        let surface_context = interface_context(300.0);
        let contact_bindings = [binding(&cards, &surface, &surface_context, FRACTION)];
        for sign in [1.0, -1.0] {
            let mut thermal = initial(&curves, 999.0, None);
            let mut actual = contacted_circuit(&cards, &surface, sign);
            let (mut q, mut flux, mut old) = (3.0, 0.0, [0.0; 2]);
            for interval in 1..=3 {
                let (j, theta, drops) = reference(0.1, q, flux, old, r, rc);
                let next = layout
                    .step_response_contact_circuit_candidate(
                        &thermal,
                        interval as f64 * 0.1,
                        &actual,
                        &[],
                        &[],
                        &contact_bindings,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                let electrical = &next.report().electrical;
                let native = &electrical.circuit;
                near(native.spatial_current_a, sign * j);
                near(native.spatial_voltage_v, sign * drops.iter().sum::<f64>());
                assert!(electrical.field_solves > 1);
                assert!(next.report().circuit_trials >= electrical.field_solves);
                let contact = &native.spatial.contacts[0];
                assert_eq!(contact.cells, [CellId(30), CellId(10)]);
                assert_eq!(contact.heat_fraction_to_first, FRACTION);
                near(contact.response.current_a, -sign * j);
                near(contact.response.voltage_jump_v[0], -sign * drops[2]);
                near(contact.response.power_w, j * drops[2]);
                near(
                    native.spatial.cell_power_w.iter().sum::<f64>() + contact.response.power_w,
                    native.heat_j / 0.1,
                );
                near(next.report().deposited_joule_j, native.heat_j);
                for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
                    near(
                        next.thermal().temperature_k(cell).unwrap(),
                        300.0 + theta[i],
                    );
                    for triangle in 2 * i..2 * i + 2 {
                        near(native.spatial.cell_power_w[triangle], 0.5 * j * drops[i]);
                        near(
                            electrical.evaluated[triangle]
                                .conductor
                                .material()
                                .query()
                                .axes()["operating"],
                            j,
                        );
                    }
                }
                let source = next.circuit().contacts()[0].source.resolved();
                near(source.state_coordinate("Ta").unwrap(), 300.0 + theta[0]);
                near(source.state_coordinate("Tb").unwrap(), 300.0 + theta[1]);
                near(
                    source.state_coordinate("Tm").unwrap(),
                    300.0 + 0.5 * (theta[0] + theta[1]),
                );
                near(source.state_coordinate("pressure").unwrap(), 1e5);
                q -= 0.1 * j;
                flux = 2.0 * L * j - flux;
                near(next.circuit().charges_c()[0], sign * q);
                near(next.circuit().fluxes_wb()[0], sign * flux);
                near(
                    next.circuit().charges_c()[0].powi(2) / (2.0 * C)
                        + next.circuit().fluxes_wb()[0].powi(2) / (2.0 * L)
                        + CP[0] * theta[0]
                        + CP[1] * theta[1],
                    9.0 / (2.0 * C),
                );
                near(next.report().total_energy_residual_j, 0.0);
                if interval == 1 && sign == 1.0 {
                    currents.push(j);
                }
                old = theta;
                (thermal, actual, _) = next.into_parts();
            }
        }
    }
    assert!(currents[0] > currents[1]);
    assert!(currents[0] > currents[2]);
}

#[test]
fn contact_owner_source_budget_and_late_cancel_refusals_preserve_all_histories() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], BETA[0], coordinate, false);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], BETA[1], coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let thermal = initial(&curves, 999.0, None);
    let context = point(300.0, coordinate);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let layout = network(&cells);
    let surface = interface(&cards, 0.5, CONTACT_SLOPE, 1000.0);
    let surface_context = interface_context(300.0);
    let contacts = [binding(&cards, &surface, &surface_context, FRACTION)];
    let actual = contacted_circuit(&cards, &surface, 1.0);
    assert!(
        layout
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
            .is_err()
    );
    assert!(
        layout
            .step_response_candidate(
                &thermal,
                0.1,
                &contacted_problem(&cards, &surface),
                ConductionDrive::Current {
                    positive: &[5, 6],
                    negative: &[0, 3],
                    current_a: 1.0
                },
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    let mut budget = config();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_response_contact_circuit_candidate(
            &thermal,
            0.1,
            &actual,
            &[],
            &[],
            &contacts,
            limits(),
            field_controls(),
            controls(),
            budget,
            &[0.0; 2],
            || false,
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let run = |cancel: &mut dyn FnMut() -> bool| {
        layout.step_response_contact_circuit_candidate(
            &thermal,
            0.1,
            &actual,
            &[],
            &[],
            &contacts,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            cancel,
        )
    };
    let mut polls = 0;
    let good = run(&mut || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 1);
    let mut count = 0;
    assert!(matches!(
        run(&mut || {
            count += 1;
            count >= polls
        }),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let replay = run(&mut || false).unwrap();
    assert_eq!(replay.circuit().charges_c(), good.circuit().charges_c());
    assert_eq!(replay.circuit().fluxes_wb(), good.circuit().fluxes_wb());
    assert_eq!(replay.circuit().contacts(), good.circuit().contacts());
    assert_eq!(actual.time_s(), 0.0);
    assert_eq!(actual.charges_c(), &[3.0]);
    assert_eq!(actual.fluxes_wb(), &[0.0]);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(30)).unwrap(), 300.0);

    let narrow = interface(&cards, 0.5, CONTACT_SLOPE, 300.00001);
    let narrow_contacts = [binding(&cards, &narrow, &surface_context, FRACTION)];
    let narrow_circuit = contacted_circuit(&cards, &narrow, 1.0);
    assert!(
        layout
            .step_response_contact_circuit_candidate(
                &thermal,
                0.1,
                &narrow_circuit,
                &[],
                &[],
                &narrow_contacts,
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    assert_eq!(narrow_circuit.time_s(), 0.0);
    assert_eq!(narrow_circuit.charges_c(), &[3.0]);
    assert_eq!(
        narrow_circuit.contacts()[0]
            .source
            .resolved()
            .state_coordinate("Tm"),
        Some(300.0)
    );
}

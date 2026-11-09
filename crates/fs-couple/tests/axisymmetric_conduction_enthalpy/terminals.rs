//! G1/G3/G4: simultaneous ideal electrodes in the original thermal residual.
//! Manufactured body balances provide independent current/work references.

use super::contacts::{binding, contact_mesh, interface, interface_context, specimen};
use super::*;
use fs_couple::circuit_enthalpy::conduction::axisymmetric::terminals::{
    AxisymmetricTerminalControl, AxisymmetricTerminalDrive,
};
use fs_em::conduction::axisymmetric::{
    AxisymmetricConductanceControls, AxisymmetricMixedConductanceControls,
};

const TERMINALS: [&[usize]; 3] = [&[0, 1], &[2, 3, 4, 5], &[6, 7]];

fn response_controls() -> AxisymmetricConductanceControls {
    AxisymmetricConductanceControls {
        maximum_terminals: 4,
        maximum_field_solves: 3,
        maximum_response_values: 512,
        maximum_dense_work: 2048,
        reference_voltage_v: 1.0,
        coefficient_tolerance_s: 1e-8,
    }
}
fn mixed_controls() -> AxisymmetricMixedConductanceControls {
    AxisymmetricMixedConductanceControls {
        maximum_current_terminals: 3,
        maximum_response_values: 256,
        maximum_dense_work: 256,
    }
}
fn floating_specimen(cards: &[MaterialCard; 2], angle: f64) -> AxisymmetricConduction {
    let conductors = [0, 0, 1, 1]
        .map(|i| {
            resolve_conductor_state_with_basis(
                &cards[i],
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap()
        })
        .to_vec();
    // The mesh pieces connect electrically ONLY through their shared ideal
    // middle electrode. Its total current is zero, while each face conducts.
    AxisymmetricConduction::new(contact_mesh(angle), conductors, limits(), &mut || false).unwrap()
}

#[test]
fn floating_electrode_transfers_current_with_refreshed_material_and_fourier_heat() {
    let alpha = 0.002;
    let (a, ca) = material(1.0, 1.0, alpha);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let network = network(&cells);
    let dt = 0.5;
    // C1=A, C2=2*A, G=A, j=2; divide both BE balances by swept area A.
    // d1=.5*[4*(1+alpha*d1)+d2-d1], 2*d2=.5*[12+d1-d2].
    let a11 = 1.0 + dt - 4.0 * dt * alpha;
    let a22 = 2.0 + dt;
    let det = a11 * a22 - dt * dt;
    let d1 = (2.0 * a22 + dt * 6.0) / det;
    let d2 = (a11 * 6.0 + dt * 2.0) / det;
    let rho_a = 1.0 + alpha * d1;
    let middle_v = 2.0 * rho_a;
    let top_v = middle_v + 6.0;
    for angle in [0.7, 2.0, core::f64::consts::TAU] {
        let area = angle / 2.0;
        let thermal = initial(&curves, angle);
        let original = floating_specimen(&cards, angle);
        let currents = [-2.0 * area, 0.0, 2.0 * area];
        let voltage = [0.0, middle_v, top_v];
        // Reversed current-control order and a huge common voltage reference
        // must not change physical drops, local heat or source work.
        let mixed_currents = [2.0 * area, 0.0];
        let endpoint_voltages = [0.0, top_v];
        for control in [
            AxisymmetricTerminalControl::Currents(&currents),
            AxisymmetricTerminalControl::Mixed {
                current_terminals: &[2, 1],
                currents_a: &mixed_currents,
                voltages_v: &[1e250],
                controls: mixed_controls(),
            },
            AxisymmetricTerminalControl::Mixed {
                current_terminals: &[1],
                currents_a: &[0.0],
                voltages_v: &endpoint_voltages,
                controls: mixed_controls(),
            },
            AxisymmetricTerminalControl::Voltages(&voltage),
        ] {
            let next = network
                .step_axisymmetric_terminal_candidate(
                    &thermal,
                    dt,
                    &original,
                    AxisymmetricTerminalDrive {
                        terminals: &TERMINALS,
                        control,
                        response: response_controls(),
                    },
                    &[],
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + d1,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + d2,
            );
            let step = &next.report().electrical;
            near(step.voltage_v[2], middle_v);
            near(step.voltage_v[4], middle_v);
            near(step.voltage_v[6], top_v);
            let reaction = |vertices: &[usize]| {
                vertices
                    .iter()
                    .map(|&v| step.nodal_boundary_current_a[v].unwrap())
                    .sum::<f64>()
            };
            near(reaction(&[2, 3]), 2.0 * area);
            near(reaction(&[4, 5]), -2.0 * area);
            near(reaction(TERMINALS[1]), 0.0);
            for (triangle, &power) in step.cell_power_w.iter().enumerate() {
                let rho = if triangle < 2 { rho_a } else { 3.0 };
                near(
                    power,
                    4.0 * rho * original.mesh().face_measure(triangle).unwrap(),
                );
                near(step.current_density_a_m2[triangle][0], 0.0);
                near(step.current_density_a_m2[triangle][1], -2.0);
            }
            near(step.supplied_work_j, dt * 2.0 * area * top_v);
            near(step.heat_j, area * (d1 + 2.0 * d2));
            near(next.report().deposited_joule_j, step.heat_j);
            near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
            assert!(
                next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j
            );
            assert!(next.report().circuit_trials > 2);
            assert_eq!(next.report().circuit_trials % 2, 0);
            assert!(step.contacts.is_empty());
            let continued = network
                .step_axisymmetric_terminal_candidate(
                    next.thermal(),
                    2.0 * dt,
                    next.problem(),
                    AxisymmetricTerminalDrive {
                        terminals: &TERMINALS,
                        control: AxisymmetricTerminalControl::Currents(&currents),
                        response: response_controls(),
                    },
                    &[],
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            assert!(continued.report().electrical.voltage_v[6] > top_v);
        }
        // Change which terminal fixes the reporting gauge while keeping the
        // same physical current drive, including a floating reference terminal.
        let reordered = [TERMINALS[1], TERMINALS[2], TERMINALS[0]];
        let next = network
            .step_axisymmetric_terminal_candidate(
                &thermal,
                dt,
                &original,
                AxisymmetricTerminalDrive {
                    terminals: &reordered,
                    control: AxisymmetricTerminalControl::Currents(&[0.0, 2.0 * area, -2.0 * area]),
                    response: response_controls(),
                },
                &[],
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        near(next.report().electrical.voltage_v[2], 0.0);
        near(next.report().electrical.voltage_v[0], -middle_v);
        near(next.report().electrical.voltage_v[6], 6.0);
        near(
            next.report().electrical.supplied_work_j,
            area * (d1 + 2.0 * d2),
        );
        near(
            next.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + d1,
        );
        near(
            next.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + d2,
        );
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn terminal_matrix_feedback_preserves_explicit_contact_heat_and_source_refresh() {
    let angle = 2.0;
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, angle);
    let interface = interface(&cards, 0.5, 0.003, 1000.0);
    let original = specimen(&cards, &interface, angle);
    let context = point(300.0);
    let surface = interface_context(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let network = network(&cells);
    let bindings = [binding(&cards, &interface, &surface, 0.25)];
    let direct = network
        .step_axisymmetric_contact_candidate(
            &thermal,
            0.5,
            &original,
            ConductionDrive::Current {
                positive: &[6, 7],
                negative: &[0, 1],
                current_a: 2.0,
            },
            &bindings,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let terminals: [&[usize]; 2] = [&[0, 1], &[6, 7]];
    let endpoint = [0.0, direct.report().electrical.voltage_v[6]];
    for control in [
        AxisymmetricTerminalControl::Currents(&[-2.0, 2.0]),
        AxisymmetricTerminalControl::Mixed {
            current_terminals: &[1],
            currents_a: &[2.0],
            voltages_v: &[1e250],
            controls: mixed_controls(),
        },
        AxisymmetricTerminalControl::Voltages(&endpoint),
    ] {
        let next = network
            .step_axisymmetric_terminal_candidate(
                &thermal,
                0.5,
                &original,
                AxisymmetricTerminalDrive {
                    terminals: &terminals,
                    control,
                    response: response_controls(),
                },
                &bindings,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        for cell in [CellId(30), CellId(10)] {
            near(
                next.thermal().temperature_k(cell).unwrap(),
                direct.thermal().temperature_k(cell).unwrap(),
            );
        }
        let actual = &next.report().electrical;
        let expected = &direct.report().electrical;
        near(actual.supplied_work_j, expected.supplied_work_j);
        near(
            actual.contacts[0].response.power_w,
            expected.contacts[0].response.power_w,
        );
        near(
            actual.contacts[0].response.current_a,
            expected.contacts[0].response.current_a,
        );
        assert_eq!(actual.contacts[0].cells, expected.contacts[0].cells);
        near(actual.contacts[0].heat_fraction_to_first, 0.25);
        let refreshed = next.problem().contacts()[0].source.resolved();
        let target = direct.problem().contacts()[0].source.resolved();
        for axis in ["Ta", "Tb", "Tm", "pressure"] {
            near(
                refreshed.state_coordinate(axis).unwrap(),
                target.state_coordinate(axis).unwrap(),
            );
        }
        assert_eq!(refreshed.card_identity(), target.card_identity());
        near(next.report().deposited_joule_j, expected.heat_j);
    }
    assert!(
        network
            .step_axisymmetric_terminal_candidate(
                &thermal,
                0.5,
                &original,
                AxisymmetricTerminalDrive {
                    terminals: &terminals,
                    control: AxisymmetricTerminalControl::Currents(&[-2.0, 2.0]),
                    response: response_controls(),
                },
                &[],
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn terminal_basis_budgets_and_late_refusals_leave_the_original_owners_reusable() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 2.0);
    let original = floating_specimen(&cards, 2.0);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let network = network(&cells);
    let drive = AxisymmetricTerminalDrive {
        terminals: &TERMINALS,
        control: AxisymmetricTerminalControl::Currents(&[-2.0, 0.0, 2.0]),
        response: response_controls(),
    };
    let run = |drive, controls, cancel: &mut dyn FnMut() -> bool| {
        network.step_axisymmetric_terminal_candidate(
            &thermal,
            0.5,
            &original,
            drive,
            &[],
            limits(),
            field_controls(),
            controls,
            &[0.0; 2],
            cancel,
        )
    };
    let energy = thermal.energy_j();
    let mut polls = 0;
    let accepted = run(drive, config(), &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut too_few = config();
    too_few.maximum_circuit_trials = 1;
    assert!(matches!(
        run(drive, too_few, &mut || false),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    too_few.maximum_circuit_trials = accepted.report().circuit_trials - 1;
    assert!(matches!(
        run(drive, too_few, &mut || false),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut late = 0;
    assert!(matches!(
        run(drive, config(), &mut || {
            late += 1;
            late >= polls - 1
        }),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut unbalanced = drive;
    unbalanced.control = AxisymmetricTerminalControl::Currents(&[-2.0, 0.0, 3.0]);
    assert!(run(unbalanced, config(), &mut || false).is_err());
    let mut repeated = drive;
    repeated.control = AxisymmetricTerminalControl::Mixed {
        current_terminals: &[2, 2],
        currents_a: &[2.0, 0.0],
        voltages_v: &[0.0],
        controls: mixed_controls(),
    };
    assert!(run(repeated, config(), &mut || false).is_err());
    let mut small_basis = drive;
    small_basis.response.maximum_field_solves = 1;
    assert!(run(small_basis, config(), &mut || false).is_err());
    let mut exact = config();
    exact.maximum_circuit_trials = accepted.report().circuit_trials;
    let retry = run(drive, exact, &mut || false).unwrap();
    assert_eq!(
        retry.report().electrical.voltage_v,
        accepted.report().electrical.voltage_v
    );
    assert_eq!(
        retry.thermal().energy_j().to_bits(),
        accepted.thermal().energy_j().to_bits()
    );
    assert_eq!(thermal.energy_j().to_bits(), energy.to_bits());
    assert_eq!(thermal.time_s(), 0.0);
    for source in original.conductors() {
        near(source.material().query().axes()["T"], 300.0);
    }
}

//! G1/G3: independent midpoint RC equations with separate bulk and contact laws.
use super::*;
use fs_couple::material_circuit::conductivity::contact::{
    ContactAxis, ContactCoordinate, ContactDirection, ContactResponse,
};

const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);
const CURRENT: Dims = Dims([0, 0, 0, 0, 1, 0]);

fn root(voltage: f64, bulk_slope: f64, contacts: impl Fn(f64) -> f64) -> f64 {
    let (mut lo, mut hi) = (0.0, 4.0);
    for _ in 0..80 {
        let i = 0.5 * (lo + hi);
        let bulk_voltage = 2.0 * i / (1.0 + (1.0 + 4.0 * bulk_slope * i).sqrt());
        if bulk_voltage + i * (contacts(i) + DT / (2.0 * CAPACITANCE)) > voltage {
            hi = i;
        } else { lo = i; }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_contact_voltage_current_and_density_sources_close_with_bulk_at_the_same_midpoint() {
    let bulk_slope = 0.2;
    let source = card(bulk_slope, false, 8.0, 16.0);
    let area = 0.5;
    let r0 = 0.4 / area;
    for (coordinate, dims) in [(ContactCoordinate::VoltageMagnitude, VOLTAGE),
        (ContactCoordinate::Current, CURRENT), (ContactCoordinate::CurrentDensityMagnitude, J)] {
        let axes = [ContactAxis { name: "drive", coordinate }];
        for direction in [ContactDirection::AlongBranch, ContactDirection::AgainstBranch] {
            for sign in [-1.0, 1.0] {
                let mut endpoints = Vec::new();
                for beta in [0.0, 0.04] {
                    let interface = contact_source::response_card(&source, &source, "drive", dims,
                        0.4, 0.0, beta, 10.0);
                    let contact = contact_source::response_contact(7, 1, &interface, area, "drive", dims);
                    let original = circuit(&source, sign * 2.0)
                        .rebind_contacts_candidate(vec![contact], &[]).unwrap();
                    let responses = [ContactResponse { id: 7, card: &interface, direction, axes: &axes }];
                    let result = original.conductivity_with_contact_response_candidate(
                        original.resistors().to_vec(), original.contacts().to_vec(), &bindings(&source),
                        Some(&responses), &[], DT, controls(), solve(), || false,
                    ).unwrap();
                    let orientation = if direction == ContactDirection::AlongBranch { sign } else { -sign };
                    let resistance = |i: f64| match coordinate {
                        ContactCoordinate::VoltageMagnitude => r0 / (1.0 - r0 * beta * i),
                        ContactCoordinate::Current => r0 * (1.0 + beta * orientation * i),
                        ContactCoordinate::CurrentDensityMagnitude => r0 * (1.0 + beta * i / area),
                        _ => unreachable!(),
                    };
                    let i = root(2.0, bulk_slope, resistance);
                    let report = result.step();
                    let sample = &report.evaluated_contacts[0];
                    let native = &report.electrical.contacts[0];
                    near(native.current_a, sign * i, 5e-11);
                    near(native.voltage_v, sign * i * resistance(i), 5e-11);
                    near(sample.point.current_a, orientation * i, 5e-11);
                    near(sample.point.current_density_a_m2, orientation * i / area, 5e-11);
                    near(sample.point.voltage_v, orientation * i * resistance(i), 5e-11);
                    near(sample.contact.resistance_ohm(), resistance(i), 5e-11);
                    assert_eq!(sample.point.direction, direction);
                    assert_eq!(sample.contact.source().resolved().card_identity(), interface.content_hash());
                    assert!(sample.resistance_defect_ohm.abs() < 1e-11);
                    assert!(sample.power_defect_w.abs() < controls().power_tolerance_w);
                    near(report.electrical.heat_j,
                        report.electrical.resistors[0].heat_j + native.heat_j, 5e-11);
                    near(report.electrical.stored_after_j + report.electrical.heat_j, 1.0, 5e-11);
                    near(result.circuit().charges_c()[0], sign * (1.0 - DT * i), 5e-11);
                    assert_eq!(report.source_evaluations, 2 * report.iterations);
                    endpoints.push(result.circuit().charges_c()[0]);
                    let next = result.circuit();
                    let continued = next.conductivity_with_contact_response_candidate(
                        next.resistors().to_vec(), next.contacts().to_vec(), &bindings(&source),
                        Some(&responses), &[], DT, controls(), solve(), || false,
                    ).unwrap();
                    assert!(continued.circuit().charges_c()[0].abs() < next.charges_c()[0].abs());
                }
                assert!((endpoints[0] - endpoints[1]).abs() > 1e-5,
                    "a source-only contact-law change alters actual charge evolution");
            }
        }
    }
}

#[test]
fn g1_each_series_contact_uses_its_own_drop_and_area_without_combining_losses() {
    let source = card(0.0, false, 8.0, 16.0);
    let interfaces = [contact_source::response_card(&source, &source, "V", VOLTAGE, 0.3, 0.0, 0.04, 10.0),
        contact_source::response_card(&source, &source, "J", J, 0.4, 0.0, 0.03, 10.0)];
    let contacts = vec![contact_source::response_contact(9, 1, &interfaces[1], 0.5, "J", J),
        contact_source::response_contact(7, 1, &interfaces[0], 1.0, "V", VOLTAGE)];
    let original = circuit(&source, 2.0).rebind_contacts_candidate(contacts, &[]).unwrap();
    let voltage_axes = [ContactAxis { name: "V", coordinate: ContactCoordinate::Voltage }];
    let density_axes = [ContactAxis { name: "J", coordinate: ContactCoordinate::CurrentDensity }];
    let responses = [ContactResponse { id: 9, card: &interfaces[1], direction: ContactDirection::AlongBranch, axes: &density_axes },
        ContactResponse { id: 7, card: &interfaces[0], direction: ContactDirection::AlongBranch, axes: &voltage_axes }];
    let result = original.conductivity_with_contact_response_candidate(
        original.resistors().to_vec(), original.contacts().to_vec(), &bindings(&source),
        Some(&responses), &[], DT, controls(), solve(), || false,
    ).unwrap();
    let resistances = |i: f64| [0.3 / (1.0 - 0.3 * 0.04 * i), 0.8 * (1.0 + 0.03 * i / 0.5)];
    let i = root(2.0, 0.0, |i| resistances(i).iter().sum());
    for ((sample, native), expected) in result.step().evaluated_contacts.iter()
        .zip(&result.step().electrical.contacts).zip(resistances(i)) {
        near(sample.contact.resistance_ohm(), expected, 5e-11);
        near(native.voltage_v, i * expected, 5e-11);
        near(native.heat_j, DT * i * i * expected, 5e-11);
        assert_eq!(sample.point.id, native.contact.id());
    }
    near(result.step().electrical.stored_after_j + result.step().electrical.heat_j, 1.0, 5e-11);
}

#[test]
fn g4_contact_source_binding_limits_budget_and_late_cancel_preserve_accepted_state() {
    let source = card(0.2, false, 8.0, 16.0);
    let interface = contact_source::response_card(&source, &source, "V", VOLTAGE, 0.4, 0.0, 0.04, 10.0);
    let original = circuit(&source, 2.0).rebind_contacts_candidate(vec![
        contact_source::response_contact(7, 1, &interface, 1.0, "V", VOLTAGE),
    ], &[]).unwrap();
    let saved = original.clone();
    let axes = [ContactAxis { name: "V", coordinate: ContactCoordinate::Voltage }];
    let run = |card: &fs_matdb::InterfaceSystemCard, bindings: &[ContactAxis<'_>], controls,
        cancel: &mut dyn FnMut() -> bool| {
        let responses = [ContactResponse { id: 7, card, direction: ContactDirection::AlongBranch, axes: bindings }];
        original.conductivity_with_contact_response_candidate(original.resistors().to_vec(),
            original.contacts().to_vec(), &super::bindings(&source), Some(&responses), &[], DT,
            super::controls(), controls, cancel)
    };
    assert!(run(&interface, &[], solve(), &mut || false).is_err());
    assert!(run(&interface, &[ContactAxis { name: "V", coordinate: ContactCoordinate::Current }],
        solve(), &mut || false).is_err());
    assert!(run(&interface, &[axes[0], axes[0]], solve(), &mut || false).is_err());
    let wrong = contact_source::response_card(&source, &source, "V", VOLTAGE, 0.6, 0.0, 0.04, 10.0);
    assert!(run(&wrong, &axes, solve(), &mut || false).is_err());
    assert!(matches!(run(&interface, &axes, ConductivitySolve { maximum_iterations: 1, ..solve() },
        &mut || false), Err(MaterialCircuitError::ConductivityConvergence { .. })));
    let mut polls = 0;
    run(&interface, &axes, solve(), &mut || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(matches!(run(&interface, &axes, solve(), &mut || { count += 1; count == polls }),
        Err(MaterialCircuitError::Cancelled)));
    assert_eq!(original, saved);

    let limited = contact_source::response_card(&source, &source, "V", VOLTAGE, 0.4, 0.0, 0.04, 0.01);
    let bounded = circuit(&source, 2.0).rebind_contacts_candidate(vec![
        contact_source::response_contact(7, 1, &limited, 1.0, "V", VOLTAGE),
    ], &[]).unwrap();
    let before = bounded.clone();
    let response = [ContactResponse { id: 7, card: &limited, direction: ContactDirection::AlongBranch, axes: &axes }];
    assert!(bounded.conductivity_with_contact_response_candidate(bounded.resistors().to_vec(),
        bounded.contacts().to_vec(), &bindings(&source), Some(&response), &[], DT,
        controls(), solve(), || false).is_err());
    assert_eq!(bounded, before);
}

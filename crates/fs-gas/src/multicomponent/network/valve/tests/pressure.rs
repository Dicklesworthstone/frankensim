//! Differential-pressure actuation through the original nozzle and storage owners.
use super::*;

fn law(cracking: f64, full_open: f64) -> PressureOpening {
    PressureOpening::new(Pressure::new(cracking), Pressure::new(full_open)).unwrap()
}

#[test]
fn pressure_actuated_aperture_has_exact_thresholds_and_rejects_invalid_laws() {
    let opening = law(0.25, 1.25);
    assert_eq!(opening.cracking_pressure(), Pressure::new(0.25));
    assert_eq!(opening.full_open_pressure(), Pressure::new(1.25));
    for (difference, expected) in [
        (-2.0, 0.0),
        (0.0, 0.0),
        (0.25, 0.0),
        (0.5, 0.25),
        (0.75, 0.5),
        (1.25, 1.0),
        (2.0, 1.0),
        (0.25, 0.0),
    ] {
        assert_eq!(
            opening.fraction(Pressure::new(difference)).unwrap(),
            expected
        );
    }
    for (cracking, full_open) in [
        (-1.0, 1.0),
        (0.0, 0.0),
        (1.0, 0.5),
        (1.0, 1.0),
        (f64::NAN, 1.0),
        (0.0, f64::NAN),
        (0.0, f64::INFINITY),
        (f64::NEG_INFINITY, 1.0),
    ] {
        assert!(PressureOpening::new(Pressure::new(cracking), Pressure::new(full_open)).is_err());
    }
    for difference in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(opening.fraction(Pressure::new(difference)).is_err());
    }
    assert!(
        law(0.0, f64::MAX)
            .fraction(Pressure::new(f64::from_bits(1)))
            .is_err()
    );
}

#[test]
fn pressure_actuated_flow_matches_choked_and_unchoked_nozzle_with_exact_reverse_closure() {
    let gas = caloric_gas(false, 0.0);
    let high = gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap();
    let opening = law(0.25, 1.25);
    let valve = Orifice::new(Area::new(0.1), 0.7, 0.6, ValveDirection::Bidirectional)
        .unwrap()
        .with_pressure_opening(opening);
    assert_eq!(valve.pressure_opening(), Some(opening));
    assert_eq!(valve.opening(), 0.6);
    close(valve.effective_area().value(), 0.042, 1e-17);
    for (backpressure, fraction) in [(1.9, 0.0), (1.5, 0.25), (1.0, 0.75), (0.2, 1.0)] {
        let low = gas
            .conserved(&[0.0, backpressure], 0.0, backpressure)
            .unwrap();
        let forward = valve
            .flow(&gas, high.as_cell(), low.as_cell(), || false)
            .unwrap();
        let reverse = valve
            .flow(&gas, low.as_cell(), high.as_cell(), || false)
            .unwrap();
        if fraction == 0.0 {
            assert_eq!(forward, ValveFlow::zero(2));
            assert_eq!(reverse, ValveFlow::zero(2));
            continue;
        }
        // The threshold changes area only. The actual pressure ratio still
        // sets the original compressible nozzle's sonic transition and energy.
        close(
            forward.mass_kg_s,
            gamma_rate(0.042 * fraction, 2.0, 2.0, backpressure),
            2e-11,
        );
        close(forward.energy_w / forward.mass_kg_s, 7.0 - 0.0025, 2e-12);
        assert_eq!(forward.species_kg_s[1], 0.0);
        assert_eq!(forward.choked, backpressure <= 1.0);
        assert_eq!(reverse.mass_kg_s, -forward.mass_kg_s);
        assert_eq!(reverse.energy_w, -forward.energy_w);
        for direction in [ValveDirection::FirstToSecond, ValveDirection::SecondToFirst] {
            let check = Orifice::new(Area::new(0.1), 0.7, 0.6, direction)
                .unwrap()
                .with_pressure_opening(opening);
            let (allowed, refused) = if direction == ValveDirection::FirstToSecond {
                (
                    [high.as_cell(), low.as_cell()],
                    [low.as_cell(), high.as_cell()],
                )
            } else {
                (
                    [low.as_cell(), high.as_cell()],
                    [high.as_cell(), low.as_cell()],
                )
            };
            close(
                check
                    .flow(&gas, allowed[0], allowed[1], || false)
                    .unwrap()
                    .mass_kg_s
                    .abs(),
                forward.mass_kg_s,
                2e-14,
            );
            assert_eq!(
                check.flow(&gas, refused[0], refused[1], || false).unwrap(),
                ValveFlow::zero(2)
            );
        }
    }
    // Exact closure at the source-decoded cracking differential, independently
    // of the temperature inversion's final floating-point pressure rounding.
    let low = gas.conserved(&[0.0, 1.75], 0.0, 1.75).unwrap();
    let difference = gas.primitive(high.as_cell(), || false).unwrap().pressure
        - gas.primitive(low.as_cell(), || false).unwrap().pressure;
    let at_crack = valve.with_pressure_opening(law(difference, difference + 1.0));
    assert_eq!(
        at_crack
            .flow(&gas, high.as_cell(), low.as_cell(), || false)
            .unwrap(),
        ValveFlow::zero(2)
    );
}

#[test]
fn pressure_actuated_flow_preserves_admission_and_refuses_area_underflow() {
    let gas = caloric_gas(false, 0.0);
    let high = gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap();
    let low = gas.conserved(&[0.0, 1.5], 0.0, 1.5).unwrap();
    let tiny = orifice(f64::from_bits(1)).with_pressure_opening(law(0.25, 1.25));
    assert_eq!(
        tiny.flow(&gas, high.as_cell(), low.as_cell(), || false),
        Err(ValveError::Invalid(
            "nonzero pressure-actuated orifice area underflowed"
        ))
    );
    let closed = orifice(0.1).with_pressure_opening(law(1.0, 2.0));
    assert_eq!(
        closed.flow(&gas, high.as_cell(), low.as_cell(), || true),
        Err(ValveError::Cancelled)
    );
    let moving = gas.conserved(&[1.0, 0.0], 0.1, 2.0).unwrap();
    assert!(
        closed
            .flow(&gas, moving.as_cell(), low.as_cell(), || false)
            .is_err()
    );
}

#[test]
fn pressure_actuated_finite_transfer_obeys_blowdown_endpoint_and_cracking_bound() {
    let gas = caloric_gas(false, 0.0);
    let initial = tanks(&gas);
    let valve = orifice(0.1).with_pressure_opening(law(0.4, 1.2));
    // Both finite tanks contain pure A, so donor/receiver pressure follows
    // independently from gamma-law storage at their different fixed volumes.
    let limiting_mass = 0.2 * (1.0 - (4.5_f64 / 6.0).powf(1.0 / 1.4));
    let mut previous = 0.0;
    for duration in [0.001, 0.1, 10.0, 100.0] {
        let update = pair::solve(
            &gas,
            [&initial[0], &initial[1]],
            valve,
            duration,
            ValveControl::default(),
            &mut || false,
        )
        .unwrap();
        let mass = update.exchange.mass_to_second_kg();
        assert!(mass > previous && mass < limiting_mass);
        previous = mass;
        let remaining = 1.0 - mass / 0.2;
        let donor = gas
            .primitive(cell(&update.states[0][..4], 2), || false)
            .unwrap();
        let receiver = gas
            .primitive(cell(&update.states[1][..4], 2), || false)
            .unwrap();
        let donor_pressure = 2.0 * remaining.powf(1.4);
        let receiver_pressure = 0.1 + 4.0 * (1.0 - remaining.powf(1.4));
        close(donor.pressure, donor_pressure, 3e-11);
        close(receiver.pressure, receiver_pressure, 3e-11);
        let difference = donor.pressure - receiver.pressure;
        assert!(difference > 0.4);
        let fraction = ((difference - 0.4) / 0.8).min(1.0);
        let endpoint = gamma_rate(
            0.1 * fraction,
            donor.pressure,
            donor.temperature,
            receiver.pressure,
        );
        close(mass, duration * endpoint, 3e-10);
        close(
            update.exchange.energy_to_second_j + 0.0025 * mass,
            1.0 - remaining.powf(1.4),
            3e-11,
        );
        let mut energy_change = 0.0;
        for side in 0..2 {
            let next = cell(&update.states[side][..4], 2);
            let sign = if side == 0 { -1.0 } else { 1.0 };
            for k in 0..2 {
                close(
                    initial[side].volume_m3()
                        * (next.partial_densities[k] - initial[side].state().partial_densities[k]),
                    sign * update.exchange.species_to_second_kg[k],
                    3e-14,
                );
            }
            energy_change += initial[side].volume_m3()
                * (gas.absolute_energy_density(next).unwrap()
                    - gas.absolute_energy_density(initial[side].state()).unwrap());
        }
        close(energy_change, 0.0, 3e-13);
    }
    // Even a deliberately loose mass tolerance cannot accept a trial that has
    // already crossed below cracking pressure and has no endpoint flow.
    let valve = orifice(0.1).with_pressure_opening(law(1.2, 1.8));
    let loose = ValveControl {
        absolute_mass_tolerance_kg: 1.0,
        ..ValveControl::default()
    };
    let update = pair::solve(
        &gas,
        [&initial[0], &initial[1]],
        valve,
        10.0,
        loose,
        &mut || false,
    )
    .unwrap();
    let pressure = |side: usize| {
        gas.primitive(cell(&update.states[side][..4], 2), || false)
            .unwrap()
            .pressure
    };
    assert!(pressure(0) - pressure(1) > 1.2);
}

fn source_owner() -> (MixtureNetwork1d, NetworkValveOperator) {
    let gas = caloric_gas(true, 0.0);
    let network = network(&gas);
    let valve = Orifice::new(Area::new(0.01), 0.8, 1.0, ValveDirection::FirstToSecond)
        .unwrap()
        .with_pressure_opening(law(0.2, 1.2));
    let operator = NetworkValveOperator::new(
        &network,
        vec![ValveLink {
            first: 0,
            second: 1,
            orifice: valve,
        }],
        || false,
    )
    .unwrap();
    (network, operator)
}

#[test]
fn pressure_actuated_network_source_updates_real_species_energy_and_rolls_back_refusal() {
    let (mut reference, mut reference_operator) = source_owner();
    let old = reference.plenums().to_vec();
    let old_time = reference.time();
    let mut polls = 0;
    let report = reference_operator
        .advance_for(
            &mut reference,
            Time::new(0.02),
            ValveSolve::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let exchange = reference_operator.accepted_exchanges().unwrap()[0];
    assert!(exchange.mass_to_second_kg() > 0.0);
    assert!(reference.plenums()[1].state().partial_densities[0] > 0.0);
    assert_eq!(reference.time(), old_time);
    for side in 0..2 {
        let new = &reference.plenums()[side];
        let sign = if side == 0 { -1.0 } else { 1.0 };
        for k in 0..2 {
            close(
                new.volume_m3()
                    * (new.state().partial_densities[k] - old[side].state().partial_densities[k]),
                sign * exchange.species_to_second_kg[k],
                3e-14,
            );
        }
        close(
            new.volume_m3() * (new.state().energy - old[side].state().energy),
            sign * exchange.energy_to_second_j,
            3e-13,
        );
    }
    close(report.absolute_energy_residual_j, 0.0, 3e-12);
    for stop in [1, polls / 2, polls] {
        let (mut network, mut operator) = source_owner();
        let original = network.plenums().to_vec();
        let mut count = 0;
        assert_eq!(
            operator.advance_for(&mut network, Time::new(0.02), ValveSolve::default(), || {
                count += 1;
                count >= stop
            }),
            Err(ValveError::Cancelled)
        );
        assert_eq!(network.plenums(), original);
        assert_eq!(operator.last_report(), None);
        assert_eq!(operator.accepted_exchanges(), None);
        let limited = ValveSolve {
            maximum_iterations: 1,
            ..ValveSolve::default()
        };
        assert!(
            operator
                .advance_for(&mut network, Time::new(0.02), limited, || false)
                .is_err()
        );
        assert_eq!(network.plenums(), original);
        operator
            .advance_for(&mut network, Time::new(0.02), ValveSolve::default(), || {
                false
            })
            .unwrap();
        assert_eq!(network.plenums(), reference.plenums());
        assert_eq!(
            operator.accepted_exchanges(),
            reference_operator.accepted_exchanges()
        );
    }
}

#[test]
fn pressure_actuated_valve_drives_original_network_waves_and_shared_energy_balance() {
    let (network, operator) = source_owner();
    let initial = network.totals();
    let mut owner = ValvedNetwork1d::new(network, operator.links().to_vec(), || false).unwrap();
    for _ in 0..4 {
        let report = owner
            .step(Time::new(0.002), ValveControl::default(), || false)
            .unwrap();
        assert!(owner.accepted_valve_exchanges().unwrap()[0].mass_to_second_kg() > 0.0);
        close(report.absolute_energy_residual_j, 0.0, 4e-12);
    }
    assert!(owner.plenums()[1].state().partial_densities[0] > 0.0);
    assert!(owner.branches()[1].totals().species_masses[0] > 0.0);
    for k in 0..2 {
        close(
            owner.totals().species_masses[k],
            initial.species_masses[k],
            3e-12,
        );
    }
    close(
        owner.totals().absolute_energy(owner.gas()).unwrap(),
        initial.absolute_energy(owner.gas()).unwrap(),
        3e-11,
    );
}

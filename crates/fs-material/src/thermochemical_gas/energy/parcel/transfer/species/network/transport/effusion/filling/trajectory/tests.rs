//! Adaptive exact-vacuum startup through the original finite-pair owners.
use super::super::{
    step_candidate as first_fill,
    tests::{AREA, controls as physical, empty, interface, link, pure, solve},
};
use super::*;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn adaptive(tolerance: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: 2.0,
        minimum_step_s: 1e-9,
        maximum_step_s: 2.0,
        maximum_attempts: 12_000,
        maximum_accepted_intervals: 2_000,
        relative_tolerance: tolerance,
        temperature_tolerance_k: 1e-7,
        pressure_tolerance_pa: 1e-7,
        species_mass_tolerance_kg: tolerance * 1e-6,
        stream_energy_tolerance_j: tolerance,
    }
}

fn coordinates(donor: &ClosedGasParcel, receiver: &ClosedGasParcel) -> [f64; 6] {
    [donor, receiver]
        .map(|gas| {
            [
                gas.mass_kg(),
                gas.state().state().temperature,
                gas.state().state().pressure,
            ]
        })
        .concat()
        .try_into()
        .unwrap()
}

fn relative_error(a: [f64; 6], b: [f64; 6]) -> f64 {
    a.into_iter()
        .zip(b)
        .map(|(x, y)| (x / y - 1.0).abs())
        .fold(0.0, f64::max)
}

// The independent temporal reference uses the original first-fill operation
// exactly once, followed only by ordinary reciprocal native transactions.
fn native_reference(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    bound: EffusionFillLink<'_>,
    duration_s: f64,
    steps: u32,
) -> BTreeMap<u64, ClosedGasParcel> {
    let dt = duration_s / f64::from(steps);
    let (next_donor, next_receiver, _) = first_fill(
        donor,
        receiver,
        bound,
        dt,
        solve(),
        solve(),
        physical(),
        || false,
    )
    .unwrap()
    .into_parts();
    let EffusionReceiver::Filled(next_receiver) = next_receiver else {
        panic!("positive reference aperture must fill")
    };
    let mut current = BTreeMap::from([(bound.donor, next_donor), (bound.receiver, next_receiver)]);
    let native_link = [EffusionNetworkLink {
        id: bound.id,
        parcels: [bound.donor, bound.receiver],
        interface: bound.interface,
        area_m2: bound.area_m2,
    }];
    let paths = GasTransportPaths {
        effusion: &native_link,
        ..Default::default()
    };
    let endpoints = BTreeMap::from([(bound.donor, solve()), (bound.receiver, solve())]);
    for _ in 1..steps {
        current = reservoir::step_candidate(
            &current,
            paths,
            &[],
            &[],
            dt,
            &endpoints,
            physical(),
            || false,
        )
        .unwrap()
        .into_parts()
        .0;
    }
    current
}

#[test]
fn g1_exact_vacuum_startup_refines_with_immediate_backflow_and_complete_closed_balances() {
    let donor = pure(2.5, 0.0, 1e12);
    let saved_donor = donor.clone();
    let receiver = empty(0.002);
    let saved_receiver = receiver.clone();
    let law = interface(&donor, 1.0, 1e-8);
    let bound = link(&law, AREA);
    let duration = 2.0;
    let coarse_reference = native_reference(&donor, &receiver, bound, duration, 1024);
    let reference = native_reference(&donor, &receiver, bound, duration, 2048);
    let oracle = coordinates(&reference[&bound.donor], &reference[&bound.receiver]);
    let reference_difference = relative_error(
        coordinates(
            &coarse_reference[&bound.donor],
            &coarse_reference[&bound.receiver],
        ),
        oracle,
    );
    assert!(reference_difference < 2e-5, "{reference_difference}");
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [2e-3, 2e-4, 2e-5] {
        let result = advance_candidate(
            &donor,
            &receiver,
            bound,
            duration,
            solve(),
            solve(),
            physical(),
            adaptive(tolerance),
            || false,
        )
        .unwrap();
        let filled = result.receiver().gas().unwrap();
        errors.push(relative_error(coordinates(result.donor(), filled), oracle));
        let report = result.step();
        work.push(report.attempted_steps);
        let backflow = report.first_reciprocal.as_ref().unwrap();
        assert!(backflow.exchange.fluxes.iter().any(|row| {
            row.transfer.donor == bound.receiver
                && row.transfer.receiver == bound.donor
                && row.transfer.mass_kg > 0.0
        }));
        assert!(report.first_fill.duration_s < duration);
        assert!(report.continuation.is_some());
        assert_eq!(report.duration_s, duration);
        assert!(report.maximum_accepted_error <= 1.0);
        assert!(report.accepted_intervals > 1);
        assert!(report.attempted_steps >= 3 * report.accepted_intervals);
        let retained_duration = report.first_fill.duration_s
            + backflow.duration_s
            + report
                .continuation
                .as_ref()
                .map_or(0.0, |step| step.duration_s);
        assert!((retained_duration - duration).abs() < 1e-13);

        let nd = result.donor().mass_kg() / 0.028;
        let nr = filled.mass_kg() / 0.028;
        let n0 = donor.mass_kg() / 0.028;
        let td = result.donor().state().state().temperature;
        let tr = filled.state().state().temperature;
        // Independent constant-Cv specialization: no formation datum appears.
        let energy = 1.5 * R * (nd * (td - 400.0) + nr * (tr - 400.0));
        assert!(energy.abs() <= physical().exchange.energy_tolerance_j);
        assert!(
            (result.donor().mass_kg() + filled.mass_kg() - donor.mass_kg()).abs()
                <= physical().exchange.mass_tolerance_kg
        );
        let entropy = R
            * (nd * (1.5 * (td / 400.0).ln() - (nd / n0).ln())
                + nr * (1.5 * (tr / 400.0).ln()
                    - ((nr / receiver.volume_m3()) / (n0 / donor.volume_m3())).ln()));
        assert!(entropy > 0.0);
        assert!((report.balance.entropy_change_j_k - entropy).abs() < 1e-11);
        assert!(report.balance.energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);
        assert!(
            report.balance.maximum_species_mass_residual_kg
                <= physical().exchange.mass_tolerance_kg
        );
    }
    assert!(errors[2] < 5e-4, "{errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.8 * pair[0], "{errors:?}");
    }
    assert!(work[2] > work[0]);
    assert_eq!(donor, saved_donor);
    assert_eq!(receiver, saved_receiver);
}

#[test]
fn g0_inactive_adaptive_window_retains_an_exact_evacuated_receiver() {
    let donor = pure(2.5, 0.0, 0.0);
    let receiver = empty(0.002);
    for (transmission, area) in [(1.0, 0.0), (0.0, AREA)] {
        let law = interface(&donor, transmission, 1e-8);
        let controls = EffusionAdvanceControls {
            initial_step_s: 0.1,
            ..adaptive(1e-4)
        };
        let result = advance_candidate(
            &donor,
            &receiver,
            link(&law, area),
            3.0,
            solve(),
            solve(),
            physical(),
            controls,
            || false,
        )
        .unwrap();
        assert_eq!(result.donor(), &donor);
        assert_eq!(
            result.receiver(),
            &EffusionReceiver::Evacuated(receiver.clone())
        );
        let report = result.step();
        assert_eq!(report.first_fill.duration_s, 3.0);
        assert_eq!(report.attempted_steps, 2);
        assert_eq!(report.accepted_intervals, 1);
        assert_eq!(report.maximum_accepted_error, 0.0);
        assert!(report.first_reciprocal.is_none());
        assert!(report.continuation.is_none());
        assert_eq!(report.balance.energy_residual_j, 0.0);
        assert_eq!(report.balance.entropy_change_j_k, 0.0);
    }
}

#[test]
fn g3_adaptive_startup_limits_entropy_losing_overfill_before_reciprocal_continuation() {
    let donor = pure(2.5, 0.0, 0.0);
    let receiver = empty(1e-9);
    let law = interface(&donor, 1.0, 1e-8);
    let bound = link(&law, AREA);
    assert!(matches!(
        first_fill(
            &donor,
            &receiver,
            bound,
            0.001,
            solve(),
            solve(),
            physical(),
            || false
        ),
        Err(EffusionFillError::EntropyDecrease { .. })
    ));
    let result = advance_candidate(
        &donor,
        &receiver,
        bound,
        0.001,
        solve(),
        solve(),
        physical(),
        adaptive(2e-3),
        || false,
    )
    .unwrap();
    assert!(result.receiver().gas().is_some());
    // The physical initial-pressure ceiling may avoid the coarse overfill
    // before any evaluated trial needs rejection. The retained interval is
    // still genuinely shorter and continues through both finite gas owners.
    assert!(result.step().first_fill.duration_s < 0.0005);
    assert!(result.step().first_reciprocal.is_some());
    assert!(result.step().continuation.is_some());
    assert!(result.step().balance.entropy_change_j_k >= -physical().entropy_tolerance_j_k);
}

#[test]
fn g4_initial_pressure_ceiling_respects_the_authored_minimum_step_without_a_pressure_floor() {
    let donor = pure(2.5, 0.0, 0.0);
    let saved_donor = donor.clone();
    let receiver = empty(1e-9);
    let saved_receiver = receiver.clone();
    let law = interface(&donor, 1.0, 1e-8);
    let minimum = 1e-4;
    let error = advance_candidate(
        &donor,
        &receiver,
        link(&law, AREA),
        0.001,
        solve(),
        solve(),
        physical(),
        EffusionAdvanceControls {
            minimum_step_s: minimum,
            ..adaptive(2e-3)
        },
        || false,
    )
    .unwrap_err();
    let Error::InitialPressureLimit {
        maximum_step_s,
        minimum_step_s,
    } = error
    else {
        panic!("expected explicit initial-pressure timestep refusal: {error:?}")
    };
    assert_eq!(minimum_step_s, minimum);
    assert!(maximum_step_s > 0.0 && maximum_step_s < minimum);
    // Independent ideal-gas ceiling from the admitted initial molar beam and
    // the authored receiver temperature bracket. No endpoint pressure is set.
    let molar_rate = AREA * 200.0 / (2.0 * core::f64::consts::PI * 0.028 * R * 400.0).sqrt();
    let expected = 200.0 * receiver.volume_m3() / (molar_rate * R * solve().maximum_temperature_k);
    assert!((maximum_step_s / expected - 1.0).abs() < 1e-12);
    assert_eq!(donor, saved_donor);
    assert_eq!(receiver, saved_receiver);
}

#[test]
fn g4_initial_and_continuation_work_limits_and_late_cancellation_publish_neither_owner() {
    let donor = pure(2.5, 0.0, 0.0);
    let saved_donor = donor.clone();
    let receiver = empty(0.002);
    let saved_receiver = receiver.clone();
    let law = interface(&donor, 1.0, 1e-8);
    let bound = link(&law, AREA);
    assert!(matches!(
        advance_candidate(
            &donor,
            &receiver,
            bound,
            2.0,
            solve(),
            solve(),
            physical(),
            EffusionAdvanceControls {
                maximum_attempts: 3,
                ..adaptive(1e-6)
            },
            || false
        ),
        Err(Error::Advance(EffusionAdvanceError::WorkBudget {
            elapsed_s: 0.0,
            attempted_steps: 3,
            accepted_intervals: 0,
        }))
    ));
    assert!(matches!(
        advance_candidate(
            &donor,
            &receiver,
            bound,
            2.0,
            solve(),
            solve(),
            physical(),
            EffusionAdvanceControls {
                maximum_accepted_intervals: 1,
                ..adaptive(1e-4)
            },
            || false
        ),
        Err(Error::Advance(EffusionAdvanceError::WorkBudget {
            accepted_intervals: 1,
            ..
        }))
    ));
    let mut polls = 0;
    let full = advance_candidate(
        &donor,
        &receiver,
        bound,
        2.0,
        solve(),
        solve(),
        physical(),
        adaptive(1e-4),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(full.step().continuation.is_some());
    let mut seen = 0;
    let error = advance_candidate(
        &donor,
        &receiver,
        bound,
        2.0,
        solve(),
        solve(),
        physical(),
        adaptive(1e-4),
        || {
            seen += 1;
            seen >= polls
        },
    )
    .unwrap_err();
    assert_eq!(error, Error::Cancelled);
    assert_eq!(donor, saved_donor);
    assert_eq!(receiver, saved_receiver);
}

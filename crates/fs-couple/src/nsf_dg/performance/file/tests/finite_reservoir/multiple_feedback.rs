//! One finite actuated source, a cork duct and two jointly reacting escapes.
use super::*;
use crate::nsf_dg::performance::{NsfDgPerformance, ReservoirObservation};

const OBSERVATION: &str = "reservoir_observation aperture 1 baffled -0.08 0 0\n";
const SECOND_OBSERVATION: &str = "reservoir_observation aperture 2 baffled -0.08 0 0\n";
const FEEDBACKS: &str =
    "reservoir_feedback aperture 1 -20 20 16 16\nreservoir_feedback aperture 2 -20 20 16 16\n";
const SOLVE: &str =
    "reservoir_feedback_solve 12 1e-5 1e-8 0 4 8 1e-5 0.5 0.9 1.5 1e-4 0.5 1e-6\n";

fn multiple(text: String) -> String {
    text.replace("reservoir_apertures 2", "reservoir_apertures 3")
        .replace("reservoir_cell 5 2 1\n",
            "reservoir_cell 5 2 1\nreservoir_aperture 1 high reference 1\nreservoir_cell 9 2 1\n")
        .replace(OBSERVATION, &format!("{OBSERVATION}{SECOND_OBSERVATION}"))
        .replace("reservoir_feedback -20 20 16 16\n", &format!("{FEEDBACKS}{SOLVE}"))
        .replace("boundary finite-reservoir 0 1\nhold 64",
            "boundary finite-reservoir-actuation 0 0.05 1\nhold 64 -2e-7 0.02")
}

macro_rules! original_state {
    ($owner:expr) => {{
        let system = &$owner;
        (
            (
                system.zone().state_hash(),
                system.zone().cumulative_ledger().clone(),
                system.zone().last_reservoir_flows().to_vec(),
                system.zone().last_reservoir_exchanges().to_vec(),
                (0..3).map(|index| system.zone().reservoir_aperture_at(index).unwrap().state)
                    .collect::<Vec<_>>(),
            ),
            (
                system.load().state().clone(),
                system.load().terminal().unwrap().pole_states().to_vec(),
                system.load().last_independent_mouths_frame().cloned(),
                system.upstream_load().unwrap().state().clone(),
            ),
            (
                system.finite_reservoir().cloned(),
                system.finite_reservoir_inventory(),
                system.accepted_steps(),
            ),
        )
    }};
}

#[test]
fn actuated_finite_source_cork_and_two_escapes_replay_on_the_physical_clock() {
    let run = |block| {
        let file = build(true, block, multiple).unwrap();
        let rate = file.info().sample_rate_hz;
        let mut pressure = vec![0.0; file.info().samples as usize];
        let mut renderer = file.into_renderer();
        renderer.validate_sample_rate(rate).unwrap();
        for chunk in pressure.chunks_mut(block) {
            renderer.block(chunk).unwrap();
        }
        (pressure, renderer)
    };
    let (whole, a) = run(13);
    let (parts, b) = run(5);
    assert_eq!(whole, parts);
    assert!(whole.iter().any(|pressure| pressure.abs() > 1e-12));
    assert_eq!(original_state!(a.system()), original_state!(b.system()));
    assert_eq!(a.system().finite_reservoir_aperture(), Some(0));
    assert_eq!(a.system().reservoir_radiations().iter()
        .map(|(index, _)| *index).collect::<Vec<_>>(), [1, 2]);
    assert!(a.system().reservoir_radiation().is_none());
    let frame = a.last_frame().unwrap();
    assert!(frame.upstream.is_some());
    assert!(frame.reservoir_radiation.is_none());
    assert_eq!(frame.reservoir_radiations.iter()
        .map(|mouth| mouth.resolved_exchange.aperture).collect::<Vec<_>>(), [1, 2]);
    let source = frame.finite_reservoir.as_ref().unwrap();
    let motion = source.actuation.as_ref().unwrap();
    let dt = a.system().load().propagator().dt();
    assert_eq!(source.source_exchange.aperture, 0);
    assert_eq!(source.source_exchange, frame.reservoir_exchanges[0]);
    assert_eq!(motion.step.supplied_heat_j, 0.02 * dt);
    assert!(motion.step.work_by_gas_j < 0.0);
    assert!((source.final_volume_m3 - source.volume_m3 + 2e-7 * dt).abs()
        <= 4.0 * f64::EPSILON * source.volume_m3);
    assert_eq!(source.final_volume_m3, a.system().finite_reservoir().unwrap().volume_m3());
    assert!(source.energy_residual_j.abs() <= 1e-9);
    assert!(source.mass_residual_kg.abs() <= 1e-16);
}

#[test]
fn numerical_preview_budget_cancellation_and_late_source_failure_roll_back_all_owners() {
    let renderer = build(true, 1, multiple).unwrap().into_renderer();
    let drive = renderer.drive_history()[0];
    let mut system = renderer.system;
    let limits = system.limits;
    let gate = CancelGate::new_clock_free();
    let cold = original_state!(system);
    system.limits.maximum_evaluations = 3;
    let failure = system.step_with_drive(drive, &gate).unwrap_err();
    assert!(matches!(&failure, NsfDgError::Convergence { evaluations: 2, .. }),
        "the constructor and finite-difference base preview consume the budget: {failure:?}");
    assert_eq!(original_state!(system), cold);
    system.limits = limits;
    system.step_with_drive(drive, &gate).unwrap();
    let accepted = original_state!(system);

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(system.step_with_drive(drive, &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(original_state!(system), accepted);

    // This valid finite control fails the post-exchange chamber-volume gate,
    // after the nonlinear fluid and both duct candidates have been prepared.
    let impossible_compression = NsfDgDrive::FiniteReservoirActuation {
        aperture: 0,
        volume_rate_m3_s: -1.0,
        heat_power_w: 0.02,
        maximum_relative_volume_change: 0.05,
    };
    assert!(system.step_with_drive(impossible_compression, &gate).is_err());
    assert_eq!(original_state!(system), accepted);
    system.step_with_drive(drive, &gate).unwrap();

    let renderer = build(true, 1, multiple).unwrap().into_renderer();
    let mut uninterrupted = renderer.system;
    uninterrupted.step_with_drive(drive, &gate).unwrap();
    uninterrupted.step_with_drive(drive, &gate).unwrap();
    assert_eq!(original_state!(system), original_state!(uninterrupted),
        "failed previews must not alter the next accepted state or nonlinear warm start");
}

#[test]
fn each_passive_owner_keeps_its_original_placement_and_excludes_the_finite_source() {
    for replacement in [
        "reservoir_feedback aperture 1 -20 20 16 16\nreservoir_feedback aperture 1 -20 20 16 16\n",
        "reservoir_feedback aperture 0 -20 20 16 16\nreservoir_feedback aperture 2 -20 20 16 16\n",
    ] {
        assert!(build(true, 1, |text| multiple(text)
            .replace(FEEDBACKS, replacement)
            .replace(OBSERVATION,
                &format!("reservoir_observation aperture 0 baffled 0.06 0 0\n{OBSERVATION}")))
            .is_err());
    }
    for (origin_m, angle_rad, accepted) in [
        ([-0.08, 0.0], 0.0, true),
        ([-0.09, 0.0], 0.0, false),
        ([-0.08, 0.0], 0.1, false),
    ] {
        // Build a fresh observer from the same retained physical owner so a
        // duplicate-observer guard cannot mask the per-aperture placement gate.
        let renderer = build(true, 1, multiple).unwrap().into_renderer();
        let bare = NsfDgPerformance::new_resampled_with_drives(
            renderer.system,
            renderer.observation,
            renderer.sample_rate_hz,
            renderer.boundary,
            renderer.output_samples,
            1,
            renderer.gate,
        ).unwrap();
        assert_eq!(bare.with_reservoir_observation_at(
            2, ReservoirObservation { origin_m, angle_rad }
        ).is_ok(), accepted);
    }
}

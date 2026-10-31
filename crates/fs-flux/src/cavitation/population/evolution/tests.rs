use super::*;
use crate::cavitation::{
    BubbleGas, BubbleLimits, BubblePopulationClass, KellerMiksisBubble, KellerMiksisLimits,
};
use fs_material::cavitation::CavitationState;

const AMBIENT: f64 = 100_000.0;
const DURATION: f64 = 8e-6;

fn liquid() -> CavitationState {
    CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap()
}
fn bubble(radius: f64, maximum_radius_factor: f64) -> KellerMiksisBubble {
    let material = liquid();
    KellerMiksisBubble::try_new(
        material,
        BubbleGas::equilibrium(material, radius, AMBIENT, 1.4).unwrap(),
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: radius * 0.5,
                maximum_radius_m: radius * maximum_radius_factor,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap()
}
fn volume(radius: f64) -> f64 {
    4.0 * core::f64::consts::PI * radius.powi(3) / 3.0
}
fn class(radius: f64, fraction: f64, maximum_radius_factor: f64) -> BubblePopulationClass {
    BubblePopulationClass::try_new(
        bubble(radius, maximum_radius_factor),
        fraction / volume(radius),
    )
    .unwrap()
}
fn population() -> BubblePopulation {
    BubblePopulation::try_new(
        liquid(),
        vec![class(50e-6, 3e-5, 2.0), class(120e-6, 2e-5, 2.0)],
        2,
    )
    .unwrap()
}
fn limits() -> PopulationEvolutionLimits {
    PopulationEvolutionLimits {
        maximum_void_fraction: 0.001,
    }
}
fn drive(time: f64) -> BubblePressureSample {
    let omega = 5e4;
    BubblePressureSample {
        absolute_pressure_pa: AMBIENT - 3000.0 * (omega * time).sin(),
        time_derivative_pa_s: -3000.0 * omega * (omega * time).cos(),
    }
}
fn numerics() -> PopulationEvolutionNumerics {
    PopulationEvolutionNumerics {
        per_class: BubbleNumerics {
            relative_tolerance: 1e-10,
            radius_absolute_tolerance_m: 1e-16,
            velocity_absolute_tolerance_m_s: 1e-10,
            initial_step_s: DURATION / 100.0,
            minimum_step_s: DURATION * 1e-10,
            maximum_step_s: DURATION / 20.0,
            maximum_attempts: 10_000,
        },
        maximum_total_attempts: 100_000,
    }
}
fn owner() -> EvolvingBubblePopulation {
    EvolvingBubblePopulation::at_reference(population(), 0.0, limits(), drive(0.0)).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

// The independent oracle invokes each existing native owner directly, retaining
// its actual radius, gas inventory and numerical controls. It also reconstructs
// the population source from spherical geometry instead of summing new outputs.
#[test]
fn original_class_trajectories_and_actual_volume_derivatives_are_preserved() {
    let mut evolving = owner();
    let original = evolving.population().clone();
    let mut expected_states = evolving.states().to_vec();
    let controls = numerics();
    let mut accepted = 0;
    let mut rejected = 0;
    let mut rhs = 0;
    let mut expected_fraction = 0.0;
    let mut expected_rate = 0.0;
    let mut expected_acceleration = 0.0;
    let mut expected_number_radius = 0.0;
    for (class, state) in original.classes().iter().zip(&mut expected_states) {
        let report = class
            .bubble()
            .advance_controlled(state, DURATION, &controls.per_class, &drive, &mut || false)
            .unwrap();
        accepted += report.accepted_steps;
        rejected += report.rejected_steps;
        rhs += report.rhs_evaluations;
        let number = class.number_density_per_m3();
        let radius = state.radius_m();
        let velocity = state.wall_velocity_m_s();
        let acceleration = report.observation.wall_acceleration_m_s2;
        expected_fraction += number * volume(radius);
        expected_rate += number * 4.0 * core::f64::consts::PI * radius.powi(2) * velocity;
        expected_acceleration += number
            * 4.0
            * core::f64::consts::PI
            * (2.0 * radius * velocity.powi(2) + radius.powi(2) * acceleration);
        expected_number_radius += number * radius;
    }
    let step = evolving.advance(DURATION, &controls, &drive).unwrap();
    assert_eq!(evolving.states(), expected_states);
    assert_eq!(evolving.population(), &original);
    assert_eq!(evolving.time_s(), DURATION);
    assert!(
        evolving
            .states()
            .iter()
            .all(|state| state.time_s() == DURATION)
    );
    assert_eq!(step.advanced_classes, 2);
    assert_eq!(step.accepted_steps, accepted);
    assert_eq!(step.rejected_steps, rejected);
    assert_eq!(step.rhs_evaluations, rhs);
    assert_eq!(step.observation.imposed_pressure, drive(DURATION));
    assert_eq!(
        step.observation.number_density_per_m3,
        original.number_density_per_m3()
    );
    close(step.observation.void_fraction, expected_fraction, 2e-15);
    close(
        step.observation.void_fraction_rate_per_s,
        expected_rate,
        2e-15,
    );
    close(
        step.observation.void_fraction_acceleration_per_s2,
        expected_acceleration,
        3e-15,
    );
    close(
        step.observation.number_radius_moment_per_m2,
        expected_number_radius,
        2e-15,
    );
    assert_ne!(
        evolving.states()[0].radius_m() / original.classes()[0].bubble().gas().reference_radius_m(),
        evolving.states()[1].radius_m() / original.classes()[1].bubble().gas().reference_radius_m(),
    );
    close(evolving.maximum_void_fraction_envelope(), 8.0 * 5e-5, 2e-15);
}

#[test]
fn splitting_a_class_preserves_nonlinear_motion_and_geometric_source() {
    let original = class(50e-6, 3e-5, 2.0);
    let split =
        BubblePopulationClass::try_new(*original.bubble(), original.number_density_per_m3() / 8.0)
            .unwrap();
    let merged = BubblePopulation::try_new(liquid(), vec![original], 1).unwrap();
    let split = BubblePopulation::try_new(liquid(), vec![split; 8], 8).unwrap();
    let mut merged =
        EvolvingBubblePopulation::at_reference(merged, 0.0, limits(), drive(0.0)).unwrap();
    let mut split =
        EvolvingBubblePopulation::at_reference(split, 0.0, limits(), drive(0.0)).unwrap();
    let a = merged.advance(DURATION, &numerics(), &drive).unwrap();
    let b = split.advance(DURATION, &numerics(), &drive).unwrap();
    for state in split.states() {
        assert_eq!(*state, merged.states()[0]);
    }
    assert_eq!(
        a.observation.number_density_per_m3,
        b.observation.number_density_per_m3
    );
    assert_eq!(a.observation.void_fraction, b.observation.void_fraction);
    assert_eq!(
        a.observation.void_fraction_rate_per_s,
        b.observation.void_fraction_rate_per_s
    );
    assert_eq!(
        a.observation.void_fraction_acceleration_per_s2,
        b.observation.void_fraction_acceleration_per_s2
    );
    assert_eq!(
        a.observation.number_radius_moment_per_m2,
        b.observation.number_radius_moment_per_m2
    );
    assert_eq!(b.accepted_steps, 8 * a.accepted_steps);
    assert_eq!(b.rejected_steps, 8 * a.rejected_steps);
}

#[test]
fn dropped_complete_proposal_keeps_every_accepted_state_and_allows_identical_retry() {
    let mut evolving = owner();
    let original = evolving.population().clone();
    let initial = evolving.states().to_vec();
    let proposal = evolving
        .propose_advance_controlled(DURATION, &numerics(), &drive, &mut || false)
        .unwrap();
    assert_eq!(proposal.population(), &original);
    assert!(
        proposal
            .states()
            .iter()
            .all(|state| state.time_s() == DURATION)
    );
    let proposed_states = proposal.states().to_vec();
    let proposed_step = *proposal.step();
    assert_eq!(proposal.observation(), &proposed_step.observation);
    drop(proposal);
    assert_eq!(evolving.states(), initial);
    assert_eq!(evolving.time_s(), 0.0);
    assert_eq!(evolving.population(), &original);
    let committed = evolving
        .propose_advance_controlled(DURATION, &numerics(), &drive, &mut || false)
        .unwrap()
        .commit();
    assert_eq!(committed, proposed_step);
    assert_eq!(evolving.states(), proposed_states);
    assert_eq!(evolving.time_s(), DURATION);
}

#[test]
fn cancellation_before_during_and_after_all_classes_is_atomic_and_replayable() {
    let mut reference = owner();
    let mut polls = 0;
    let expected = reference
        .advance_controlled(DURATION, &numerics(), &drive, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 20);
    for stop in [1, polls / 2, polls] {
        let mut evolving = owner();
        let initial = evolving.states().to_vec();
        let original = evolving.population().clone();
        let mut count = 0;
        let result = evolving.advance_controlled(DURATION, &numerics(), &drive, &mut || {
            count += 1;
            count == stop
        });
        assert_eq!(result.unwrap_err(), PopulationEvolutionError::Cancelled);
        assert_eq!(evolving.states(), initial);
        assert_eq!(evolving.time_s(), 0.0);
        assert_eq!(evolving.population(), &original);
        let retry = evolving.advance(DURATION, &numerics(), &drive).unwrap();
        assert_eq!(retry, expected);
        assert_eq!(evolving.states(), reference.states());
    }
}

#[test]
fn shared_work_limit_counts_all_native_trials_and_rolls_back_completed_classes() {
    let mut evolving = owner();
    let initial = evolving.states().to_vec();
    let mut controls = numerics();
    let mut first = initial[0];
    let first_step = evolving.population().classes()[0]
        .bubble()
        .advance(&mut first, DURATION, &controls.per_class, &drive)
        .unwrap();
    let first_work = first_step.accepted_steps + first_step.rejected_steps;
    controls.maximum_total_attempts = first_work + 1;
    let pressure_calls = core::cell::Cell::new(0usize);
    let result = evolving.advance(DURATION, &controls, &|time| {
        pressure_calls.set(pressure_calls.get() + 1);
        drive(time)
    });
    assert_eq!(
        result.unwrap_err(),
        PopulationEvolutionError::WorkBudgetExceeded {
            class_index: 1,
            consumed_attempts: controls.maximum_total_attempts,
            maximum_total_attempts: controls.maximum_total_attempts,
        }
    );
    assert!(pressure_calls.get() <= 2 * 2 + 6 * controls.maximum_total_attempts + 1);
    assert_eq!(evolving.states(), initial);
    assert_eq!(evolving.time_s(), 0.0);
    assert!(first.time_s() > evolving.states()[0].time_s());
}

#[test]
fn a_late_class_domain_failure_keeps_the_previously_completed_class_unpublished() {
    let original = BubblePopulation::try_new(
        liquid(),
        vec![class(50e-6, 3e-5, 2.0), class(120e-6, 2e-5, 1.0 + 1e-8)],
        2,
    )
    .unwrap();
    let mut evolving =
        EvolvingBubblePopulation::at_reference(original.clone(), 0.0, limits(), drive(0.0))
            .unwrap();
    let initial = evolving.states().to_vec();
    let error = evolving.advance(DURATION, &numerics(), &drive).unwrap_err();
    assert!(
        matches!(
            error,
            PopulationEvolutionError::Class {
                class_index: 1,
                cause: CavitationError::RadiusOutsideDomain { .. },
            }
        ),
        "{error:?}"
    );
    assert_eq!(evolving.states(), initial);
    assert_eq!(evolving.time_s(), 0.0);
    assert_eq!(evolving.population(), &original);
}

#[test]
fn absent_populations_have_zero_source_and_can_advance_without_a_trial_budget() {
    let empty = BubblePopulation::try_new(liquid(), vec![], 0).unwrap();
    let ordinary = bubble(100e-6, 2.0);
    let unstable = KellerMiksisBubble::try_new(
        liquid(),
        BubbleGas::try_new(100e-6, 0.0, 1.0).unwrap(),
        ordinary.compressibility_limits(),
    )
    .unwrap();
    assert!(unstable.linearization().is_err());
    let absent = BubblePopulation::try_new(
        liquid(),
        vec![BubblePopulationClass::try_new(unstable, 0.0).unwrap()],
        1,
    )
    .unwrap();
    let mut controls = numerics();
    controls.maximum_total_attempts = 0;
    for population in [empty, absent] {
        let mut evolving =
            EvolvingBubblePopulation::at_reference(population, 0.0, limits(), drive(0.0)).unwrap();
        let step = evolving.advance(DURATION, &controls, &drive).unwrap();
        assert_eq!(evolving.time_s(), DURATION);
        assert!(
            evolving
                .states()
                .iter()
                .all(|state| state.time_s() == DURATION)
        );
        assert_eq!(step.advanced_classes, 0);
        assert_eq!(step.accepted_steps, 0);
        assert_eq!(step.rejected_steps, 0);
        assert_eq!(step.rhs_evaluations, 0);
        assert_eq!(step.observation.occupied_classes, 0);
        assert_eq!(step.observation.number_density_per_m3, 0.0);
        assert_eq!(step.observation.void_fraction, 0.0);
        assert_eq!(step.observation.void_fraction_rate_per_s, 0.0);
        assert_eq!(step.observation.void_fraction_acceleration_per_s2, 0.0);
        assert_eq!(step.observation.number_radius_moment_per_m2, 0.0);
    }
    let mut occupied = owner();
    let step = occupied.advance(0.0, &controls, &drive).unwrap();
    assert_eq!(step.accepted_steps, 0);
    assert_eq!(step.observation, occupied.observe(drive(0.0)).unwrap());
    assert!(matches!(
        occupied.advance(DURATION, &controls, &drive),
        Err(PopulationEvolutionError::WorkBudgetExceeded {
            class_index: 0,
            consumed_attempts: 0,
            maximum_total_attempts: 0,
        })
    ));
}

#[test]
fn entire_radius_envelope_and_exact_shared_clock_are_admission_conditions() {
    let broad = BubblePopulation::try_new(liquid(), vec![class(50e-6, 3e-5, 20.0)], 1).unwrap();
    assert!(broad.void_fraction() < limits().maximum_void_fraction);
    assert!(matches!(
        EvolvingBubblePopulation::at_reference(broad, 0.0, limits(), drive(0.0)),
        Err(PopulationEvolutionError::DilutionEnvelopeExceeded { .. })
    ));
    let states = vec![
        BubbleState::try_new(0.0, 50e-6, 0.0).unwrap(),
        BubbleState::try_new(1e-20, 120e-6, 0.0).unwrap(),
    ];
    assert!(matches!(
        EvolvingBubblePopulation::try_new(population(), 0.0, states, limits(), drive(0.0)),
        Err(PopulationEvolutionError::ClockMismatch { class_index: 1, .. })
    ));
    assert!(matches!(
        EvolvingBubblePopulation::try_new(population(), 0.0, vec![], limits(), drive(0.0)),
        Err(PopulationEvolutionError::StateCountMismatch {
            expected: 2,
            actual: 0
        })
    ));
    for maximum_void_fraction in [0.0, 0.010_000_001, f64::NAN] {
        assert!(
            EvolvingBubblePopulation::at_reference(
                population(),
                0.0,
                PopulationEvolutionLimits {
                    maximum_void_fraction
                },
                drive(0.0),
            )
            .is_err()
        );
    }
    let mut evolving = owner();
    let initial = evolving.states().to_vec();
    assert!(evolving.advance(-1.0, &numerics(), &drive).is_err());
    let mut invalid = numerics();
    invalid.maximum_total_attempts = usize::MAX;
    assert!(evolving.advance(DURATION, &invalid, &drive).is_err());
    assert_eq!(evolving.states(), initial);
}

#[test]
fn source_derivatives_match_independent_local_geometry_at_nonzero_wall_velocity() {
    let population = population();
    let time = 3e-4;
    let states = vec![
        BubbleState::try_new(time, 1.03 * 50e-6, 0.35).unwrap(),
        BubbleState::try_new(time, 0.98 * 120e-6, -0.12).unwrap(),
    ];
    let pressure = BubblePressureSample::try_new(AMBIENT, 0.0).unwrap();
    let evolving = EvolvingBubblePopulation::try_new(
        population.clone(),
        time,
        states.clone(),
        limits(),
        pressure,
    )
    .unwrap();
    let source = evolving.observe(pressure).unwrap();
    let h = 1e-9;
    let mut minus = 0.0;
    let mut center = 0.0;
    let mut plus = 0.0;
    for (class, state) in population.classes().iter().zip(states) {
        let acceleration = class
            .bubble()
            .observe(&state, pressure)
            .unwrap()
            .wall_acceleration_m_s2;
        let r = state.radius_m();
        let u = state.wall_velocity_m_s();
        let n = class.number_density_per_m3();
        minus += n * volume(r - u * h + 0.5 * acceleration * h * h);
        center += n * volume(r);
        plus += n * volume(r + u * h + 0.5 * acceleration * h * h);
    }
    close(source.void_fraction, center, 2e-15);
    close(
        source.void_fraction_rate_per_s,
        (plus - minus) / (2.0 * h),
        2e-7,
    );
    close(
        source.void_fraction_acceleration_per_s2,
        (plus - 2.0 * center + minus) / (h * h),
        5e-5,
    );
}

#[test]
fn occupied_class_underflow_refuses_without_rejecting_zero_or_cancelled_derivatives() {
    // Idealized inviscid, surface-free material makes the reference state
    // exactly balanced, isolating weighted-source representability.
    let material = CavitationState::try_new(300.0, 1000.0, 0.0, 0.0, 0.0, 1500.0).unwrap();
    let model = KellerMiksisBubble::try_new(
        material,
        BubbleGas::equilibrium(material, 1.0, AMBIENT, 1.4).unwrap(),
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.25,
                maximum_radius_m: 2.0,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap();
    let population = |number| {
        BubblePopulation::try_new(
            material,
            vec![BubblePopulationClass::try_new(model, number).unwrap()],
            1,
        )
        .unwrap()
    };
    let static_pressure = BubblePressureSample::try_new(AMBIENT, 0.0).unwrap();
    let minimum_number = f64::from_bits(1);
    let reference = EvolvingBubblePopulation::at_reference(
        population(minimum_number),
        0.0,
        limits(),
        static_pressure,
    )
    .unwrap();
    let zero = reference.observe(static_pressure).unwrap();
    assert!(zero.void_fraction > 0.0);
    assert!(zero.number_radius_moment_per_m2 > 0.0);
    assert_eq!(zero.void_fraction_rate_per_s, 0.0);
    assert_eq!(zero.void_fraction_acceleration_per_s2, 0.0);

    let cases = [
        // Actual occupied volume vanishes even though its reference is valid.
        (minimum_number, 0.4, 0.0, static_pressure, 0),
        // Volume still rounds positive, but the n*R moment vanishes.
        (minimum_number, 0.5, 0.0, static_pressure, 1),
        (1e-300, 1.0, 1e-100, static_pressure, 2),
        (
            1e-300,
            1.0,
            0.0,
            BubblePressureSample::try_new(AMBIENT, 1e-200).unwrap(),
            3,
        ),
    ];
    for (number, radius, velocity, pressure, component) in cases {
        let state = BubbleState::try_new(0.0, radius, velocity).unwrap();
        let native = model.observe(&state, pressure).unwrap();
        match component {
            0 => assert_eq!(number * native.volume_m3, 0.0),
            1 => {
                assert!(number * native.volume_m3 > 0.0);
                assert_eq!(number * radius, 0.0);
            }
            2 => {
                assert_ne!(native.volume_rate_m3_s, 0.0);
                assert_eq!(number * native.volume_rate_m3_s, 0.0);
            }
            3 => {
                assert_eq!(native.volume_rate_m3_s, 0.0);
                assert_ne!(native.volume_acceleration_m3_s2, 0.0);
                assert_eq!(number * native.volume_acceleration_m3_s2, 0.0);
            }
            _ => unreachable!(),
        }
        assert!(matches!(
            EvolvingBubblePopulation::try_new(
                population(number),
                0.0,
                vec![state],
                limits(),
                pressure,
            ),
            Err(PopulationEvolutionError::InvalidInput(
                "unrepresentable weighted bubble source",
            )),
        ));
    }

    // Cancellation across two individually representable rates is physical.
    let class = BubblePopulationClass::try_new(model, 1e-20).unwrap();
    let balanced = EvolvingBubblePopulation::try_new(
        BubblePopulation::try_new(material, vec![class, class], 2).unwrap(),
        0.0,
        vec![
            BubbleState::try_new(0.0, 1.0, 1e-4).unwrap(),
            BubbleState::try_new(0.0, 1.0, -1e-4).unwrap(),
        ],
        limits(),
        static_pressure,
    )
    .unwrap();
    assert_eq!(
        balanced
            .observe(static_pressure)
            .unwrap()
            .void_fraction_rate_per_s,
        0.0,
    );
}

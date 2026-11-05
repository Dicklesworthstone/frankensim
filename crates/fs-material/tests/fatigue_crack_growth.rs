//! G0/G1/G3/G4: bounded Paris crack evolution and transactional load blocks.
//! Analytical, synthetic fixtures verify implementation; no measured crack
//! growth, material transferability or component service life is claimed.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGrowthError, CrackGrowthLimits, CrackGrowthState, ParisCrackGrowth,
    ParisDomain, ParisLaw, ParisParameters, StressCycle,
};

const INITIAL_LENGTH_M: f64 = 0.004;
const REFERENCE_RATE_M_PER_CYCLE: f64 = 2.0e-7;
const GEOMETRY_FACTOR: f64 = 1.12;
const STRESS_RANGE_PA: f64 = 100.0e6;

fn scope() -> FatigueScope {
    FatigueScope {
        source_reference: "synthetic analytical Paris fixture; not measured".into(),
        material: "synthetic crack-growth coupon".into(),
        process: "declared coupon process".into(),
        surface_and_geometry: "constant-Y mode-I crack".into(),
        environment: "fixed synthetic environment".into(),
        load_basis: "tensile normal stress range in Pa; full cycles".into(),
        temperature_k: 293.15,
    }
}

fn reference_delta_k() -> f64 {
    GEOMETRY_FACTOR * STRESS_RANGE_PA * (std::f64::consts::PI * INITIAL_LENGTH_M).sqrt()
}

fn parameters(exponent: f64) -> ParisParameters {
    ParisParameters {
        reference_delta_k_pa_sqrt_m: reference_delta_k(),
        reference_rate_m_per_cycle: REFERENCE_RATE_M_PER_CYCLE,
        exponent,
    }
}

fn domain() -> ParisDomain {
    ParisDomain {
        minimum_delta_k_pa_sqrt_m: 0.25 * reference_delta_k(),
        maximum_delta_k_pa_sqrt_m: 8.0 * reference_delta_k(),
        minimum_load_ratio: 0.0,
        maximum_load_ratio: 0.8,
        maximum_k_pa_sqrt_m: 10.0 * reference_delta_k(),
    }
}

fn law(exponent: f64) -> Arc<ParisLaw> {
    Arc::new(ParisLaw::new(scope(), parameters(exponent), domain()).unwrap())
}

fn geometry() -> ConstantGeometryCrack {
    ConstantGeometryCrack {
        geometry_factor: GEOMETRY_FACTOR,
        minimum_length_m: 0.001,
        maximum_length_m: 0.1,
    }
}

fn limits() -> CrackGrowthLimits {
    CrackGrowthLimits {
        maximum_blocks: 100,
        maximum_cycles: 100_000.0,
    }
}

fn run(exponent: f64) -> ParisCrackGrowth {
    ParisCrackGrowth::new(law(exponent), geometry(), INITIAL_LENGTH_M, limits()).unwrap()
}

fn cycle(cycles: f64) -> StressCycle {
    StressCycle {
        minimum_pa: 50.0e6,
        maximum_pa: 150.0e6,
        cycles,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3.0e-12 * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:e} != {expected:e}"
    );
}

fn state_bits(state: CrackGrowthState) -> (u64, u64, u64) {
    (
        state.crack_length_m.to_bits(),
        state.cycles.to_bits(),
        state.blocks,
    )
}

#[test]
fn g1_constant_geometry_matches_three_independent_closed_forms() {
    let cycles = 2_000.0;
    let growth_fraction = REFERENCE_RATE_M_PER_CYCLE * cycles / INITIAL_LENGTH_M;
    // At a0, the declared stress range gives exactly the reference Delta K.
    // Integrating da/dN = r0 * (a/a0)^(m/2) gives these independent forms.
    let cases = [
        (
            1.0,
            INITIAL_LENGTH_M * (1.0 + growth_fraction / 2.0).powi(2),
        ),
        (2.0, INITIAL_LENGTH_M * growth_fraction.exp()),
        (4.0, INITIAL_LENGTH_M / (1.0 - growth_fraction)),
    ];
    for (exponent, expected_length) in cases {
        let mut evolution = run(exponent);
        let before = evolution.state();
        let step = evolution.advance(cycle(cycles)).unwrap();
        assert_eq!(state_bits(step.before), state_bits(before));
        assert_eq!(state_bits(step.after), state_bits(evolution.state()));
        close(step.after.crack_length_m, expected_length);
        assert_eq!(step.after.cycles, cycles);
        assert_eq!(step.after.blocks, 1);
        close(step.initial_delta_k_pa_sqrt_m, reference_delta_k());
        close(
            step.final_delta_k_pa_sqrt_m,
            reference_delta_k() * (expected_length / INITIAL_LENGTH_M).sqrt(),
        );
        close(step.initial_maximum_k_pa_sqrt_m, 1.5 * reference_delta_k());
        close(
            step.final_maximum_k_pa_sqrt_m,
            1.5 * reference_delta_k() * (expected_length / INITIAL_LENGTH_M).sqrt(),
        );
    }
}

#[test]
fn g1_near_exponent_two_retains_the_continuous_limit_and_its_first_correction() {
    let cycles = 10_000.0;
    let t = REFERENCE_RATE_M_PER_CYCLE * cycles / INITIAL_LENGTH_M;
    for exponent in [2.0 - 1.0e-8, 2.0, 2.0 + 1.0e-8] {
        let q = 1.0 - 0.5 * exponent;
        // Independently expand the integral near m=2. Here |q*t| < 3e-9;
        // the omitted term in log(a/a0) is below 3e-27. Simply treating all
        // nearby exponents as two would lose the measurable first correction.
        let log_ratio = t - 0.5 * q * t * t + q * q * t * t * t / 3.0;
        let expected = INITIAL_LENGTH_M * log_ratio.exp();
        let mut evolution = run(exponent);
        evolution.advance(cycle(cycles)).unwrap();
        close(evolution.state().crack_length_m, expected);
    }
}

#[test]
fn g1_paris_uses_stress_intensity_range_and_full_cycles() {
    let model = law(2.0);
    close(
        model.growth_rate(reference_delta_k(), 1.0 / 3.0).unwrap(),
        REFERENCE_RATE_M_PER_CYCLE,
    );
    close(
        model
            .growth_rate(0.5 * reference_delta_k(), 1.0 / 3.0)
            .unwrap(),
        0.25 * REFERENCE_RATE_M_PER_CYCLE,
    );
    let cycles = 1_000.0;
    let mut original = run(2.0);
    let mut half_range = run(2.0);
    original.advance(cycle(cycles)).unwrap();
    half_range
        .advance(StressCycle {
            minimum_pa: 25.0e6,
            maximum_pa: 75.0e6,
            cycles,
        })
        .unwrap();
    let growth_fraction = REFERENCE_RATE_M_PER_CYCLE * cycles / INITIAL_LENGTH_M;
    close(
        original.state().crack_length_m,
        INITIAL_LENGTH_M * growth_fraction.exp(),
    );
    close(
        half_range.state().crack_length_m,
        INITIAL_LENGTH_M * (0.25 * growth_fraction).exp(),
    );
    assert!(original.state().crack_length_m > half_range.state().crack_length_m);
}

#[test]
fn g3_variable_load_blocks_and_callback_partitions_preserve_the_same_evolution() {
    let blocks = [
        cycle(250.0),
        StressCycle {
            minimum_pa: 25.0e6,
            maximum_pa: 75.0e6,
            cycles: 700.0,
        },
        StressCycle {
            minimum_pa: 75.0e6,
            maximum_pa: 225.0e6,
            cycles: 450.0,
        },
    ];
    let mut whole = run(2.0);
    let result = whole.advance_block(&blocks, &mut || false).unwrap();
    assert_eq!(state_bits(result), state_bits(whole.state()));
    let equivalent_cycles = 250.0 + 0.25 * 700.0 + 2.25 * 450.0;
    close(
        result.crack_length_m,
        INITIAL_LENGTH_M
            * (REFERENCE_RATE_M_PER_CYCLE * equivalent_cycles / INITIAL_LENGTH_M).exp(),
    );
    assert_eq!(result.cycles, 1_400.0);
    assert_eq!(result.blocks, 3);

    let mut sequential = run(2.0);
    for block in blocks {
        sequential.advance(block).unwrap();
    }
    assert_eq!(state_bits(result), state_bits(sequential.state()));
    for split in 1..blocks.len() {
        let mut partitioned = run(2.0);
        partitioned
            .advance_block(&blocks[..split], &mut || false)
            .unwrap();
        partitioned
            .advance_block(&blocks[split..], &mut || false)
            .unwrap();
        assert_eq!(state_bits(result), state_bits(partitioned.state()));
    }

    // Physically subdividing a constant load adds block boundaries but has the
    // same analytical crack evolution, within floating-point roundoff.
    let mut single = run(4.0);
    single.advance(cycle(2_000.0)).unwrap();
    let mut subdivided = run(4.0);
    for count in [200.0, 500.0, 1_300.0] {
        subdivided.advance(cycle(count)).unwrap();
    }
    close(
        subdivided.state().crack_length_m,
        single.state().crack_length_m,
    );
    assert_eq!(subdivided.state().cycles, single.state().cycles);
}

#[test]
fn g0_zero_range_does_not_grow_but_static_intensity_still_limits_admission() {
    let mut evolution = run(2.0);
    for stress in [0.0, 75.0e6] {
        let step = evolution
            .advance(StressCycle {
                minimum_pa: stress,
                maximum_pa: stress,
                cycles: 250.0,
            })
            .unwrap();
        assert_eq!(step.after.crack_length_m, INITIAL_LENGTH_M);
        assert_eq!(step.initial_delta_k_pa_sqrt_m, 0.0);
        assert_eq!(step.final_delta_k_pa_sqrt_m, 0.0);
    }
    assert_eq!(evolution.state().cycles, 500.0);
    assert_eq!(evolution.state().blocks, 2);
    let before = evolution.state();
    assert!(
        evolution
            .advance(StressCycle {
                minimum_pa: 1.0e12,
                maximum_pa: 1.0e12,
                cycles: 1.0,
            })
            .is_err()
    );
    assert_eq!(state_bits(evolution.state()), state_bits(before));
}

#[test]
fn g0_missing_support_and_compression_refuse_without_inventing_an_endurance_limit() {
    let mut evolution = run(2.0);
    evolution.advance(cycle(100.0)).unwrap();
    let before = evolution.state();
    let unsupported = [
        StressCycle {
            minimum_pa: 5.0e6,
            maximum_pa: 15.0e6,
            cycles: 100.0,
        },
        StressCycle {
            minimum_pa: -10.0e6,
            maximum_pa: 150.0e6,
            cycles: 100.0,
        },
        StressCycle {
            minimum_pa: 149.0e6,
            maximum_pa: 150.0e6,
            cycles: 100.0,
        },
        StressCycle {
            minimum_pa: 150.0e6,
            maximum_pa: 50.0e6,
            cycles: 100.0,
        },
        StressCycle {
            minimum_pa: 50.0e6,
            maximum_pa: 150.0e6,
            cycles: f64::NAN,
        },
    ];
    for block in unsupported {
        assert!(evolution.advance(block).is_err());
        assert_eq!(state_bits(evolution.state()), state_bits(before));
    }
    assert!(
        evolution
            .law()
            .growth_rate(0.1 * reference_delta_k(), 1.0 / 3.0)
            .is_err()
    );
    assert!(
        evolution
            .law()
            .growth_rate(reference_delta_k(), -0.1)
            .is_err()
    );
}

#[test]
fn g0_growth_must_remain_inside_geometry_and_intensity_domains_without_clipping() {
    let mut bounded_geometry = geometry();
    bounded_geometry.maximum_length_m = 1.05 * INITIAL_LENGTH_M;
    let mut geometry_run =
        ParisCrackGrowth::new(law(2.0), bounded_geometry, INITIAL_LENGTH_M, limits()).unwrap();
    let before = geometry_run.state();
    assert!(geometry_run.advance(cycle(2_000.0)).is_err());
    assert_eq!(state_bits(geometry_run.state()), state_bits(before));

    let mut range_domain = domain();
    range_domain.maximum_delta_k_pa_sqrt_m = 1.04 * reference_delta_k();
    let range_law = Arc::new(ParisLaw::new(scope(), parameters(2.0), range_domain).unwrap());
    let mut range_run =
        ParisCrackGrowth::new(range_law, geometry(), INITIAL_LENGTH_M, limits()).unwrap();
    let before = range_run.state();
    assert!(range_run.advance(cycle(2_000.0)).is_err());
    assert_eq!(state_bits(range_run.state()), state_bits(before));

    let mut static_domain = domain();
    static_domain.maximum_k_pa_sqrt_m = 1.55 * reference_delta_k();
    let static_law = Arc::new(ParisLaw::new(scope(), parameters(2.0), static_domain).unwrap());
    let mut static_run =
        ParisCrackGrowth::new(static_law, geometry(), INITIAL_LENGTH_M, limits()).unwrap();
    let before = static_run.state();
    assert!(static_run.advance(cycle(2_000.0)).is_err());
    assert_eq!(state_bits(static_run.state()), state_bits(before));
}

#[test]
fn g0_calibration_scope_and_evolution_budgets_are_explicit() {
    let model = law(2.0);
    assert_eq!(model.scope(), &scope());
    assert_eq!(model.parameters().exponent, 2.0);
    assert_eq!(model.domain().maximum_load_ratio, 0.8);
    let evolution = ParisCrackGrowth::new(model, geometry(), INITIAL_LENGTH_M, limits()).unwrap();
    assert_eq!(evolution.law().scope(), &scope());
    assert_eq!(evolution.geometry().geometry_factor, GEOMETRY_FACTOR);

    let mut missing_source = scope();
    missing_source.source_reference.clear();
    assert!(ParisLaw::new(missing_source, parameters(2.0), domain()).is_err());
    let mut invalid_rate = parameters(2.0);
    invalid_rate.reference_rate_m_per_cycle = 0.0;
    assert!(ParisLaw::new(scope(), invalid_rate, domain()).is_err());
    assert!(ParisLaw::new(scope(), parameters(0.0), domain()).is_err());
    let mut reversed_domain = domain();
    reversed_domain.maximum_delta_k_pa_sqrt_m = 0.5 * reversed_domain.minimum_delta_k_pa_sqrt_m;
    assert!(ParisLaw::new(scope(), parameters(2.0), reversed_domain).is_err());
    assert!(ParisCrackGrowth::new(law(2.0), geometry(), 0.0, limits()).is_err());

    let mut capped = ParisCrackGrowth::new(
        law(2.0),
        geometry(),
        INITIAL_LENGTH_M,
        CrackGrowthLimits {
            maximum_blocks: 1,
            maximum_cycles: 500.0,
        },
    )
    .unwrap();
    let before = capped.state();
    assert!(capped.advance(cycle(501.0)).is_err());
    assert_eq!(state_bits(capped.state()), state_bits(before));
    capped.advance(cycle(100.0)).unwrap();
    let before = capped.state();
    assert!(capped.advance(cycle(100.0)).is_err());
    assert_eq!(state_bits(capped.state()), state_bits(before));
}

#[test]
fn g4_late_failure_and_cancellation_publish_no_partial_crack_growth() {
    let mut evolution = run(2.0);
    evolution.advance(cycle(100.0)).unwrap();
    let before = evolution.state();
    assert!(
        evolution
            .advance_block(
                &[
                    cycle(100.0),
                    StressCycle {
                        minimum_pa: -1.0,
                        ..cycle(100.0)
                    }
                ],
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(state_bits(evolution.state()), state_bits(before));

    let mut polls = 0;
    let mut reference = run(2.0);
    reference.advance(cycle(100.0)).unwrap();
    reference
        .advance_block(&[cycle(100.0), cycle(200.0)], &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(
        polls > 1,
        "cancellation must cover more than initial admission"
    );
    for stop_at in [1, polls] {
        let mut attempt = run(2.0);
        attempt.advance(cycle(100.0)).unwrap();
        let accepted = attempt.state();
        let mut calls = 0;
        let result = attempt.advance_block(&[cycle(100.0), cycle(200.0)], &mut || {
            calls += 1;
            calls == stop_at
        });
        assert!(matches!(result, Err(CrackGrowthError::Cancelled)));
        assert_eq!(state_bits(attempt.state()), state_bits(accepted));
    }
    let mut cancelled = || true;
    assert!(matches!(
        evolution.advance_with_cancel(cycle(100.0), &mut cancelled),
        Err(CrackGrowthError::Cancelled)
    ));
    assert_eq!(state_bits(evolution.state()), state_bits(before));
}

#[test]
fn g1_forecast_inverts_three_independent_integrals_from_the_accepted_crack() {
    let target_length = 1.5 * INITIAL_LENGTH_M;
    for exponent in [1.0, 2.0, 4.0] {
        let mut evolution = run(exponent);
        evolution.advance(cycle(250.0)).unwrap();
        let before = evolution.state();
        let current = before.crack_length_m;
        // Independently integrate dN = da / (r_ref (a/a_ref)^(m/2)).
        let expected_cycles = if exponent == 1.0 {
            2.0 * INITIAL_LENGTH_M.sqrt() / REFERENCE_RATE_M_PER_CYCLE * (target_length - current)
                / (target_length.sqrt() + current.sqrt())
        } else if exponent == 2.0 {
            INITIAL_LENGTH_M / REFERENCE_RATE_M_PER_CYCLE * (target_length / current).ln()
        } else {
            INITIAL_LENGTH_M.powi(2) / REFERENCE_RATE_M_PER_CYCLE * (target_length - current)
                / (current * target_length)
        };
        let forecast = evolution
            .cycles_to_length(50.0e6, 150.0e6, target_length)
            .unwrap();
        assert_eq!(state_bits(evolution.state()), state_bits(before));
        assert_eq!(state_bits(forecast.before), state_bits(before));
        assert_eq!(forecast.target_length_m, target_length);
        close(forecast.additional_cycles, expected_cycles);
        close(
            forecast.initial_delta_k_pa_sqrt_m,
            reference_delta_k() * (current / INITIAL_LENGTH_M).sqrt(),
        );
        close(
            forecast.target_delta_k_pa_sqrt_m,
            reference_delta_k() * (target_length / INITIAL_LENGTH_M).sqrt(),
        );
        close(
            forecast.initial_maximum_k_pa_sqrt_m,
            1.5 * reference_delta_k() * (current / INITIAL_LENGTH_M).sqrt(),
        );
        close(
            forecast.target_maximum_k_pa_sqrt_m,
            1.5 * reference_delta_k() * (target_length / INITIAL_LENGTH_M).sqrt(),
        );
        let reached = evolution
            .advance(cycle(forecast.additional_cycles))
            .unwrap();
        close(reached.after.crack_length_m, target_length);
        assert_eq!(
            reached.after.cycles,
            before.cycles + forecast.additional_cycles
        );
        assert_eq!(reached.after.blocks, before.blocks + 1);
    }
}

#[test]
fn g1_forecast_preserves_the_near_two_correction_and_adjacent_float_target() {
    let target_length = 1.5 * INITIAL_LENGTH_M;
    let log_ratio = (target_length / INITIAL_LENGTH_M).ln();
    for exponent in [2.0 - 1.0e-8, 2.0, 2.0 + 1.0e-8] {
        let mut evolution = run(exponent);
        let before = evolution.state();
        let q = 1.0 - 0.5 * exponent;
        // Independent Taylor integral about m=2. The omitted term in the
        // bracket is below 2e-28, while the first correction is measurable.
        let expected = INITIAL_LENGTH_M / REFERENCE_RATE_M_PER_CYCLE
            * (log_ratio + 0.5 * q * log_ratio.powi(2) + q * q * log_ratio.powi(3) / 6.0);
        let forecast = evolution
            .cycles_to_length(50.0e6, 150.0e6, target_length)
            .unwrap();
        close(forecast.additional_cycles, expected);
        assert_eq!(state_bits(evolution.state()), state_bits(before));
        evolution
            .advance(cycle(forecast.additional_cycles))
            .unwrap();
        close(evolution.state().crack_length_m, target_length);
    }

    let mut evolution = run(2.0);
    let before = evolution.state();
    let adjacent_target = f64::from_bits(INITIAL_LENGTH_M.to_bits() + 1);
    let increment = adjacent_target - INITIAL_LENGTH_M;
    let relative_increment = increment / INITIAL_LENGTH_M;
    // For one representable length increment, forming target/current first
    // loses several percent of the requested cycle count. The independent
    // log(1+x) series keeps the actual, exactly subtracted length increment.
    let expected = increment / REFERENCE_RATE_M_PER_CYCLE
        * (1.0 - 0.5 * relative_increment + relative_increment.powi(2) / 3.0);
    let forecast = evolution
        .cycles_to_length(50.0e6, 150.0e6, adjacent_target)
        .unwrap();
    assert!(forecast.additional_cycles > 0.0);
    close(forecast.additional_cycles, expected);
    assert_eq!(state_bits(evolution.state()), state_bits(before));
    evolution
        .advance(cycle(forecast.additional_cycles))
        .unwrap();
    assert_eq!(
        evolution.state().crack_length_m.to_bits(),
        adjacent_target.to_bits()
    );
}

#[test]
fn g0_forecast_refuses_invalid_targets_and_the_entire_unsupported_growth_interval() {
    let evolution = run(2.0);
    let before = evolution.state();
    for target in [
        INITIAL_LENGTH_M,
        0.5 * INITIAL_LENGTH_M,
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
    ] {
        assert!(evolution.cycles_to_length(50.0e6, 150.0e6, target).is_err());
        assert_eq!(state_bits(evolution.state()), state_bits(before));
    }
    for (minimum, maximum) in [
        (0.0, 0.0),
        (75.0e6, 75.0e6),
        (-1.0, 150.0e6),
        (150.0e6, 50.0e6),
        (f64::NAN, 150.0e6),
        (50.0e6, f64::INFINITY),
        (5.0e6, 15.0e6),
        (149.0e6, 150.0e6),
        (50.0e6, 1.0e12),
    ] {
        assert!(
            evolution
                .cycles_to_length(minimum, maximum, 1.1 * INITIAL_LENGTH_M)
                .is_err()
        );
        assert_eq!(state_bits(evolution.state()), state_bits(before));
    }

    let mut bounded_geometry = geometry();
    bounded_geometry.maximum_length_m = 1.05 * INITIAL_LENGTH_M;
    let geometry_run =
        ParisCrackGrowth::new(law(2.0), bounded_geometry, INITIAL_LENGTH_M, limits()).unwrap();
    assert!(matches!(
        geometry_run.cycles_to_length(50.0e6, 150.0e6, 1.1 * INITIAL_LENGTH_M),
        Err(CrackGrowthError::CrackLengthOutsideDomain { .. })
    ));
    assert_eq!(state_bits(geometry_run.state()), state_bits(before));

    let mut range_domain = domain();
    range_domain.maximum_delta_k_pa_sqrt_m = 1.04 * reference_delta_k();
    let range_law = Arc::new(ParisLaw::new(scope(), parameters(2.0), range_domain).unwrap());
    let range_run =
        ParisCrackGrowth::new(range_law, geometry(), INITIAL_LENGTH_M, limits()).unwrap();
    assert!(matches!(
        range_run.cycles_to_length(50.0e6, 150.0e6, 1.1 * INITIAL_LENGTH_M),
        Err(CrackGrowthError::DeltaKOutsideDomain { .. })
    ));
    assert_eq!(state_bits(range_run.state()), state_bits(before));

    let mut static_domain = domain();
    static_domain.maximum_k_pa_sqrt_m = 1.55 * reference_delta_k();
    let static_law = Arc::new(ParisLaw::new(scope(), parameters(2.0), static_domain).unwrap());
    let static_run =
        ParisCrackGrowth::new(static_law, geometry(), INITIAL_LENGTH_M, limits()).unwrap();
    assert!(matches!(
        static_run.cycles_to_length(50.0e6, 150.0e6, 1.1 * INITIAL_LENGTH_M),
        Err(CrackGrowthError::MaximumKExceeded { .. })
    ));
    assert_eq!(state_bits(static_run.state()), state_bits(before));
}

#[test]
fn g4_forecast_respects_remaining_cycle_block_and_float_resolution_budgets() {
    let mut exposure_limited = ParisCrackGrowth::new(
        law(2.0),
        geometry(),
        INITIAL_LENGTH_M,
        CrackGrowthLimits {
            maximum_blocks: 10,
            maximum_cycles: 1_000.0,
        },
    )
    .unwrap();
    exposure_limited.advance(cycle(600.0)).unwrap();
    let before = exposure_limited.state();
    let fractional_rate = REFERENCE_RATE_M_PER_CYCLE / INITIAL_LENGTH_M;
    let target_after_500 = before.crack_length_m * (fractional_rate * 500.0).exp();
    assert!(matches!(
        exposure_limited.cycles_to_length(50.0e6, 150.0e6, target_after_500),
        Err(CrackGrowthError::Budget)
    ));
    assert_eq!(state_bits(exposure_limited.state()), state_bits(before));
    let target_after_200 = before.crack_length_m * (fractional_rate * 200.0).exp();
    let admitted = exposure_limited
        .cycles_to_length(50.0e6, 150.0e6, target_after_200)
        .unwrap();
    close(admitted.additional_cycles, 200.0);
    assert_eq!(state_bits(exposure_limited.state()), state_bits(before));

    let mut block_limited = ParisCrackGrowth::new(
        law(2.0),
        geometry(),
        INITIAL_LENGTH_M,
        CrackGrowthLimits {
            maximum_blocks: 1,
            maximum_cycles: 100_000.0,
        },
    )
    .unwrap();
    block_limited.advance(cycle(100.0)).unwrap();
    let before = block_limited.state();
    assert!(matches!(
        block_limited.cycles_to_length(50.0e6, 150.0e6, 1.1 * before.crack_length_m),
        Err(CrackGrowthError::Budget)
    ));
    assert_eq!(state_bits(block_limited.state()), state_bits(before));

    let mut resolution_limited = ParisCrackGrowth::new(
        law(2.0),
        geometry(),
        INITIAL_LENGTH_M,
        CrackGrowthLimits {
            maximum_blocks: 10,
            maximum_cycles: (1_u64 << 54) as f64,
        },
    )
    .unwrap();
    resolution_limited
        .advance(StressCycle {
            minimum_pa: 0.0,
            maximum_pa: 0.0,
            cycles: (1_u64 << 53) as f64,
        })
        .unwrap();
    let before = resolution_limited.state();
    let adjacent_target = f64::from_bits(INITIAL_LENGTH_M.to_bits() + 1);
    assert!(matches!(
        resolution_limited.cycles_to_length(50.0e6, 150.0e6, adjacent_target),
        Err(CrackGrowthError::Budget)
    ));
    assert_eq!(state_bits(resolution_limited.state()), state_bits(before));
}

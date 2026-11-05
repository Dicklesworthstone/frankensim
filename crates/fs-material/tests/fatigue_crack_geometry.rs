//! G0/G1/G3/G4: evolving crack geometry against independent closed forms.
//! These coefficients and geometry tables are synthetic implementation fixtures,
//! not measured crack-growth data, component calibrations, or service-life claims.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGeometry, CrackGeometryPoint, CrackGeometryTableContext,
    CrackGrowthError, CrackGrowthLimits, CrackGrowthState, ParisCrackGrowth, ParisDomain, ParisLaw,
    ParisParameters, StressCycle, TabulatedCrackGeometry,
};

const A0: f64 = 1.0e-3;
const RATE: f64 = 1.0e-6;
const MINIMUM_PA: f64 = 25.0e6;
const MAXIMUM_PA: f64 = 75.0e6;

fn reference_k() -> f64 {
    (MAXIMUM_PA - MINIMUM_PA) * (std::f64::consts::PI * A0).sqrt()
}

fn context() -> CrackGeometryTableContext {
    CrackGeometryTableContext {
        source_reference: "synthetic analytical geometry; not measured".into(),
        crack_length_convention: "a is the declared half crack length in metres".into(),
        nominal_stress_basis: "far-field normal engineering stress in Pa".into(),
    }
}

fn domain() -> ParisDomain {
    ParisDomain {
        minimum_delta_k_pa_sqrt_m: 0.1 * reference_k(),
        maximum_delta_k_pa_sqrt_m: 10.0 * reference_k(),
        minimum_load_ratio: 0.0,
        maximum_load_ratio: 0.8,
        maximum_k_pa_sqrt_m: 15.0 * reference_k(),
    }
}

fn law(exponent: f64, fit: ParisDomain) -> Arc<ParisLaw> {
    Arc::new(
        ParisLaw::new(
            FatigueScope {
                source_reference: "synthetic normalized Paris law; not measured".into(),
                material: "synthetic coupon".into(),
                process: "synthetic fixed process".into(),
                surface_and_geometry: "declared mode-I nominal-stress geometry".into(),
                environment: "synthetic fixed environment".into(),
                load_basis: "full tensile stress-intensity range; full cycles".into(),
                temperature_k: 300.0,
            },
            ParisParameters {
                reference_delta_k_pa_sqrt_m: reference_k(),
                reference_rate_m_per_cycle: RATE,
                exponent,
            },
            fit,
        )
        .unwrap(),
    )
}

fn points(values: &[(f64, f64)]) -> Vec<CrackGeometryPoint> {
    values
        .iter()
        .map(|&(crack_length_m, geometry_factor)| CrackGeometryPoint {
            crack_length_m,
            geometry_factor,
        })
        .collect()
}

fn table(values: &[(f64, f64)]) -> Arc<TabulatedCrackGeometry> {
    Arc::new(TabulatedCrackGeometry::new(context(), points(values), 8).unwrap())
}

fn limits() -> CrackGrowthLimits {
    CrackGrowthLimits {
        maximum_blocks: 100,
        maximum_cycles: 100_000.0,
    }
}

fn run(geometry: CrackGeometry) -> ParisCrackGrowth {
    ParisCrackGrowth::new(law(2.0, domain()), geometry, A0, limits()).unwrap()
}

fn mixed_geometry() -> CrackGeometry {
    // Y=sqrt(a/a0) until 2a0; thereafter Y=2sqrt(a0/a). Thus the
    // rate is r0(a/a0)^2 first, then exactly 4r0 across two intervals.
    CrackGeometry::Tabulated(table(&[
        (A0, 1.0),
        (2.0 * A0, 2.0_f64.sqrt()),
        (3.0 * A0, 2.0 / 3.0_f64.sqrt()),
        (4.0 * A0, 1.0),
    ]))
}

fn cycle(cycles: f64) -> StressCycle {
    StressCycle {
        minimum_pa: MINIMUM_PA,
        maximum_pa: MAXIMUM_PA,
        cycles,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4.0e-12 * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:e} != {expected:e}"
    );
}

fn bits(state: CrackGrowthState) -> (u64, u64, u64) {
    (
        state.crack_length_m.to_bits(),
        state.cycles.to_bits(),
        state.blocks,
    )
}

#[test]
fn g1_sqrt_geometry_integrates_quadratic_crack_rate_and_retains_its_source() {
    let geometry = table(&[(A0, 1.0), (4.0 * A0, 2.0)]);
    assert_eq!(geometry.points().len(), 2);
    assert_eq!(geometry.points()[1].crack_length_m, 4.0 * A0);
    assert_eq!(
        geometry.context().source_reference,
        context().source_reference
    );
    assert_eq!(
        geometry.context().crack_length_convention,
        context().crack_length_convention
    );
    assert_eq!(
        geometry.context().nominal_stress_basis,
        context().nominal_stress_basis
    );
    close(geometry.factor_at(2.0 * A0).unwrap(), 2.0_f64.sqrt());
    let mut evolution = run(CrackGeometry::Tabulated(geometry));
    assert_eq!(evolution.geometry().length_bounds(), (A0, 4.0 * A0));
    let cycles = 250.0;
    // Y proportional to sqrt(a) makes Delta K proportional to a. With m=2,
    // da/dN=r0(a/a0)^2, whose independently integrated solution is rational.
    let expected = A0 / (1.0 - RATE * cycles / A0);
    let step = evolution.advance(cycle(cycles)).unwrap();
    close(step.after.crack_length_m, expected);
    close(step.initial_delta_k_pa_sqrt_m, reference_k());
    close(step.final_delta_k_pa_sqrt_m, reference_k() * expected / A0);
    close(
        step.final_maximum_k_pa_sqrt_m,
        1.5 * reference_k() * expected / A0,
    );
    assert_eq!(step.after.cycles, cycles);
    assert_eq!(step.after.blocks, 1);
}

#[test]
fn g1_inverse_sqrt_geometry_keeps_intensity_constant_and_growth_linear() {
    let mut evolution = run(CrackGeometry::Tabulated(table(&[
        (A0, 1.0),
        (2.0 * A0, 1.0 / 2.0_f64.sqrt()),
        (4.0 * A0, 0.5),
    ])));
    // Here Y sqrt(a) is constant: no exponential length law is applicable.
    let cycles = 1_500.0;
    let step = evolution.advance(cycle(cycles)).unwrap();
    close(step.after.crack_length_m, A0 + RATE * cycles);
    close(step.initial_delta_k_pa_sqrt_m, reference_k());
    close(step.final_delta_k_pa_sqrt_m, reference_k());
    close(step.final_maximum_k_pa_sqrt_m, 1.5 * reference_k());
    let accepted = evolution.state();
    let target = 3.5 * A0;
    let forecast = evolution
        .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, target)
        .unwrap();
    close(
        forecast.additional_cycles,
        (target - accepted.crack_length_m) / RATE,
    );
    assert_eq!(bits(evolution.state()), bits(accepted));
    evolution
        .advance(cycle(forecast.additional_cycles))
        .unwrap();
    close(evolution.state().crack_length_m, target);
}

#[test]
fn g3_piecewise_crossings_match_partition_restart_and_inverse_forecast() {
    let target = 3.5 * A0;
    // First integral: a0/r0 * (1-a0/(2a0)). The remaining 1.5a0 grows at
    // 4r0, so this three-interval exposure is exactly 500+375=875 cycles.
    let expected_cycles = 0.5 * A0 / RATE + 1.5 * A0 / (4.0 * RATE);
    let mut whole = run(mixed_geometry());
    let before = whole.state();
    let forecast = whole
        .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, target)
        .unwrap();
    close(forecast.additional_cycles, expected_cycles);
    assert_eq!(bits(whole.state()), bits(before));
    whole.advance(cycle(expected_cycles)).unwrap();
    close(whole.state().crack_length_m, target);

    let mut partitioned = run(mixed_geometry());
    for count in [125.0, 375.0, 100.0, 275.0] {
        partitioned.advance(cycle(count)).unwrap();
    }
    close(partitioned.state().crack_length_m, target);
    assert_eq!(partitioned.state().cycles, whole.state().cycles);

    let mut prefix = run(mixed_geometry());
    prefix.advance(cycle(250.0)).unwrap();
    let accepted = prefix.state();
    let mut resumed = prefix.clone();
    let remaining = resumed
        .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, target)
        .unwrap();
    close(remaining.additional_cycles, expected_cycles - 250.0);
    assert_eq!(bits(resumed.state()), bits(accepted));
    resumed.advance(cycle(remaining.additional_cycles)).unwrap();
    prefix.advance(cycle(remaining.additional_cycles)).unwrap();
    assert_eq!(bits(resumed.state()), bits(prefix.state()));
    close(resumed.state().crack_length_m, target);
}

#[test]
fn g1_constant_tables_match_constant_geometry_for_three_paris_exponents() {
    let geometry = CrackGeometry::Tabulated(table(&[
        (A0, 1.0),
        (2.0 * A0, 1.0),
        (3.0 * A0, 1.0),
        (4.0 * A0, 1.0),
    ]));
    for exponent in [1.0, 2.0, 4.0] {
        let model = law(exponent, domain());
        let mut constant = ParisCrackGrowth::new(
            model.clone(),
            ConstantGeometryCrack {
                geometry_factor: 1.0,
                minimum_length_m: A0,
                maximum_length_m: 4.0 * A0,
            },
            A0,
            limits(),
        )
        .unwrap();
        let mut tabulated = ParisCrackGrowth::new(model, geometry.clone(), A0, limits()).unwrap();
        let target = 3.5 * A0;
        let forecast = constant
            .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, target)
            .unwrap();
        let table_forecast = tabulated
            .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, target)
            .unwrap();
        close(table_forecast.additional_cycles, forecast.additional_cycles);
        let a = constant.advance(cycle(forecast.additional_cycles)).unwrap();
        let b = tabulated
            .advance(cycle(forecast.additional_cycles))
            .unwrap();
        close(a.after.crack_length_m, target);
        close(b.after.crack_length_m, target);
        close(b.final_delta_k_pa_sqrt_m, a.final_delta_k_pa_sqrt_m);
        assert_eq!(a.after.cycles, b.after.cycles);
        assert_eq!(a.after.blocks, b.after.blocks);
    }
}

#[test]
fn g1_forecasts_land_exactly_on_each_original_knot_without_extending_support() {
    for prior_cycles in [0.0, 50.0] {
        for target in [2.0 * A0, 3.0 * A0, 4.0 * A0] {
            let mut evolution = run(mixed_geometry());
            if prior_cycles > 0.0 {
                evolution.advance(cycle(prior_cycles)).unwrap();
            }
            let before = evolution.state();
            let forecast = evolution
                .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, target)
                .unwrap();
            if target == 4.0 * A0 {
                let excess = f64::from_bits(forecast.additional_cycles.to_bits() + 1);
                assert!(matches!(
                    evolution.advance(cycle(excess)),
                    Err(CrackGrowthError::GeometrySupportExhausted { .. })
                ));
                assert_eq!(bits(evolution.state()), bits(before));
            }
            evolution
                .advance(cycle(forecast.additional_cycles))
                .unwrap();
            assert_eq!(evolution.state().crack_length_m.to_bits(), target.to_bits());
        }
    }
}

#[test]
fn g0_interior_intensity_peak_refuses_despite_supported_endpoints() {
    let geometry = CrackGeometry::Tabulated(table(&[
        (A0, 1.0),
        (2.0 * A0, 2.0_f64.sqrt()),
        (4.0 * A0, 0.5),
    ]));
    // Endpoint Delta K is Kref at both a0 and 4a0; at 2a0 it is 2Kref.
    // Endpoint checks alone therefore miss both kinds of interior violation.
    for length in [A0, 4.0 * A0] {
        let delta_k = geometry.factor_at(length).unwrap()
            * (MAXIMUM_PA - MINIMUM_PA)
            * (std::f64::consts::PI * length).sqrt();
        close(delta_k, reference_k());
    }
    let mut maximum_k = domain();
    maximum_k.maximum_k_pa_sqrt_m = 2.25 * reference_k();
    let mut fit_range = domain();
    fit_range.maximum_delta_k_pa_sqrt_m = 1.5 * reference_k();
    for restricted in [maximum_k, fit_range] {
        let mut evolution =
            ParisCrackGrowth::new(law(2.0, restricted), geometry.clone(), A0, limits()).unwrap();
        let before = evolution.state();
        let refusal = evolution
            .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, 4.0 * A0)
            .unwrap_err();
        assert!(matches!(
            refusal,
            CrackGrowthError::MaximumKExceeded { .. }
                | CrackGrowthError::DeltaKOutsideDomain { .. }
        ));
        // The second segment has da/dN=16*r0*a0^2/a^2. Its exposure from
        // 2a0 to 4a0 is (7/6)*a0/r0, after the first segment's a0/(2r0).
        let refusal = evolution
            .advance(cycle((0.5 + 7.0 / 6.0) * A0 / RATE))
            .unwrap_err();
        assert!(matches!(
            refusal,
            CrackGrowthError::MaximumKExceeded { .. }
                | CrackGrowthError::DeltaKOutsideDomain { .. }
        ));
        assert_eq!(bits(evolution.state()), bits(before));
    }
}

#[test]
fn g0_geometry_requires_source_ordered_positive_points_and_bounded_support() {
    let valid = [(A0, 1.0), (2.0 * A0, 2.0)];
    for missing in 0..3 {
        let mut incomplete = context();
        match missing {
            0 => incomplete.source_reference = "  ".into(),
            1 => incomplete.crack_length_convention.clear(),
            _ => incomplete.nominal_stress_basis.clear(),
        }
        assert!(TabulatedCrackGeometry::new(incomplete, points(&valid), 2).is_err());
    }
    assert!(TabulatedCrackGeometry::new(context(), Vec::new(), 2).is_err());
    assert!(TabulatedCrackGeometry::new(context(), points(&valid[..1]), 2).is_err());
    assert!(TabulatedCrackGeometry::new(context(), points(&valid), 1).is_err());
    for invalid in [
        [(A0, 1.0), (A0, 2.0)],
        [(2.0 * A0, 1.0), (A0, 2.0)],
        [(0.0, 1.0), (A0, 2.0)],
        [(A0, 1.0), (f64::INFINITY, 2.0)],
        [(f64::NAN, 1.0), (2.0 * A0, 2.0)],
        [(A0, 0.0), (2.0 * A0, 2.0)],
        [(A0, -1.0), (2.0 * A0, 2.0)],
        [(A0, f64::NAN), (2.0 * A0, 2.0)],
    ] {
        assert!(TabulatedCrackGeometry::new(context(), points(&invalid), 2).is_err());
    }
    let supported = table(&valid);
    for outside in [0.5 * A0, 2.1 * A0, f64::NAN] {
        assert!(supported.factor_at(outside).is_err());
        assert!(
            ParisCrackGrowth::new(
                law(2.0, domain()),
                CrackGeometry::Tabulated(supported.clone()),
                outside,
                limits(),
            )
            .is_err()
        );
    }
}

#[test]
fn g0_exhausted_geometry_reports_unconsumed_cycles_without_partial_growth() {
    let mut evolution = run(CrackGeometry::Tabulated(table(&[
        (A0, 1.0),
        (2.0 * A0, 1.0 / 2.0_f64.sqrt()),
    ])));
    let before = evolution.state();
    let error = evolution.advance(cycle(1_500.0)).unwrap_err();
    match error {
        CrackGrowthError::GeometrySupportExhausted {
            length_m,
            remaining_cycles,
        } => {
            close(length_m, 2.0 * A0);
            close(remaining_cycles, 500.0);
        }
        other => panic!("expected explicit geometry exhaustion, received {other:?}"),
    }
    assert_eq!(bits(evolution.state()), bits(before));
    assert!(
        evolution
            .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, 2.5 * A0)
            .is_err()
    );
    assert_eq!(bits(evolution.state()), bits(before));
    // The squared-root table values and evaluated reference intensity are
    // rounded floats, so the supported endpoint's model count may differ by
    // a few ULPs from the independent mathematical value 1000. Use the model's
    // actual endpoint count instead of extending its support with a tolerance.
    let endpoint_cycles = evolution
        .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, 2.0 * A0)
        .unwrap()
        .additional_cycles;
    close(endpoint_cycles, 1_000.0);
    let excess = f64::from_bits(endpoint_cycles.to_bits() + 1);
    assert!(matches!(
        evolution.advance(cycle(excess)),
        Err(CrackGrowthError::GeometrySupportExhausted { .. })
    ));
    assert_eq!(bits(evolution.state()), bits(before));
    evolution.advance(cycle(endpoint_cycles)).unwrap();
    close(evolution.state().crack_length_m, 2.0 * A0);
    let endpoint = evolution.state();
    assert!(evolution.advance(cycle(1.0)).is_err());
    assert_eq!(bits(evolution.state()), bits(endpoint));
}

#[test]
fn g4_cancellation_during_geometry_traversal_preserves_the_accepted_state() {
    let mut accepted = run(mixed_geometry());
    accepted.advance(cycle(50.0)).unwrap();
    let before = accepted.state();
    let target = 3.5 * A0;
    let remaining = 825.0;
    let mut complete = accepted.clone();
    let mut polls = 0;
    complete
        .advance_with_cancel(cycle(remaining), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(
        polls >= 3,
        "crossed intervals require interior cancellation polling"
    );
    let mut attempt = accepted.clone();
    let mut calls = 0;
    assert!(matches!(
        attempt.advance_with_cancel(cycle(remaining), &mut || {
            calls += 1;
            calls == polls - 1
        }),
        Err(CrackGrowthError::Cancelled)
    ));
    assert!(calls > 1);
    assert_eq!(bits(attempt.state()), bits(before));
    attempt.advance(cycle(remaining)).unwrap();
    assert_eq!(bits(attempt.state()), bits(complete.state()));

    let mut forecast_polls = 0;
    let forecast = accepted
        .cycles_to_length_with_cancel(MINIMUM_PA, MAXIMUM_PA, target, &mut || {
            forecast_polls += 1;
            false
        })
        .unwrap();
    close(forecast.additional_cycles, remaining);
    assert!(forecast_polls >= 3);
    let mut calls = 0;
    assert!(matches!(
        accepted.cycles_to_length_with_cancel(MINIMUM_PA, MAXIMUM_PA, target, &mut || {
            calls += 1;
            calls == forecast_polls - 1
        }),
        Err(CrackGrowthError::Cancelled)
    ));
    assert!(calls > 1);
    assert_eq!(bits(accepted.state()), bits(before));
}

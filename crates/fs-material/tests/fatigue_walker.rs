//! Normalized tensile Walker growth against independent analytical relations.
//! All coefficients and geometry tables are synthetic implementation fixtures.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGeometryPoint, CrackGeometryTableContext, CrackGrowthError,
    CrackGrowthLimits, ParisCrackGrowth, ParisDomain, ParisLaw, ParisParameters, StressCycle,
    TabulatedCrackGeometry, WalkerParameters,
};

const A0: f64 = 1.0e-3;
const REFERENCE_K: f64 = 1.0e6;
const RATE: f64 = 1.0e-7;

fn scope() -> FatigueScope {
    FatigueScope {
        source_reference: "synthetic normalized Walker fit; not measured".into(),
        material: "synthetic coupon".into(),
        process: "synthetic fixed process".into(),
        surface_and_geometry: "declared tensile mode-I geometry".into(),
        environment: "synthetic fixed environment".into(),
        load_basis: "full stress-intensity range; full cycles".into(),
        temperature_k: 300.0,
    }
}

fn parameters(exponent: f64) -> ParisParameters {
    ParisParameters {
        reference_delta_k_pa_sqrt_m: REFERENCE_K,
        reference_rate_m_per_cycle: RATE,
        exponent,
    }
}

fn domain() -> ParisDomain {
    ParisDomain {
        minimum_delta_k_pa_sqrt_m: 0.1 * REFERENCE_K,
        maximum_delta_k_pa_sqrt_m: 100.0 * REFERENCE_K,
        minimum_load_ratio: 0.0,
        maximum_load_ratio: 0.9,
        maximum_k_pa_sqrt_m: 1_000.0 * REFERENCE_K,
    }
}

fn walker(reference_load_ratio: f64, gamma: f64, exponent: f64) -> Arc<ParisLaw> {
    Arc::new(
        ParisLaw::new_walker(
            scope(),
            parameters(exponent),
            domain(),
            WalkerParameters {
                reference_load_ratio,
                gamma,
            },
        )
        .unwrap(),
    )
}

fn limits() -> CrackGrowthLimits {
    CrackGrowthLimits {
        maximum_blocks: 100,
        maximum_cycles: 100_000.0,
    }
}

fn evolution(law: Arc<ParisLaw>) -> ParisCrackGrowth {
    ParisCrackGrowth::new(
        law,
        ConstantGeometryCrack {
            geometry_factor: 1.0,
            minimum_length_m: A0,
            maximum_length_m: 10.0 * A0,
        },
        A0,
        limits(),
    )
    .unwrap()
}

fn cycle(ratio: f64, cycles: f64) -> StressCycle {
    let range = REFERENCE_K / (std::f64::consts::PI * A0).sqrt();
    let maximum_pa = range / (1.0 - ratio);
    StressCycle {
        minimum_pa: ratio * maximum_pa,
        maximum_pa,
        cycles,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4.0e-12 * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn reference_pair_is_exact_and_equal_ranges_follow_the_original_walker_equation() {
    let law = walker(0.2, 0.4, 3.0);
    assert_eq!(law.walker().unwrap().reference_load_ratio, 0.2);
    assert_eq!(law.walker().unwrap().gamma, 0.4);
    assert_eq!(law.scope(), &scope());
    assert!(law.source().is_none());
    assert_eq!(
        law.growth_rate(REFERENCE_K, 0.2).unwrap().to_bits(),
        RATE.to_bits()
    );
    // Independently use AFGROW's K_max form, not the implementation's log shift.
    let reference_equivalent = (1.0_f64 - 0.2).powf(0.4) * REFERENCE_K / (1.0 - 0.2);
    for ratio in [0.0, 0.2, 0.6] {
        let delta_k = 2.0 * REFERENCE_K;
        let maximum_k = delta_k / (1.0 - ratio);
        let equivalent = (1.0_f64 - ratio).powf(0.4) * maximum_k;
        let expected = RATE * (equivalent / reference_equivalent).powi(3);
        close(law.growth_rate(delta_k, ratio).unwrap(), expected);
    }
    assert!(law.growth_rate(REFERENCE_K, 0.6).unwrap() > RATE);
    assert!(law.growth_rate(REFERENCE_K, 0.0).unwrap() < RATE);
}

#[test]
fn gamma_limits_preserve_paris_bits_and_recover_maximum_intensity_scaling() {
    let plain = Arc::new(ParisLaw::new(scope(), parameters(2.0), domain()).unwrap());
    let gamma_one = walker(0.25, 1.0, 2.0);
    assert!(plain.walker().is_none());
    for delta_k in [0.5 * REFERENCE_K, REFERENCE_K, 2.0 * REFERENCE_K] {
        for ratio in [0.0, 0.25, 0.6] {
            assert_eq!(
                plain.growth_rate(delta_k, ratio).unwrap().to_bits(),
                gamma_one.growth_rate(delta_k, ratio).unwrap().to_bits()
            );
        }
    }
    let mut plain_growth = evolution(plain);
    let mut gamma_one_growth = evolution(gamma_one);
    for block in [cycle(0.0, 250.0), cycle(0.5, 750.0)] {
        assert_eq!(
            plain_growth.advance(block).unwrap(),
            gamma_one_growth.advance(block).unwrap()
        );
    }
    let gamma_zero = walker(0.25, 0.0, 2.0);
    let maximum_k = 2.0 * REFERENCE_K;
    let expected = RATE * (maximum_k / (REFERENCE_K / 0.75)).powi(2);
    for ratio in [0.0, 0.25, 0.5, 0.75] {
        close(
            gamma_zero
                .growth_rate(maximum_k * (1.0 - ratio), ratio)
                .unwrap(),
            expected,
        );
    }
    let near_reference = walker(0.8, 0.5, 2.0);
    let ratio = 0.8 + 1.0e-8;
    close(
        near_reference.growth_rate(REFERENCE_K, ratio).unwrap(),
        RATE * ((1.0 - 0.8) / (1.0 - ratio)),
    );
}

#[test]
fn constant_geometry_growth_and_remaining_exposure_have_the_walker_time_scale() {
    let mut growth = evolution(walker(0.25, 0.5, 2.0));
    let block = cycle(0.5, 2_000.0);
    // For n=2 and gamma=1/2, R=.5 versus R_ref=.25 multiplies rate by1.5.
    let coefficient = 1.5 * RATE / A0;
    let expected = A0 * (coefficient * block.cycles).exp();
    let step = growth.advance(block).unwrap();
    close(step.after.crack_length_m, expected);
    close(step.initial_delta_k_pa_sqrt_m, REFERENCE_K);
    close(
        step.final_maximum_k_pa_sqrt_m,
        2.0 * REFERENCE_K * (expected / A0).sqrt(),
    );
    let before = growth.state();
    let target = 2.5 * A0;
    let forecast = growth
        .cycles_to_length(block.minimum_pa, block.maximum_pa, target)
        .unwrap();
    close(
        forecast.additional_cycles,
        (target / A0).ln() / coefficient - block.cycles,
    );
    assert_eq!(growth.state(), before);
    growth
        .advance(StressCycle {
            cycles: forecast.additional_cycles,
            ..block
        })
        .unwrap();
    close(growth.state().crack_length_m, target);
}

#[test]
fn evolving_geometry_uses_the_walker_rate_across_every_knot_and_exact_final_support() {
    // First Y=sqrt(a/a0), then Y=2sqrt(a0/a): quadratic, then constant rate.
    let points = [
        (A0, 1.0),
        (2.0 * A0, 2.0_f64.sqrt()),
        (3.0 * A0, 2.0 / 3.0_f64.sqrt()),
        (4.0 * A0, 1.0),
    ]
    .into_iter()
    .map(|(crack_length_m, geometry_factor)| CrackGeometryPoint {
        crack_length_m,
        geometry_factor,
    })
    .collect();
    let table = TabulatedCrackGeometry::new(
        CrackGeometryTableContext {
            source_reference: "synthetic exact piecewise geometry".into(),
            crack_length_convention: "physical half crack length in metres".into(),
            nominal_stress_basis: "remote tensile stress in Pa".into(),
        },
        points,
        4,
    )
    .unwrap();
    let mut growth = ParisCrackGrowth::new(walker(0.25, 0.5, 2.0), table, A0, limits()).unwrap();
    let block = cycle(0.5, 1.0);
    let forecast = growth
        .cycles_to_length(block.minimum_pa, block.maximum_pa, 4.0 * A0)
        .unwrap();
    // First interval costs a0/(2*1.5*r0), remaining2a0 at rate4*1.5*r0.
    close(forecast.additional_cycles, A0 / (1.5 * RATE));
    let step = growth
        .advance(StressCycle {
            cycles: forecast.additional_cycles,
            ..block
        })
        .unwrap();
    assert_eq!(step.after.crack_length_m, 4.0 * A0);
    close(step.final_delta_k_pa_sqrt_m, 2.0 * REFERENCE_K);
}

#[test]
fn changing_ratio_blocks_preserve_partition_restart_and_atomic_rollback() {
    let law = walker(0.25, 0.5, 2.0);
    let blocks = [
        cycle(0.0, 1_000.0),
        cycle(0.5, 2_000.0),
        cycle(0.75, 1_000.0),
    ];
    let mut whole = evolution(law.clone());
    whole.advance_block(&blocks, &mut || false).unwrap();
    let expected = A0 * (RATE / A0 * (0.75 * 1_000.0 + 1.5 * 2_000.0 + 3.0 * 1_000.0)).exp();
    close(whole.state().crack_length_m, expected);
    let mut prefix = evolution(law.clone());
    prefix.advance(blocks[0]).unwrap();
    let mut restart = prefix.clone();
    prefix.advance_block(&blocks[1..], &mut || false).unwrap();
    for block in &blocks[1..] {
        restart.advance(*block).unwrap();
    }
    assert_eq!(prefix.state(), restart.state());
    assert_eq!(prefix.state(), whole.state());
    let mut partitioned = evolution(law);
    for block in blocks {
        for _ in 0..2 {
            partitioned
                .advance(StressCycle {
                    cycles: block.cycles / 2.0,
                    ..block
                })
                .unwrap();
        }
    }
    close(partitioned.state().crack_length_m, expected);
    assert_eq!(partitioned.state().cycles, whole.state().cycles);
    let before = prefix.state();
    let unsupported = cycle(0.95, 1.0);
    assert!(matches!(
        prefix.advance_block(&[blocks[0], unsupported], &mut || false),
        Err(CrackGrowthError::LoadRatioOutsideDomain { .. })
    ));
    assert_eq!(prefix.state(), before);
    let mut polls = 0;
    assert!(matches!(
        prefix.advance_block(&blocks, &mut || {
            polls += 1;
            polls >= 5
        }),
        Err(CrackGrowthError::Cancelled)
    ));
    assert_eq!(prefix.state(), before);
}

#[test]
fn walker_calibration_and_actual_loading_support_cannot_be_bypassed() {
    for gamma in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            ParisLaw::new_walker(
                scope(),
                parameters(2.0),
                domain(),
                WalkerParameters {
                    reference_load_ratio: 0.25,
                    gamma
                }
            )
            .is_err()
        );
    }
    for reference_load_ratio in [-0.1, 0.95, 1.0, f64::NAN] {
        assert!(
            ParisLaw::new_walker(
                scope(),
                parameters(2.0),
                domain(),
                WalkerParameters {
                    reference_load_ratio,
                    gamma: 0.5
                }
            )
            .is_err()
        );
    }
    let mut restricted = domain();
    restricted.maximum_k_pa_sqrt_m = 2.0 * REFERENCE_K;
    assert!(matches!(
        ParisLaw::new_walker(
            scope(),
            parameters(2.0),
            restricted,
            WalkerParameters {
                reference_load_ratio: 0.5,
                gamma: 0.5
            }
        ),
        Err(CrackGrowthError::MaximumKExceeded { .. })
    ));
    let law = walker(0.0, 0.0, 2.0);
    // A Walker-equivalent range above the lower fit bound must not admit an
    // actual Delta-K below it. R and upper K also keep their original meaning.
    assert!(matches!(
        law.growth_rate(0.09 * REFERENCE_K, 0.5),
        Err(CrackGrowthError::DeltaKOutsideDomain { .. })
    ));
    assert!(matches!(
        law.growth_rate(REFERENCE_K, 0.95),
        Err(CrackGrowthError::LoadRatioOutsideDomain { .. })
    ));
    let mut static_growth = evolution(law);
    let block = cycle(0.5, 1.0);
    let step = static_growth
        .advance(StressCycle {
            minimum_pa: block.maximum_pa,
            ..block
        })
        .unwrap();
    assert_eq!(step.after.crack_length_m, A0);
    assert_eq!(step.initial_delta_k_pa_sqrt_m, 0.0);
    assert!(
        static_growth
            .cycles_to_length(block.maximum_pa, block.maximum_pa, 2.0 * A0)
            .is_err()
    );
}

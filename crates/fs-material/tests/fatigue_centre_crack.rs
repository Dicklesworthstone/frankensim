//! Handbook centre-crack geometry and its effect on bounded Paris evolution.
//! The Paris coefficients are synthetic, not measured material or service life.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGeometry, CrackGrowthError, CrackGrowthLimits,
    FiniteWidthCentreCrack, ParisCrackGrowth, ParisDomain, ParisLaw, ParisParameters, StressCycle,
};

const WIDTH: f64 = 0.1;
const INITIAL: f64 = 0.005;
const TARGET: f64 = 0.03;
const MINIMUM_PA: f64 = 1.0e6;
const MAXIMUM_PA: f64 = 11.0e6;
const REFERENCE_K: f64 = 1.0e6;
const RATE: f64 = 1.0e-8;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3.0e-14 * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}

fn lengths(intervals: usize) -> Vec<f64> {
    (0..=intervals)
        .map(|i| INITIAL + (TARGET - INITIAL) * (i as f64 / intervals as f64))
        .collect()
}

fn evolution(geometry: impl Into<CrackGeometry>) -> ParisCrackGrowth {
    let law = ParisLaw::new(
        FatigueScope {
            source_reference: "synthetic Paris coefficients for geometry verification".into(),
            material: "synthetic elastic plate".into(),
            process: "synthetic fixed process".into(),
            surface_and_geometry: "centred through crack under remote uniform tension".into(),
            environment: "synthetic fixed environment".into(),
            load_basis: "full tensile stress-intensity range; full cycles".into(),
            temperature_k: 300.0,
        },
        ParisParameters {
            reference_delta_k_pa_sqrt_m: REFERENCE_K,
            reference_rate_m_per_cycle: RATE,
            exponent: 3.0,
        },
        ParisDomain {
            minimum_delta_k_pa_sqrt_m: 1.0,
            maximum_delta_k_pa_sqrt_m: 1.0e8,
            minimum_load_ratio: 0.0,
            maximum_load_ratio: 0.8,
            maximum_k_pa_sqrt_m: 1.0e9,
        },
    )
    .unwrap();
    ParisCrackGrowth::new(
        Arc::new(law),
        geometry,
        INITIAL,
        CrackGrowthLimits {
            maximum_blocks: 10,
            maximum_cycles: 1.0e9,
        },
    )
    .unwrap()
}

#[test]
fn handbook_factor_uses_full_width_half_length_and_is_scale_invariant() {
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    assert_eq!(model.plate_width_m(), WIDTH);
    assert_eq!(model.maximum_half_length_m(), 0.4 * WIDTH);
    // a/W=1/4 means full-crack/full-width=1/2, and Y=2^(1/4).
    let expected = 2.0_f64.sqrt().sqrt();
    close(model.factor_at(0.25 * WIDTH).unwrap(), expected);
    let scaled = FiniteWidthCentreCrack::new(100.0).unwrap();
    close(scaled.factor_at(25.0).unwrap(), expected);
    close(
        FiniteWidthCentreCrack::new(1.0e12)
            .unwrap()
            .factor_at(0.025)
            .unwrap(),
        1.0,
    );
    let requested = vec![INITIAL, 0.25 * WIDTH, model.maximum_half_length_m()];
    let table = model.tabulate(requested.clone(), 3).unwrap();
    for (point, length) in table.points().iter().zip(requested) {
        assert_eq!(point.crack_length_m, length);
        assert_eq!(point.geometry_factor, model.factor_at(length).unwrap());
    }
    assert!(table.context().source_reference.contains("AFGROW"));
    assert!(table.context().source_reference.contains("NASA TN D-5140"));
    assert!(
        table
            .context()
            .source_reference
            .contains(&format!("{WIDTH:.17e}"))
    );
    assert!(
        table
            .context()
            .crack_length_convention
            .contains("half-length")
    );
    assert!(
        table
            .context()
            .nominal_stress_basis
            .contains("gross-section")
    );
}

#[test]
fn support_original_point_order_work_limit_and_cancellation_are_enforced() {
    for width in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(FiniteWidthCentreCrack::new(width).is_err());
    }
    assert!(FiniteWidthCentreCrack::new(f64::from_bits(1)).is_err());
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    let maximum = model.maximum_half_length_m();
    assert!(model.factor_at(maximum).is_ok());
    for length in [0.0, -INITIAL, f64::NAN, f64::INFINITY, maximum.next_up()] {
        assert!(model.factor_at(length).is_err());
    }
    for points in [
        vec![],
        vec![INITIAL],
        vec![TARGET, INITIAL],
        vec![INITIAL, INITIAL],
        vec![INITIAL, maximum.next_up()],
    ] {
        assert!(model.tabulate(points, 3).is_err());
    }
    assert!(matches!(
        model.tabulate(vec![INITIAL, TARGET], 1),
        Err(CrackGrowthError::Budget)
    ));
    assert!(matches!(
        model.tabulate_with_cancel(vec![INITIAL, TARGET], 2, &mut || true),
        Err(CrackGrowthError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        model.tabulate_with_cancel(lengths(8), 9, &mut || {
            polls += 1;
            polls >= 4
        }),
        Err(CrackGrowthError::Cancelled)
    ));
    assert!(model.tabulate(vec![INITIAL, TARGET], 2).is_ok());
}

fn direct_secant_cycles() -> f64 {
    // Composite Simpson integration of dN/da using the unsampled handbook Y.
    // This reference uses platform math, not geometry interpolation or the
    // production analytical power-law integral.
    let reciprocal_rate = |a: f64| {
        let y = (std::f64::consts::PI * a / WIDTH).cos().sqrt().recip();
        let delta_k = (MAXIMUM_PA - MINIMUM_PA) * (std::f64::consts::PI * a).sqrt() * y;
        (RATE * (delta_k / REFERENCE_K).powi(3)).recip()
    };
    let panels = 16_384;
    let step = (TARGET - INITIAL) / panels as f64;
    let mut sum = reciprocal_rate(INITIAL) + reciprocal_rate(TARGET);
    for i in 1..panels {
        let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += weight * reciprocal_rate(INITIAL + i as f64 * step);
    }
    sum * step / 3.0
}

#[test]
fn refined_table_life_converges_to_independent_direct_secant_integration() {
    let reference_cycles = direct_secant_cycles();
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    let mut previous_error = f64::INFINITY;
    for intervals in [2, 8, 32, 128] {
        let table = model.tabulate(lengths(intervals), intervals + 1).unwrap();
        let forecast = evolution(table)
            .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, TARGET)
            .unwrap();
        let error = (forecast.additional_cycles / reference_cycles - 1.0).abs();
        assert!(error < previous_error, "refinement error {error:e}");
        previous_error = error;
    }
    assert!(
        previous_error < 2.0e-5,
        "relative life error {previous_error:e}"
    );
}

#[test]
fn finite_width_increases_growth_and_reduces_cycles_to_the_same_length() {
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    let mut finite = evolution(model.tabulate(lengths(128), 129).unwrap());
    let mut infinite = evolution(ConstantGeometryCrack {
        geometry_factor: 1.0,
        minimum_length_m: INITIAL,
        maximum_length_m: TARGET,
    });
    let finite_life = finite
        .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, TARGET)
        .unwrap()
        .additional_cycles;
    let infinite_life = infinite
        .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, TARGET)
        .unwrap()
        .additional_cycles;
    assert!(finite_life < infinite_life);
    let cycle = StressCycle {
        minimum_pa: MINIMUM_PA,
        maximum_pa: MAXIMUM_PA,
        cycles: finite_life / 4.0,
    };
    let finite_step = finite.advance(cycle).unwrap();
    let infinite_step = infinite.advance(cycle).unwrap();
    assert!(finite_step.after.crack_length_m > infinite_step.after.crack_length_m);
    assert!(finite_step.final_delta_k_pa_sqrt_m > infinite_step.final_delta_k_pa_sqrt_m);
}

#[test]
fn automatic_table_obeys_the_relative_factor_envelope_between_knots() {
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    let minimum = WIDTH * 0.01;
    let maximum = model.maximum_half_length_m();
    let tolerance = 1.0e-4;
    let table = model
        .tabulate_to_tolerance(minimum, maximum, tolerance, 512)
        .unwrap();
    assert_eq!(table.length_bounds(), (minimum, maximum));
    assert!(table.points().len() < 512);
    assert_eq!(
        table,
        model
            .tabulate_to_tolerance(minimum, maximum, tolerance, 512)
            .unwrap()
    );
    assert!(
        table
            .context()
            .source_reference
            .contains(&format!("relative Y tolerance={tolerance:.17e}"))
    );
    for interval in table.points().windows(2) {
        let left = interval[0].crack_length_m;
        let right = interval[1].crack_length_m;
        assert!(left < right);
        for i in 0..=64 {
            let a = match i {
                0 => left,
                64 => right,
                _ => left * (right / left).powf(i as f64 / 64.0),
            };
            let direct = (std::f64::consts::PI * a / WIDTH).cos().sqrt().recip();
            let relative = table.factor_at(a).unwrap() / direct - 1.0;
            assert!(relative >= -3.0e-14, "lost convex envelope: {relative:e}");
            assert!(relative <= tolerance + 3.0e-14, "factor error {relative:e}");
        }
    }
}

#[test]
fn tighter_factor_tolerance_improves_the_independent_paris_life_integral() {
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    let reference = direct_secant_cycles();
    let mut previous_error = f64::INFINITY;
    let mut previous_points = 0;
    for tolerance in [1.0e-2, 1.0e-4, 1.0e-6] {
        let table = model
            .tabulate_to_tolerance(INITIAL, TARGET, tolerance, 1025)
            .unwrap();
        assert!(table.points().len() > previous_points);
        previous_points = table.points().len();
        let cycles = evolution(table)
            .cycles_to_length(MINIMUM_PA, MAXIMUM_PA, TARGET)
            .unwrap()
            .additional_cycles;
        // For this fixed synthetic m=3 fit alone, the Y envelope bounds
        // the inverse-rate integral between (1+tolerance)^-3 and 1.
        assert!(cycles >= reference / (1.0 + tolerance).powi(3) * (1.0 - 1.0e-11));
        assert!(cycles <= reference * (1.0 + 1.0e-11));
        let error = (cycles / reference - 1.0).abs();
        assert!(error < previous_error);
        previous_error = error;
    }
    assert!(previous_error < 3.0e-6);
}

#[test]
fn automatic_sampling_preserves_scaling_and_resolves_narrow_supported_intervals() {
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    let table = model
        .tabulate_to_tolerance(INITIAL, TARGET, 1.0e-4, 512)
        .unwrap();
    for scale in [1.0e-3, 1.0e3] {
        let scaled = FiniteWidthCentreCrack::new(WIDTH * scale)
            .unwrap()
            .tabulate_to_tolerance(INITIAL * scale, TARGET * scale, 1.0e-4, 512)
            .unwrap();
        assert_eq!(scaled.points().len(), table.points().len());
        for (point, original) in scaled.points().iter().zip(table.points()) {
            close(point.crack_length_m / scale, original.crack_length_m);
            close(point.geometry_factor, original.geometry_factor);
        }
    }
    let maximum = model.maximum_half_length_m();
    let adjacent = maximum.next_down();
    let narrow = model
        .tabulate_to_tolerance(adjacent, maximum, 1.0e-30, 2)
        .unwrap();
    assert_eq!(narrow.length_bounds(), (adjacent, maximum));
    assert_eq!(narrow.points().len(), 2);
    assert!(matches!(
        model.tabulate_to_tolerance(adjacent, maximum, 1.0e-40, 3),
        Err(CrackGrowthError::Arithmetic)
    ));
    assert_eq!(
        model
            .tabulate_to_tolerance(WIDTH * 1.0e-12, WIDTH * 2.0e-12, 1.0e-12, 2)
            .unwrap()
            .points()
            .len(),
        2
    );
}

#[test]
fn automatic_sampling_refuses_bad_domains_exhausted_budgets_and_cancellation() {
    let model = FiniteWidthCentreCrack::new(WIDTH).unwrap();
    for tolerance in [0.0, -1.0e-3, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .tabulate_to_tolerance(INITIAL, TARGET, tolerance, 64)
                .is_err()
        );
    }
    for (minimum, maximum) in [
        (INITIAL, INITIAL),
        (TARGET, INITIAL),
        (0.0, TARGET),
        (INITIAL, model.maximum_half_length_m().next_up()),
    ] {
        assert!(
            model
                .tabulate_to_tolerance(minimum, maximum, 1.0e-3, 64)
                .is_err()
        );
    }
    for budget in [1, 2, 8] {
        let mut polls = 0;
        let result =
            model.tabulate_to_tolerance_with_cancel(INITIAL, TARGET, 1.0e-8, budget, &mut || {
                polls += 1;
                false
            });
        assert!(matches!(result, Err(CrackGrowthError::Budget)));
        assert!(polls <= 2 * budget + 1, "refinement exceeded bounded work");
    }
    for cancel_after in [1, 5] {
        let mut polls = 0;
        let result =
            model.tabulate_to_tolerance_with_cancel(INITIAL, TARGET, 1.0e-6, 1025, &mut || {
                polls += 1;
                polls >= cancel_after
            });
        assert!(matches!(result, Err(CrackGrowthError::Cancelled)));
    }
}

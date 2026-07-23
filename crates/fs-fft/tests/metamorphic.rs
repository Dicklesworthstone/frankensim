//! Gauntlet G3 relations for the production FFT plan.
//!
//! Existing oracle, round-trip, and golden tests remain the primary numeric
//! and bit-contract pins; this battery adds shrinkable declared relations.

use fs_fft::{C64, Fft};
use fs_propcheck::metamorphic::{
    RelationCase, RelationObservation, Tolerance, check_relation, unit_rescaling,
};

fn minimum_complex_margin(
    tolerance: Tolerance,
    reference: &[C64],
    candidate: &[C64],
) -> RelationObservation {
    let margin = reference
        .iter()
        .zip(candidate)
        .flat_map(|(reference, candidate)| {
            [
                tolerance
                    .evaluate_scalar(reference.re, candidate.re)
                    .margin(),
                tolerance
                    .evaluate_scalar(reference.im, candidate.im)
                    .margin(),
            ]
        })
        .fold(f64::MAX, f64::min);
    RelationObservation::new(
        margin,
        "every complex FFT component obeys the declared signal rescaling",
    )
}

#[test]
fn g3_forward_fft_is_equivariant_under_signal_rescaling() {
    const N: usize = 8;
    let plan = Fft::new(N);
    let operator = |input: &Vec<(f64, f64)>| {
        let mut data = vec![C64::default(); N];
        for (slot, &(re, im)) in data.iter_mut().zip(input) {
            *slot = C64::new(re, im);
        }
        let mut scratch = vec![C64::default(); N];
        plan.forward(&mut data, &mut scratch);
        data
    };
    let relation = unit_rescaling(
        "forward-signal-scale-equivariance",
        Tolerance::AbsoluteRelative {
            max_abs: 1.0e-13,
            max_relative: 1.0e-13,
        },
        |input: &Vec<(f64, f64)>, &exponent: &i64| {
            let scale = fs_math::det::powi(2.0, exponent as i32);
            input
                .iter()
                .map(|&(re, im)| (re * scale, im * scale))
                .collect()
        },
        |base: &Vec<C64>, transformed: &Vec<C64>, &exponent: &i64, tolerance| {
            let scale = fs_math::det::powi(2.0, exponent as i32);
            let expected: Vec<C64> = base
                .iter()
                .map(|value| C64::new(value.re * scale, value.im * scale))
                .collect();
            minimum_complex_margin(tolerance, &expected, transformed)
        },
    );

    check_relation(
        "fs-fft::Fft::forward",
        0x2ACE_0002,
        384,
        |stream| {
            let exponent = match stream.int_in(0, 5) {
                0 => -3,
                1 => -2,
                2 => -1,
                3 => 1,
                4 => 2,
                _ => 3,
            };
            RelationCase::new(
                (0..N)
                    .map(|_| (stream.f64_in(-32.0, 32.0), stream.f64_in(-32.0, 32.0)))
                    .collect(),
                exponent,
            )
        },
        &operator,
        &relation,
    );
}

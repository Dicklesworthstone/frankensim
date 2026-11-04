//! Exact production interpolation/convolution leaf, also directly rustc-testable.
use super::*;

fn response(taps: &[f64], first: usize, omega: f64) -> [f64; 2] {
    taps.iter()
        .enumerate()
        .fold([0.0; 2], |[re, im], (i, value)| {
            let phase = -omega * (first + i) as f64;
            [re + value * phase.cos(), im + value * phase.sin()]
        })
}
fn product(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

#[test]
fn original_fractional_delay_complex_accuracy_survives_shared_extraction() {
    let w = window();
    for j in 0..=64 {
        let fraction = j as f64 / 65.0;
        let h = taps(fraction, &w);
        assert!((h.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        for bin in 0..=90 {
            let omega = TAU * 0.45 * bin as f64 / 90.0;
            let actual = response(&h, 0, omega);
            let phase = -omega * (HALF as f64 + fraction);
            assert!((actual[0] - phase.cos()).hypot(actual[1] - phase.sin()) < 4e-4);
        }
    }
}

#[test]
fn composition_retains_each_path_phase_and_signed_source_weight() {
    let sources = [vec![0.6, 0.3, 0.1], vec![0.2, 0.5, 0.3]];
    let delays = [41.25, 83.75];
    let weights = [1.4, -0.75];
    let mut budget = 10_000;
    let paths = [
        compose(&sources[0], delays[0], 1024, &mut budget).unwrap(),
        compose(&sources[1], delays[1], 1024, &mut budget).unwrap(),
    ];
    for bin in 0..=90 {
        let omega = TAU * 0.45 * bin as f64 / 90.0;
        let (mut actual, mut exact) = ([0.0; 2], [0.0; 2]);
        for i in 0..2 {
            let observed = response(&paths[i].1, paths[i].0, omega);
            let phase = -omega * (delays[i] + HALF as f64);
            let independent = product(response(&sources[i], 0, omega), [phase.cos(), phase.sin()]);
            for a in 0..2 {
                actual[a] += weights[i] * observed[a];
                exact[a] += weights[i] * independent[a];
            }
        }
        assert!((actual[0] - exact[0]).hypot(actual[1] - exact[1]) < 9e-4);
    }
    // The actual convolution theorem, independently of ideal sinc accuracy.
    let omega = 0.731;
    for i in 0..2 {
        let expected = product(
            response(&sources[i], 0, omega),
            response(
                &taps(delays[i].fract(), &window()),
                delays[i].floor() as usize,
                omega,
            ),
        );
        let actual = response(&paths[i].1, paths[i].0, omega);
        assert!((actual[0] - expected[0]).hypot(actual[1] - expected[1]) < 2e-14);
    }
}

#[test]
fn a_path_response_is_not_replaced_by_a_gain_or_dc_normalization() {
    let mut budget = 1000;
    let (first, h) = compose(&[0.2, 0.1], 10.0, 256, &mut budget).unwrap();
    assert_eq!(first, 10);
    assert!((h.iter().sum::<f64>() - 0.3).abs() < 1e-15);
    assert_eq!(h[HALF], 0.2);
    assert_eq!(h[HALF + 1], 0.1);
    assert!(h[..HALF].iter().all(|x| *x == 0.0));
    assert_eq!(budget, 1000 - 2 * TAPS);
}

#[test]
fn residual_travel_and_common_padding_count_filter_delay_once() {
    // Two finite causal material responses with different intrinsic delay.
    // The second needs negative residual travel; common padding preserves
    // their physical relative timing without ever reading a future sample.
    let laws = [
        vec![0.0, 0.0, 0.75, 0.25],
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.25],
    ];
    let physical = [5.25, 5.75];
    let own_delay = [2.25, 8.25];
    let padding = 3.0;
    let mut budget = 10_000;
    for i in 0..2 {
        let (first, h) = compose(
            &laws[i],
            physical[i] - own_delay[i] + padding,
            256,
            &mut budget,
        )
        .unwrap();
        let gain: f64 = h.iter().sum();
        let moment: f64 = h
            .iter()
            .enumerate()
            .map(|(j, x)| (first + j) as f64 * x)
            .sum();
        assert!((moment / gain - (physical[i] + padding + HALF as f64)).abs() < 1e-4);
        assert!(first < 256 && first + h.len() <= 256);
    }
}

#[test]
fn work_history_and_invalid_path_refusals_precede_publication() {
    for (h, delay, frames) in [
        (vec![], 1.0, 256),
        (vec![f64::NAN], 1.0, 256),
        (vec![1.0], -0.1, 256),
        (vec![1.0], f64::INFINITY, 256),
        (vec![1.0], 240.0, 256),
        (vec![1.0], 1.0, 0),
    ] {
        let mut budget = 1000;
        assert!(compose(&h, delay, frames, &mut budget).is_err());
        assert_eq!(budget, 1000);
    }
    let mut budget = TAPS - 1;
    assert!(compose(&[1.0], 1.0, 256, &mut budget).is_err());
    assert_eq!(budget, TAPS - 1);
}

#[test]
fn lossless_path_is_the_original_fractional_kernel_bit_for_bit() {
    for delay in [1.0_f64, 15.25, 64.5, 150.999] {
        let mut budget = 1000;
        let (first, h) = compose(&[1.0], delay, 1024, &mut budget).unwrap();
        assert_eq!(first, delay.floor() as usize);
        assert_eq!(h, taps(delay - delay.floor(), &window()));
    }
}

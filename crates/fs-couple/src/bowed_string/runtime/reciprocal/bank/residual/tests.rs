//! Residual flexibility against closed forms and against the exact
//! point-supported string (MagentaBadger 7592, gate 1). The exact roots come
//! from shooting the piecewise-sine solution here, with std trigonometry; they
//! share nothing with the residual owner. Predictions: sharp/pred-r3.md.

use super::*;

const SPAN_M: f64 = 0.383;
const A_TENSION_N: f64 = 55.486_491_050_429_706;
const DENSITY_KG_M: f64 = 1e-3;
const BRIDGE: f64 = 0.856_396_866_840_731_1;
const BRIDGE_N_M: f64 = 2e5;

fn card(bending: f64, modes: usize) -> BowedStringCard {
    BowedStringCard { length_m: SPAN_M, tension_n: A_TENSION_N, linear_density_kg_m: DENSITY_KG_M,
        bending_stiffness_n_m2: bending, viscous_bending_n_m2_s: 0.0, mode_count: modes,
        zetas: vec![0.0; modes], sample_rate_hz: 48_000 }
}

#[test]
fn static_green_is_the_complete_modal_sum() {
    // With EI the omitted tail falls as 1/N^3: 40 000 sines leave about 2e-10 of it.
    let beam = card(2e-5, 40_000);
    for (x, y) in [(0.1, 0.1), (0.1, BRIDGE), (0.5, 0.52), (BRIDGE, BRIDGE), (0.93, 0.07)] {
        let exact = static_green_m_n(&beam, x, y);
        let sum = retained_compliance_m_n(&beam, x, y);
        assert!((exact-sum).abs() <= 1e-8*exact.abs(), "x {x} y {y}: {exact:e} vs {sum:e}");
        assert_eq!(exact.to_bits(), static_green_m_n(&beam, y, x).to_bits());
        assert_eq!(sum.to_bits(), retained_compliance_m_n(&beam, y, x).to_bits());
    }
    // The taut string's own diagonal tail is 4/(pi^2 N) of it.
    let string = card(0.0, 1_000_000);
    let exact = static_green_m_n(&string, 0.3, 0.3);
    let textbook = 0.3*0.7*SPAN_M/A_TENSION_N;
    assert!((exact-textbook).abs() <= 4e-15*textbook, "{exact:e} vs {textbook:e}");
    let sum = retained_compliance_m_n(&string, 0.3, 0.3);
    assert!(exact-sum > 0.0 && exact-sum <= 1e-6*exact, "{exact:e} vs {sum:e}");
}

#[test]
fn residual_is_positive_symmetric_and_refuses_supports_off_the_span() {
    let a = card(1.6e-5, 25);
    let stations = [0.09, 0.2, BRIDGE];
    for &x in &stations {
        assert!(residual_compliance_m_n(&a, x, x).unwrap() > 0.0);
        for &y in &stations {
            assert_eq!(residual_compliance_m_n(&a, x, y).unwrap().to_bits(),
                residual_compliance_m_n(&a, y, x).unwrap().to_bits());
        }
    }
    for (x, y) in [(0.0, 0.5), (0.5, 1.0), (-0.1, 0.5), (0.5, f64::NAN)] {
        assert!(residual_compliance_m_n(&a, x, y).is_err());
    }
    let mut slack = a.clone();
    slack.tension_n = 0.0;
    assert!(residual_compliance_m_n(&slack, 0.5, 0.5).is_err());
    assert!(condense(-1e-9, BRIDGE_N_M, 0.0).is_err());
    // k G = 5: the spring and dashpot see 1/6 and 1/36 of themselves.
    let condensed = condense(2.5e-5, BRIDGE_N_M, 3.0).unwrap();
    assert!((condensed.stiffness_n_m*6.0/BRIDGE_N_M-1.0).abs() <= 1e-15);
    assert!((condensed.damping_n_s_m*36.0/3.0-1.0).abs() <= 1e-15);
}

// y(L) for wavenumber kw: y = sin(kw x) from the nut, each spring K to ground
// adding the slope jump K y / T. Springs are in ascending station order.
fn shoot(kw: f64, springs: &[(f64, f64)]) -> f64 {
    let (mut y, mut dy, mut from) = (0.0_f64, kw, 0.0);
    for &(station, k) in springs {
        let x = station*SPAN_M;
        let (s, c) = (kw*(x-from)).sin_cos();
        (y, dy) = (y*c + dy/kw*s, -y*kw*s + dy*c);
        dy += k*y/A_TENSION_N;
        from = x;
    }
    let (s, c) = (kw*(SPAN_M-from)).sin_cos();
    y*c + dy/kw*s
}

// The lowest `count` exact frequencies [Hz], by a fine wavenumber scan and bisection.
fn exact(springs: &[(f64, f64)], count: usize) -> Vec<f64> {
    let step = 1e-3;
    let mut roots = Vec::new();
    let (mut lo, mut f_lo) = (step, shoot(step, springs));
    while roots.len() < count {
        let hi = lo+step;
        let f_hi = shoot(hi, springs);
        if f_lo == 0.0 || f_lo.signum() != f_hi.signum() {
            let (mut a, mut b, f_a) = (lo, hi, f_lo);
            for _ in 0..200 {
                let mid = 0.5*(a+b);
                if mid <= a || mid >= b { break; }
                if shoot(mid, springs).signum() == f_a.signum() { a = mid; } else { b = mid; }
            }
            roots.push(0.5*(a+b));
        }
        (lo, f_lo) = (hi, f_hi);
        assert!(lo < 5e3, "scan found only {} roots", roots.len());
    }
    roots.iter().map(|kw| kw*(A_TENSION_N/DENSITY_KG_M).sqrt()/(2.0*core::f64::consts::PI)).collect()
}

// Mass-normalized sine shape of `card` at `station`.
fn shape(card: &BowedStringCard, station: f64) -> Vec<f64> {
    let scale = (2.0/(card.linear_density_kg_m*card.length_m)).sqrt();
    (1..=card.mode_count).map(|n| scale*(n as f64*core::f64::consts::PI*station).sin()).collect()
}

// Ascending frequencies [Hz] of the retained sines plus the supports, each a
// rank-one stiffness k c c^T.
fn retained(card: &BowedStringCard, supports: &[(f64, Vec<f64>)]) -> Vec<f64> {
    let n = card.mode_count;
    let mut stiffness = vec![0.0; n*n];
    for m in 0..n {
        let k = (m+1) as f64*core::f64::consts::PI/card.length_m;
        stiffness[m*n+m] = (card.tension_n*k*k + card.bending_stiffness_n_m2*k*k*k*k)/card.linear_density_kg_m;
    }
    for (k, column) in supports {
        for i in 0..n {
            for j in 0..n { stiffness[i*n+j] += k*column[i]*column[j]; }
        }
    }
    let (squared, _) = fs_la::eigen::jacobi_eigh(&stiffness, n);
    squared.iter().map(|w2| w2.sqrt()/(2.0*core::f64::consts::PI)).collect()
}

fn cents(f: f64, reference: f64) -> f64 { 1200.0*(f/reference).log2() }

#[test]
fn residual_flexibility_restores_the_exact_point_supported_string() {
    // The taut A string at N 25 with B4's or D5's stop (its summed K at that
    // station, to ground), alone and beside the 2e5 N/m bridge. Gate: f1..f4
    // within 1 c of exact and at most a quarter of the uncorrected error.
    let a = card(0.0, 25);
    for (note, semitones) in [("B4", 2.0_f64), ("D5", 5.0)] {
        let stop = BRIDGE*(1.0-(-semitones/12.0).exp2());
        for bridge_n_m in [0.0, BRIDGE_N_M] {
            for stop_n_m in [1e5, 1e6, 1e7] {
                let mut springs = vec![(stop, stop_n_m)];
                if bridge_n_m > 0.0 { springs.push((BRIDGE, bridge_n_m)); }
                let truth = exact(&springs, 10);
                let (phi_s, phi_b) = (shape(&a, stop), shape(&a, BRIDGE));
                let plain = retained(&a, &[(stop_n_m, phi_s.clone()), (bridge_n_m, phi_b.clone())]);

                let bridge = ResidualBridge { card: &a, station_fraction: BRIDGE, axis: [1.0, 0.0],
                    condensed: condense(residual_compliance_m_n(&a, BRIDGE, BRIDGE).unwrap(), bridge_n_m, 0.0).unwrap() };
                let (cross, compliance) = stop_residual(&[bridge],
                    &[ResidualStop { string: 0, station_fraction: stop, direction: [1.0, 0.0] }]).unwrap();
                let column: Vec<f64> = phi_s.iter().zip(&phi_b).map(|(s, b)| s-cross[0]*b).collect();
                let corrected = retained(&a, &[(stop_n_m/(1.0+stop_n_m*compliance[0]), column),
                    (bridge.condensed.stiffness_n_m, phi_b)]);

                let before: Vec<f64> = plain.iter().zip(&truth).map(|(f, e)| cents(*f, *e)).collect();
                let after: Vec<f64> = corrected.iter().zip(&truth).map(|(f, e)| cents(*f, *e)).collect();
                eprintln!("{note} bridge {bridge_n_m:e} K {stop_n_m:e}: uncorrected f1..4 {:+.3?}; corrected f1..4 \
                    {:+.3?}; corrected n5..10 (ungated) {:+.2?}", &before[..4], &after[..4], &after[4..10]);
                for (n, (after, before)) in after.iter().zip(&before).take(4).enumerate() {
                    assert!(after.abs() <= 1.0 && after.abs() <= 0.25*before.abs(),
                        "{note} bridge {bridge_n_m:e} K {stop_n_m:e} f{}: {after:+.3} c corrected, {before:+.3} c uncorrected",
                        n+1);
                }
            }
        }
    }
}

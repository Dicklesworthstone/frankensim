//! Cold-path parts of the live kernel against the reference facilities they
//! stand in for (integration test: the library is built WITHOUT `cfg(test)`,
//! because the reference modules' own unit tests need the full example tree).
use fs_live_piano::dense::{eigh_gen, exp_neg_small, spd_inverse};
use fs_live_piano::fastfelt::Site;
use fs_live_piano::refmods::{felt, geometry, steinway_scale};
use fs_material::Uniaxial;

// Actual repeated-contact operands, with independent 100-digit Decimal answers.
// The live adapter must not depend on the reference implementation to detect
// a shared arithmetic mistake. These do not establish musical fidelity.
type CapturedFeltCase = (f64, f64, f64, f64);

fn captured_felt_cases() -> (fs_material::WoolFelt, felt::State, [CapturedFeltCase; 5]) {
    let law = fs_material::WoolFelt::new(10784971.283384772, 0.2, 2.5116355500000003,
        3.2116355500000004, 0.25, 0.8).unwrap();
    let state = felt::State { eps_max: 0.05313385293183455, sig_max: 386343.9999062737 };
    let a = 0.0003418749043622657;
    let b2 = 0.00034187507285155427;
    let b32 = 0.00034187507285085426;
    (law, state, [
        (a, b2, 6.911_181_590_113_403, 45597.24015008384),
        (b2, a, 6.911_181_590_113_403, 45597.21688026879),
        (a, b32, 6.911181590081485, 45597.24014989049),
        (b32, a, 6.911181590081485, 45597.21688017211),
        (a, a, 6.91117390747077, 45597.19361045862),
    ])
}

#[test]
fn captured_felt_quadrature_force_matches_decimal() {
    let (law, state, cases) = captured_felt_cases();
    for (start, end, expected, _) in cases {
        let (force, _) = Site::new(&law, &state, start,
            0.007413793103448275, 3.3333333333333335e-5).force(end);
        assert!((force - expected).abs() < 1e-12,
            "{start:.17e}->{end:.17e}: force={force:.17e}, reference={expected:.17e}");
    }
}

#[test]
fn captured_felt_quadrature_derivative_matches_decimal() {
    let (law, state, cases) = captured_felt_cases();
    for (start, end, _, expected) in cases {
        let (_, slope) = Site::new(&law, &state, start,
            0.007413793103448275, 3.3333333333333335e-5).force(end);
        assert!((slope - expected).abs() < 1e-8,
            "{start:.17e}->{end:.17e}: slope={slope:.17e}, reference={expected:.17e}");
    }
}

#[test]
fn prepared_felt_cubic_boundaries_match_analytic_moments() {
    for crushed in [false, true] {
        let law = fs_material::WoolFelt::new(4.0, 1.0, if crushed { 2.0 } else { 3.0 },
            3.0, if crushed { 0.25 } else { 0.0 }, 2.0).unwrap();
        let state = if crushed { felt::State { eps_max: 1.0, sig_max: 4.0 } }
            else { law.initial_state() };
        let origin = if crushed { 0.25 } else { 0.0 };
        let coefficient = if crushed { 4.0 / (0.75 * 0.75 * 0.75) } else { 4.0 };
        let a = 0.5_f64;
        let series_edge = |t: f64| origin + (a - origin) * (1.0 + t) / (1.0 - t);
        for b in [a, a.next_up(), a * (1.0 + 0.999e-7), a * (1.0 + 1.001e-7),
            series_edge((1.0 - 1e-8) / 24.0), series_edge((1.0 + 1e-8) / 24.0)] {
            for (start, end) in [(a, b), (b, a)] {
                let (x, y) = (start - origin, end - origin);
                let force = coefficient * (x*x*x + x*x*y + x*y*y + y*y*y) / 4.0;
                let slope = coefficient * (x*x + 2.0*x*y + 3.0*y*y) / 4.0;
                let (actual_force, actual_slope) = Site::new(&law, &state, start, 1.0, 1.0).force(end);
                assert!((actual_force - force).abs() < 1e-12,
                    "crushed={crushed} {start:.17e}->{end:.17e}: force {actual_force} vs {force}");
                assert!((actual_slope - slope).abs() < 1e-11,
                    "crushed={crushed} {start:.17e}->{end:.17e}: slope {actual_slope} vs {slope}");
            }
        }
    }
}

#[test]
fn captured_felt_quadrature_committed_maximum_matches_cubic() {
    let law = fs_material::WoolFelt::new(4.0, 1.0, 2.0, 3.0, 0.25, 2.0).unwrap();
    let state = felt::State { eps_max: 1.0, sig_max: 4.0 };
    let (start, end) = (1.0, 0.9999998);
    let (x, y) = (start - 0.25, end - 0.25);
    let coefficient = 4.0 / (0.75 * 0.75 * 0.75);
    let force = coefficient * (x*x*x + x*x*y + x*y*y + y*y*y) / 4.0;
    let slope = coefficient * (x*x + 2.0*x*y + 3.0*y*y) / 4.0;
    // Starting exactly at the maximum and unloading never crosses branches.
    // Keep both answers in the failure output so neither hides the other.
    let actual = [
        Site::new(&law, &state, start, 1.0, 1.0).force(end),
        felt::average(&law, &state, start, end, 1.0, 1.0),
    ];
    assert!(actual.iter().all(|(f, df)| (f - force).abs() < 1e-12 && (df - slope).abs() < 1e-11),
        "analytic ({force:.17e}, {slope:.17e}); prepared/reference {actual:?}");
}

#[test]
fn prepared_felt_tiny_positive_width_retains_tangent() {
    let law = fs_material::WoolFelt::new(4.0, 1.0, 2.0, 2.0, 0.0, 2.0).unwrap();
    let a = 1e-200_f64;
    let b = a.next_up();
    let (_, slope) = Site::new(&law, &law.initial_state(), a, 1.0, 1.0).force(b);
    let expected = 4.0 * (a + 2.0*b) / 3.0;
    assert!(slope.is_finite() && slope > 0.0, "lost positive tangent: {slope:e}");
    assert!((slope / expected - 1.0).abs() < 2e-14,
        "tiny-interval slope {slope:e}, analytic reference {expected:e}");
}


fn pencil(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut k = vec![0.0; n * n];
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let g = ((i * 7 + j * 3) % 11) as f64 * 0.01;
            k[i * n + j] += g * 1.0e4;
            k[j * n + i] += g * 1.0e4;
            m[i * n + j] += g * 0.005;
            m[j * n + i] += g * 0.005;
        }
        k[i * n + i] += 1.0e6 * (i + 1) as f64;
        m[i * n + i] += 1.0;
    }
    (k, m)
}

#[test]
fn generalized_modes_match_the_reference_eigen_facility() {
    for n in [1usize, 4, 25, 43] {
        let (k, m) = pencil(n);
        let (lambda, phi) = eigh_gen(&k, &m, n).unwrap();
        let reference = fs_modal::eigh_gen_dense(&k, &m, n).unwrap();
        for j in 0..n {
            assert!((lambda[j] / reference[j].lambda - 1.0).abs() < 1.0e-12, "mode {j}");
            // K phi = lambda M phi and Phi^T M Phi = I.
            for i in 0..n {
                let kp: f64 = (0..n).map(|a| k[i * n + a] * phi[a * n + j]).sum();
                let mp: f64 = (0..n).map(|a| m[i * n + a] * phi[a * n + j]).sum();
                assert!((kp - lambda[j] * mp).abs() < 1.0e-8 * lambda[j].abs());
            }
            for other in 0..n {
                let dotm: f64 = (0..n).map(|a| phi[a * n + j]
                    * (0..n).map(|b| m[a * n + b] * phi[b * n + other]).sum::<f64>()).sum();
                assert!((dotm - if other == j { 1.0 } else { 0.0 }).abs() < 1.0e-12);
            }
        }
    }
}

#[test]
fn spd_inverse_and_series_exponential_are_what_they_claim() {
    let n = 25;
    let (_, m) = pencil(n);
    let inverse = spd_inverse(&m, n, "test").unwrap();
    for i in 0..n {
        for j in 0..n {
            let value: f64 = (0..n).map(|a| m[i * n + a] * inverse[a * n + j]).sum();
            assert!((value - if i == j { 1.0 } else { 0.0 }).abs() < 1.0e-12);
        }
    }
    let x: Vec<f64> = m.iter().map(|v| 1.0e-3 * v).collect();
    let e = exp_neg_small(&x, n).unwrap();
    // Reference: spectral exponential through fs-modal.
    let eye: Vec<f64> = (0..n * n).map(|i| if i / n == i % n { 1.0 } else { 0.0 }).collect();
    let modes = fs_modal::eigh_gen_dense(&x, &eye, n).unwrap();
    for i in 0..n {
        for j in 0..n {
            let value: f64 = modes.iter().map(|p| (-p.lambda).exp() * p.phi[i] * p.phi[j]).sum();
            assert!((value - e[i * n + j]).abs() < 1.0e-14);
        }
    }
    assert!(spd_inverse(&[0.0; 4], 2, "zero").is_err());
    assert!(exp_neg_small(&[1.0; 4], 2).is_err());
}


#[test]
fn prepared_force_and_solve_reproduce_the_reference_felt_adapter() {
    let scale = geometry::read_scale(fs_live_piano::kernel::SCALE_CSV).unwrap();
    let mut worst_force = 0.0_f64;
    let mut worst_solve = 0.0_f64;
    for &midi in &[21u8, 40, 60, 84, 108] {
        let c = scale[usize::from(midi - 21)];
        let (law, _) = steinway_scale::hammer_material(&c).unwrap();
        let (t, area) = (c.felt_thickness_m, c.felt_area_m2 / c.unison as f64);
        let mut state = law.initial_state();
        for commit in [0.0, 0.05, 0.18, 0.4] {
            state = law.update_state(commit, &state);
            for &start in &[-0.02 * t, 0.0, 0.01 * t, 0.1 * t, 0.3 * t] {
                let site = Site::new(&law, &state, start, t, area);
                for i in 0..60 {
                    let end = (-0.05 + 0.012 * i as f64) * t;
                    let (f, df) = site.force(end);
                    let (rf, rdf) = felt::average(&law, &state, start, end, t, area);
                    worst_force = worst_force.max((f - rf).abs() / (1.0e-9 + rf.abs()));
                    assert!((df - rdf).abs() <= 1.0e-8 * (1.0 + rdf.abs()), "tangent {df} {rdf}");
                }
                for &free in &[-1.0e-4, 0.0, 2.0e-5, 2.0e-4, 8.0e-4] {
                    for &compliance in &[2.0e-8, 3.0e-6] {
                        let reference = felt::solve_with_rate(&law, &state, start, free, compliance, t, area, 0.0, 1.0);
                        let fast = site.solve(free, compliance, 0.0);
                        match (reference, fast) {
                            (Ok(a), Ok(b)) => {
                                // Both satisfy |end + C F - free| <= 1e-13 m.
                                worst_solve = worst_solve.max((a - b).abs() * compliance);
                                let warm = site.solve(free, compliance, a * 1.01).unwrap();
                                worst_solve = worst_solve.max((a - warm).abs() * compliance);
                            }
                            (Err(_), Err(_)) => {}
                            other => panic!("solver admission differs: {other:?}"),
                        }
                    }
                }
            }
        }
    }
    assert!(worst_force < 1.0e-12, "relative force deviation {worst_force:e}");
    assert!(worst_solve <= 2.0e-13, "overlap deviation {worst_solve:e} m");
}

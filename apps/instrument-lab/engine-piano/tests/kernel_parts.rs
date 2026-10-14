//! Cold-path parts of the live kernel against the reference facilities they
//! stand in for (integration test: the library is built WITHOUT `cfg(test)`,
//! because the reference modules' own unit tests need the full example tree).
use fs_live_piano::dense::{eigh_gen, exp_neg_small, spd_inverse};
use fs_live_piano::fastfelt::Site;
use fs_live_piano::refmods::{felt, geometry, steinway_scale};
use fs_material::Uniaxial;


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
    assert!(spd_inverse(&vec![0.0; 4], 2, "zero").is_err());
    assert!(exp_neg_small(&vec![1.0; 4], 2).is_err());
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

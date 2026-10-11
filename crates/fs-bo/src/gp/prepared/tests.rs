use super::*;
use crate::gp::{Kernel, Matern};

fn fixture(family: Matern) -> (Gp, Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
    let mut x: Vec<_> = (0..12).map(|i| vec![i as f64 / 11.0, (i * i % 7) as f64 / 7.0]).collect();
    x[3] = x[2].clone(); // noisy replicate, not a second latent variable
    let y: Vec<_> = x.iter().map(|p| p[0] * p[0] - 0.7 * p[1]).collect();
    let noise: Vec<_> = (0..12).map(|i| 0.01 + 0.003 * i as f64).collect();
    let kernel = Kernel { family, signal: 1.3, lengthscales: vec![0.37, 0.81] };
    (Gp::try_fit_diag(&x, &y, kernel, &noise).unwrap(), x, y, noise)
}

// Independent dense reference uses full C^-1 solves, not the production
// L^-1 projection/dot route, and factors the ENTIRE joined covariance.
fn dense(gp: &Gp, train: &[Vec<f64>], y: &[f64], noise: &[f64],
    points: &[Vec<f64>], jitter: f64) -> (Vec<f64>, Vec<f64>, Cholesky)
{
    let n = train.len(); let m = points.len();
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i { c[i * n + j] = gp.kernel.eval(&train[i], &train[j]); }
        c[i * n + i] += noise[i];
    }
    let c = cholesky(&c, n).unwrap();
    let mut alpha = y.to_vec(); c.solve(&mut alpha);
    let cross: Vec<Vec<f64>> = points.iter().map(|p| train.iter().map(|t| gp.kernel.eval(p, t)).collect()).collect();
    let solved: Vec<_> = cross.iter().map(|row| { let mut r = row.clone(); c.solve(&mut r); r }).collect();
    let mean = cross.iter().map(|r| r.iter().zip(&alpha).map(|(a, b)| a * b).sum()).collect();
    let mut covariance = vec![0.0; m * m];
    for i in 0..m {
        for j in 0..=i {
            let v = gp.kernel.eval(&points[i], &points[j])
                - cross[i].iter().zip(&solved[j]).map(|(a, b)| a * b).sum::<f64>();
            let v = if i == j { v.max(0.0) + jitter } else { v };
            covariance[i * m + j] = v; covariance[j * m + i] = v;
        }
    }
    let factor = cholesky(&covariance, m).unwrap();
    (mean, covariance, factor)
}

#[test]
fn g1_extended_covariances_match_independent_dense_solves_across_factor_panels() {
    for family in [Matern::Half, Matern::ThreeHalves, Matern::FiveHalves] {
        let (gp, x, y, noise) = fixture(family);
        for b in [1, 7, 33, 65] {
            let baseline: Vec<_> = (0..b).map(|i| vec![(i + 1) as f64 / (b + 1) as f64, 0.23]).collect();
            let candidates = vec![vec![0.13, 0.67], vec![0.81, 0.37], vec![1.2, -0.4]];
            let prepared = PreparedPosterior::new(&gp, &baseline, 1e-6).unwrap();
            let extension = prepared.extend(&candidates).unwrap();
            let mut points = baseline.clone(); points.extend(candidates);
            let m = points.len();
            let (mean, covariance, factor) = dense(&gp, &x, &y, &noise, &points, 1e-6);
            assert_eq!(extension.lower_rows.len(), 3 * m, "candidate memory only, no baseline copy");
            assert_eq!(prepared.diagonal_jitter(), 1e-6);
            let l = |i: usize, j: usize| if i < b { prepared.factor().l(i, j) }
                else { extension.lower_rows[(i - b) * m + j] };
            for i in 0..m {
                let got_mean = if i < b { prepared.mean()[i] } else { extension.mean[i - b] };
                assert!((got_mean - mean[i]).abs() < 2e-12);
                for j in 0..=i {
                    let rebuilt: f64 = (0..=j).map(|k| l(i, k) * l(j, k)).sum();
                    assert!((rebuilt - covariance[i * m + j]).abs() < 2e-11,
                        "{family:?} b={b} at {i},{j}: {rebuilt} != {}", covariance[i * m + j]);
                    assert!((l(i, j) - factor.l(i, j)).abs() < 2e-8);
                }
            }
        }
    }
}

#[test]
fn g3_prefix_means_and_draws_never_change_across_candidate_probes() {
    let (gp, _, _, _) = fixture(Matern::FiveHalves);
    let baseline = vec![vec![0.1, 0.2], vec![0.8, 0.9]];
    let prepared = PreparedPosterior::new(&gp, &baseline, 1e-8).unwrap();
    let original = (prepared.mean().to_vec(), prepared.factor().l(1, 0), prepared.factor().l(1, 1));
    let candidates = vec![vec![0.10001, 0.2], vec![0.79, 0.88]];
    let a = prepared.extend(&candidates).unwrap();
    let b = prepared.extend(&[vec![-5.0, 4.0]]).unwrap();
    assert!(b.mean[0].is_finite());
    assert_eq!(a, prepared.extend(&candidates).unwrap());
    assert_eq!(original, (prepared.mean().to_vec(), prepared.factor().l(1, 0), prepared.factor().l(1, 1)));
}

#[test]
fn g0_duplicate_and_malformed_points_and_jitter_refuse() {
    let (gp, _, _, _) = fixture(Matern::Half);
    for jitter in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(matches!(PreparedPosterior::new(&gp, &[vec![0.0, 0.0]], jitter), Err(PreparedError::InvalidInput)));
    }
    assert!(matches!(PreparedPosterior::new(&gp, &[vec![0.0, 0.0], vec![-0.0, 0.0]], 1e-8),
        Err(PreparedError::DuplicatePoint)));
    let p = PreparedPosterior::new(&gp, &[vec![0.0, 0.0]], 1e-8).unwrap();
    for bad in [vec![], vec![vec![0.1]], vec![vec![0.1, f64::NAN]]] {
        assert_eq!(p.extend(&bad), Err(PreparedError::InvalidInput));
    }
    assert_eq!(p.extend(&[vec![-0.0, 0.0]]), Err(PreparedError::DuplicatePoint));
    assert_eq!(p.extend(&[vec![0.3, 0.4], vec![0.3, 0.4]]), Err(PreparedError::DuplicatePoint));
    assert_eq!(entries(usize::MAX, 2), Err(PreparedError::ShapeOverflow));
}

#[test]
fn g4_cancelled_extensions_retain_the_reusable_baseline() {
    let (gp, _, _, _) = fixture(Matern::ThreeHalves);
    let base = vec![vec![0.0, 0.2], vec![0.7, 0.9]];
    assert!(matches!(PreparedPosterior::new_controlled(&gp, &base, 1e-8, &mut || false),
        Err(PreparedError::Cancelled)));
    let p = PreparedPosterior::new(&gp, &base, 1e-8).unwrap();
    let candidates = vec![vec![0.3, 0.1], vec![0.8, 0.7]];
    let expected = p.extend(&candidates).unwrap();
    let mut checks = 0;
    assert_eq!(p.extend_controlled(&candidates, &mut || { checks += 1; checks < 4 }),
        Err(PreparedError::Cancelled));
    assert_eq!(checks, 4); // second candidate factor row
    assert_eq!(p.extend(&candidates).unwrap(), expected);
}

#[test]
fn g0_singular_nominal_posterior_requires_explicit_regularization() {
    let kernel = Kernel { family: Matern::Half, signal: 1.0, lengthscales: vec![1.0] };
    let gp = Gp::fit(&[vec![0.0]], &[0.0], kernel, 0.0);
    assert!(matches!(PreparedPosterior::new(&gp, &[vec![0.0]], 0.0),
        Err(PreparedError::NotPositiveDefinite)));
    let p = PreparedPosterior::new(&gp, &[vec![0.0]], 1e-8).unwrap();
    assert_eq!(p.diagonal_jitter(), 1e-8);
    assert_eq!(gp.noise, 0.0, "no new observation noise");
    assert!(p.extend(&[vec![1e-10]]).is_ok());
}

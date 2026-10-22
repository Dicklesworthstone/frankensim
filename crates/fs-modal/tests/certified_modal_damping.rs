//! IRP M0.2-W5: certified modal damping against dense truth.
//!
//! The truth oracle is independent of the certifier. The dense truth comes
//! from `fs_modal::quadratic_eigenvalues`, refined per pole by inverse
//! iteration and the two-sided Rayleigh functional `xᵀQ(λ)x = 0`. Q is
//! complex symmetric, so x is also the left eigenvector and the refinement
//! converges quadratically. The raw dense values can sit ~1e-12 relative
//! off. That is comparable to tight certified radii, so the test reports the
//! correction it applied. The refined values carry an allowance of
//! `TRUTH_ULPS·ε·max|λ|`.

use core::f64::consts::PI;
use fs_math::c64::C64;
use fs_math::det;
use fs_mms::{RefinementLadder, fit_order};
use fs_modal::{
    CertifiedModalDamping, MAX_DAMPING_CLUSTER, ModalDampingBlock, ModalDampingBudget, ModalError,
    certify_modal_damping, certify_modal_damping_windowed, certify_modal_damping_windowed_far_rows,
    damping_window_neighbours, quadratic_eigenvalues,
};
use fs_rand::StreamKey;

const TRUTH_ULPS: f64 = 8.0;

/// So loose that nothing merges for budget reasons.
fn loose() -> ModalDampingBudget {
    ModalDampingBudget {
        max_rel_decay_error: 1e9,
        max_freq_error_cents: 1e9,
        max_cluster_size: MAX_DAMPING_CLUSTER,
    }
}

fn budget() -> ModalDampingBudget {
    ModalDampingBudget {
        max_rel_decay_error: 1e-3,
        max_freq_error_cents: 0.01,
        max_cluster_size: MAX_DAMPING_CLUSTER,
    }
}

/// Solve `A x = b` by complex Gaussian elimination with partial pivoting.
fn solve(mut a: Vec<C64>, mut b: Vec<C64>, n: usize) -> Vec<C64> {
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&i, &j| a[i * n + col].abs().total_cmp(&a[j * n + col].abs()))
            .unwrap();
        if pivot != col {
            for k in 0..n {
                a.swap(col * n + k, pivot * n + k);
            }
            b.swap(col, pivot);
        }
        if a[col * n + col].abs() == 0.0 {
            // Exactly singular at an exact eigenvalue: any tiny pivot keeps
            // inverse iteration going.
            a[col * n + col] = C64::from_re(f64::MIN_POSITIVE.sqrt());
        }
        let p = a[col * n + col];
        for i in col + 1..n {
            let f = a[i * n + col] / p;
            for k in col..n {
                a[i * n + k] = a[i * n + k] - f * a[col * n + k];
            }
            b[i] = b[i] - f * b[col];
        }
    }
    let mut x = vec![C64::ZERO; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        for k in i + 1..n {
            sum = sum - a[i * n + k] * x[k];
        }
        x[i] = sum / a[i * n + i];
    }
    x
}

/// One pole refined from `seed`: inverse iteration on `Q(λ)` and the root of
/// `xᵀQ(λ)x` nearest the current value.
fn refine(omega: &[f64], c: &[f64], seed: C64) -> C64 {
    let n = omega.len();
    let mut lambda = seed;
    let mut x = vec![C64::ONE; n];
    for _ in 0..8 {
        let mut q = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                q[n * i + j] = lambda * C64::from_re(c[n * i + j]);
            }
            q[n * i + i] = q[n * i + i] + lambda * lambda + C64::from_re(omega[i] * omega[i]);
        }
        x = solve(q, x, n);
        let scale = x.iter().map(|v| v.abs()).fold(0.0, f64::max);
        for v in &mut x {
            *v = v.scale(1.0 / scale);
        }
        let (mut a2, mut a1, mut a0) = (C64::ZERO, C64::ZERO, C64::ZERO);
        for i in 0..n {
            a2 = a2 + x[i] * x[i];
            a0 = a0 + x[i] * x[i] * C64::from_re(omega[i] * omega[i]);
            for j in 0..n {
                a1 = a1 + x[i] * C64::from_re(c[n * i + j]) * x[j];
            }
        }
        let root = (a1 * a1 - C64::from_re(4.0) * a2 * a0).sqrt();
        let two_a2 = C64::from_re(2.0) * a2;
        let (r1, r2) = (
            (C64::ZERO - a1 + root) / two_a2,
            (C64::ZERO - a1 - root) / two_a2,
        );
        let next = if (r1 - lambda).abs() <= (r2 - lambda).abs() {
            r1
        } else {
            r2
        };
        let step = (next - lambda).abs();
        lambda = next;
        if !(step > 4.0 * f64::EPSILON * lambda.abs()) {
            break;
        }
    }
    lambda
}

/// The upper-half-plane poles (refined), the oracle allowance, and the
/// largest correction the refinement applied to the raw dense values.
fn truth(omega: &[f64], c: &[f64]) -> (Vec<C64>, f64, f64) {
    let n = omega.len();
    let mut identity = vec![0.0; n * n];
    let mut stiffness = vec![0.0; n * n];
    for i in 0..n {
        identity[n * i + i] = 1.0;
        stiffness[n * i + i] = omega[i] * omega[i];
    }
    let raw = quadratic_eigenvalues(&identity, c, &stiffness, n).expect("dense truth");
    let mut poles = Vec::with_capacity(n);
    let mut correction = 0.0_f64;
    for seed in raw.into_iter().filter(|p| p.im > 0.0) {
        let refined = refine(omega, c, seed);
        // A refinement that wanders off belongs to a neighbour; keep the seed.
        if (refined - seed).abs() <= 1e-6 * seed.abs() {
            correction = correction.max((refined - seed).abs() / seed.abs());
            poles.push(refined);
        } else {
            poles.push(seed);
        }
    }
    let largest = poles.iter().map(|p| p.abs()).fold(0.0, f64::max);
    (poles, TRUTH_ULPS * f64::EPSILON * largest, correction)
}

/// Every certified pole with its radius, in block order.
fn certified_poles(result: &CertifiedModalDamping) -> Vec<(C64, f64)> {
    let mut out = Vec::new();
    for block in &result.blocks {
        match block {
            ModalDampingBlock::Diagonal {
                eigenvalue, radius, ..
            } => out.push((*eigenvalue, *radius)),
            ModalDampingBlock::Full {
                eigenvalues, radii, ..
            } => out.extend(eigenvalues.iter().copied().zip(radii.iter().copied())),
        }
    }
    out
}

fn nearest(poles: &[C64], z: C64) -> f64 {
    poles
        .iter()
        .map(|&p| (p - z).abs())
        .fold(f64::INFINITY, f64::min)
}

/// `C = S (κG + (1 − κ)I) S`: `G` is a random Gram matrix normalized to a
/// unit diagonal (so |G_ij| ≤ 1; κ = 1 is strongly nonproportional) and
/// `S = diag(√(2ζω))`, giving `C_ii = 2ζ_iω_i` up to rounding. A convex
/// combination of PSD matrices, so PSD by construction.
fn random_damping(
    omega: &[f64],
    zeta: &[f64],
    rank: usize,
    coupling: f64,
    s: &mut fs_rand::Stream,
) -> Vec<f64> {
    let n = omega.len();
    let b: Vec<f64> = (0..n * rank).map(|_| 2.0 * s.next_f64() - 1.0).collect();
    let mut g = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0_f64;
            for k in 0..rank {
                sum = b[i * rank + k].mul_add(b[j * rank + k], sum);
            }
            g[n * i + j] = sum;
        }
    }
    let diag: Vec<f64> = (0..n).map(|i| det::sqrt(g[n * i + i])).collect();
    let scale: Vec<f64> = (0..n)
        .map(|i| det::sqrt(2.0 * zeta[i] * omega[i]))
        .collect();
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let unit = if i == j {
                1.0
            } else {
                coupling * g[n * i + j] / (diag[i] * diag[j])
            };
            c[n * i + j] = scale[i] * unit * scale[j];
        }
    }
    // Exact symmetry: rounding of the two products may differ.
    for i in 0..n {
        for j in 0..i {
            c[n * i + j] = c[n * j + i];
        }
    }
    c
}

fn log_uniform(s: &mut fs_rand::Stream, lo: f64, hi: f64) -> f64 {
    det::exp(det::ln(lo) + s.next_f64() * (det::ln(hi) - det::ln(lo)))
}

#[test]
fn g0_certified_radius_encloses_dense_truth() {
    let started = std::time::Instant::now();
    let (mut violations, mut checked) = (0usize, 0usize);
    let mut ratios = Vec::new();
    let (mut diagonal, mut full, mut largest_block) = (0usize, 0usize, 0usize);
    let mut partial_blocks = 0usize;
    let mut max_correction = 0.0_f64;
    let mut worst = (0.0_f64, 0usize);
    for trial in 0..4096_u32 {
        let mut s = StreamKey {
            seed: 0x4d30_3257_3561,
            kernel: 5,
            tile: trial,
        }
        .stream();
        let n = 4 + (s.next_below(45) as usize);
        let mut omega: Vec<f64> = (0..n)
            .map(|_| log_uniform(&mut s, 2.0 * PI * 50.0, 2.0 * PI * 5000.0))
            .collect();
        if s.next_f64() < 0.25 {
            // A planted near-degenerate pair.
            omega[1] = omega[0] * (1.0 + 1e-4 * s.next_f64());
        }
        let zeta: Vec<f64> = (0..n).map(|_| log_uniform(&mut s, 1e-4, 5e-2)).collect();
        let rank = 1 + s.next_below(n as u64) as usize;
        // From nearly proportional to fully nonproportional, so both the
        // singleton and the block paths carry weight.
        let kappa = log_uniform(&mut s, 1e-4, 1.0);
        let c = random_damping(&omega, &zeta, rank, kappa, &mut s);
        let result = certify_modal_damping(&omega, &c, &budget()).expect("certified");
        let (poles, slack, correction) = truth(&omega, &c);
        max_correction = max_correction.max(correction);
        let certified = certified_poles(&result);
        assert_eq!(certified.len(), n, "one upper pole per mode");
        for &(mu, r) in &certified {
            checked += 1;
            let d = nearest(&poles, mu);
            if d > r + slack {
                violations += 1;
            }
            if r > 100.0 * slack {
                // The oracle's own allowance is not charged to the radius.
                let ratio = (d - slack).max(0.0) / r;
                if ratio > worst.0 {
                    worst = (ratio, trial as usize);
                }
                ratios.push(ratio);
            }
        }
        for block in &result.blocks {
            match block {
                ModalDampingBlock::Diagonal { .. } => diagonal += 1,
                ModalDampingBlock::Full { modes, .. } => {
                    full += 1;
                    largest_block = largest_block.max(modes.len());
                    if modes.len() < n {
                        partial_blocks += 1;
                    }
                }
            }
        }
        if trial < 8 {
            assert_eq!(
                certify_modal_damping(&omega, &c, &budget()).unwrap(),
                result,
                "bitwise replay"
            );
        }
    }
    ratios.sort_by(f64::total_cmp);
    let median = ratios[ratios.len() / 2];
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g0_certified_radius_encloses_dense_truth\",\"tier\":\"G0\",\"inputs\":{{\"trials\":4096,\"n\":[4,48],\"zeta\":[1e-4,5e-2],\"planted_pair_fraction\":0.25}},\"provenance\":\"fs-rand seed 0x4d3032573561 kernel 5\",\"observed\":{{\"poles_checked\":{checked},\"violations\":{violations},\"ratio_samples\":{},\"max_enclosure_ratio\":{:e},\"worst_trial\":{},\"median_enclosure_ratio\":{median:e},\"diagonal_blocks\":{diagonal},\"full_blocks\":{full},\"full_blocks_with_complement\":{partial_blocks},\"largest_block\":{largest_block}}},\"tolerance\":{{\"truth_ulps\":{TRUTH_ULPS},\"raw_dense_truth_max_rel_correction\":{max_correction:e}}},\"verdict\":\"{}\",\"elapsed_s\":{:.1}}}",
        ratios.len(),
        worst.0,
        worst.1,
        if violations == 0 && worst.0 <= 1.0 {
            "pass"
        } else {
            "fail"
        },
        started.elapsed().as_secs_f64()
    );
    assert_eq!(violations, 0, "certified radii must enclose dense truth");
    // Coverage, not a quality claim: both certificate paths must be live.
    assert!(
        diagonal >= checked / 5,
        "diagonal path barely exercised: {diagonal}"
    );
    assert!(
        partial_blocks >= 100,
        "block-vs-complement path barely exercised"
    );
    assert!(worst.0 <= 1.0, "enclosure ratio {} > 1", worst.0);
}

#[test]
fn g0_rouche_radius_is_second_order() {
    // Well-separated modes, all singletons: the radius scales as ‖F‖².
    let omega = [100.0, 173.0, 251.0, 337.0, 419.0, 503.0];
    let n = omega.len();
    let zeta = [2e-3, 3e-3, 1e-3, 4e-3, 2.5e-3, 1.5e-3];
    let mut s = StreamKey {
        seed: 0x4d30_3257_3561,
        kernel: 6,
        tile: 0,
    }
    .stream();
    let base = random_damping(&omega, &zeta, n, 1.0, &mut s);
    let scales = [1.0, 0.5, 0.25, 0.125];
    let mut radii = Vec::new();
    for &t in &scales {
        let c: Vec<f64> = (0..n * n)
            .map(|k| if k / n == k % n { base[k] } else { t * base[k] })
            .collect();
        // The order of the radius, not the budget, is under test.
        let result = certify_modal_damping(&omega, &c, &loose()).unwrap();
        let (poles, slack, _) = truth(&omega, &c);
        let mut rung = Vec::new();
        for block in &result.blocks {
            let ModalDampingBlock::Diagonal {
                eigenvalue, radius, ..
            } = block
            else {
                panic!("t = {t}: expected all singletons, got {block:?}");
            };
            assert!(nearest(&poles, *eigenvalue) <= radius + slack);
            rung.push(*radius);
        }
        radii.push(rung);
    }
    let mut orders = Vec::new();
    for mode in 0..n {
        let errors: Vec<f64> = radii.iter().map(|r| r[mode]).collect();
        let fit = fit_order(&RefinementLadder::new(scales.to_vec(), errors).unwrap());
        orders.push(fit.observed);
    }
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g0_rouche_radius_is_second_order\",\"tier\":\"G0\",\"rungs\":{{\"f_scale\":{scales:?}}},\"observed\":{{\"order_per_mode\":{orders:?}}},\"tolerance\":{{\"order\":2.0,\"band\":0.2}}}}"
    );
    for (mode, order) in orders.iter().enumerate() {
        assert!((order - 2.0).abs() <= 0.2, "mode {mode}: order {order}");
    }
}

#[test]
fn g1_clustered_pair_triggers_full_block_fallback() {
    // Two modes split by 0.3 %, coupled through damping as strongly as PSD
    // allows; a third, well-separated mode stays diagonal.
    let omega = [1000.0, 1003.0, 1700.0];
    let (d0, d1, d2) = (4.0, 0.1, 3.0);
    let c12 = 0.6;
    let c = vec![
        d0, c12, 0.05, //
        c12, d1, 0.02, //
        0.05, 0.02, d2,
    ];
    let result = certify_modal_damping(&omega, &c, &budget()).unwrap();
    let (poles, slack, _) = truth(&omega, &c);
    let mut pair = None;
    for block in &result.blocks {
        match block {
            ModalDampingBlock::Full {
                modes,
                eigenvalues,
                radii,
                ..
            } => {
                assert_eq!(modes, &vec![0, 1]);
                for (&mu, &r) in eigenvalues.iter().zip(radii) {
                    let d = nearest(&poles, mu);
                    assert!(d <= r + slack, "block pole off by {d} > {r}");
                    assert!(
                        r <= 1e-3 * mu.re.abs(),
                        "block radius {r} breaks the decay budget"
                    );
                }
                pair = Some(eigenvalues.clone());
            }
            ModalDampingBlock::Diagonal { mode, radius, .. } => {
                assert_eq!(*mode, 2);
                assert!(*radius <= 1e-3 * 0.5 * d2);
            }
        }
    }
    let pair = pair.expect("the near-degenerate pair needs a full block");
    // Forcing the diagonal: a budget so loose that nothing merges.
    let forced = certify_modal_damping(&omega, &c, &loose()).unwrap();
    let mut worst_decay = 0.0_f64;
    for block in &forced.blocks {
        if let ModalDampingBlock::Diagonal {
            mode, eigenvalue, ..
        } = block
            && *mode < 2
        {
            let truth_pole = poles
                .iter()
                .copied()
                .filter(|p| p.im > 0.0)
                .min_by(|a, b| {
                    (*a - *eigenvalue)
                        .abs()
                        .total_cmp(&(*b - *eigenvalue).abs())
                })
                .unwrap();
            worst_decay =
                worst_decay.max((truth_pole.re - eigenvalue.re).abs() / truth_pole.re.abs());
        }
    }
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g1_clustered_pair_triggers_full_block_fallback\",\"tier\":\"G1\",\"inputs\":{{\"fixture\":\"synthetic 3-mode, 0.3% split\"}},\"observed\":{{\"block_poles\":\"{pair:?}\",\"forced_diagonal_rel_decay_error\":{worst_decay:e}}},\"tolerance\":{{\"decay\":1e-3,\"cents\":0.01}}}}"
    );
    assert!(
        worst_decay > 1e-2,
        "forcing the diagonal must violate the 1e-3 decay budget measurably: {worst_decay:e}"
    );
}

#[test]
fn g0_cluster_cap_and_cap_plus_one() {
    // Every mode coupled to every other inside a 1e-6 relative band: the
    // minimal certified cluster is everything.
    let system = |n: usize| {
        let omega: Vec<f64> = (0..n).map(|i| 1000.0 * (1.0 + 1e-6 * i as f64)).collect();
        let mut c = vec![0.3; n * n];
        for i in 0..n {
            c[n * i + i] = 2.0 + 1e-3 * i as f64;
        }
        (omega, c)
    };
    let (omega, c) = system(MAX_DAMPING_CLUSTER);
    let result = certify_modal_damping(&omega, &c, &budget()).expect("64-mode cluster");
    assert_eq!(result.blocks.len(), 1);
    assert!(
        matches!(&result.blocks[0], ModalDampingBlock::Full { modes, .. } if modes.len() == 64)
    );
    let (omega, c) = system(MAX_DAMPING_CLUSTER + 1);
    match certify_modal_damping(&omega, &c, &budget()) {
        Err(ModalError::DampingClusterTooLarge { modes, cap }) => {
            assert_eq!(cap, MAX_DAMPING_CLUSTER);
            assert!(modes.len() > MAX_DAMPING_CLUSTER, "{}", modes.len());
        }
        other => panic!("65 coupled modes must refuse: {other:?}"),
    }
    for cap in [0, MAX_DAMPING_CLUSTER + 1] {
        let bad = ModalDampingBudget {
            max_cluster_size: cap,
            ..budget()
        };
        assert!(matches!(
            certify_modal_damping(&[1.0], &[0.0], &bad),
            Err(ModalError::DampingInput { .. })
        ));
    }
}

#[test]
fn g0_budget_inputs_are_required_and_finite() {
    let omega = [10.0, 20.0];
    let c = [0.1, 0.01, 0.01, 0.2];
    for value in [0.0, -1e-3, f64::NAN, f64::INFINITY] {
        for field in 0..2 {
            let mut b = budget();
            if field == 0 {
                b.max_rel_decay_error = value;
            } else {
                b.max_freq_error_cents = value;
            }
            assert!(
                matches!(
                    certify_modal_damping(&omega, &c, &b),
                    Err(ModalError::DampingInput { .. })
                ),
                "budget field {field} = {value} must refuse"
            );
        }
    }
    let refuses = |omega: &[f64], c: &[f64]| {
        matches!(
            certify_modal_damping(omega, c, &budget()),
            Err(ModalError::DampingInput { .. })
        )
    };
    assert!(refuses(&[10.0, 0.0], &c), "zero frequency");
    assert!(refuses(&[10.0, f64::NAN], &c), "non-finite frequency");
    assert!(refuses(&omega, &c[..3]), "shape");
    assert!(refuses(&omega, &[0.1, 0.01, 0.02, 0.2]), "asymmetric");
    assert!(refuses(&omega, &[-0.1, 0.0, 0.0, 0.2]), "negative diagonal");
    assert!(
        refuses(&omega, &[0.1, f64::NAN, f64::NAN, 0.2]),
        "non-finite"
    );
    assert!(refuses(&[], &[]), "empty");
    assert!(matches!(
        certify_modal_damping(&omega, &[20.0, 0.0, 0.0, 0.2], &budget()),
        Err(ModalError::DampingNotUnderdamped { modes }) if modes == vec![0]
    ));
}

#[test]
fn g0_bessel_far_coupling_bound_dominates_true_norm() {
    // A complete 40-mode space; the certifier sees the lowest 24 through a
    // window, and Bessel's inequality bounds what it does not read.
    let full_n = 40;
    let mut s = StreamKey {
        seed: 0x4d30_3257_3561,
        kernel: 7,
        tile: 0,
    }
    .stream();
    let mut omega_full: Vec<f64> = (0..full_n)
        .map(|i| 100.0 * (1.0 + i as f64) + 30.0 * s.next_f64())
        .collect();
    omega_full.sort_by(f64::total_cmp);
    let zeta: Vec<f64> = (0..full_n)
        .map(|_| log_uniform(&mut s, 1e-4, 1e-3))
        .collect();
    let c_full = random_damping(&omega_full, &zeta, 6, 1.0, &mut s);
    let column = |j: usize| {
        det::sqrt(
            (0..full_n)
                .map(|k| c_full[full_n * j + k] * c_full[full_n * j + k])
                .sum(),
        )
    };
    let mut worst_margin = f64::INFINITY;
    let mut floor = 0.0_f64;
    for (retained, window) in [(24usize, 350.0), (24, 120.0), (full_n, 1e6)] {
        let omega = &omega_full[..retained];
        let bounds: Vec<f64> = (0..retained).map(column).collect();
        let near = |i: usize, j: usize| c_full[full_n * i + j];
        let result =
            certify_modal_damping_windowed(omega, &near, &bounds, window, &loose()).unwrap();
        assert_eq!(result.far_coupling.window_rad_s, Some(window));
        for j in 0..retained {
            let true_far = det::sqrt(
                (0..retained)
                    .filter(|&k| k != j && (omega[j] - omega[k]).abs() >= window)
                    .map(|k| c_full[full_n * j + k] * c_full[full_n * j + k])
                    .sum(),
            );
            let bound = result.far_coupling.far_row_norm[j];
            assert!(
                bound >= true_far,
                "mode {j}: far bound {bound} < true {true_far}"
            );
            worst_margin = worst_margin.min(bound - true_far);
            if retained == full_n {
                // A difference of squares keeps about half the digits: the
                // floor is ~√ε·‖C e_j‖, so 1e-10 relative is out of reach.
                floor = floor.max(bound / column(j));
                assert!(
                    bound <= 1e-6 * column(j),
                    "complete space, full window: bound {bound} should vanish"
                );
            }
        }
    }
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g0_bessel_far_coupling_bound_dominates_true_norm\",\"tier\":\"G0\",\"inputs\":{{\"fixture\":\"synthetic complete 40-mode space, 24 retained\"}},\"observed\":{{\"min_bound_minus_true\":{worst_margin:e},\"complete_space_relative_floor\":{floor:e}}}}}"
    );
}

#[test]
fn g0_windowed_certificate_encloses_dense_truth() {
    let (mut checked, mut violations, mut certified_trials) = (0usize, 0usize, 0usize);
    for trial in 0..128_u32 {
        let mut s = StreamKey {
            seed: 0x4d30_3257_3561,
            kernel: 8,
            tile: trial,
        }
        .stream();
        let n = 6 + s.next_below(20) as usize;
        let mut omega: Vec<f64> = (0..n)
            .map(|_| log_uniform(&mut s, 300.0, 30_000.0))
            .collect();
        omega.sort_by(f64::total_cmp);
        let zeta: Vec<f64> = (0..n).map(|_| log_uniform(&mut s, 1e-4, 1e-2)).collect();
        let rank = 1 + s.next_below(n as u64) as usize;
        let c = random_damping(&omega, &zeta, rank, 0.05, &mut s);
        let bounds: Vec<f64> = (0..n)
            .map(|j| det::sqrt((0..n).map(|k| c[n * j + k] * c[n * j + k]).sum()))
            .map(|b| b * (1.0 + 1e-12))
            .collect();
        let window = log_uniform(&mut s, 50.0, 5000.0);
        let near = |i: usize, j: usize| c[n * i + j];
        let moderate = ModalDampingBudget {
            max_rel_decay_error: 0.1,
            max_freq_error_cents: 1.0,
            max_cluster_size: MAX_DAMPING_CLUSTER,
        };
        let result = match certify_modal_damping_windowed(&omega, &near, &bounds, window, &moderate)
        {
            Ok(result) => {
                certified_trials += 1;
                result
            }
            Err(ModalError::DampingBlockSpansWindow { .. }) => continue,
            Err(e) => panic!("trial {trial}: {e}"),
        };
        let (poles, slack, _) = truth(&omega, &c);
        for (mu, r) in certified_poles(&result) {
            checked += 1;
            if nearest(&poles, mu) > r + slack {
                violations += 1;
            }
        }
    }
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g0_windowed_certificate_encloses_dense_truth\",\"tier\":\"G0\",\"observed\":{{\"certified_trials\":{certified_trials},\"of\":128,\"poles_checked\":{checked},\"violations\":{violations}}}}}"
    );
    assert!(
        certified_trials >= 64,
        "only {certified_trials}/128 trials certified"
    );
    assert_eq!(violations, 0);
}

#[test]
fn g0_window_neighbours_are_exactly_the_entries_read() {
    // A caller bounding the unread entries itself must exclude exactly what
    // the certifier reads; record every read and compare.
    let mut s = StreamKey {
        seed: 0x4d30_3257_3561,
        kernel: 9,
        tile: 0,
    }
    .stream();
    let n = 30;
    let mut omega: Vec<f64> = (0..n).map(|_| log_uniform(&mut s, 300.0, 3000.0)).collect();
    // A tie and a gap of exactly the window: both sides of `<` are exercised.
    omega[7] = omega[3];
    omega[4] = 1000.0;
    omega[9] = 1100.0;
    let zeta: Vec<f64> = (0..n).map(|_| log_uniform(&mut s, 1e-4, 1e-3)).collect();
    let c = random_damping(&omega, &zeta, 4, 1e-3, &mut s);
    let window = 100.0;
    let neighbours = damping_window_neighbours(&omega, window).unwrap();
    for (i, row) in neighbours.iter().enumerate() {
        assert!(row.windows(2).all(|w| w[0] < w[1]), "ascending");
        for k in 0..n {
            let inside = k != i && (omega[i] - omega[k]).abs() < window;
            assert_eq!(row.contains(&k), inside, "mode {i} vs {k}");
            assert_eq!(neighbours[k].contains(&i), row.contains(&k), "symmetric");
        }
    }
    assert!(
        !neighbours[4].contains(&9),
        "a gap of exactly the window is far"
    );
    assert!(neighbours[3].contains(&7), "a tie is near");
    let read = std::cell::RefCell::new(std::collections::BTreeSet::new());
    let near = |i: usize, j: usize| {
        read.borrow_mut().insert((i, j));
        c[n * i + j]
    };
    let far: Vec<f64> = (0..n)
        .map(|j| {
            det::sqrt(
                (0..n)
                    .filter(|&k| k != j && !neighbours[j].contains(&k))
                    .map(|k| c[n * j + k] * c[n * j + k])
                    .sum(),
            )
        })
        .collect();
    certify_modal_damping_windowed_far_rows(&omega, &near, &far, window, &loose()).unwrap();
    let expected: std::collections::BTreeSet<(usize, usize)> = (0..n)
        .flat_map(|i| {
            neighbours[i]
                .iter()
                .map(move |&k| (i, k))
                .chain(std::iter::once((i, i)))
        })
        .collect();
    assert_eq!(*read.borrow(), expected);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            damping_window_neighbours(&omega, bad),
            Err(ModalError::DampingInput { .. })
        ));
        assert!(matches!(
            damping_window_neighbours(&[1.0, bad], 1.0),
            Err(ModalError::DampingInput { .. })
        ));
    }
}

#[test]
fn g0_far_row_bound_path_encloses_dense_truth_and_refuses_bad_bounds() {
    // The same trials as the windowed column-norm path, with each row's far
    // norm summed directly: no difference of squares, and a far bound that
    // vanishes exactly when the window covers everything.
    let (mut checked, mut violations, mut certified_trials) = (0usize, 0usize, 0usize);
    let mut full_window_far = 0.0_f64;
    for trial in 0..128_u32 {
        let mut s = StreamKey {
            seed: 0x4d30_3257_3561,
            kernel: 8,
            tile: trial,
        }
        .stream();
        let n = 6 + s.next_below(20) as usize;
        let mut omega: Vec<f64> = (0..n)
            .map(|_| log_uniform(&mut s, 300.0, 30_000.0))
            .collect();
        omega.sort_by(f64::total_cmp);
        let zeta: Vec<f64> = (0..n).map(|_| log_uniform(&mut s, 1e-4, 1e-2)).collect();
        let rank = 1 + s.next_below(n as u64) as usize;
        let c = random_damping(&omega, &zeta, rank, 0.05, &mut s);
        let window = log_uniform(&mut s, 50.0, 5000.0);
        let moderate = ModalDampingBudget {
            max_rel_decay_error: 0.1,
            max_freq_error_cents: 1.0,
            max_cluster_size: MAX_DAMPING_CLUSTER,
        };
        let near = |i: usize, j: usize| c[n * i + j];
        for w in [window, 2.0 * (omega[n - 1] - omega[0])] {
            let neighbours = damping_window_neighbours(&omega, w).unwrap();
            let far: Vec<f64> = (0..n)
                .map(|j| {
                    det::sqrt(
                        (0..n)
                            .filter(|&k| k != j && !neighbours[j].contains(&k))
                            .map(|k| c[n * j + k] * c[n * j + k])
                            .sum(),
                    )
                })
                .collect();
            let result =
                match certify_modal_damping_windowed_far_rows(&omega, &near, &far, w, &moderate) {
                    Ok(result) => result,
                    Err(ModalError::DampingBlockSpansWindow { .. }) if w == window => continue,
                    Err(e) => panic!("trial {trial}, window {w}: {e}"),
                };
            assert_eq!(result.far_coupling.far_row_norm, far, "bounds pass through");
            if w != window {
                certified_trials += 1;
                full_window_far = full_window_far.max(far.iter().copied().fold(0.0, f64::max));
            }
            let (poles, slack, _) = truth(&omega, &c);
            for (mu, r) in certified_poles(&result) {
                checked += 1;
                if nearest(&poles, mu) > r + slack {
                    violations += 1;
                }
            }
        }
    }
    let omega = [10.0, 20.0];
    let c = [0.1, 0.01, 0.01, 0.2];
    let near = |i: usize, j: usize| c[2 * i + j];
    for far in [
        vec![0.0],
        vec![-1e-3, 0.0],
        vec![f64::NAN, 0.0],
        vec![0.0, f64::INFINITY],
    ] {
        assert!(
            matches!(
                certify_modal_damping_windowed_far_rows(&omega, &near, &far, 5.0, &budget()),
                Err(ModalError::DampingInput { .. })
            ),
            "far bound {far:?} must refuse"
        );
    }
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g0_far_row_bound_path_encloses_dense_truth_and_refuses_bad_bounds\",\"tier\":\"G0\",\"observed\":{{\"full_window_certified_trials\":{certified_trials},\"of\":128,\"poles_checked\":{checked},\"violations\":{violations},\"full_window_max_far_bound\":{full_window_far:e}}}}}"
    );
    assert_eq!(certified_trials, 128, "the full window spans no block");
    assert_eq!(
        full_window_far, 0.0,
        "nothing is unread when the window covers everything"
    );
    assert_eq!(violations, 0);
}

/// `Σ aᵢbᵢ` as an unevaluated pair `hi + lo`, about ε² relative to `Σ|aᵢbᵢ|`:
/// every product and partial sum keeps its exact rounding error (Dot2,
/// Ogita, Rump and Oishi 2005).
fn dot2(a: &[f64], b: &[f64]) -> (f64, f64) {
    let (mut sum, mut carry) = (0.0_f64, 0.0_f64);
    for (&x, &y) in a.iter().zip(b) {
        let p = x * y;
        let product_error = x.mul_add(y, -p);
        let t = sum + p;
        let z = t - sum;
        carry += (sum - (t - z)) + (p - z) + product_error;
        sum = t;
    }
    (sum, carry)
}

#[test]
fn g0_residual_far_row_bound_encloses_exact_far_norm() {
    // The residual route must still be an UPPER bound. Each trial builds a
    // random lumped space: a diagonal mass spanning six decades (rotary-like
    // DOFs), an SPD stiffness and a PSD C_FE. Modes come from the dense
    // generalized solve, all of them or a lowest subset retained. In a
    // quarter of the trials the retained shapes are perturbed by 1e-14..1e-6,
    // so the charged basis defect carries weight. Truth is every unformed
    // retained entry φ_kᵀC_FEφ_j, computed in compensated arithmetic on the
    // same f64 inputs.
    let started = std::time::Instant::now();
    let (mut checked, mut far_violations, mut column_violations) = (0usize, 0usize, 0usize);
    let (mut worst, mut worst_trial) = (0.0_f64, 0u32);
    let mut tight_ratios = Vec::new();
    let (mut complete_trials, mut perturbed_trials, mut max_defect) = (0usize, 0usize, 0.0_f64);
    for trial in 0..4096_u32 {
        let mut s = StreamKey {
            seed: 0x4d30_3257_3561,
            kernel: 10,
            tile: trial,
        }
        .stream();
        let dofs = 4 + s.next_below(41) as usize;
        let mass: Vec<f64> = (0..dofs).map(|_| log_uniform(&mut s, 1e-6, 1.0)).collect();
        let a: Vec<f64> = (0..dofs * dofs).map(|_| 2.0 * s.next_f64() - 1.0).collect();
        let (mut k, mut m) = (vec![0.0; dofs * dofs], vec![0.0; dofs * dofs]);
        for i in 0..dofs {
            m[dofs * i + i] = mass[i];
            for j in 0..=i {
                let sum: f64 = (0..dofs).map(|r| a[dofs * i + r] * a[dofs * j + r]).sum();
                k[dofs * i + j] = sum;
                k[dofs * j + i] = sum;
            }
            k[dofs * i + i] += 1e-3 * dofs as f64;
        }
        let pairs = fs_modal::eigh_gen_dense(&k, &m, dofs).expect("dense modes");
        let complete = s.next_f64() < 0.5;
        let n = if complete {
            dofs
        } else {
            2 + s.next_below(dofs as u64 - 1) as usize
        };
        complete_trials += usize::from(n == dofs);
        let omega: Vec<f64> = pairs[..n].iter().map(|p| det::sqrt(p.lambda)).collect();
        let mut phi: Vec<Vec<f64>> = pairs[..n].iter().map(|p| p.phi.clone()).collect();
        if s.next_f64() < 0.25 {
            perturbed_trials += 1;
            let level = det::exp(-det::ln(10.0) * (6.0 + 8.0 * s.next_f64()));
            for v in &mut phi {
                for x in v.iter_mut() {
                    *x *= 1.0 + level * (2.0 * s.next_f64() - 1.0);
                }
            }
        }
        let rank = 1 + s.next_below(dofs as u64) as usize;
        let b: Vec<f64> = (0..dofs * rank).map(|_| 2.0 * s.next_f64() - 1.0).collect();
        let mut c_fe = vec![0.0; dofs * dofs];
        for i in 0..dofs {
            for j in 0..=i {
                let sum: f64 = (0..rank).map(|r| b[rank * i + r] * b[rank * j + r]).sum();
                c_fe[dofs * i + j] = sum;
                c_fe[dofs * j + i] = sum;
            }
        }
        // x_j = C_FE φ_j in f64, and its rounding: |δx| ≤ γ_dofs |C_FE||φ_j|.
        let gamma = dofs as f64 * f64::EPSILON / (1.0 - dofs as f64 * f64::EPSILON);
        let applied: Vec<Vec<f64>> = phi
            .iter()
            .map(|p| {
                (0..dofs)
                    .map(|r| (0..dofs).map(|q| c_fe[dofs * r + q] * p[q]).sum())
                    .collect()
            })
            .collect();
        let applied_error: Vec<f64> = phi
            .iter()
            .map(|p| {
                let size_sq: f64 = (0..dofs)
                    .map(|r| {
                        let size: f64 = (0..dofs).map(|q| (c_fe[dofs * r + q] * p[q]).abs()).sum();
                        size * size / mass[r]
                    })
                    .sum();
                1.01 * gamma * det::sqrt(size_sq)
            })
            .collect();
        let window = log_uniform(&mut s, 1e-3 * omega[n - 1], 2.0 * omega[n - 1]);
        let bounds = fs_modal::residual_far_row_bounds(
            &omega,
            window,
            &mass,
            &phi,
            &applied,
            &applied_error,
        )
        .expect("bounds");
        max_defect = max_defect.max(bounds.basis_defect);
        let neighbours = damping_window_neighbours(&omega, window).unwrap();
        // Exact x_j as hi + lo pairs, then exact retained entries C_kj.
        let exact: Vec<Vec<f64>> = phi
            .iter()
            .map(|p| {
                let (mut hi, mut lo) = (vec![0.0; dofs], vec![0.0; dofs]);
                for r in 0..dofs {
                    (hi[r], lo[r]) = dot2(&c_fe[dofs * r..dofs * (r + 1)], p);
                }
                hi.extend(lo);
                hi
            })
            .collect();
        for j in 0..n {
            let entries: Vec<f64> = (0..n)
                .map(|kk| {
                    let doubled: Vec<f64> = phi[kk].iter().chain(&phi[kk]).copied().collect();
                    let (hi, lo) = dot2(&doubled, &exact[j]);
                    hi + lo
                })
                .collect();
            let far: Vec<f64> = (0..n)
                .filter(|&kk| kk != j && neighbours[j].binary_search(&kk).is_err())
                .map(|kk| entries[kk])
                .collect();
            let (hi, lo) = dot2(&far, &far);
            let true_far = det::sqrt(hi + lo);
            let (hi, lo) = dot2(&entries, &entries);
            let true_column = det::sqrt(hi + lo);
            checked += 1;
            let bound = bounds.far_row_norm[j];
            if true_far > bound {
                far_violations += 1;
            }
            if true_column > bounds.column_norm[j] {
                column_violations += 1;
            }
            let ratio = true_far / bound;
            if ratio > worst {
                (worst, worst_trial) = (ratio, trial);
            }
            // Complete spaces are where the bound can be tight.
            if n == dofs && true_far > 0.0 {
                tight_ratios.push(ratio);
            }
        }
    }
    tight_ratios.sort_by(f64::total_cmp);
    let median = tight_ratios
        .get(tight_ratios.len() / 2)
        .copied()
        .unwrap_or(0.0);
    println!(
        "{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g0_residual_far_row_bound_encloses_exact_far_norm\",\"tier\":\"G0\",\"inputs\":{{\"trials\":4096,\"dofs\":[4,44],\"mass_decades\":6,\"complete_trials\":{complete_trials},\"perturbed_trials\":{perturbed_trials}}},\"provenance\":\"fs-rand seed 0x4d3032573561 kernel 10\",\"observed\":{{\"rows_checked\":{checked},\"far_violations\":{far_violations},\"column_violations\":{column_violations},\"max_true_over_bound\":{worst:e},\"worst_trial\":{worst_trial},\"complete_space_median_true_over_bound\":{median:e},\"complete_space_samples\":{},\"max_basis_defect\":{max_defect:e}}},\"tolerance\":{{\"true_over_bound\":1.0}},\"verdict\":\"{}\",\"elapsed_s\":{:.1}}}",
        tight_ratios.len(),
        if far_violations == 0 && column_violations == 0 {
            "pass"
        } else {
            "fail"
        },
        started.elapsed().as_secs_f64()
    );
    assert_eq!(
        far_violations, 0,
        "far-row bounds must enclose the exact far norm"
    );
    assert_eq!(
        column_violations, 0,
        "column bounds must enclose the exact column norm"
    );
    // Coverage, not a quality claim: both kinds of trial must be live.
    assert!(complete_trials >= 1000 && perturbed_trials >= 500);
}

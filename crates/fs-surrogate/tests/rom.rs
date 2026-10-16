//! G0/G1/G2 battery for the Koopman/DMD and POD-Galerkin+DEIM reduced models:
//! exact spectral recovery on linear systems, exact finite-section Koopman
//! recovery on an invariant-subspace system, conformal forecast coverage on
//! held-out nonlinear trajectories (certify-or-escalate wired), and a
//! hyper-reduced Allen–Cahn ROM against its full-order model.

use fs_surrogate::{
    Decision, DmdRank, Forecaster, GalerkinDeimRom, certify_or_escalate, deim, dmd, edmd,
    forecast_bands,
};

/// Deterministic LCG in [0, 1).
fn lcg(s: &mut u64) -> f64 {
    *s = s
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    ((*s >> 11) as f64) / (1u64 << 53) as f64
}

/// Orthonormal n×k columns by Gram–Schmidt on LCG vectors.
fn orthonormal(n: usize, k: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut s = seed;
    let mut out: Vec<Vec<f64>> = Vec::new();
    while out.len() < k {
        let mut v: Vec<f64> = (0..n).map(|_| lcg(&mut s) - 0.5).collect();
        for u in &out {
            let c: f64 = u.iter().zip(&v).map(|(a, b)| a * b).sum();
            for (vi, ui) in v.iter_mut().zip(u) {
                *vi -= c * ui;
            }
        }
        let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        out.push(v.iter().map(|x| x / nrm).collect());
    }
    out
}

#[test]
fn dmd_g0_recovers_embedded_linear_spectrum() {
    // Latent 4-D dynamics: a damped rotation (0.95 e^{±0.3i}) plus 0.8, 0.6,
    // embedded in ℝ⁴⁰ by an orthonormal map.
    let (r, th) = (0.95f64, 0.3f64);
    let a = [
        [r * th.cos(), -r * th.sin(), 0.0, 0.0],
        [r * th.sin(), r * th.cos(), 0.0, 0.0],
        [0.0, 0.0, 0.8, 0.0],
        [0.0, 0.0, 0.0, 0.6],
    ];
    let q = orthonormal(40, 4, 7);
    let mut z = [1.0, -0.5, 0.7, 0.3];
    let mut snaps = Vec::new();
    for _ in 0..30 {
        snaps.push(
            (0..40)
                .map(|i| (0..4).map(|k| q[k][i] * z[k]).sum())
                .collect::<Vec<f64>>(),
        );
        let zn: Vec<f64> = (0..4)
            .map(|i| (0..4).map(|j| a[i][j] * z[j]).sum())
            .collect();
        z.copy_from_slice(&zn);
    }
    let model = dmd(&snaps, 0.1, DmdRank::Fixed(4)).unwrap();
    assert_eq!(model.rank(), 4);
    let ev = model.eigenvalues();
    assert!((ev[0].0.hypot(ev[0].1) - 0.95).abs() < 1e-9, "{ev:?}");
    assert!((ev[0].1.abs() - r * th.sin()).abs() < 1e-9);
    assert!((ev[2].0 - 0.8).abs() < 1e-9 && (ev[3].0 - 0.6).abs() < 1e-9);
    assert!(model.spectral_radius() < 1.0);
    assert!(model.fit_residual() < 1e-10, "{}", model.fit_residual());
    // Continuous spectrum: growth ln(0.95)/0.1, frequency 0.3/0.1.
    let cs = model.continuous_spectrum();
    assert!((cs[0].0 - 0.95f64.ln() / 0.1).abs() < 1e-8);
    assert!((cs[0].1.abs() - 3.0).abs() < 1e-8);
    // Forecast reproduces the snapshots.
    let f = model.forecast(&snaps[0], 20);
    for (h, fh) in f.iter().enumerate() {
        let e: f64 = fh
            .iter()
            .zip(&snaps[h + 1])
            .map(|(p, t)| (p - t).abs())
            .fold(0.0, f64::max);
        assert!(e < 1e-9, "horizon {h}: {e}");
    }
}

#[test]
fn edmd_g0_finite_koopman_invariant_system_is_exact() {
    // x⁺ = λx, y⁺ = μy + c x²: span{x, y, x²} is Koopman-invariant, so a
    // degree-2 dictionary reproduces the dynamics EXACTLY and the spectrum
    // contains λ, μ, λ².
    let (lam, mu, c) = (0.9f64, 0.5f64, 0.3f64);
    let mut pairs = Vec::new();
    let mut s = 11u64;
    for _ in 0..60 {
        let x = 2.0 * lcg(&mut s) - 1.0;
        let y = 2.0 * lcg(&mut s) - 1.0;
        pairs.push((vec![x, y], vec![lam * x, mu * y + c * x * x]));
    }
    let model = edmd(&pairs, 2, 0.0, 1.0).unwrap();
    assert_eq!(model.dictionary_size(), 6);
    let x0 = vec![0.7, -0.4];
    let traj = model.forecast(&x0, 10);
    let (mut x, mut y) = (x0[0], x0[1]);
    for p in &traj {
        let (xn, yn) = (lam * x, mu * y + c * x * x);
        x = xn;
        y = yn;
        assert!(
            (p[0] - x).abs() < 1e-10 && (p[1] - y).abs() < 1e-10,
            "{p:?} vs ({x}, {y})"
        );
    }
    let ev = model.eigenvalues();
    for target in [lam, mu, lam * lam] {
        assert!(
            ev.iter()
                .any(|(re, im)| (re - target).abs() < 1e-8 && im.abs() < 1e-8),
            "missing eigenvalue {target}: {ev:?}"
        );
    }
}

/// Nonlinear pitch dynamics θ̈ = −kθ(1 − θ²/θs²) − dθ̇ sampled by RK4.
fn pitch_traj(theta0: f64, omega0: f64, dt: f64, steps: usize) -> Vec<Vec<f64>> {
    let (k, d, ths) = (2.5f64, 0.8f64, 0.35f64);
    let f = |s: [f64; 2]| {
        [
            s[1],
            -k * s[0] * (1.0 - s[0] * s[0] / (ths * ths)) - d * s[1],
        ]
    };
    let mut s = [theta0, omega0];
    let mut out = vec![s.to_vec()];
    let sub = 10;
    let h = dt / f64::from(sub);
    for _ in 0..steps {
        for _ in 0..sub {
            let k1 = f(s);
            let k2 = f([s[0] + 0.5 * h * k1[0], s[1] + 0.5 * h * k1[1]]);
            let k3 = f([s[0] + 0.5 * h * k2[0], s[1] + 0.5 * h * k2[1]]);
            let k4 = f([s[0] + h * k3[0], s[1] + h * k3[1]]);
            for i in 0..2 {
                s[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
        }
        out.push(s.to_vec());
    }
    out
}

#[test]
fn koopman_g2_conformal_forecast_bands_cover_fresh_trajectories() {
    let dt = 0.05;
    let horizon = 20;
    let mut s = 0x5eed_u64;
    let mut draw = |n: usize| -> Vec<Vec<Vec<f64>>> {
        (0..n)
            .map(|_| {
                let th = 0.25 * (2.0 * lcg(&mut s) - 1.0);
                let om = 0.4 * (2.0 * lcg(&mut s) - 1.0);
                pitch_traj(th, om, dt, horizon)
            })
            .collect()
    };
    let train = draw(40);
    // Simultaneous bands at α/H = 0.005 need ≥ 199 calibration trajectories
    // for a finite order statistic.
    let calib = draw(240);
    let fresh = draw(200);
    let pairs: Vec<(Vec<f64>, Vec<f64>)> = train
        .iter()
        .flat_map(|t| t.windows(2).map(|w| (w[0].clone(), w[1].clone())))
        .collect();
    let lin = edmd(&pairs, 1, 1e-12, dt).unwrap();
    let cubic = edmd(&pairs, 3, 1e-12, dt).unwrap();
    // The cubic dictionary captures the stall nonlinearity: smaller error.
    let err = |m: &dyn Forecaster| -> f64 {
        fresh
            .iter()
            .map(|t| {
                let f = m.forecast(&t[0], horizon);
                f.iter()
                    .zip(&t[1..])
                    .map(|(p, q)| (p[0] - q[0]).hypot(p[1] - q[1]))
                    .fold(0.0, f64::max)
            })
            .sum::<f64>()
            / fresh.len() as f64
    };
    let (e_lin, e_cub) = (err(&lin), err(&cubic));
    assert!(
        e_cub < 0.5 * e_lin,
        "cubic EDMD {e_cub} not better than linear {e_lin}"
    );
    // The linear fit's spectrum is the damped pitch pair (|λ| < 1); the cubic
    // dictionary's finite section may carry spurious eigenvalues (EDMD
    // spectral pollution — the monomial span is not Koopman-invariant), so
    // its spectrum is NOT asserted: forecasts and their bands are the claim.
    let lin_ev = lin.eigenvalues();
    assert!(
        lin_ev.iter().all(|(re, im)| re.hypot(*im) < 1.0),
        "{lin_ev:?}"
    );
    // Simultaneous 90% bands: fresh coverage at the declared level.
    let bands = forecast_bands(&cubic, &calib, horizon, 0.1, true);
    assert!(bands.max_half_width().is_finite(), "{bands:?}");
    let covered = fresh
        .iter()
        .filter(|t| bands.covers(&cubic.forecast(&t[0], horizon), &t[1..]))
        .count();
    let cov = covered as f64 / fresh.len() as f64;
    assert!(cov >= 0.85, "simultaneous coverage {cov} below 0.9 − slack");
    // Certify-or-escalate on the whole-horizon band.
    let widest = bands.bands.iter().copied().fold(bands.bands[0], |a, b| {
        if b.half_width > a.half_width { b } else { a }
    });
    assert!(matches!(
        certify_or_escalate(&widest, true, widest.half_width * 2.0),
        Decision::UseSurrogate { .. }
    ));
    assert!(matches!(
        certify_or_escalate(&widest, true, widest.half_width * 0.5),
        Decision::Escalate { .. }
    ));
}

/// Allen–Cahn ẋ = ν Δ_h x + x − x³ on (0, 1), homogeneous Dirichlet.
struct AllenCahn {
    n: usize,
    nu: f64,
}

impl AllenCahn {
    fn lap(&self, x: &[f64]) -> Vec<f64> {
        let h = 1.0 / (self.n as f64 + 1.0);
        let c = self.nu / (h * h);
        (0..self.n)
            .map(|i| {
                let l = if i > 0 { x[i - 1] } else { 0.0 };
                let r = if i + 1 < self.n { x[i + 1] } else { 0.0 };
                c * (l - 2.0 * x[i] + r) + x[i]
            })
            .collect()
    }
    fn rhs(&self, x: &[f64]) -> Vec<f64> {
        self.lap(x)
            .iter()
            .zip(x)
            .map(|(l, xi)| l - xi * xi * xi)
            .collect()
    }
    fn run(&self, x0: &[f64], dt: f64, steps: usize) -> Vec<Vec<f64>> {
        let mut x = x0.to_vec();
        let mut out = vec![x.clone()];
        for _ in 0..steps {
            let k1 = self.rhs(&x);
            let k2 = self.rhs(
                &x.iter()
                    .zip(&k1)
                    .map(|(a, k)| a + 0.5 * dt * k)
                    .collect::<Vec<_>>(),
            );
            let k3 = self.rhs(
                &x.iter()
                    .zip(&k2)
                    .map(|(a, k)| a + 0.5 * dt * k)
                    .collect::<Vec<_>>(),
            );
            let k4 = self.rhs(
                &x.iter()
                    .zip(&k3)
                    .map(|(a, k)| a + dt * k)
                    .collect::<Vec<_>>(),
            );
            for i in 0..x.len() {
                x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            out.push(x.clone());
        }
        out
    }
}

fn initial(n: usize, a: f64, b: f64) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let s = (i as f64 + 1.0) / (n as f64 + 1.0);
            a * (std::f64::consts::PI * s).sin() + b * (3.0 * std::f64::consts::PI * s).sin()
        })
        .collect()
}

#[test]
fn deim_g1_allen_cahn_rom_tracks_full_model_with_few_samples() {
    // ν/h² ≈ 202: the stiffest Laplacian mode (≈ −808) sits inside RK4's
    // stability interval (≈ −2.78/dt = −1392) at dt = 2e-3.
    let model = AllenCahn { n: 200, nu: 0.005 };
    let (dt, steps) = (2e-3, 500);
    // Training: two initial conditions of the family.
    // Every 5th state of each training run (the ROM needs the manifold, not
    // every step).
    let mut snaps: Vec<Vec<f64>> = model
        .run(&initial(200, 0.3, 0.2), dt, steps)
        .into_iter()
        .step_by(5)
        .collect();
    snaps.extend(
        model
            .run(&initial(200, 0.1, -0.15), dt, steps)
            .into_iter()
            .step_by(5),
    );
    let rom = GalerkinDeimRom::new(
        &snaps,
        1.0 - 1e-9,
        16,
        1.0 - 1e-9,
        16,
        |x| model.lap(x),
        &vec![0.0; 200],
        |x| -x * x * x,
    )
    .unwrap();
    assert!(rom.state_rank() <= 16 && rom.deim().rank() <= 16);
    // Online: an UNSEEN initial condition inside the family.
    let x0 = initial(200, 0.2, 0.05);
    let full = model.run(&x0, dt, steps);
    let red = rom.integrate(&x0, dt, steps);
    let mut worst = 0.0f64;
    for (zf, xf) in red.iter().zip(&full) {
        let x = rom.lift(zf);
        let e = x
            .iter()
            .zip(xf)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        let nrm = xf.iter().map(|b| b * b).sum::<f64>().sqrt();
        worst = worst.max(e / nrm.max(1e-12));
    }
    assert!(worst < 1e-2, "worst relative ROM error {worst}");
    // DEIM's a-priori bound holds on every training nonlinearity snapshot.
    let d = rom.deim();
    for s in snaps.iter().step_by(5) {
        let g: Vec<f64> = s.iter().map(|x| -x * x * x).collect();
        let (e, bnd) = (d.interpolation_error(&g), d.error_bound(&g));
        assert!(
            e <= bnd * (1.0 + 1e-9) + 1e-14,
            "DEIM error {e} exceeds bound {bnd}"
        );
    }
    // Distinct sample rows, far fewer than n.
    let mut idx = d.indices().to_vec();
    idx.sort_unstable();
    idx.dedup();
    assert_eq!(idx.len(), d.rank());
    assert!(d.rank() * 10 <= 200);
}

#[test]
fn deim_g5_replays_bit_identically() {
    let snaps: Vec<Vec<f64>> = (0..12)
        .map(|k| {
            (0..50)
                .map(|i| (f64::from(i * (k + 1)) * 0.07).sin())
                .collect()
        })
        .collect();
    let a = deim(&snaps, 0.999_999, 8).unwrap();
    let b = deim(&snaps, 0.999_999, 8).unwrap();
    assert_eq!(a, b);
    assert_eq!(a.error_constant().to_bits(), b.error_constant().to_bits());
}

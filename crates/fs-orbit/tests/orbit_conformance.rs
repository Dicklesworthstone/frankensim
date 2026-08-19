//! fs-orbit conformance battery (music bead
//! `frankensim-music-v8-root-3ez8g.11.1`): ORACLES WITH NO MUSIC IN
//! THEM, exactly as the I09 charter demands.
//!
//! - ob-001: conservative Duffing backbone vs the fs-nlmodal pinned
//!   first-order law (the amplitude-pinned HB mode).
//! - ob-002: van der Pol — the weak-mu perturbation laws (amplitude
//!   -> 2, `T = 2 pi (1 + mu^2/16)`) and the mu = 1 orbit against
//!   TWO INDEPENDENT published sources:
//!   Amore, "Computing the solutions of the van der Pol equation to
//!   arbitrary precision" (arXiv:2111.12198): T = 6.663286859323130,
//!   amplitude = 2.008619860874843; Gasull, Giacomini & Grau,
//!   "Proving the existence of numerically detected planar limit
//!   cycles" (arXiv:1602.00113): T ~ 6.6632866.
//! - ob-003: pseudo-arclength continuation TRAVERSES both folds of
//!   the forced-Duffing response S-curve; fold frequencies match the
//!   independent first-harmonic scalar law.
//! - ob-004: shooting-vs-HB cross-validation on the same van der Pol
//!   orbit (two independent methods agreeing is the artifact
//!   detector) + Floquet stability from the monodromy.
//! - ob-005: refusals by name, including the TorusSuspected named
//!   no-claim on a constructed two-frequency problem.
//! - ob-006: bitwise determinism.

use fs_orbit::{
    ContinuableProblem, ContinuationBudget, HbAnchor, HbBudget, OrbitError, OrbitProblem,
    ShootBudget, continue_branch, solve_hb, solve_shooting,
};

const TAU: f64 = core::f64::consts::TAU;

/// Duffing: `x'' + 2 zeta x' + x + eps x^3 = f cos(omega t)`.
struct Duffing {
    zeta: f64,
    eps: f64,
    force: f64,
    omega: f64,
}

impl OrbitProblem for Duffing {
    fn dim(&self) -> usize {
        2
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[1];
        out[1] = self.force * fs_math::det::cos(self.omega * t)
            - 2.0 * self.zeta * x[1]
            - x[0]
            - self.eps * x[0] * x[0] * x[0];
    }
    fn autonomous(&self) -> bool {
        self.force == 0.0
    }
}

impl ContinuableProblem for Duffing {
    fn set_parameter(&mut self, lambda: f64) {
        self.omega = lambda;
    }
}

/// Van der Pol: `x'' - mu (1 - x^2) x' + x = 0`.
struct VanDerPol {
    mu: f64,
}

impl OrbitProblem for VanDerPol {
    fn dim(&self) -> usize {
        2
    }
    fn island(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[1];
        out[1] = self.mu * (1.0 - x[0] * x[0]) * x[1] - x[0];
    }
    fn autonomous(&self) -> bool {
        true
    }
}

#[test]
fn ob_001_duffing_backbone_matches_the_pinned_law() {
    let eps = 0.1;
    let problem = Duffing {
        zeta: 0.0,
        eps,
        force: 0.0,
        omega: 1.0,
    };
    let budget = HbBudget::default();
    let mut worst = 0.0f64;
    for amplitude in [0.2f64, 0.5, 0.8] {
        let orbit = solve_hb(
            &problem,
            HbAnchor::Backbone {
                amplitude,
                omega_guess: 1.0,
            },
            amplitude,
            &budget,
        )
        .expect("backbone point");
        let pinned = fs_nlmodal::duffing_backbone(1.0, eps, amplitude);
        // The pinned law is FIRST order; the HB answer carries the
        // higher orders, so the band is the eps^2 a^4 term's scale.
        let band = 0.02 * eps * eps * amplitude.powi(4) + 1.0e-6;
        let dev = (orbit.omega - pinned).abs();
        worst = worst.max(dev / pinned);
        assert!(
            dev < band.max(2.0e-4),
            "a={amplitude}: HB omega {:.8} vs pinned {pinned:.8}",
            orbit.omega
        );
        println!(
            "{{\"suite\":\"fs-orbit\",\"case\":\"ob-001-backbone\",\"amplitude\":{amplitude},\
             \"omega_hb\":{:.8},\"omega_pinned\":{pinned:.8},\"newton_iters\":{},\
             \"residual\":{:.3e}}}",
            orbit.omega,
            orbit.residual_trace.len(),
            orbit.residual
        );
    }
    assert!(worst < 5.0e-4, "worst backbone deviation {worst:.2e}");
}

#[test]
fn ob_002_van_der_pol_matches_published_values() {
    // Weak mu: the textbook perturbation laws.
    let weak = VanDerPol { mu: 0.2 };
    let budget = HbBudget {
        harmonics: 12,
        ..HbBudget::default()
    };
    let orbit = solve_hb(
        &weak,
        HbAnchor::Autonomous { omega_guess: 1.0 },
        2.0,
        &budget,
    )
    .expect("weak vdP");
    let period = TAU / orbit.omega;
    let t_law = TAU * (1.0 + 0.2f64 * 0.2 / 16.0);
    assert!(
        (period / t_law - 1.0).abs() < 1.0e-3,
        "weak-mu period {period:.6} vs law {t_law:.6}"
    );
    let amp = orbit.peak(0);
    assert!((amp - 2.0).abs() < 0.02, "weak-mu amplitude {amp:.4}");
    // mu = 1 against the two published sources.
    let strong = VanDerPol { mu: 1.0 };
    let budget = HbBudget {
        harmonics: 17,
        max_newton: 60,
        ..HbBudget::default()
    };
    let orbit = solve_hb(
        &strong,
        HbAnchor::Autonomous { omega_guess: 0.94 },
        2.0,
        &budget,
    )
    .expect("mu=1 vdP");
    let period = TAU / orbit.omega;
    let t_amore = 6.663_286_859_323_130;
    let t_gasull = 6.663_286_6;
    assert!(
        (period / t_amore - 1.0).abs() < 1.0e-4,
        "mu=1 period {period:.8} vs Amore {t_amore}"
    );
    assert!(
        (period / t_gasull - 1.0).abs() < 1.0e-4,
        "mu=1 period {period:.8} vs Gasull-Giacomini-Grau {t_gasull}"
    );
    let amp = orbit.peak(0);
    let a_amore = 2.008_619_860_874_843;
    assert!(
        (amp / a_amore - 1.0).abs() < 2.0e-3,
        "mu=1 amplitude {amp:.6} vs Amore {a_amore}"
    );
    println!(
        "{{\"suite\":\"fs-orbit\",\"case\":\"ob-002-vdp\",\"mu\":1.0,\"period\":{period:.9},\
         \"published\":[{t_amore},{t_gasull}],\"amplitude\":{amp:.6},\
         \"published_amplitude\":{a_amore},\"harmonics\":17,\"residual\":{:.3e}}}",
        orbit.residual
    );
}

/// First-harmonic scalar law of the forced Duffing:
/// `F^2 = r^2 ((1 - w^2 + 3 eps r^2 / 4)^2 + (2 zeta w)^2)`. For a
/// response amplitude `r` the two `w^2` roots are analytic; folds are
/// the extrema of `w(r)` — computed here by a fine scan, giving an
/// INDEPENDENT fold-location oracle.
fn duffing_fold_omegas(zeta: f64, eps: f64, force: f64) -> (f64, f64) {
    let mut upper_fold = f64::NAN;
    let mut lower_fold = f64::NAN;
    let mut prev_w_hi = f64::NAN;
    let mut rising = true;
    let mut r = 0.02;
    let mut prev_r_ok = false;
    while r < 6.0 {
        let a = 1.0 + 0.75 * eps * r * r;
        // w^2 roots of (a - w^2)^2 + 4 zeta^2 w^2 = (F/r)^2.
        let c = force / r;
        // (w^2)^2 - 2(a - 2 zeta^2) w^2 + a^2 - c^2 = 0.
        let b_half = a - 2.0 * zeta * zeta;
        let disc = b_half * b_half - (a * a - c * c);
        if disc >= 0.0 {
            let sq = disc.sqrt();
            let w2_hi = b_half + sq;
            if w2_hi > 0.0 {
                let w_hi = w2_hi.sqrt();
                if prev_r_ok {
                    let now_rising = w_hi > prev_w_hi;
                    if rising && !now_rising {
                        upper_fold = prev_w_hi;
                    }
                    if !rising && now_rising {
                        lower_fold = prev_w_hi;
                    }
                    rising = now_rising;
                }
                prev_w_hi = w_hi;
                prev_r_ok = true;
            }
        }
        r += 1.0e-4;
    }
    (lower_fold, upper_fold)
}

#[test]
fn ob_003_continuation_traverses_the_duffing_folds() {
    let (zeta, eps, force) = (0.05, 0.15, 0.3);
    let mut problem = Duffing {
        zeta,
        eps,
        force,
        omega: 0.6,
    };
    let budget = ContinuationBudget {
        max_steps: 400,
        initial_step: 0.02,
        min_step: 1.0e-5,
        hb: HbBudget {
            harmonics: 7,
            max_newton: 25,
            tolerance: 1.0e-9,
        },
    };
    let points = continue_branch(&mut problem, |lam| lam, 0.6, 1.8, 0.35, &budget)
        .expect("branch traversal");
    assert!(points.len() > 10, "branch points: {}", points.len());
    // Fold traversal: the parameter sequence reverses direction at
    // least twice (up the low branch, back through the upper fold,
    // forward again through the lower fold).
    let mut reversals = 0usize;
    let mut fold_lambdas = Vec::new();
    for w in points.windows(3) {
        let d1 = w[1].lambda - w[0].lambda;
        let d2 = w[2].lambda - w[1].lambda;
        if d1 * d2 < 0.0 {
            reversals += 1;
            fold_lambdas.push(w[1].lambda);
        }
    }
    assert!(
        reversals >= 2,
        "expected two fold reversals, saw {reversals}"
    );
    // Independent scalar-law fold locations (first-harmonic; the HB
    // branch carries harmonics, so the band is authored at 3%).
    let (fold_lo, fold_hi) = duffing_fold_omegas(zeta, eps, force);
    let seen_hi = fold_lambdas[0];
    let seen_lo = fold_lambdas[1];
    assert!(
        (seen_hi / fold_hi - 1.0).abs() < 0.03,
        "upper fold {seen_hi:.4} vs scalar law {fold_hi:.4}"
    );
    assert!(
        (seen_lo / fold_lo - 1.0).abs() < 0.03,
        "lower fold {seen_lo:.4} vs scalar law {fold_lo:.4}"
    );
    // The branch visits both response levels around the fold pair.
    let max_amp = points
        .iter()
        .map(|p| p.orbit.first_harmonic_amplitude(0))
        .fold(0.0f64, f64::max);
    let min_amp = points
        .iter()
        .map(|p| p.orbit.first_harmonic_amplitude(0))
        .fold(f64::INFINITY, f64::min);
    assert!(
        max_amp > 3.0 * min_amp,
        "S-curve span {min_amp:.3}..{max_amp:.3}"
    );
    println!(
        "{{\"suite\":\"fs-orbit\",\"case\":\"ob-003-folds\",\"points\":{},\
         \"reversals\":{reversals},\"fold_hi_hb\":{seen_hi:.4},\"fold_hi_law\":{fold_hi:.4},\
         \"fold_lo_hb\":{seen_lo:.4},\"fold_lo_law\":{fold_lo:.4},\
         \"amp_span\":[{min_amp:.3},{max_amp:.3}],\
         \"step_sizes\":[{:.4},{:.4}]}}",
        points.len(),
        points[1].step,
        points[points.len() - 1].step
    );
}

#[test]
fn ob_004_shooting_cross_validates_hb_and_classifies_stability() {
    let problem = VanDerPol { mu: 1.0 };
    let hb = solve_hb(
        &problem,
        HbAnchor::Autonomous { omega_guess: 0.94 },
        2.0,
        &HbBudget {
            harmonics: 17,
            max_newton: 60,
            ..HbBudget::default()
        },
    )
    .expect("HB orbit");
    // Seed shooting FROM the HB orbit (phase 0 point).
    let seed = vec![hb.sample(0, 0.0), hb.sample(1, 0.0)];
    let shot = solve_shooting(&problem, &seed, TAU / hb.omega, &ShootBudget::default())
        .expect("shooting orbit");
    let t_hb = TAU / hb.omega;
    let rel = (shot.period / t_hb - 1.0).abs();
    assert!(
        rel < 5.0e-5,
        "shooting period {:.9} vs HB {:.9} (rel {rel:.2e})",
        shot.period,
        t_hb
    );
    // Floquet: one trivial multiplier at +1, the other inside the
    // unit circle (the vdP cycle is stable).
    let mut trivial = false;
    let mut stable = false;
    for &(re, im) in &shot.multipliers {
        let mag = (re * re + im * im).sqrt();
        if (re - 1.0).abs() < 1.0e-2 && im.abs() < 1.0e-2 {
            trivial = true;
        } else if mag < 0.9 {
            stable = true;
        }
    }
    assert!(trivial && stable, "multipliers {:?}", shot.multipliers);
    println!(
        "{{\"suite\":\"fs-orbit\",\"case\":\"ob-004-cross\",\"t_hb\":{t_hb:.9},\
         \"t_shoot\":{:.9},\"rel\":{rel:.2e},\"multipliers\":{:?},\
         \"shoot_iters\":{}}}",
        shot.period,
        shot.multipliers,
        shot.residual_trace.len()
    );
}

/// Two independent rotations at incommensurate frequencies: a
/// constructed monodromy with a complex unit-circle pair — the
/// TorusSuspected detector's fixture.
struct TwoCenters;

impl OrbitProblem for TwoCenters {
    fn dim(&self) -> usize {
        4
    }
    fn island(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[1];
        out[1] = -x[0];
        let w2 = core::f64::consts::SQRT_2;
        out[2] = w2 * x[3];
        out[3] = -w2 * x[2];
    }
    fn autonomous(&self) -> bool {
        false
    }
}

#[test]
fn ob_005_refusals_fire_by_name() {
    // Bad parameters.
    let vdp = VanDerPol { mu: 1.0 };
    assert!(matches!(
        solve_hb(
            &vdp,
            HbAnchor::Autonomous { omega_guess: -1.0 },
            2.0,
            &HbBudget::default()
        ),
        Err(OrbitError::BadParameter { .. })
    ));
    // Budget exhaustion is a typed stall with a residual trace.
    match solve_hb(
        &vdp,
        HbAnchor::Autonomous { omega_guess: 0.94 },
        2.0,
        &HbBudget {
            harmonics: 17,
            max_newton: 2,
            tolerance: 1.0e-14,
        },
    ) {
        Err(OrbitError::NewtonStalled { trace, .. }) => {
            assert_eq!(trace.len(), 2, "trace discloses per-iteration residuals");
        }
        other => panic!("expected NewtonStalled, got {other:?}"),
    }
    // The named quasi-periodic no-claim: the two-center fixture
    // converges as a periodic orbit of the first pair (the second
    // pair seeded at rest) and its monodromy carries the
    // incommensurate rotation — TorusSuspected by name.
    let torus = solve_shooting(
        &TwoCenters,
        &[1.0, 0.0, 0.0, 0.0],
        TAU,
        &ShootBudget::default(),
    );
    assert!(
        matches!(torus, Err(OrbitError::TorusSuspected { .. })),
        "expected TorusSuspected, got {torus:?}"
    );
    println!("{{\"suite\":\"fs-orbit\",\"case\":\"ob-005-refusals\",\"verdict\":\"pass\"}}");
}

#[test]
fn ob_006_bitwise_determinism() {
    let problem = VanDerPol { mu: 0.7 };
    let budget = HbBudget {
        harmonics: 11,
        ..HbBudget::default()
    };
    let a = solve_hb(
        &problem,
        HbAnchor::Autonomous { omega_guess: 0.97 },
        2.0,
        &budget,
    )
    .expect("first");
    let b = solve_hb(
        &problem,
        HbAnchor::Autonomous { omega_guess: 0.97 },
        2.0,
        &budget,
    )
    .expect("second");
    assert!(a.omega.to_bits() == b.omega.to_bits());
    for (ca, cb) in a.coeffs.iter().zip(&b.coeffs) {
        for (x, y) in ca.iter().zip(cb) {
            assert!(x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits());
        }
    }
    assert_eq!(a.residual_trace.len(), b.residual_trace.len());
    println!("{{\"suite\":\"fs-orbit\",\"case\":\"ob-006-determinism\",\"verdict\":\"pass\"}}");
}

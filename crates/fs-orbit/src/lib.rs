//! fs-orbit — Periodic-orbit service: harmonic balance, shooting,
//! continuation (music bead `frankensim-music-v8-root-3ez8g.11.1`;
//! the first implementation of the fs-vmanifest I09 slot).
//!
//! Layer: L2. See CONTRACT.md for invariants, error model, determinism
//! class, cancellation behavior, and no-claim boundaries.
//!
//! THE PROBLEM SHAPE is island-plus-linear-port, never
//! instrument-shaped: a state vector whose dynamics carry a nonlinear
//! part evaluated pointwise in TIME (the island — a reed law, a device
//! card, a cubic spring) and a linear part specified per HARMONIC in
//! the FREQUENCY domain (the port — `s I` for a plain ODE, a TMM
//! impedance `Z(n omega)` for an acoustic load). Harmonic balance
//! closes the loop by alternating-frequency/time (AFT) evaluation at a
//! fixed truncation `N` (the truncation is DISCLOSED structure,
//! X-Struct). Shooting is the independent second method — the
//! artifact detector, exactly like VF-vs-Loewner — and its monodromy
//! matrix yields Floquet multipliers through the fs-la complex
//! eigensolver. Pseudo-arclength continuation traverses folds;
//! thresholds and slot maps ARE fold structures, so continuation is
//! the product, not a luxury.
//!
//! Determinism: fixed sample grids, fixed iteration caps, direct
//! trigonometric AFT sums through `fs_math::det`, no RNG, canonical
//! unknown ordering (DC, then per-harmonic Re/Im pairs, state-major).
//! Quasi-periodic drift (torus birth) is a NAMED refusal in v1, never
//! a wrong answer.

use fs_la::eigen_complex::{EigFailure, eig};
use fs_la::factor::lu;
use fs_math::c64::C64;
use fs_math::det;

const TAU: f64 = core::f64::consts::TAU;

/// Typed refusals.
#[derive(Debug, Clone, PartialEq)]
pub enum OrbitError {
    /// A parameter is non-finite or out of range.
    BadParameter {
        /// Which parameter refused.
        what: &'static str,
    },
    /// Newton exhausted its iteration budget without meeting the
    /// tolerance. The residual trace is disclosed for diagnosis.
    NewtonStalled {
        /// Final residual norm.
        residual: f64,
        /// Iterations spent.
        iterations: usize,
        /// Per-iteration residual norms.
        trace: Vec<f64>,
    },
    /// The HB/shooting Jacobian became singular.
    SingularJacobian,
    /// Continuation exhausted its step budget or its minimum step.
    ContinuationExhausted {
        /// Steps completed before exhaustion.
        steps: usize,
    },
    /// A non-trivial Floquet multiplier pair sits on the unit circle
    /// with nonzero angle: a torus (quasi-periodic drift) is suspected.
    /// v1 NAMES this and refuses — no quasi-periodic claim exists.
    TorusSuspected {
        /// The offending multiplier.
        multiplier: (f64, f64),
    },
    /// Newton converged onto the trivial equilibrium instead of an
    /// orbit — refused by name (a zero solution satisfies every
    /// autonomous balance; claiming it as a cycle would be silent
    /// garbage).
    TrivialCollapse,
    /// Eigensolver refusal from fs-la, forwarded.
    Eigen(EigFailure),
}

impl core::fmt::Display for OrbitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            OrbitError::BadParameter { what } => write!(f, "FS-ORBIT-PARAM: {what}"),
            OrbitError::NewtonStalled {
                residual,
                iterations,
                ..
            } => write!(
                f,
                "FS-ORBIT-STALL: residual {residual:.3e} after {iterations} iterations"
            ),
            OrbitError::SingularJacobian => write!(f, "FS-ORBIT-SINGULAR: Jacobian factorization"),
            OrbitError::ContinuationExhausted { steps } => {
                write!(f, "FS-ORBIT-CONTINUATION: exhausted after {steps} steps")
            }
            OrbitError::TorusSuspected { multiplier } => write!(
                f,
                "FS-ORBIT-TORUS: non-trivial unit-circle multiplier ({}, {}) — \
                 quasi-periodic drift is a named no-claim in v1",
                multiplier.0, multiplier.1
            ),
            OrbitError::TrivialCollapse => {
                write!(
                    f,
                    "FS-ORBIT-TRIVIAL: converged to the equilibrium, not an orbit"
                )
            }
            OrbitError::Eigen(e) => write!(f, "FS-ORBIT-EIG: {e:?}"),
        }
    }
}

impl std::error::Error for OrbitError {}

/// The island-plus-linear-port problem shape.
pub trait OrbitProblem {
    /// State dimension `d`.
    fn dim(&self) -> usize;

    /// The NONLINEAR island: the full right-hand side `g(t, x)` of
    /// `M x' = g` evaluated pointwise in time (AFT samples this).
    /// For a plain ODE this is the whole vector field.
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]);

    /// The LINEAR port at Laplace argument `s = i n omega`, row-major
    /// `d x d`: the harmonic-domain operator applied to `X_n` on the
    /// left of the balance `P(s) X_n = G_n`. Default `s I` (a plain
    /// first-order ODE). A TMM impedance load overrides this.
    fn port(&self, s: C64) -> Vec<C64> {
        let d = self.dim();
        let mut m = vec![C64::new(0.0, 0.0); d * d];
        for i in 0..d {
            m[i * d + i] = s;
        }
        m
    }

    /// Whether the problem is autonomous (the period is then an
    /// unknown and a phase anchor closes the system).
    fn autonomous(&self) -> bool;
}

/// How the phase/frequency degrees of freedom are closed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HbAnchor {
    /// Forced problem at a known angular frequency.
    Forced {
        /// The forcing angular frequency [rad/s].
        omega: f64,
    },
    /// Autonomous orbit: `omega` is an unknown; the phase is anchored
    /// by `Im X_1[0] = 0`.
    Autonomous {
        /// Initial guess for the angular frequency.
        omega_guess: f64,
    },
    /// Conservative BACKBONE point: amplitude of `X_1[0]` is pinned
    /// (`Re X_1[0] = amplitude`, `Im X_1[0] = 0`), `omega` is the
    /// unknown, and the redundant phase equation is dropped.
    Backbone {
        /// Pinned first-harmonic amplitude of state 0.
        amplitude: f64,
        /// Initial guess for the angular frequency.
        omega_guess: f64,
    },
}

/// Fixed solver budgets (deterministic; exhaustion is a typed refusal).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HbBudget {
    /// Harmonic truncation `N` (DC + N complex harmonics).
    pub harmonics: usize,
    /// Newton iteration cap.
    pub max_newton: usize,
    /// Relative residual tolerance.
    pub tolerance: f64,
}

impl Default for HbBudget {
    fn default() -> Self {
        HbBudget {
            harmonics: 9,
            max_newton: 40,
            tolerance: 1.0e-10,
        }
    }
}

/// A converged harmonic-balance orbit.
#[derive(Debug, Clone)]
pub struct HbOrbit {
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Coefficients: `coeffs[k][n]` is `X_n` of state `k` (`n = 0` is
    /// DC, real part only meaningful there).
    pub coeffs: Vec<Vec<C64>>,
    /// Per-iteration residual norms (the disclosed log).
    pub residual_trace: Vec<f64>,
    /// Final residual norm.
    pub residual: f64,
}

impl HbOrbit {
    /// Synthesize the state at phase `theta` in `[0, tau)`.
    #[must_use]
    pub fn sample(&self, k: usize, theta: f64) -> f64 {
        let mut x = self.coeffs[k][0].re;
        for (n, c) in self.coeffs[k].iter().enumerate().skip(1) {
            let a = n as f64 * theta;
            x += 2.0 * (c.re * det::cos(a) - c.im * det::sin(a));
        }
        x
    }

    /// First-harmonic amplitude of state `k` (`2 |X_1|`).
    #[must_use]
    pub fn first_harmonic_amplitude(&self, k: usize) -> f64 {
        2.0 * self.coeffs[k][1].abs()
    }

    /// Peak of state `k` over a fine phase grid (256 points).
    #[must_use]
    pub fn peak(&self, k: usize) -> f64 {
        let mut peak = f64::NEG_INFINITY;
        for j in 0..256i32 {
            peak = peak.max(self.sample(k, TAU * f64::from(j) / 256.0));
        }
        peak
    }
}

/// Internal: canonical packing of unknowns.
///
/// `u = [X_0[0..d] (Re), {Re X_n[0..d], Im X_n[0..d]}_{n=1..N}, omega?]`
struct Packing {
    d: usize,
    n_harm: usize,
    omega_slot: Option<usize>,
}

impl Packing {
    fn coeff_len(&self) -> usize {
        self.d * (2 * self.n_harm + 1)
    }
    fn len(&self) -> usize {
        self.coeff_len() + usize::from(self.omega_slot.is_some())
    }
    fn re_index(&self, n: usize, k: usize) -> usize {
        if n == 0 {
            k
        } else {
            self.d + (n - 1) * 2 * self.d + k
        }
    }
    fn im_index(&self, n: usize, k: usize) -> usize {
        debug_assert!(n >= 1);
        self.d + (n - 1) * 2 * self.d + self.d + k
    }
    fn unpack(&self, u: &[f64]) -> (Vec<Vec<C64>>, f64) {
        let mut coeffs = vec![vec![C64::new(0.0, 0.0); self.n_harm + 1]; self.d];
        for k in 0..self.d {
            coeffs[k][0] = C64::new(u[self.re_index(0, k)], 0.0);
            for n in 1..=self.n_harm {
                coeffs[k][n] = C64::new(u[self.re_index(n, k)], u[self.im_index(n, k)]);
            }
        }
        let omega = self.omega_slot.map_or(f64::NAN, |s| u[s]);
        (coeffs, omega)
    }
}

/// The HB residual: for each harmonic `n`, `P(i n omega) X_n - G_n`,
/// with `G` the AFT transform of the island over `m` samples.
fn hb_residual<P: OrbitProblem>(
    problem: &P,
    pack: &Packing,
    u: &[f64],
    omega_fixed: f64,
    out: &mut [f64],
) {
    let d = pack.d;
    let n_harm = pack.n_harm;
    let omega = pack.omega_slot.map_or(omega_fixed, |s| u[s]);
    let m = (4 * n_harm + 4).next_power_of_two();
    // Synthesize samples.
    let mut xs = vec![0.0f64; m * d];
    for j in 0..m {
        let theta = TAU * j as f64 / m as f64;
        for k in 0..d {
            let mut x = u[pack.re_index(0, k)];
            for n in 1..=n_harm {
                let a = n as f64 * theta;
                x += 2.0
                    * (u[pack.re_index(n, k)] * det::cos(a) - u[pack.im_index(n, k)] * det::sin(a));
            }
            xs[j * d + k] = x;
        }
    }
    // Island evaluation per sample.
    let mut gs = vec![0.0f64; m * d];
    let mut buf = vec![0.0f64; d];
    let period = if omega > 0.0 { TAU / omega } else { f64::NAN };
    for j in 0..m {
        let t = period * j as f64 / m as f64;
        problem.island(t, &xs[j * d..(j + 1) * d], &mut buf);
        gs[j * d..(j + 1) * d].copy_from_slice(&buf);
    }
    // Analysis: G_n = (1/m) sum_j g_j e^{-i n theta_j}.
    // Balance rows in the canonical packing order.
    for k in 0..d {
        // DC row.
        let mut g0 = 0.0;
        for j in 0..m {
            g0 += gs[j * d + k];
        }
        g0 /= m as f64;
        // P(0) X_0 - G_0 (real part; the DC port row).
        let p0 = problem.port(C64::new(0.0, 0.0));
        let mut px = C64::new(0.0, 0.0);
        for kk in 0..d {
            px = px + p0[k * d + kk].scale(u[pack.re_index(0, kk)]);
        }
        out[pack.re_index(0, k)] = px.re - g0;
    }
    for n in 1..=n_harm {
        let s = C64::new(0.0, n as f64 * omega);
        let p = problem.port(s);
        for k in 0..d {
            let mut gn = C64::new(0.0, 0.0);
            for j in 0..m {
                let a = n as f64 * TAU * j as f64 / m as f64;
                gn = gn + C64::new(det::cos(a), -det::sin(a)).scale(gs[j * d + k]);
            }
            gn = gn.scale(1.0 / m as f64);
            let mut px = C64::new(0.0, 0.0);
            for kk in 0..d {
                let xkk = C64::new(u[pack.re_index(n, kk)], u[pack.im_index(n, kk)]);
                px = px + p[k * d + kk] * xkk;
            }
            out[pack.re_index(n, k)] = px.re - gn.re;
            out[pack.im_index(n, k)] = px.im - gn.im;
        }
    }
}

fn norm(v: &[f64]) -> f64 {
    det::sqrt(v.iter().map(|x| x * x).sum::<f64>())
}

/// Solve a harmonic-balance orbit.
///
/// `guess` seeds the first harmonic of state 0 (`Re X_1[0]`); richer
/// seeding goes through [`solve_hb_seeded`].
///
/// # Errors
/// Typed [`OrbitError`] (budgets, singular Jacobian, bad parameters).
pub fn solve_hb<P: OrbitProblem>(
    problem: &P,
    anchor: HbAnchor,
    guess_amplitude: f64,
    budget: &HbBudget,
) -> Result<HbOrbit, OrbitError> {
    let d = problem.dim();
    let pack = Packing {
        d,
        n_harm: budget.harmonics,
        omega_slot: match anchor {
            HbAnchor::Forced { .. } => None,
            _ => Some(d * (2 * budget.harmonics + 1)),
        },
    };
    let mut u = vec![0.0f64; pack.len()];
    u[pack.re_index(1, 0)] = 0.5 * guess_amplitude;
    // Kinematically consistent seed for first-order forms of
    // second-order problems: state 1 as the derivative of state 0
    // (`x = a cos -> v = -a w sin`). Harmless for other shapes (it is
    // only a seed).
    if d >= 2 {
        let w0 = match anchor {
            HbAnchor::Forced { omega } => omega,
            HbAnchor::Autonomous { omega_guess } | HbAnchor::Backbone { omega_guess, .. } => {
                omega_guess
            }
        };
        if w0.is_finite() && w0 > 0.0 {
            u[pack.im_index(1, 1)] = 0.5 * guess_amplitude * w0;
        }
    }
    match anchor {
        HbAnchor::Forced { omega } => {
            if !(omega > 0.0 && omega.is_finite()) {
                return Err(OrbitError::BadParameter {
                    what: "forcing frequency must be positive and finite",
                });
            }
        }
        HbAnchor::Autonomous { omega_guess } | HbAnchor::Backbone { omega_guess, .. } => {
            if !(omega_guess > 0.0 && omega_guess.is_finite()) {
                return Err(OrbitError::BadParameter {
                    what: "frequency guess must be positive and finite",
                });
            }
            u[pack.omega_slot.expect("slot")] = omega_guess;
        }
    }
    if let HbAnchor::Backbone { amplitude, .. } = anchor {
        u[pack.re_index(1, 0)] = 0.5 * amplitude;
    }
    solve_hb_seeded(problem, anchor, u, budget)
}

/// [`solve_hb`] from a full packed seed (continuation uses this).
///
/// # Errors
/// Typed [`OrbitError`].
#[allow(clippy::too_many_lines)] // one coherent masked-Newton stage
pub fn solve_hb_seeded<P: OrbitProblem>(
    problem: &P,
    anchor: HbAnchor,
    mut u: Vec<f64>,
    budget: &HbBudget,
) -> Result<HbOrbit, OrbitError> {
    let d = problem.dim();
    if d == 0 {
        return Err(OrbitError::BadParameter {
            what: "state dimension must be positive",
        });
    }
    let pack = Packing {
        d,
        n_harm: budget.harmonics,
        omega_slot: match anchor {
            HbAnchor::Forced { .. } => None,
            _ => Some(d * (2 * budget.harmonics + 1)),
        },
    };
    if u.len() != pack.len() {
        return Err(OrbitError::BadParameter {
            what: "seed length does not match the packing",
        });
    }
    let omega_fixed = match anchor {
        HbAnchor::Forced { omega } => omega,
        _ => f64::NAN,
    };
    // Masks: fixed unknown indices + dropped equation indices.
    // Backbone closes with an APPENDED amplitude-norm equation
    // (|X_1[0]| = a/2) rather than pinned components: every balance
    // row stays enforced and the phase stays free (the conservative
    // family's phase nullspace is handled by the regularized step; a
    // converged point satisfies ALL equations at the pinned
    // amplitude, so omega is phase-independent).
    let amplitude_row: Option<f64> = match anchor {
        HbAnchor::Backbone { amplitude, .. } => Some(amplitude),
        _ => None,
    };
    let (fixed, dropped): (Vec<usize>, Vec<usize>) = match anchor {
        HbAnchor::Forced { .. } | HbAnchor::Backbone { .. } => (vec![], vec![]),
        HbAnchor::Autonomous { .. } => (vec![pack.im_index(1, 0)], vec![]),
    };
    if let HbAnchor::Autonomous { .. } = anchor {
        u[pack.im_index(1, 0)] = 0.0;
    }
    let free: Vec<usize> = (0..pack.len()).filter(|i| !fixed.contains(i)).collect();
    let kept: Vec<usize> = (0..pack.coeff_len())
        .filter(|i| !dropped.contains(i))
        .collect();
    let extra = usize::from(amplitude_row.is_some());
    if free.len() != kept.len() + extra {
        return Err(OrbitError::BadParameter {
            what: "anchor produced a non-square system",
        });
    }
    let nn = free.len();
    let mut res = vec![0.0f64; pack.coeff_len()];
    let mut trace = Vec::new();
    let amp_residual = |u: &[f64]| -> f64 {
        amplitude_row.map_or(0.0, |a| {
            let re = u[pack.re_index(1, 0)];
            let im = u[pack.im_index(1, 0)];
            re * re + im * im - 0.25 * a * a
        })
    };
    for iter in 0..budget.max_newton {
        hb_residual(problem, &pack, &u, omega_fixed, &mut res);
        let mut r: Vec<f64> = kept.iter().map(|&i| res[i]).collect();
        if amplitude_row.is_some() {
            r.push(amp_residual(&u));
        }
        let rn = norm(&r);
        trace.push(rn);
        // Scale over the COEFFICIENTS only: the omega slot is in
        // rad/s and would inflate the relative gate (measured: a
        // sub-threshold reed "converged" at 7.7e-7 because omega
        // ~1600 stretched the scale).
        let scale = norm(&u[..pack.coeff_len()]).max(1.0);
        if rn < budget.tolerance * scale {
            // Trivial-collapse guard: an all-zero harmonic content is
            // the equilibrium, not an orbit.
            let first_harm: f64 = (0..d)
                .map(|k| {
                    let re = u[pack.re_index(1, k)];
                    let im = u[pack.im_index(1, k)];
                    det::sqrt(re * re + im * im)
                })
                .sum();
            if first_harm < 1.0e-9 {
                return Err(OrbitError::TrivialCollapse);
            }
            let (coeffs, omega_u) = pack.unpack(&u);
            let omega = if pack.omega_slot.is_some() {
                omega_u
            } else {
                omega_fixed
            };
            return Ok(HbOrbit {
                omega,
                coeffs,
                residual_trace: trace,
                residual: rn,
            });
        }
        // Forward-difference Jacobian over the free unknowns.
        let mut jac = vec![0.0f64; nn * nn];
        let mut res2 = vec![0.0f64; pack.coeff_len()];
        let amp0 = amp_residual(&u);
        for (col, &ui) in free.iter().enumerate() {
            let h = 1.0e-7 * u[ui].abs().max(1.0);
            let saved = u[ui];
            u[ui] = saved + h;
            hb_residual(problem, &pack, &u, omega_fixed, &mut res2);
            if amplitude_row.is_some() {
                jac[kept.len() * nn + col] = (amp_residual(&u) - amp0) / h;
            }
            u[ui] = saved;
            for (row, &ri) in kept.iter().enumerate() {
                jac[row * nn + col] = (res2[ri] - res[ri]) / h;
            }
        }
        // Conservative families carry a parity nullspace; a fixed
        // Tikhonov jitter keeps the step deterministic and the
        // converged residual is still the only accepted evidence.
        let factorized = if let Ok(f) = lu(&jac, nn) {
            f
        } else {
            let scale_j = jac.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let delta = 1.0e-10 * scale_j.max(1.0);
            for i in 0..nn {
                jac[i * nn + i] += delta;
            }
            match lu(&jac, nn) {
                Ok(f) => f,
                Err(_) => return Err(OrbitError::SingularJacobian),
            }
        };
        let mut step: Vec<f64> = r.clone();
        factorized.solve(&mut step);
        // Deterministic backtracking: halve the step (up to 8 times)
        // until the residual does not grow — the undamped Newton
        // overshoots island problems into the equilibrium basin
        // (measured on the reed gate).
        let mut best_alpha = 1.0f64;
        let mut best_norm = f64::INFINITY;
        let mut trial = vec![0.0f64; pack.coeff_len()];
        let mut alpha = 1.0f64;
        for _ in 0..8 {
            let mut u_try = u.clone();
            for (col, &ui) in free.iter().enumerate() {
                u_try[ui] -= alpha * step[col];
            }
            hb_residual(problem, &pack, &u_try, omega_fixed, &mut trial);
            let mut tn: f64 = kept.iter().map(|&i| trial[i] * trial[i]).sum();
            if amplitude_row.is_some() {
                let ar = amp_residual(&u_try);
                tn += ar * ar;
            }
            let tn = det::sqrt(tn);
            if tn < best_norm {
                best_norm = tn;
                best_alpha = alpha;
            }
            if tn < rn {
                break;
            }
            alpha *= 0.5;
        }
        for (col, &ui) in free.iter().enumerate() {
            u[ui] -= best_alpha * step[col];
        }
        let _ = iter;
    }
    Err(OrbitError::NewtonStalled {
        residual: *trace.last().unwrap_or(&f64::NAN),
        iterations: budget.max_newton,
        trace,
    })
}

// ---------------------------------------------------------------------
// Shooting + Floquet.

/// Shooting budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShootBudget {
    /// RK4 steps per period (fixed, deterministic).
    pub steps_per_period: usize,
    /// Newton iteration cap.
    pub max_newton: usize,
    /// Relative residual tolerance.
    pub tolerance: f64,
}

impl Default for ShootBudget {
    fn default() -> Self {
        ShootBudget {
            steps_per_period: 4096,
            max_newton: 30,
            tolerance: 1.0e-10,
        }
    }
}

/// A converged shooting orbit.
#[derive(Debug, Clone)]
pub struct ShootOrbit {
    /// Period [s].
    pub period: f64,
    /// A point on the orbit.
    pub x0: Vec<f64>,
    /// Floquet multipliers of the monodromy matrix.
    pub multipliers: Vec<(f64, f64)>,
    /// Per-iteration residual norms.
    pub residual_trace: Vec<f64>,
}

/// Fixed-step RK4 flow over `[0, t_end]`.
fn flow<P: OrbitProblem>(problem: &P, x0: &[f64], t_end: f64, steps: usize) -> Vec<f64> {
    let d = problem.dim();
    let h = t_end / steps as f64;
    let mut x = x0.to_vec();
    let (mut k1, mut k2, mut k3, mut k4) = (vec![0.0; d], vec![0.0; d], vec![0.0; d], vec![0.0; d]);
    let mut tmp = vec![0.0; d];
    for j in 0..steps {
        let t = h * j as f64;
        problem.island(t, &x, &mut k1);
        for i in 0..d {
            tmp[i] = 0.5f64.mul_add(h * k1[i], x[i]);
        }
        problem.island(t + 0.5 * h, &tmp, &mut k2);
        for i in 0..d {
            tmp[i] = 0.5f64.mul_add(h * k2[i], x[i]);
        }
        problem.island(t + 0.5 * h, &tmp, &mut k3);
        for i in 0..d {
            tmp[i] = h.mul_add(k3[i], x[i]);
        }
        problem.island(t + h, &tmp, &mut k4);
        for i in 0..d {
            x[i] += h / 6.0 * (k4[i] + k1[i] + 2.0 * (k2[i] + k3[i]));
        }
    }
    x
}

/// Shooting: Newton on the period map. Autonomous problems anchor
/// `x0[0]` at the seed value and treat the period as an unknown. On
/// convergence the monodromy matrix is differenced and its Floquet
/// multipliers extracted; a non-trivial unit-circle pair with nonzero
/// angle refuses as [`OrbitError::TorusSuspected`].
///
/// # Errors
/// Typed [`OrbitError`].
#[allow(clippy::too_many_lines)] // one coherent shooting stage
pub fn solve_shooting<P: OrbitProblem>(
    problem: &P,
    seed: &[f64],
    period_guess: f64,
    budget: &ShootBudget,
) -> Result<ShootOrbit, OrbitError> {
    let d = problem.dim();
    if seed.len() != d {
        return Err(OrbitError::BadParameter {
            what: "seed length must equal the state dimension",
        });
    }
    if !(period_guess > 0.0 && period_guess.is_finite()) {
        return Err(OrbitError::BadParameter {
            what: "period guess must be positive and finite",
        });
    }
    let auto = problem.autonomous();
    // Unknowns: forced -> x0 (period fixed); autonomous -> x0 with
    // x0[0] frozen at the seed, plus the period.
    let nn = d;
    let mut x0 = seed.to_vec();
    let mut period = period_guess;
    let mut trace = Vec::new();
    for _ in 0..budget.max_newton {
        let xt = flow(problem, &x0, period, budget.steps_per_period);
        let r: Vec<f64> = (0..d).map(|i| xt[i] - x0[i]).collect();
        let rn = norm(&r);
        trace.push(rn);
        let scale = norm(&x0).max(1.0);
        if rn < budget.tolerance * scale {
            // Equilibrium guard: a fixed point satisfies the period
            // map for EVERY period.
            let mut fx = vec![0.0f64; d];
            problem.island(0.0, &x0, &mut fx);
            if norm(&fx) < 1.0e-9 * scale {
                return Err(OrbitError::TrivialCollapse);
            }
            // Monodromy by forward differences.
            let mut mono = vec![0.0f64; d * d];
            for col in 0..d {
                let h = 1.0e-7 * x0[col].abs().max(1.0e-3);
                let mut xp = x0.clone();
                xp[col] += h;
                let xph = flow(problem, &xp, period, budget.steps_per_period);
                for row in 0..d {
                    mono[row * d + col] = (xph[row] - xt[row]) / h;
                }
            }
            let mc: Vec<C64> = mono.iter().map(|&v| C64::new(v, 0.0)).collect();
            let eigs = eig(&mc, d).map_err(OrbitError::Eigen)?;
            let multipliers: Vec<(f64, f64)> = eigs.iter().map(|e| (e.re, e.im)).collect();
            // Torus detection: a NON-trivial multiplier on the unit
            // circle away from the real axis.
            for &(re, im) in &multipliers {
                let mag = det::sqrt(re * re + im * im);
                let trivial = auto && (re - 1.0).abs() < 1.0e-3 && im.abs() < 1.0e-3;
                if !trivial && (mag - 1.0).abs() < 1.0e-4 && im.abs() > 1.0e-3 {
                    return Err(OrbitError::TorusSuspected {
                        multiplier: (re, im),
                    });
                }
            }
            return Ok(ShootOrbit {
                period,
                x0,
                multipliers,
                residual_trace: trace,
            });
        }
        // Jacobian of the period-map residual.
        let mut jac = vec![0.0f64; nn * nn];
        for col in 0..nn {
            let (mut xp, mut pp) = (x0.clone(), period);
            let h;
            if auto && col == 0 {
                // Column 0 differentiates the PERIOD (x0[0] anchored).
                h = 1.0e-7 * period.abs().max(1.0e-3);
                pp += h;
            } else {
                h = 1.0e-7 * x0[col].abs().max(1.0e-3);
                xp[col] += h;
            }
            let xph = flow(problem, &xp, pp, budget.steps_per_period);
            let xt0 = flow(problem, &x0, period, budget.steps_per_period);
            for row in 0..d {
                let base = xt0[row] - x0[row];
                let pert = xph[row] - xp[row];
                jac[row * nn + col] = (pert - base) / h;
            }
        }
        let Ok(factorized) = lu(&jac, nn) else {
            return Err(OrbitError::SingularJacobian);
        };
        let mut step: Vec<f64> = r.clone();
        factorized.solve(&mut step);
        let mut best = (1.0f64, f64::INFINITY);
        let mut alpha = 1.0f64;
        for _ in 0..6 {
            let (mut xp, mut pp) = (x0.clone(), period);
            for col in 0..nn {
                if auto && col == 0 {
                    pp -= alpha * step[col];
                } else {
                    xp[col] -= alpha * step[col];
                }
            }
            if pp > 0.0 && pp.is_finite() {
                let xt_try = flow(problem, &xp, pp, budget.steps_per_period);
                let tn = det::sqrt(
                    (0..d)
                        .map(|i| (xt_try[i] - xp[i]) * (xt_try[i] - xp[i]))
                        .sum::<f64>(),
                );
                if tn < best.1 {
                    best = (alpha, tn);
                }
                if tn < rn {
                    break;
                }
            }
            alpha *= 0.5;
        }
        for col in 0..nn {
            if auto && col == 0 {
                period -= best.0 * step[col];
            } else {
                x0[col] -= best.0 * step[col];
            }
        }
        if !(period > 0.0 && period.is_finite()) {
            return Err(OrbitError::BadParameter {
                what: "period left the positive reals during Newton",
            });
        }
    }
    Err(OrbitError::NewtonStalled {
        residual: *trace.last().unwrap_or(&f64::NAN),
        iterations: budget.max_newton,
        trace,
    })
}

// ---------------------------------------------------------------------
// Pseudo-arclength continuation.

/// A problem family over one physical parameter.
pub trait ContinuableProblem: OrbitProblem {
    /// Set the continuation parameter (blowing pressure, forcing
    /// frequency, lip tension, ...).
    fn set_parameter(&mut self, lambda: f64);
}

/// Continuation budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContinuationBudget {
    /// Maximum accepted steps.
    pub max_steps: usize,
    /// Initial parameter step.
    pub initial_step: f64,
    /// Minimum step before exhaustion refuses.
    pub min_step: f64,
    /// HB budget per corrector solve.
    pub hb: HbBudget,
}

/// One accepted continuation point.
#[derive(Debug, Clone)]
pub struct ContinuationPoint {
    /// Parameter value.
    pub lambda: f64,
    /// The orbit at this parameter.
    pub orbit: HbOrbit,
    /// Arclength step size that reached this point.
    pub step: f64,
}

/// TRUE pseudo-arclength continuation of a FORCED orbit branch in one
/// parameter: the corrector solves the bordered system (HB residuals
/// plus the arclength row `t . (z - z_pred) = 0`) for
/// `z = (coefficients, lambda)` JOINTLY, so folds in `lambda` are
/// TRAVERSED, not stalled at. Secant-tangent predictor, deterministic
/// step halving/growth, fixed budgets. v1 scope (disclosed): forced
/// anchors only — the parameter enters through
/// [`ContinuableProblem::set_parameter`] and doubles as the forcing
/// frequency when `omega_is_lambda` is set (the response-curve shape).
///
/// # Errors
/// Typed [`OrbitError`].
#[allow(clippy::too_many_lines)] // one coherent bordered-Newton stage
pub fn continue_branch<P: ContinuableProblem>(
    problem: &mut P,
    omega_of: impl Fn(f64) -> f64,
    lambda0: f64,
    lambda1: f64,
    guess_amplitude: f64,
    budget: &ContinuationBudget,
) -> Result<Vec<ContinuationPoint>, OrbitError> {
    if !(budget.initial_step > 0.0 && budget.initial_step.is_finite()) {
        return Err(OrbitError::BadParameter {
            what: "initial step must be positive and finite",
        });
    }
    let d = problem.dim();
    let pack = Packing {
        d,
        n_harm: budget.hb.harmonics,
        omega_slot: None,
    };
    let nc = pack.coeff_len();
    let nz = nc + 1; // coefficients + lambda
    // First point: plain HB at lambda0.
    problem.set_parameter(lambda0);
    let first = solve_hb(
        problem,
        HbAnchor::Forced {
            omega: omega_of(lambda0),
        },
        guess_amplitude,
        &budget.hb,
    )?;
    let mut z = pack_orbit(&first, &HbAnchor::Forced { omega: 1.0 });
    z.push(lambda0);
    let mut points = vec![ContinuationPoint {
        lambda: lambda0,
        orbit: first,
        step: 0.0,
    }];
    let dir = (lambda1 - lambda0).signum();
    // Initial tangent: pure-parameter march.
    let mut tangent = vec![0.0f64; nz];
    tangent[nc] = dir;
    let mut step = budget.initial_step;
    let mut residual = vec![0.0f64; nc];
    let eval = |problem: &mut P, z: &[f64], out: &mut [f64]| {
        problem.set_parameter(z[nc]);
        hb_residual(problem, &pack, &z[..nc], omega_of(z[nc]), out);
    };
    for _ in 0..budget.max_steps {
        // Reached the far end of the parameter window?
        let lam = z[nc];
        if (dir > 0.0 && lam >= lambda1) || (dir < 0.0 && lam <= lambda1) {
            return Ok(points);
        }
        // Predictor.
        let mut zp: Vec<f64> = z.iter().zip(&tangent).map(|(a, t)| a + step * t).collect();
        // Corrector: bordered Newton.
        let mut converged = false;
        for _ in 0..budget.hb.max_newton {
            eval(problem, &zp, &mut residual);
            let mut r = residual.clone();
            let arc: f64 = zp
                .iter()
                .zip(z.iter())
                .zip(&tangent)
                .map(|((a, b), t)| (a - b) * t)
                .sum::<f64>()
                - step;
            r.push(arc);
            let rn = norm(&r);
            let scale = norm(&zp[..nc]).max(1.0);
            if rn < budget.hb.tolerance * scale {
                converged = true;
                break;
            }
            let mut jac = vec![0.0f64; nz * nz];
            let mut r2 = vec![0.0f64; nc];
            for col in 0..nz {
                let h = 1.0e-7 * zp[col].abs().max(1.0);
                let saved = zp[col];
                zp[col] = saved + h;
                eval(problem, &zp, &mut r2);
                zp[col] = saved;
                for row in 0..nc {
                    jac[row * nz + col] = (r2[row] - residual[row]) / h;
                }
                jac[nc * nz + col] = tangent[col];
            }
            eval(problem, &zp, &mut residual);
            let Ok(factorized) = lu(&jac, nz) else {
                return Err(OrbitError::SingularJacobian);
            };
            let mut delta = r;
            factorized.solve(&mut delta);
            for col in 0..nz {
                zp[col] -= delta[col];
            }
        }
        if converged {
            let new_tangent: Vec<f64> = {
                let diff: Vec<f64> = zp.iter().zip(z.iter()).map(|(a, b)| a - b).collect();
                let nrm = norm(&diff).max(1.0e-300);
                diff.iter().map(|v| v / nrm).collect()
            };
            tangent = new_tangent;
            z = zp;
            let (coeffs, _) = pack.unpack(&z[..nc]);
            points.push(ContinuationPoint {
                lambda: z[nc],
                orbit: HbOrbit {
                    omega: omega_of(z[nc]),
                    coeffs,
                    residual_trace: vec![],
                    residual: 0.0,
                },
                step,
            });
            step = (step * 1.3).min(budget.initial_step * 4.0);
        } else if step > budget.min_step {
            step *= 0.5;
        } else {
            return Err(OrbitError::ContinuationExhausted {
                steps: points.len(),
            });
        }
    }
    Err(OrbitError::ContinuationExhausted {
        steps: points.len(),
    })
}

/// Repack an orbit as a seed vector for the given anchor.
fn pack_orbit(orbit: &HbOrbit, anchor: &HbAnchor) -> Vec<f64> {
    let d = orbit.coeffs.len();
    let n_harm = orbit.coeffs[0].len() - 1;
    let pack = Packing {
        d,
        n_harm,
        omega_slot: match anchor {
            HbAnchor::Forced { .. } => None,
            _ => Some(d * (2 * n_harm + 1)),
        },
    };
    let mut u = vec![0.0f64; pack.len()];
    for k in 0..d {
        u[pack.re_index(0, k)] = orbit.coeffs[k][0].re;
        for n in 1..=n_harm {
            u[pack.re_index(n, k)] = orbit.coeffs[k][n].re;
            u[pack.im_index(n, k)] = orbit.coeffs[k][n].im;
        }
    }
    if let Some(slot) = pack.omega_slot {
        u[slot] = orbit.omega;
    }
    u
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

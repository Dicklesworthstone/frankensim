//! Prewarped bilinear discretization of a fitted model to runtime
//! filter forms: factored section cascade (first-order + biquad
//! sections summed in PARALLEL, mirroring the pole-residue structure)
//! and discrete state-space.
//!
//! Bilinear map `s = K (z - 1)/(z + 1)` with `K = omega_pw /
//! tan(omega_pw * T / 2)` (prewarp: the continuous and discrete
//! responses agree EXACTLY at `omega_pw`; `K -> 2/T` as `omega_pw ->
//! 0`). The improper `s*e` term maps to `e*K (z-1)/(z+1)` — a proper
//! first-order z-section, so improper continuous models are exactly
//! representable after discretization.
//!
//! Sections are PARALLEL (summed), not cascaded products: parallel
//! form preserves the per-resonance structure the pole-residue model
//! carries, keeps coefficient sensitivity local to each resonance, and
//! makes the quantization report per-section attributable.

use crate::model::{PoleTerm, RationalModel};
use fs_math::c64::C64;
use fs_math::det;

/// One second-order (or degenerate first-order) parallel section:
/// `(b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Biquad {
    /// Numerator coefficients.
    pub b: [f64; 3],
    /// Denominator coefficients (`a0 == 1` implied).
    pub a: [f64; 2],
}

impl Biquad {
    /// Frequency response at `z = e^{i*omega*t_s}`.
    #[must_use]
    pub fn eval(&self, omega: f64, t_s: f64) -> C64 {
        let zi = C64::new(det::cos(omega * t_s), -det::sin(omega * t_s));
        let zi2 = zi * zi;
        let num = C64::from_re(self.b[0]) + zi.scale(self.b[1]) + zi2.scale(self.b[2]);
        let den = C64::ONE + zi.scale(self.a[0]) + zi2.scale(self.a[1]);
        num * den.recip()
    }

    /// The same response with coefficients rounded through f32 — the
    /// quantization-sensitivity probe.
    #[must_use]
    pub fn eval_f32_quantized(&self, omega: f64, t_s: f64) -> C64 {
        let q = Biquad {
            b: [
                f64::from(self.b[0] as f32),
                f64::from(self.b[1] as f32),
                f64::from(self.b[2] as f32),
            ],
            a: [f64::from(self.a[0] as f32), f64::from(self.a[1] as f32)],
        };
        q.eval(omega, t_s)
    }

    /// Section poles inside the unit circle?
    #[must_use]
    pub fn is_stable(&self) -> bool {
        // Jury criterion for 1 + a1 z^-1 + a2 z^-2.
        let (a1, a2) = (self.a[0], self.a[1]);
        a2.abs() < 1.0 && (a1.abs() < 1.0 + a2)
    }
}

/// A parallel bank of sections plus the sample interval.
#[derive(Debug, Clone, PartialEq)]
pub struct DigitalFilter {
    /// Parallel sections (their responses SUM).
    pub sections: Vec<Biquad>,
    /// Constant (direct) term.
    pub direct: f64,
    /// Sample interval [s].
    pub t_s: f64,
    /// Prewarp frequency actually used [rad/s] (0 = none).
    pub prewarp: f64,
}

impl DigitalFilter {
    /// Total response at angular frequency `omega`.
    #[must_use]
    pub fn eval(&self, omega: f64) -> C64 {
        let mut acc = C64::from_re(self.direct);
        for s in &self.sections {
            acc = acc + s.eval(omega, self.t_s);
        }
        acc
    }

    /// Total response with every section's coefficients f32-quantized.
    #[must_use]
    pub fn eval_f32_quantized(&self, omega: f64) -> C64 {
        let mut acc = C64::from_re(f64::from(self.direct as f32));
        for s in &self.sections {
            acc = acc + s.eval_f32_quantized(omega, self.t_s);
        }
        acc
    }

    /// All sections stable? The exact lossless differentiator the
    /// improper `e` term maps to (`b1 = -b0`, pole exactly at
    /// `z = -1`) is marginally stable BY DESIGN and is admitted.
    #[must_use]
    pub fn is_stable(&self) -> bool {
        self.sections
            .iter()
            .all(|s| s.is_stable() || is_lossless_differentiator(s))
    }
}

fn is_lossless_differentiator(s: &Biquad) -> bool {
    // Exact bitwise structural match against the section `bilinear`
    // emits for the e-term — an identity check, not a numeric
    // comparison, so strict equality is the correct operator.
    #[allow(clippy::float_cmp)]
    fn exact(s: &Biquad) -> bool {
        s.b[1] == -s.b[0] && s.b[2] == 0.0 && s.a == [1.0, 0.0]
    }
    exact(s)
}

/// Typed discretization failure.
#[derive(Debug, Clone, PartialEq)]
pub enum DiscretizeError {
    /// Sample rate too low: a pole (or the prewarp point) sits at or
    /// beyond Nyquist.
    BeyondNyquist {
        /// Offending frequency [rad/s].
        omega: f64,
        /// Nyquist [rad/s].
        nyquist: f64,
    },
    /// Non-positive sample interval.
    BadSampleInterval,
}

impl core::fmt::Display for DiscretizeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DiscretizeError::BeyondNyquist { omega, nyquist } => {
                write!(
                    f,
                    "frequency {omega} rad/s at/beyond nyquist {nyquist} rad/s"
                )
            }
            DiscretizeError::BadSampleInterval => write!(f, "sample interval must be positive"),
        }
    }
}

impl std::error::Error for DiscretizeError {}

/// Bilinear-transform the model at sample interval `t_s` with prewarp
/// at `omega_pw` (pass 0 for the unwarped `K = 2/T` map).
///
/// # Errors
/// [`DiscretizeError`] when the sample rate cannot represent the model
/// (pole resonance or prewarp at/beyond Nyquist) or `t_s <= 0`.
pub fn bilinear(
    model: &RationalModel,
    t_s: f64,
    omega_pw: f64,
) -> Result<DigitalFilter, DiscretizeError> {
    if t_s <= 0.0 || t_s.is_nan() {
        return Err(DiscretizeError::BadSampleInterval);
    }
    let nyquist = core::f64::consts::PI / t_s;
    if omega_pw >= nyquist {
        return Err(DiscretizeError::BeyondNyquist {
            omega: omega_pw,
            nyquist,
        });
    }
    for t in &model.terms {
        let w = match t {
            PoleTerm::Real { pole, .. } => pole.abs(),
            PoleTerm::Pair { pole, .. } => pole.abs(),
        };
        if w >= nyquist {
            return Err(DiscretizeError::BeyondNyquist { omega: w, nyquist });
        }
    }
    let k = if omega_pw > 0.0 {
        omega_pw / det::tan(omega_pw * t_s / 2.0)
    } else {
        2.0 / t_s
    };
    let mut sections = Vec::new();
    for t in &model.terms {
        match *t {
            PoleTerm::Real { pole, residue } => {
                // r/(s-p), s = K(z-1)/(z+1):
                //   r (1 + z^-1) / ((K - p) + (-K - p) z^-1)
                let a0 = k - pole;
                sections.push(Biquad {
                    b: [residue / a0, residue / a0, 0.0],
                    a: [(-k - pole) / a0, 0.0],
                });
            }
            PoleTerm::Pair { pole, residue } => {
                // Pair as a real second-order s-section:
                //   (beta1 s + beta0) / (s^2 + alpha1 s + alpha0)
                let (a_re, b_im) = (pole.re, pole.im);
                let (rho, sigma) = (residue.re, residue.im);
                let beta1 = 2.0 * rho;
                let beta0 = -2.0 * (rho * a_re + sigma * b_im);
                let alpha1 = -2.0 * a_re;
                let alpha0 = a_re * a_re + b_im * b_im;
                // Substitute and clear (z+1)^2:
                //   num = beta1 K (z^2 - 1) + beta0 (z+1)^2
                //   den = K^2 (z-1)^2 + alpha1 K (z^2 - 1) + alpha0 (z+1)^2
                let k2 = k * k;
                let n0 = beta1 * k + beta0; // z^2
                let n1 = 2.0 * beta0; // z^1
                let n2 = -beta1 * k + beta0; // z^0
                let d0 = k2 + alpha1 * k + alpha0;
                let d1 = -2.0 * k2 + 2.0 * alpha0;
                let d2 = k2 - alpha1 * k + alpha0;
                sections.push(Biquad {
                    b: [n0 / d0, n1 / d0, n2 / d0],
                    a: [d1 / d0, d2 / d0],
                });
            }
        }
    }
    if model.e != 0.0 {
        // e*s -> e K (z - 1)/(z + 1): first-order, pole at z = -1
        // nudged is NOT needed — bilinear puts it exactly at z=-1 which
        // rings at Nyquist; standard practice keeps it exact (the
        // section is lossless) and the band tolerance test covers it.
        sections.push(Biquad {
            b: [model.e * k, -model.e * k, 0.0],
            a: [1.0, 0.0],
        });
    }
    Ok(DigitalFilter {
        sections,
        direct: model.d,
        t_s,
        prewarp: omega_pw,
    })
}

/// Discrete state-space by the same bilinear map (Tustin):
/// `Ad = (K I - A)^{-1} (K I + A)`, `Bd = sqrt(2 K) (K I - A)^{-1} B`,
/// `Cd = sqrt(2 K) C (K I - A)^{-1}`, `Dd = D + C (K I - A)^{-1} B`.
/// The improper `e` term is NOT representable in proper discrete
/// state-space and is returned as the leftover coefficient for the
/// caller to realize as the extra first-order section (as
/// [`bilinear`] does) — a named seam, not a silent drop.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscreteStateSpace {
    /// State dimension.
    pub n: usize,
    /// Row-major discrete state matrix.
    pub a: Vec<f64>,
    /// Input map.
    pub b: Vec<f64>,
    /// Output map.
    pub c: Vec<f64>,
    /// Direct term.
    pub d: f64,
    /// UNREALIZED improper coefficient (see type docs).
    pub e_leftover: f64,
    /// Sample interval.
    pub t_s: f64,
}

/// Tustin state-space discretization (see [`DiscreteStateSpace`]).
///
/// # Errors
/// As [`bilinear`], plus singular `(K I - A)` (a pole exactly at `K`,
/// impossible for stable models since `K > 0`).
pub fn bilinear_state_space(
    model: &RationalModel,
    t_s: f64,
    omega_pw: f64,
) -> Result<DiscreteStateSpace, DiscretizeError> {
    if t_s <= 0.0 || t_s.is_nan() {
        return Err(DiscretizeError::BadSampleInterval);
    }
    let nyquist = core::f64::consts::PI / t_s;
    if omega_pw >= nyquist {
        return Err(DiscretizeError::BeyondNyquist {
            omega: omega_pw,
            nyquist,
        });
    }
    let k = if omega_pw > 0.0 {
        omega_pw / det::tan(omega_pw * t_s / 2.0)
    } else {
        2.0 / t_s
    };
    let ss = model.state_space();
    let n = ss.n;
    // M = K I - A, factored once.
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = -ss.a[i * n + j];
        }
        m[i * n + i] += k;
    }
    let fact = fs_la::factor::lu(&m, n)
        .map_err(|_| DiscretizeError::BeyondNyquist { omega: k, nyquist })?;
    // X = M^{-1} (K I + A) column by column; Y = M^{-1} B.
    let mut ad = vec![0.0f64; n * n];
    for col in 0..n {
        let mut rhs: Vec<f64> = (0..n).map(|row| ss.a[row * n + col]).collect();
        rhs[col] += k;
        fact.solve(&mut rhs);
        for row in 0..n {
            ad[row * n + col] = rhs[row];
        }
    }
    let mut y = ss.b.clone();
    fact.solve(&mut y);
    let g = det::sqrt(2.0 * k);
    let bd: Vec<f64> = y.iter().map(|&v| g * v).collect();
    // Cd = g * C M^{-1}  (solve M^T z = C^T).
    let mut z = ss.c.clone();
    fact.solve_transpose(&mut z);
    let cd: Vec<f64> = z.iter().map(|&v| g * v).collect();
    let mut dd = ss.d;
    for (ci, yi) in ss.c.iter().zip(&y) {
        dd += ci * yi;
    }
    Ok(DiscreteStateSpace {
        n,
        a: ad,
        b: bd,
        c: cd,
        d: dd,
        e_leftover: ss.e,
        t_s,
    })
}

impl DiscreteStateSpace {
    /// Frequency response `Cd (z I - Ad)^{-1} Bd + Dd` at `z =
    /// e^{i*omega*t_s}` (the `e_leftover` term is the caller's
    /// section).
    ///
    /// # Errors
    /// Singular `(z I - Ad)` (resonance exactly on the unit circle).
    pub fn eval(&self, omega: f64) -> Result<C64, fs_la::eigen_complex::EigFailure> {
        let n = self.n;
        let z = C64::new(det::cos(omega * self.t_s), det::sin(omega * self.t_s));
        let mut m = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                m[i * n + j] = C64::from_re(-self.a[i * n + j]);
            }
            m[i * n + i] = m[i * n + i] + z;
        }
        let lu = fs_la::eigen_complex::lu_complex(&m, n)?;
        let mut x: Vec<C64> = self.b.iter().map(|&v| C64::from_re(v)).collect();
        lu.solve(&mut x);
        let mut acc = C64::from_re(self.d);
        for (ci, xi) in self.c.iter().zip(&x) {
            acc = acc + xi.scale(*ci);
        }
        Ok(acc)
    }
}

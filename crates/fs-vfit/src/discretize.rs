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
    /// A public state-space field does not match the declared state
    /// dimension.
    InvalidRuntimeDimensions {
        /// Field whose length is inconsistent.
        field: &'static str,
        /// Required number of scalar entries.
        expected: usize,
        /// Supplied number of scalar entries.
        actual: usize,
    },
    /// Squaring the declared state dimension overflowed `usize`.
    RuntimeDimensionOverflow {
        /// Declared state dimension.
        n: usize,
    },
    /// A runtime input, coefficient, state coordinate, or computed result is
    /// not finite.
    NonFiniteRuntimeValue {
        /// Value source.
        field: &'static str,
        /// Coordinate for vector/matrix fields; absent for scalars.
        index: Option<usize>,
    },
    /// The proper state-space step cannot silently omit the Tustin
    /// differentiator represented by `e_leftover`.
    UnrealizedImproperTerm {
        /// Coefficient requiring a separate differentiator section.
        e_leftover: f64,
    },
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
            DiscretizeError::InvalidRuntimeDimensions {
                field,
                expected,
                actual,
            } => write!(
                f,
                "runtime state-space field {field} has length {actual}, expected {expected}"
            ),
            DiscretizeError::RuntimeDimensionOverflow { n } => {
                write!(f, "runtime state dimension {n} cannot be squared")
            }
            DiscretizeError::NonFiniteRuntimeValue { field, index } => match index {
                Some(index) => write!(f, "runtime value {field}[{index}] must be finite"),
                None => write!(f, "runtime value {field} must be finite"),
            },
            DiscretizeError::UnrealizedImproperTerm { e_leftover } => write!(
                f,
                "runtime state-space step cannot realize e_leftover={e_leftover}; use the separate differentiator section"
            ),
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

/// Mutable, allocation-free-after-construction runtime for one proper
/// [`DiscreteStateSpace`] realization.
///
/// The lifetime binds the state to the exact realization used to construct it;
/// callers cannot accidentally step that state with a different realization.
/// A nonzero [`DiscreteStateSpace::e_leftover`] is refused at construction
/// because silently dropping its separate differentiator section would change
/// the transfer function.
#[derive(Debug)]
pub struct DiscreteStateSpaceRuntime<'realization> {
    realization: &'realization DiscreteStateSpace,
    state: Vec<f64>,
    next_state: Vec<f64>,
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
    // Same per-pole refusal as the section route (review finding: a
    // resonance at/beyond Nyquist must not alias silently through the
    // Tustin route either).
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
    /// Construct a zero-state runtime bound to this proper realization.
    ///
    /// # Errors
    /// Returns [`DiscretizeError`] if public realization fields are malformed
    /// or nonfinite, or if `e_leftover` requires the separate differentiator
    /// section.
    pub fn try_runtime(&self) -> Result<DiscreteStateSpaceRuntime<'_>, DiscretizeError> {
        self.validate_runtime_realization()?;
        Ok(DiscreteStateSpaceRuntime {
            realization: self,
            state: vec![0.0; self.n],
            next_state: vec![0.0; self.n],
        })
    }

    fn validate_runtime_realization(&self) -> Result<(), DiscretizeError> {
        let matrix_len = self
            .n
            .checked_mul(self.n)
            .ok_or(DiscretizeError::RuntimeDimensionOverflow { n: self.n })?;
        require_len("a", self.a.len(), matrix_len)?;
        require_len("b", self.b.len(), self.n)?;
        require_len("c", self.c.len(), self.n)?;
        require_finite("a", &self.a)?;
        require_finite("b", &self.b)?;
        require_finite("c", &self.c)?;
        require_finite_scalar("d", self.d)?;
        require_finite_scalar("e_leftover", self.e_leftover)?;
        require_finite_scalar("t_s", self.t_s)?;
        if self.t_s <= 0.0 {
            return Err(DiscretizeError::BadSampleInterval);
        }
        if self.e_leftover != 0.0 {
            return Err(DiscretizeError::UnrealizedImproperTerm {
                e_leftover: self.e_leftover,
            });
        }
        Ok(())
    }

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

impl DiscreteStateSpaceRuntime<'_> {
    /// Current realization state in canonical coordinate order.
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.state
    }

    /// Advance one scalar sample using `y = Cx + D u`, followed by
    /// `x_next = A x + B u`.
    ///
    /// Arithmetic order is fixed by row and coordinate index. If finite inputs
    /// overflow, the call refuses without committing the partially computed
    /// next state.
    ///
    /// # Errors
    /// [`DiscretizeError::NonFiniteRuntimeValue`] for a nonfinite input or
    /// computed output/state coordinate.
    pub fn step(&mut self, input: f64) -> Result<f64, DiscretizeError> {
        require_finite_scalar("input", input)?;

        let mut output = 0.0;
        for index in 0..self.realization.n {
            output += self.realization.c[index] * self.state[index];
        }
        output += self.realization.d * input;
        require_finite_scalar("output", output)?;

        for row in 0..self.realization.n {
            let mut value = 0.0;
            for column in 0..self.realization.n {
                value += self.realization.a[row * self.realization.n + column] * self.state[column];
            }
            value += self.realization.b[row] * input;
            if !value.is_finite() {
                return Err(DiscretizeError::NonFiniteRuntimeValue {
                    field: "next_state",
                    index: Some(row),
                });
            }
            self.next_state[row] = value;
        }
        core::mem::swap(&mut self.state, &mut self.next_state);
        Ok(output)
    }
}

fn require_len(field: &'static str, actual: usize, expected: usize) -> Result<(), DiscretizeError> {
    if actual == expected {
        Ok(())
    } else {
        Err(DiscretizeError::InvalidRuntimeDimensions {
            field,
            expected,
            actual,
        })
    }
}

fn require_finite(field: &'static str, values: &[f64]) -> Result<(), DiscretizeError> {
    if let Some(index) = values.iter().position(|value| !value.is_finite()) {
        Err(DiscretizeError::NonFiniteRuntimeValue {
            field,
            index: Some(index),
        })
    } else {
        Ok(())
    }
}

fn require_finite_scalar(field: &'static str, value: f64) -> Result<(), DiscretizeError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(DiscretizeError::NonFiniteRuntimeValue { field, index: None })
    }
}

#[cfg(test)]
mod runtime_tests {
    use super::*;

    fn first_order() -> DiscreteStateSpace {
        DiscreteStateSpace {
            n: 1,
            a: vec![0.5],
            b: vec![1.0],
            c: vec![2.0],
            d: 3.0,
            e_leftover: 0.0,
            t_s: 0.01,
        }
    }

    #[test]
    fn g0_runtime_step_matches_analytic_first_order_recurrence() {
        let system = first_order();
        let mut runtime = system.try_runtime().expect("proper realization");
        assert_eq!(runtime.step(1.0).unwrap(), 3.0);
        assert_eq!(runtime.state(), [1.0]);
        assert_eq!(runtime.step(2.0).unwrap(), 8.0);
        assert_eq!(runtime.state(), [2.5]);
    }

    #[test]
    fn g0_runtime_replay_is_bitwise() {
        let system = first_order();
        let mut left = system.try_runtime().unwrap();
        let mut right = system.try_runtime().unwrap();
        for input in [1.0, -0.25, 0.0, 2.5, -4.0] {
            assert_eq!(
                left.step(input).unwrap().to_bits(),
                right.step(input).unwrap().to_bits()
            );
            assert_eq!(left.state(), right.state());
        }
    }

    #[test]
    fn g0_zero_state_and_input_remain_zero() {
        let mut system = first_order();
        system.d = 0.0;
        let mut runtime = system.try_runtime().unwrap();
        assert_eq!(runtime.step(0.0).unwrap().to_bits(), 0.0f64.to_bits());
        assert_eq!(runtime.state(), [0.0]);
    }

    #[test]
    fn g0_runtime_refuses_malformed_nonfinite_and_improper_inputs() {
        let mut malformed = first_order();
        malformed.a.clear();
        assert!(matches!(
            malformed.try_runtime(),
            Err(DiscretizeError::InvalidRuntimeDimensions { field: "a", .. })
        ));
        let mut improper = first_order();
        improper.e_leftover = 0.25;
        assert!(matches!(
            improper.try_runtime(),
            Err(DiscretizeError::UnrealizedImproperTerm { .. })
        ));
        let system = first_order();
        let mut runtime = system.try_runtime().unwrap();
        assert!(matches!(
            runtime.step(f64::NAN),
            Err(DiscretizeError::NonFiniteRuntimeValue { field: "input", .. })
        ));
        assert_eq!(runtime.state(), [0.0]);
    }

    #[test]
    fn g0_impulse_response_agrees_with_frequency_evaluation() {
        let system = first_order();
        let omega = 37.0;
        let expected = system.eval(omega).unwrap();
        let mut runtime = system.try_runtime().unwrap();
        let mut measured = C64::ZERO;
        for sample in 0..96 {
            let output = runtime.step(if sample == 0 { 1.0 } else { 0.0 }).unwrap();
            let phase = -(sample as f64) * omega * system.t_s;
            measured = measured + C64::new(det::cos(phase), det::sin(phase)).scale(output);
        }
        assert!((measured - expected).abs() < 1.0e-12);
    }
}

//! # fs-nlmodal — geometric nonlinearity in modal coordinates
//!
//! Thin-structure geometric nonlinearity (von Karman plates,
//! Kirchhoff-Carrier strings) reduced to N coupled modal oscillators
//! with a QUARTIC SUM-OF-SQUARES potential:
//!
//! `H = 1/2 sum_k (p_k^2 + w_k^2 q_k^2) + 1/4 sum_j c_j (q^T E_j q)^2`
//!
//! with `c_j >= 0` and `E_j` symmetric — the Airy-stress-eliminated von
//! Karman form (Ducceschi-Touze / Bilbao modal formulation). Because
//! the nonlinear potential is a positive sum of squares, `H` is
//! bounded below and the whole system is a PORT-HAMILTONIAN system
//! with non-quadratic storage: time stepping, damping, striking, and
//! the original-energy ledger come from fs-phs discrete gradients.
//! The opt-in [`sav`] path instead uses exact linear propagation and one
//! non-iterative auxiliary-variable solve. Its conserved numerical energy
//! and its discrepancy from the physical quartic energy are reported separately.
//!
//! Constructors:
//! - [`von_karman_ss_plate`] — simply-supported rectangular plate with
//!   ANALYTIC sine modes for displacement AND Airy stress (both are
//!   biharmonic eigenfunctions on the SS rectangle); coupling
//!   integrals by exact trigonometric products, with integer-mode
//!   selection rules and no spatial quadrature aliasing.
//! - [`von_karman_sampled_plate`] — isotropic membrane elasticity on a
//!   *sampled* displacement grid (FE, clamped). The tensor counterpart
//!   [`von_karman_sampled_plate_with_membrane`] retains directional normal,
//!   shear and rotated-axis elasticity in the same Airy basis.
//! - [`kirchhoff_carrier_string`] — the 1D sibling: one stress "mode"
//!   whose coupling matrix is exactly `diag(k^2 pi^2/L^2)` — tension
//!   modulation and pitch glide nearly free. Fixed-fixed sines.
//! - [`kirchhoff_carrier_moving_end`] — the same channel on
//!   free-fixed cosines `κ = (k−½)π/L`, `φ(0) ≠ 0`. A cable on a
//!   moving support, not a second instrument law.
//! - [`longitudinal::LongitudinalString`] — finite longitudinal inertia
//!   coupled to one or two transverse slope fields. Unretained axial
//!   strain channels are statically equilibrated; the zero-mode limit
//!   recovers the same Kirchhoff--Carrier mean strain.
//!
//! Honest scope: moderate-rotation von Karman regime; no damage,
//! plasticity, or wrinkling. In-plane edges remain movable (Airy
//! sine); out-of-plane shape is whatever grid the caller samples.

use fs_math::det;
use fs_phs::{
    ConstitutiveDiscreteGradient, ContinuousStorageTangent, PhsError, PortHamiltonian, Storage,
};

mod membrane;
pub use membrane::{VkSampledPlateParams, von_karman_sampled_plate_with_membrane};

pub mod backbone;
pub mod component;
pub mod loaded;
pub mod longitudinal;
pub mod parameter;
pub mod preload;
pub mod sav;
pub use backbone::{
    DuffingBackboneDerivative, DuffingBackboneDirection, DuffingBackboneError,
    DuffingKinematicSpectrum, DuffingSpectrum, MAX_DUFFING_HARMONICS,
    duffing_backbone_directional_derivative, duffing_backbone_exact,
    duffing_kinematic_spectrum_exact, duffing_spectrum_exact,
};
pub use parameter::{SosModalParameterDirection, StressChannelParameterDirection};
pub use preload::{
    PreloadedModalStorage, PreloadedSpectrum, PreloadedSpectrumBudget, PreloadedSpectrumError,
    PreloadedSpectrumGuard,
};

/// One Airy-stress coupling channel: coefficient `c >= 0` and the
/// symmetric matrix `E` (row-major `n x n`) of `q^T E q`.
#[derive(Debug, Clone)]
pub struct StressChannel {
    /// Positive channel coefficient `c_j` (absorbs `E h / (2 xi^4)`).
    pub coefficient: f64,
    /// Symmetric coupling matrix, row-major `n x n`.
    pub coupling: Vec<f64>,
}

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum NlModalError {
    /// Bad geometry/material parameter.
    Parameter {
        /// Which one.
        what: &'static str,
    },
    /// A coupling matrix failed its symmetry certificate.
    AsymmetricCoupling {
        /// Worst residual relative to the matrix scale.
        residual: f64,
    },
    /// The two independent quadrature orders disagreed beyond
    /// tolerance — the tensor is not converged.
    QuadratureMismatch {
        /// Worst relative disagreement.
        residual: f64,
    },
    /// Underlying pHS admission failure.
    Phs(PhsError),
}

impl core::fmt::Display for NlModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            NlModalError::Parameter { what } => write!(f, "bad parameter: {what}"),
            NlModalError::AsymmetricCoupling { residual } => {
                write!(f, "coupling matrix asymmetry {residual:.3e}")
            }
            NlModalError::QuadratureMismatch { residual } => {
                write!(f, "quadrature orders disagree by {residual:.3e}")
            }
            NlModalError::Phs(e) => write!(f, "pHS admission failed: {e}"),
        }
    }
}

impl std::error::Error for NlModalError {}

/// The quartic sum-of-squares storage
/// `H = 1/2 sum (p^2 + w^2 q^2) + 1/4 sum_j c_j (q^T E_j q)^2`.
///
/// State layout matches fs-phs `modal_bank`: `[q_0, p_0, q_1, p_1,
/// ...]`. The gradient is the EXACT analytic gradient of the coded
/// `H` — force-vs-energy divergence (the classic hand-derived-tensor
/// bug class) is therefore impossible through this type; the battery
/// demonstrates that under Gonzalez stepping such divergence is a
/// TRAJECTORY error, never an energy error (the architectural
/// finding pinned in the tests).
pub struct SosModalStorage {
    /// Modal angular frequencies [rad/s].
    pub omegas: Vec<f64>,
    /// Stress channels.
    pub channels: Vec<StressChannel>,
}

impl SosModalStorage {
    fn n_modes(&self) -> usize {
        self.omegas.len()
    }

    /// Admit the public frequency and quartic-channel data before wrapping
    /// this storage in a multiport or full-viscosity system. This is the same
    /// storage admission used by [`assemble`]; it changes no coefficients.
    ///
    /// # Errors
    /// Empty/nonfinite/nonpositive frequencies, malformed or nonfinite
    /// channels, negative quartic weights, or an asymmetric coupling matrix.
    pub fn validate(&self) -> Result<(), NlModalError> {
        let n = self.n_modes();
        let Some(square) = n.checked_mul(n) else {
            return Err(NlModalError::Parameter {
                what: "coupling matrix size overflow",
            });
        };
        if n == 0 || self.omegas.iter().any(|w| !w.is_finite() || *w <= 0.0) {
            return Err(NlModalError::Parameter {
                what: "modal frequency must be positive and finite",
            });
        }
        for ch in &self.channels {
            if ch.coefficient < 0.0 || !ch.coefficient.is_finite() {
                return Err(NlModalError::Parameter {
                    what: "channel coefficient must be non-negative and finite",
                });
            }
            if ch.coupling.iter().any(|v| !v.is_finite()) {
                return Err(NlModalError::Parameter {
                    what: "coupling entries must be finite",
                });
            }
            if ch.coupling.len() != square {
                return Err(NlModalError::Parameter {
                    what: "coupling matrix size",
                });
            }
            let scale = ch
                .coupling
                .iter()
                .fold(f64::MIN_POSITIVE, |a, &v| a.max(v.abs()));
            let mut worst = 0.0f64;
            for p in 0..n {
                for q in 0..p {
                    worst = worst.max((ch.coupling[p * n + q] - ch.coupling[q * n + p]).abs());
                }
            }
            if worst > SYM_TOL * scale {
                return Err(NlModalError::AsymmetricCoupling {
                    residual: worst / scale,
                });
            }
        }
        Ok(())
    }

    /// Exact conservative frequency of one invariant displacement axis.
    ///
    /// The physical peak displacement is in the original mass-normalized
    /// coordinates. This uses the SAME source stiffness and stress channels;
    /// every active channel must leave the selected axis invariant. In
    /// particular, every Kirchhoff-Carrier axis is admitted, while a plate mode
    /// that drives another mode refuses instead of silently dropping coupling.
    /// Channel-by-channel invariance is sufficient, not necessary: cancellations
    /// between distinct noninvariant channels are not inferred here.
    ///
    /// The result can seed the frequency of a nonlinear periodic-orbit solve.
    /// Its turning amplitude is not the solver's first-harmonic amplitude, and
    /// its existence does not imply transverse or forced-response stability.
    /// No state, model coefficient, damping or assembled dynamics is changed.
    pub fn invariant_axis_backbone(
        &self,
        mode: usize,
        amplitude: f64,
    ) -> Result<f64, NlModalError> {
        self.validate()?;
        backbone::invariant_axis_backbone(
            &self.omegas,
            self.channels
                .iter()
                .map(|channel| (channel.coefficient, channel.coupling.as_slice())),
            mode,
            amplitude,
        )
        .map_err(|error| NlModalError::Parameter {
            what: if error == DuffingBackboneError::NonInvariantAxis {
                "backbone axis couples to other displacement coordinates"
            } else {
                "backbone mode, amplitude or scalar frequency is not representable"
            },
        })
    }

    /// Exact local frequency change per unit of physical turning amplitude.
    ///
    /// This differentiates the same source-bound conservative orbit as
    /// [`Self::invariant_axis_backbone`], including finite-amplitude hardening.
    /// The derivative uses the original mass-normalized displacement units;
    /// it is not a weak-amplitude approximation or a forced-peak derivative.
    /// Every active stress channel must preserve the selected axis. No source
    /// coefficients, state, damping, or coupled dynamics are changed.
    pub fn invariant_axis_backbone_amplitude_derivative(
        &self,
        mode: usize,
        amplitude: f64,
    ) -> Result<DuffingBackboneDerivative, NlModalError> {
        self.validate()?;
        backbone::invariant_axis_backbone_amplitude_derivative(
            &self.omegas,
            self.channels
                .iter()
                .map(|channel| (channel.coefficient, channel.coupling.as_slice())),
            mode,
            amplitude,
        )
        .map_err(|error| NlModalError::Parameter {
            what: if error == DuffingBackboneError::NonInvariantAxis {
                "backbone axis couples to other displacement coordinates"
            } else {
                "backbone mode, amplitude, frequency or derivative is not representable"
            },
        })
    }

    /// Nonlinear harmonic amplitudes of an original invariant displacement axis.
    ///
    /// This shares [`Self::invariant_axis_backbone`]'s source coefficients and
    /// channel-by-channel invariance admission. `amplitude` is a physical
    /// turning displacement; the returned fundamental need not equal it.
    /// `odd_harmonics` retains the explicitly requested odd orders, with a
    /// separate omitted-displacement bound. Increasing it never renormalizes
    /// existing coefficients. See [`duffing_spectrum_exact`] for the Fourier
    /// convention and the conservative reference/solver-seed scope.
    pub fn invariant_axis_spectrum(
        &self,
        mode: usize,
        amplitude: f64,
        odd_harmonics: usize,
    ) -> Result<DuffingSpectrum, NlModalError> {
        self.validate()?;
        let result = backbone::invariant_axis_parameters(
            &self.omegas,
            self.channels
                .iter()
                .map(|channel| (channel.coefficient, channel.coupling.as_slice())),
            mode,
        )
        .and_then(|(omega, beta)| duffing_spectrum_exact(omega, beta, amplitude, odd_harmonics));
        result.map_err(|error| NlModalError::Parameter {
            what: if error == DuffingBackboneError::NonInvariantAxis {
                "spectrum axis couples to other displacement coordinates"
            } else {
                "spectrum mode, amplitude, harmonic count or scalar spectrum is not representable"
            },
        })
    }

    /// Physical displacement, velocity and acceleration harmonics of one axis.
    ///
    /// This uses the same original coefficients and channel-by-channel
    /// invariance rule as [`Self::invariant_axis_spectrum`]. Derivative
    /// coefficients retain their physical phases and units; separate tail
    /// estimates expose the greater bandwidth required for velocity and force.
    /// Products are scaled before rounding, so a representable derivative can
    /// survive even when its displacement harmonic is unrepresentably small.
    /// An unrepresentable derivative coefficient or tail refuses this opt-in
    /// query without narrowing the displacement-only query's range.
    pub fn invariant_axis_kinematic_spectrum(
        &self,
        mode: usize,
        amplitude: f64,
        odd_harmonics: usize,
    ) -> Result<DuffingKinematicSpectrum, NlModalError> {
        self.validate()?;
        backbone::invariant_axis_parameters(
            &self.omegas,
            self.channels
                .iter()
                .map(|channel| (channel.coefficient, channel.coupling.as_slice())),
            mode,
        )
        .and_then(|(omega, beta)| {
            duffing_kinematic_spectrum_exact(omega, beta, amplitude, odd_harmonics)
        })
        .map_err(|error| NlModalError::Parameter {
            what: if error == DuffingBackboneError::NonInvariantAxis {
                "kinematic spectrum axis couples to other displacement coordinates"
            } else {
                "kinematic spectrum parameters, coefficients or tails are not representable"
            },
        })
    }

    /// Exact continuous effort Jacobian applied to a physical state direction.
    ///
    /// Both inputs and the output use the original interleaved
    /// `[q_0, p_0, q_1, p_1, ...]` mass-normalized coordinates. For a symmetric
    /// stress channel, differentiating `c (q^T E q) E q` gives
    /// `c [(q^T E q) E v + 2 (q^T E v) E q]`. The linear displacement rows add
    /// `omega^2 v_q`; momentum rows are exactly `v_p`. Off-diagonal modal
    /// interactions, negative channel entries and all stress channels remain.
    ///
    /// This differentiates the ORIGINAL [`Storage::gradient`], including any
    /// roundoff-sized coupling asymmetry accepted by [`Self::validate`]: the
    /// quadratic form's directional derivative is `v^T E q + q^T E v` rather
    /// than assuming these two represented contractions are identical. For an
    /// exactly symmetric channel this is the Hamiltonian Hessian action.
    ///
    /// No matrix, finite-difference state, or temporary vector is allocated.
    /// Two matrix-vector passes per channel keep memory constant and work
    /// `O(channel_count * mode_count^2)`. Traversal order is fixed. Original
    /// model admission is rechecked because its frequency/channel fields are
    /// public. [`Storage::continuous_tangent`] exposes this same action after
    /// boxing, composition and assembly, so periodic/variational solvers can
    /// discover it without a second model or a solver dependency here.
    ///
    /// This is a continuous constitutive tangent, not the endpoint Jacobian of
    /// a discrete-gradient time step. It establishes neither harmonic
    /// truncation accuracy nor nonlinear orbit stability by itself.
    ///
    /// # Errors
    /// Refuses invalid storage, incomplete/nonfinite inputs, a wrong output
    /// length, or unrepresentable intermediate/output arithmetic. Shape/input
    /// admission leaves `out` unchanged; after an arithmetic error its contents
    /// are unspecified.
    pub fn hessian_action(
        &self,
        state: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> Result<(), NlModalError> {
        self.validate()?;
        let n = self.n_modes();
        let dimension = n.checked_mul(2).ok_or(NlModalError::Parameter {
            what: "modal tangent dimension overflow",
        })?;
        if state.len() != dimension
            || direction.len() != dimension
            || out.len() != dimension
            || state.iter().chain(direction).any(|v| !v.is_finite())
        {
            return Err(NlModalError::Parameter {
                what: "modal tangent needs complete finite state/direction and matching output",
            });
        }
        for (k, omega) in self.omegas.iter().enumerate() {
            out[2 * k] = omega * omega * direction[2 * k];
            out[2 * k + 1] = direction[2 * k + 1];
        }
        for channel in &self.channels {
            let mut form = 0.0;
            let mut form_direction = 0.0;
            for k in 0..n {
                let mut eq = 0.0;
                let mut ev = 0.0;
                for j in 0..n {
                    let entry = channel.coupling[k * n + j];
                    eq += entry * state[2 * j];
                    ev += entry * direction[2 * j];
                }
                form += state[2 * k] * eq;
                form_direction += direction[2 * k] * eq + state[2 * k] * ev;
            }
            if !form.is_finite() || !form_direction.is_finite() {
                return Err(NlModalError::Parameter {
                    what: "modal tangent channel contraction overflow",
                });
            }
            for k in 0..n {
                let mut eq = 0.0;
                let mut ev = 0.0;
                for j in 0..n {
                    let entry = channel.coupling[k * n + j];
                    eq += entry * state[2 * j];
                    ev += entry * direction[2 * j];
                }
                out[2 * k] += channel.coefficient * (form * ev + form_direction * eq);
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(NlModalError::Parameter {
                what: "modal tangent output overflow",
            });
        }
        Ok(())
    }

    /// `q^T E q` for channel `ch` given the interleaved state.
    fn quadratic_form(&self, ch: &StressChannel, x: &[f64]) -> f64 {
        let n = self.n_modes();
        let mut acc = 0.0;
        for p in 0..n {
            let mut row = 0.0;
            for q in 0..n {
                row += ch.coupling[p * n + q] * x[2 * q];
            }
            acc += x[2 * p] * row;
        }
        acc
    }

    // The public source permits roundoff-sized asymmetry in E. Its coded
    // effort then differs from the exact energy gradient: retain Gonzalez for
    // that case instead of advertising a constitutive rule with a false
    // coincidence limit. Native string/plate constructors are exactly symmetric.
    fn has_exact_symmetric_channels(&self) -> bool {
        let n = self.n_modes();
        let Some(square) = n.checked_mul(n) else {
            return false;
        };
        self.channels.iter().all(|channel| {
            channel.coupling.len() == square
                && (0..n).all(|i| {
                    (0..i).all(|j| channel.coupling[i * n + j] == channel.coupling[j * n + i])
                })
        })
    }

    fn valid_discrete_endpoints(&self, a: &[f64], b: &[f64], output_len: usize) -> bool {
        self.n_modes().checked_mul(2).is_some_and(|dimension| {
            a.len() == dimension
                && b.len() == dimension
                && output_len == dimension
                && a.iter().chain(b).all(|value| value.is_finite())
        }) && self.validate().is_ok()
            && self.has_exact_symmetric_channels()
    }
}

impl ContinuousStorageTangent for SosModalStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        self.hessian_action(state, direction, out).is_ok()
    }
}

impl Storage for SosModalStorage {
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }

    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        self.has_exact_symmetric_channels()
            .then_some(self as &dyn ConstitutiveDiscreteGradient)
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], midpoint: &[f64]) -> Option<f64> {
        // Let d=b-a and o=(a-m)+(b-m), so the exact mean is m+o/2.
        // For a symmetric channel define D=d.E.d, O=o.E.o, C=d.E.o,
        // M=m.E.d, N=m.E.o, S=m.E.m. Its complete remainder is
        // c*((M+C/2)*(N+(D+O)/4)+S*C/2). At an exact midpoint this
        // reduces to c*D*M/4. Retain o: a represented midpoint can round
        // to an endpoint, where even quadratic storage has a remainder.
        let n = self.n_modes();
        let mut remainder = 0.0;
        for (k, omega) in self.omegas.iter().enumerate() {
            let dq = b[2 * k] - a[2 * k];
            let dp = b[2 * k + 1] - a[2 * k + 1];
            let oq = (a[2 * k] - midpoint[2 * k]) + (b[2 * k] - midpoint[2 * k]);
            let op = (a[2 * k + 1] - midpoint[2 * k + 1]) + (b[2 * k + 1] - midpoint[2 * k + 1]);
            remainder += 0.5 * (omega * omega * dq * oq + dp * op);
        }
        for ch in &self.channels {
            if ch.coefficient == 0.0 {
                continue;
            }
            let mut delta_form = 0.0;
            let mut midpoint_delta_form = 0.0;
            let mut offset_form = 0.0;
            let mut delta_offset_form = 0.0;
            let mut midpoint_offset_form = 0.0;
            let mut midpoint_form = 0.0;
            let mut asymmetry_work = 0.0;
            for k in 0..n {
                let dk = b[2 * k] - a[2 * k];
                let ok = (a[2 * k] - midpoint[2 * k]) + (b[2 * k] - midpoint[2 * k]);
                let diagonal = ch.coupling[k * n + k];
                delta_form += diagonal * dk * dk;
                midpoint_delta_form += diagonal * midpoint[2 * k] * dk;
                offset_form += diagonal * ok * ok;
                delta_offset_form += diagonal * dk * ok;
                midpoint_offset_form += diagonal * midpoint[2 * k] * ok;
                midpoint_form += diagonal * midpoint[2 * k] * midpoint[2 * k];
                for j in 0..k {
                    let dj = b[2 * j] - a[2 * j];
                    let oj = (a[2 * j] - midpoint[2 * j]) + (b[2 * j] - midpoint[2 * j]);
                    let forward = ch.coupling[k * n + j];
                    let reverse = ch.coupling[j * n + k];
                    let symmetric = f64::midpoint(forward, reverse);
                    delta_form += 2.0 * symmetric * dk * dj;
                    midpoint_delta_form +=
                        symmetric * (midpoint[2 * k] * dj + midpoint[2 * j] * dk);
                    offset_form += 2.0 * symmetric * ok * oj;
                    delta_offset_form += symmetric * (dk * oj + dj * ok);
                    midpoint_offset_form +=
                        symmetric * (midpoint[2 * k] * oj + midpoint[2 * j] * ok);
                    midpoint_form += 2.0 * symmetric * midpoint[2 * k] * midpoint[2 * j];
                    asymmetry_work +=
                        (forward - reverse) * (midpoint[2 * k] * dj - midpoint[2 * j] * dk);
                }
            }
            remainder += ch.coefficient
                * ((midpoint_delta_form + 0.5 * delta_offset_form)
                    * (midpoint_offset_form + 0.25 * (delta_form + offset_form))
                    + 0.5 * midpoint_form * delta_offset_form);
            // Admission allows roundoff-sized asymmetry. Account for the
            // coded gradient's corresponding mismatch without subtracting
            // nearly equal work terms; exactly symmetric channels skip it.
            if asymmetry_work != 0.0 {
                remainder +=
                    0.5 * ch.coefficient * self.quadratic_form(ch, midpoint) * asymmetry_work;
            }
        }
        Some(remainder)
    }

    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        let n = 2 * self.n_modes();
        let mut matrix = vec![0.0; n * n];
        for (k, w) in self.omegas.iter().enumerate() {
            matrix[2 * k * n + 2 * k] = w * w;
            matrix[(2 * k + 1) * n + 2 * k + 1] = 1.0;
        }
        Some(fs_phs::linear_exact::QuadraticPart {
            matrix,
            exact: self
                .channels
                .iter()
                .all(|channel| channel.coefficient == 0.0),
        })
    }

    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let n = self.n_modes();
        let mut h = 0.0;
        for k in 0..n {
            let (q, p) = (x[2 * k], x[2 * k + 1]);
            h += f64::midpoint(p * p, self.omegas[k] * self.omegas[k] * q * q);
        }
        for ch in &self.channels {
            let s = self.quadratic_form(ch, x);
            h += 0.25 * ch.coefficient * s * s;
        }
        h
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let n = self.n_modes();
        for k in 0..n {
            out[2 * k] = self.omegas[k] * self.omegas[k] * x[2 * k];
            out[2 * k + 1] = x[2 * k + 1];
        }
        for ch in &self.channels {
            let s = self.quadratic_form(ch, x);
            // d/dq_k [c/4 s^2] = c s (E q)_k (E symmetric).
            for k in 0..n {
                let mut row = 0.0;
                for q in 0..n {
                    row += ch.coupling[k * n + q] * x[2 * q];
                }
                out[2 * k] += ch.coefficient * s * row;
            }
        }
    }
}

/// Exact average of the original cubic effort along the endpoint segment.
/// Each stress channel keeps its own elastic work; no quartic correction is
/// applied to momentum, a channel's null directions, or unrelated storage.
/// Evaluation and the endpoint action use constant temporary memory and fixed
/// matrix traversal, without numerical quadrature or finite differences.
impl ConstitutiveDiscreteGradient for SosModalStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], midpoint: &[f64], out: &mut [f64]) {
        if !self.valid_discrete_endpoints(a, b, out.len())
            || midpoint.len() != a.len()
            || midpoint.iter().any(|value| !value.is_finite())
        {
            out.fill(f64::NAN);
            return;
        }
        if a == b {
            self.gradient(a, out);
            return;
        }
        let n = self.n_modes();
        for (k, omega) in self.omegas.iter().enumerate() {
            out[2 * k] = omega * omega * midpoint[2 * k];
            out[2 * k + 1] = midpoint[2 * k + 1];
        }
        for channel in &self.channels {
            if channel.coefficient == 0.0 {
                continue;
            }
            let mut form = 0.0;
            let mut delta_form = 0.0;
            let mut mixed_form = 0.0;
            for k in 0..n {
                let mut em = 0.0;
                let mut ed = 0.0;
                for j in 0..n {
                    let entry = channel.coupling[k * n + j];
                    em += entry * midpoint[2 * j];
                    ed += entry * (b[2 * j] - a[2 * j]);
                }
                form += midpoint[2 * k] * em;
                delta_form += (b[2 * k] - a[2 * k]) * ed;
                mixed_form += midpoint[2 * k] * ed;
            }
            // Integral_0^1 c (q(t).E.q(t)) E.q(t) dt, with q(t)=m+(t-1/2)d:
            // c*((m.E.m+d.E.d/12) E.m + (m.E.d/6) E.d).
            let mean_factor = channel.coefficient * (form + delta_form / 12.0);
            let delta_factor = channel.coefficient * (mixed_form / 6.0);
            for k in 0..n {
                let mut em = 0.0;
                let mut ed = 0.0;
                for j in 0..n {
                    let entry = channel.coupling[k * n + j];
                    em += entry * midpoint[2 * j];
                    ed += entry * (b[2 * j] - a[2 * j]);
                }
                out[2 * k] += mean_factor * em + delta_factor * ed;
            }
        }
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        if !self.valid_discrete_endpoints(a, b, out.len())
            || direction.len() != a.len()
            || direction.iter().any(|value| !value.is_finite())
        {
            return false;
        }
        let n = self.n_modes();
        for (k, omega) in self.omegas.iter().enumerate() {
            out[2 * k] = 0.5 * omega * omega * direction[2 * k];
            out[2 * k + 1] = 0.5 * direction[2 * k + 1];
        }
        for channel in &self.channels {
            if channel.coefficient == 0.0 {
                continue;
            }
            let mut form = 0.0;
            let mut delta_form = 0.0;
            let mut mixed_form = 0.0;
            let mut mean_direction_form = 0.0;
            let mut delta_direction_form = 0.0;
            for k in 0..n {
                let mut em = 0.0;
                let mut ed = 0.0;
                let mut ev = 0.0;
                for j in 0..n {
                    let entry = channel.coupling[k * n + j];
                    em += entry * f64::midpoint(a[2 * j], b[2 * j]);
                    ed += entry * (b[2 * j] - a[2 * j]);
                    ev += entry * direction[2 * j];
                }
                let mk = f64::midpoint(a[2 * k], b[2 * k]);
                let dk = b[2 * k] - a[2 * k];
                form += mk * em;
                delta_form += dk * ed;
                mixed_form += mk * ed;
                mean_direction_form += mk * ev;
                delta_direction_form += dk * ev;
            }
            // D_b of the same polynomial, with D_b m[v]=v/2 and D_b d[v]=v.
            // This is also Integral_0^1 t Hessian(V)(a+t*d)[v] dt;
            // it is not the continuous tangent at a guessed midpoint.
            let direction_factor = 0.5 * (form + delta_form / 12.0) + mixed_form / 6.0;
            let mean_factor = mean_direction_form + delta_direction_form / 6.0;
            let delta_factor = mean_direction_form / 6.0 + delta_direction_form / 12.0;
            for k in 0..n {
                let mut em = 0.0;
                let mut ed = 0.0;
                let mut ev = 0.0;
                for j in 0..n {
                    let entry = channel.coupling[k * n + j];
                    em += entry * f64::midpoint(a[2 * j], b[2 * j]);
                    ed += entry * (b[2 * j] - a[2 * j]);
                    ev += entry * direction[2 * j];
                }
                out[2 * k] += channel.coefficient
                    * (direction_factor * ev + mean_factor * em + delta_factor * ed);
            }
        }
        out.iter().all(|value| value.is_finite())
    }
}

const SYM_TOL: f64 = 1.0e-10;

/// Assemble the pHS: symplectic pair blocks, per-mode viscous damping
/// `R = diag(0, 2 zeta_k w_k)` (zeta from the caller — the
/// visco-damping facility's per-mode output slots in here; this crate
/// does not invent a second damping representation), and a single
/// strike port weighted by the mode shapes at the strike point.
///
/// # Errors
/// [`NlModalError`] on length mismatches, non-positive frequencies,
/// negative damping, asymmetric couplings, or pHS admission failure.
pub fn assemble(
    storage: SosModalStorage,
    zetas: &[f64],
    strike_weights: &[f64],
) -> Result<PortHamiltonian, NlModalError> {
    let n = storage.n_modes();
    if zetas.len() != n || strike_weights.len() != n {
        return Err(NlModalError::Parameter {
            what: "zetas/strike_weights length vs mode count",
        });
    }
    storage.validate()?;
    for &z in zetas {
        if z.is_nan() || z < 0.0 {
            return Err(NlModalError::Parameter {
                what: "damping ratio must be non-negative",
            });
        }
    }
    let mut g = vec![0.0; 2 * n];
    for k in 0..n {
        g[2 * k + 1] = strike_weights[k];
    }
    let omegas = storage.omegas.clone();
    assemble_storage(n, &omegas, zetas, 1, g, Box::new(storage))
}

/// Assemble a modal pHS around already-built storage.
///
/// Same symplectic `J` and viscous `R` as [`assemble`], but the
/// storage may already wrap contact, a reduction, or any other
/// `fs-phs::Storage`. `g` is row-major `(2 n_modes) × m_ports`.
///
/// # Errors
/// [`NlModalError`] on length mismatches, non-physical frequencies
/// or damping, or pHS admission failure.
pub fn assemble_storage(
    n_modes: usize,
    omegas: &[f64],
    zetas: &[f64],
    m_ports: usize,
    g: Vec<f64>,
    storage: Box<dyn Storage>,
) -> Result<PortHamiltonian, NlModalError> {
    assemble_storage_damping(
        n_modes,
        omegas,
        ModalDamping::Ratios(zetas),
        m_ports,
        g,
        storage,
    )
}

/// Assemble a modal pHS with the complete projected viscous operator.
///
/// `modal_damping_per_s` is the symmetric positive-semidefinite, row-major
/// `n_modes × n_modes` matrix `C` in the SAME mass-normalized basis as
/// `storage`, `omegas` and `g`. The equations are `q_dot = p` and
/// `p_dot = -dV/dq - C p + G_p u`, with instantaneous viscous power
/// `p^T C p`. Local material or attachment loss generally couples modes,
/// including modes with nearly equal frequencies; replacing C by its
/// diagonal changes both the transient and the undamped motions.
///
/// Every supplied entry is retained in the momentum block of R. Rank-deficient
/// damping, zero rows and coupled lossless directions require no artificial
/// loss floor. The existing fs-phs owner performs its symmetry/PSD admission
/// and time integration. C already contains ALL caller-owned viscous loss;
/// no extra per-mode ratios are added here. Constitutive memory still belongs
/// to an explicit relaxation-state model, not this instantaneous operator.
///
/// # Errors
/// [`NlModalError`] for malformed/nonfinite input, nonpositive frequencies,
/// or the existing pHS symmetry/positive-semidefiniteness refusal.
pub fn assemble_storage_with_damping(
    n_modes: usize,
    omegas: &[f64],
    modal_damping_per_s: &[f64],
    m_ports: usize,
    g: Vec<f64>,
    storage: Box<dyn Storage>,
) -> Result<PortHamiltonian, NlModalError> {
    assemble_storage_damping(
        n_modes,
        omegas,
        ModalDamping::Matrix(modal_damping_per_s),
        m_ports,
        g,
        storage,
    )
}

enum ModalDamping<'a> {
    Ratios(&'a [f64]),
    Matrix(&'a [f64]),
}

fn assemble_storage_damping(
    n_modes: usize,
    omegas: &[f64],
    damping: ModalDamping<'_>,
    m_ports: usize,
    g: Vec<f64>,
    storage: Box<dyn Storage>,
) -> Result<PortHamiltonian, NlModalError> {
    let Some(dim) = n_modes.checked_mul(2) else {
        return Err(NlModalError::Parameter {
            what: "modal state dimension overflow",
        });
    };
    let Some(square) = dim.checked_mul(dim) else {
        return Err(NlModalError::Parameter {
            what: "modal matrix dimension overflow",
        });
    };
    if n_modes == 0 || omegas.len() != n_modes {
        return Err(NlModalError::Parameter {
            what: "positive mode count and matching omegas required",
        });
    }
    if m_ports == 0 || dim.checked_mul(m_ports) != Some(g.len()) {
        return Err(NlModalError::Parameter {
            what: "G shape vs (2 n_modes) × m_ports",
        });
    }
    if g.iter().any(|value| !value.is_finite()) {
        return Err(NlModalError::Parameter {
            what: "port projections must be finite",
        });
    }
    for &w in omegas {
        if w <= 0.0 || !w.is_finite() {
            return Err(NlModalError::Parameter {
                what: "modal frequency must be positive and finite",
            });
        }
    }
    match damping {
        ModalDamping::Ratios(zetas) => {
            if zetas.len() != n_modes {
                return Err(NlModalError::Parameter {
                    what: "zetas length vs mode count",
                });
            }
            if zetas.iter().zip(omegas).any(|(&zeta, &omega)| {
                zeta < 0.0 || !zeta.is_finite() || !(2.0 * zeta * omega).is_finite()
            }) {
                return Err(NlModalError::Parameter {
                    what: "damping ratios and rates must be non-negative and finite",
                });
            }
        }
        ModalDamping::Matrix(matrix) => {
            if matrix.len() != n_modes * n_modes {
                return Err(NlModalError::Parameter {
                    what: "damping matrix size vs mode count",
                });
            }
            if matrix.iter().any(|value| !value.is_finite()) {
                return Err(NlModalError::Parameter {
                    what: "damping matrix entries must be finite",
                });
            }
        }
    }
    let mut j = vec![0.0; square];
    let mut r = vec![0.0; square];
    for k in 0..n_modes {
        let (qi, pi) = (2 * k, 2 * k + 1);
        j[qi * dim + pi] = 1.0;
        j[pi * dim + qi] = -1.0;
        match damping {
            ModalDamping::Ratios(zetas) => r[pi * dim + pi] = 2.0 * zetas[k] * omegas[k],
            ModalDamping::Matrix(matrix) => {
                for l in 0..n_modes {
                    r[pi * dim + 2 * l + 1] = matrix[k * n_modes + l];
                }
            }
        }
    }
    PortHamiltonian::new(dim, m_ports, j, r, g, storage).map_err(NlModalError::Phs)
}

// ---------------------------------------------------------------------
// Von Karman simply-supported rectangular plate
// ---------------------------------------------------------------------

/// Von Karman SS plate parameters (isotropic, in-plane movable edges).
#[derive(Debug, Clone, Copy)]
pub struct VkPlateParams {
    /// Side lengths [m].
    pub lx: f64,
    /// Side lengths [m].
    pub ly: f64,
    /// Thickness [m].
    pub h: f64,
    /// Young's modulus [Pa].
    pub young: f64,
    /// Poisson ratio.
    pub nu: f64,
    /// Density [kg/m^3].
    pub rho: f64,
    /// Isotropic in-plane pretension [N/m]. Zero is the unloaded
    /// plate. Raises linear frequencies via `ω² = (D k⁴ + T k²)/(ρ h)`.
    pub pretension_n_m: f64,
}

/// One sine mode index pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SineMode {
    /// x half-waves.
    pub m: usize,
    /// y half-waves.
    pub n: usize,
}

/// The assembled von Karman modal model plus its certification data.
pub struct VkModel {
    /// The storage (frequencies + channels), ready for [`assemble`].
    pub storage: SosModalStorage,
    /// Displacement mode list (index-aligned with the storage).
    pub disp_modes: Vec<SineMode>,
    /// Retained Airy stress basis. Tensor membrane elasticity can mix these
    /// modes when factoring the positive complementary energy into channels.
    pub stress_modes: Vec<SineMode>,
    /// Worst relative disagreement between the two spatial approximations
    /// for sampled modes. Zero for [`von_karman_ss_plate`], whose integrals
    /// are analytic and do not have a quadrature truncation error. Zero does
    /// not certify the separate displacement/stress modal truncations.
    pub quadrature_residual: f64,
}

/// Gauss-Legendre nodes/weights on [0, 1] (order-n rule; nodes by
/// Newton on the Legendre recurrence — deterministic).
#[cfg(test)]
fn gauss_01(order: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = Vec::with_capacity(order);
    let mut weights = Vec::with_capacity(order);
    let nf = order as f64;
    for i in 0..order {
        // Initial guess (Chebyshev-like), then Newton.
        let mut x = det::cos(core::f64::consts::PI * (i as f64 + 0.75) / (nf + 0.5));
        for _ in 0..60 {
            // Legendre P_n(x) and derivative by recurrence.
            let (mut p0, mut p1) = (1.0f64, x);
            for k in 2..=order {
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            let dp = nf * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1.0e-15 {
                break;
            }
        }
        // Recompute derivative at the root for the weight.
        let (mut p0, mut p1) = (1.0f64, x);
        for k in 2..=order {
            let kf = k as f64;
            let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
            p0 = p1;
            p1 = p2;
        }
        let dp = nf * (x * p1 - p0) / (x * x - 1.0);
        nodes.push(f64::midpoint(x, 1.0));
        weights.push(1.0 / ((1.0 - x * x) * dp * dp));
    }
    (nodes, weights)
}

/// The von Karman bracket integrand `Psi_j L(Phi_p, Phi_q)` for sine
/// modes on the SS rectangle, integrated by the given quadrature.
/// `L(f, g) = f_xx g_yy + f_yy g_xx - 2 f_xy g_xy`. Shapes are UNIT
/// sine products; normalization factors are applied by the caller.
#[cfg(test)]
fn coupling_integral(
    lx: f64,
    ly: f64,
    stress: SineMode,
    a: SineMode,
    b: SineMode,
    nodes: &[f64],
    weights: &[f64],
) -> f64 {
    // Reflection about either midline changes the integrand's sign when the
    // three half-wave indices sum to an even number. Enforce that exact
    // selection rule instead of injecting roundoff-sized forbidden channels.
    if (stress.m ^ a.m ^ b.m) & 1 == 0 || (stress.n ^ a.n ^ b.n) & 1 == 0 {
        return 0.0;
    }
    // Sine products make the Hessian bracket separable. With S=sin(s*pi*x),
    // the two needed integrals are integral S*sin(a*pi*x)*sin(b*pi*x) and
    // integral S*cos(a*pi*x)*cos(b*pi*x), on the unit interval. Each tensor
    // quadrature therefore costs O(Q), rather than evaluating a Q-by-Q grid.
    // Both independent quadrature orders remain owned by the constructor.
    let axis = |s: usize, a: usize, b: usize| {
        let ks = s as f64 * core::f64::consts::PI;
        let ka = a as f64 * core::f64::consts::PI;
        let kb = b as f64 * core::f64::consts::PI;
        let (mut sine, mut cosine) = (0.0, 0.0);
        for (&x, &weight) in nodes.iter().zip(weights) {
            let (sa, ca) = det::sin_cos(ka * x);
            let (sb, cb) = det::sin_cos(kb * x);
            let weighted_stress = weight * det::sin(ks * x);
            sine += weighted_stress * sa * sb;
            cosine += weighted_stress * ca * cb;
        }
        (sine, cosine)
    };
    let (sine_x, cosine_x) = axis(stress.m, a.m, b.m);
    let (sine_y, cosine_y) = axis(stress.n, a.n, b.n);
    let kx = |m: usize| m as f64 * core::f64::consts::PI / lx;
    let ky = |n: usize| n as f64 * core::f64::consts::PI / ly;
    let (pax, pay) = (kx(a.m), ky(a.n));
    let (pbx, pby) = (kx(b.m), ky(b.n));
    let normal = (pax * pax * pby * pby + pay * pay * pbx * pbx) * sine_x * sine_y;
    let shear = 2.0 * pax * pay * pbx * pby * cosine_x * cosine_y;
    (normal - shear) * lx * ly
}

/// Exact integral of three integer sine modes on [0, 1]. Product-to-sum
/// gives -4*s*a*b / (pi*(s+a+b)*(s+a-b)*(s-a+b)*(s-a-b)) when their
/// sum is odd, and zero otherwise. Keeping the signed differences in
/// integer arithmetic avoids losing a near-resonant difference at high
/// mode numbers; the factored expression avoids cancellation of four
/// nearly equal reciprocal frequencies for well-separated modes.
fn triple_sine_integral(s: usize, a: usize, b: usize) -> f64 {
    if (s ^ a ^ b) & 1 == 0 {
        return 0.0;
    }
    let mut indices = [s, a, b];
    indices.sort_unstable();
    let [s, a, b] = indices.map(|index| index as i128);
    let denominator =
        (s + a + b) as f64 * (s + a - b) as f64 * (s - a + b) as f64 * (s - a - b) as f64;
    -4.0 * (s as f64) * (a as f64) * (b as f64) / (core::f64::consts::PI * denominator)
}

/// Exact separated Hessian bracket. In one axis,
/// integral sin(s)cos(a)cos(b) / integral sin(s)sin(a)sin(b)
/// is (a^2+b^2-s^2)/(2*a*b). Substitution cancels the four modal
/// wavenumbers in the shear term. The remaining integer polynomial is
/// evaluated exactly whenever i128 can hold it, retaining weak couplings
/// that would lose precision when the normal and shear terms nearly cancel.
fn analytic_coupling_integral(lx: f64, ly: f64, s: SineMode, a: SineMode, b: SineMode) -> f64 {
    if (s.m ^ a.m ^ b.m) & 1 == 0 || (s.n ^ a.n ^ b.n) & 1 == 0 {
        return 0.0;
    }
    let square = |x: usize| (x as i128) * (x as i128);
    let (ax, ay, bx, by, sx, sy) = (
        square(a.m),
        square(a.n),
        square(b.m),
        square(b.n),
        square(s.m),
        square(s.n),
    );
    let (dx, dy) = (ax + bx - sx, ay + by - sy);
    let numerator = ax
        .checked_mul(by)
        .and_then(|v| ay.checked_mul(bx).and_then(|w| v.checked_add(w)))
        .and_then(|v| v.checked_mul(2))
        .and_then(|v| dx.checked_mul(dy).and_then(|w| v.checked_sub(w)))
        .map_or_else(
            || {
                2.0 * ((ax as f64) * (by as f64) + (ay as f64) * (bx as f64))
                    - (dx as f64) * (dy as f64)
            },
            |v| v as f64,
        );
    let pi = core::f64::consts::PI;
    let scale = (pi / lx) * (pi / ly) * pi * pi;
    0.5 * numerator
        * triple_sine_integral(s.m, a.m, b.m)
        * triple_sine_integral(s.n, a.n, b.n)
        * scale
}

/// Parameter/mode-list validation for [`von_karman_ss_plate`].
fn validate_vk_inputs(
    params: &VkPlateParams,
    disp_modes: &[SineMode],
    stress_modes: &[SineMode],
) -> Result<(), NlModalError> {
    let VkPlateParams {
        lx,
        ly,
        h,
        young,
        nu,
        rho,
        pretension_n_m,
    } = *params;
    if !(lx > 0.0
        && ly > 0.0
        && h > 0.0
        && young > 0.0
        && rho > 0.0
        && pretension_n_m >= 0.0
        && (0.0..0.5).contains(&nu))
    {
        return Err(NlModalError::Parameter {
            what: "plate parameters",
        });
    }
    if [lx, ly, h, young, nu, rho, pretension_n_m]
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(NlModalError::Parameter {
            what: "finite plate parameters",
        });
    }
    validate_vk_modes(disp_modes, stress_modes)
}

fn validate_vk_modes(
    disp_modes: &[SineMode],
    stress_modes: &[SineMode],
) -> Result<(), NlModalError> {
    if disp_modes.is_empty() || stress_modes.is_empty() {
        return Err(NlModalError::Parameter {
            what: "empty mode list",
        });
    }
    for md in disp_modes.iter().chain(stress_modes) {
        if md.m == 0 || md.n == 0 {
            return Err(NlModalError::Parameter {
                what: "sine mode indices start at 1",
            });
        }
        // Exact integer-to-f64 conversion of each half-wave index is needed
        // by the frequencies and analytic overlap. This bound also keeps
        // all sums of index squares representable in the i128 calculation.
        if md.m as u128 > (1_u128 << 53) || md.n as u128 > (1_u128 << 53) {
            return Err(NlModalError::Parameter {
                what: "sine mode indices must be exactly representable in f64",
            });
        }
    }
    // Duplicates silently double-count a physical mode: refuse in each list.
    // Sort a cold copy, retaining the authored basis order. A quadratic scan
    // would dominate analytic assembly for many Airy modes and few motions.
    for list in [disp_modes, stress_modes] {
        let mut ordered = list.to_vec();
        ordered.sort_unstable_by_key(|mode| (mode.m, mode.n));
        if ordered.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(NlModalError::Parameter {
                what: "duplicate mode in list",
            });
        }
    }
    Ok(())
}

/// Build the von Karman modal model for a simply-supported rectangle.
///
/// Displacement modes are MASS-NORMALIZED sine products (`rho h
/// integral Phi^2 = 1`); stress modes are UNIT-normalized sine
/// products (`integral Psi^2 = 1`), both biharmonic eigenfunctions of
/// the SS rectangle. Channel coefficients follow the Airy elimination
/// `c_j = E h / (2 xi_j^4)` with `xi_j^4 = ((m pi/lx)^2 +
/// (n pi/ly)^2)^2`.
///
/// The stress-mode count is a SEPARATE truncation from the
/// displacement count (the literature uses more stress modes); both
/// lists are explicit inputs. Integer trigonometric products give the
/// coupling tensor analytically, including exact reflection-parity zeros.
/// Assembly costs O(stress_modes * disp_modes²), independently of the
/// largest half-wave index; it does not sample and alias spatial harmonics.
/// The separate modal truncations still need physical convergence checks.
///
/// # Errors
/// [`NlModalError`] on bad parameters, unrepresentable mode indices, or
/// nonfinite assembled storage.
pub fn von_karman_ss_plate(
    params: &VkPlateParams,
    disp_modes: &[SineMode],
    stress_modes: &[SineMode],
) -> Result<VkModel, NlModalError> {
    let VkPlateParams {
        lx,
        ly,
        h,
        young,
        nu,
        rho,
        pretension_n_m,
    } = *params;
    validate_vk_inputs(params, disp_modes, stress_modes)?;
    let d_bend = young * h * h * h / (12.0 * (1.0 - nu * nu));
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = disp_modes
        .iter()
        .map(|md| {
            let k2 = (md.m as f64 * pi / lx) * (md.m as f64 * pi / lx)
                + (md.n as f64 * pi / ly) * (md.n as f64 * pi / ly);
            det::sqrt((d_bend * k2 * k2 + pretension_n_m * k2) / (rho * h))
        })
        .collect();
    // Normalizations: Phi = phi_norm * sin sin with rho h * (lx ly /4)
    // * phi_norm^2 = 1; Psi = psi_norm * sin sin with (lx ly/4) *
    // psi_norm^2 = 1.
    let phi_norm = det::sqrt(4.0 / (rho * h * lx * ly));
    let psi_norm = det::sqrt(4.0 / (lx * ly));
    let nq = disp_modes.len();
    let matrix_len = nq.checked_mul(nq).ok_or(NlModalError::Parameter {
        what: "plate coupling matrix dimensions",
    })?;
    let mut channels = Vec::with_capacity(stress_modes.len());
    for &sm in stress_modes {
        let xi4 = {
            let k2 = (sm.m as f64 * pi / lx) * (sm.m as f64 * pi / lx)
                + (sm.n as f64 * pi / ly) * (sm.n as f64 * pi / ly);
            k2 * k2
        };
        let mut e = vec![0.0; matrix_len];
        for p in 0..nq {
            for q in 0..=p {
                let raw = analytic_coupling_integral(lx, ly, sm, disp_modes[p], disp_modes[q]);
                let value = psi_norm * phi_norm * phi_norm * raw;
                e[p * nq + q] = value;
                e[q * nq + p] = value;
            }
        }
        channels.push(StressChannel {
            coefficient: young * h / (2.0 * xi4),
            coupling: e,
        });
    }
    let storage = SosModalStorage { omegas, channels };
    storage.validate()?;
    Ok(VkModel {
        storage,
        disp_modes: disp_modes.to_vec(),
        stress_modes: stress_modes.to_vec(),
        quadrature_residual: 0.0,
    })
}

/// A displacement (or Airy) field sampled on a regular `[0, Lx] × [0, Ly]`
/// grid. FE eigenpairs, clamped shapes, and orthotropic modes enter
/// here — the constructor does not assume sine products.
#[derive(Debug, Clone)]
pub struct SampledPlateMode {
    /// Linear angular frequency [rad/s] (from the bending pencil).
    pub omega: f64,
    /// Samples of `w` (or `ψ`), `ny * nx`, x-fast.
    pub w: Vec<f64>,
    /// Sample count along x (must be ≥ 5).
    pub nx: usize,
    /// Sample count along y (must be ≥ 5).
    pub ny: usize,
}

/// Von Karman storage from *sampled* displacement modes.
///
/// Second derivatives use central finite differences at every interior
/// grid point. The Airy sine trace is zero at the physical boundary, so
/// boundary Hessians are not needed. Prepared displacement Hessians are
/// reused across all mode pairs and stress channels.
/// Airy stress uses the same sine products as the SS constructor
/// (in-plane movable), so clamped/orthotropic bending is FE while
/// the membrane channel stays the analytic Airy basis. This constructor
/// uses isotropic membrane elasticity. Use
/// [`von_karman_sampled_plate_with_membrane`] for directional stiffness,
/// independent shear stiffness, or a rotated material axis.
///
/// # Errors
/// [`NlModalError`] on bad grids, non-physical parameters, Airy indices
/// at or above the grid Nyquist limit, or an
/// FD-stencil disagreement above `0.9` of the channel scale
/// (FE Hessians are first-order; the sine reprint test is the
/// accuracy pin).
pub fn von_karman_sampled_plate(
    params: &VkPlateParams,
    disp: &[SampledPlateMode],
    stress_modes: &[SineMode],
) -> Result<VkModel, NlModalError> {
    validate_vk_inputs(params, &dummy_disp_ids(disp.len()), stress_modes)?;
    let mut model = sampled_plate_couplings(
        params.lx,
        params.ly,
        params.rho * params.h,
        disp,
        stress_modes,
    )?;
    for (channel, mode) in model.storage.channels.iter_mut().zip(stress_modes) {
        let k2 = (mode.m as f64 * core::f64::consts::PI / params.lx).powi(2)
            + (mode.n as f64 * core::f64::consts::PI / params.ly).powi(2);
        channel.coefficient = params.young * params.h / (2.0 * (k2 * k2).max(1.0e-30));
    }
    model.storage.validate()?;
    Ok(model)
}

// Assemble geometry once; constitutive callers supply the complementary
// membrane metric. Unit coefficients here are replaced before publication.
#[allow(clippy::items_after_statements)]
fn sampled_plate_couplings(
    lx: f64,
    ly: f64,
    areal_density: f64,
    disp: &[SampledPlateMode],
    stress_modes: &[SineMode],
) -> Result<VkModel, NlModalError> {
    if [lx, ly, areal_density]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(NlModalError::Parameter {
            what: "finite positive sampled plate dimensions and areal density",
        });
    }
    validate_vk_modes(&dummy_disp_ids(disp.len()), stress_modes)?;
    if disp.is_empty() {
        return Err(NlModalError::Parameter {
            what: "empty sampled displacement list",
        });
    }
    let nx = disp[0].nx;
    let ny = disp[0].ny;
    if nx < 5 || ny < 5 {
        return Err(NlModalError::Parameter {
            what: "sampled grid must be at least 5×5",
        });
    }
    let Some(points) = nx.checked_mul(ny) else {
        return Err(NlModalError::Parameter {
            what: "sampled grid size overflow",
        });
    };
    for mode in disp {
        if mode.nx != nx
            || mode.ny != ny
            || mode.w.len() != points
            || !(mode.omega > 0.0 && mode.omega.is_finite())
            || mode.w.iter().any(|v| !v.is_finite())
        {
            return Err(NlModalError::Parameter {
                what: "sampled modes must share a finite positive grid",
            });
        }
    }
    if stress_modes
        .iter()
        .any(|mode| mode.m >= nx - 1 || mode.n >= ny - 1)
    {
        return Err(NlModalError::Parameter {
            what: "sampled Airy mode must be below the grid Nyquist limit",
        });
    }
    let dx = lx / (nx - 1) as f64;
    let dy = ly / (ny - 1) as f64;
    let area_el = dx * dy;
    let mut hessians = Vec::with_capacity(disp.len());
    let mut omegas = Vec::with_capacity(disp.len());
    for mode in disp {
        let energy: f64 = mode.w.iter().map(|v| v * v).sum::<f64>() * area_el * areal_density;
        if !(energy > 0.0 && energy.is_finite()) {
            return Err(NlModalError::Parameter {
                what: "sampled mode must have finite positive L2 mass",
            });
        }
        let scale = det::sqrt(1.0 / energy);
        let field: Vec<_> = mode.w.iter().map(|v| v * scale).collect();
        hessians.push(sampled_hessians(&field, nx, ny, dx, dy));
        omegas.push(mode.omega);
    }
    type RawChannel = (Vec<(usize, usize, f64, f64)>, f64);
    let nq = hessians.len();
    let pi = core::f64::consts::PI;
    let psi_norm = det::sqrt(4.0 / (lx * ly));
    let kmin2 = (pi / lx) * (pi / lx) + (pi / ly) * (pi / ly);
    // Mass-normalized φ ~ 1/√(ρ h A); a nonzero raw integral is ~ k⁴/(ρ h).
    let char_scale = kmin2 * kmin2 / areal_density.max(1.0e-30);
    let mut channels = Vec::with_capacity(stress_modes.len());
    let mut worst_rel = 0.0f64;
    let mut global_scale = f64::MIN_POSITIVE;
    let mut raw: Vec<RawChannel> = Vec::new();
    for &sm in stress_modes {
        let psi = sample_sine(sm, lx, ly, nx, ny);
        let mut pairs = Vec::new();
        let mut scale = f64::MIN_POSITIVE;
        for p in 0..nq {
            for q in 0..=p {
                let (a2, a4) = sampled_coupling(&psi, &hessians[p], &hessians[q], area_el);
                scale = scale.max(a2.abs()).max(a4.abs());
                pairs.push((p, q, a2, a4));
            }
        }
        global_scale = global_scale.max(scale);
        raw.push((pairs, scale));
    }
    for (pairs, scale) in raw {
        let judge = scale.max(1.0e-12 * global_scale).max(1.0e-8 * char_scale);
        let mut e = vec![0.0; nq * nq];
        for (p, q, a2, a4) in pairs {
            worst_rel = worst_rel.max((a2 - a4).abs() / judge);
            // 2nd-order FD is the production value: DKT/FE w is only C⁰,
            // so the 4th-order stencil is not more accurate there. The
            // 4th-order residual is a convergence probe, not a Gauss
            // certificate.
            let v = psi_norm * a2;
            e[p * nq + q] = v;
            e[q * nq + p] = v;
        }
        channels.push(StressChannel {
            coefficient: 1.0,
            coupling: e,
        });
    }
    if worst_rel > 0.9 {
        return Err(NlModalError::QuadratureMismatch {
            residual: worst_rel,
        });
    }
    let storage = SosModalStorage { omegas, channels };
    storage.validate()?;
    Ok(VkModel {
        storage,
        disp_modes: dummy_disp_ids(nq),
        stress_modes: stress_modes.to_vec(),
        quadrature_residual: worst_rel,
    })
}

fn dummy_disp_ids(n: usize) -> Vec<SineMode> {
    (1..=n).map(|m| SineMode { m, n: 1 }).collect()
}

fn sample_sine(mode: SineMode, lx: f64, ly: f64, nx: usize, ny: usize) -> Vec<f64> {
    let mut out = vec![0.0; nx * ny];
    let pi = core::f64::consts::PI;
    for j in 0..ny {
        let y = ly * j as f64 / (ny - 1) as f64;
        for i in 0..nx {
            let x = lx * i as f64 / (nx - 1) as f64;
            out[j * nx + i] =
                det::sin(mode.m as f64 * pi * x / lx) * det::sin(mode.n as f64 * pi * y / ly);
        }
    }
    out
}

fn fd_d2(w: &[f64], nx: usize, i: usize, j: usize, axis: usize, h: f64, order: u8) -> f64 {
    let at = |ii: i32, jj: i32| -> f64 {
        let ii = ii.clamp(0, nx as i32 - 1) as usize;
        let jj = jj.clamp(0, (w.len() / nx) as i32 - 1) as usize;
        w[jj * nx + ii]
    };
    let (ii, jj) = (i as i32, j as i32);
    let (di, dj) = if axis == 0 { (1, 0) } else { (0, 1) };
    if order >= 4 && i >= 2 && i + 2 < nx && j >= 2 && j + 2 < w.len() / nx {
        (-at(ii + 2 * di, jj + 2 * dj) + 16.0 * at(ii + di, jj + dj) - 30.0 * at(ii, jj)
            + 16.0 * at(ii - di, jj - dj)
            - at(ii - 2 * di, jj - 2 * dj))
            / (12.0 * h * h)
    } else {
        (at(ii + di, jj + dj) - 2.0 * at(ii, jj) + at(ii - di, jj - dj)) / (h * h)
    }
}

fn fd_dxy(w: &[f64], nx: usize, i: usize, j: usize, dx: f64, dy: f64) -> f64 {
    let ny = w.len() / nx;
    let at = |ii: i32, jj: i32| -> f64 {
        let ii = ii.clamp(0, nx as i32 - 1) as usize;
        let jj = jj.clamp(0, ny as i32 - 1) as usize;
        w[jj * nx + ii]
    };
    let (ii, jj) = (i as i32, j as i32);
    (at(ii + 1, jj + 1) - at(ii + 1, jj - 1) - at(ii - 1, jj + 1) + at(ii - 1, jj - 1))
        / (4.0 * dx * dy)
}

// [xx(second order), yy(second order), xy, xx(fourth order), yy(fourth order)].
// The mixed derivative is shared by the production and diagnostic stencils.
fn sampled_hessians(w: &[f64], nx: usize, ny: usize, dx: f64, dy: f64) -> Vec<[f64; 5]> {
    let mut hessians = vec![[0.0; 5]; nx * ny];
    // Every interior node has a valid second-order central stencil,
    // including the first ring next to the boundary. Omitting that ring
    // biases the harmonic overlaps, particularly the mixed-curvature term.
    // Use second-order differences there for both orders; the Airy trace
    // is zero on the boundary itself, where the Hessian is left zero.
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            hessians[j * nx + i] = [
                fd_d2(w, nx, i, j, 0, dx, 2),
                fd_d2(w, nx, i, j, 1, dy, 2),
                fd_dxy(w, nx, i, j, dx, dy),
                fd_d2(w, nx, i, j, 0, dx, 4),
                fd_d2(w, nx, i, j, 1, dy, 4),
            ];
        }
    }
    hessians
}

fn sampled_coupling(psi: &[f64], a: &[[f64; 5]], b: &[[f64; 5]], area_el: f64) -> (f64, f64) {
    let mut a2 = 0.0;
    let mut a4 = 0.0;
    for ((&psi, a), b) in psi.iter().zip(a).zip(b) {
        let mixed = 2.0 * a[2] * b[2];
        a2 += psi * (a[0] * b[1] + a[1] * b[0] - mixed);
        a4 += psi * (a[3] * b[4] + a[4] * b[3] - mixed);
    }
    (a2 * area_el, a4 * area_el)
}

// ---------------------------------------------------------------------
// Kirchhoff-Carrier string
// ---------------------------------------------------------------------

/// Kirchhoff-Carrier string parameters.
#[derive(Debug, Clone, Copy)]
pub struct KcStringParams {
    /// Speaking length [m].
    pub length: f64,
    /// Tension [N].
    pub tension: f64,
    /// Linear density [kg/m].
    pub lin_density: f64,
    /// Longitudinal stiffness `E A` [N].
    pub ea: f64,
}

/// Kirchhoff-Carrier modal model: `n_modes` sine string modes with the
/// single tension-modulation channel
/// `H_nl = (E A / (8 L)) * (sum_k (k pi / L)^2 Q_k^2)^2` where `Q_k`
/// is the PHYSICAL modal displacement; in mass-normalized coordinates
/// (`q_k = sqrt(mu L / 2) Q_k`) the channel matrix is
/// `E[k,k] = (k pi / L)^2 * 2/(mu L)` with coefficient `E A L / 8`
/// — the exact Kirchhoff-Carrier averaged-tension form.
///
/// # Errors
/// [`NlModalError::Parameter`] on non-physical inputs.
pub fn kirchhoff_carrier_string(
    params: &KcStringParams,
    n_modes: usize,
) -> Result<SosModalStorage, NlModalError> {
    let KcStringParams {
        length,
        tension,
        lin_density,
        ea,
    } = *params;
    if !(length > 0.0 && tension > 0.0 && lin_density > 0.0 && ea >= 0.0) || n_modes == 0 {
        return Err(NlModalError::Parameter {
            what: "string parameters",
        });
    }
    let c = det::sqrt(tension / lin_density);
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = (1..=n_modes).map(|k| k as f64 * pi * c / length).collect();
    // Tension modulation: T(t) = T0 + (EA/2L) * integral (w')^2 dx;
    // with w = sum Q_k sin(k pi x / L): integral (w')^2 = (L/2) sum
    // (k pi/L)^2 Q_k^2. The added potential is (EA/(8L)) * (that
    // integral without L/2... ) — assembled so that dH/dQ reproduces
    // the KC modal force  (EA/(4L)) * S * (k pi/L)^2 * L Q_k with
    // S = sum (j pi/L)^2 Q_j^2 / ... The mass-normalized channel:
    let e: Vec<f64> = {
        let mut m = vec![0.0; n_modes * n_modes];
        for k in 1..=n_modes {
            let kk = k as f64 * pi / length;
            m[(k - 1) * n_modes + (k - 1)] = kk * kk * 2.0 / (lin_density * length);
        }
        m
    };
    let coefficient = ea * length / 8.0;
    Ok(SosModalStorage {
        omegas,
        channels: vec![StressChannel {
            coefficient,
            coupling: e,
        }],
    })
}

/// Kirchhoff–Carrier on a free-fixed (moving-end) basis.
///
/// Eigenfunctions `φ_k(x) = √(2/μL) cos(κ_k x)` with
/// `κ_k = (k−½)π/L`. `φ_k(0) ≠ 0`, so the attachment port can
/// Dirac-join a body. The averaged-tension channel is the same
/// SOS form as [`kirchhoff_carrier_string`]:
/// `E[k,k] = κ_k² · 2/(μ L)`, coefficient `E A L / 8`. Fixed-fixed
/// sines have `φ(0)=0` and cannot carry this join.
///
/// # Errors
/// [`NlModalError::Parameter`] on non-physical inputs.
pub fn kirchhoff_carrier_moving_end(
    params: &KcStringParams,
    n_modes: usize,
) -> Result<SosModalStorage, NlModalError> {
    let KcStringParams {
        length,
        tension,
        lin_density,
        ea,
    } = *params;
    if !(length > 0.0 && tension > 0.0 && lin_density > 0.0 && ea >= 0.0) || n_modes == 0 {
        return Err(NlModalError::Parameter {
            what: "string parameters",
        });
    }
    let c = det::sqrt(tension / lin_density);
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = (0..n_modes)
        .map(|k| (k as f64 + 0.5) * pi * c / length)
        .collect();
    let e: Vec<f64> = {
        let mut m = vec![0.0; n_modes * n_modes];
        for k in 0..n_modes {
            let kappa = (k as f64 + 0.5) * pi / length;
            m[k * n_modes + k] = kappa * kappa * 2.0 / (lin_density * length);
        }
        m
    };
    let coefficient = ea * length / 8.0;
    Ok(SosModalStorage {
        omegas,
        channels: vec![StressChannel {
            coefficient,
            coupling: e,
        }],
    })
}

/// Linear angular frequency [rad/s] of pinned sine mode `k` (1-based).
///
/// Uses the same Euler–Bernoulli dispersion as
/// [`prestressed_beam_omega_at_wavenumber`], with `κ = k π/L`.
/// Zero tension retains the pure bending spectrum; `ei == 0` gives the
/// flexible taut-string limit. Negative tension retains the linear Euler
/// compression boundary: zero stiffness gives zero frequency, and an unstable
/// mode gives NaN. Invalid geometry or data also returns NaN.
#[must_use]
pub fn prestressed_beam_omega(
    length: f64,
    tension: f64,
    lin_density: f64,
    ei: f64,
    k: usize,
) -> f64 {
    if !length.is_finite() || length <= 0.0 || k == 0 {
        return f64::NAN;
    }
    prestressed_beam_omega_at_wavenumber(
        k as f64 * core::f64::consts::PI / length,
        tension,
        lin_density,
        ei,
    )
}

/// Euler–Bernoulli angular frequency [rad/s] at spatial wavenumber `κ` [rad/m].
///
/// `ω² = (T κ² + EI κ⁴)/μ`, directly from
/// `μ u_tt = T u_xx - EI u_xxxx`. The caller owns boundary conditions and
/// selects its spatial spectrum; this function does not infer support types.
/// It is shared by pinned strings and the moving-end Galerkin reduction.
///
/// Combining the tension and bending frequencies with a scaled norm avoids
/// dividing by tension or squaring large frequencies. Thus the zero-tension
/// beam and small-stiffness string limits remain continuous. Compressed modes
/// retain the signed stiffness sum and the linear Euler threshold. Invalid
/// geometry, density or bending stiffness and nonfinite inputs return NaN;
/// unrepresentable arithmetic is not clamped. Both restoring terms may be zero.
///
/// This slender, small-slope model omits rotary inertia and transverse shear:
/// its formal high-wavenumber speeds are not a broadband material-speed claim.
/// The linear compression limit does not model postbuckling or nonlinear stability.
/// See Ducceschi and Bilbao, JASA 140 (2016), doi:10.1121/1.4962553.
#[must_use]
pub fn prestressed_beam_omega_at_wavenumber(
    wave_number: f64,
    tension: f64,
    lin_density: f64,
    ei: f64,
) -> f64 {
    if !wave_number.is_finite()
        || wave_number < 0.0
        || !tension.is_finite()
        || !lin_density.is_finite()
        || lin_density <= 0.0
        || !ei.is_finite()
        || ei < 0.0
    {
        return f64::NAN;
    }
    if wave_number == 0.0 {
        return 0.0;
    }
    if tension < 0.0 {
        let modal_load = if ei == 0.0 {
            tension
        } else {
            tension + (ei * wave_number) * wave_number
        };
        return wave_number * (det::sqrt(modal_load) / det::sqrt(lin_density));
    }
    // Separate square roots avoid overflowing T/μ or EI/μ while their
    // square roots are still representable.
    let root_mass = det::sqrt(lin_density);
    let omega_tension = wave_number * (det::sqrt(tension) / root_mass);
    if ei == 0.0 {
        return omega_tension;
    }
    let omega_bending = wave_number * (wave_number * (det::sqrt(ei) / root_mass));
    let scale = omega_tension.max(omega_bending);
    if scale == 0.0 || !scale.is_finite() {
        return scale;
    }
    let ratio = omega_tension.min(omega_bending) / scale;
    scale * det::sqrt(1.0 + ratio * ratio)
}

#[cfg(test)]
mod beam_dispersion_tests {
    use super::{prestressed_beam_omega, prestressed_beam_omega_at_wavenumber};

    fn close(actual: f64, expected: f64) {
        assert!(actual.is_finite(), "nonfinite dispersion: {actual}");
        assert!(
            (actual - expected).abs() <= 4.0e-14 * expected.abs(),
            "dispersion {actual} != {expected}"
        );
    }

    // G1: independently substitute sine eigenfunctions in μ u_tt = T u_xx - EI u_xxxx.
    #[test]
    fn g1_modes_follow_the_tension_plus_bending_operator() {
        let (length, tension, mass, ei): (f64, f64, f64, f64) = (0.65, 80.0, 0.003, 0.006);
        for mode in [1, 2, 7, 32, 1024] {
            let kappa = mode as f64 * core::f64::consts::PI / length;
            let expected = ((tension * kappa * kappa + ei * kappa.powi(4)) / mass).sqrt();
            close(
                prestressed_beam_omega(length, tension, mass, ei, mode),
                expected,
            );
        }
    }

    #[test]
    fn g1_string_and_zero_tension_beam_limits_are_continuous() {
        let (length, tension, mass, ei): (f64, f64, f64, f64) = (0.8, 60.0, 0.005, 0.02);
        for mode in [1, 2, 10, 80] {
            let kappa = mode as f64 * core::f64::consts::PI / length;
            close(
                prestressed_beam_omega(length, tension, mass, 0.0, mode),
                kappa * (tension / mass).sqrt(),
            );
            let bending = kappa * kappa * (ei / mass).sqrt();
            for small_tension in [0.0, f64::from_bits(1), f64::MIN_POSITIVE, 1.0e-200] {
                close(
                    prestressed_beam_omega(length, small_tension, mass, ei, mode),
                    bending,
                );
            }
            for small_ei in [f64::from_bits(1), f64::MIN_POSITIVE, 1.0e-200] {
                close(
                    prestressed_beam_omega(length, tension, mass, small_ei, mode),
                    kappa * (tension / mass).sqrt(),
                );
            }
        }
    }

    #[test]
    fn g1_phase_and_group_speeds_have_the_physical_dispersion() {
        let kappa = 7.0;
        let h = 1.0e-4;
        let mass = 0.02;
        for (tension, ei) in [(80.0, 0.0), (0.0, 0.01), (80.0, 0.01)] {
            let omega = |k| prestressed_beam_omega_at_wavenumber(k, tension, mass, ei);
            let phase_speed = omega(kappa) / kappa;
            let group_speed = (omega(kappa + h) - omega(kappa - h)) / (2.0 * h);
            let expected_ratio =
                (tension + 2.0 * ei * kappa * kappa) / (tension + ei * kappa * kappa);
            assert!((group_speed / phase_speed - expected_ratio).abs() < 2.0e-9);
        }
    }

    #[test]
    fn g0_large_finite_frequencies_do_not_square_or_divide_by_tension() {
        // The squared frequency and both T/μ and EI/μ overflow; ω remains finite.
        close(
            prestressed_beam_omega_at_wavenumber(1.0e50, 1.0e200, 1.0e-200, 1.0e200),
            1.0e300,
        );
        assert_eq!(
            prestressed_beam_omega_at_wavenumber(3.0, 0.0, 1.0, 0.0),
            0.0
        );
        assert_eq!(
            prestressed_beam_omega_at_wavenumber(0.0, 1.0, 1.0, 1.0),
            0.0
        );
    }

    #[test]
    fn g0_invalid_stiffness_never_becomes_a_flexible_string() {
        for ei in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(prestressed_beam_omega(1.0, 1.0, 1.0, ei, 1).is_nan());
        }
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(prestressed_beam_omega(invalid, 1.0, 1.0, 1.0, 1).is_nan());
            assert!(prestressed_beam_omega(1.0, 1.0, invalid, 1.0, 1).is_nan());
            assert!(prestressed_beam_omega_at_wavenumber(invalid, 1.0, 1.0, 1.0).is_nan());
        }
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(prestressed_beam_omega(1.0, invalid, 1.0, 1.0, 1).is_nan());
        }
        assert!(prestressed_beam_omega(0.0, 1.0, 1.0, 1.0, 1).is_nan());
        assert!(prestressed_beam_omega(1.0, 1.0, 0.0, 1.0, 1).is_nan());
        assert!(prestressed_beam_omega(1.0, 1.0, 1.0, 1.0, 0).is_nan());
    }
}

/// First-order Duffing backbone: for
/// `qddot + w0^2 q + beta q^3 = 0`, the amplitude-dependent frequency
/// is `w(a) = w0 (1 + 3 beta a^2 / (8 w0^2))` — the analytic
/// perturbation pin the batteries measure against. This is valid only for
/// `|beta| a^2 / w0^2 << 1`; use [`duffing_backbone_exact`] for the full
/// conservative scalar orbit at a physical turning-point amplitude.
#[must_use]
pub fn duffing_backbone(omega0: f64, beta: f64, amplitude: f64) -> f64 {
    omega0 * (1.0 + 3.0 * beta * amplitude * amplitude / (8.0 * omega0 * omega0))
}

/// The effective single-mode cubic coefficient `beta_k` of mode `k`
/// (the `q_k^3` force coefficient with every other mode at rest):
/// `beta_k = sum_j c_j E_j[k,k]^2`.
#[must_use]
pub fn single_mode_beta(storage: &SosModalStorage, k: usize) -> f64 {
    let n = storage.omegas.len();
    storage
        .channels
        .iter()
        .map(|ch| {
            let e = ch.coupling[k * n + k];
            ch.coefficient * e * e
        })
        .sum()
}

#[cfg(test)]
mod fd_probe {
    use super::{
        SineMode, SosModalStorage, StressChannel, analytic_coupling_integral, coupling_integral,
        fd_coupling, fd_d2, gauss_01, sample_sine, triple_sine_integral,
    };
    use fs_phs::Storage;

    #[test]
    fn g1_separable_bracket_matches_independent_area_quadrature() {
        // Evaluate the physical Hessian bracket point by point, without the
        // production separation or parity shortcut and with native trig.
        let direct = |lx: f64,
                      ly: f64,
                      s: SineMode,
                      a: SineMode,
                      b: SineMode,
                      nodes: &[f64],
                      weights: &[f64]| {
            let wave = |m: SineMode| {
                (
                    m.m as f64 * core::f64::consts::PI / lx,
                    m.n as f64 * core::f64::consts::PI / ly,
                )
            };
            let (sx, sy) = wave(s);
            let (ax, ay) = wave(a);
            let (bx, by) = wave(b);
            let mut result = 0.0;
            for (&xi, &wx) in nodes.iter().zip(weights) {
                for (&eta, &wy) in nodes.iter().zip(weights) {
                    let (x, y) = (lx * xi, ly * eta);
                    let psi = (sx * x).sin() * (sy * y).sin();
                    let av = (ax * x).sin() * (ay * y).sin();
                    let bv = (bx * x).sin() * (by * y).sin();
                    let axy = ax * ay * (ax * x).cos() * (ay * y).cos();
                    let bxy = bx * by * (bx * x).cos() * (by * y).cos();
                    result += wx
                        * wy
                        * psi
                        * ((-ax * ax * av) * (-by * by * bv) + (-ay * ay * av) * (-bx * bx * bv)
                            - 2.0 * axy * bxy);
                }
            }
            result * lx * ly
        };
        let mode = |m, n| SineMode { m, n };
        let triples = [
            (mode(1, 1), mode(1, 1), mode(1, 1)),
            (mode(3, 5), mode(6, 2), mode(4, 6)),
            (mode(4, 4), mode(4, 6), mode(1, 1)),
            (mode(3, 5), mode(13, 9), mode(1, 1)),
            (mode(4, 6), mode(1, 3), mode(6, 8)),
            (mode(2, 1), mode(1, 1), mode(1, 1)),
        ];
        for order in [96, 109] {
            let (nodes, weights) = gauss_01(order);
            for (lx, ly) in [(0.7, 0.4), (2.3, 0.9)] {
                for (s, a, b) in triples {
                    let actual = coupling_integral(lx, ly, s, a, b, &nodes, &weights);
                    let expected = direct(lx, ly, s, a, b, &nodes, &weights);
                    assert!(
                        (actual - expected).abs() < 2.0e-10 * expected.abs().max(1.0),
                        "order {order}, {s:?}, {a:?}, {b:?}: {actual} vs {expected}"
                    );
                    let swapped = coupling_integral(lx, ly, s, b, a, &nodes, &weights);
                    assert!((actual - swapped).abs() < 2.0e-12 * actual.abs().max(1.0));
                    let analytic = analytic_coupling_integral(lx, ly, s, a, b);
                    assert!(
                        (analytic - expected).abs() < 2.0e-10 * expected.abs().max(1.0),
                        "analytic bracket: {s:?}, {a:?}, {b:?}: {analytic} vs {expected}"
                    );
                    assert_eq!(analytic, analytic_coupling_integral(lx, ly, s, b, a));
                }
            }
        }
    }

    #[test]
    fn exact_spatial_harmonics_retain_far_overlap_and_near_resonance() {
        // For a=b=1, integrate sin(s*pi*x)*(1-cos(2*pi*x))/2.
        // Combining the antiderivatives before floating point gives this
        // independent reduced identity, including its sign for s>2.
        for s in [1_usize, 3, 257, 1_000_001] {
            let sf = s as f64;
            let expected = -4.0 / (core::f64::consts::PI * sf * (sf * sf - 4.0));
            let actual = triple_sine_integral(s, 1, 1);
            assert!((actual / expected - 1.0).abs() < 4.0 * f64::EPSILON);
            assert_eq!(actual, triple_sine_integral(1, s, 1));
            assert_eq!(actual, triple_sine_integral(1, 1, s));
        }
        // A sum-frequency neighbour has a denominator factor of exactly one.
        // Mode differences must be formed before conversion to floating point.
        let n = 1_000_000_usize;
        let actual = triple_sine_integral(2 * n + 1, n, n);
        let nf = n as f64;
        let expected =
            -4.0 * nf * nf / (core::f64::consts::PI * (4.0 * nf + 1.0) * (2.0 * nf + 1.0));
        assert!((actual / expected - 1.0).abs() < 4.0 * f64::EPSILON);
        assert_eq!(triple_sine_integral(2 * n, n, n), 0.0);
    }

    #[test]
    fn quadratic_metadata_preserves_the_live_stress_remainder() {
        let mut storage = SosModalStorage {
            omegas: vec![2.0, 3.0],
            channels: vec![StressChannel {
                coefficient: 7.0,
                coupling: vec![1.0, 0.0, 0.0, 2.0],
            }],
        };
        let part = storage.quadratic_part().unwrap();
        assert!(!part.exact);
        assert_eq!(
            part.matrix,
            vec![
                4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0
            ]
        );
        let x = [1.0, 2.0, 3.0, 4.0];
        let linear = 0.5 * (4.0 + 4.0 + 81.0 + 16.0);
        assert_eq!(
            storage.hamiltonian(&x),
            linear + 0.25 * 7.0 * 19.0_f64.powi(2)
        );
        storage.channels[0].coefficient = 0.0;
        assert!(storage.quadratic_part().unwrap().exact);
        assert_eq!(storage.hamiltonian(&x), linear);
    }

    #[test]
    fn sine_second_derivative_orders_agree_at_center() {
        let (nx, ny) = (21, 17);
        let (lx, ly) = (0.4, 0.3);
        let dx = lx / (nx - 1) as f64;
        let dy = ly / (ny - 1) as f64;
        let w = sample_sine(SineMode { m: 1, n: 1 }, lx, ly, nx, ny);
        let i = nx / 2;
        let j = ny / 2;
        let a2 = fd_d2(&w, nx, i, j, 0, dx, 2);
        let a4 = fd_d2(&w, nx, i, j, 0, dx, 4);
        let want = -(core::f64::consts::PI / lx).powi(2)
            * (core::f64::consts::PI * i as f64 / (nx - 1) as f64).sin()
            * (core::f64::consts::PI * j as f64 / (ny - 1) as f64).sin();
        assert!((a4 - want).abs() / want.abs() < 0.02);
        assert!((a2 - a4).abs() / a4.abs() < 0.05);
        let psi = sample_sine(SineMode { m: 1, n: 1 }, lx, ly, nx, ny);
        let c2 = fd_coupling(&psi, &w, &w, nx, ny, dx, dy, 2);
        let c4 = fd_coupling(&psi, &w, &w, nx, ny, dx, dy, 4);
        assert!(
            c4.abs() > 1.0e-6 && (c2 - c4).abs() / c4.abs() < 0.05,
            "interior coupling 2={c2} 4={c4}"
        );
    }
}

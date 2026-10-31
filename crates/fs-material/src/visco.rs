//! Viscoelastic damping constitutive models (bead
//! frankensim-fsim-visco-damping-ybc75, musical-acoustics program) — the
//! decisive instrument physics: damping sets decay times, brightness-vs-
//! time, and the audible difference between spruce and plywood.
//!
//! FOUR TIERS, one runtime path:
//! - Tier 0 [`RayleighDamping`]: classical `C = αM + βK` and constant modal
//!   damping — trivial, but it must exist so quick studies are not forced
//!   into Prony.
//! - Tier a [`FractionalZener`]: the FITTING-side canonical model. Wood,
//!   polymers, and felt show nearly constant loss factor over decades of
//!   frequency; integer-order Prony needs many terms to mimic that, the
//!   4-parameter fractional Zener captures it directly. Evaluated
//!   generically over [`fs_ad::Real`], so exact AD parameter tangents come
//!   from the same code path (and the f64 instantiation routes through
//!   `fs_math::det` — no platform libm in solver state).
//! - Tier b [`GeneralizedMaxwell`] (Prony series): the RUNTIME model.
//!   Unconditionally stable, convolution-free exact-exponential recursive
//!   update for piecewise-linear strain, with a per-step dissipation
//!   ledger (work = stored + dissipated, closure tested at 1e-12).
//! - [`lower_to_prony`]: the bridge — deterministic log-spaced relaxation
//!   ladder + Lawson–Hanson NNLS fit of the fractional model, returning a
//!   [`LoweredModel`] carrying its VALIDITY BAND and the measured supremum
//!   relative modulus error over a dense in-band verification grid.
//!   Evaluation outside the band REFUSES (`FS-MAT-VISCO-OUT-OF-BAND`) —
//!   the certificate is a measured bound on the stated band, never a
//!   global claim.
//! - [`ThermoelasticZener`]: Zener's closed-form thermoelastic loss for
//!   transverse beam/plate vibration — the dominant damping in metal bars
//!   (vibraphone) at low modes, nearly free analytically.
//!
//! Modal consumption: `ζ_k = η(ω_k)/2` ([`loss_factor_to_zeta`]); the
//! resulting complex eigenvalues match the exact quadratic roots for light
//! damping (tested).

use crate::MaterialError;
use fs_ad::Real;
use fs_ad::dual::gradient;
use fs_math::det;

/// Passive full-tensor Prony (generalized Maxwell) viscoelastic law in Mandel coordinates.
pub mod tensor;
/// Zener thermoelastic damping of solid circular sections (exact thermal-mode series).
pub mod circular_ted;
/// Bounded two-exponent fractional law, including measured violin-varnish spectra.
pub mod modified_zener;
pub mod source;
pub mod stress_control;
pub use modified_zener::ModifiedFractionalZener;
pub use tensor::{Mandel6, PronyBranch, TensorProny, TensorPronyState};

/// Typed refusals of the viscoelastic layer. Display strings carry stable
/// `FS-MAT-VISCO-*` codes.
#[derive(Debug, Clone, PartialEq)]
pub enum ViscoError {
    /// Parameters outside the admissible region.
    Parameters {
        /// Diagnosis.
        what: &'static str,
    },
    /// A lowered model was evaluated outside its certified band.
    OutOfBand {
        /// Requested angular frequency [rad/s].
        omega: f64,
        /// The certified band [rad/s].
        band: (f64, f64),
    },
    /// The NNLS fit could not reach the authored tolerance.
    FitTolerance {
        /// Measured supremum relative error over the verification grid.
        achieved: f64,
        /// The authored gate.
        required: f64,
    },
    /// The nonlinear fractional-Zener fit exhausted its iteration budget
    /// without meeting the convergence criteria.
    FitDiverged {
        /// Final relative residual norm.
        residual: f64,
        /// Iterations consumed.
        iterations: usize,
    },
}

impl core::fmt::Display for ViscoError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ViscoError::Parameters { what } => {
                write!(f, "FS-MAT-VISCO-PARAMETERS: {what}")
            }
            ViscoError::OutOfBand { omega, band } => write!(
                f,
                "FS-MAT-VISCO-OUT-OF-BAND: omega {omega:.3e} outside certified band \
                 [{:.3e}, {:.3e}] rad/s",
                band.0, band.1
            ),
            ViscoError::FitTolerance { achieved, required } => write!(
                f,
                "FS-MAT-VISCO-FIT-TOLERANCE: sup relative modulus error {achieved:.3e} \
                 exceeds the authored gate {required:.3e}"
            ),
            ViscoError::FitDiverged {
                residual,
                iterations,
            } => write!(
                f,
                "FS-MAT-VISCO-FIT-DIVERGED: relative residual {residual:.3e} after \
                 {iterations} iterations without meeting the convergence criteria"
            ),
        }
    }
}

impl std::error::Error for ViscoError {}

impl From<ViscoError> for MaterialError {
    fn from(e: ViscoError) -> MaterialError {
        MaterialError::Parameters {
            what: e.to_string(),
        }
    }
}

/// Modal damping ratio from a loss factor: `ζ = η/2` (valid for light
/// damping, the structural-dynamics convention).
#[must_use]
pub fn loss_factor_to_zeta(eta: f64) -> f64 {
    0.5 * eta
}

// ---------------------------------------------------------------------------
// Tier 0: Rayleigh
// ---------------------------------------------------------------------------

/// Classical Rayleigh damping `C = αM + βK`; `ζ(ω) = α/(2ω) + βω/2`.
/// The matrix itself is assembled by the consumer holding M and K (fs-time's
/// generalized-α takes C directly); this type owns the coefficients and the
/// per-mode damping ratio.
#[derive(Debug, Clone, Copy)]
pub struct RayleighDamping {
    /// Mass-proportional coefficient α [1/s].
    pub alpha: f64,
    /// Stiffness-proportional coefficient β [s].
    pub beta: f64,
}

impl RayleighDamping {
    /// Admission: both coefficients finite and non-negative.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`].
    pub fn new(alpha: f64, beta: f64) -> Result<RayleighDamping, ViscoError> {
        if !(alpha.is_finite() && beta.is_finite() && alpha >= 0.0 && beta >= 0.0) {
            return Err(ViscoError::Parameters {
                what: "Rayleigh coefficients must be finite and non-negative",
            });
        }
        Ok(RayleighDamping { alpha, beta })
    }

    /// Damping ratio at angular frequency ω.
    #[must_use]
    pub fn zeta_at(&self, omega: f64) -> f64 {
        0.5 * (self.alpha / omega + self.beta * omega)
    }

    /// Coefficients matching prescribed ratios at two frequencies (the
    /// classical two-point fit). A passive fit must match BOTH pins: an
    /// incompatible pair refuses instead of clipping a negative coefficient
    /// and silently changing the requested resonance damping.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`] on invalid pins, an incompatible passive
    /// fit, or unrepresentable coefficients.
    pub fn from_two_points(
        omega1: f64,
        zeta1: f64,
        omega2: f64,
        zeta2: f64,
    ) -> Result<RayleighDamping, ViscoError> {
        if !(omega1.is_finite()
            && omega2.is_finite()
            && zeta1.is_finite()
            && zeta2.is_finite()
            && omega1 > 0.0
            && omega2 > 0.0
            && omega1.to_bits() != omega2.to_bits()
            && zeta1 >= 0.0
            && zeta2 >= 0.0)
        {
            return Err(ViscoError::Parameters {
                what: "two-point Rayleigh fit needs finite distinct positive frequencies and nonnegative damping",
            });
        }
        let (lo, zlo, hi, zhi) = if omega1 < omega2 {
            (omega1, zeta1, omega2, zeta2)
        } else {
            (omega2, zeta2, omega1, zeta1)
        };
        // Solve for alpha/lo and beta*hi using r=lo/hi. No squared
        // dimensional frequency appears, and (hi-lo)/hi retains the gap
        // even when r rounds close to one.
        let r = lo / hi;
        let denominator = ((hi - lo) / hi) * (1.0 + r);
        let mass_part = zlo - r * zhi;
        let stiffness_part = zhi - r * zlo;
        // Only roundoff at a pure-mass/pure-stiffness boundary is allowed.
        // A materially negative coefficient cannot represent passive damping
        // at both supplied pins and must not be replaced with zero.
        let mass_roundoff = 8.0 * f64::EPSILON * zlo.max(r * zhi);
        let stiffness_roundoff = 8.0 * f64::EPSILON * zhi.max(r * zlo);
        if mass_part < -mass_roundoff || stiffness_part < -stiffness_roundoff {
            return Err(ViscoError::Parameters {
                what: "requested damping ratios have no nonnegative two-point Rayleigh fit",
            });
        }
        let alpha = (2.0 * (mass_part.max(0.0) / denominator)) * lo;
        let beta = (2.0 * (stiffness_part.max(0.0) / denominator)) / hi;
        let model = RayleighDamping::new(alpha, beta)?;
        // Check the original equations, including representability of a very
        // small coefficient. A clipped roundoff term cannot erase a real pin.
        for (omega, requested) in [(lo, zlo), (hi, zhi)] {
            let actual = model.zeta_at(omega);
            let allowance = 64.0 * f64::EPSILON * requested;
            if !actual.is_finite() || (actual - requested).abs() > allowance {
                return Err(ViscoError::Parameters {
                    what: "two-point Rayleigh coefficients do not reproduce the requested damping",
                });
            }
        }
        Ok(model)
    }
}

// ---------------------------------------------------------------------------
// Tier a: fractional Zener (fitting canonical)
// ---------------------------------------------------------------------------

/// Four-parameter fractional Zener solid:
/// `E*(ω) = (E0 + E∞·(iωτ)^α) / (1 + (iωτ)^α)`, `0 < α ≤ 1`,
/// `0 < E0 < E∞`, `τ > 0`. At α = 1 it reduces exactly to the classical
/// Zener (standard linear solid); small α gives the nearly-frequency-
/// constant loss factor that wood and felt exhibit.
#[derive(Debug, Clone, Copy)]
pub struct FractionalZener {
    /// Equilibrium (low-frequency) modulus E0 [Pa].
    pub e0: f64,
    /// Instantaneous (high-frequency) modulus E∞ [Pa].
    pub e_inf: f64,
    /// Fractional order α ∈ (0, 1].
    pub alpha: f64,
    /// Relaxation time τ [s].
    pub tau: f64,
}

impl FractionalZener {
    /// Admission of the physical parameter region.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`].
    pub fn new(e0: f64, e_inf: f64, alpha: f64, tau: f64) -> Result<FractionalZener, ViscoError> {
        let finite = e0.is_finite() && e_inf.is_finite() && alpha.is_finite() && tau.is_finite();
        if !finite || e0 <= 0.0 || e_inf <= e0 || alpha <= 0.0 || alpha > 1.0 || tau <= 0.0 {
            return Err(ViscoError::Parameters {
                what: "fractional Zener needs 0 < E0 < Einf, 0 < alpha <= 1, tau > 0",
            });
        }
        Ok(FractionalZener {
            e0,
            e_inf,
            alpha,
            tau,
        })
    }

    /// Storage and loss modulus `(E′, E″)` at angular frequency ω — the
    /// [`Real`]-generic core shared by the f64 path and the AD tangents.
    pub fn modulus_parts<R: Real>(e0: R, e_inf: R, alpha: R, tau: R, omega: R) -> (R, R) {
        // (iωτ)^α = (ωτ)^α · (cos(απ/2) + i·sin(απ/2)).
        let half_pi = R::from_f64(core::f64::consts::FRAC_PI_2);
        let x = (alpha * (omega * tau).ln()).exp();
        let c = (alpha * half_pi).cos();
        let s = (alpha * half_pi).sin();
        let xr = x * c;
        let xi = x * s;
        // (E0 + E∞·(xr + i·xi)) / (1 + xr + i·xi).
        let nr = e0 + e_inf * xr;
        let ni = e_inf * xi;
        let dr = R::one() + xr;
        let di = xi;
        let den = dr * dr + di * di;
        let ep = (nr * dr + ni * di) / den;
        let epp = (ni * dr - nr * di) / den;
        (ep, epp)
    }

    /// `(E′, E″)` at ω [rad/s] (deterministic `fs_math::det` elementary
    /// functions via the f64 [`Real`] instantiation).
    #[must_use]
    pub fn modulus(&self, omega: f64) -> (f64, f64) {
        FractionalZener::modulus_parts(self.e0, self.e_inf, self.alpha, self.tau, omega)
    }

    /// Loss factor `η(ω) = E″/E′`.
    #[must_use]
    pub fn loss_factor(&self, omega: f64) -> f64 {
        let (ep, epp) = self.modulus(omega);
        epp / ep
    }

    /// Exact AD gradients of `(E′, E″)` with respect to
    /// `[E0, E∞, α, τ]` at fixed ω (forward-mode duals through the same
    /// generic code path as the value).
    #[must_use]
    pub fn modulus_gradients(&self, omega: f64) -> ((f64, [f64; 4]), (f64, [f64; 4])) {
        let x = [self.e0, self.e_inf, self.alpha, self.tau];
        let storage = gradient(x, |v| {
            let w = fs_ad::dual::Dual64::<4>::constant(omega);
            FractionalZener::modulus_parts(v[0], v[1], v[2], v[3], w).0
        });
        let loss = gradient(x, |v| {
            let w = fs_ad::dual::Dual64::<4>::constant(omega);
            FractionalZener::modulus_parts(v[0], v[1], v[2], v[3], w).1
        });
        (storage, loss)
    }
}

// ---------------------------------------------------------------------------
// Tier b: generalized Maxwell (Prony) — the runtime model
// ---------------------------------------------------------------------------

/// Generalized Maxwell (Prony) series:
/// `E*(ω) = E∞ + Σ_j E_j·(iωτ_j)/(1 + iωτ_j)`.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralizedMaxwell {
    /// Equilibrium modulus [Pa].
    pub e_inf: f64,
    /// `(E_j, τ_j)` branch pairs [Pa, s].
    pub terms: Vec<(f64, f64)>,
}

impl GeneralizedMaxwell {
    /// Admission: positive equilibrium modulus, non-negative branch moduli,
    /// positive relaxation times.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`].
    pub fn new(e_inf: f64, terms: Vec<(f64, f64)>) -> Result<GeneralizedMaxwell, ViscoError> {
        if !(e_inf > 0.0 && e_inf.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "equilibrium modulus must be positive",
            });
        }
        for &(e, tau) in &terms {
            if !(e >= 0.0 && e.is_finite() && tau > 0.0 && tau.is_finite()) {
                return Err(ViscoError::Parameters {
                    what: "Prony terms need E_j >= 0 and tau_j > 0",
                });
            }
        }
        Ok(GeneralizedMaxwell { e_inf, terms })
    }

    /// One Maxwell branch that matches loss factor `η` at `omega`
    /// with relaxation time `τ = 1/ω`.
    ///
    /// At `ωτ = 1`, `η = E₁ / (2 E∞ + E₁)`, so
    /// `E₁ = 2 η E∞ / (1 − η)`. This is the integer-order runtime
    /// stand-in for a single authored `η(ω₀)`, not a wood spectrum.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`] if `e_inf` is not positive, `omega`
    /// is not positive, or `η` is not in `(0, 1)`.
    pub fn matching_loss(
        e_inf: f64,
        omega: f64,
        eta: f64,
    ) -> Result<GeneralizedMaxwell, ViscoError> {
        if !(e_inf > 0.0 && e_inf.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "equilibrium modulus must be positive",
            });
        }
        if !(omega > 0.0 && omega.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "match frequency must be positive",
            });
        }
        if !(eta > 0.0 && eta < 1.0 && eta.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "loss factor must lie in (0, 1)",
            });
        }
        let e1 = 2.0 * eta * e_inf / (1.0 - eta);
        GeneralizedMaxwell::new(e_inf, vec![(e1, 1.0 / omega)])
    }

    /// `(E′, E″)` at ω [rad/s].
    #[must_use]
    pub fn modulus(&self, omega: f64) -> (f64, f64) {
        let mut ep = self.e_inf;
        let mut epp = 0.0;
        for &(e, tau) in &self.terms {
            let wt = omega * tau;
            let den = 1.0 + wt * wt;
            ep += e * wt * wt / den;
            epp += e * wt / den;
        }
        (ep, epp)
    }

    /// Loss factor `η(ω) = E″/E′`.
    #[must_use]
    pub fn loss_factor(&self, omega: f64) -> f64 {
        let (ep, epp) = self.modulus(omega);
        epp / ep
    }

    /// Modal damping ratios `ζ_k = η(ω_k)/2` from this Prony series.
    ///
    /// Runtime bridge into a modal bank: one authored
    /// [`Self::matching_loss`] (or a lowered fractional model) becomes
    /// a frequency-dependent `ζ` vector. It is not Rayleigh `αM+βK`.
    #[must_use]
    pub fn modal_zetas(&self, omegas: &[f64]) -> Vec<f64> {
        omegas
            .iter()
            .map(|&w| loss_factor_to_zeta(self.loss_factor(w)))
            .collect()
    }

    /// Fresh time-stepping state (all internal variables relaxed, ledger
    /// zeroed).
    #[must_use]
    pub fn state(&self) -> PronyState {
        PronyState {
            q: vec![0.0; self.terms.len()],
            strain: 0.0,
            stress: 0.0,
            work: 0.0,
            dissipated: 0.0,
        }
    }

    /// Advance the SAME exact piecewise-linear-strain constitutive history.
    /// Work integrates the actual (exponential, not linearly interpolated)
    /// stress, and heat integrates `sum(E_j q_j(t)^2 / tau_j)`. Endpoint stress
    /// and the original branch-strain recursion are retained. In particular,
    /// a reversal cannot manufacture negative heat through trapezoidal work.
    /// `dt = 0` is an elastic strain jump with exactly zero viscous loss.
    ///
    /// Energies are densities [J/m³], not joules. This fixed-temperature law
    /// does not itself advance a thermal field or change source parameters.
    pub fn step(&self, state: &mut PronyState, strain_new: f64, dt: f64) -> f64 {
        assert!(
            dt.is_finite() && dt >= 0.0,
            "Prony step requires a finite non-negative dt"
        );
        assert_eq!(state.q.len(), self.terms.len(), "Prony state must match its branches");
        let d_eps = strain_new - state.strain;
        let mut mean_stress = self.e_inf * (state.strain + 0.5 * d_eps);
        let mut stress = self.e_inf * strain_new;
        let mut dissipated = 0.0;
        for (q, &(e, tau)) in state.q.iter_mut().zip(&self.terms) {
            let (next, mean, loss) = prony_ramp(*q, d_eps, dt / tau);
            *q = next;
            mean_stress += e * mean;
            stress += e * next;
            dissipated += e * loss;
        }
        state.strain = strain_new;
        state.stress = stress;
        state.work += mean_stress * d_eps;
        state.dissipated += dissipated;
        stress
    }
}

/// Internal-variable state and energy ledger for [`GeneralizedMaxwell`]
/// time stepping.
#[derive(Debug, Clone)]
pub struct PronyState {
    /// Internal (viscous-spring) strains, one per branch.
    pub q: Vec<f64>,
    /// Current total strain.
    pub strain: f64,
    /// Current stress.
    pub stress: f64,
    /// Accumulated exact ramp work density [J/m³], relative to initialization.
    pub work: f64,
    /// Accumulated physical dashpot heat density [J/m³]. Integrated independently
    /// of work and stored energy; nonnegative and nondecreasing for admitted data.
    pub dissipated: f64,
}

// Exact scalar Maxwell ramp; r = dt/tau. Returns endpoint strain, mean
// strain, and integral(q²/tau dt). Complete the square in the analytic
// integral rather than subtracting nearly equal work/storage quantities:
// A*(q0+h*delta)² + R*delta², A=(1-exp(-2r))/2,
// h=tanh(r/2)/r, R=(r-2*tanh(r/2))/r² >= 0.
pub(crate) fn prony_ramp(q0: f64, delta: f64, r: f64) -> (f64, f64, f64) {
    if r == 0.0 {
        return (q0 + delta, q0 + 0.5 * delta, 0.0);
    }
    let a = det::exp(-r);
    let s = -det::expm1(-r);
    let b = s / r;
    let h = b / (1.0 + a);
    let (c, remainder) = if r < 0.1 {
        // Analytic series, not a loss clamp. Avoid cancellation in both
        // (1-b)/r and (1-2h)/r at acoustic dt << material relaxation time.
        let c = 0.5 + r * (-1.0 / 6.0 + r * (1.0 / 24.0
            + r * (-1.0 / 120.0 + r * (1.0 / 720.0
            + r * (-1.0 / 5040.0 + r * (1.0 / 40320.0
            + r * (-1.0 / 362880.0 + r / 3628800.0)))))));
        let r2 = r * r;
        let remainder = r * (1.0 / 12.0 + r2 * (-1.0 / 120.0
            + r2 * (17.0 / 20160.0 + r2 * (-31.0 / 362880.0
            + r2 * (691.0 / 79833600.0)))));
        (c, remainder)
    } else {
        // Ordered divisions also cover an overflowing dt/tau: the fully
        // relaxed branch releases its old stored energy, with no NaN.
        ((1.0 - b) / r, (1.0 - 2.0 * h) / r)
    };
    let shifted = q0 + h * delta;
    let loss = (0.5 * s * (1.0 + a)) * shifted * shifted
        + (remainder * delta) * delta;
    (a * q0 + b * delta, b * q0 + c * delta, loss)
}

// ---------------------------------------------------------------------------
// The bridge: fractional → Prony lowering with a measured certificate
// ---------------------------------------------------------------------------

/// A Prony model lowered from a fractional Zener, carrying its certified
/// validity band and the MEASURED supremum relative complex-modulus error
/// over the in-band verification grid. The certificate is exactly that
/// measurement — a dense-grid bound on the stated band, not an analytic
/// global proof (no-claim boundary).
#[derive(Debug, Clone)]
pub struct LoweredModel {
    /// The runtime Prony model.
    pub model: GeneralizedMaxwell,
    /// Certified band [rad/s].
    pub band: (f64, f64),
    /// Measured `sup |E*_prony − E*_frac| / |E*_frac|` over the grid.
    pub sup_rel_err: f64,
    /// Measured supremum relative storage-modulus error on the same grid.
    pub sup_rel_storage_err: f64,
    /// Measured supremum relative loss-modulus error on the same grid.
    /// This is normalized by loss, so weak damping cannot hide behind storage.
    pub sup_rel_loss_err: f64,
    /// Number of verification grid points.
    pub verification_points: usize,
}

impl LoweredModel {
    /// `(E′, E″)` at ω, REFUSING outside the certified band.
    ///
    /// # Errors
    /// [`ViscoError::OutOfBand`].
    pub fn modulus_checked(&self, omega: f64) -> Result<(f64, f64), ViscoError> {
        if !(omega >= self.band.0 && omega <= self.band.1) {
            return Err(ViscoError::OutOfBand {
                omega,
                band: self.band,
            });
        }
        Ok(self.model.modulus(omega))
    }

    /// Loss factor at ω with the band refusal.
    ///
    /// # Errors
    /// [`ViscoError::OutOfBand`].
    pub fn loss_factor_checked(&self, omega: f64) -> Result<f64, ViscoError> {
        let (ep, epp) = self.modulus_checked(omega)?;
        Ok(epp / ep)
    }

    /// Read-only diagnostic: supremum relative storage-modulus error
    /// `max_ω |E′_prony(ω) − E′_frac(ω)| / E′_frac(ω)` over the verification grid.
    #[must_use]
    pub fn sup_rel_storage_error(&self, fz: &FractionalZener) -> f64 {
        let mut worst = 0.0f64;
        let log_lo = det::ln(self.band.0);
        let log_hi = det::ln(self.band.1);
        let n = self.verification_points.max(2);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let w = det::exp(log_lo + t * (log_hi - log_lo));
            let (fe, _) = fz.modulus(w);
            let (ge, _) = self.model.modulus(w);
            worst = worst.max((ge - fe).abs() / fe.abs().max(1e-300));
        }
        worst
    }

    /// Read-only diagnostic: supremum relative loss-modulus error
    /// `max_ω |E″_prony(ω) − E″_frac(ω)| / E″_frac(ω)` over the verification grid.
    #[must_use]
    pub fn sup_rel_loss_error(&self, fz: &FractionalZener) -> f64 {
        let mut worst = 0.0f64;
        let log_lo = det::ln(self.band.0);
        let log_hi = det::ln(self.band.1);
        let n = self.verification_points.max(2);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let w = det::exp(log_lo + t * (log_hi - log_lo));
            let (_, fi) = fz.modulus(w);
            let (_, gi) = self.model.modulus(w);
            worst = worst.max((gi - fi).abs() / fi.abs().max(1e-300));
        }
        worst
    }
}

/// Dense-grid supremum relative complex-modulus error of `gm` against `fz`
/// over `[w_lo, w_hi]` (log-spaced `points`).
#[must_use]
pub fn sup_rel_modulus_error(
    fz: &FractionalZener,
    gm: &GeneralizedMaxwell,
    w_lo: f64,
    w_hi: f64,
    points: usize,
) -> f64 {
    let mut worst = 0.0f64;
    let log_lo = det::ln(w_lo);
    let log_hi = det::ln(w_hi);
    for i in 0..points {
        let t = i as f64 / (points - 1) as f64;
        let w = det::exp(log_lo + t * (log_hi - log_lo));
        let (fe, fi) = fz.modulus(w);
        let (ge, gi) = gm.modulus(w);
        let num = det::sqrt((ge - fe) * (ge - fe) + (gi - fi) * (gi - fi));
        let den = det::sqrt(fe * fe + fi * fi);
        worst = worst.max(num / den);
    }
    worst
}

/// Lower a fractional Zener to an `n_terms` Prony series certified on the
/// band `[f_lo, f_hi]` Hz: deterministic log-spaced relaxation-time ladder
/// spanning the band with margin, Lawson–Hanson NNLS on storage+loss
/// samples (E∞ pinned to the fractional model's E0 so the ω→0 limits agree
/// exactly), then a 512-point verification sweep. Refuses if the measured
/// supremum error exceeds `tol`.
///
/// # Errors
/// [`ViscoError::Parameters`] on a bad band or term count;
/// [`ViscoError::FitTolerance`] when the certificate misses `tol`.
#[allow(clippy::too_many_lines)] // ladder + NNLS + certification, one pipeline
pub fn lower_to_prony(
    fz: &FractionalZener,
    f_lo: f64,
    f_hi: f64,
    n_terms: usize,
    tol: f64,
) -> Result<LoweredModel, ViscoError> {
    lower_scalar_to_prony(
        |omega| Ok(fz.modulus(omega)),
        fz.e0,
        (f_lo, f_hi),
        n_terms,
        10.0,
        ScalarLoweringTolerance {
            complex: Some(tol),
            storage: None,
            loss: None,
        },
    )
}

/// Component gates are separate from the original complex-error API. A small
/// complex residual alone says little about accuracy of a weak loss modulus.
struct ScalarLoweringTolerance {
    complex: Option<f64>,
    storage: Option<f64>,
    loss: Option<f64>,
}

/// One positive-Prony solver for both fractional parameterizations. The target
/// supplies its own applicability checks; no sampled value is extrapolated.
#[allow(clippy::too_many_lines)]
fn lower_scalar_to_prony(
    target: impl Fn(f64) -> Result<(f64, f64), ViscoError>,
    equilibrium: f64,
    (f_lo, f_hi): (f64, f64),
    n_terms: usize,
    margin: f64,
    tolerance: ScalarLoweringTolerance,
) -> Result<LoweredModel, ViscoError> {
    if !(f_lo > 0.0 && f_hi > f_lo && f_lo.is_finite() && f_hi.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "lowering band needs 0 < f_lo < f_hi",
        });
    }
    if n_terms == 0 || n_terms > 64 {
        return Err(ViscoError::Parameters {
            what: "term count must be in 1..=64",
        });
    }
    if !(equilibrium > 0.0 && equilibrium.is_finite() && margin >= 1.0 && margin.is_finite())
        || [tolerance.complex, tolerance.storage, tolerance.loss]
            .into_iter()
            .flatten()
            .any(|value| !(value > 0.0 && value.is_finite()))
    {
        return Err(ViscoError::Parameters {
            what: "lowering requires finite positive equilibrium and tolerances, margin >= 1",
        });
    }
    let two_pi = 2.0 * core::f64::consts::PI;
    let w_lo = two_pi * f_lo;
    let w_hi = two_pi * f_hi;
    let tau_hi = margin / w_lo;
    let tau_lo = (1.0 / w_hi) / margin;
    if !(w_lo > 0.0 && w_hi.is_finite() && tau_lo > 0.0 && tau_hi.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "lowering band and relaxation ladder must remain finite and positive",
        });
    }
    let taus: Vec<f64> = (0..n_terms)
        .map(|j| {
            let t = if n_terms == 1 {
                0.5
            } else {
                j as f64 / (n_terms - 1) as f64
            };
            det::exp(det::ln(tau_lo) + t * (det::ln(tau_hi) - det::ln(tau_lo)))
        })
        .collect();
    // Sample matrix: rows = [E′(ω_i) − E0 ; E″(ω_i)], columns = branches.
    let m_samples = (8 * n_terms).max(48);
    let omegas: Vec<f64> = (0..m_samples)
        .map(|i| lowering_grid_frequency(w_lo, w_hi, i, m_samples))
        .collect();
    let rows = 2 * m_samples;
    let mut a = vec![0.0f64; rows * n_terms];
    let mut rhs = vec![0.0f64; rows];
    for (i, &w) in omegas.iter().enumerate() {
        let (fe, fi) = target(w)?;
        if !(fe > 0.0 && fe.is_finite() && fi >= 0.0 && fi.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "lowering target must have finite positive storage and nonnegative loss",
            });
        }
        // Weight storage and loss separately. Scaling both by |E*| allows a
        // nearly elastic target to lose all damping while appearing accurate.
        // The floor keeps zero/sub-roundoff loss rows numerically bounded;
        // the final loss gate still measures loss itself, not this floor.
        let storage_scale = fe;
        let loss_scale = fi.max(fe * 1.0e-12).max(1.0e-300);
        rhs[i] = (fe - equilibrium) / storage_scale;
        rhs[m_samples + i] = fi / loss_scale;
        for (j, &tau) in taus.iter().enumerate() {
            let wt = w * tau;
            let (storage, loss) = if wt <= 1.0 {
                let den = 1.0 + wt * wt;
                (wt * wt / den, wt / den)
            } else {
                let inv = 1.0 / wt;
                let den = 1.0 + inv * inv;
                (1.0 / den, inv / den)
            };
            a[i * n_terms + j] = storage / storage_scale;
            a[(m_samples + i) * n_terms + j] = loss / loss_scale;
        }
    }
    // NNLS uses an absolute gradient threshold. Unit-norm columns make its
    // decisions independent of whether moduli are Pa, MPa, or another scale.
    let mut norms = vec![0.0; n_terms];
    for j in 0..n_terms {
        let scale = (0..rows)
            .map(|i| a[i * n_terms + j].abs())
            .fold(0.0_f64, f64::max);
        let norm = scale
            * det::sqrt(
                (0..rows)
                    .map(|i| {
                        let value = a[i * n_terms + j] / scale;
                        value * value
                    })
                    .sum(),
            );
        if !(norm > 0.0 && norm.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "lowering matrix has a nonrepresentable column norm",
            });
        }
        norms[j] = norm;
        for i in 0..rows {
            a[i * n_terms + j] /= norm;
        }
    }
    let mut coeffs = nnls(&a, &rhs, rows, n_terms);
    for (coefficient, norm) in coeffs.iter_mut().zip(norms) {
        *coefficient /= norm;
        if !coefficient.is_finite() || *coefficient < 0.0 {
            return Err(ViscoError::Parameters {
                what: "lowering returned a nonfinite or negative Prony modulus",
            });
        }
    }
    let terms: Vec<(f64, f64)> = coeffs
        .iter()
        .zip(&taus)
        .filter(|pair| *pair.0 > 0.0)
        .map(|(&e, &tau)| (e, tau))
        .collect();
    let model = GeneralizedMaxwell::new(equilibrium, terms)?;
    let verification_points = 512;
    let mut errors = [0.0_f64; 3];
    for i in 0..verification_points {
        let w = lowering_grid_frequency(w_lo, w_hi, i, verification_points);
        let (fe, fi) = target(w)?;
        let (ge, gi) = model.modulus(w);
        let storage_error = (ge - fe) / fe;
        let loss_error = (gi - fi) / fi.abs().max(1.0e-300);
        let imaginary_relative = (gi - fi) / fe;
        let loss_ratio = fi / fe;
        let complex_error = det::sqrt(
            (storage_error * storage_error + imaginary_relative * imaginary_relative)
                / (1.0 + loss_ratio * loss_ratio),
        );
        for (worst, error) in
            errors
                .iter_mut()
                .zip([complex_error, storage_error.abs(), loss_error.abs()])
        {
            if !error.is_finite() {
                return Err(ViscoError::FitTolerance {
                    achieved: f64::INFINITY,
                    required: tolerance
                        .loss
                        .or(tolerance.storage)
                        .or(tolerance.complex)
                        .unwrap_or(0.0),
                });
            }
            *worst = worst.max(error);
        }
    }
    for (achieved, required) in
        errors
            .into_iter()
            .zip([tolerance.complex, tolerance.storage, tolerance.loss])
    {
        if let Some(required) = required {
            if achieved > required {
                return Err(ViscoError::FitTolerance { achieved, required });
            }
        }
    }
    Ok(LoweredModel {
        model,
        band: (w_lo, w_hi),
        sup_rel_err: errors[0],
        sup_rel_storage_err: errors[1],
        sup_rel_loss_err: errors[2],
        verification_points,
    })
}

fn lowering_grid_frequency(lo: f64, hi: f64, index: usize, points: usize) -> f64 {
    if index == 0 {
        lo
    } else if index + 1 == points {
        hi
    } else {
        let fraction = index as f64 / (points - 1) as f64;
        det::exp(det::ln(lo) + fraction * (det::ln(hi) - det::ln(lo))).clamp(lo, hi)
    }
}

/// Lawson–Hanson non-negative least squares (dense, deterministic): solve
/// `min ‖A x − b‖` s.t. `x ≥ 0`, `A` row-major `rows x cols`. Sizes here are
/// tiny (≤ 64 unknowns). A column enters while its gradient exceeds 1e-12,
/// an absolute threshold, so callers scale `A`'s columns to unit norm.
#[must_use]
pub fn nnls(a: &[f64], b: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut passive = vec![false; cols];
    let mut x = vec![0.0f64; cols];
    for _outer in 0..(2 * cols + 4) {
        // Gradient of ½‖Ax−b‖²: w = Aᵀ(b − Ax).
        let mut r = b.to_vec();
        for (i, ri) in r.iter_mut().enumerate() {
            for j in 0..cols {
                *ri -= a[i * cols + j] * x[j];
            }
        }
        let mut best = usize::MAX;
        let mut best_w = 1e-12;
        for j in 0..cols {
            if !passive[j] {
                let mut wj = 0.0;
                for i in 0..rows {
                    wj += a[i * cols + j] * r[i];
                }
                if wj > best_w {
                    best_w = wj;
                    best = j;
                }
            }
        }
        if best == usize::MAX {
            break; // KKT satisfied
        }
        passive[best] = true;
        // Inner loop: solve on the passive set; demote negatives.
        loop {
            let idx: Vec<usize> = (0..cols).filter(|&j| passive[j]).collect();
            let z = ls_solve_subset(a, b, rows, cols, &idx);
            if z.iter().all(|&v| v > 0.0) {
                for (k, &j) in idx.iter().enumerate() {
                    x[j] = z[k];
                }
                for (j, slot) in x.iter_mut().enumerate() {
                    if !passive[j] {
                        *slot = 0.0;
                    }
                }
                break;
            }
            // Step toward z until the first variable hits zero; demote it.
            let mut alpha = 1.0f64;
            let mut drop_j = usize::MAX;
            for (k, &j) in idx.iter().enumerate() {
                if z[k] <= 0.0 {
                    let step = x[j] / (x[j] - z[k]);
                    if step < alpha {
                        alpha = step;
                        drop_j = j;
                    }
                }
            }
            for (k, &j) in idx.iter().enumerate() {
                x[j] += alpha * (z[k] - x[j]);
            }
            if drop_j == usize::MAX {
                break;
            }
            passive[drop_j] = false;
            x[drop_j] = 0.0;
        }
    }
    x
}

/// Least squares on a column subset via normal equations + Gaussian
/// elimination with partial pivoting (sizes ≤ 64).
fn ls_solve_subset(a: &[f64], b: &[f64], rows: usize, cols: usize, idx: &[usize]) -> Vec<f64> {
    let k = idx.len();
    let mut ata = vec![0.0f64; k * k];
    let mut atb = vec![0.0f64; k];
    for (p, &jp) in idx.iter().enumerate() {
        for (q, &jq) in idx.iter().enumerate() {
            let mut acc = 0.0;
            for i in 0..rows {
                acc += a[i * cols + jp] * a[i * cols + jq];
            }
            ata[p * k + q] = acc;
        }
        let mut acc = 0.0;
        for i in 0..rows {
            acc += a[i * cols + jp] * b[i];
        }
        atb[p] = acc;
    }
    // Gaussian elimination with partial pivoting.
    for col in 0..k {
        let mut piv = col;
        for r in col + 1..k {
            if ata[r * k + col].abs() > ata[piv * k + col].abs() {
                piv = r;
            }
        }
        if piv != col {
            for c in 0..k {
                ata.swap(col * k + c, piv * k + c);
            }
            atb.swap(col, piv);
        }
        let d = ata[col * k + col];
        if d.abs() < 1e-300 {
            continue;
        }
        for r in col + 1..k {
            let f = ata[r * k + col] / d;
            for c in col..k {
                ata[r * k + c] -= f * ata[col * k + c];
            }
            atb[r] -= f * atb[col];
        }
    }
    let mut z = vec![0.0f64; k];
    for r in (0..k).rev() {
        let mut acc = atb[r];
        for c in r + 1..k {
            acc -= ata[r * k + c] * z[c];
        }
        let d = ata[r * k + r];
        z[r] = if d.abs() < 1e-300 { 0.0 } else { acc / d };
    }
    z
}

// ---------------------------------------------------------------------------
// Fitting front-end: complex-modulus data → fractional Zener
// ---------------------------------------------------------------------------

/// Result of [`fit_fractional_zener`]: the admitted model plus the
/// per-iteration relative residual norms (the fit's replayable evidence).
#[derive(Debug, Clone)]
pub struct ZenerFit {
    /// The fitted, admission-checked model.
    pub model: FractionalZener,
    /// Per-iteration relative residual norms.
    pub residual_history: Vec<f64>,
}

/// Fit a [`FractionalZener`] to measured complex-modulus samples
/// `(ω [rad/s], E′ [Pa], E″ [Pa])` by damped Gauss–Newton
/// (Levenberg–Marquardt) over an unconstrained reparameterization that
/// enforces the admissible region by construction:
/// `E0 = exp(θ0)`, `E∞ = E0 + exp(θ1)`, `α = logistic(θ2)`, `τ = exp(θ3)`.
/// Residuals are relative to each sample's complex-modulus magnitude, and
/// the Jacobian comes from the exact AD tangents
/// ([`FractionalZener::modulus_gradients`]) chain-ruled through the
/// transform — no finite differences anywhere.
///
/// Deterministic: fixed iteration policy, fixed damping schedule, no
/// randomness.
///
/// `tol_rel` is the caller-AUTHORED residual budget: measured data cannot
/// fit below its own noise floor, so the caller states what residual level
/// counts as a successful fit (budgets-first, like every other gate in
/// this workspace). Convergence is stationarity — an accepted step with
/// `‖δθ‖ < 1e-13`, a damping factor no improving step can beat, or an
/// early residual below `1e-12` — and a stationary point is accepted only
/// when its relative residual norm is within `tol_rel`; otherwise the fit
/// REFUSES rather than returning a model that cannot represent the data
/// (e.g. active/negative-loss samples fed to a passive model).
///
/// # Errors
/// [`ViscoError::Parameters`] for fewer than 4 samples, non-finite or
/// non-positive data, or a non-finite/non-positive `tol_rel`;
/// [`ViscoError::FitDiverged`] when the fit goes stationary or exhausts
/// `max_iters` above the authored budget.
pub fn fit_fractional_zener(
    samples: &[(f64, f64, f64)],
    init: &FractionalZener,
    max_iters: usize,
    tol_rel: f64,
) -> Result<ZenerFit, ViscoError> {
    if !(tol_rel.is_finite() && tol_rel > 0.0) {
        return Err(ViscoError::Parameters {
            what: "the residual budget tol_rel must be finite and positive",
        });
    }
    if samples.len() < 4 {
        return Err(ViscoError::Parameters {
            what: "the 4-parameter fit needs at least 4 complex-modulus samples",
        });
    }
    for &(w, ep, epp) in samples {
        if !(w.is_finite() && ep.is_finite() && epp.is_finite() && w > 0.0 && ep > 0.0) {
            return Err(ViscoError::Parameters {
                what: "samples need finite omega > 0, E' > 0, finite E''",
            });
        }
    }

    // θ from the admitted initial model (its invariants guarantee the
    // transforms below are well-defined).
    let logit = |p: f64| det::ln(p / (1.0 - p));
    let mut theta = [
        det::ln(init.e0),
        det::ln(init.e_inf - init.e0),
        // α = 1 sits on the closed admission boundary but at logit's pole;
        // nudge inside so the transform stays finite.
        logit(init.alpha.min(1.0 - 1e-9)),
        det::ln(init.tau),
    ];
    // Relative residual 2-norm and (residuals, Jacobian) at θ.
    let m = samples.len() * 2;
    let mut jac = Vec::new();
    let (mut residuals, mut norm) = zener_fit_residuals(samples, &theta, Some(&mut jac), m);
    let mut history = vec![norm];
    let mut lambda = 1e-3f64;
    for _iteration in 0..max_iters {
        if norm < 1e-12 {
            let (e0, e_inf, alpha, tau) = unpack_zener_fit_parameters(&theta);
            let model = FractionalZener::new(e0, e_inf, alpha, tau)?;
            return Ok(ZenerFit {
                model,
                residual_history: history,
            });
        }
        // Normal equations with Marquardt scaling on the diagonal.
        let mut ata = [0.0f64; 16];
        let mut atb = [0.0f64; 4];
        for (row, &ri) in residuals.iter().enumerate() {
            for c in 0..4 {
                atb[c] -= jac[row * 4 + c] * ri;
                for c2 in 0..4 {
                    ata[c * 4 + c2] += jac[row * 4 + c] * jac[row * 4 + c2];
                }
            }
        }
        let mut damped = ata;
        for d in 0..4 {
            damped[d * 4 + d] += lambda * ata[d * 4 + d].max(1e-30);
        }
        let delta = solve4(&damped, &atb);
        let candidate = [
            theta[0] + delta[0],
            theta[1] + delta[1],
            theta[2] + delta[2],
            theta[3] + delta[3],
        ];
        let (cand_res, cand_norm) = zener_fit_residuals(samples, &candidate, None, m);
        if cand_norm < norm {
            theta = candidate;
            lambda = (lambda * 0.5).max(1e-12);
            let step = det::sqrt(
                delta[0] * delta[0]
                    + delta[1] * delta[1]
                    + delta[2] * delta[2]
                    + delta[3] * delta[3],
            );
            let (_, refreshed_norm) = zener_fit_residuals(samples, &theta, Some(&mut jac), m);
            residuals = cand_res;
            norm = refreshed_norm;
            history.push(norm);
            if step < 1e-13 {
                break;
            }
        } else {
            lambda *= 4.0;
            if lambda > 1e12 {
                // No improving step exists at any damping: numerically
                // stationary. Acceptance is decided by the budget below.
                break;
            }
        }
    }
    if norm <= tol_rel {
        let (e0, e_inf, alpha, tau) = unpack_zener_fit_parameters(&theta);
        let model = FractionalZener::new(e0, e_inf, alpha, tau)?;
        return Ok(ZenerFit {
            model,
            residual_history: history,
        });
    }
    Err(ViscoError::FitDiverged {
        residual: norm,
        iterations: history.len() - 1,
    })
}

fn unpack_zener_fit_parameters(theta: &[f64; 4]) -> (f64, f64, f64, f64) {
    let e0 = det::exp(theta[0]);
    let e_inf = e0 + det::exp(theta[1]);
    let alpha = 1.0 / (1.0 + det::exp(-theta[2]));
    let tau = det::exp(theta[3]);
    (e0, e_inf, alpha, tau)
}

fn zener_fit_residuals(
    samples: &[(f64, f64, f64)],
    theta: &[f64; 4],
    jac: Option<&mut Vec<f64>>,
    m: usize,
) -> (Vec<f64>, f64) {
    let (e0, e_inf, alpha, tau) = unpack_zener_fit_parameters(theta);
    // Chain rule of the transform: rows are raw params, cols are θ.
    let dalpha = alpha * (1.0 - alpha);
    let chain = [
        [e0, 0.0, 0.0, 0.0],
        [e0, e_inf - e0, 0.0, 0.0],
        [0.0, 0.0, dalpha, 0.0],
        [0.0, 0.0, 0.0, tau],
    ];
    let mut r = Vec::with_capacity(m);
    let mut jrows = jac;
    if let Some(j) = jrows.as_deref_mut() {
        j.clear();
        j.reserve(m * 4);
    }
    let mut sq = 0.0f64;
    for &(w, ep_d, epp_d) in samples {
        let mag = det::sqrt(ep_d * ep_d + epp_d * epp_d);
        let (ep, epp) = FractionalZener::modulus_parts(e0, e_inf, alpha, tau, w);
        let r0 = (ep - ep_d) / mag;
        let r1 = (epp - epp_d) / mag;
        sq += r0 * r0 + r1 * r1;
        r.push(r0);
        r.push(r1);
        if let Some(j) = jrows.as_deref_mut() {
            let probe = FractionalZener {
                e0,
                e_inf,
                alpha,
                tau,
            };
            let ((_, gp), (_, gl)) = probe.modulus_gradients(w);
            for grad in [gp, gl] {
                for c in 0..4 {
                    let mut acc = 0.0;
                    for (raw, chain_row) in chain.iter().enumerate() {
                        acc += grad[raw] * chain_row[c];
                    }
                    j.push(acc / mag);
                }
            }
        }
    }
    (r, det::sqrt(sq / m as f64))
}

/// Dense 4×4 solve by Gaussian elimination with partial pivoting
/// (deterministic; the LM normal equations are tiny and well-damped).
fn solve4(a: &[f64; 16], b: &[f64; 4]) -> [f64; 4] {
    let mut m = *a;
    let mut rhs = *b;
    let mut perm = [0usize, 1, 2, 3];
    for col in 0..4 {
        let mut pivot = col;
        for row in col + 1..4 {
            if m[perm[row] * 4 + col].abs() > m[perm[pivot] * 4 + col].abs() {
                pivot = row;
            }
        }
        perm.swap(col, pivot);
        let d = m[perm[col] * 4 + col];
        if d.abs() < 1e-300 {
            continue;
        }
        for row in col + 1..4 {
            let f = m[perm[row] * 4 + col] / d;
            for c in col..4 {
                m[perm[row] * 4 + c] -= f * m[perm[col] * 4 + c];
            }
            rhs[perm[row]] -= f * rhs[perm[col]];
        }
    }
    let mut x = [0.0f64; 4];
    for row in (0..4).rev() {
        let mut acc = rhs[perm[row]];
        for c in row + 1..4 {
            acc -= m[perm[row] * 4 + c] * x[c];
        }
        let d = m[perm[row] * 4 + row];
        x[row] = if d.abs() < 1e-300 { 0.0 } else { acc / d };
    }
    x
}

// ---------------------------------------------------------------------------
// Thermoelastic damping (Zener's closed form)
// ---------------------------------------------------------------------------

/// Zener's thermoelastic damping for transverse vibration of a thin
/// beam/plate: `η(ω) = Δ·ωτ/(1 + (ωτ)²)` with relaxation strength
/// `Δ = Eα²T₀/(ρc_p)` and thermal time `τ = h²ρc_p/(π²κ)`. Dominant in
/// vibraphone-bar metals at low modes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoelasticZener {
    /// Young's modulus E [Pa].
    pub e: f64,
    /// Linear thermal expansion α_T [1/K].
    pub alpha_t: f64,
    /// Absolute temperature T₀ [K].
    pub t0: f64,
    /// Density ρ [kg/m³].
    pub rho: f64,
    /// Specific heat c_p [J/(kg·K)].
    pub cp: f64,
    /// Thermal conductivity κ [W/(m·K)].
    pub conductivity: f64,
}

impl ThermoelasticZener {
    /// Relaxation strength `Δ = Eα²T₀/(ρc_p)`.
    #[must_use]
    pub fn relaxation_strength(&self) -> f64 {
        self.e * self.alpha_t * self.alpha_t * self.t0 / (self.rho * self.cp)
    }

    /// Thermal relaxation time for thickness h: `τ = h²ρc_p/(π²κ)`.
    #[must_use]
    pub fn relaxation_time(&self, thickness: f64) -> f64 {
        let pi = core::f64::consts::PI;
        thickness * thickness * self.rho * self.cp / (pi * pi * self.conductivity)
    }

    /// Loss factor at angular frequency ω for thickness h.
    #[must_use]
    pub fn loss_factor(&self, omega: f64, thickness: f64) -> f64 {
        let tau = self.relaxation_time(thickness);
        let wt = omega * tau;
        self.relaxation_strength() * wt / (1.0 + wt * wt)
    }

    /// Handbook aluminum (6061-class) at temperature `t0` [K].
    ///
    /// # Errors
    /// [`ViscoError::Parameters`] if `t0` is not positive and finite.
    pub fn aluminum(t0: f64) -> Result<Self, ViscoError> {
        if !(t0 > 0.0 && t0.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "thermoelastic T0 must be positive and finite",
            });
        }
        Ok(Self {
            e: 69.0e9,
            alpha_t: 23.1e-6,
            t0,
            rho: 2_700.0,
            cp: 897.0,
            conductivity: 167.0,
        })
    }

    /// Handbook structural steel at temperature `t0` [K].
    ///
    /// # Errors
    /// [`ViscoError::Parameters`] if `t0` is not positive and finite.
    pub fn structural_steel(t0: f64) -> Result<Self, ViscoError> {
        if !(t0 > 0.0 && t0.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "thermoelastic T0 must be positive and finite",
            });
        }
        Ok(Self {
            e: 200.0e9,
            alpha_t: 11.7e-6,
            t0,
            rho: 7_850.0,
            cp: 486.0,
            conductivity: 51.9,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    #[test]
    fn single_maxwell_closed_forms_over_six_decades() {
        // One-branch Maxwell + equilibrium spring: E′, E″ closed forms and
        // the loss peak at ωτ = 1 with η_max = E1/(2√(E∞(E∞+E1))).
        let e_inf = 2.0e9;
        let e1 = 1.0e9;
        let tau = 1e-3;
        let gm = GeneralizedMaxwell::new(e_inf, vec![(e1, tau)]).expect("gm");
        for k in -3..=3 {
            let w = det::exp(det::ln(10.0) * f64::from(k)) / tau;
            let wt = w * tau;
            let (ep, epp) = gm.modulus(w);
            let ep_want = e_inf + e1 * wt * wt / (1.0 + wt * wt);
            let epp_want = e1 * wt / (1.0 + wt * wt);
            assert!((ep - ep_want).abs() < 1e-6 * ep_want);
            assert!((epp - epp_want).abs() < 1e-6 * epp_want.max(1.0));
        }
        // η peaks at ωτ = √(E∞/(E∞+E1)) — NOT at ωτ = 1, which is where
        // E″ peaks — with η_max = E1/(2√(E∞(E∞+E1))).
        let wt_star = det::sqrt(e_inf / (e_inf + e1));
        let eta_peak = gm.loss_factor(wt_star / tau);
        let eta_max_want = e1 / (2.0 * det::sqrt(e_inf * (e_inf + e1)));
        assert!(
            (eta_peak - eta_max_want).abs() < 1e-12 * eta_max_want,
            "peak loss {eta_peak} vs closed form {eta_max_want}"
        );
        // And it is a maximum: neighbors are lower.
        assert!(gm.loss_factor(1.5 * wt_star / tau) < eta_peak);
        assert!(gm.loss_factor(0.6 * wt_star / tau) < eta_peak);
    }

    #[test]
    fn fractional_zener_alpha_one_reduces_to_classical_zener_exactly() {
        let fz = FractionalZener::new(1.0e9, 3.0e9, 1.0, 2e-4).expect("fz");
        let gm = GeneralizedMaxwell::new(1.0e9, vec![(2.0e9, 2e-4)]).expect("gm");
        for k in -30..=30 {
            let w = det::exp(0.3 * f64::from(k)) / 2e-4;
            let (fe, fi) = fz.modulus(w);
            let (ge, gi) = gm.modulus(w);
            assert!((fe - ge).abs() < 1e-6 * ge.abs(), "E' at {w}");
            assert!((fi - gi).abs() < 1e-6 * gi.abs().max(1.0), "E'' at {w}");
        }
    }

    #[test]
    fn prony_recursion_matches_direct_convolution() {
        // Random piecewise-linear strain path: the recursive update must
        // match the O(n²) direct convolution of the same discretization at
        // 1e-12 relative.
        let gm = GeneralizedMaxwell::new(1.5e9, vec![(8.0e8, 3e-4), (4.0e8, 5e-3)]).expect("gm");
        let dt = 1e-4;
        let n = 400;
        let mut seed = 0x5EEDu64;
        let strains: Vec<f64> = (0..n).map(|_| lcg(&mut seed) * 1e-3).collect();
        let mut state = gm.state();
        let mut recursive = Vec::with_capacity(n);
        for &eps in &strains {
            recursive.push(gm.step(&mut state, eps, dt));
        }
        // Direct convolution: q_j(t_n) = Σ_m a^(n-m)·b·Δε_m exactly mirrors
        // the recursion algebra but is evaluated independently per step.
        let mut prev = 0.0;
        let deltas: Vec<f64> = strains
            .iter()
            .map(|&e| {
                let d = e - prev;
                prev = e;
                d
            })
            .collect();
        for (nstep, &eps) in strains.iter().enumerate() {
            let mut stress = gm.e_inf * eps;
            for &(e, tau) in &gm.terms {
                let a = det::exp(-dt / tau);
                let b = tau * (1.0 - a) / dt;
                let mut q = 0.0;
                for (m, &d) in deltas[..=nstep].iter().enumerate() {
                    // a^(n-m) via repeated multiplication is O(n²) total —
                    // an INDEPENDENT evaluation order from the recursion.
                    let mut pw = 1.0;
                    for _ in 0..(nstep - m) {
                        pw *= a;
                    }
                    q += pw * b * d;
                }
                stress += e * q;
            }
            let scale = stress.abs().max(1e3);
            assert!(
                (stress - recursive[nstep]).abs() / scale < 1e-12,
                "step {nstep}: convolution {stress} vs recursion {}",
                recursive[nstep]
            );
        }
    }

    #[test]
    fn dissipation_ledger_closes_and_mutation_is_caught() {
        let gm = GeneralizedMaxwell::new(1.0e9, vec![(5.0e8, 1e-3)]).expect("gm");
        let dt = 5e-5;
        let mut state = gm.state();
        let mut seed = 0xD15Cu64;
        let mut eps = 0.0;
        let mut last_d = 0.0;
        for step in 0..2000 {
            eps += lcg(&mut seed) * 2e-4;
            gm.step(&mut state, eps, dt);
            assert!(
                state.dissipated >= last_d - 1e-9 * state.work.abs().max(1.0),
                "dissipation must be monotone at step {step}"
            );
            last_d = state.dissipated;
        }
        assert!(state.dissipated > 0.0, "cyclic loading must dissipate");
        // Heat is now integrated independently of work and storage. The
        // check above still requires monotonicity. MUTATION: flip the
        // internal-variable update sign (the bead's named control) — the
        // ledger must go NEGATIVE, i.e. the mutant fabricates energy.
        let mut q = 0.0f64;
        let (e1, tau) = gm.terms[0];
        let mut eps_m = 0.0f64;
        let mut prev_eps = 0.0f64;
        let mut work = 0.0;
        let mut stress_old = 0.0;
        let mut seed2 = 0xD15Cu64;
        let mut min_d = f64::INFINITY;
        for _ in 0..2000 {
            eps_m += lcg(&mut seed2) * 2e-4;
            let a = det::exp(-dt / tau);
            let b = tau * (1.0 - a) / dt;
            q = a * q - b * (eps_m - prev_eps); // SIGN FLIPPED
            let stress = gm.e_inf * eps_m + e1 * q;
            work += 0.5 * (stress_old + stress) * (eps_m - prev_eps);
            stress_old = stress;
            prev_eps = eps_m;
            let stored = 0.5 * gm.e_inf * eps_m * eps_m + 0.5 * e1 * q * q;
            min_d = min_d.min(work - stored);
        }
        assert!(
            min_d < -1e-3,
            "sign-flipped internal update must fabricate energy (min ledger {min_d})"
        );
    }

    #[test]
    fn lowering_reproduces_constant_loss_plateau_with_eight_terms() {
        // α = 0.35 gives the nearly-frequency-constant loss factor wood
        // shows; 8 terms over 4 decades must certify at ≤ 2%.
        let fz = FractionalZener::new(8.0e9, 16.0e9, 0.35, 1e-3).expect("fz");
        let lowered = lower_to_prony(&fz, 20.0, 2.0e5, 8, 0.02).expect("lowering");
        assert!(lowered.model.terms.len() <= 8);
        assert!(
            lowered.sup_rel_err <= 0.02,
            "certificate: {}",
            lowered.sup_rel_err
        );
        // Plateau reproduced: lowered η stays within 10% of the fractional
        // η across the band interior.
        for k in 0..=20 {
            let w = det::exp(
                det::ln(2.0 * core::f64::consts::PI * 20.0)
                    + (f64::from(k) / 20.0)
                        * (det::ln(2.0 * core::f64::consts::PI * 2.0e5)
                            - det::ln(2.0 * core::f64::consts::PI * 20.0)),
            );
            // exp(ln(...)) can land one ulp outside the stored band edge.
            let w = w.clamp(lowered.band.0, lowered.band.1);
            let eta_f = fz.loss_factor(w);
            let eta_l = lowered.loss_factor_checked(w).expect("in band");
            assert!(
                (eta_l - eta_f).abs() / eta_f < 0.10,
                "plateau at w={w}: fractional {eta_f} vs lowered {eta_l}"
            );
        }
        // Out-of-band evaluation refuses with the named code.
        let err = lowered.loss_factor_checked(1.0).unwrap_err();
        assert!(err.to_string().contains("FS-MAT-VISCO-OUT-OF-BAND"));
        // MUTATION: dropping a Prony term must break the certificate.
        let mut truncated = lowered.model.clone();
        truncated.terms.pop();
        let err = sup_rel_modulus_error(&fz, &truncated, lowered.band.0, lowered.band.1, 256);
        assert!(
            err > lowered.sup_rel_err * 2.0,
            "dropped term must degrade the fit: {err} vs {}",
            lowered.sup_rel_err
        );
        println!(
            "{{\"suite\":\"fs-material-visco\",\"case\":\"fractional-lowering\",\"terms\":{},\"sup_rel_err\":{:.3e},\"verdict\":\"pass\"}}",
            lowered.model.terms.len(),
            lowered.sup_rel_err
        );
    }

    #[test]
    fn modal_injection_matches_quadratic_roots_for_light_damping() {
        // Single mode ω0, loss η: ζ = η/2; the quadratic
        // λ² + ηω0λ + ω0² = 0 has roots −ζω0 ± iω0√(1−ζ²).
        let omega0 = 2.0 * core::f64::consts::PI * 440.0;
        let eta = 0.02;
        let zeta = loss_factor_to_zeta(eta);
        let re_want = -zeta * omega0;
        let im_want = omega0 * det::sqrt(1.0 - zeta * zeta);
        // Solve the quadratic directly.
        let b = eta * omega0;
        let c = omega0 * omega0;
        let disc = b * b - 4.0 * c;
        assert!(disc < 0.0);
        let re = -b / 2.0;
        let im = det::sqrt(-disc) / 2.0;
        assert!((re - re_want).abs() < 1e-9 * omega0);
        assert!((im - im_want).abs() < 1e-9 * omega0);
    }

    #[test]
    fn maxwell_matching_loss_hits_the_authored_eta() {
        let gm = GeneralizedMaxwell::matching_loss(2.0e9, 1.0e3, 0.04).expect("match");
        let eta = gm.loss_factor(1.0e3);
        assert!((eta - 0.04).abs() < 1.0e-9, "eta={eta}");
        assert!(GeneralizedMaxwell::matching_loss(2.0e9, 1.0e3, 1.5).is_err());
        let z = gm.modal_zetas(&[1.0e3, 4.0e3]);
        assert!((z[0] - 0.02).abs() < 1.0e-9);
        assert!(
            (z[1] - z[0]).abs() > 1.0e-4,
            "Prony ζ must move with ω (got {} vs {})",
            z[1],
            z[0]
        );
    }

    #[test]
    fn handbook_metals_are_admitted() {
        let al = ThermoelasticZener::aluminum(293.0).expect("al");
        let st = ThermoelasticZener::structural_steel(293.0).expect("st");
        assert!(al.loss_factor(1.0e3, 0.002) > 0.0);
        assert!(st.loss_factor(1.0e3, 0.002) > 0.0);
        assert!(ThermoelasticZener::aluminum(0.0).is_err());
    }

    #[test]
    fn thermoelastic_zener_peak_and_aluminum_magnitude() {
        // Aluminum: E=70 GPa, α=23e-6/K, T=293K, ρ=2700, cp=900, κ=237.
        let ted = ThermoelasticZener {
            e: 70e9,
            alpha_t: 23e-6,
            t0: 293.0,
            rho: 2700.0,
            cp: 900.0,
            conductivity: 237.0,
        };
        let h = 0.01; // 10 mm vibraphone-bar scale
        let tau = ted.relaxation_time(h);
        let delta = ted.relaxation_strength();
        // Peak exactly at ωτ = 1 with η = Δ/2.
        let eta_peak = ted.loss_factor(1.0 / tau, h);
        assert!((eta_peak - delta / 2.0).abs() < 1e-15);
        // Aluminum relaxation strength ~4.5e-3 ⇒ peak loss ~2e-3.
        assert!(delta > 1e-3 && delta < 1e-2, "Delta = {delta}");
        // Symmetry in log frequency: η(ωτ = x) = η(ωτ = 1/x).
        for x in [0.1, 0.3, 3.0, 10.0] {
            let lo = ted.loss_factor(x / tau, h);
            let hi = ted.loss_factor(1.0 / (x * tau), h);
            assert!((lo - hi).abs() < 1e-15 * lo.max(hi));
        }
    }

    #[test]
    fn ad_gradients_match_finite_differences() {
        let fz = FractionalZener::new(5.0e9, 12.0e9, 0.4, 3e-4).expect("fz");
        let omega = 2.0 * core::f64::consts::PI * 1000.0;
        let ((ep, gep), (epp, gepp)) = fz.modulus_gradients(omega);
        let (ep0, epp0) = fz.modulus(omega);
        assert!((ep - ep0).abs() < 1e-6 * ep0.abs());
        assert!((epp - epp0).abs() < 1e-6 * epp0.abs());
        let x = [fz.e0, fz.e_inf, fz.alpha, fz.tau];
        for p in 0..4 {
            let h = 1e-6 * x[p].abs();
            let mut xp = x;
            xp[p] += h;
            let mut xm = x;
            xm[p] -= h;
            let fp = FractionalZener::modulus_parts(xp[0], xp[1], xp[2], xp[3], omega);
            let fm = FractionalZener::modulus_parts(xm[0], xm[1], xm[2], xm[3], omega);
            let fd_ep = (fp.0 - fm.0) / (2.0 * h);
            let fd_epp = (fp.1 - fm.1) / (2.0 * h);
            let scale_e = gep[p].abs().max(fd_ep.abs()).max(1e-12);
            let scale_l = gepp[p].abs().max(fd_epp.abs()).max(1e-12);
            assert!(
                (gep[p] - fd_ep).abs() / scale_e < 1e-5,
                "dE'/dp{p}: AD {} vs FD {fd_ep}",
                gep[p]
            );
            assert!(
                (gepp[p] - fd_epp).abs() / scale_l < 1e-5,
                "dE''/dp{p}: AD {} vs FD {fd_epp}",
                gepp[p]
            );
        }
    }

    #[test]
    fn rayleigh_tier_and_refusals() {
        let r = RayleighDamping::from_two_points(100.0, 0.01, 1000.0, 0.02).expect("fit");
        assert!((r.zeta_at(100.0) - 0.01).abs() < 1e-12);
        assert!((r.zeta_at(1000.0) - 0.02).abs() < 1e-12);
        assert!(RayleighDamping::new(-1.0, 0.0).is_err());
        assert!(FractionalZener::new(2.0e9, 1.0e9, 0.5, 1e-3).is_err()); // E∞ < E0
        assert!(FractionalZener::new(1.0e9, 2.0e9, 1.5, 1e-3).is_err()); // α > 1
        assert!(GeneralizedMaxwell::new(1.0e9, vec![(1.0, -1.0)]).is_err());
        let err = lower_to_prony(
            &FractionalZener::new(1.0e9, 2.0e9, 0.9, 1e-3).expect("fz"),
            100.0,
            10.0,
            4,
            0.1,
        )
        .unwrap_err();
        assert!(err.to_string().contains("FS-MAT-VISCO-PARAMETERS"));
    }

    #[test]
    fn rayleigh_fit_preserves_both_resonances_and_pin_order() {
        let reference = RayleighDamping::new(0.8, 2.0e-5).unwrap();
        let (w1, w2) = (80.0, 1800.0);
        let (z1, z2) = (reference.zeta_at(w1), reference.zeta_at(w2));
        let fit = RayleighDamping::from_two_points(w1, z1, w2, z2).unwrap();
        let reversed = RayleighDamping::from_two_points(w2, z2, w1, z1).unwrap();
        assert_eq!(fit.alpha.to_bits(), reversed.alpha.to_bits());
        assert_eq!(fit.beta.to_bits(), reversed.beta.to_bits());
        for omega in [50.0, w1, 300.0, w2, 10000.0] {
            assert!((fit.zeta_at(omega) / reference.zeta_at(omega) - 1.0).abs() < 2.0e-14);
        }
        // A mass-normalized oscillator's displacement at its natural
        // frequency is |F|/(2*zeta*omega^2). Check the actual physical
        // resonant compliance, not just signs of fitted coefficients.
        for (omega, zeta) in [(w1, z1), (w2, z2)] {
            let damping = fit.alpha + fit.beta * omega * omega;
            let compliance = 1.0 / (omega * damping);
            let expected = 1.0 / (2.0 * zeta * omega * omega);
            assert!((compliance / expected - 1.0).abs() < 2.0e-14);
        }
    }

    #[test]
    fn rayleigh_fit_refuses_incompatible_passive_loss() {
        // Both clipped fits used to succeed while missing one requested
        // resonant damping by more than a factor of ten.
        for (z1, z2) in [(0.1, 0.001), (0.001, 0.1), (0.0, 0.01), (0.01, 0.0)] {
            assert!(RayleighDamping::from_two_points(100.0, z1, 1000.0, z2).is_err());
            assert!(RayleighDamping::from_two_points(1000.0, z2, 100.0, z1).is_err());
        }
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
            assert!(RayleighDamping::from_two_points(invalid, 0.01, 1000.0, 0.02).is_err());
            assert!(RayleighDamping::from_two_points(100.0, invalid, 1000.0, 0.02).is_err());
            assert!(RayleighDamping::from_two_points(100.0, 0.01, invalid, 0.02).is_err());
            assert!(RayleighDamping::from_two_points(100.0, 0.01, 1000.0, invalid).is_err());
        }
    }

    #[test]
    fn rayleigh_fit_is_covariant_under_frequency_units() {
        let reference = RayleighDamping::from_two_points(100.0, 0.01, 1000.0, 0.02).unwrap();
        for scale in [1.0e-200, 1.0, 1.0e200] {
            let fit = RayleighDamping::from_two_points(
                100.0 * scale, 0.01, 1000.0 * scale, 0.02,
            ).unwrap();
            assert!((fit.alpha / scale / reference.alpha - 1.0).abs() < 2.0e-14);
            assert!((fit.beta * scale / reference.beta - 1.0).abs() < 2.0e-14);
            assert!((fit.zeta_at(100.0 * scale) - 0.01).abs() < 2.0e-16);
            assert!((fit.zeta_at(1000.0 * scale) - 0.02).abs() < 4.0e-16);
        }
    }

    #[test]
    fn rayleigh_fit_retains_pure_mass_stiffness_and_lossless_limits() {
        for reference in [
            RayleighDamping::new(8.0, 0.0).unwrap(),
            RayleighDamping::new(0.0, 0.125).unwrap(),
            RayleighDamping::new(0.0, 0.0).unwrap(),
        ] {
            let fit = RayleighDamping::from_two_points(
                2.0, reference.zeta_at(2.0), 16.0, reference.zeta_at(16.0),
            ).unwrap();
            assert_eq!(fit.alpha, reference.alpha);
            assert_eq!(fit.beta, reference.beta);
        }
    }

    #[test]
    fn lowering_is_bitwise_deterministic() {
        let fz = FractionalZener::new(8.0e9, 16.0e9, 0.35, 1e-3).expect("fz");
        let a = lower_to_prony(&fz, 20.0, 2.0e5, 8, 0.02).expect("a");
        let b = lower_to_prony(&fz, 20.0, 2.0e5, 8, 0.02).expect("b");
        assert_eq!(a.model.terms.len(), b.model.terms.len());
        for ((ea, ta), (eb, tb)) in a.model.terms.iter().zip(&b.model.terms) {
            assert_eq!(ea.to_bits(), eb.to_bits());
            assert_eq!(ta.to_bits(), tb.to_bits());
        }
        assert_eq!(a.sup_rel_err.to_bits(), b.sup_rel_err.to_bits());
    }

    #[test]
    fn prony_step_zero_dt_is_the_exact_elastic_jump_and_small_dt_stays_accurate() {
        let model = GeneralizedMaxwell::new(1.0e9, vec![(5.0e8, 10.0)]).expect("model");
        // dt = 0: instantaneous strain jump, branch strain follows
        // elastically (b = 1), nothing NaN.
        let mut state = model.state();
        let stress = model.step(&mut state, 1.0e-3, 0.0);
        assert!(stress.is_finite() && state.q[0].is_finite());
        assert!((state.q[0] - 1.0e-3).abs() < 1.0e-18);
        assert!((stress - (1.0e9 + 5.0e8) * 1.0e-3).abs() < 1.0e-3);
        // dt << tau: b must match the series limit 1 - h/2 + h^2/6 to
        // machine precision (the direct 1 - exp(-h) form loses ~six
        // digits at h = 2e-8 and returns NaN-free but wrong updates).
        let tau = 10.0;
        let dt = 2.0e-7; // h = 2e-8
        let h = dt / tau;
        let b_reference = 1.0 - h / 2.0 + h * h / 6.0;
        let b_impl = -tau * det::expm1(-dt / tau) / dt;
        assert!(((b_impl - b_reference) / b_reference).abs() < 1.0e-14);
    }

    #[test]
    fn fit_recovers_ground_truth_from_an_offset_start() {
        let truth = FractionalZener::new(9.0e9, 1.5e10, 0.35, 2.0e-4).expect("truth");
        let mut samples = Vec::new();
        for k in 0..24 {
            let w = det::exp(
                det::ln(1.0e1) + (f64::from(k) / 23.0) * (det::ln(1.0e7) - det::ln(1.0e1)),
            );
            let (ep, epp) = truth.modulus(w);
            samples.push((w, ep, epp));
        }
        // Every parameter starts far off (E0 /3, E∞ ×2.7, α ×1.7, τ /20).
        let init = FractionalZener::new(3.0e9, 4.0e10, 0.6, 1.0e-5).expect("init");
        let fit = fit_fractional_zener(&samples, &init, 200, 1.0e-10).expect("fit converges");
        assert!(fit.residual_history.len() >= 2, "history is evidence");
        let last = *fit.residual_history.last().expect("nonempty");
        assert!(
            last < 1.0e-10,
            "clean data fits to the numerical floor: {last}"
        );
        for (got, want) in [
            (fit.model.e0, truth.e0),
            (fit.model.e_inf, truth.e_inf),
            (fit.model.alpha, truth.alpha),
            (fit.model.tau, truth.tau),
        ] {
            assert!(
                ((got - want) / want).abs() < 1.0e-6,
                "parameter recovery: {got} vs {want}"
            );
        }
        // MUTATION: negated loss data describes an active material no
        // passive fractional Zener can represent; the fit must refuse by
        // name rather than return a silently wrong model.
        let mut corrupted = samples.clone();
        for s in &mut corrupted {
            s.2 = -s.2;
        }
        let err = fit_fractional_zener(&corrupted, &init, 200, 1.0e-3).unwrap_err();
        assert!(
            err.to_string().contains("FS-MAT-VISCO-FIT-DIVERGED"),
            "{err}"
        );
    }

    #[test]
    fn fit_refuses_degenerate_input_by_name() {
        let init = FractionalZener::new(1.0e9, 2.0e9, 0.5, 1.0e-3).expect("init");
        // Under-determined: fewer samples than parameters.
        let err = fit_fractional_zener(&[(1.0, 1.0e9, 1.0e7)], &init, 10, 1.0e-3).unwrap_err();
        assert!(err.to_string().contains("FS-MAT-VISCO-PARAMETERS"), "{err}");
        // Non-physical storage modulus.
        let bad = [
            (1.0, -1.0, 0.0),
            (2.0, 1.0, 0.0),
            (3.0, 1.0, 0.0),
            (4.0, 1.0, 0.0),
        ];
        let err = fit_fractional_zener(&bad, &init, 10, 1.0e-3).unwrap_err();
        assert!(err.to_string().contains("FS-MAT-VISCO-PARAMETERS"), "{err}");
    }

    #[test]
    fn fitted_and_lowered_models_keep_eta_nonnegative_and_relaxation_monotone() {
        // Property sweep: for admitted fractional models and their Prony
        // lowerings, η(ω) ≥ 0 across the certified band, every lowered
        // weight is non-negative, and the relaxation modulus
        // G(t) = E∞ + Σ E_j·exp(−t/τ_j) is non-increasing (G′ ≤ 0 follows
        // analytically from E_j ≥ 0; the grid check pins the
        // implementation to that analysis).
        for &(alpha, tau) in &[(0.2, 1.0e-4), (0.5, 1.0e-3), (0.9, 5.0e-3)] {
            let fz = FractionalZener::new(5.0e9, 1.2e10, alpha, tau).expect("fz");
            let lowered = lower_to_prony(&fz, 20.0, 2.0e4, 10, 0.05).expect("lowering");
            for k in 0..=64 {
                let t = f64::from(k) / 64.0;
                let w = det::exp(
                    det::ln(lowered.band.0)
                        + t * (det::ln(lowered.band.1) - det::ln(lowered.band.0)),
                )
                .clamp(lowered.band.0, lowered.band.1);
                assert!(fz.loss_factor(w) >= 0.0, "fractional eta at w={w}");
                assert!(
                    lowered.loss_factor_checked(w).expect("in band") >= 0.0,
                    "lowered eta at w={w}"
                );
            }
            for &(e, tau_j) in &lowered.model.terms {
                assert!(e >= 0.0 && tau_j > 0.0, "admitted Prony weights");
            }
            let g = |t: f64| {
                lowered.model.e_inf
                    + lowered
                        .model
                        .terms
                        .iter()
                        .map(|&(e, tj)| e * det::exp(-t / tj))
                        .sum::<f64>()
            };
            let mut prev = g(0.0);
            for k in 1..=100 {
                let t = 1.0e-7 * det::exp(f64::from(k) * 0.15);
                let cur = g(t);
                assert!(
                    cur <= prev * (1.0 + 1.0e-15),
                    "relaxation must be monotone: G({t}) = {cur} > {prev}"
                );
                prev = cur;
            }
        }
    }
}

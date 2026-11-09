//! Harmonic response of a retained mechanical pencil, under `exp(-i omega t)`.
//!
//! Solve `(K - omega^2 M - i omega C) q = F` without discarding any projected
//! mass, stiffness, or velocity-coupling entries. The real velocity operator
//! may contain both symmetric viscous damping and skew gyroscopic coupling.
//! Resonance, antiresonance, phase and
//! transfer between ports are properties of this forced system, not just its
//! free-decay eigenvalues. Forces are peak phasors, so cycle-averaged power is
//! `Re(F^H v)/2`, with `v = -i omega q`.
//!
//! Preparation uses the quadratic-pencil owner's exact binary coordinate/time
//! scaling, then fs-la complex LU. A prepared frequency reuses its factor for
//! multiple excitations. Exactly diagonal pencils use scalar division. Singular
//! lossless poles refuse; no fictitious damping or stiffness is added. At DC a
//! free coordinate has no unique displacement and refuses too.
//!
//! This is a particular periodic response of the supplied finite-dimensional
//! pencil. It does not certify modal truncation, radiation, material validity,
//! the startup transient, or a steady state for a lossless forced resonance.
//! PSD admission uses a diagonally equilibrated Jacobi check with a rounding
//! allowance; it is not an interval passivity certificate and never modifies
//! the source matrix.

use fs_exec::Cx;
use fs_la::eigen::{admit_jacobi_eigh, jacobi_eigh};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_blocked_with_cancel};
use fs_la::factor::cholesky;
use fs_math::c64::C64;

use crate::{ModalError, quadratic_fraction, quadratic_power_scale};

mod impedance;
use impedance::PreparedPortImpedance;
mod compact;
use compact::CompactPortFactor;
mod refinement;
pub use impedance::{HarmonicPortControls, HarmonicPortResponse, admit_passive_port_impedance};
mod frequency;
pub use frequency::{
    HarmonicFrequencyControls, HarmonicFrequencyCurvatureControls,
    HarmonicFrequencyCurvatureResponse, HarmonicFrequencyResponse,
    HarmonicPortFrequencyCurvatureResponse, HarmonicPortFrequencyResponse,
    harmonic_mean_power_curvature_with_cancel, harmonic_mean_power_parameter_jet_with_cancel,
    harmonic_quadratic_power_curvature_with_cancel,
    harmonic_quadratic_power_derivative_with_cancel,
    harmonic_quadratic_power_parameter_jet_with_cancel,
};
mod parameter;
pub use parameter::{
    HarmonicFrequencyParameterControls, HarmonicFrequencyParameterResponse,
    HarmonicParameterControls, HarmonicParameterResponse, HarmonicPortFrequencyParameterDirection,
    HarmonicPortFrequencyParameterResponse, HarmonicPortParameterDirection,
    HarmonicPortParameterResponse, harmonic_mean_power_derivative_with_cancel,
};
mod interference;
pub use interference::{
    HarmonicEnsembleControls, HarmonicFactorStatistics, HarmonicForceFactor, HarmonicForceTone,
    HarmonicImpedanceSample, HarmonicInterferenceControls, HarmonicInterferenceLoad,
    HarmonicInterferenceMoments, HarmonicReceiver, HarmonicReceiverControlEvaluation,
    HarmonicReceiverControlLimits, HarmonicReceiverControls, HarmonicReceiverOptimumError,
    HarmonicReceiverOptimumLimits, HarmonicReceiverOptimumRequest, HarmonicReceiverOptimumSolution,
    HarmonicReceiverPairMoments, HarmonicReceiverQuantity, MAX_ENSEMBLE_FORCE_TONES,
    MAX_INTERFERENCE_FACTORS, MAX_INTERFERENCE_RECEIVERS, MAX_INTERFERENCE_TONES,
    MAX_RECEIVER_CONTROL_SOURCES, PreparedHarmonicInterference,
    PreparedHarmonicInterferenceEnsemble, PreparedHarmonicReceiverControl,
    PreparedHarmonicReceiverOptimum, PreparedHarmonicReceivers,
};
mod supports;
pub use supports::{PreparedSupportedHarmonicResponse, SupportedHarmonicResponse};
/// Local driven-response extrema and half-magnitude-squared linewidths.
pub mod resonance;

/// Explicit size, numerical payload and solve-residual limits.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicResponseControls {
    /// Maximum retained coordinates, checked before allocation.
    pub maximum_coordinates: usize,
    /// Conservative peak bytes of owned numerical arrays, including preparation
    /// and one solve. Caller input, allocator overhead and thread stacks are excluded.
    /// Dense owners admit `144*n*n + 128*n`; exactly diagonal, unloaded modal
    /// owners admit `272*n` without allocating dense coefficient matrices.
    pub maximum_payload_bytes: usize,
    /// Maximum componentwise relative backward error of the dynamic solve.
    /// This is not a forward-error bound near an ill-conditioned resonance.
    pub maximum_relative_residual: f64,
}

/// Phasor motion and independent physical power diagnostics.
#[derive(Clone, Debug)]
pub struct HarmonicResponse {
    /// Peak generalized displacement amplitudes, in the original coordinates.
    pub displacement: Vec<C64>,
    /// Peak generalized velocity amplitudes, `-i omega displacement`.
    pub velocity: Vec<C64>,
    /// Maximum rowwise backward error `|Aq-F| / (|A||q|+|F|)`.
    pub relative_residual: f64,
    /// `Re(F^H v)/2` [W] for physical force/displacement coordinates.
    pub average_input_power_w: f64,
    /// `Re(v^H D v)/2` [W], with the complete `D=(C+C^T)/2`. The skew
    /// gyroscopic part does no work and is excluded before accumulation.
    /// Small signed roundoff is retained, never clipped into a passivity claim.
    pub average_viscous_power_w: f64,
    /// Physical load-port flow, effort and absorbed work, when a load is supplied.
    pub port_load: Option<HarmonicPortResponse>,
    /// Input power minus structural viscous and supplied load-port power [W].
    /// Absorbed boundary work is not automatically radiation or material heat.
    pub power_balance_residual_w: f64,
}

/// Transactional preparation/solve refusal; no source data is mutated.
#[derive(Clone, Debug, PartialEq)]
pub enum HarmonicResponseError {
    /// Invalid dimensions, domains, or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The physical observation window or its phase/product is unrepresentable.
    Observation(fs_fft::harmonic_window::WindowError),
    /// Cancellation observed before publication.
    Cancelled,
    /// Explicit count or numerical-payload budget exceeded.
    Budget {
        /// Resource being admitted.
        quantity: &'static str,
        /// Requested resource count.
        requested: usize,
        /// Caller-authored ceiling.
        cap: usize,
    },
    /// A harmonic factorization or same-factor frequency tangent exceeds its work cap.
    WorkBudget {
        /// Requested work score.
        requested: u64,
        /// Caller-authored maximum work score.
        cap: u64,
    },
    /// The common finite/symmetric damping admission rejected the source.
    Damping(ModalError),
    /// The mass matrix failed positive-definite admission.
    MassNotPositiveDefinite,
    /// Materially negative stiffness/damping direction after scaled PSD admission.
    NotPositiveSemidefinite {
        /// Failed coefficient matrix.
        matrix: &'static str,
        /// Minimum eigenvalue of its diagonally equilibrated image, when computed.
        minimum_scaled_eigenvalue: f64,
    },
    /// No unique harmonic response exists at this frequency. In particular,
    /// a lossless resonance cannot be assigned a finite steady-state amplitude.
    Singular {
        /// Physical driving frequency [rad/s].
        angular_frequency_rad_s: f64,
        /// Failing elimination pivot; compact loads may index their bordered core.
        pivot: usize,
    },
    /// A solve did not meet the caller's backward-error budget.
    Residual {
        /// Measured componentwise relative residual.
        measured: f64,
        /// Caller-authored ceiling.
        limit: f64,
    },
}

impl core::fmt::Display for HarmonicResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-HARMONIC: {self:?}")
    }
}
impl std::error::Error for HarmonicResponseError {}

/// One immutable dynamic-stiffness factor at one physical frequency.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicResponse {
    n: usize,
    layout: MatrixLayout,
    angular_frequency_rad_s: f64,
    scaled_frequency: f64,
    time_exponent: i32,
    coordinate_exponents: Vec<i32>,
    operator: Vec<C64>,
    scaled_mass: Vec<f64>,
    scaled_damping: Vec<f64>,
    damping: Vec<f64>,
    factor: Option<LuComplex>,
    compact_port_factor: Option<CompactPortFactor>,
    port_load: Option<PreparedPortImpedance>,
    maximum_relative_residual: f64,
}

// Structural coefficient storage. A compact physical port factor is separate.
// A small off-diagonal velocity entry is never a structural sparsity test.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MatrixLayout {
    Dense,
    Diagonal,
}

impl MatrixLayout {
    fn columns(self, row: usize, n: usize) -> core::ops::Range<usize> {
        match self {
            Self::Dense => 0..n,
            Self::Diagonal => row..row + 1,
        }
    }

    fn index(self, row: usize, column: usize, n: usize) -> usize {
        match self {
            Self::Dense => row * n + column,
            Self::Diagonal => row,
        }
    }
}

impl PreparedHarmonicResponse {
    /// Prepare a complete real row-major `M,C,K` pencil. `M` and `K` must be
    /// exactly symmetric; `M` must be positive definite, while `K` and the
    /// symmetric viscous part `D=(C+C^T)/2` pass numerical PSD admission.
    /// The skew part `G=(C-C^T)/2` is retained as gyroscopic coupling, including
    /// its sign: it changes resonances and reciprocity but dissipates no energy.
    /// Rigid coordinates and zero damping are permitted. Frequency
    /// is nonnegative [rad/s]. No mass whitening or diagonal projection is applied.
    ///
    /// Polls cancellation between dense admission stages, per assembly row,
    /// and within fs-la LU. Cholesky/Jacobi admission remains one indivisible
    /// bounded-size call; wall-time cancellation during it is not claimed.
    #[allow(clippy::too_many_arguments)] // complete physical pencil and explicit controls
    pub fn try_new(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        Self::prepare(
            mass,
            damping,
            stiffness,
            n,
            angular_frequency_rad_s,
            controls,
            cx,
            None,
            None,
            MatrixLayout::Dense,
            true,
        )
    }

    /// Prepare `(K - omega^2 M - i omega C - i omega B Z B^T) q = F`.
    /// `B` is the real row-major `n x ports` force map; the conjugate flow
    /// is `B^T v`. The full complex `Z` maps that flow to opposing load effort
    /// under `exp(-i omega t)`. All mutual and reactive entries are retained.
    ///
    /// The finite impedance's Hermitian part must pass pointwise numerical PSD
    /// admission. This is not a broadband causality/passivity certificate. The
    /// caller owns the source frequency band, physical port units and any fit.
    /// Load and structure are assembled before factorization: a physical load
    /// can resolve an otherwise singular bare resonance without inverting it.
    #[allow(clippy::too_many_arguments)] // physical pencil, ports and separate explicit budgets
    pub fn try_with_port_impedance(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        force_map: &[f64],
        ports: usize,
        impedance: &[C64],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        admit_size(n, controls)?;
        let load = PreparedPortImpedance::new(n, force_map, ports, impedance, port_controls, cx)?;
        Self::prepare(
            mass,
            damping,
            stiffness,
            n,
            angular_frequency_rad_s,
            controls,
            cx,
            None,
            Some(load),
            MatrixLayout::Dense,
            true,
        )
    }

    #[allow(clippy::too_many_arguments)] // same pencil, with optional exact modal frequencies
    fn prepare(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
        modal_frequencies: Option<&[f64]>,
        port_load: Option<PreparedPortImpedance>,
        layout: MatrixLayout,
        factorize: bool,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        admit_layout(n, controls, layout)?;
        let entries = match layout {
            MatrixLayout::Dense => n.checked_mul(n),
            MatrixLayout::Diagonal => Some(n),
        };
        if entries != Some(mass.len())
            || damping.len() != mass.len()
            || stiffness.len() != mass.len()
        {
            return Err(invalid("M,C,K dimensions must match their storage layout"));
        }
        if !angular_frequency_rad_s.is_finite() || angular_frequency_rad_s < 0.0 {
            return Err(invalid("finite nonnegative driving frequency required"));
        }
        if layout == MatrixLayout::Dense {
            symmetric(mass, n)?;
            symmetric(stiffness, n)?;
            validate_velocity_coupling(damping, n)?;
        } else if mass
            .iter()
            .chain(damping)
            .chain(stiffness)
            .any(|v| !v.is_finite())
            || damping.iter().chain(stiffness).any(|&v| v < 0.0)
        {
            return Err(invalid(
                "finite nonnegative diagonal modal coefficients required",
            ));
        }
        if (0..n).any(|i| mass[layout.index(i, i, n)] <= 0.0) {
            return Err(HarmonicResponseError::MassNotPositiveDefinite);
        }
        let mut coordinate_exponents = Vec::with_capacity(n);
        for i in 0..n {
            coordinate_exponents.push(
                quadratic_fraction(mass[layout.index(i, i, n)])
                    .1
                    .div_euclid(2),
            );
        }
        let mut time_exponent = if angular_frequency_rad_s == 0.0 {
            i32::MIN
        } else {
            quadratic_fraction(angular_frequency_rad_s).1
        };
        for i in 0..n {
            poll(cx)?;
            for j in layout.columns(i, n) {
                let index = layout.index(i, j, n);
                let coordinate = coordinate_exponents[i] + coordinate_exponents[j];
                if damping[index] != 0.0 {
                    time_exponent =
                        time_exponent.max(quadratic_fraction(damping[index]).1 - coordinate);
                }
                if stiffness[index] != 0.0 {
                    time_exponent = time_exponent.max(
                        (quadratic_fraction(stiffness[index]).1 - coordinate + 1).div_euclid(2),
                    );
                }
            }
        }
        // A finite velocity impedance contributes exactly zero at DC. In
        // particular it must not rescale an unrelated static stiffness solve.
        let port_scaling = if angular_frequency_rad_s == 0.0 {
            None
        } else {
            port_load
                .as_ref()
                .map(|load| load.scaling(&coordinate_exponents, cx))
                .transpose()?
        };
        if let Some(scaling) = &port_scaling {
            time_exponent = time_exponent.max(scaling.minimum_time_exponent);
        }
        if time_exponent == i32::MIN {
            time_exponent = 0;
        }
        let transform = |a: &[f64], order: i32| -> Result<Vec<f64>, HarmonicResponseError> {
            let mut transformed = Vec::with_capacity(a.len());
            for (index, &v) in a.iter().enumerate() {
                let (row, column) = match layout {
                    MatrixLayout::Dense => (index / n, index % n),
                    MatrixLayout::Diagonal => (index, index),
                };
                transformed.push(scale(
                    v,
                    -coordinate_exponents[row]
                        - coordinate_exponents[column]
                        - order * time_exponent,
                )?);
            }
            Ok(transformed)
        };
        let m = transform(mass, 0)?;
        let c = transform(damping, 1)?;
        let k = transform(stiffness, 2)?;
        poll(cx)?;
        if layout == MatrixLayout::Dense {
            cholesky(&m, n).map_err(|_| HarmonicResponseError::MassNotPositiveDefinite)?;
            require_psd(&c, n, "damping", cx)?;
            require_psd(&k, n, "stiffness", cx)?;
        }
        let mu = scale(angular_frequency_rad_s, -time_exponent)?;
        let mut operator = vec![C64::ZERO; mass.len()];
        for i in 0..n {
            poll(cx)?;
            for j in layout.columns(i, n) {
                let index = layout.index(i, j, n);
                // Compensate the first multiplication in mu*(mu*M), so a
                // near-resonant detuning is not lost in two rounded products.
                let product = mu * m[index];
                let error = mu.mul_add(m[index], -product);
                let mut real = (-mu).mul_add(product, k[index]) - mu * error;
                if i == j
                    && let Some(frequencies) = modal_frequencies
                {
                    // The modal API supplies omega, not a rounded omega^2.
                    // Factor its detuning so exact resonance stays exact.
                    real = scale(frequencies[i] - angular_frequency_rad_s, -time_exponent)?
                        * (scale(frequencies[i], -time_exponent)? + mu);
                }
                operator[index] = C64::new(real, -mu * c[index]);
                if !finite(operator[index]) {
                    return Err(invalid("dynamic stiffness is not representable"));
                }
            }
        }
        let compact_port_factor = match (&port_load, &port_scaling, layout) {
            (Some(load), Some(scaling), MatrixLayout::Diagonal) => Some(CompactPortFactor::new(
                &operator,
                load,
                &coordinate_exponents,
                time_exponent,
                mu,
                angular_frequency_rad_s,
                scaling,
                controls.maximum_relative_residual,
                cx,
            )?),
            _ => None,
        };
        if compact_port_factor.is_none()
            && let (Some(load), Some(scaling)) = (&port_load, &port_scaling)
        {
            load.add_to_operator(
                &mut operator,
                &coordinate_exponents,
                time_exponent,
                mu,
                scaling,
                cx,
            )?;
        }
        let mut prepared = Self {
            n,
            layout,
            angular_frequency_rad_s,
            scaled_frequency: mu,
            time_exponent,
            coordinate_exponents,
            operator,
            scaled_mass: m,
            scaled_damping: c,
            damping: damping.to_vec(),
            factor: None,
            compact_port_factor,
            port_load,
            maximum_relative_residual: controls.maximum_relative_residual,
        };
        if factorize {
            prepared.factorize(cx)?;
        }
        poll(cx)?;
        Ok(prepared)
    }

    // Supports first admit the complete physical pencil, then factor only its
    // free principal block. A rigid full-body mode must not prevent a supported
    // static solve, and the reaction rows must retain the original operator.
    fn factorize(&mut self, cx: &Cx<'_>) -> Result<(), HarmonicResponseError> {
        poll(cx)?;
        if self.compact_port_factor.is_some() {
            return Ok(());
        }
        let n = self.n;
        let diagonal = self.layout == MatrixLayout::Diagonal
            || (0..n).all(|i| (0..n).all(|j| i == j || self.operator[i * n + j] == C64::ZERO));
        if diagonal {
            for i in 0..n {
                if self.operator[self.layout.index(i, i, n)] == C64::ZERO {
                    return Err(HarmonicResponseError::Singular {
                        angular_frequency_rad_s: self.angular_frequency_rad_s,
                        pivot: i,
                    });
                }
            }
        } else {
            self.factor = Some(
                lu_complex_blocked_with_cancel(&self.operator, n, 1, || cx.checkpoint().is_err())
                    .map_err(|e| lu_error(e, self.angular_frequency_rad_s))?,
            );
        }
        poll(cx)
    }

    /// Convenience for mass-normalized modes: `M=I`, `K=diag(omega_i^2)`.
    /// The full supplied velocity matrix is retained, including viscous
    /// off-diagonals and skew gyroscopic terms, under the same PSD admission.
    /// Factored frequency detuning preserves exact resonance and narrow peaks.
    /// Exactly diagonal damping uses linear owned storage, solves and frequency
    /// tangents after inspecting the caller's complete damping matrix.
    pub fn try_modal(
        natural_frequencies_rad_s: &[f64],
        damping_per_s: &[f64],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        Self::prepare_modal(
            natural_frequencies_rad_s,
            damping_per_s,
            angular_frequency_rad_s,
            controls,
            cx,
            None,
        )
    }

    /// Prepare independent mass-normalized modes from two length-n vectors.
    /// Each damping entry is the finite nonnegative coefficient c [1/s] in
    /// q'' + c*q' + omega_i^2*q = F, not a damping ratio. The authored model
    /// has no off-diagonal damping or port load; no matrix is projected here.
    ///
    /// Input, preparation, owned storage, solves and frequency tangents all
    /// use O(n) space/work. The same `272*n` preparation-payload admission and
    /// frequency controls apply as for an exactly diagonal `try_modal` owner.
    /// Coordinate and payload bounds precede allocation; each mode polls
    /// cancellation. Zero frequencies/damping and singular-response refusal,
    /// representable squared frequencies, and factored detuning are shared
    /// with `try_modal`.
    pub fn try_diagonal_modal(
        natural_frequencies_rad_s: &[f64],
        damping_rates_per_s: &[f64],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        Self::prepare_modal_coefficients(
            natural_frequencies_rad_s,
            damping_rates_per_s,
            angular_frequency_rad_s,
            controls,
            cx,
            None,
            MatrixLayout::Diagonal,
        )
    }

    /// Independent modal vectors with a complete mutual physical port load.
    /// Damping entries are viscous coefficients [1/s], as in `try_diagonal_modal`.
    /// The full complex impedance and power-conjugate force map retain the
    /// same admission, physical units and slope requirements as the dense API.
    ///
    /// A small bordered LU retains exact and near bare poles before elimination;
    /// neither a bare zero diagonal nor the impedance is inverted. All modes
    /// and couplings remain. Load controls bound the retained core and streamed
    /// audits. Storage is O(n*p + k*k), preparation O(n*p*p + k*k*k), and a
    /// factor solve O(n*p + k*k), where k is ports plus retained pivot modes.
    /// The unchanged componentwise residual streams O(n*n*p) original entries.
    /// A difficult cluster refuses its budget or residual rather than dropping
    /// a mode. At DC the finite load vanishes from the operator but remains in
    /// the analytic frequency derivative.
    #[allow(clippy::too_many_arguments)] // physical vectors, ports and distinct resource controls
    pub fn try_diagonal_modal_with_port_impedance(
        natural_frequencies_rad_s: &[f64],
        damping_rates_per_s: &[f64],
        force_map: &[f64],
        ports: usize,
        impedance: &[C64],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        let n = natural_frequencies_rad_s.len();
        admit_layout(n, controls, MatrixLayout::Diagonal)?;
        if damping_rates_per_s.len() != n {
            return Err(invalid(
                "one modal damping rate per natural frequency required",
            ));
        }
        let load = PreparedPortImpedance::new(n, force_map, ports, impedance, port_controls, cx)?;
        Self::prepare_modal_coefficients(
            natural_frequencies_rad_s,
            damping_rates_per_s,
            angular_frequency_rad_s,
            controls,
            cx,
            Some(load),
            MatrixLayout::Diagonal,
        )
    }

    /// Mass-normalized modal convenience with the same complete physical load.
    /// Retains factored modal detuning at narrow and exact bare resonances.
    #[allow(clippy::too_many_arguments)] // original modal source and physical ports
    pub fn try_modal_with_port_impedance(
        natural_frequencies_rad_s: &[f64],
        damping_per_s: &[f64],
        force_map: &[f64],
        ports: usize,
        impedance: &[C64],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        admit_size(natural_frequencies_rad_s.len(), controls)?;
        let load = PreparedPortImpedance::new(
            natural_frequencies_rad_s.len(),
            force_map,
            ports,
            impedance,
            port_controls,
            cx,
        )?;
        Self::prepare_modal(
            natural_frequencies_rad_s,
            damping_per_s,
            angular_frequency_rad_s,
            controls,
            cx,
            Some(load),
        )
    }

    fn prepare_modal(
        natural_frequencies_rad_s: &[f64],
        damping_per_s: &[f64],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
        port_load: Option<PreparedPortImpedance>,
    ) -> Result<Self, HarmonicResponseError> {
        let n = natural_frequencies_rad_s.len();
        poll(cx)?;
        // Bound coordinates and the minimum possible payload before scanning
        // the full caller-owned damping input. No numerical coupling threshold.
        admit_layout(n, controls, MatrixLayout::Diagonal)?;
        validate_velocity_coupling(damping_per_s, n)?;
        let mut layout = MatrixLayout::Diagonal;
        for i in 0..n {
            poll(cx)?;
            if port_load.is_some() || (0..n).any(|j| i != j && damping_per_s[i * n + j] != 0.0) {
                layout = MatrixLayout::Dense;
                break;
            }
        }
        let mut diagonal_damping = if layout == MatrixLayout::Diagonal {
            Vec::with_capacity(n)
        } else {
            Vec::new()
        };
        if layout == MatrixLayout::Diagonal {
            for i in 0..n {
                poll(cx)?;
                diagonal_damping.push(damping_per_s[i * n + i]);
            }
        }
        Self::prepare_modal_coefficients(
            natural_frequencies_rad_s,
            if layout == MatrixLayout::Diagonal {
                &diagonal_damping
            } else {
                damping_per_s
            },
            angular_frequency_rad_s,
            controls,
            cx,
            port_load,
            layout,
        )
    }

    // Both authored diagonal vectors and the full-matrix modal convenience
    // share frequency admission, M/K assembly and the same physical solver.
    fn prepare_modal_coefficients(
        natural_frequencies_rad_s: &[f64],
        damping_per_s: &[f64],
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
        port_load: Option<PreparedPortImpedance>,
        layout: MatrixLayout,
    ) -> Result<Self, HarmonicResponseError> {
        let n = natural_frequencies_rad_s.len();
        poll(cx)?;
        admit_layout(n, controls, layout)?;
        let entries = if layout == MatrixLayout::Diagonal {
            n
        } else {
            n * n
        };
        if damping_per_s.len() != entries {
            return Err(invalid(
                "modal damping dimensions must match its storage layout",
            ));
        }
        let mut m = vec![0.0; entries];
        let mut k = vec![0.0; entries];
        for (i, &w) in natural_frequencies_rad_s.iter().enumerate() {
            poll(cx)?;
            if !w.is_finite() || w < 0.0 || !(w * w).is_finite() || w != 0.0 && w * w == 0.0 {
                return Err(invalid(
                    "finite nonnegative modal frequencies with representable stiffness required",
                ));
            }
            if layout == MatrixLayout::Diagonal
                && (!damping_per_s[i].is_finite() || damping_per_s[i] < 0.0)
            {
                return Err(invalid("finite nonnegative modal damping rates required"));
            }
            let index = layout.index(i, i, n);
            m[index] = 1.0;
            k[index] = w * w;
        }
        Self::prepare(
            &m,
            damping_per_s,
            &k,
            n,
            angular_frequency_rad_s,
            controls,
            cx,
            Some(natural_frequencies_rad_s),
            port_load,
            layout,
            true,
        )
    }

    /// Solve a new peak-force phasor using the prepared frequency's factor.
    /// Generalized forces and outputs remain in the source pencil's coordinates.
    /// Up to two residual corrections reuse that factor when needed to meet
    /// the unchanged caller-authored componentwise backward-error threshold.
    /// Returns no partial response on cancellation, overflow or residual refusal.
    pub fn solve(
        &self,
        force: &[C64],
        cx: &Cx<'_>,
    ) -> Result<HarmonicResponse, HarmonicResponseError> {
        let (solution, force_exponent, relative_residual) = self.solve_force(force, cx)?;
        self.observe(force, &solution, force_exponent, relative_residual, cx)
    }

    fn solve_force(
        &self,
        force: &[C64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<C64>, i32, f64), HarmonicResponseError> {
        poll(cx)?;
        if force.len() != self.n || force.iter().any(|z| !finite(*z)) {
            return Err(invalid("one finite force phasor per coordinate required"));
        }
        // Normalize each RHS as a whole, independently of the source units or
        // response size. Output scaling restores the exact physical exponents.
        let mut force_exponent = i32::MIN;
        for (i, f) in force.iter().enumerate() {
            for v in [f.re, f.im] {
                if v != 0.0 {
                    force_exponent =
                        force_exponent.max(quadratic_fraction(v).1 - self.coordinate_exponents[i]);
                }
            }
        }
        if force_exponent == i32::MIN {
            force_exponent = 0;
        }
        let mut rhs = Vec::with_capacity(self.n);
        for (i, &f) in force.iter().enumerate() {
            rhs.push(scale_complex(
                f,
                -self.coordinate_exponents[i] - force_exponent,
            )?);
        }
        let (solution, residual) = self.solve_normalized(&rhs, cx)?;
        Ok((solution, force_exponent, residual))
    }

    fn solve_normalized(
        &self,
        rhs: &[C64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<C64>, f64), HarmonicResponseError> {
        poll(cx)?;
        let mut solution = rhs.to_vec();
        self.solve_factored(&mut solution, cx)?;
        let mut relative_residual = backward_error(self, &solution, rhs, cx)?;
        for _ in 0..2 {
            if relative_residual <= self.maximum_relative_residual
                || self.factor.is_none() && self.compact_port_factor.is_none()
            {
                break;
            }
            refinement::correct(self, &mut solution, rhs, cx)?;
            relative_residual = backward_error(self, &solution, rhs, cx)?;
        }
        if relative_residual > self.maximum_relative_residual {
            return Err(HarmonicResponseError::Residual {
                measured: relative_residual,
                limit: self.maximum_relative_residual,
            });
        }
        poll(cx)?;
        Ok((solution, relative_residual))
    }

    // The immutable LU/compact/scalar dispatch serves primal, tangent and
    // bounded residual-correction solves without any new factorization.
    fn solve_factored(
        &self,
        solution: &mut [C64],
        cx: &Cx<'_>,
    ) -> Result<(), HarmonicResponseError> {
        poll(cx)?;
        if let Some(compact) = &self.compact_port_factor {
            compact.solve(&self.operator, solution, cx)?;
        } else if let Some(factor) = &self.factor {
            factor
                .solve_with_cancel(solution, || cx.checkpoint().is_err())
                .map_err(|e| lu_error(e, self.angular_frequency_rad_s))?;
        } else {
            for (i, value) in solution.iter_mut().enumerate() {
                poll(cx)?;
                *value = *value / self.operator[self.layout.index(i, i, self.n)];
            }
        }
        if solution.iter().any(|z| !finite(*z)) {
            return Err(invalid("harmonic solution is not representable"));
        }
        poll(cx)
    }

    fn observe(
        &self,
        force: &[C64],
        solution: &[C64],
        force_exponent: i32,
        relative_residual: f64,
        cx: &Cx<'_>,
    ) -> Result<HarmonicResponse, HarmonicResponseError> {
        let mut displacement = Vec::with_capacity(self.n);
        let mut velocity = Vec::with_capacity(self.n);
        for (i, &value) in solution.iter().enumerate() {
            poll(cx)?;
            let exponent = force_exponent - self.coordinate_exponents[i];
            displacement.push(scale_complex(value, exponent - 2 * self.time_exponent)?);
            velocity.push(scale_complex(
                C64::new(value.im, -value.re).scale(self.scaled_frequency),
                exponent - self.time_exponent,
            )?);
        }
        let mut input = 0.0;
        let mut dissipation = 0.0;
        for i in 0..self.n {
            poll(cx)?;
            input += 0.5 * real_inner(force[i], velocity[i]);
            let mut cv = C64::ZERO;
            for j in self.layout.columns(i, self.n) {
                cv = cv + velocity[j].scale(self.viscous_entry(i, j));
            }
            dissipation += 0.5 * real_inner(velocity[i], cv);
        }
        let port_load = self
            .port_load
            .as_ref()
            .map(|load| load.response(&velocity, cx))
            .transpose()?;
        let load_power = port_load
            .as_ref()
            .map_or(0.0, |load| load.average_absorbed_power_w);
        let balance = input - dissipation - load_power;
        if ![input, dissipation, balance]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("physical harmonic power is not representable"));
        }
        poll(cx)?;
        Ok(HarmonicResponse {
            displacement,
            velocity,
            relative_residual,
            average_input_power_w: input,
            average_viscous_power_w: dissipation,
            port_load,
            power_balance_residual_w: balance,
        })
    }

    fn operator_columns(&self, row: usize) -> core::ops::Range<usize> {
        if self.compact_port_factor.is_some() {
            0..self.n
        } else {
            self.layout.columns(row, self.n)
        }
    }

    fn operator_entry(&self, row: usize, column: usize) -> Result<C64, HarmonicResponseError> {
        if let Some(compact) = &self.compact_port_factor {
            compact.entry(&self.operator, row, column)
        } else {
            Ok(self.operator[self.layout.index(row, column, self.n)])
        }
    }

    fn viscous_entry(&self, row: usize, column: usize) -> f64 {
        match self.layout {
            MatrixLayout::Dense => symmetric_entry(&self.damping, self.n, row, column),
            MatrixLayout::Diagonal => self.damping[row],
        }
    }
}

fn admit_size(n: usize, controls: HarmonicResponseControls) -> Result<(), HarmonicResponseError> {
    admit_layout(n, controls, MatrixLayout::Dense)
}
fn admit_layout(
    n: usize,
    controls: HarmonicResponseControls,
    layout: MatrixLayout,
) -> Result<(), HarmonicResponseError> {
    if n == 0
        || !controls.maximum_relative_residual.is_finite()
        || controls.maximum_relative_residual <= 0.0
        || controls.maximum_relative_residual >= 1.0
    {
        return Err(invalid(
            "nonempty pencil and relative residual budget in (0,1) required",
        ));
    }
    if n > controls.maximum_coordinates {
        return Err(HarmonicResponseError::Budget {
            quantity: "coordinates",
            requested: n,
            cap: controls.maximum_coordinates,
        });
    }
    let bytes = match layout {
        MatrixLayout::Dense => n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(144))
            .and_then(|x| n.checked_mul(128).and_then(|y| x.checked_add(y))),
        // Compact owner: 44*n bytes (three real diagonals, one complex
        // diagonal, coordinate exponents). Modal input construction and its
        // scaled stiffness peak at another 32*n; one primal response peaks
        // at another 48*n. Exact-capacity vectors keep both below 272*n.
        MatrixLayout::Diagonal => n.checked_mul(272),
    }
    .ok_or(invalid("harmonic payload extent overflow"))?;
    if bytes > controls.maximum_payload_bytes {
        return Err(HarmonicResponseError::Budget {
            quantity: "numerical payload bytes",
            requested: bytes,
            cap: controls.maximum_payload_bytes,
        });
    }
    if layout == MatrixLayout::Dense {
        admit_jacobi_eigh(n)
            .map_err(|_| invalid("dense PSD admission extent exceeds fs-la work cap"))?;
    }
    Ok(())
}
fn validate_velocity_coupling(damping: &[f64], n: usize) -> Result<(), HarmonicResponseError> {
    if n.checked_mul(n) != Some(damping.len()) {
        return crate::damping::validate_dense_damping(damping, n)
            .map_err(HarmonicResponseError::Damping);
    }
    if damping.iter().any(|v| !v.is_finite()) {
        return Err(invalid("finite velocity-coupling entries required"));
    }
    // Retain the original symmetric material-damping refusal semantics.
    // A nonsymmetric source is admitted through its viscous part later,
    // without modifying the full operator used for motion and tangents.
    if (0..n).all(|i| (0..n).all(|j| damping[i * n + j] == damping[j * n + i])) {
        crate::damping::validate_dense_damping(damping, n)
            .map_err(HarmonicResponseError::Damping)?;
    }
    Ok(())
}

fn symmetric(a: &[f64], n: usize) -> Result<(), HarmonicResponseError> {
    for i in 0..n {
        for j in 0..n {
            if !a[i * n + j].is_finite() || a[i * n + j] != a[j * n + i] {
                return Err(invalid(
                    "finite exactly symmetric mass and stiffness required",
                ));
            }
        }
    }
    Ok(())
}
fn require_psd(
    a: &[f64],
    n: usize,
    matrix: &'static str,
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    poll(cx)?;
    let refusal = |value| HarmonicResponseError::NotPositiveSemidefinite {
        matrix,
        minimum_scaled_eigenvalue: value,
    };
    let mut roots = Vec::with_capacity(n);
    for i in 0..n {
        // A zero PSD diagonal requires an exactly zero row. Preserve true
        // nullspaces; neither a global norm nor diagonal jitter can hide these.
        if a[i * n + i] < 0.0
            || a[i * n + i] == 0.0 && (0..n).any(|j| symmetric_entry(a, n, i, j) != 0.0)
        {
            return Err(refusal(a[i * n + i]));
        }
        roots.push(a[i * n + i].sqrt());
    }
    if (0..n).all(|i| (0..n).all(|j| i == j || symmetric_entry(a, n, i, j) == 0.0)) {
        return Ok(());
    }
    let mut normalized = vec![0.0; n * n];
    let allowance = 64.0 * n as f64 * f64::EPSILON;
    for i in 0..n {
        poll(cx)?;
        if roots[i] == 0.0 {
            continue;
        }
        normalized[i * n + i] = 1.0;
        for j in 0..i {
            if roots[j] == 0.0 {
                continue;
            }
            let value = (symmetric_entry(a, n, i, j) / roots[i]) / roots[j];
            if !value.is_finite() || value.abs() > 1.0 + allowance {
                return Err(refusal(1.0 - value.abs()));
            }
            normalized[i * n + j] = value;
            normalized[j * n + i] = value;
        }
    }
    let (values, _) = jacobi_eigh(&normalized, n);
    poll(cx)?;
    let minimum = values.iter().copied().fold(f64::INFINITY, f64::min);
    if !minimum.is_finite() || minimum < -allowance {
        return Err(refusal(minimum));
    }
    Ok(())
}
// The full C stays in the motion operator; only its symmetric part enters
// passivity admission and physical heat. Equal pairs retain their exact bits,
// including subnormals, and the alternate sum avoids overflow for large pairs.
fn symmetric_entry(a: &[f64], n: usize, i: usize, j: usize) -> f64 {
    let x = a[i * n + j];
    let y = a[j * n + i];
    if x == y {
        x
    } else {
        let sum = x + y;
        if sum.is_finite() {
            0.5 * sum
        } else {
            0.5 * x + 0.5 * y
        }
    }
}
fn backward_error(
    prepared: &PreparedHarmonicResponse,
    x: &[C64],
    b: &[C64],
    cx: &Cx<'_>,
) -> Result<f64, HarmonicResponseError> {
    let scale = x
        .iter()
        .chain(b)
        .fold(0.0_f64, |s, v| s.max(v.re.abs()).max(v.im.abs()));
    if scale == 0.0 {
        return Ok(0.0);
    }
    let normalized = |v: C64| C64::new(v.re / scale, v.im / scale);
    let mut maximum = 0.0_f64;
    for i in 0..prepared.n {
        poll(cx)?;
        let mut residual = -normalized(b[i]);
        let mut denominator = normalized(b[i]).abs();
        for j in prepared.operator_columns(i) {
            let value = normalized(x[j]);
            let coefficient = prepared.operator_entry(i, j)?;
            residual = residual + coefficient * value;
            denominator += coefficient.abs() * value.abs();
        }
        let error = if denominator == 0.0 {
            0.0
        } else {
            residual.abs() / denominator
        };
        if !error.is_finite() {
            return Err(invalid("harmonic residual is not representable"));
        }
        maximum = maximum.max(error);
    }
    Ok(maximum)
}
fn real_inner(a: C64, b: C64) -> f64 {
    a.re.mul_add(b.re, a.im * b.im)
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn invalid(what: &'static str) -> HarmonicResponseError {
    HarmonicResponseError::Invalid(what)
}
fn poll(cx: &Cx<'_>) -> Result<(), HarmonicResponseError> {
    cx.checkpoint()
        .map_err(|_| HarmonicResponseError::Cancelled)
}
fn scale(value: f64, exponent: i32) -> Result<f64, HarmonicResponseError> {
    let scaled = quadratic_power_scale(value, exponent);
    if !scaled.is_finite() || value != 0.0 && scaled == 0.0 {
        Err(invalid("binary scaled response is not representable"))
    } else {
        Ok(scaled)
    }
}
fn scale_complex(z: C64, exponent: i32) -> Result<C64, HarmonicResponseError> {
    Ok(C64::new(scale(z.re, exponent)?, scale(z.im, exponent)?))
}
fn lu_error(error: LuComplexError, angular_frequency_rad_s: f64) -> HarmonicResponseError {
    match error {
        LuComplexError::Cancelled => HarmonicResponseError::Cancelled,
        LuComplexError::Singular { index } => HarmonicResponseError::Singular {
            angular_frequency_rad_s,
            pivot: index,
        },
    }
}

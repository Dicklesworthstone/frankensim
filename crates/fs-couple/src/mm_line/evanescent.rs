//! Passive causal higher-mode ports below their local cutoffs.
//!
//! For a lossless uniform circular guide, the independently matched radial
//! mode has `Z(s) = z0 s/sqrt(s^2 + wc^2)` in the causal `exp(+i omega t)`
//! convention. It has the positive spectral representation
//! `Z(s) = (2 z0 wc/pi) integral_0^(pi/2) s/(wc^2+s^2 cos^2(theta)) dtheta`.
//! Midpoint Gauss-Chebyshev quadrature gives a series chain of parallel LC
//! cells: `L_j=z0/(N wc)`, `C_j=N cos^2(theta_j)/(z0 wc)`,
//! `theta_j=(j+1/2)pi/(2N)`. Every coefficient follows from the original gas,
//! radius and radial eigenvalue; there is no fit, damping or gain correction.
//!
//! The shared physical-cell midpoint owner closes the *current* wave sample.
//! Its identity is `delta H = dt (a^2-b^2)/z0`, for arbitrary real inputs.
//! Its physical frequency accuracy is qualified only on the authored band
//! below cutoff. Finite LC quadrature does not reproduce propagating radiation
//! above cutoff. The receipts bound both quadrature and Tustin warping; this is
//! not an exact-ZOH or an all-frequency square-root realization. No propagation
//! delay is hidden in processing latency. The higher modes here are matched
//! at the far end, so their input impedance is independent of cylinder length.
//!
//! Quadrature identity: for `r=omega/wc<1`, let
//! `q=(r/(1+sqrt(1-r^2)))^2`. The relative impedance error is exactly
//! `2 q^(2N)/(1+q^(2N))`. The Mobius map to reflection bounds its absolute
//! complex error by `2 r/sqrt(1-r^2)` times that error. Midpoint substitutes
//! `omega_tilde=2 tan(omega dt/2)/dt`. Since the exact reflection is
//! `-exp(-2i asin(r))`, the additional complex error is
//! `2 sin(asin(omega_tilde/wc)-asin(omega/wc))`. All these bounds increase
//! over the qualified band. They are exact-arithmetic approximation bounds,
//! not interval certificates for floating-point roundoff.

use super::{MmComboReceipt, MmLineBank, MmLineConfig, MmLineError, MmLoad};
use fs_duct::modal::{MAX_MODES, j1_roots};
use fs_duct::{Duct, DuctError, LossModel, Segment, Termination};
use fs_material::gas::GasState;
use fs_math::det;
use fs_vfit::discretize::{DelayedFilter, DiscretizeError};
use fs_vfit::relaxation::{
    MAX_PARALLEL_CELLS, ParallelCellFrame, ParallelCellLoad, ParallelCellSpec, ParallelLcCell,
};
use fs_vfit::scattering::ScatteringResponse;
use fs_vfit::waveguide::WaveguideError;

/// Authored approximation domain and refinement level; no hidden defaults.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EvanescentConfig {
    /// Positive spectral quadrature order, at most [`MAX_PARALLEL_CELLS`].
    pub cells_per_mode: usize,
    /// Upper angular frequency [rad/s]; the warped band must stay below cutoff.
    pub band_max_rad_s: f64,
    /// Maximum absolute complex reflection error on the band. This bounds
    /// quadrature plus time warping in exact arithmetic, excluding roundoff.
    pub max_reflection_error: f64,
}

/// Source-derived bounds for one retained radial mode (plane mode excluded).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EvanescentReceipt {
    /// Original radial mode index, starting at one.
    pub mode: usize,
    /// Physical local cutoff `c*j1_root/radius` [rad/s].
    pub cutoff_rad_s: f64,
    /// Relative impedance quadrature error, maximized over the warped band.
    pub impedance_relative_error_bound: f64,
    /// Absolute complex reflection error caused by finite quadrature.
    pub quadrature_reflection_error_bound: f64,
    /// Absolute complex reflection error caused by midpoint frequency warping.
    pub time_reflection_error_bound: f64,
    /// Sum of the two reflection bounds on the authored physical band.
    pub reflection_error_bound: f64,
}

/// Typed refusals; an inadmissible sample never advances accepted history.
#[derive(Debug)]
pub enum EvanescentError {
    /// Invalid or unsupported geometry, domain, clock or sample.
    Invalid(&'static str),
    /// The original duct authority refused the physical source.
    Duct(DuctError),
    /// The existing plane-mode FIR realization refused.
    Plane(MmLineError),
    /// The existing plane filter refused a sample or pending response.
    Discrete(DiscretizeError),
    /// The physical cell owner refused a coefficient, state or work value.
    Cell(WaveguideError),
    /// The authored resolution and clock cannot meet the requested band error.
    Accuracy {
        /// Original radial mode index.
        mode: usize,
        /// Analytic complex reflection error bound.
        bound: f64,
        /// Caller-authored maximum.
        tolerance: f64,
    },
}

impl core::fmt::Display for EvanescentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "evanescent port: {what}"),
            Self::Duct(error) => write!(f, "evanescent source: {error}"),
            Self::Plane(error) => write!(f, "uniform modal plane: {error}"),
            Self::Discrete(error) => write!(f, "uniform modal plane: {error}"),
            Self::Cell(error) => write!(f, "evanescent storage: {error}"),
            Self::Accuracy {
                mode,
                bound,
                tolerance,
            } => write!(
                f,
                "evanescent mode {mode}: reflection error bound {bound:e} exceeds {tolerance:e}"
            ),
        }
    }
}
impl core::error::Error for EvanescentError {}

/// Independently matched higher ports; each uses the same positive real `z0`
/// wave reference as [`super::MatrixFirLine`]. At most seven by sixteen cells.
/// Hot stepping is allocation-free and retains uncorrected storage/work frames.
pub struct MatchedEvanescentPorts {
    gas: GasState,
    radius: f64,
    z0: f64,
    dt: f64,
    config: EvanescentConfig,
    receipts: Vec<EvanescentReceipt>,
    loads: Vec<ParallelCellLoad>,
    pending_loads: Vec<ParallelCellLoad>,
    incoming: Vec<f64>,
    pending_incoming: Vec<f64>,
    frames: Vec<ParallelCellFrame>,
    pending_frames: Vec<ParallelCellFrame>,
}

impl MatchedEvanescentPorts {
    /// Realize modes `1..n_modes` from the original gas and circular radius.
    /// `n_modes=1` retains no higher ports. Construction is bounded by the
    /// duct mode cap and the existing physical-cell capacity.
    ///
    /// # Errors
    /// Invalid source/clock/order, warped band reaching a cutoff, or an error
    /// bound exceeding the authored tolerance. No best-effort fallback exists.
    pub fn new(
        gas: &GasState,
        radius: f64,
        n_modes: usize,
        dt: f64,
        config: EvanescentConfig,
    ) -> Result<Self, EvanescentError> {
        if n_modes == 0
            || n_modes > MAX_MODES
            || config.cells_per_mode == 0
            || config.cells_per_mode > MAX_PARALLEL_CELLS
            || !dt.is_finite()
            || dt <= 0.0
            || !config.band_max_rad_s.is_finite()
            || config.band_max_rad_s <= 0.0
            || !config.max_reflection_error.is_finite()
            || config.max_reflection_error <= 0.0
            || config.max_reflection_error > 2.0
        {
            return Err(EvanescentError::Invalid(
                "invalid mode/order, clock or approximation band",
            ));
        }
        let half_band = 0.5 * config.band_max_rad_s * dt;
        if !(half_band > 0.0 && half_band < core::f64::consts::FRAC_PI_2) {
            return Err(EvanescentError::Invalid(
                "the authored band must be below Nyquist",
            ));
        }
        let warped = config.band_max_rad_s * (det::tan(half_band) / half_band);
        let wave = fs_duct::segment_wave(gas, radius, config.band_max_rad_s, LossModel::Lossless)
            .map_err(EvanescentError::Duct)?;
        let z0 = wave.characteristic_impedance.re;
        if !warped.is_finite() || !z0.is_finite() || z0 <= 0.0 {
            return Err(EvanescentError::Invalid(
                "source impedance or warped band is not finite",
            ));
        }
        let roots = j1_roots(n_modes - 1).map_err(EvanescentError::Duct)?;
        let mut loads = Vec::with_capacity(roots.len());
        let mut receipts = Vec::with_capacity(roots.len());
        let n = config.cells_per_mode;
        for (index, root) in roots.into_iter().enumerate() {
            let cutoff = gas.sound_speed * root / radius;
            if !cutoff.is_finite() || cutoff <= 0.0 || warped >= cutoff {
                return Err(EvanescentError::Invalid(
                    "the warped band reaches a higher-mode cutoff",
                ));
            }
            let (r, t) = (config.band_max_rad_s / cutoff, warped / cutoff);
            let (sr, st) = (
                det::sqrt((1.0 - r) * (1.0 + r)),
                det::sqrt((1.0 - t) * (1.0 + t)),
            );
            let q = (t / (1.0 + st)).powi(2);
            let qn = det::powi(q, (2 * n) as i32);
            let relative = 2.0 * qn / (1.0 + qn);
            let quadrature = 2.0 * (t / st) * relative;
            let time = 2.0 * (t * sr - r * st).abs();
            let bound = quadrature + time;
            if !bound.is_finite() || bound > config.max_reflection_error {
                return Err(EvanescentError::Accuracy {
                    mode: index + 1,
                    bound,
                    tolerance: config.max_reflection_error,
                });
            }
            let mut cells = [ParallelLcCell {
                inertance_pa_s2_m3: 0.0,
                compliance_m3_pa: 0.0,
            }; MAX_PARALLEL_CELLS];
            for (j, cell) in cells[..n].iter_mut().enumerate() {
                let theta = (j as f64 + 0.5) * core::f64::consts::PI / (2 * n) as f64;
                let cosine = det::cos(theta);
                *cell = ParallelLcCell {
                    inertance_pa_s2_m3: z0 / (n as f64 * cutoff),
                    compliance_m3_pa: n as f64 * cosine * cosine / (z0 * cutoff),
                };
            }
            let spec = ParallelCellSpec::lossless_lc(&cells[..n]).map_err(EvanescentError::Cell)?;
            let load = ParallelCellLoad::new(spec, dt).map_err(EvanescentError::Cell)?;
            if !(load.affine_pressure().0 + z0).is_finite() {
                return Err(EvanescentError::Invalid(
                    "wave closure is not representable",
                ));
            }
            loads.push(load);
            receipts.push(EvanescentReceipt {
                mode: index + 1,
                cutoff_rad_s: cutoff,
                impedance_relative_error_bound: relative,
                quadrature_reflection_error_bound: quadrature,
                time_reflection_error_bound: time,
                reflection_error_bound: bound,
            });
        }
        let count = loads.len();
        Ok(Self {
            gas: *gas,
            radius,
            z0,
            dt,
            config,
            receipts,
            pending_loads: loads.clone(),
            loads,
            incoming: vec![0.0; count],
            pending_incoming: vec![0.0; count],
            frames: vec![ParallelCellFrame::default(); count],
            pending_frames: vec![ParallelCellFrame::default(); count],
        })
    }

    /// Original physical gas; no reidentification occurs.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
    /// Original circular radius [m].
    #[must_use]
    pub const fn radius_m(&self) -> f64 {
        self.radius
    }
    /// Shared real characteristic-wave impedance [Pa s/m^3].
    #[must_use]
    pub const fn zc0(&self) -> f64 {
        self.z0
    }
    /// Actual sample interval [s].
    #[must_use]
    pub const fn dt(&self) -> f64 {
        self.dt
    }
    /// Authored approximation settings, unchanged.
    #[must_use]
    pub const fn config(&self) -> EvanescentConfig {
        self.config
    }
    /// Source-derived per-mode approximation bounds.
    #[must_use]
    pub fn receipts(&self) -> &[EvanescentReceipt] {
        &self.receipts
    }
    /// The actual LC declaration for higher-port index `i` (radial mode `i+1`).
    #[must_use]
    pub fn cell_spec(&self, i: usize) -> Option<&ParallelCellSpec> {
        self.loads.get(i).map(ParallelCellLoad::spec)
    }
    /// Last accepted higher-mode waves, in radial mode order starting at one.
    #[must_use]
    pub fn incoming(&self) -> &[f64] {
        &self.incoming
    }
    /// Last accepted physical cell observations, without energy correction.
    #[must_use]
    pub fn frames(&self) -> &[ParallelCellFrame] {
        &self.frames
    }
    /// Actual reactive storage [J], summed over all higher modes.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.loads
            .iter()
            .map(ParallelCellLoad::stored_energy_j)
            .sum()
    }

    /// Current-sample affine incoming wave at higher-port index `i` (mode i+1).
    /// An implicit source solves with this before accepting the same sample.
    pub fn boundary_response(&self, i: usize) -> Result<ScatteringResponse, EvanescentError> {
        let load = self
            .loads
            .get(i)
            .ok_or(EvanescentError::Invalid("higher-port index out of range"))?;
        let (direct, history) = load.affine_pressure();
        let divisor = direct + self.z0;
        let response = ScatteringResponse {
            direct: (direct - self.z0) / divisor,
            history: (self.z0 / divisor) * history,
        };
        if !response.direct.is_finite() || !response.history.is_finite() {
            return Err(EvanescentError::Invalid(
                "pending higher-mode wave is not finite",
            ));
        }
        Ok(response)
    }

    fn prepare(&mut self, outgoing: &[f64]) -> Result<(), EvanescentError> {
        if outgoing.len() != self.loads.len() || outgoing.iter().any(|value| !value.is_finite()) {
            return Err(EvanescentError::Invalid(
                "one finite outgoing wave per higher mode is required",
            ));
        }
        let mut total_storage = 0.0;
        for (i, &a) in outgoing.iter().enumerate() {
            let (direct, history) = self.loads[i].affine_pressure();
            let flow = (2.0 * a - history) / (direct + self.z0);
            let b = a - self.z0 * flow;
            if !b.is_finite() {
                return Err(EvanescentError::Invalid(
                    "higher-mode wave closure is not finite",
                ));
            }
            self.pending_loads[i] = self.loads[i];
            self.pending_frames[i] = self.pending_loads[i]
                .step(flow)
                .map_err(EvanescentError::Cell)?;
            total_storage += self.pending_frames[i].stored_energy_j;
            if !total_storage.is_finite() {
                return Err(EvanescentError::Invalid(
                    "total modal storage is not finite",
                ));
            }
            self.pending_incoming[i] = b;
        }
        Ok(())
    }

    fn accept(&mut self) {
        self.loads.copy_from_slice(&self.pending_loads);
        self.frames.copy_from_slice(&self.pending_frames);
        self.incoming.copy_from_slice(&self.pending_incoming);
    }

    /// Advance every higher port on this sample. Storage can release energy:
    /// reflected amplitudes are never clipped. No hot-path allocation.
    ///
    /// # Errors
    /// Wrong count or nonfinite/unrepresentable waves, states, storage or work.
    /// Every accepted state and observation is preserved on failure.
    pub fn push(&mut self, outgoing: &[f64]) -> Result<&[f64], EvanescentError> {
        self.prepare(outgoing)?;
        self.accept();
        Ok(&self.incoming)
    }
}

/// One uniform lossless cylinder with its original plane-mode FIR and passive
/// causal independently matched higher ports. The geometry has no mode
/// conversion. Scalar radiation/table loads retain precisely their disclosed
/// matched-higher-mode approximation; this is not a full aperture matrix.
/// The plane FIR keeps its original finite-grid admission receipt. Only the
/// higher ports have the exact LC storage/work identity.
pub struct UniformModalLine {
    plane: DelayedFilter,
    plane_receipt: MmComboReceipt,
    higher: MatchedEvanescentPorts,
    incoming: Vec<f64>,
}

impl UniformModalLine {
    /// Explicit opt-in to the qualified LC realization. No existing FIR route
    /// silently substitutes this approximation. Original plane taps and delay
    /// remain owned by the admitted one-mode [`MmLineBank`]'s filter.
    ///
    /// # Errors
    /// Requires exactly one lossless cylinder and scalar radiation or a strict
    /// table. Closed, pressure-release, full-modal, matched-plane and diagnostic
    /// fallback loads refuse: their higher-port boundary is a different problem.
    pub fn new(
        duct: &Duct,
        gas: &GasState,
        load: MmLoad<'_>,
        config: &MmLineConfig,
        approximation: EvanescentConfig,
    ) -> Result<Self, EvanescentError> {
        let [Segment::Cylinder { radius, length }] = duct.segments.as_slice() else {
            return Err(EvanescentError::Invalid(
                "a single uniform cylinder is required",
            ));
        };
        if !matches!(config.loss, LossModel::Lossless) || !length.is_finite() || *length <= 0.0 {
            return Err(EvanescentError::Invalid(
                "a positive-length lossless cylinder is required",
            ));
        }
        if !matches!(
            load,
            MmLoad::Tabulated(_)
                | MmLoad::Analytic(
                    Termination::UnflangedOpen
                        | Termination::UnflangedThinWall { .. }
                        | Termination::FlangedOpen
                )
        ) {
            return Err(EvanescentError::Invalid(
                "the source must independently match higher modes",
            ));
        }
        let higher = MatchedEvanescentPorts::new(
            gas,
            *radius,
            config.n_modes,
            1.0 / f64::from(config.sample_rate_hz),
            approximation,
        )?;
        let plane_config = MmLineConfig {
            n_modes: 1,
            ..*config
        };
        let mut bank = MmLineBank::new(
            core::slice::from_ref(duct),
            &["uniform cylinder"],
            gas,
            load,
            &plane_config,
        )
        .map_err(EvanescentError::Plane)?;
        Ok(Self {
            plane: bank
                .lines
                .pop()
                .expect("one admitted duct produces one line"),
            plane_receipt: bank
                .receipts
                .pop()
                .expect("one admitted duct produces one receipt"),
            higher,
            incoming: vec![0.0; config.n_modes],
        })
    }

    /// The original plane-mode FIR admission, including its physical delay.
    #[must_use]
    pub const fn plane_receipt(&self) -> &MmComboReceipt {
        &self.plane_receipt
    }
    /// Higher-mode source, cells, approximation bounds, storage and work.
    #[must_use]
    pub const fn higher_ports(&self) -> &MatchedEvanescentPorts {
        &self.higher
    }
    /// Last accepted vector, including plane mode at index zero.
    #[must_use]
    pub fn incoming(&self) -> &[f64] {
        &self.incoming
    }
    /// Current-sample affine incoming response for original radial `mode`.
    pub fn boundary_response(&self, mode: usize) -> Result<ScatteringResponse, EvanescentError> {
        if mode == 0 {
            self.plane
                .boundary_response()
                .map_err(EvanescentError::Discrete)
        } else {
            self.higher.boundary_response(mode - 1)
        }
    }

    /// Step the original plane owner and every higher port once, on one clock.
    /// No extra latency and no allocation. Refusals preserve all accepted states.
    pub fn push(&mut self, outgoing: &[f64]) -> Result<&[f64], EvanescentError> {
        if outgoing.len() != self.incoming.len() || outgoing.iter().any(|value| !value.is_finite())
        {
            return Err(EvanescentError::Invalid(
                "one finite outgoing wave per retained mode is required",
            ));
        }
        self.higher.prepare(&outgoing[1..])?;
        // The plain FIR owner validates its entire output before its sole commit.
        // Higher states are still only candidates if that validation refuses.
        let plane = self
            .plane
            .push(outgoing[0])
            .map_err(EvanescentError::Discrete)?;
        self.higher.accept();
        self.incoming[0] = plane;
        self.incoming[1..].copy_from_slice(self.higher.incoming());
        Ok(&self.incoming)
    }
}

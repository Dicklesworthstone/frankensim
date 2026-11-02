//! Homogeneous Fourier synthesis from the complete Mann spectral tensor.
//!
//! Each mode has two solenoidal velocity polarizations `B`, with
//! `B B^T = w Phi(k)`. Independent unit-variance Gaussian cosine/sine
//! quadratures give covariance `w Phi(k) cos(k.r)` at every location, not
//! only a fitted one-point stress. Log-radial strata and randomly shifted
//! low-discrepancy directions approximate the explicitly retained finite
//! wavenumber band. This is a finite quadrature of the homogeneous model;
//! neither the omitted bands nor quadrature error are silently renormalized.
//!
//! Four exact-discrete OU paths evolve the two polarizations and quadratures.
//! Their correlation times are caller inputs, separate from the Mann eddy
//! lifetime used in the spatial tensor. Advection is explicit. Sampling
//! evaluates precomputed polarizations without allocating or drawing noise.
//!
//! **There is no wall boundary condition.** Spatial coordinates can have
//! either sign, including negative z; vertical velocity generally does not
//! vanish at z = 0. Do not substitute this field for a wall-compatible field
//! near terrain without a separately justified boundary treatment.
//!
//! Spatial synthesis follows the spectral-tensor construction of Mann,
//! *Wind field simulation*, Probabilistic Engineering Mechanics 13 (1998),
//! 269–282, <https://doi.org/10.1016/S0266-8920(97)00036-2>.

use crate::mann::{MIN_KL, MannParams, mann_tensor};
use crate::ou::{MAX_ADVANCE_TICKS, OuCheckpoint, OuMode, StationaryOuPath};
use crate::{FlowSample, MAX_MODES, Refusal, TICK_HZ, refuse};
use fs_math::det;
use fs_rand::StreamKey;

/// Exact spatial observation filters on the represented Fourier field.
pub mod filter;
use filter::MannBoxFilter;

/// Independent Philox namespace for homogeneous Mann quadrature geometry.
pub const MANN_FIELD_KERNEL: u32 = 0x4d41_4e4e;

const PATH_KEYS: [u64; 4] = [
    0x6d61_6e6e_6331_5f31,
    0x6d61_6e6e_7331_5f31,
    0x6d61_6e6e_6332_5f31,
    0x6d61_6e6e_7332_5f31,
];

/// Explicit finite-band quadrature. Wavenumber is in radians per metre;
/// the bounds are dimensionless `|k| L`, using the target's length scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannQuadrature {
    /// Number of wavevectors (each carries two polarizations), 1..=256.
    pub n_modes: usize,
    /// Lower resolved `|k| L`, at least [`MIN_KL`].
    pub min_kl: f64,
    /// Upper resolved `|k| L`, finite and greater than `min_kl`.
    pub max_kl: f64,
}

#[derive(Clone, Debug, PartialEq)]
struct SpectralMode {
    k: [f64; 3],
    weight: f64,
    polarization: [[f64; 3]; 2],
}

/// Immutable, seed-specific quadrature of the full homogeneous Mann tensor.
/// Geometry and amplitude-path seeds are separate, so independent velocity
/// realizations can share precisely the same represented spatial statistics.
#[derive(Clone, Debug, PartialEq)]
pub struct MannSpectralModel {
    geometry_seed: u64,
    parameters: MannParams,
    quadrature: MannQuadrature,
    advection_mps: [f64; 3],
    modes: Vec<SpectralMode>,
}

fn invalid(message: &str) -> Refusal {
    refuse(
        "mann-field-invalid",
        message.into(),
        "use a finite admitted target, 1..=256 modes and explicit representable kL bounds",
    )
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn radical_inverse(mut n: usize, base: usize) -> f64 {
    let mut value = 0.0;
    let mut factor = 1.0 / base as f64;
    while n != 0 {
        value += (n % base) as f64 * factor;
        n /= base;
        factor /= base as f64;
    }
    value
}

fn shifted_unit(value: f64, shift: f64) -> f64 {
    let sum = value + shift;
    if sum >= 1.0 { sum - 1.0 } else { sum }
}

// Factor in k's transverse plane instead of attempting a singular 3x3
// Cholesky. A pivoted 2x2 factor retains cross-correlations and enforces
// incompressibility structurally in both independent polarizations.
fn polarizations(k: [f64; 3], phi: [[f64; 3]; 3], weight: f64) -> Result<[[f64; 3]; 2], Refusal> {
    let norm = det::sqrt(dot(k, k));
    let n = k.map(|v| v / norm);
    let axis = (0..3)
        .min_by(|&i, &j| n[i].abs().total_cmp(&n[j].abs()))
        .expect("three axes");
    let mut reference = [0.0; 3];
    reference[axis] = 1.0;
    let transverse = cross(n, reference);
    let transverse_norm = det::sqrt(dot(transverse, transverse));
    let mut e0 = transverse.map(|v| v / transverse_norm);
    let mut e1 = cross(n, e0);
    let scale = phi[0][0].max(phi[1][1]).max(phi[2][2]);
    if !scale.is_finite() || scale <= 0.0 || !phi.iter().flatten().all(|v| v.is_finite()) {
        return Err(invalid("spectral tensor is not representable"));
    }
    let project = |a: [f64; 3], b: [f64; 3]| -> f64 {
        (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| a[i] * (phi[i][j] / scale) * b[j])
                    .sum::<f64>()
            })
            .sum()
    };
    let mut c00 = project(e0, e0);
    let mut c11 = project(e1, e1);
    let c01 = project(e0, e1);
    if c11 > c00 {
        core::mem::swap(&mut c00, &mut c11);
        core::mem::swap(&mut e0, &mut e1);
    }
    if !c00.is_finite() || c00 <= 0.0 || !c11.is_finite() || !c01.is_finite() {
        return Err(invalid(
            "spectral tensor has no positive transverse variance",
        ));
    }
    let l00 = det::sqrt(c00);
    let l10 = c01 / l00;
    let remainder = c11 - l10 * l10;
    if remainder < -1e-10 * c00 {
        return Err(invalid("spectral tensor is not positive semidefinite"));
    }
    let l11 = det::sqrt(remainder.max(0.0));
    let amplitude = det::sqrt(scale) * det::sqrt(weight);
    let result = [
        core::array::from_fn(|i| amplitude * (e0[i] * l00 + e1[i] * l10)),
        core::array::from_fn(|i| amplitude * e1[i] * l11),
    ];
    if !result.iter().flatten().all(|v| v.is_finite()) {
        return Err(invalid(
            "weighted spectral polarization is not representable",
        ));
    }
    Ok(result)
}

impl MannSpectralModel {
    /// Build a bounded finite-band quadrature of `Phi_ij(k) d^3k`.
    /// Log-radius is stratified; two independently shifted low-discrepancy
    /// sequences supply uniform solid angle. The full-sphere weight is
    /// `4 pi k^3 delta(log k) / N`. A real cosine/sine pair represents the
    /// Hermitian +/-k pair, so no additional factor of two is applied.
    ///
    /// # Errors
    /// Mann-target admission errors or `mann-field-invalid` for invalid
    /// quadrature, advection, or unrepresentable spectral coefficients.
    pub fn build(
        geometry_seed: u64,
        parameters: MannParams,
        quadrature: MannQuadrature,
        advection_mps: [f64; 3],
    ) -> Result<Self, Refusal> {
        parameters.admit()?;
        if !(1..=MAX_MODES).contains(&quadrature.n_modes)
            || !quadrature.min_kl.is_finite()
            || quadrature.min_kl < MIN_KL
            || !quadrature.max_kl.is_finite()
            || quadrature.max_kl <= quadrature.min_kl
            || !advection_mps.iter().all(|v| v.is_finite())
        {
            return Err(invalid("invalid mode count, kL band or advection"));
        }
        let mut stream = StreamKey {
            seed: geometry_seed,
            kernel: MANN_FIELD_KERNEL,
            tile: 0,
        }
        .stream();
        let shifts = [stream.next_f64(), stream.next_f64(), stream.next_f64()];
        let log_min = det::ln(quadrature.min_kl);
        let log_span = det::ln(quadrature.max_kl) - log_min;
        let mut modes = Vec::with_capacity(quadrature.n_modes);
        for index in 0..quadrature.n_modes {
            let radial = (index as f64 + shifts[0]) / quadrature.n_modes as f64;
            let radius = det::exp(log_min + radial * log_span) / parameters.length_m;
            let mu = 2.0 * shifted_unit(radical_inverse(index + 1, 2), shifts[1]) - 1.0;
            let azimuth =
                core::f64::consts::TAU * shifted_unit(radical_inverse(index + 1, 3), shifts[2]);
            let planar = det::sqrt((1.0 - mu * mu).max(0.0));
            let k = [
                radius * planar * det::cos(azimuth),
                radius * planar * det::sin(azimuth),
                radius * mu,
            ];
            let weight = 4.0 * core::f64::consts::PI * log_span * radius * radius * radius
                / quadrature.n_modes as f64;
            if !weight.is_finite() || weight <= 0.0 || !k.iter().all(|v| v.is_finite()) {
                return Err(invalid("wavenumber quadrature is not representable"));
            }
            let phi = mann_tensor(&parameters, k)?;
            modes.push(SpectralMode {
                k,
                weight,
                polarization: polarizations(k, phi, weight)?,
            });
        }
        Ok(Self {
            geometry_seed,
            parameters,
            quadrature,
            advection_mps,
            modes,
        })
    }

    /// Number of wavevectors; each has two independent polarizations.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.modes.len()
    }

    /// Explicit resolved band and quadrature size. No variance correction
    /// compensates for unresolved low or high wavenumbers.
    #[must_use]
    pub fn quadrature(&self) -> MannQuadrature {
        self.quadrature
    }

    /// Spatial target used for every mode, including all cross-components.
    #[must_use]
    pub fn parameters(&self) -> MannParams {
        self.parameters
    }

    /// Exact represented equal-time covariance `E[u_i(x) u_j(x+r)]`.
    /// At zero separation this gives the retained Reynolds-stress tensor;
    /// at nonzero separation it reports the actual finite-mode coherence.
    /// This is not an error estimate against the continuum Mann integral.
    pub fn covariance(&self, separation_m: [f64; 3]) -> Result<[[f64; 3]; 3], Refusal> {
        finite_vector(separation_m)?;
        let mut covariance = [[0.0; 3]; 3];
        for mode in &self.modes {
            let phase = dot(mode.k, separation_m);
            if !phase.is_finite() {
                return Err(invalid("covariance separation is not representable"));
            }
            accumulate_covariance(mode, det::cos(phase), &mut covariance);
        }
        finite_covariance(covariance)
    }
}

fn finite_vector(value: [f64; 3]) -> Result<(), Refusal> {
    if !value.iter().all(|v| v.is_finite()) {
        return Err(refuse(
            "non-finite-input",
            format!("{value:?}"),
            "supply finite coordinates",
        ));
    }
    Ok(())
}

fn accumulate_covariance(mode: &SpectralMode, multiplier: f64, covariance: &mut [[f64; 3]; 3]) {
    for i in 0..3 {
        for j in 0..3 {
            covariance[i][j] += multiplier
                * (mode.polarization[0][i] * mode.polarization[0][j]
                    + mode.polarization[1][i] * mode.polarization[1][j]);
        }
    }
}

fn finite_covariance(covariance: [[f64; 3]; 3]) -> Result<[[f64; 3]; 3], Refusal> {
    if !covariance.iter().flatten().all(|v| v.is_finite()) {
        return Err(invalid("represented covariance is not finite"));
    }
    Ok(covariance)
}

/// Native prescribed homogeneous turbulence with explicit sequential time.
/// There is no mean wind in the returned velocity and no wall at z = 0.
#[derive(Clone, Debug, PartialEq)]
pub struct HomogeneousMannField {
    model: MannSpectralModel,
    evolution_seed: u64,
    correlation_times_s: Vec<f64>,
    paths: [StationaryOuPath; 4],
}

/// Complete model inputs and four OU states. Resume deterministically
/// rebuilds the quadrature and tensor factors; no external model is needed.
#[derive(Clone, Debug, PartialEq)]
pub struct MannFieldCheckpoint {
    /// Seed defining the quadrature geometry.
    pub geometry_seed: u64,
    /// Full spatial spectral target.
    pub parameters: MannParams,
    /// Explicit resolved wavenumber band and mode count.
    pub quadrature: MannQuadrature,
    /// Uniform phase-convection velocity [m/s].
    pub advection_mps: [f64; 3],
    /// Independent amplitude-path seed.
    pub evolution_seed: u64,
    /// One temporal decorrelation time [s] per wavevector.
    pub correlation_times_s: Vec<f64>,
    /// Cosine/sine states of polarization 1, then polarization 2.
    pub paths: [OuCheckpoint; 4],
}

fn ou_modes(model: &MannSpectralModel, times: &[f64]) -> Result<Vec<OuMode>, Refusal> {
    if times.len() != model.mode_count() {
        return Err(invalid("one correlation time is required per wavevector"));
    }
    times
        .iter()
        .map(|&tau| OuMode::from_correlation_time(1.0, tau))
        .collect()
}

impl HomogeneousMannField {
    /// Draw a stationary field at the shared canonical OU anchor. Each
    /// polarization has independent unit-variance cosine and sine paths;
    /// using both retains spatial homogeneity and the complete modal tensor.
    pub fn at_anchor(
        model: MannSpectralModel,
        evolution_seed: u64,
        correlation_times_s: Vec<f64>,
    ) -> Result<Self, Refusal> {
        let modes = ou_modes(&model, &correlation_times_s)?;
        let make = |index: usize| {
            StationaryOuPath::stationary_at_anchor(evolution_seed ^ PATH_KEYS[index], modes.clone())
        };
        let paths = [make(0)?, make(1)?, make(2)?, make(3)?];
        Ok(Self {
            model,
            evolution_seed,
            correlation_times_s,
            paths,
        })
    }

    /// Immutable represented spatial model, available for covariance queries.
    #[must_use]
    pub fn model(&self) -> &MannSpectralModel {
        &self.model
    }

    /// Current signed physical tick, including canonical prehistory.
    #[must_use]
    pub fn tick(&self) -> i64 {
        self.paths[0].tick()
    }

    /// Copy only the evolving state from an identical model, reusing all
    /// storage. Useful for accepted/trial fields and atomic consumer steps;
    /// spectral coefficients are neither cloned nor recomputed.
    ///
    /// # Errors
    /// `mann-field-invalid` if geometry, target, cutoffs, advection, amplitude
    /// seed or any OU correlation time differs. Refusal leaves `self` intact.
    pub fn copy_temporal_state_from(&mut self, source: &Self) -> Result<(), Refusal> {
        if self.model != source.model
            || self.evolution_seed != source.evolution_seed
            || self.correlation_times_s != source.correlation_times_s
        {
            return Err(invalid(
                "temporal state belongs to a different spatial or OU model",
            ));
        }
        for (target, source) in self.paths.iter_mut().zip(&source.paths) {
            target.copy_temporal_state_from(source);
        }
        Ok(())
    }

    /// Advance all quadratures together. Backward or excessive requests
    /// refuse before changing any path; samples never advance the clock.
    pub fn advance_to(&mut self, target_tick: i64) -> Result<(), Refusal> {
        if target_tick < self.tick() {
            return Err(refuse(
                "ou-advance-backwards",
                "target precedes the field clock".into(),
                "resume an earlier checkpoint",
            ));
        }
        if target_tick
            .checked_sub(self.tick())
            .is_none_or(|span| span > MAX_ADVANCE_TICKS)
        {
            return Err(refuse(
                "ou-advance-span-exceeded",
                "advance exceeds the bounded span".into(),
                "advance in bounded chunks",
            ));
        }
        for path in &mut self.paths {
            path.advance_to(target_tick)?;
        }
        Ok(())
    }

    /// Velocity and exact analytic gradient at the current physical time.
    /// Every polarization is perpendicular to its wavevector, so the
    /// instantaneous gradient trace vanishes to floating-point roundoff.
    /// Coordinates can have either sign; z = 0 has no special meaning.
    /// This operation draws no randomness and allocates no storage.
    pub fn sample(&self, point_m: [f64; 3]) -> Result<FlowSample, Refusal> {
        self.sample_with_gains(point_m, None)
    }

    /// Exact average of velocity and its analytic gradient over a prepared
    /// box, face or line centered at `point_m`. The same scalar gain acts on
    /// both transverse polarizations, preserving solenoidality. This filters
    /// velocity, not nonlinear momentum or energy fluxes.
    pub fn sample_filtered(
        &self,
        point_m: [f64; 3],
        filter: &MannBoxFilter,
    ) -> Result<FlowSample, Refusal> {
        self.sample_with_gains(point_m, filter.gains(&self.model)?)
    }

    fn sample_with_gains(
        &self,
        point_m: [f64; 3],
        gains: Option<&[f64]>,
    ) -> Result<FlowSample, Refusal> {
        finite_vector(point_m)?;
        let time = self.tick() as f64 / TICK_HZ;
        let position = core::array::from_fn(|i| point_m[i] - self.model.advection_mps[i] * time);
        let mut sample = FlowSample {
            u: [0.0; 3],
            grad: [[0.0; 3]; 3],
        };
        for (index, mode) in self.model.modes.iter().enumerate() {
            let phase = dot(mode.k, position);
            if !phase.is_finite() {
                return Err(invalid("sample phase is not representable"));
            }
            let (s, c) = (det::sin(phase), det::cos(phase));
            let a0 = self.paths[0].amplitudes()[index];
            let b0 = self.paths[1].amplitudes()[index];
            let a1 = self.paths[2].amplitudes()[index];
            let b1 = self.paths[3].amplitudes()[index];
            let mut value = [a0 * c + b0 * s, a1 * c + b1 * s];
            let mut derivative = [-a0 * s + b0 * c, -a1 * s + b1 * c];
            if let Some(gains) = gains {
                for v in value.iter_mut().chain(&mut derivative) {
                    *v *= gains[index];
                }
            }
            for i in 0..3 {
                sample.u[i] +=
                    mode.polarization[0][i] * value[0] + mode.polarization[1][i] * value[1];
                let gradient = mode.polarization[0][i] * derivative[0]
                    + mode.polarization[1][i] * derivative[1];
                for j in 0..3 {
                    sample.grad[i][j] += gradient * mode.k[j];
                }
            }
        }
        if !sample
            .u
            .iter()
            .chain(sample.grad.iter().flatten())
            .all(|v| v.is_finite())
        {
            return Err(invalid(
                "sample coordinates or velocity are not representable",
            ));
        }
        Ok(sample)
    }

    /// Exact represented covariance `E[u_i(x,t) u_j(x+r,t+lag)]`.
    /// Every modal tensor is multiplied by its advected phase and OU decay.
    /// This exposes the retained space-time model without sampling noise.
    pub fn covariance(&self, separation_m: [f64; 3], lag_s: f64) -> Result<[[f64; 3]; 3], Refusal> {
        self.covariance_with_gains(separation_m, lag_s, None)
    }

    /// Exact represented covariance when the same spatial filter is applied
    /// at both observations. Each modal tensor receives the SQUARED gain;
    /// attenuation is retained, never renormalized to the pointwise variance.
    pub fn covariance_filtered(
        &self,
        separation_m: [f64; 3],
        lag_s: f64,
        filter: &MannBoxFilter,
    ) -> Result<[[f64; 3]; 3], Refusal> {
        self.covariance_with_gains(separation_m, lag_s, filter.gains(&self.model)?)
    }

    fn covariance_with_gains(
        &self,
        separation_m: [f64; 3],
        lag_s: f64,
        gains: Option<&[f64]>,
    ) -> Result<[[f64; 3]; 3], Refusal> {
        finite_vector(separation_m)?;
        if !lag_s.is_finite() {
            return Err(invalid("covariance lag must be finite"));
        }
        let displacement =
            core::array::from_fn(|i| separation_m[i] - self.model.advection_mps[i] * lag_s);
        let mut covariance = [[0.0; 3]; 3];
        for (index, (mode, tau)) in self
            .model
            .modes
            .iter()
            .zip(&self.correlation_times_s)
            .enumerate()
        {
            let phase = dot(mode.k, displacement);
            if !phase.is_finite() {
                return Err(invalid("space-time covariance phase is not representable"));
            }
            let mut multiplier = det::cos(phase) * det::exp(-lag_s.abs() / tau);
            if let Some(gains) = gains {
                multiplier *= gains[index] * gains[index];
            }
            accumulate_covariance(mode, multiplier, &mut covariance);
        }
        finite_covariance(covariance)
    }

    /// Retain every model input and quadrature path for exact restoration.
    #[must_use]
    pub fn checkpoint(&self) -> MannFieldCheckpoint {
        MannFieldCheckpoint {
            geometry_seed: self.model.geometry_seed,
            parameters: self.model.parameters,
            quadrature: self.model.quadrature,
            advection_mps: self.model.advection_mps,
            evolution_seed: self.evolution_seed,
            correlation_times_s: self.correlation_times_s.clone(),
            paths: core::array::from_fn(|i| self.paths[i].checkpoint()),
        }
    }

    /// Rebuild the admitted finite-band model and restore all four paths.
    /// A mismatched clock/key, malformed model or nonfinite state refuses.
    pub fn resume(checkpoint: MannFieldCheckpoint) -> Result<Self, Refusal> {
        let model = MannSpectralModel::build(
            checkpoint.geometry_seed,
            checkpoint.parameters,
            checkpoint.quadrature,
            checkpoint.advection_mps,
        )?;
        let modes = ou_modes(&model, &checkpoint.correlation_times_s)?;
        let tick = checkpoint.paths[0].tick;
        if checkpoint
            .paths
            .iter()
            .enumerate()
            .any(|(i, p)| p.tick != tick || p.seed != checkpoint.evolution_seed ^ PATH_KEYS[i])
        {
            return Err(invalid(
                "checkpoint quadratures have inconsistent clocks or stream keys",
            ));
        }
        let [c0, s0, c1, s1] = checkpoint.paths;
        let paths = [
            StationaryOuPath::resume(c0, modes.clone())?,
            StationaryOuPath::resume(s0, modes.clone())?,
            StationaryOuPath::resume(c1, modes.clone())?,
            StationaryOuPath::resume(s1, modes)?,
        ];
        Ok(Self {
            model,
            evolution_seed: checkpoint.evolution_seed,
            correlation_times_s: checkpoint.correlation_times_s,
            paths,
        })
    }
}

#[cfg(test)]
mod tests;

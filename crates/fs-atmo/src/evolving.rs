//! Advected, decorrelating turbulence on the wall-compatible spatial basis.
//!
//! Each vector-potential mode uses independent stationary OU quadratures:
//! `c(t) cos(theta) + s(t) sin(theta)`, with `Var(c) = Var(s) = 1/2` and
//! `theta = kx (x - U_adv t) + ky y + phi`. A single random multiplier of
//! a frozen mode would give position-dependent variance; the quadrature pair
//! preserves horizontal homogeneity and the frozen mode's phase-averaged
//! covariance, including its cross-component stress. Modal autocovariances
//! decay as `exp(-abs(dt)/tau)` in addition to mean advection.
//!
//! Time evolution acts on the potential before taking its curl. Consequently
//! every instantaneous field remains divergence-free and has zero normal
//! velocity at the wall. Sampling never advances time or draws random numbers.
//! This is a finite-mode Gaussian OU construction, not a turbulence PDE or a
//! validation of the basis's spectrum against atmospheric measurements. See
//! the Fourier OU construction in *A spatio-temporal random synthetic turbulent
//! velocity field: The underlying Gaussian structure*, JFM (2026),
//! <https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/5686250725BCEC52F0D1D8E8E48C68A5>.

use crate::ou::{MAX_ADVANCE_TICKS, OuCheckpoint, OuMode, StationaryOuPath};
use crate::{
    AirScenario, AirState, Atmosphere, FlatSiteLogLaw, FlowSample, MAX_MODES, Refusal, TICK_HZ,
    TurbulenceField, refuse,
};
use fs_math::det;

// Separate Philox keys without changing existing OU stream semantics.
const SINE_SEED_XOR: u64 = 0x7369_6e65_5f4f_555f;

/// Wall-compatible spatial field with explicit, checkpointed temporal state.
/// The spatial coefficients are fixed; each modal quadrature has the caller's
/// declared correlation time in seconds. No time is inferred from the frozen
/// builder's unvalidated spectral envelope.
#[derive(Clone, Debug, PartialEq)]
pub struct EvolvingTurbulenceField {
    basis: TurbulenceField,
    correlation_times_s: Vec<f64>,
    cosine: StationaryOuPath,
    sine: StationaryOuPath,
}

/// Self-contained field state, including the spatial basis and temporal model.
/// Resuming requires no external correlation times or advection speed.
#[derive(Clone, Debug, PartialEq)]
pub struct EvolvingTurbulenceCheckpoint {
    /// Wall-compatible spatial modes and their advection speed.
    pub basis: TurbulenceField,
    /// One correlation time [s] per spatial mode.
    pub correlation_times_s: Vec<f64>,
    /// Cosine-quadrature state.
    pub cosine: OuCheckpoint,
    /// Sine-quadrature state, at the same tick as `cosine`.
    pub sine: OuCheckpoint,
}

fn admit_model(basis: &TurbulenceField, times: &[f64]) -> Result<Vec<OuMode>, Refusal> {
    if basis.modes.is_empty()
        || basis.modes.len() > MAX_MODES
        || times.len() != basis.modes.len()
        || !basis.u_adv_mps.is_finite()
        || basis.modes.iter().any(|m| {
            ![m.kx, m.ky, m.kh, m.phi, m.a[0], m.a[1], m.a[2]]
                .iter()
                .all(|v| v.is_finite())
        })
    {
        return Err(refuse(
            "evolving-params-invalid",
            format!(
                "{} finite spatial modes require the same number of correlation times, got {}",
                basis.modes.len(),
                times.len()
            ),
            "supply one positive finite correlation time per admitted spatial mode",
        ));
    }
    times
        .iter()
        .map(|&tau| OuMode::from_correlation_time(core::f64::consts::FRAC_1_SQRT_2, tau))
        .collect()
}

fn admit_point(point: [f64; 3]) -> Result<(), Refusal> {
    if !point.iter().all(|v| v.is_finite()) {
        return Err(refuse(
            "non-finite-input",
            format!("{point:?}"),
            "sample finite positions",
        ));
    }
    if point[2] < 0.0 {
        return Err(refuse(
            "below-surface-query",
            format!(
                "height {} m is below the aerodynamic ground plane",
                point[2]
            ),
            "sample at or above h = 0",
        ));
    }
    Ok(())
}

impl EvolvingTurbulenceField {
    /// Start both quadratures in their stationary distributions at
    /// [`crate::ou::STATIONARY_ANCHOR_TICK`]. Build the spatial basis with
    /// [`TurbulenceField::build`] or [`crate::fit::build_fitted_field`].
    ///
    /// # Errors
    /// `evolving-params-invalid` for a malformed basis or mode-count mismatch;
    /// `ou-params-invalid` for a non-positive or non-finite correlation time.
    pub fn at_anchor(
        basis: TurbulenceField,
        seed: u64,
        correlation_times_s: Vec<f64>,
    ) -> Result<Self, Refusal> {
        let modes = admit_model(&basis, &correlation_times_s)?;
        Ok(Self {
            basis,
            correlation_times_s,
            cosine: StationaryOuPath::stationary_at_anchor(seed, modes.clone())?,
            sine: StationaryOuPath::stationary_at_anchor(seed ^ SINE_SEED_XOR, modes)?,
        })
    }

    /// Physical tick of every sample, including negative pre-roll ticks.
    #[must_use]
    pub fn tick(&self) -> i64 {
        self.cosine.tick()
    }

    /// Advance both quadratures. This is the only operation that changes
    /// physical time. Invalid requests leave both paths intact.
    ///
    /// # Errors
    /// `ou-advance-backwards` or `ou-advance-span-exceeded`, as for
    /// [`StationaryOuPath::advance_to`].
    pub fn advance_to(&mut self, target_tick: i64) -> Result<(), Refusal> {
        if target_tick < self.tick() {
            return Err(refuse(
                "ou-advance-backwards",
                format!("target {target_tick} behind current {}", self.tick()),
                "resume an earlier checkpoint",
            ));
        }
        if target_tick
            .checked_sub(self.tick())
            .is_none_or(|span| span > MAX_ADVANCE_TICKS)
        {
            return Err(refuse(
                "ou-advance-span-exceeded",
                format!(
                    "advance from {} to {target_tick} exceeds the span cap",
                    self.tick()
                ),
                "advance in bounded chunks",
            ));
        }
        // The immutable model and clocks match; preflight establishes every
        // possible advance refusal before either path moves.
        self.cosine.advance_to(target_tick)?;
        self.sine.advance_to(target_tick)?;
        Ok(())
    }

    fn sample_admitted(&self, point: [f64; 3]) -> FlowSample {
        let [x, y, h] = point;
        let t = self.tick() as f64 / TICK_HZ;
        let mut sample = FlowSample {
            u: [0.0; 3],
            grad: [[0.0; 3]; 3],
        };
        for ((m, &a), &b) in self
            .basis
            .modes
            .iter()
            .zip(self.cosine.amplitudes())
            .zip(self.sine.amplitudes())
        {
            let theta = m.kx * (x - self.basis.u_adv_mps * t) + m.ky * y + m.phi;
            let (s, c) = (det::sin(theta), det::cos(theta));
            // Preserve the derivative pair: modulate the potential, not
            // individual velocity components after projection.
            m.accumulate(h, a * s - b * c, a * c + b * s, &mut sample);
        }
        sample
    }

    /// Velocity and exact analytic gradient at the current tick. This query
    /// allocates nothing and consumes no randomness; observer count and query
    /// order cannot affect the flow.
    ///
    /// # Errors
    /// `non-finite-input` or `below-surface-query`.
    pub fn sample(&self, x_m: f64, y_m: f64, h_m: f64) -> Result<FlowSample, Refusal> {
        let point = [x_m, y_m, h_m];
        admit_point(point)?;
        Ok(self.sample_admitted(point))
    }

    /// Fill caller-owned storage at one tick, without allocations. All points
    /// are admitted before writing, so a refusal leaves output intact.
    ///
    /// # Errors
    /// `sample-count-mismatch`, `non-finite-input`, or `below-surface-query`.
    pub fn sample_batch(
        &self,
        points: &[[f64; 3]],
        output: &mut [FlowSample],
    ) -> Result<(), Refusal> {
        if points.len() != output.len() {
            return Err(refuse(
                "sample-count-mismatch",
                format!("{} points, {} output slots", points.len(), output.len()),
                "provide one output slot per point",
            ));
        }
        for &point in points {
            admit_point(point)?;
        }
        for (&point, sample) in points.iter().zip(output) {
            *sample = self.sample_admitted(point);
        }
        Ok(())
    }

    /// Snapshot the spatial model and both sequential paths.
    #[must_use]
    pub fn checkpoint(&self) -> EvolvingTurbulenceCheckpoint {
        EvolvingTurbulenceCheckpoint {
            basis: self.basis.clone(),
            correlation_times_s: self.correlation_times_s.clone(),
            cosine: self.cosine.checkpoint(),
            sine: self.sine.checkpoint(),
        }
    }

    /// Reconstruct a checkpoint without a pre-roll or new random draws.
    ///
    /// # Errors
    /// `evolving-checkpoint-invalid` if the quadrature clocks or stream keys
    /// disagree; model and OU-checkpoint admission errors pass through.
    pub fn resume(checkpoint: EvolvingTurbulenceCheckpoint) -> Result<Self, Refusal> {
        let modes = admit_model(&checkpoint.basis, &checkpoint.correlation_times_s)?;
        if checkpoint.cosine.tick != checkpoint.sine.tick
            || checkpoint.sine.seed != checkpoint.cosine.seed ^ SINE_SEED_XOR
        {
            return Err(refuse(
                "evolving-checkpoint-invalid",
                "modal quadratures need matched clocks and independent stream keys".into(),
                "restore both quadratures from the same checkpoint",
            ));
        }
        Ok(Self {
            basis: checkpoint.basis,
            correlation_times_s: checkpoint.correlation_times_s,
            cosine: StationaryOuPath::resume(checkpoint.cosine, modes.clone())?,
            sine: StationaryOuPath::resume(checkpoint.sine, modes)?,
        })
    }
}

/// Mean wind, evolving turbulence, and one shared air-property scenario.
#[derive(Clone, Debug, PartialEq)]
pub struct EvolvingAtmosphere {
    /// Scenario-level mean-wind law.
    pub mean: FlatSiteLogLaw,
    /// Explicit sequential turbulence state.
    pub turbulence: EvolvingTurbulenceField,
    /// Air properties used by dynamic pressure and Reynolds number.
    pub air: AirScenario,
}

impl Atmosphere {
    /// Enable modal decorrelation explicitly. Frozen atmosphere APIs and
    /// their deterministic realizations retain their existing behavior.
    ///
    /// # Errors
    /// Mean-law and [`EvolvingTurbulenceField::at_anchor`] admission errors.
    pub fn into_evolving(
        self,
        seed: u64,
        correlation_times_s: Vec<f64>,
    ) -> Result<EvolvingAtmosphere, Refusal> {
        self.mean.admit()?;
        Ok(EvolvingAtmosphere {
            mean: self.mean,
            turbulence: EvolvingTurbulenceField::at_anchor(
                self.turbulence,
                seed,
                correlation_times_s,
            )?,
            air: self.air,
        })
    }
}

impl EvolvingAtmosphere {
    /// Full air state at the turbulence state's current physical tick.
    ///
    /// # Errors
    /// Mean-law or spatial-query admission errors.
    pub fn sample_air_state(&self, x_m: f64, y_m: f64, h_m: f64) -> Result<AirState, Refusal> {
        let flow = self.sample_flow(x_m, y_m, h_m)?;
        Ok(AirState {
            velocity_mps: flow.u,
            grad: flow.grad,
            rho_kg_m3: self.air.rho_kg_m3,
            mu_kg_m_s: self.air.mu_kg_m_s,
            temperature_k: self.air.temperature_k,
            pressure_pa: self.air.pressure_pa,
            provenance: self.air.provenance,
        })
    }

    /// Sample the existing mean and current turbulence without imposing constant
    /// air properties. This is read-only and consumes no random samples; a
    /// higher-layer material owner may attach a position-dependent gas state.
    pub fn sample_flow(&self, x_m: f64, y_m: f64, h_m: f64) -> Result<FlowSample, Refusal> {
        self.mean.admit()?;
        let turbulent = self.turbulence.sample(x_m, y_m, h_m)?;
        let mut velocity = turbulent.u;
        velocity[0] += self.mean.speed(h_m);
        let mut grad = turbulent.grad;
        grad[0][2] += self.mean.dspeed_dh(h_m);
        Ok(FlowSample { u: velocity, grad })
    }
}

#[cfg(test)]
mod tests;

//! Force-driven plate images in the existing sample-accurate renderer.
//!
//! `thin_plate` owns geometry/material reduction, exact-ZOH mechanics and the
//! baffled compact-monopole observation. This module only hosts those bodies.
//! A force is the total force on the reduction's declared footprint; its signed
//! participation and radiation area are distinct conjugate/output projections.
//! Zero force releases the input, never the stored vibration. No reaction load,
//! frequency, normalization, contact law or material data is invented here.
//!
//! Pressure is the existing endpoint-acceleration observation, not an imaginary
//! narrow-band transfer applied to displacement (which would radiate a static
//! deflection). Compact radiation remains an approximation: no retardation,
//! atmospheric propagation, directivity or whole-fluid energy claim is added.

use super::RenderError;
use crate::acoustic_realize::AcousticRealizeError;
use crate::modal_acoustic_time::MAX_TIME_DOMAIN_ACOUSTIC_MODES;
use crate::thin_plate::{CompactBody, certified_radiators};
use fs_scenario::ThinPlate;

/// Explicit callback admission; these ceilings do not rescale sound or forces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateVoiceConfig {
    /// Fixed output/mechanical sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Maximum retained modes; checked before reducing a rectangular plate.
    pub max_modes: usize,
    /// Fraction of Nyquist allowed for every retained natural frequency.
    pub nyquist_guard_fraction: f64,
    /// Surrounding fluid density in the compact pressure law [kg/m^3].
    pub density_kg_m3: f64,
    /// Fixed observation distance [m].
    pub listener_m: f64,
    /// Maximum absolute total footprint force [N].
    pub maximum_abs_force_n: f64,
    /// Maximum absolute sum of compact modal pressures [Pa].
    pub maximum_abs_pressure_pa: f64,
}

impl PlateVoiceConfig {
    fn validate(self) -> Result<(), RenderError> {
        if self.sample_rate_hz == 0 || self.max_modes == 0
            || self.max_modes > MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || !self.nyquist_guard_fraction.is_finite()
            || !(0.0 < self.nyquist_guard_fraction && self.nyquist_guard_fraction < 1.0)
            || [self.density_kg_m3, self.listener_m, self.maximum_abs_force_n,
                self.maximum_abs_pressure_pa].iter().any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid("plate voice requires explicit positive finite clock, mode, fluid, force and pressure budgets"));
        }
        Ok(())
    }
}

/// One independent linear plate, driven through the reduction's force footprint.
/// Bodies and their state are retained across every callback/control boundary.
/// Construction does not attach a piston-reaction filter or silently choose one.
/// Bodies supplied with a filter retain it and its existing approximation scope.
pub struct CompactPlateVoice {
    bodies: Vec<CompactBody>,
    config: PlateVoiceConfig,
    force_n: f64,
    sample_index: u64,
    poisoned: bool,
}

impl CompactPlateVoice {
    /// Reduce actual rectangular geometry/materials through the existing DKT
    /// and modal owners, then admit the resulting runtime. The existing reducer
    /// owns its mesh/window choices. The nonlinear plate description refuses;
    /// it is not silently converted to a linear oscillator bank.
    ///
    /// # Errors
    /// Invalid budgets/force, mode count, nonlinear model or reducer refusal.
    pub fn from_plate(
        plate: ThinPlate,
        initial_force_n: f64,
        config: PlateVoiceConfig,
    ) -> Result<Self, RenderError> {
        config.validate()?;
        if !initial_force_n.is_finite() || initial_force_n.abs() > config.maximum_abs_force_n {
            return Err(invalid("plate footprint force exceeds its finite newton budget"));
        }
        if plate.geometric_nonlinearity || plate.n_modes == 0 || plate.n_modes > config.max_modes {
            return Err(invalid("plate voice requires a linear plate within its retained-mode budget"));
        }
        Self::from_radiators(certified_radiators(plate).map_err(RenderError::Voice)?, initial_force_n, config)
    }

    /// Host already reduced bodies, including those returned by
    /// `thin_plate::certified_chart_radiators` for a regional-material mesh and
    /// force footprint. Signed areas/participations and any existing vibration
    /// survive unchanged. The clock counts samples since this wrapper was made.
    /// A body list is not proof of source validity or modal completeness.
    ///
    /// # Errors
    /// Empty/oversized banks, invalid coefficients, Nyquist or force admission.
    pub fn from_radiators(
        bodies: Vec<CompactBody>,
        initial_force_n: f64,
        config: PlateVoiceConfig,
    ) -> Result<Self, RenderError> {
        config.validate()?;
        if bodies.is_empty() || bodies.len() > config.max_modes {
            return Err(invalid("plate body count must be nonzero and within max_modes"));
        }
        let maximum_omega = core::f64::consts::PI * f64::from(config.sample_rate_hz)
            * config.nyquist_guard_fraction;
        for body in &bodies {
            if !body.mass_kg.is_finite() || body.mass_kg <= 0.0
                || !body.omega.is_finite() || body.omega <= 0.0 || body.omega > maximum_omega
                || !body.zeta.is_finite() || body.zeta < 0.0
                || !body.area_m2.is_finite() || !body.drive_participation.is_finite()
                || !body.volume_velocity().is_finite()
            {
                return Err(invalid("plate bodies need finite physical coefficients and frequencies below the Nyquist guard"));
            }
        }
        let voice = Self { bodies, config, force_n: initial_force_n, sample_index: 0, poisoned: false };
        voice.validate_force(initial_force_n)?;
        Ok(voice)
    }

    /// Read the retained geometric reduction without mutating runtime state.
    #[must_use]
    pub fn bodies(&self) -> &[CompactBody] { &self.bodies }

    /// Number of completed samples since wrapping these bodies.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.sample_index }

    /// Exact mechanical period used by the shared body stepper [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { f64::from(self.config.sample_rate_hz).recip() }

    /// Validate a new total footprint force before any batch is applied.
    pub fn validate_force(&self, force_n: f64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if !force_n.is_finite() || force_n.abs() > self.config.maximum_abs_force_n
            || self.bodies.iter().any(|b| !(force_n * b.drive_participation).is_finite())
        {
            return Err(invalid("plate footprint force exceeds its finite newton budget"));
        }
        Ok(())
    }

    // RenderContext validates the whole batch before calling this mutation.
    pub(super) fn set_force_admitted(&mut self, force_n: f64) { self.force_n = force_n; }

    /// Advance existing body mechanics, sum existing compact pressures in body
    /// order, and retain every state. Refusals leave a partial callback and poison
    /// this voice. Clock/empty-input admission occurs before any state changes.
    /// The successful no-filter loop allocates no buffers; measured callback
    /// timing/allocation and loaded-filter performance remain separate evidence.
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if out.is_empty() { return Err(RenderError::EmptyBlock); }
        let count = u64::try_from(out.len()).map_err(|_| invalid("plate sample count overflow"))?;
        self.sample_index.checked_add(count).ok_or_else(|| invalid("plate sample clock overflow"))?;
        let dt = self.sample_period_s();
        for slot in out {
            let mut pressure = 0.0_f64;
            for body in &mut self.bodies {
                let force = self.force_n * body.drive_participation;
                let sample = match body.drive_and_radiate(force, dt, self.config.density_kg_m3, self.config.listener_m) {
                    Ok(sample) => sample,
                    Err(error) => {
                        self.poisoned = true;
                        return Err(RenderError::Voice(error));
                    }
                };
                pressure += sample;
            }
            if !pressure.is_finite() || pressure.abs() > self.config.maximum_abs_pressure_pa {
                self.poisoned = true;
                return Err(RenderError::Voice(AcousticRealizeError::InvalidDescription {
                    what: "compact plate pressure exceeds the declared finite output budget",
                }));
            }
            *slot = pressure;
            self.sample_index += 1;
        }
        Ok(())
    }
}

fn invalid(what: &'static str) -> RenderError { RenderError::Control { what } }

/// Bounded geometry/material/force-file input for the existing plate reduction.
pub mod file;
/// Nonlinear shell/film impact with physical striker and hysteretic felt storage.
pub mod impact;

//! Original oscillator mechanics on the native accepted acoustic-heat timeline.
//!
//! The material owner supplies exact spatial bulk/shear work and advances the
//! existing thermal history. This adapter retains the original WallPin and its
//! carrier response; it owns no thermal observer, integrator or heat account.
//! Its carrier is quasi-steady with frozen original material/support laws.
//! Envelope changes do not resolve acoustic startup or oscillator ring-up.

use super::{OscillatorBackedViscoelasticLayers, SupportedViscoelasticResponse};
use super::super::ElasticSupportError;
use fs_material::elastic::acoustic::viscoelastic::thermal::{
    AcousticBackingThermalCell, AcousticHeatingBudget, AcousticHeatingError,
    AcousticPowerSegment, AcousticThermalBin, SupportedViscoelasticHeating,
};
use fs_math::c64::C64;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit, StepConfig};

/// Original mechanical backing plus the original material's accepted heat owner.
///
/// The native source retains physical layer bins, caloric bindings, geometry,
/// mass, area, temperature bands and accepted time. A distinct backing receiver
/// is explicit; without it the resistance work leaves the thermal system.
/// Oscillator kinetic/spring storage is a separate harmonic observation.
pub struct OscillatorBackedViscoelasticHeating<'s, 'a> {
    mechanics: &'s OscillatorBackedViscoelasticLayers<'a>,
    heating: SupportedViscoelasticHeating<'s, 'a>,
}

impl<'a> OscillatorBackedViscoelasticLayers<'a> {
    /// Bind this original carrier solution to existing thermal depth cells.
    ///
    /// The peak pressure and area feed the same native supported coating used
    /// by `drive`. Positive cell work already includes the phasor factor 1/2.
    /// The native thermal owner chooses accepted durations and credits that
    /// work once. No pin storage or power-balance residual becomes heat.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_heating(
        &self,
        incident_peak_pressure_pa: C64,
        area_m2: f64,
        bins: &[AcousticThermalBin],
        thermal: &DepositionSimulation<'_>,
        budget: AcousticHeatingBudget,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<OscillatorBackedViscoelasticHeating<'_, 'a>, AcousticHeatingError> {
        let heating = SupportedViscoelasticHeating::try_new(
            self.layers(),
            incident_peak_pressure_pa,
            area_m2,
            bins,
            thermal,
            budget,
            checkpoint,
        )?;
        Ok(OscillatorBackedViscoelasticHeating {
            mechanics: self,
            heating,
        })
    }
}

impl<'s, 'a> OscillatorBackedViscoelasticHeating<'s, 'a> {
    /// Original mass, spring, resistance, carrier frequency and material layers.
    pub const fn mechanics(&self) -> &'s OscillatorBackedViscoelasticLayers<'a> {
        self.mechanics
    }

    /// Native spatial source and accepted accounting history; no copied account.
    pub const fn heating(&self) -> &SupportedViscoelasticHeating<'s, 'a> {
        &self.heating
    }

    /// Explicitly receive the resistance work in one distinct original cell.
    /// This forwards native geometry, caloric, temperature and no-double-credit
    /// admission. Reactive support storage remains outside the thermal sources.
    pub fn with_backing_receiver(
        mut self,
        thermal: &DepositionSimulation<'_>,
        receiver: AcousticBackingThermalCell,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, AcousticHeatingError> {
        self.heating = self.heating.with_backing_receiver(thermal, receiver, checkpoint)?;
        Ok(self)
    }

    /// Bind the native quasi-steady power envelope before its first advance.
    /// The material owner keeps jump times, affine exposure and source budgets.
    pub fn with_power_envelope(
        mut self,
        thermal: &DepositionSimulation<'_>,
        segments: &[AcousticPowerSegment],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, AcousticHeatingError> {
        self.heating = self.heating.with_power_envelope(thermal, segments, checkpoint)?;
        Ok(self)
    }

    /// Original harmonic response at the AUTHORED PEAK carrier pressure.
    ///
    /// This is independent of the last envelope multiplier and is not a causal
    /// support state or an accepted thermal-storage increment. It retains the
    /// actual rear velocity, displacement and original pin kinetic/spring
    /// storage under the coating's exp(+i omega t) convention.
    pub fn peak_response(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SupportedViscoelasticResponse<'_, 'a>, ElasticSupportError> {
        self.mechanics.drive(
            self.heating.incident_peak_pressure_pa(),
            self.heating.area_m2(),
            cancelled,
        )
    }

    /// Advance the native thermal candidate and its acoustic account together.
    /// The original heat owner retains conduction, sources, events, budgets,
    /// rollback and exact accepted durations. Inputs remain usable on refusal.
    pub fn advance_to<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'t>, Self, EnergyAudit), AcousticHeatingError> {
        let (thermal, heating, audit) =
            self.heating.advance_to(thermal, target_s, config, checkpoint)?;
        Ok((thermal, Self { mechanics: self.mechanics, heating }, audit))
    }
}

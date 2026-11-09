//! Original supported-coating heat on the shared accepted thermal timeline.
//!
//! Only independently integrated intrinsic bulk/shear work heats the coating.
//! Backing resistance work either leaves this accounting boundary or heats one
//! explicitly selected distinct native cell. Backing spring and kinetic energy
//! are never thermal sources. Original acoustic coefficients stay frozen inside
//! the caller's declared temperature bands.

use super::{
    AcousticBackingThermalCell, AcousticHeating, AcousticHeatingBudget, AcousticHeatingError,
    AcousticThermalBin, BoundBacking, BoundCell, Result, source,
};
use crate::elastic::acoustic::viscoelastic::layers::{
    SupportedViscoelasticLayers, SupportedViscoelasticPower,
};
use fs_math::c64::C64;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Accepted-time heating from one physical incident fluid, original supported
/// viscoelastic layers and their actual passive backing. No second fluid or
/// transmitted acoustic heat channel is introduced.
pub type SupportedViscoelasticHeating<'s, 'a> =
    AcousticHeating<&'s SupportedViscoelasticLayers<'a>, C64, SupportedViscoelasticPower>;

impl<'s, 'a> AcousticHeating<&'s SupportedViscoelasticLayers<'a>, C64, SupportedViscoelasticPower> {
    /// Bind every positive-thickness coating layer to a complete physical depth
    /// partition. Original interval strain work supplies each bin's heat.
    ///
    /// Empty and all-zero-thickness coatings admit zero bins. Their backing
    /// work is still accounted and may be routed to an explicit receiver.
    /// Without a receiver, backing work is exported and does not heat a cell.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        stack: &'s SupportedViscoelasticLayers<'a>,
        incident_peak_pressure_pa: C64,
        area_m2: f64,
        bins: &[AcousticThermalBin],
        thermal: &DepositionSimulation<'_>,
        budget: AcousticHeatingBudget,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Self> {
        let prepared = source::prepare(
            source::NativeStack::Supported(stack),
            [incident_peak_pressure_pa, C64::ZERO],
            area_m2,
            bins,
            thermal,
            budget,
            &mut checkpoint,
        )?;
        let source::NativePower::Supported(peak_power) = prepared.power else {
            return Err(AcousticHeatingError::Invalid("acoustic source kind disagrees"));
        };
        Ok(Self::from_prepared(
            stack, incident_peak_pressure_pa, peak_power, area_m2, thermal, budget, prepared,
        ))
    }

    /// Route all original backing resistance work to one explicitly authored
    /// native thermal cell, distinct from every coating bin. The native cell's
    /// original material, mass, geometry and caloric binding remain authoritative.
    /// Bind once before the first advance; later source/material changes refuse.
    ///
    /// The temperature band admits the caller's frozen support approximation.
    /// It does not infer temperature-dependent impedance, source refresh,
    /// mechanical storage release, or heat transfer through an interface.
    pub fn with_backing_receiver(
        mut self,
        thermal: &DepositionSimulation<'_>,
        binding: AcousticBackingThermalCell,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<Self> {
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        self.no_pending()?;
        if self.energy.duration_s != 0.0 || self.backing_receiver.is_some() {
            return Err(AcousticHeatingError::Invalid(
                "bind the backing receiver once before the first advance",
            ));
        }
        if self.bins.iter().any(|bin| bin.power.binding.cell == binding.cell) {
            return Err(AcousticHeatingError::Invalid(
                "backing receiver must be distinct from every physical coating bin",
            ));
        }
        let cell = BoundCell::new(binding.cell, binding.temperature_band_k, thermal)?;
        checkpoint()?;
        self.backing_receiver = Some(BoundBacking { binding, cell });
        Ok(self)
    }

    /// Exact explicitly selected receiver; None means backing work is exported.
    pub const fn backing_receiver(&self) -> Option<AcousticBackingThermalCell> {
        match &self.backing_receiver {
            Some(receiver) => Some(receiver.binding),
            None => None,
        }
    }
}

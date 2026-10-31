//! Frozen-phase sound propagation from an accepted finite water inventory.
//!
//! The cell evolves its equilibrium phase masses through the original IF97
//! energy/volume inverse. This observation freezes those masses over an acoustic
//! perturbation and reuses the original single-phase sound speeds or Wood
//! compliance. It is an explicit acoustic closure, not the fully equilibrated
//! boiling sound speed or a bubble-resonance/phase-transfer model.
//!
//! The observer borrows one accepted cell state and preserves its physical
//! clock, material endpoints and inversion residuals. In particular it does
//! not replace the original source density with M/V to hide a source recovery
//! residual. Long wavelengths, common pressure/velocity, adiabatic individual
//! phases and negligible interphase heat/mass transfer remain caller conditions.

use super::WaterCellState;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics, FrozenHomogeneousMixture};

/// Explicit frozen-phase acoustic view of one accepted water-cell snapshot.
///
/// This implements the existing gas/liquid acoustic interface, so the native
/// duct propagation and other acoustic consumers use the original material
/// state. Construct a new view after accepting a mass or energy transfer.
#[derive(Clone, Copy, Debug)]
pub struct FrozenWaterCellAcoustics<'a> {
    source: &'a WaterCellState,
    acoustic: FluidAcousticState,
    mixture: Option<FrozenHomogeneousMixture>,
}

impl<'a> FrozenWaterCellAcoustics<'a> {
    /// Read retained IF97 source endpoints; no new pressure/temperature solve or
    /// guessed vapor EOS is performed by observation.
    pub fn new(source: &'a WaterCellState) -> Result<Self, MaterialError> {
        let thermodynamic = source.source_state();
        let (acoustic, mixture) = if let Some(single) = thermodynamic.single_phase_source() {
            (single.acoustic_state()?, None)
        } else if let Some(pair) = thermodynamic.saturation_pair() {
            let mixture = FrozenHomogeneousMixture::from_mass_fraction(
                pair.liquid().acoustic_state()?,
                pair.vapor().acoustic_state()?,
                thermodynamic.vapor_mass_fraction(),
            )?;
            (mixture.state(), Some(mixture))
        } else {
            return Err(MaterialError::State {
                what: "water acoustic observation requires the retained IF97 phase source".into(),
            });
        };
        Ok(Self {
            source,
            acoustic,
            mixture,
        })
    }

    /// Complete original accepted inventory, clock and caloric/source residuals.
    pub const fn source(&self) -> &'a WaterCellState {
        self.source
    }

    /// Original phase acoustic coefficients and actual volume fraction when
    /// the IF97 state is represented by its retained coexistence pair.
    pub const fn mixture(&self) -> Option<FrozenHomogeneousMixture> {
        self.mixture
    }

    /// Frozen-phase material response; this read does not advance the cell.
    pub const fn state(&self) -> FluidAcousticState {
        self.acoustic
    }
}

impl FluidAcoustics for FrozenWaterCellAcoustics<'_> {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

#[cfg(test)]
#[path = "acoustics/tests.rs"]
mod tests;

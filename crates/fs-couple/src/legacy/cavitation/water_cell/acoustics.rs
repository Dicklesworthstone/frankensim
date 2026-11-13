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
use fs_material::phase::water::acoustics::HomogeneousEquilibriumWater;
use fs_thermochem::water::if97::regions::If97FluidSaturationPair;

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
    density_representation_difference_kg_m3: f64,
}

impl<'a> FrozenWaterCellAcoustics<'a> {
    /// Read retained IF97 source endpoints; no new pressure/temperature solve or
    /// guessed vapor EOS is performed by observation.
    pub fn new(source: &'a WaterCellState) -> Result<Self, MaterialError> {
        let thermodynamic = source.source_state();
        let (acoustic, mixture) = if let Some(single) = thermodynamic.single_phase_source() {
            (single.acoustic_state()?, None)
        } else if let Some(pair) = thermodynamic.regional_saturation_pair() {
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
        let (acoustic, density_representation_difference_kg_m3) =
            retain_source_density(source, acoustic)?;
        Ok(Self {
            source,
            acoustic,
            mixture,
            density_representation_difference_kg_m3,
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

    /// Frozen-phase response with exactly the retained thermal source density;
    /// this read does not advance the cell or conceal its volume residual.
    pub const fn state(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Thermal 1/v minus the independently summed acoustic mixture density.
    /// This is representation roundoff, not a change to conserved water mass.
    pub const fn density_representation_difference_kg_m3(&self) -> f64 {
        self.density_representation_difference_kg_m3
    }
}

impl FluidAcoustics for FrozenWaterCellAcoustics<'_> {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

/// Explicit fast-equilibration acoustic view of one accepted water inventory.
///
/// Unlike the frozen-phase observer, this lets liquid and vapor exchange heat
/// and mass during an adiabatic perturbation. The original IF97 energy/volume
/// tangent supplies its compressibility. Both are distinct asymptotic closures;
/// no frequency interpolation or relaxation time is inferred. Use this view
/// only when equilibration is fast and the perturbation stays on its source
/// branch. Exact coexistence endpoints refuse the ambiguous two-sided tangent.
#[derive(Clone, Copy, Debug)]
pub struct EquilibriumWaterCellAcoustics<'a> {
    source: &'a WaterCellState,
    acoustic: FluidAcousticState,
    equilibrium: Option<HomogeneousEquilibriumWater>,
    regional_equilibrium: Option<HomogeneousEquilibriumWater<If97FluidSaturationPair>>,
    density_representation_difference_kg_m3: f64,
}

impl<'a> EquilibriumWaterCellAcoustics<'a> {
    pub fn new(source: &'a WaterCellState) -> Result<Self, MaterialError> {
        let thermodynamic = source.source_state();
        let (acoustic, equilibrium, regional_equilibrium) = if let Some(single) =
            thermodynamic.single_phase_source()
        {
            (single.acoustic_state()?, None, None)
        } else if let Some(pair) = thermodynamic.saturation_pair() {
            let equilibrium =
                HomogeneousEquilibriumWater::try_new(pair, thermodynamic.vapor_mass_fraction())?;
            (equilibrium.state(), Some(equilibrium), None)
        } else if let Some(pair) = thermodynamic.regional_saturation_pair() {
            let equilibrium = HomogeneousEquilibriumWater::try_from_regional(
                pair,
                thermodynamic.vapor_mass_fraction(),
            )?;
            (equilibrium.state(), None, Some(equilibrium))
        } else {
            return Err(MaterialError::State {
                what: "equilibrium acoustic observation requires the retained IF97 phase source"
                    .into(),
            });
        };
        let (acoustic, density_representation_difference_kg_m3) =
            retain_source_density(source, acoustic)?;
        Ok(Self {
            source,
            acoustic,
            equilibrium,
            regional_equilibrium,
            density_representation_difference_kg_m3,
        })
    }

    pub const fn source(&self) -> &'a WaterCellState {
        self.source
    }
    pub const fn state(&self) -> FluidAcousticState {
        self.acoustic
    }
    /// Original Gibbs-only equilibrium source, tangent and source-fit defects.
    /// Regional coexistence charts expose their original pair separately.
    pub const fn equilibrium(&self) -> Option<HomogeneousEquilibriumWater> {
        self.equilibrium
    }
    /// Original regional coexistence source, including retained Region 3 roots
    /// and the actual pressure/latent-energy residuals used by its tangent.
    pub const fn regional_equilibrium(
        &self,
    ) -> Option<HomogeneousEquilibriumWater<If97FluidSaturationPair>> {
        self.regional_equilibrium
    }
    /// Thermal 1/v minus the selected tangent's independently summed density.
    pub const fn density_representation_difference_kg_m3(&self) -> f64 {
        self.density_representation_difference_kg_m3
    }
}

impl FluidAcoustics for EquilibriumWaterCellAcoustics<'_> {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

// Preserve thermal storage density and the selected physical bulk modulus.
// The original mixture/tangent remains accessible for independent diagnostics.
fn retain_source_density(
    source: &WaterCellState,
    selected: FluidAcousticState,
) -> Result<(FluidAcousticState, f64), MaterialError> {
    let density = source.source_state().bulk_density_kg_m3();
    let difference = density - selected.density_kg_m3();
    let acoustic = if density.to_bits() == selected.density_kg_m3().to_bits() {
        selected
    } else {
        FluidAcousticState::from_adiabatic_bulk_modulus(
            density,
            selected.adiabatic_bulk_modulus_pa(),
        )?
    };
    Ok((acoustic, difference))
}

#[cfg(test)]
#[path = "acoustics/tests.rs"]
mod tests;

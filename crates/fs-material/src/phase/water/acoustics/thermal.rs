//! Acoustic response of an already solved regional water thermal state.
//!
//! The thermal owner retains the original pressure, caloric state, phase
//! inventory and IF97 roots. This adapter performs no temperature/pressure
//! inverse or material query. Coexistence requires an explicit acoustic
//! timescale closure; a thermal equilibrium state alone does not select one.
//!
//! The exported density is exactly the thermal owner's `1/v`. The selected
//! source modulus is unchanged. Independently equivalent phase-volume sums
//! can differ by rounding; the original response and density difference remain
//! inspectable. No physical mass, enthalpy, phase fraction or source is changed.
//! This is a local nondispersive response, not finite-rate phase transfer,
//! acoustic attenuation, bubble dynamics or a spatial wave solver.

use super::super::regional::RegionalWaterState;
use super::HomogeneousEquilibriumWater;
use crate::MaterialError;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics, FrozenHomogeneousMixture};
use fs_thermochem::water::if97::regions::{If97FluidSaturationPair, If97FluidState};

/// Caller-selected response of a liquid/vapor coexistence state. There is no
/// default: equilibrium thermal storage does not imply fast acoustic exchange.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoexistenceAcousticClosure {
    /// Fixed phase masses; each phase responds adiabatically and both share
    /// pressure and velocity. The original Wood-compliance owner is retained.
    FrozenPhaseMasses,
    /// Fast phase/thermal equilibrium and one velocity. The original
    /// energy-conserving saturation tangent and its source defects are retained.
    HomogeneousEquilibrium,
}

/// The actual response owner, before choosing the thermal density representation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RegionalWaterAcousticResponse {
    /// Original homogeneous IF97 state, also used at an exact frozen endpoint.
    SinglePhase(If97FluidState),
    /// Original frozen phase responses, mass quality and both volume fractions.
    FrozenPhaseMasses(FrozenHomogeneousMixture),
    /// Original saturation pair, quality, tangent and unmodified source defects.
    HomogeneousEquilibrium(HomogeneousEquilibriumWater<If97FluidSaturationPair>),
}

/// An immutable acoustic view retaining the complete accepted thermal state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionalWaterAcoustics {
    source: RegionalWaterState,
    closure: CoexistenceAcousticClosure,
    response: RegionalWaterAcousticResponse,
    acoustic: FluidAcousticState,
    density_representation_difference_kg_m3: f64,
}

impl RegionalWaterAcoustics {
    /// Attach sound to an existing thermal/boiling result without resolving it
    /// again. The constructor allocates no storage on successful evaluation.
    ///
    /// A homogeneous state always uses its own IF97 response, independently of
    /// the coexistence choice. At exact quality zero/one, frozen phase masses
    /// select the original remaining phase. Equilibrium acoustics refuses these
    /// endpoints: infinitesimal condensation/boiling and single-phase changes
    /// have different one-sided tangents. Interior qualities are never clipped
    /// or switched to another closure when their selected response refuses.
    ///
    /// The selected modulus is combined with the ORIGINAL thermal density.
    /// When both density representations coincide, the acoustic snapshot is
    /// returned unchanged, including the exact single-phase sound speed.
    /// Otherwise c and impedance are derived from that density and the same K;
    /// the signed rounding difference and original response remain available.
    pub fn try_new(
        source: RegionalWaterState,
        closure: CoexistenceAcousticClosure,
    ) -> Result<Self, MaterialError> {
        let response = if let Some(single) = source.single_phase_source() {
            RegionalWaterAcousticResponse::SinglePhase(single)
        } else {
            let pair = *source.law().saturation().ok_or_else(|| MaterialError::State {
                what: "thermal water coexistence requires its original saturation pair".into(),
            })?;
            let quality = source.vapor_quality().ok_or_else(|| MaterialError::State {
                what: "thermal water coexistence requires its original vapor mass quality".into(),
            })?;
            match closure {
                CoexistenceAcousticClosure::FrozenPhaseMasses
                    if quality == 0.0 || quality == 1.0 =>
                {
                    RegionalWaterAcousticResponse::SinglePhase(if quality == 0.0 {
                        pair.liquid()
                    } else {
                        pair.vapor()
                    })
                }
                CoexistenceAcousticClosure::FrozenPhaseMasses => {
                    RegionalWaterAcousticResponse::FrozenPhaseMasses(
                        FrozenHomogeneousMixture::from_mass_fraction(
                            pair.liquid().acoustic_state()?,
                            pair.vapor().acoustic_state()?,
                            quality,
                        )?,
                    )
                }
                CoexistenceAcousticClosure::HomogeneousEquilibrium => {
                    RegionalWaterAcousticResponse::HomogeneousEquilibrium(
                        HomogeneousEquilibriumWater::try_from_regional(pair, quality)?,
                    )
                }
            }
        };
        let selected = match response {
            RegionalWaterAcousticResponse::SinglePhase(single) => single.acoustic_state()?,
            RegionalWaterAcousticResponse::FrozenPhaseMasses(frozen) => frozen.state(),
            RegionalWaterAcousticResponse::HomogeneousEquilibrium(equilibrium) => {
                equilibrium.state()
            }
        };
        let density = source.bulk_density_kg_m3();
        let acoustic = if density.to_bits() == selected.density_kg_m3().to_bits() {
            selected
        } else {
            FluidAcousticState::from_adiabatic_bulk_modulus(
                density,
                selected.adiabatic_bulk_modulus_pa(),
            )?
        };
        Ok(Self {
            source,
            closure,
            response,
            acoustic,
            density_representation_difference_kg_m3: density - selected.density_kg_m3(),
        })
    }

    /// Complete thermal state, original source roots and caloric/pressure defects.
    pub const fn source(&self) -> &RegionalWaterState {
        &self.source
    }

    /// Explicit coexistence policy, retained even for a homogeneous state.
    pub const fn coexistence_closure(&self) -> CoexistenceAcousticClosure {
        self.closure
    }

    /// Original selected response; equilibrium tangents and source matching
    /// defects remain here. Its independently summed density may differ from
    /// `state()` by the reported representation roundoff.
    pub const fn response(&self) -> &RegionalWaterAcousticResponse {
        &self.response
    }

    /// Coherent acoustic coefficients using exactly the thermal storage density.
    pub const fn state(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Thermal `1/v` minus the selected response's independently formed density,
    /// kg/m³. This is numerical representation roundoff, never transported mass.
    pub const fn density_representation_difference_kg_m3(&self) -> f64 {
        self.density_representation_difference_kg_m3
    }
}

impl FluidAcoustics for RegionalWaterAcoustics {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

#[cfg(test)]
#[path = "thermal/tests.rs"]
mod tests;

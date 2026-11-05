//! Prescribed uniform temperature changes of the existing diffuse-interface field.
//!
//! Cooling/heating changes the original pair's chi, chemical potentials and
//! Maxwell-Stefan mobility. It does not flash cells into bulk equilibrium, seed
//! phases, discard interface history or advance the clock. Geometry, site
//! concentration, kappa and authored diffusivity stay fixed; temperature-dependent
//! diffusivity or thermal expansion require separately qualified source models.
//!
//! For the original chi = A + B/T law, fixed-composition mixing entropy is
//! independent of temperature. Consequently the external change of field free
//! energy is exactly -delta(T) * S_mix; unchanged gradient energy cancels before
//! summation. This is not heat input or dissipation. Original pure-liquid calorics
//! can supply the separate sensible heat/entropy account through [`caloric`].

use super::{
    BinaryMiscibilityNetwork, BinaryMiscibilityState, Face, Owner, Result, Sum, face_mobility,
    finite, poll,
};
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::solution::gradient::BinarySquareGradientV1;
use std::sync::Arc;

pub mod caloric;

/// External thermodynamic change at unchanged composition and physical time.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityTemperatureChange {
    pub temperature_before: Temperature,
    pub temperature_after: Temperature,
    pub chi_before: f64,
    pub chi_after: f64,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    /// Direct increment, preserving small changes beside large interface energy.
    /// Either sign is admissible; this quantity is neither supplied heat nor work.
    pub external_free_energy_change_j: f64,
    /// Homogeneous mixing entropy only, including the original excess entropy.
    /// Pure-component thermal entropy requires qualified caloric sources.
    pub mixing_entropy_j_per_k: f64,
}

/// All-or-nothing new thermal owner and the unchanged conserved spatial state.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityTemperatureCandidate {
    network: BinaryMiscibilityNetwork,
    state: BinaryMiscibilityState,
    change: BinaryMiscibilityTemperatureChange,
}
impl BinaryMiscibilityTemperatureCandidate {
    pub const fn network(&self) -> &BinaryMiscibilityNetwork {
        &self.network
    }
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub const fn change(&self) -> &BinaryMiscibilityTemperatureChange {
        &self.change
    }
    pub fn into_parts(
        self,
    ) -> (
        BinaryMiscibilityNetwork,
        BinaryMiscibilityState,
        BinaryMiscibilityTemperatureChange,
    ) {
        (self.network, self.state, self.change)
    }
}

impl BinaryMiscibilityNetwork {
    /// Change the externally imposed uniform temperature without reinitializing
    /// the field. Both inventories, compensated low words, independent fractions
    /// and elapsed time survive bit for bit. Subsequent transport must use the
    /// returned network; prepared flow can be explicitly rebound with
    /// [`Self::rebind_advection`].
    ///
    /// The original model's full temperature domain and positive-kappa spinodal
    /// admission still apply. Every mobility is recomputed from original c, D,
    /// area and distance, so temperature cycles never accumulate scaling error.
    /// Unchanged geometry/CSR arrays are shared without sorting or reconstruction.
    /// Work is O(cells+faces), within the original transaction workspace budget.
    pub fn temperature_candidate(
        &self,
        initial: &BinaryMiscibilityState,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityTemperatureCandidate> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        let original = self.model();
        let model = BinarySquareGradientV1::new(
            original.thermodynamics().clone(),
            temperature,
            AmountConcentration::new(original.molar_concentration_mol_m3()),
            original.gradient_coefficient_j_m(),
        )?;
        let network = if &model == original {
            self.clone()
        } else {
            let mut faces = Vec::with_capacity(self.owner.faces.len());
            for face in self.owner.faces.iter() {
                poll(&mut cancelled)?;
                let physical = face.physical;
                faces.push(Face {
                    physical,
                    gradient_j: face.gradient_j,
                    mobility_mol2_j_s: face_mobility(
                        &model,
                        physical.area_m2 / physical.distance_m,
                        physical.diffusivity_m2_s,
                    )?,
                });
            }
            Self {
                owner: Arc::new(Owner {
                    model,
                    cells: Arc::clone(&self.owner.cells),
                    capacities_mol: Arc::clone(&self.owner.capacities_mol),
                    faces: faces.into(),
                    adjacent_offsets: Arc::clone(&self.owner.adjacent_offsets),
                    adjacent_faces: Arc::clone(&self.owner.adjacent_faces),
                    gradient_diagonal_j: Arc::clone(&self.owner.gradient_diagonal_j),
                    workspace_bytes: self.owner.workspace_bytes,
                }),
            }
        };
        let mut entropy = Sum::default();
        for (&fractions, &capacity) in initial
            .fractions
            .iter()
            .zip(self.owner.capacities_mol.iter())
        {
            poll(&mut cancelled)?;
            entropy.add(finite(
                capacity
                    * original
                        .evaluate_composition(fractions)?
                        .mixing_entropy_j_per_mol_k(),
                "field mixing entropy",
            )?)?;
        }
        let mixing_entropy_j_per_k = entropy.value()?;
        let delta_temperature = finite(
            temperature.value() - original.temperature().value(),
            "prescribed temperature change",
        )?;
        let external_free_energy_change_j = finite(
            -delta_temperature * mixing_entropy_j_per_k,
            "external field free-energy increment",
        )?;
        let free_energy_before_j = self.free_energy_j(initial, &mut cancelled)?;
        let mut state = initial.clone();
        state.owner = Arc::clone(&network.owner);
        let free_energy_after_j = network.free_energy_j(&state, &mut cancelled)?;
        let change = BinaryMiscibilityTemperatureChange {
            temperature_before: original.temperature(),
            temperature_after: temperature,
            chi_before: original.chi(),
            chi_after: network.model().chi(),
            free_energy_before_j,
            free_energy_after_j,
            external_free_energy_change_j,
            mixing_entropy_j_per_k,
        };
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityTemperatureCandidate {
            network,
            state,
            change,
        })
    }
}

#[cfg(test)]
mod tests;

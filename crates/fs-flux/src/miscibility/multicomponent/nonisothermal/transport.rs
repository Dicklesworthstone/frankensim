//! Atomic prescribed stirring and molecular/heat transport over one interval.
//!
//! Advection and diffusion both advance over the same physical duration. This
//! Lie composition is first order; donor CFL subcycling alone does not estimate
//! the temporal error of the combined generator. A final original-state audit
//! applies one complete-interval species, enthalpy and entropy allowance.

use super::advection::{
    MulticomponentNonisothermalAdvectionControls, MulticomponentNonisothermalAdvectionStep,
};
use super::*;
use crate::miscibility::multicomponent::advection::MulticomponentMiscibilityAdvection;

/// Accuracy control of the complete stirred species/temperature evolution.
pub mod adaptive;

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalTransportStep {
    pub duration_s: f64,
    pub advection: MulticomponentNonisothermalAdvectionStep,
    pub diffusion: MulticomponentNonisothermalStep,
    /// Canonical physical face then species, signed left to right [mol].
    pub transferred_mol: Vec<Vec<f64>>,
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalTransportCandidate {
    field: MulticomponentNonisothermalField,
    step: MulticomponentNonisothermalTransportStep,
}
impl MulticomponentNonisothermalTransportCandidate {
    pub const fn field(&self) -> &MulticomponentNonisothermalField {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentNonisothermalTransportStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentNonisothermalField,
        MulticomponentNonisothermalTransportStep,
    ) {
        (self.field, self.step)
    }
}

impl MulticomponentNonisothermalField {
    /// Advect all original species/parcel enthalpy, then exchange molecular
    /// species and heat over one duration. Every stage retains the same source
    /// owner, fixed volumes, authored thermal face data and amount low words.
    pub fn transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection_controls: MulticomponentNonisothermalAdvectionControls,
        mut controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentNonisothermalTransportCandidate> {
        poll(&mut cancelled)?;
        let (mut advected, advection) = self
            .advect_candidate(
                flow,
                duration_s,
                advection_controls,
                controls,
                &mut cancelled,
            )?
            .into_parts();
        // Composition of generators changes the state twice but advances one
        // physical clock, including continuation from nonzero elapsed time.
        advected.state.elapsed_s = self.elapsed_s();
        let (field, diffusion) = advected
            .step_candidate(duration_s, controls, &mut cancelled)?
            .into_parts();
        controls.amount_balance_tolerance_mol = controls
            .amount_balance_tolerance_mol
            .min(advection_controls.amount_balance_tolerance_mol);
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let faces = self.network().faces().len();
        let n = self.network().species_count();
        let mut transferred_mol = Vec::with_capacity(faces);
        let mut transferred_energy_j = Vec::with_capacity(faces);
        for face in 0..faces {
            poll(&mut cancelled)?;
            let mut amounts = Vec::with_capacity(n);
            for species in 0..n {
                let mut amount = Sum::default();
                amount.add(advection.transferred_mol[face][species])?;
                amount.add(diffusion.transferred_mol[face][species])?;
                amounts.push(amount.value()?);
            }
            transferred_mol.push(amounts);
            let mut heat = Sum::default();
            heat.add(advection.transferred_energy_j[face])?;
            heat.add(diffusion.transferred_energy_j[face])?;
            transferred_energy_j.push(heat.value()?);
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentNonisothermalTransportCandidate {
            field,
            step: MulticomponentNonisothermalTransportStep {
                duration_s,
                advection,
                diffusion,
                transferred_mol,
                transferred_energy_j,
                component_residual_mol,
                maximum_capacity_residual_mol,
                energy_residual_j,
                entropy_change_j_per_k,
            },
        })
    }
}

#[cfg(test)]
mod tests;

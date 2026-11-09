//! Atomic composition of local-temperature advection and molecular/heat exchange.
//!
//! Both generators advance over the SAME interval. Their Lie composition is
//! first order; CFL subcycling alone does not control coupled temporal error.
//! The final original-state audit applies one energy/entropy/species allowance.

use super::advection::BinaryNonisothermalAdvectionStep;
use super::*;
use crate::miscibility::advection::{
    BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls,
};

/// Original stage receipts plus independently checked whole-interval balances.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalTransportStep {
    pub duration_s: f64,
    pub advection: BinaryNonisothermalAdvectionStep,
    pub diffusion: BinaryNonisothermalStep,
    pub transferred_component_moles: Vec<[f64; 2]>,
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

/// An accepted complete interval; no intermediate field escapes on failure.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalTransportCandidate {
    field: BinaryNonisothermalField,
    step: BinaryNonisothermalTransportStep,
}
impl BinaryNonisothermalTransportCandidate {
    pub const fn field(&self) -> &BinaryNonisothermalField {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalTransportStep {
        &self.step
    }
    pub fn into_parts(self) -> (BinaryNonisothermalField, BinaryNonisothermalTransportStep) {
        (self.field, self.step)
    }
}

impl BinaryNonisothermalField {
    /// Advect original species/enthalpy, then diffuse species/heat over one
    /// physical duration. Retains local temperatures, original material
    /// diffusivity and compensated inventories through both operators.
    pub fn transport_candidate(
        &self,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        mut controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryNonisothermalTransportCandidate> {
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
        // Two generators share one physical clock, including after continuation.
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
        let faces = self.owner.faces.len();
        let mut transferred_component_moles = Vec::with_capacity(faces);
        let mut transferred_energy_j = Vec::with_capacity(faces);
        for face in 0..faces {
            poll(&mut cancelled)?;
            let mut amount = [Sum::default(); 2];
            for (species, sign) in [1.0, -1.0].into_iter().enumerate() {
                amount[species].add(advection.transferred_component_moles[face][species])?;
                amount[species].add(sign * diffusion.transferred_a_mol[face])?;
            }
            transferred_component_moles.push([amount[0].value()?, amount[1].value()?]);
            let mut heat = Sum::default();
            heat.add(advection.transferred_energy_j[face])?;
            heat.add(diffusion.transferred_energy_j[face])?;
            transferred_energy_j.push(heat.value()?);
        }
        poll(&mut cancelled)?;
        Ok(BinaryNonisothermalTransportCandidate {
            field,
            step: BinaryNonisothermalTransportStep {
                duration_s,
                advection,
                diffusion,
                transferred_component_moles,
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

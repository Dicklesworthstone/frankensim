//! Atomic first-order prescribed stirring and full Maxwell–Stefan diffusion.
//!
//! Donor advection and coupled molecular diffusion each act over the same
//! physical interval. The second operator reads the complete first endpoint;
//! only the private trial clock is reset. Coupled temporal accuracy therefore
//! requires refinement of this Lie composition, independently of donor CFL.
//! Both stages and the complete convex-law endpoint retain their own original
//! conservation/free-energy gates. No chemical energy change is called heat or
//! stirring work.

/// Temporal refinement of the complete prescribed-flow and diffusion step.
pub mod adaptive;

use super::advection::{
    MulticomponentMiscibilityAdvection, MulticomponentMiscibilityAdvectionControls,
    MulticomponentMiscibilityAdvectionStep,
};
use super::{
    MulticomponentMiscibilityField, MulticomponentMiscibilitySolve, MulticomponentMiscibilityStep,
    Result, Sum, admit_controls, poll,
};

/// Complete stage reports and independently audited original-to-final balances.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityTransportStep {
    pub duration_s: f64,
    pub advection: MulticomponentMiscibilityAdvectionStep,
    pub diffusion: MulticomponentMiscibilityStep,
    /// Both stage contributions in canonical face/species order [mol].
    pub transferred_mol: Vec<Vec<f64>>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_change_j: f64,
    /// Direct whole-step change minus both stage changes; a roundoff account.
    pub stage_energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityTransportCandidate {
    field: MulticomponentMiscibilityField,
    step: MulticomponentMiscibilityTransportStep,
}
impl MulticomponentMiscibilityTransportCandidate {
    pub const fn field(&self) -> &MulticomponentMiscibilityField {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentMiscibilityTransportStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentMiscibilityField,
        MulticomponentMiscibilityTransportStep,
    ) {
        (self.field, self.step)
    }
}

impl MulticomponentMiscibilityField {
    /// Advect then diffuse over one physical interval, preserving every
    /// component and its low word. The stricter stage amount/free-energy
    /// allowances apply once again to the complete interval; they are never
    /// multiplied by faces, donor substeps or operator count.
    pub fn transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection_controls: MulticomponentMiscibilityAdvectionControls,
        diffusion_controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityTransportCandidate> {
        poll(&mut cancelled)?;
        admit_controls(diffusion_controls)?;
        let (mut advected, advection) = self
            .advect_candidate(flow, duration_s, advection_controls, &mut cancelled)?
            .into_parts();
        advected.elapsed_s = self.elapsed_s;
        let (field, diffusion) = advected
            .step_candidate(duration_s, diffusion_controls, &mut cancelled)?
            .into_parts();
        let (component_residual_mol, maximum_capacity_residual_mol, free_energy_change_j) = self.balances(
            &field,
            MulticomponentMiscibilitySolve {
                amount_balance_tolerance_mol: advection_controls.amount_balance_tolerance_mol
                    .min(diffusion_controls.amount_balance_tolerance_mol),
                free_energy_tolerance_j: advection_controls.free_energy_tolerance_j
                    .min(diffusion_controls.free_energy_tolerance_j),
                ..diffusion_controls
            },
            &mut cancelled,
        )?;
        let mut transferred_mol = Vec::with_capacity(self.network.faces().len());
        for (advected, diffused) in advection.transferred_mol.iter().zip(&diffusion.transferred_mol) {
            poll(&mut cancelled)?;
            let mut row = Vec::with_capacity(self.network.species_count());
            for (&advected, &diffused) in advected.iter().zip(diffused) {
                let mut transfer = Sum::default();
                transfer.add(advected)?;
                transfer.add(diffused)?;
                row.push(transfer.value()?);
            }
            transferred_mol.push(row);
        }
        let mut energy = Sum::default();
        energy.add(free_energy_change_j)?;
        energy.add(-advection.free_energy_change_j)?;
        energy.add(-diffusion.free_energy_change_j)?;
        let stage_energy_residual_j = energy.value()?;
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityTransportCandidate {
            field,
            step: MulticomponentMiscibilityTransportStep {
                duration_s,
                advection,
                diffusion,
                transferred_mol,
                component_residual_mol,
                maximum_capacity_residual_mol,
                free_energy_change_j,
                stage_energy_residual_j,
            },
        })
    }
}

#[cfg(test)]
mod tests;

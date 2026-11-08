//! Recover original physical Laurent coefficients without changing their units.

use super::{C64, response};
use fs_modal::relaxation_modes::MaxwellForceResidue;

/// Coefficients of 1/(s-s0), not response amplitudes at a real frequency.
/// For each field its residue has that field's response unit multiplied by
/// inverse seconds. For example pressure_residue_pa_per_s/(s-s0) is Pa.
#[derive(Clone, Debug)]
pub struct StructuralInteriorRelaxationForceResidue {
    pub laplace_pole_per_s: C64,
    /// ell^T f/(ell^T T'(s0)q), in the retained right mode's normalization.
    pub participation: C64,
    pub structural_displacement_residue: Vec<C64>,
    pub structural_velocity_residue: Vec<C64>,
    pub pressure_potential_residue_pa: Vec<C64>,
    /// s0*potential residue, including the real decay contribution.
    pub pressure_residue_pa_per_s: Vec<C64>,
    /// Original B^T*pressure residue, excluding explicit branch forces.
    pub fluid_force_residue: Vec<C64>,
    /// Original B*velocity residue, excluding explicit branch flow.
    pub displaced_flow_residue_m3_s2: Vec<C64>,
    /// Original stress residue for every retained branch, in source order.
    pub branch_stress_residue: Vec<C64>,
}

pub(super) fn observe(
    coupling: &[f64],
    ns: usize,
    na: usize,
    mut native: MaxwellForceResidue,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StructuralInteriorRelaxationForceResidue, response::ResponseError> {
    if cancelled() {
        return Err(response::ResponseError::Cancelled);
    }
    if ns.checked_add(na) != Some(native.displacement.len())
        || native.velocity.len() != native.displacement.len()
    {
        return Err(response::ResponseError::Invalid(
            "native residue must retain every original coordinate",
        ));
    }
    let potential = native.displacement.split_off(ns);
    let pressure = native.velocity.split_off(ns);
    let (force, flow) = response::interface(coupling, &native.velocity, &pressure, cancelled)?;
    if cancelled() {
        return Err(response::ResponseError::Cancelled);
    }
    Ok(StructuralInteriorRelaxationForceResidue {
        laplace_pole_per_s: native.laplace_frequency_per_s,
        participation: native.participation,
        structural_displacement_residue: native.displacement,
        structural_velocity_residue: native.velocity,
        pressure_potential_residue_pa: potential,
        pressure_residue_pa_per_s: pressure,
        fluid_force_residue: force,
        displaced_flow_residue_m3_s2: flow,
        branch_stress_residue: native.branch_stress,
    })
}

//! Reference material inventory -> original distributed-load integration.
use fs_exec::Cx;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_solid::linear3::{TetAssemblyBudget, TetThermalStressReport};
use fs_solid::loads3::{LoadError3, TetLoadGeometry3, VolumeLoad3};

use super::{TetrahedralMass, ThermalBodyError, checkpoint, invalid};

impl TetrahedralMass<'_> {
    /// Integrate gravity using each cell's original material density and the
    /// existing solid P1 volume-load owner. The returned sorted force vector can
    /// enter `solve_thermal_with_nodal_forces` or a constrained heat history.
    ///
    /// Acceleration is expressed in the body's reference frame [m/s^2]. Closed
    /// matter retains its reference mass at all later temperatures, so neither
    /// current density nor current deformed volume is queried here. This fixed
    /// reference-force route does not rotate gravity with a moving body.
    /// The original mesh/face-entry budget and cancellation boundaries apply.
    pub fn gravity_nodal_forces(
        &self,
        acceleration_m_s2: [f64; 3],
        budget: TetAssemblyBudget,
        maximum_face_entries: usize,
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, ThermalBodyError> {
        checkpoint(cx)?;
        if acceleration_m_s2.iter().any(|x| !x.is_finite()) {
            return Err(LoadError3::Invalid("gravity acceleration must be finite").into());
        }
        let geometry = TetLoadGeometry3::new(
            self.nodes_m(),
            self.tetrahedra(),
            budget,
            maximum_face_entries,
            cx,
        )?;
        let mut force_density = vec![[0.0; 3]; self.tetrahedra().len()];
        for (region, &element) in self
            .reference()
            .regions()
            .iter()
            .zip(self.region_elements())
        {
            checkpoint(cx)?;
            let density = region
                .material()
                .property(DENSITY_PROPERTY)
                .expect("reference material density resolved")
                .value_si();
            force_density[element] = acceleration_m_s2.map(|a| a * density);
        }
        let load = geometry.assemble(&[], VolumeLoad3::PerElement(&force_density), &[], cx)?;
        let mut forces = Vec::with_capacity(load.nodal_forces_n().len());
        for &force in load.nodal_forces_n() {
            checkpoint(cx)?;
            forces.push(force);
        }
        checkpoint(cx)?;
        Ok(forces)
    }
}

// The solid owner supplies integral B^T sigma. A mechanical load applied at a
// support must be subtracted to observe the actual support force, not silently
// discarded because that displacement component was constrained in the solve.
pub(super) fn support_reactions(
    stress: &TetThermalStressReport,
    fixed_dofs: &[usize],
    nodal_forces_n: &[(usize, f64)],
    cx: &Cx<'_>,
) -> Result<Vec<[f64; 3]>, ThermalBodyError> {
    checkpoint(cx)?;
    let mut reactions = vec![[0.0; 3]; stress.nodal_internal_force_n().len()];
    for &dof in fixed_dofs {
        checkpoint(cx)?;
        let applied = nodal_forces_n
            .binary_search_by_key(&dof, |&(key, _)| key)
            .ok()
            .map_or(0.0, |i| nodal_forces_n[i].1);
        let reaction = stress.nodal_internal_force_n()[dof / 3][dof % 3] - applied;
        if !reaction.is_finite() {
            return Err(invalid("support reaction is unrepresentable"));
        }
        reactions[dof / 3][dof % 3] = reaction;
    }
    checkpoint(cx)?;
    Ok(reactions)
}

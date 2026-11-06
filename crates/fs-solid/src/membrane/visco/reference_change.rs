//! Prescribed natural-strain events on the original membrane and material memory.
//!
//! Actual nodal positions, the original reference mesh, mass and physical clock
//! remain fixed. Each original material point changes its prescribed reference
//! and then restores its transverse plane-stress equilibrium through the native
//! zero-duration transaction. Viscous metrics are neither reset nor relaxed.
//!
//! Exact elastic volume includes the original signed pressure-volume reference
//! actuation. Reference work, transverse mechanical work and numerical splitting
//! error are separate; none becomes relaxation heat. This event may create a
//! global force imbalance. A subsequent `advance` owns the supported or inertial
//! mechanical response. No continuous thermal path or caloric feedback is inferred.

use std::sync::Arc;

use fs_material::sheet::SheetError;

use super::dynamics::FiniteViscoMembraneMotion;
use super::{
    FiniteViscoMembrane, FiniteViscoMembraneError, FiniteViscoStepControls, MembranePatch,
    ViscoelasticLawField, admit_sheet_binding, energy_gate, invalid, material_error, poll,
};

/// Physical and numerical work of one zero-duration reference event [J].
/// These are increments, separate from cumulative native material ledgers.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteViscoMembraneReferenceEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Exactly zero: no physical time or branch relaxation occurs in this event.
    pub dissipated_j: f64,
    /// Signed native input at fixed full deformation and viscous metrics.
    pub reference_work_j: f64,
    /// Signed exact-volume actuation, +p_endpoint Delta det(G), integrated
    /// over original volume. Zero for finite compliance. Total physical
    /// reference input is this plus `reference_work_j`, counted once.
    pub constraint_reference_work_j: f64,
    /// Native transverse frozen-metric loading increment. The exact-volume
    /// ledger includes constraint actuation here; this is not another physical
    /// source to add independently to the reference and transverse work.
    pub native_transverse_loading_work_j: f64,
    /// Physical full-Piola transverse endpoint work. Traction-free transverse
    /// equilibrium makes this zero within the local solve accuracy.
    pub transverse_endpoint_work_j: f64,
    /// Fixed-full-F reference input plus native loading minus storage change.
    pub material_balance_error_j: f64,
    /// Transverse endpoint work plus constraint input minus native loading.
    /// Signed quadrature/splitting error; never physical relaxation heat.
    pub endpoint_work_defect_j: f64,
    /// Sum of absolute point defects. The interval defect allowance gates this
    /// quantity so opposing regional reference events cannot hide local error.
    pub absolute_endpoint_work_defect_j: f64,
}

/// Complete locally admitted proposal at the unchanged nodal configuration.
/// Global force equilibrium is assessed by a subsequent native mechanical step.
#[derive(Debug)]
pub struct FiniteViscoMembraneReferenceChange<'a> {
    pub membrane: FiniteViscoMembrane<'a>,
    pub energy: FiniteViscoMembraneReferenceEnergy,
}

/// Prescribed reference work with unchanged velocity and original reference mass.
/// The event supplies no mechanical impulse, kinetic-energy change or heat.
#[derive(Debug)]
pub struct FiniteViscoMembraneMotionReferenceChange<'a> {
    pub motion: FiniteViscoMembraneMotion<'a>,
    pub energy: FiniteViscoMembraneReferenceEnergy,
    /// Exactly zero. Every original nodal momentum is retained, not just its sum.
    pub impulse_kg_m_s: [f64; 3],
    /// Exactly zero because both velocity and the complete reference mass survive.
    pub kinetic_energy_change_j: f64,
}

fn admit_original_patch(
    previous: &MembranePatch,
    target: &MembranePatch,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FiniteViscoMembraneError> {
    poll(cancelled)?;
    if core::ptr::eq(previous, target) {
        return Ok(());
    }
    if previous.mesh.nodes.len() != target.mesh.nodes.len()
        || previous.elements.len() != target.elements.len()
        || previous.mesh.patches.len() != target.mesh.patches.len()
        || previous.domain.minimum_area_ratio != target.domain.minimum_area_ratio
        || previous.domain.maximum_abs_strain != target.domain.maximum_abs_strain
    {
        return Err(invalid(
            "reference events require the original membrane topology and strain envelope",
        ));
    }
    for (node, (old, new)) in previous
        .mesh
        .nodes
        .iter()
        .zip(&target.mesh.nodes)
        .enumerate()
    {
        if node % 256 == 0 {
            poll(cancelled)?;
        }
        if old != new || previous.nodal_mass_kg[node] != target.nodal_mass_kg[node] {
            return Err(invalid(
                "reference events cannot change original membrane coordinates or mass",
            ));
        }
    }
    for ((old_patch, old_edges), (new_patch, new_edges)) in
        previous.mesh.patches.iter().zip(&target.mesh.patches)
    {
        poll(cancelled)?;
        if old_patch != new_patch || old_edges.len() != new_edges.len() {
            return Err(invalid("reference events retain original boundary patches"));
        }
        for (index, (old, new)) in old_edges.iter().zip(new_edges).enumerate() {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            if old != new {
                return Err(invalid("reference events retain original boundary edges"));
            }
        }
    }
    for (index, (old, new)) in previous.elements.iter().zip(&target.elements).enumerate() {
        poll(cancelled)?;
        if old.nodes != new.nodes
            || old.gradients != new.gradients
            || old.area != new.area
            || old.grain_sin_cos != new.grain_sin_cos
            || previous.grain_rad[index] != target.grain_rad[index]
            || previous.reference_areal_density_kg_m2[index]
                != target.reference_areal_density_kg_m2[index]
        {
            return Err(invalid(
                "reference events retain original triangle geometry, material axes and corner inventory",
            ));
        }
    }
    poll(cancelled)?;
    Ok(())
}

impl FiniteViscoMembrane<'_> {
    /// Change prescribed stress-free strain while retaining the actual membrane
    /// configuration, material-frame viscous metrics and physical time.
    ///
    /// Build the target patch from the original source bundles, mesh, grain,
    /// thickness, corner density and strain envelope. Its selected natural
    /// stretches must match the target laws exactly. Original equilibrium and
    /// ordered branch energies/domain remain fixed; the native material owner
    /// admits changed relaxation times and prescribed reference deformation.
    /// Full and exact-volume treatments cannot switch during this event.
    ///
    /// The complete original source bundle and selected model card are retained;
    /// a target thermal reference may carry a new admitted endpoint context.
    /// An explicit natural reference does not acquire thermal source evidence.
    /// This is a parameter/transverse-equilibrium split, not a coordinate-frame
    /// remapping or an implicitly resolved finite thermal path. Reduce the event
    /// when its separately gated quadrature defect is too large.
    ///
    /// Failure, late cancellation or a work refusal publishes no candidate.
    /// The returned body borrows only the target patch/laws; original accepted
    /// buffers are shared or copied without imposing their former borrow lifetime.
    pub fn change_reference<'next>(
        &self,
        target_patch: &'next MembranePatch,
        target_material: ViscoelasticLawField<'next>,
        controls: FiniteViscoStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteViscoMembraneReferenceChange<'next>, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        if ![
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
        ]
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0)
        {
            return Err(invalid(
                "finite nonnegative membrane reference-work allowances required",
            ));
        }
        admit_original_patch(self.patch, target_patch, &mut cancelled)?;
        if let ViscoelasticLawField::PerElement(laws) = target_material
            && laws.len() != self.states.len()
        {
            return Err(invalid(
                "one original target law per membrane triangle required",
            ));
        }
        let count = self.states.len();
        let mut states = Vec::with_capacity(count);
        let mut first_piola_pa = Vec::with_capacity(count);
        let mut pressure_pa = Vec::with_capacity(count);
        let mut resultants_n_m = Vec::with_capacity(count);
        let mut energy = FiniteViscoMembraneReferenceEnergy {
            stored_before_j: self.stored_energy_j,
            ..Default::default()
        };
        for (index, element) in target_patch.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let target = target_material.at(index);
            let sheet = admit_sheet_binding(target_patch, index, target)?;
            let previous_sheet = self.sheet(index);
            if sheet.material() != previous_sheet.material()
                || sheet.model().model_card() != previous_sheet.model().model_card()
                || sheet.reference_thickness_m() != previous_sheet.reference_thickness_m()
                || sheet.volumetric_treatment() != previous_sheet.volumetric_treatment()
            {
                return Err(material_error(
                    index,
                    SheetError::Invalid(
                        "reference events retain the complete original sheet source, model, thickness and volume treatment",
                    ),
                ));
            }
            let strain = target_patch.admit_kinematics(
                index,
                target_patch.kinematics(element, &self.displacement_m),
            )?;
            let change = target
                .change_membrane_reference_from(
                    self.material.at(index),
                    &self.states[index],
                    strain,
                    sheet.reference_thickness_m(),
                    sheet.settings(),
                )
                .map_err(|source| material_error(index, source))?;
            energy.stored_after_j += element.area * change.stored_after_j_m2;
            energy.reference_work_j += element.area * change.reference_work_j_m2;
            energy.constraint_reference_work_j +=
                element.area * change.constraint_reference_work_j_m2;
            energy.native_transverse_loading_work_j +=
                element.area * change.native_transverse_loading_work_j_m2;
            energy.transverse_endpoint_work_j +=
                element.area * change.transverse_endpoint_work_j_m2;
            energy.endpoint_work_defect_j += element.area * change.endpoint_work_defect_j_m2;
            energy.absolute_endpoint_work_defect_j +=
                element.area * change.endpoint_work_defect_j_m2.abs();
            states.push(change.state);
            first_piola_pa.push(change.first_piola_pa);
            pressure_pa.push(change.pressure_pa);
            resultants_n_m.push(change.sheet.resultants_n_m);
        }
        energy.material_balance_error_j = energy.reference_work_j
            + energy.native_transverse_loading_work_j
            - (energy.stored_after_j - energy.stored_before_j);
        energy_gate(
            "membrane reference-event material balance",
            energy.material_balance_error_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "membrane reference-event endpoint integration defect",
            energy.absolute_endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if ![
            energy.stored_after_j,
            energy.reference_work_j,
            energy.constraint_reference_work_j,
            energy.native_transverse_loading_work_j,
            energy.transverse_endpoint_work_j,
            energy.endpoint_work_defect_j,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(invalid("unrepresentable membrane reference-event work"));
        }
        poll(&mut cancelled)?;
        Ok(FiniteViscoMembraneReferenceChange {
            membrane: FiniteViscoMembrane {
                patch: target_patch,
                material: target_material,
                displacement_m: Arc::clone(&self.displacement_m),
                states: states.into(),
                first_piola_pa: first_piola_pa.into(),
                pressure_pa: pressure_pa.into(),
                resultants_n_m: resultants_n_m.into(),
                stored_energy_j: energy.stored_after_j,
                time_s: self.time_s,
            },
            energy,
        })
    }
}

impl FiniteViscoMembraneMotion<'_> {
    /// Apply the same prescribed reference event with unchanged displacement,
    /// velocity, time and original reference inertia. It creates no impulse or
    /// relaxation heat. A subsequent inertial `advance` resolves any changed
    /// material forces with these retained velocities and viscous histories.
    pub fn change_reference<'next>(
        &self,
        target_patch: &'next MembranePatch,
        target_material: ViscoelasticLawField<'next>,
        controls: FiniteViscoStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteViscoMembraneMotionReferenceChange<'next>, FiniteViscoMembraneError> {
        let change = self.membrane().change_reference(
            target_patch,
            target_material,
            controls,
            &mut cancelled,
        )?;
        let motion = FiniteViscoMembraneMotion::try_new(
            change.membrane,
            self.velocity_mps(),
            self.mass_image(),
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(FiniteViscoMembraneMotionReferenceChange {
            motion,
            energy: change.energy,
            impulse_kg_m_s: [0.0; 3],
            kinetic_energy_change_j: 0.0,
        })
    }
}

//! Source-temperature refresh inside the original rotor feedback solve.

use super::{
    Result, RotorControls, SpatialInductionRotor, SpatialRotorCandidate, SpatialRotorError,
};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    vector_potential::{VectorMagneticBoundary, induction::SpatialInductionControls},
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{conductor::ResolvedConductorState, field_storage::spatial::SpatialMagneticLaw};

impl<'m> SpatialInductionRotor<'m> {
    /// Advance the same accepted rotor with source-pinned endpoint materials.
    ///
    /// The original field owner admits the exact material identity, temperature
    /// refresh, electrical basis and unchanged conducting support. One operator
    /// preparation serves every speed trial. The rotor additionally retains
    /// homogeneous scalar conductivity and uniform scalar permeability: a
    /// nonuniform or oriented endpoint would require transported material or
    /// magnetic-interface loads outside this fixed-domain motion model.
    ///
    /// This method returns physical heat and parameter work; a finite thermal
    /// driver must own their caloric and magnetic-entropy accounts. It does not
    /// choose a temperature or deposit the reported heat a second time.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_candidate(
        &self,
        target_s: f64,
        endpoint_current_a_m2: &[[f64; 3]],
        shaft_torque_nm: f64,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<ResolvedConductorState>>,
        assembly: WeightedAssemblyLimits,
        field_controls: SpatialInductionControls,
        controls: RotorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialRotorCandidate<'m>> {
        self.step_materials_boundary_candidate(
            target_s,
            endpoint_current_a_m2,
            None,
            shaft_torque_nm,
            materials,
            conductors,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// Refresh the original material sources while an applied magnetic trace
    /// drives the finite-inertia rotor. The one endpoint operator preparation
    /// serves every speed trial; boundary supply and magnetic parameter work
    /// remain separate in the original native report.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_and_boundary_candidate(
        &self,
        target_s: f64,
        endpoint_current_a_m2: &[[f64; 3]],
        endpoint_boundary: &VectorMagneticBoundary,
        shaft_torque_nm: f64,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<ResolvedConductorState>>,
        assembly: WeightedAssemblyLimits,
        field_controls: SpatialInductionControls,
        controls: RotorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialRotorCandidate<'m>> {
        self.step_materials_boundary_candidate(
            target_s,
            endpoint_current_a_m2,
            Some(endpoint_boundary),
            shaft_torque_nm,
            materials,
            conductors,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_materials_boundary_candidate(
        &self,
        target_s: f64,
        endpoint_current_a_m2: &[[f64; 3]],
        endpoint_boundary: Option<&VectorMagneticBoundary>,
        shaft_torque_nm: f64,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<ResolvedConductorState>>,
        assembly: WeightedAssemblyLimits,
        field_controls: SpatialInductionControls,
        controls: RotorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialRotorCandidate<'m>> {
        super::poll(cancelled)?;
        crate::induction_rotor::admit_controls(
            target_s - self.field.time_s(),
            shaft_torque_nm,
            controls,
        )?;
        let prepared = self
            .field
            .prepare_materials_trial(materials, conductors, assembly, cancelled)?;
        let first_cell = self.binding.cells[0];
        let first =
            prepared.conductors()[first_cell]
                .as_ref()
                .ok_or(SpatialRotorError::Invalid(
                    "thermal rotor lost its original conductor",
                ))?;
        let problem = prepared.problem();
        let permeability = problem.permeability()[first_cell];
        if permeability.isotropic_h_m().is_none() {
            return Err(SpatialRotorError::Invalid(
                "thermal rotor requires scalar permeability",
            ));
        }
        for cell in 0..problem.mesh().tets.len() {
            super::poll(cancelled)?;
            if problem.magnetic_laws()[cell] != SpatialMagneticLaw::Isotropic
                || problem.permeability()[cell] != permeability
                || prepared.conductors()[cell]
                    .as_ref()
                    .is_some_and(|source| source != first)
            {
                return Err(SpatialRotorError::Invalid(
                    "thermal rotor requires homogeneous conductor and uniform permeability",
                ));
            }
        }
        self.step_with(
            target_s,
            endpoint_current_a_m2,
            shaft_torque_nm,
            controls,
            |velocity, checkpoint| {
                match endpoint_boundary {
                    Some(boundary) => prepared.step_with_boundary_candidate(
                        endpoint_current_a_m2,
                        velocity,
                        boundary,
                        target_s - self.field.time_s(),
                        field_controls,
                        &mut || checkpoint(),
                    ),
                    None => prepared.step_candidate(
                        endpoint_current_a_m2,
                        velocity,
                        target_s - self.field.time_s(),
                        field_controls,
                        &mut || checkpoint(),
                    ),
                }
                .map_err(Into::into)
            },
            cancelled,
        )
    }
}

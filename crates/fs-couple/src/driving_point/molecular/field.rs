//! One source-bound continuous gas field for feedback and every observed mouth.
use super::*;
use fs_duct::profile::RiccatiSteps;
use fs_material::gas_field::AxialGasField;
use fs_material::relaxation::RelaxationAbsorption;

impl MolecularDrivingPoint {
    /// Bind the original primitive temperature/composition field and molecular
    /// law to characteristic pressure feedback and physical opening flows.
    ///
    /// The same fourth-order field/geometry quadrature supplies propagation
    /// and the integral of frozen molecular slowness used for causal delays.
    /// The exact inlet state drives the exciter reference; every chimney uses
    /// its own local state. Radiation uses only the explicit exterior gas.
    /// Uniform fields retain the original homogeneous numerical path.
    ///
    /// `steps_per_span` is the original axial Magnus refinement axis, separate
    /// from wall diffusion poles and the consumer's FIR-length ladder. All
    /// original field, source, common-pressure, composition-work and internal
    /// mutual-hole validity refusals remain binding. This is a prescribed,
    /// frozen field, not an evolving transport or nonlinear gas solver.
    #[allow(clippy::too_many_arguments)]
    pub fn from_gas_field(
        geometry: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        wall: Option<MolecularWallResolution>,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        composition: ComposeOptions,
        exterior: GasState,
    ) -> Result<Self, DrivingPointError> {
        Self::from_field_bound(
            geometry,
            field,
            law,
            wall,
            None,
            steps_per_span,
            uniform_cones,
            composition,
            exterior,
        )
    }

    /// The same continuous molecular field with one source-resolved material
    /// wall shared by the bore, finite chimneys, feedback and observation.
    /// The original `LocalMaterialWall` admission and source band are retained.
    /// Its strictly positive surface mass preserves the gas-frozen wavefront.
    #[allow(clippy::too_many_arguments)]
    pub fn from_gas_field_with_material_wall(
        geometry: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        wall: Option<MolecularWallResolution>,
        material_wall: LocalMaterialWall,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        composition: ComposeOptions,
        exterior: GasState,
    ) -> Result<Self, DrivingPointError> {
        Self::from_field_bound(
            geometry,
            field,
            law,
            wall,
            Some(Arc::new(material_wall)),
            steps_per_span,
            uniform_cones,
            composition,
            exterior,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_field_bound(
        geometry: &Duct,
        field: &AxialGasField,
        law: &RelaxationAbsorption,
        wall: Option<MolecularWallResolution>,
        material_wall: Option<Arc<LocalMaterialWall>>,
        steps_per_span: RiccatiSteps,
        uniform_cones: ConeSlices,
        composition: ComposeOptions,
        exterior: GasState,
    ) -> Result<Self, DrivingPointError> {
        let inlet = field.state_at(0.0).map_err(DrivingPointError::Material)?;
        Self::admit_geometry_exterior(geometry, inlet.state(), &exterior)?;
        let material = match &material_wall {
            Some(material) => MolecularDuct::from_gas_field_with_compliant_wall(
                geometry,
                field,
                law,
                wall,
                material.clone(),
                steps_per_span,
                uniform_cones,
                composition.maximum_blocks,
            )?,
            None => MolecularDuct::from_gas_field(
                geometry,
                field,
                law,
                wall,
                steps_per_span,
                uniform_cones,
                composition.maximum_blocks,
            )?,
        };
        Self::from_prepared(geometry, material, material_wall, composition, exterior)
    }
}

#[cfg(test)]
mod tests;

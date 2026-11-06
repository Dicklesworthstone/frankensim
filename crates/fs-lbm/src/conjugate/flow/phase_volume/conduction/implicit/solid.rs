//! Finite fixed-solid storage in the same fluid/wall implicit thermal graph.

use super::*;
use crate::conjugate::flow::{
    StaggeredBoundaryPressure, ThermalVolumeMomentumInput, ThermalVolumeMomentumStep,
};
use fs_material::phase::EquilibriumPhaseState;

/// Explicit caloric inventory for one fixed-solid voxel. Its actual mass and
/// geometry are held; chart density never creates solid motion or fluid mass.
/// The same original source can be borrowed by any number of wall voxels.
#[derive(Clone, Copy)]
pub struct ThermalVolumeSolidCell<'a, 'm> {
    pub cell: usize,
    pub conductivity: &'a PreparedPhaseConductivity<'m>,
    pub mass_kg: f64,
    pub phase: EquilibriumPhaseState,
    pub net_enthalpy_increment_j: f64,
}

/// Accepted wall state, in increasing voxel order, independent of fluid routes.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalVolumeSolidEndpoint {
    pub cell: usize,
    pub mass_kg: f64,
    pub phase: EquilibriumPhaseState,
    pub conductivity_source_identity: ContentHash,
    /// Net physical heat from adjacent cells and prescribed baths, excluding
    /// this wall cell's independently authored caloric increment [J].
    pub conductive_enthalpy_increment_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConjugateThermalVolumeReport {
    /// The implicit audit covers ALL participating fluid and wall cells. Its
    /// total defect also includes the subsequent original fluid transport.
    /// Individual wall source identities are retained in their endpoints.
    pub conduction: ImplicitConductiveThermalVolumeReport,
    pub solid_authored_heat_j: f64,
    pub solid_enthalpy_change_j: f64,
    /// Signed physical interface heat, positive from finite solids into fluid.
    pub solid_to_fluid_heat_j: f64,
    /// Signed inward bath heat on explicitly bound solids [J], already included
    /// in `conduction.boundary_heat_j`. Solid storage obeys authored heat plus
    /// this boundary heat minus `solid_to_fluid_heat_j`, within the shared gate.
    pub solid_boundary_heat_j: f64,
    /// Sum of absolute solid-bath face heats, without cancellation [J].
    pub absolute_solid_boundary_heat_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConjugateThermalVolumeStep {
    pub thermal: ThermalVolumeStep,
    /// Fluid entries only; every fixed-solid entry is zero.
    pub conductive_enthalpy_increment_j: Vec<f64>,
    pub solid_endpoints: Vec<ThermalVolumeSolidEndpoint>,
    pub report: ImplicitConjugateThermalVolumeReport,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConjugateThermalVolumeMomentumStep {
    pub evolution: ThermalVolumeMomentumStep,
    pub conductive_enthalpy_increment_j: Vec<f64>,
    pub solid_endpoints: Vec<ThermalVolumeSolidEndpoint>,
    pub report: ImplicitConjugateThermalVolumeReport,
}

impl<T> ImplicitAdvance<T> {
    pub(super) fn into_conjugate(
        self,
    ) -> (
        T,
        Vec<f64>,
        Vec<ThermalVolumeSolidEndpoint>,
        ImplicitConjugateThermalVolumeReport,
    ) {
        (
            self.evolution,
            self.increments,
            self.solid_endpoints,
            ImplicitConjugateThermalVolumeReport {
                conduction: self.report,
                solid_authored_heat_j: self.solid_authored_heat_j,
                solid_enthalpy_change_j: self.solid_enthalpy_change_j,
                solid_to_fluid_heat_j: self.solid_to_fluid_heat_j,
                solid_boundary_heat_j: self.solid_boundary_heat_j,
                absolute_solid_boundary_heat_j: self.absolute_solid_boundary_heat_j,
            },
        )
    }
}

impl PreparedImplicitHeat {
    pub(super) fn into_conjugate(
        self,
    ) -> (
        Vec<f64>,
        Vec<f64>,
        Vec<ThermalVolumeSolidEndpoint>,
        ImplicitConjugateThermalVolumeReport,
    ) {
        (
            self.combined,
            self.increments,
            self.solid_endpoints,
            ImplicitConjugateThermalVolumeReport {
                conduction: self.report,
                solid_authored_heat_j: self.solid_authored_heat_j,
                solid_enthalpy_change_j: self.solid_enthalpy_change_j,
                solid_to_fluid_heat_j: self.solid_to_fluid_heat_j,
                solid_boundary_heat_j: self.solid_boundary_heat_j,
                absolute_solid_boundary_heat_j: self.absolute_solid_boundary_heat_j,
            },
        )
    }
}

impl FlowField {
    /// Evolve finite walls and fluid on ONE implicit enthalpy graph, then
    /// advance the original phase-density volume and conservative fluid routes.
    ///
    /// Adjacent participating voxels conduct through both actual half cells
    /// and any explicitly supplied fixed or source-owned face resistance. Each
    /// cell retains its source k(T,f), tangent and original latent/sensible law.
    /// Unlisted fixed solids and unselected exposed faces are insulated. Baths
    /// may heat or cool either fluid or explicitly bound solid cells through
    /// their own half-cell source law. A bath cannot also own an internal face
    /// between active bodies. Wall authored and boundary heat each enter once.
    ///
    /// Solid mass, geometry and occupancy stay fixed, including across a wall's
    /// latent interval. Source contact tangents enter the same native Newton
    /// system; external pressure and context stay with their original owner.
    /// No melting motion, thermal expansion or wall advection is inferred.
    /// Returned walls are separate inventories; fluid solid-slot mass/h remain empty.
    /// Workspace needs one entry per fluid cell plus explicitly supplied wall.
    /// The original clock, total heat budget and first-order heat/flow split
    /// remain authoritative. Empty solids/contacts preserve the existing path.
    ///
    /// # Errors
    /// Invalid ownership/contact/source, any thermal or native-flow refusal,
    /// or cancellation publishes neither a fluid nor a wall endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConjugateThermalVolumeStep, ChtError> {
        let accepted = self.implicit_conjugate_then_advance(
            domain,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            controls,
            solve,
            workspace,
            gate,
            |heated| {
                self.advance_piecewise_thermal_volume(
                    domain,
                    faces,
                    conductivity.caloric_curve(),
                    heated,
                    controls.volume,
                    gate,
                )
            },
            |thermal| thermal,
        )?;
        let (thermal, conductive_enthalpy_increment_j, solid_endpoints, report) =
            accepted.into_conjugate();
        Ok(ImplicitConjugateThermalVolumeStep {
            thermal,
            conductive_enthalpy_increment_j,
            solid_endpoints,
            report,
        })
    }

    /// The same finite-wall interval with the original accepted pressure impulse
    /// and material momentum transport. Thermal wall exchange supplies no mass
    /// or mechanical impulse. Optional reservoir pressure retains its owner.
    ///
    /// # Errors
    /// Any conjugate thermal, native momentum or cancellation refusal publishes
    /// no candidate and changes none of the supplied fluid/wall inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConjugateThermalVolumeMomentumStep, ChtError> {
        let accepted = self.implicit_conjugate_then_advance(
            domain,
            conductivity,
            baths,
            solids,
            contacts,
            input.volume,
            ConductiveThermalVolumeControls {
                volume: controls.evolution.volume,
                max_conductive_faces: controls.max_conductive_faces,
            },
            solve,
            workspace,
            gate,
            |volume| {
                let heated = ThermalVolumeMomentumInput { volume, ..input };
                match profile {
                    Some(profile) => self
                        .advance_piecewise_thermal_volume_momentum_with_pressure_profile(
                            domain,
                            faces,
                            profile,
                            conductivity.caloric_curve(),
                            heated,
                            controls.evolution,
                            gate,
                        ),
                    None => self.advance_piecewise_thermal_volume_momentum(
                        domain,
                        faces,
                        conductivity.caloric_curve(),
                        heated,
                        controls.evolution,
                        gate,
                    ),
                }
            },
            |step| &step.thermal,
        )?;
        let (evolution, conductive_enthalpy_increment_j, solid_endpoints, report) =
            accepted.into_conjugate();
        Ok(ImplicitConjugateThermalVolumeMomentumStep {
            evolution,
            conductive_enthalpy_increment_j,
            solid_endpoints,
            report,
        })
    }
}

pub(super) fn prepare(
    domain: &VoxelDomain,
    solids: &[ThermalVolumeSolidCell<'_, '_>],
    samples: &mut [Option<PhaseConductivitySample>],
    gate: &CancelGate,
) -> Result<Vec<Option<usize>>, ChtError> {
    if solids.len() > domain.cell_count() {
        return Err(refuse("finite-solid cell work budget exhausted"));
    }
    if solids.is_empty() {
        return Ok(Vec::new());
    }
    let mut at = vec![None; domain.cell_count()];
    for (index, solid) in solids.iter().enumerate() {
        poll_every(gate, index)?;
        if solid.cell >= domain.cell_count()
            || domain.is_fluid(solid.cell)
            || at[solid.cell].replace(index).is_some()
        {
            return Err(refuse(
                "finite thermal walls need unique fixed-solid voxels",
            ));
        }
        finite_positive("flow.implicit_solid_mass", solid.mass_kg)?;
        finite(
            "flow.implicit_solid_authored_heat",
            solid.net_enthalpy_increment_j,
        )?;
        admit_state(solid.conductivity, solid.phase)?;
        samples[solid.cell] = Some(
            solid
                .conductivity
                .state_at_specific_enthalpy(solid.phase.specific_enthalpy_j_kg())
                .map_err(|e| refuse(e.to_string()))?,
        );
    }
    poll(gate)?;
    Ok(at)
}

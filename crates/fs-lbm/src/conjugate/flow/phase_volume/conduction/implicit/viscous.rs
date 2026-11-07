//! One original molecular kick supplies the same implicit caloric endpoint.

mod feedback;
pub use feedback::{
    ImplicitThermorheologicalThermalVolumeMomentumStep,
    ImplicitThermoviscousThermalVolumeMomentumStep, ThermalViscosityFeedbackControls,
    ThermalViscosityFeedbackReport,
};

use super::*;
use crate::conjugate::flow::material_momentum::PreparedMaterialViscosity;
use crate::conjugate::flow::{
    MaterialMomentumViscosityInput, MaterialViscosityControls, StaggeredBoundaryPressure,
    ThermalVolumeMomentumInput, ViscousThermalVolumeMomentumStep,
};

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConductiveViscousThermalVolumeMomentumStep {
    pub evolution: ViscousThermalVolumeMomentumStep,
    pub conductive_enthalpy_increment_j: Vec<f64>,
    /// The shared thermal defect also includes molecular heat-addition roundoff.
    pub report: ImplicitConductiveThermalVolumeReport,
}

/// One molecular momentum kick with fluid and finite-wall caloric endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConjugateViscousThermalVolumeMomentumStep {
    pub evolution: ViscousThermalVolumeMomentumStep,
    /// Conductive heat supplied to fluid transport; solid entries stay zero.
    pub conductive_enthalpy_increment_j: Vec<f64>,
    pub solid_endpoints: Vec<ThermalVolumeSolidEndpoint>,
    pub report: ImplicitConjugateThermalVolumeReport,
}

impl FlowField {
    /// Prepare molecular momentum and physical heat once, solve that heat with
    /// source conduction and baths, then advance original volume and momentum.
    ///
    /// Viscosity retains its original frozen-phase sampling and actual M/V
    /// inertia. Authored and molecular heat enter the same implicit k(T,f)/latent
    /// endpoint; its accepted heat drives native flow using the prepared cell
    /// momentum. Pressure work and numerical viscous losses stay separate.
    /// This is the existing first-order split, without a thermal-viscosity
    /// constitutive iteration. The optional reservoir retains the original
    /// opposed mechanical mouths and pressure profile.
    ///
    /// # Errors
    /// Original source, viscosity, implicit heat, pressure, material, work-budget
    /// and cancellation refusals return no endpoint or changed input inventory.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conductive_thermal_volume_momentum_with_viscosity_and_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConductiveViscousThermalVolumeMomentumStep, ChtError> {
        let accepted = self
            .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_viscosity_and_baths(
                domain,
                reservoir,
                conductivity,
                baths.into(),
                &[],
                &[],
                input,
                viscosity,
                controls,
                viscosity_controls,
                solve,
                workspace,
                gate,
            )?;
        Ok(ImplicitConductiveViscousThermalVolumeMomentumStep {
            evolution: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.conductive_enthalpy_increment_j,
            report: accepted.report.conduction,
        })
    }

    /// Conduct independent shear/bulk heat in the same implicit caloric
    /// endpoint, then advance the prepared original material momentum once.
    ///
    /// Bulk coefficients [Pa s] have one slot per voxel, with fixed-solid slots
    /// ignored, and must be evaluated at the original shear source's state.
    /// The frozen predictor uses `2 mu dev(D) + zeta div(u) I`; its physical
    /// heat participates in every endpoint k(T,f) and bath-flux evaluation.
    /// This remains a first-order split without thermal-viscosity iteration.
    ///
    /// # Errors
    /// Invalid bulk response or any original molecular, implicit thermal,
    /// material or cancellation refusal publishes no endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conductive_thermal_volume_momentum_with_bulk_viscosity_and_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConductiveViscousThermalVolumeMomentumStep, ChtError> {
        let accepted = self
            .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_bulk_viscosity_and_baths(
                domain,
                reservoir,
                conductivity,
                baths.into(),
                &[],
                &[],
                input,
                viscosity,
                bulk_viscosity_pa_s,
                controls,
                viscosity_controls,
                solve,
                workspace,
                gate,
            )?;
        Ok(ImplicitConductiveViscousThermalVolumeMomentumStep {
            evolution: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.conductive_enthalpy_increment_j,
            report: accepted.report.conduction,
        })
    }

    /// Conduct authored and physical molecular heat into finite solid storage
    /// in the same implicit enthalpy graph, before native fluid transport.
    ///
    /// Viscosity is prepared once from the original fluid inventory. The solid
    /// owns its supplied held mass, caloric chart and k(T,f), and may cross a
    /// latent interval without entering mass or momentum transport. Each
    /// fluid/solid or solid/solid face has one equal-and-opposite heat transfer;
    /// explicit contact resistance belongs to that same native graph link.
    /// Molecular numerical loss and pressure work remain separate from heat.
    ///
    /// This retains the original frozen-viscosity, first-order heat/flow split.
    /// Workspace has one entry per fluid cell plus explicitly bound solid.
    /// All caloric, molecular-assembly and transport defects share the original
    /// thermal budget. Failure or cancellation returns neither endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_viscosity_and_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConjugateViscousThermalVolumeMomentumStep, ChtError> {
        self.advance_implicit_viscous_thermal_volume(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            viscosity,
            None,
            controls,
            viscosity_controls,
            solve,
            workspace,
            gate,
        )
    }

    /// Include independently supplied bulk viscosity in the same fluid/solid
    /// implicit enthalpy graph. Shear plus bulk heat is supplied once before
    /// evaluating fluid/solid conduction, contact resistance and bath exchange.
    ///
    /// Bulk coefficients follow the conductive bulk entry's original-state
    /// convention. Actual fluid inertia, held solid inventory, physical
    /// half-cell geometry and the shared thermal budget keep their native owners.
    ///
    /// # Errors
    /// Invalid bulk response or any original source, molecular, conjugate
    /// thermal, material or cancellation refusal returns neither endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_bulk_viscosity_and_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConjugateViscousThermalVolumeMomentumStep, ChtError> {
        self.advance_implicit_viscous_thermal_volume(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            viscosity,
            Some(bulk_viscosity_pa_s),
            controls,
            viscosity_controls,
            solve,
            workspace,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_implicit_viscous_thermal_volume(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: Option<&[f64]>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConjugateViscousThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let (faces, profile, mouth) = match reservoir {
            Some((faces, profile)) => (
                faces,
                Some(profile),
                Some(Self::admit_material_viscosity_pressure_mouth(
                    faces,
                    viscosity.boundaries,
                    input.boundary_velocity_m_s,
                )?),
            ),
            None => ([FlowFace::Wall; 6], None, None),
        };
        let (viscous, heat_budget) = self.prepare_thermal_volume_viscosity(
            domain,
            faces,
            mouth,
            conductivity.caloric_curve(),
            input,
            viscosity,
            bulk_viscosity_pa_s,
            controls.evolution,
            viscosity_controls,
            gate,
        )?;
        let (combined, molecular_assembly) = Self::assemble_thermal_volume_viscous_heat(
            domain,
            input.volume.net_enthalpy_increment_j,
            &viscous.report,
            viscosity_controls,
            heat_budget,
            gate,
        )?;
        let prepared = self.prepare_implicit_conjugate_heat(
            domain,
            conductivity,
            baths,
            solids,
            contacts,
            CoexistenceVolumeInput {
                net_enthalpy_increment_j: &combined,
                ..input.volume
            },
            ConductiveThermalVolumeControls {
                volume: controls.evolution.volume,
                max_conductive_faces: controls.max_conductive_faces,
            },
            solve,
            workspace,
            gate,
        )?;
        self.finish_implicit_viscous_thermal_volume(
            domain,
            faces,
            profile,
            conductivity,
            input,
            controls,
            prepared,
            viscous,
            molecular_assembly,
            0.0,
            heat_budget,
            gate,
        )
    }

    /// One material/pressure route consumes the accepted heat and momentum.
    #[allow(clippy::too_many_arguments)]
    fn finish_implicit_viscous_thermal_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        conductivity: &PreparedPhaseConductivity<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        prepared: PreparedImplicitHeat,
        viscous: PreparedMaterialViscosity,
        molecular_assembly: Dd,
        coupling_defect_j: f64,
        heat_budget: f64,
        gate: &CancelGate,
    ) -> Result<ImplicitConjugateViscousThermalVolumeMomentumStep, ChtError> {
        let heated = ThermalVolumeMomentumInput {
            momentum_kg_m_s: &viscous.momentum_kg_m_s,
            ..input
        };
        let mut accepted = self.finish_implicit_conjugate_heat(
            prepared,
            conductivity,
            input.volume,
            ConductiveThermalVolumeControls {
                volume: controls.evolution.volume,
                max_conductive_faces: controls.max_conductive_faces,
            },
            gate,
            |volume| {
                let input = ThermalVolumeMomentumInput { volume, ..heated };
                match profile {
                    Some(profile) => self
                        .advance_piecewise_thermal_volume_momentum_with_pressure_profile(
                            domain,
                            faces,
                            profile,
                            conductivity.caloric_curve(),
                            input,
                            controls.evolution,
                            gate,
                        ),
                    None => self.advance_piecewise_thermal_volume_momentum(
                        domain,
                        faces,
                        conductivity.caloric_curve(),
                        input,
                        controls.evolution,
                        gate,
                    ),
                }
            },
            |evolution| &evolution.thermal,
        )?;
        let heat_defect = molecular_assembly + Dd::from_f64(coupling_defect_j);
        let viscous_defect = (heat_defect
            + Dd::from_f64(
                accepted
                    .evolution
                    .thermal
                    .volume
                    .report
                    .absolute_enthalpy_defect_j,
            ))
        .to_f64();
        let total = (heat_defect + Dd::from_f64(accepted.report.absolute_enthalpy_defect_j))
            .to_f64();
        finite("flow.implicit_viscous_enthalpy_defect", viscous_defect)?;
        finite("flow.implicit_viscous_total_enthalpy_defect", total)?;
        if total > heat_budget || viscous_defect > heat_budget {
            return Err(refuse(
                "implicit conduction and molecular heat exceed shared enthalpy budget",
            ));
        }
        accepted.report.absolute_enthalpy_defect_j = total;
        poll(gate)?;
        let (evolution, increments, solid_endpoints, report) = accepted.into_conjugate();
        Ok(ImplicitConjugateViscousThermalVolumeMomentumStep {
            evolution: ViscousThermalVolumeMomentumStep {
                evolution,
                viscosity_reconstruction: viscous.reconstruction,
                viscosity: viscous.report,
                absolute_heat_assembly_defect_j: molecular_assembly.to_f64(),
                absolute_enthalpy_defect_j: viscous_defect,
            },
            conductive_enthalpy_increment_j: increments,
            solid_endpoints,
            report,
        })
    }
}

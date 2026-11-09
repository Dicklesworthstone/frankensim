//! Implicit spatial heat with both original physical dissipation owners.

mod feedback;
pub use feedback::{ImplicitThermorheologicalResistiveThermalVolumeMomentumStep, ImplicitThermoviscousResistiveThermalVolumeMomentumStep};

use super::*;
use crate::conjugate::{
    FvBoundary, MaterialViscosityControls, ViscousResistiveThermalVolumeMomentumStep,
};

/// One accepted molecular/drag/heat/volume step and its finite-wall endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConjugateViscousResistiveThermalVolumeMomentumStep {
    pub evolution: ViscousResistiveThermalVolumeMomentumStep,
    /// Fourier/bath/contact heat of the accepted resistance trial [J].
    pub conductive_enthalpy_increment_j: Vec<f64>,
    pub solid_endpoints: Vec<ThermalVolumeSolidEndpoint>,
    /// The caloric defect includes molecular heat assembly, implicit heat,
    /// resistance coupling and conservative material transport exactly once.
    pub report: ImplicitConjugateThermalVolumeReport,
    /// Cumulative heat work, including every rejected resistance trial.
    pub work: ImplicitResistanceHeatWork,
}

impl FlowField {
    /// Advance full molecular stress, implicit spatial heat and phase Darcy
    /// drag on the same original density chart and conservative M/Mh/P route.
    ///
    /// Original source-resolved viscosity supplies one native full-stress
    /// stage with actual M/V inertia. Its physical heat is added once. Every
    /// subsequent phase/resistance trial solves the SAME original fluid and
    /// finite-wall inventories with authored heat, that molecular heat and
    /// its physical Darcy heat guess. Source k(T,f), baths and contacts use
    /// the existing analytic implicit graph; rejected trials advance neither
    /// wall time nor fluid inventory. Only the accepted trial is transported.
    ///
    /// Molecular viscosity is frozen at the original state, conduction uses
    /// the implicit pre-advection endpoint and permeability uses the accepted
    /// transported phase. This is the existing first-order split and requires
    /// time refinement. Pressure/source work retains the original caller-owned
    /// net-enthalpy convention; numerical mechanical losses remain separate
    /// from physical heat. Exact solid fluid phases, moving skeletons and free
    /// surfaces remain excluded; bound fixed solids exchange heat only.
    ///
    /// `reservoir` supplies opposed pressure mouths, otherwise normal fluxes
    /// are closed. Viscous tangential walls, optional bulk stress and original
    /// source-density admission retain their native contracts. Workspace has
    /// one entry per fluid cell plus bound solid. Residual/Krylov ceilings cover
    /// every trial; caloric and both mechanical gates include all heat defects.
    /// A source, solve, budget or cancellation refusal publishes no endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conjugate_viscous_phase_resistance_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        bulk_viscosity_pa_s: Option<&[f64]>,
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ImplicitConjugateViscousResistiveThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let curve = conductivity.caloric_curve();
        if curve.identity() != law.curve().identity()
            || curve.material_card_identity() != law.curve().material_card_identity()
        {
            return Err(refuse(
                "conductivity and phase resistance need the same original chart and material",
            ));
        }
        let mut work = ImplicitResistanceHeatWork::default();
        let (evolution, (increments, solid_endpoints, mut report)) = self
            .advance_viscous_phase_resistance_thermal_volume_with_heat(
                domain,
                reservoir,
                law,
                input,
                boundaries,
                bulk_viscosity_pa_s,
                controls.resistance,
                viscosity_controls,
                gate,
                resolve,
                |heat| {
                    self.prepare_resistive_implicit_heat_trial(
                        domain,
                        conductivity,
                        baths,
                        solids,
                        contacts,
                        input.volume,
                        &heat,
                        controls,
                        solve,
                        &mut work,
                        workspace,
                        gate,
                    )
                },
            )?;
        report.conduction.absolute_enthalpy_defect_j = evolution.absolute_enthalpy_defect_j;
        let native = &evolution.evolution.evolution;
        for (cell, state) in native.thermal.volume.phase_states.iter().enumerate() {
            poll_every(gate, cell)?;
            if let Some(state) = state {
                admit_state(conductivity, *state)?;
            }
        }
        poll(gate)?;
        Ok(ImplicitConjugateViscousResistiveThermalVolumeMomentumStep {
            evolution,
            conductive_enthalpy_increment_j: increments,
            solid_endpoints,
            report,
            work,
        })
    }
}

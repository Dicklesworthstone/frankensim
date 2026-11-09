//! One molecular-stress and phase-resistance step on the original density chart.

use super::super::resistance::{
    ResistiveThermalVolumeMomentumControls, ResolvedResistiveThermalVolumeMomentumStep,
    admit_mobile_phase, source_resistance,
};
use super::*;
use crate::conjugate::FvBoundary;
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::transport::ResolvedPhaseTransportState;

/// Atomic molecular diffusion, mushy drag, physical heat and density evolution.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscousResistiveThermalVolumeMomentumStep {
    /// Accepted phase-dependent resistance, volume flow and conservative M/Mh/P.
    pub evolution: ResolvedResistiveThermalVolumeMomentumStep,
    /// Original source snapshots that supplied frozen molecular shear viscosity.
    /// The resistance endpoint sources remain in `evolution`; they may differ.
    pub viscosity_sources: Vec<Option<ResolvedPhaseTransportState>>,
    pub viscosity_reconstruction: DensityMomentumReconstructionReport,
    /// Physical molecular heat and separate numerical time/scatter losses.
    pub viscosity: MaterialViscosityReport,
    /// Roundoff from adding molecular heat to the authored caloric increments.
    pub absolute_heat_assembly_defect_j: f64,
    /// Assembly, caloric transport and accepted resistance coupling defects.
    pub absolute_enthalpy_defect_j: f64,
    /// The above caloric defects plus both native physical-work defects.
    /// Numerical dissipation is accounted separately, never included as heat.
    pub absolute_energy_defect_j: f64,
}

impl FlowField {
    /// Advance molecular stress AND phase-dependent Darcy resistance with their
    /// physical heat on the same density-changing material timeline.
    ///
    /// Original source-resolved mobile phases supply shear viscosity once.
    /// The existing full-stress owner kicks original M/P with actual M/V
    /// inertia. Its physical heat joins the authored NET caloric increments
    /// once, and its returned cell momentum enters the existing coupled
    /// phase/resistance/heat/volume owner. Every resistance trial starts from
    /// this same once-prepared state. The accepted native pressure and material
    /// routes then carry M, M*h and P once; neither trial heat nor impulse is
    /// accumulated across iterations.
    ///
    /// `resolve` must use the original law and current chart density, such as
    /// `resolve_phase_transport_state_for_volume_change`. The same checked
    /// initial sources supply both molecular viscosity and initial drag.
    /// Resistance is closed against its transported endpoint phase; molecular
    /// viscosity is frozen at the original state. This is a first-order split,
    /// not simultaneous endpoint shear/temperature/drag coupling. Refinement
    /// remains necessary. Resolve must not mutate external physical history.
    ///
    /// `reservoir` selects the existing opposed normal-traction mouths;
    /// otherwise all normal boundary fluxes are closed. Tangential moving
    /// walls and their work retain the native viscous owner. Optional bulk
    /// viscosity gives 2 mu dev(D) + zeta div(u) I; None retains 2 mu D.
    /// Initial, incoming, heated and transported phases must remain mobile.
    /// Exact solid locking, moving skeletons and free surfaces remain excluded.
    ///
    /// Each physical heat enters h once. Backward-Euler, reconstruction,
    /// scatter, pressure and donor losses stay numerical. Source-pressure and
    /// reservoir work keep the original caller-owned net-enthalpy convention.
    /// Caloric and both mechanical defects share the supplied energy budgets.
    /// No source, budget, solve or cancellation failure publishes a partial
    /// momentum, heated phase or transported inventory.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_viscous_phase_resistance_thermal_volume_momentum(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        bulk_viscosity_pa_s: Option<&[f64]>,
        controls: ResistiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ViscousResistiveThermalVolumeMomentumStep, ChtError> {
        let heat_budget = admit_viscous_thermal_input(domain, input, controls.thermal, gate)?;
        let (faces, profile, mouth) = match reservoir {
            Some((faces, profile)) => (
                faces,
                Some(profile),
                Some(Self::admit_material_viscosity_pressure_mouth(
                    faces,
                    boundaries,
                    input.boundary_velocity_m_s,
                )?),
            ),
            None => ([FlowFace::Wall; 6], None, None),
        };
        for phase in input.volume.boundary_inflow.into_iter().flatten() {
            admit_mobile_phase(law, &phase)?;
        }
        let mut original_sources = Vec::with_capacity(domain.cell_count());
        for cell in 0..domain.cell_count() {
            super::super::super::poll_every(gate, cell)?;
            if !domain.is_fluid(cell) {
                if input.volume.phase_states[cell].is_some() {
                    return Err(refuse(
                        "fixed solids cannot contain viscous phase inventory",
                    ));
                }
                original_sources.push(None);
                continue;
            }
            let phase = input.volume.phase_states[cell]
                .ok_or_else(|| refuse("molecular resistance needs the original mobile phase"))?;
            // Exact solids refuse before invoking a liquid-property resolver.
            let mobile = admit_mobile_phase(law, &phase)?;
            let source = resolve(cell, &phase)?;
            poll(gate)?;
            source_resistance(mobile, &source)?;
            original_sources.push(Some(source));
        }
        let borrowed_sources: Vec<Option<&dyn FluidTransport>> = original_sources
            .iter()
            .map(|source| source.as_ref().map(|s| s as &dyn FluidTransport))
            .collect();
        let (viscous, _) = self.prepare_thermal_volume_viscosity(
            domain,
            faces,
            mouth,
            law.curve(),
            input,
            MaterialMomentumViscosityInput {
                sources: &borrowed_sources,
                boundaries,
            },
            bulk_viscosity_pa_s,
            controls.thermal,
            viscosity_controls,
            gate,
        )?;
        let (combined, assembly) = Self::assemble_thermal_volume_viscous_heat(
            domain,
            input.volume.net_enthalpy_increment_j,
            &viscous.report,
            viscosity_controls,
            heat_budget,
            gate,
        )?;
        let heated = ThermalVolumeMomentumInput {
            volume: CoexistenceVolumeInput {
                net_enthalpy_increment_j: &combined,
                ..input.volume
            },
            momentum_kg_m_s: &viscous.momentum_kg_m_s,
            ..input
        };
        let mut source_at =
            |cell: usize, phase: &EquilibriumPhaseState| match &original_sources[cell] {
                Some(source) if source.phase_state() == phase => Ok(source.clone()),
                _ => resolve(cell, phase),
            };
        let evolution = match profile {
            Some(profile) => self
                .advance_resolved_phase_resistance_thermal_volume_momentum_with_pressure_profile(
                    domain,
                    faces,
                    profile,
                    law,
                    heated,
                    controls,
                    gate,
                    &mut source_at,
                )?,
            None => self.advance_resolved_phase_resistance_thermal_volume_momentum(
                domain,
                faces,
                law,
                heated,
                controls,
                gate,
                &mut source_at,
            )?,
        };
        let native = &evolution.evolution;
        let caloric = assembly
            + Dd::from_f64(native.thermal.volume.report.absolute_enthalpy_defect_j)
            + Dd::from_f64(native.coupling.absolute_heat_coupling_defect_j)
            + Dd::from_f64(native.coupling.absolute_constitutive_heat_defect_j);
        let energy = caloric
            + Dd::from_f64(viscous.report.absolute_energy_defect_j)
            + Dd::from_f64(native.resistance.absolute_energy_defect_j);
        let absolute_enthalpy_defect_j = caloric.to_f64();
        let absolute_energy_defect_j = energy.to_f64();
        finite(
            "flow.viscous_resistance_enthalpy_defect",
            absolute_enthalpy_defect_j,
        )?;
        finite(
            "flow.viscous_resistance_energy_defect",
            absolute_energy_defect_j,
        )?;
        if absolute_enthalpy_defect_j > heat_budget
            || absolute_energy_defect_j > viscosity_controls.maximum_energy_defect_j
            || absolute_energy_defect_j > controls.thermal.mechanics.forcing.maximum_energy_defect_j
        {
            return Err(refuse(
                "molecular and phase-resistance heat exceed their shared energy budgets",
            ));
        }
        poll(gate)?;
        Ok(ViscousResistiveThermalVolumeMomentumStep {
            evolution,
            viscosity_sources: original_sources,
            viscosity_reconstruction: viscous.reconstruction,
            viscosity: viscous.report,
            absolute_heat_assembly_defect_j: assembly.to_f64(),
            absolute_enthalpy_defect_j,
            absolute_energy_defect_j,
        })
    }
}

//! Thermal molecular response inside the original mushy/Darcy heat timeline.

mod rheology;
pub use rheology::ImplicitThermorheologicalResistiveThermalVolumeMomentumStep;

use super::*;
use std::cell::RefCell;

use super::super::super::viscous::feedback::{
    PreparedThermalFeedback, ThermalViscositySources, feedback_phase_at, resolve_thermal_sources,
};
use crate::conjugate::flow::phase_volume::momentum::admit_viscous_thermal_input;

/// Accepted thermal viscosity, implicit heat and transported phase resistance.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitThermoviscousResistiveThermalVolumeMomentumStep {
    pub evolution: ImplicitConjugateViscousResistiveThermalVolumeMomentumStep,
    /// Independent bulk response at the exact native pretransport caloric state.
    /// Its shear/phase source receipts are in `evolution.evolution.viscosity_sources`.
    pub bulk_viscosity_pa_s: Vec<f64>,
    /// Physical molecular heat supplied to the accepted native heat graph [J].
    /// Its non-cancelling difference from actual molecular heat stays explicit.
    pub used_molecular_heat_j: Vec<f64>,
    /// Includes every discarded Darcy/thermal trial and molecular predictor.
    pub feedback: ThermalViscosityFeedbackReport,
}

impl FlowField {
    /// Close thermal shear/bulk viscosity and physical molecular heat inside
    /// the same original-inventory implicit heat / phase-resistance timeline.
    ///
    /// For each native Darcy heat guess, the existing thermal-viscosity fixed
    /// point solves original fluid M/h and finite-solid histories. Exact native
    /// pretransport caloric samples resolve shear and independent bulk response;
    /// their native stress kick always uses ORIGINAL M/P and M/V inertia.
    /// Source density must match that thermal phase, and is never substituted
    /// for the original mass in the kick. The law must admit mobile mushy or
    /// liquid material with the same chart and resolved permeability evidence.
    ///
    /// That converged molecular momentum feeds the existing Darcy/pressure/
    /// volume owner. Darcy sources close at the transported endpoint; molecular
    /// sources close before advection. This remains an explicit first-order
    /// mechanical split, with no claim of a simultaneous unsplit endpoint.
    /// Both physical heats enter the same original thermal interval once;
    /// numerical losses, pressure work and boundary exchange remain separate.
    ///
    /// Feedback maximum_trials, solve residual/Krylov limits and molecular
    /// max_sweeps are WHOLE-CALL ceilings across all nested rejected trials.
    /// Repeated coefficients within a thermal solve reuse their original kick.
    /// The resolver owns fixed composition/pressure context and thermal laws,
    /// not dependence on unknown endpoint rate. Exact solids require their
    /// existing immobile owner; finite held walls retain their own histories.
    ///
    /// # Errors
    /// Source, phase-support, heat, work, mechanical or cancellation refusal
    /// returns neither fluid nor finite-solid endpoints. Workspace stays scratch.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_isotropic_implicit_conjugate_thermoviscous_phase_resistance_thermal_volume_momentum_with_baths(
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
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<ImplicitThermoviscousResistiveThermalVolumeMomentumStep, ChtError> {
        // The thermal-only entry retains one sequential source resolver for
        // both molecular stress and its declared phase-resistance response.
        let resolver = RefCell::new(resolve);
        self.advance_implicit_resistive_thermal_feedback(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            law,
            input,
            boundaries,
            controls,
            viscosity_controls,
            solve,
            feedback,
            workspace,
            gate,
            |cell, phase| (resolver.borrow_mut())(cell, phase).map(|(source, _)| source),
            |faces,
             mouth,
             heat,
             previous: Option<PreparedThermalFeedback<(ThermalViscositySources, ())>>,
             work| {
                let sources = resolve_thermal_sources(
                    domain,
                    law.curve(),
                    Some(law),
                    |cell| feedback_phase_at(law.curve(), input, heat, cell),
                    &mut *resolver.borrow_mut(),
                    work,
                    gate,
                )?;
                let viscosity = match previous {
                    Some(previous)
                        if previous.sources.0.shear_pa_s == sources.shear_pa_s
                            && previous.sources.0.bulk_pa_s == sources.bulk_pa_s =>
                    {
                        previous.viscosity
                    }
                    _ => self.prepare_thermoviscous_momentum(
                        domain,
                        faces,
                        mouth,
                        input,
                        boundaries,
                        ConductiveThermalVolumeMomentumControls {
                            evolution: controls.resistance.thermal,
                            max_conductive_faces: controls.max_conductive_faces,
                        },
                        viscosity_controls,
                        &sources,
                        work,
                        gate,
                    )?,
                };
                Ok(PreparedThermalFeedback {
                    viscosity,
                    sources: (sources, ()),
                })
            },
        )
        .map(|(evolution, ())| evolution)
    }

    /// Shared original-history molecular/heat/Darcy owner. A native molecular
    /// preparation returns exact sources plus an opaque accepted receipt; only
    /// the receipt belonging to the accepted thermal and Darcy trial escapes.
    /// The independent resistance callback retains its transported-phase owner.
    #[allow(
        clippy::too_many_arguments,
        clippy::too_many_lines,
        clippy::type_complexity
    )]
    fn advance_implicit_resistive_thermal_feedback<S>(
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
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        mut resolve_resistance: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
        mut prepare_molecular: impl FnMut(
            [FlowFace; 6],
            Option<usize>,
            Option<&PreparedImplicitHeat>,
            Option<PreparedThermalFeedback<(ThermalViscositySources, S)>>,
            &mut ThermalViscosityFeedbackReport,
        ) -> Result<
            PreparedThermalFeedback<(ThermalViscositySources, S)>,
            ChtError,
        >,
    ) -> Result<(ImplicitThermoviscousResistiveThermalVolumeMomentumStep, S), ChtError> {
        let heat_budget =
            admit_viscous_thermal_input(domain, input, controls.resistance.thermal, gate)?;
        if conductivity.caloric_curve().identity() != law.curve().identity()
            || conductivity.caloric_curve().material_card_identity()
                != law.curve().material_card_identity()
        {
            return Err(refuse(
                "thermal viscosity and Darcy need the same caloric chart",
            ));
        }
        let (faces, profile, mouth) = match reservoir {
            Some((faces, profile)) => (
                faces,
                Some(BoundaryPressure::Profile(profile)),
                Some(Self::admit_material_viscosity_pressure_mouth(
                    faces,
                    boundaries,
                    input.boundary_velocity_m_s,
                )?),
            ),
            None => ([FlowFace::Wall; 6], None, None),
        };
        let thermal_controls = ConductiveThermalVolumeMomentumControls {
            evolution: controls.resistance.thermal,
            max_conductive_faces: controls.max_conductive_faces,
        };
        let mut work = ThermalViscosityFeedbackReport::default();
        let mut darcy_source_evaluations = 0_usize;
        let (evolution, receipt) = self
            .advance_resolved_phase_resistance_thermal_volume_with_heat(
                domain,
                faces,
                profile,
                law,
                input,
                controls.resistance,
                gate,
                |cell, phase| {
                    darcy_source_evaluations = darcy_source_evaluations
                        .checked_add(1)
                        .ok_or_else(|| refuse("thermal viscosity source work overflow"))?;
                    resolve_resistance(cell, phase)
                },
                |heat| {
                    let trial_input = ThermalVolumeMomentumInput {
                        volume: CoexistenceVolumeInput {
                            net_enthalpy_increment_j: &heat,
                            ..input.volume
                        },
                        ..input
                    };
                    let prepared = self.prepare_implicit_thermal_feedback(
                        domain,
                        faces,
                        mouth,
                        conductivity,
                        baths,
                        solids,
                        contacts,
                        trial_input,
                        thermal_controls,
                        viscosity_controls,
                        solve,
                        feedback,
                        workspace,
                        gate,
                        &mut work,
                        &mut prepare_molecular,
                    )?;
                    let molecular_heat_defect =
                        prepared.assembly + Dd::from_f64(prepared.coupling_defect_j);
                    let caloric_defect = molecular_heat_defect
                        + Dd::from_f64(prepared.heat.report.absolute_enthalpy_defect_j);
                    let (combined, increments, solid_endpoints, report) =
                        prepared.heat.into_conjugate();
                    Ok(ResistanceHeatTrial {
                        increment_j: combined,
                        absolute_enthalpy_defect_j: caloric_defect.to_f64(),
                        momentum_kg_m_s: Some(
                            prepared.mechanical.viscosity.momentum_kg_m_s.clone(),
                        ),
                        receipt: (
                            prepared.mechanical,
                            prepared.used_molecular_heat_j,
                            prepared.assembly,
                            prepared.coupling_defect_j,
                            increments,
                            solid_endpoints,
                            report,
                        ),
                    })
                },
            )?;
        let (
            mechanical,
            used_molecular_heat_j,
            assembly,
            coupling_defect,
            increments,
            solid_endpoints,
            mut report,
        ) = receipt;
        let transformed_defect = (Dd::from_f64(coupling_defect)
            + Dd::from_f64(report.conduction.absolute_enthalpy_defect_j))
        .to_f64();
        let (sources, native_receipt) = mechanical.sources;
        let evolution = Self::finish_viscous_phase_resistance_thermal_volume(
            evolution,
            sources.original,
            mechanical.viscosity,
            assembly,
            transformed_defect,
            controls.resistance,
            viscosity_controls,
            heat_budget,
            gate,
        )?;
        // The accepted conservative transport endpoint must remain in the
        // original conductivity support too; no advection hides source escape.
        for (cell, phase) in evolution
            .evolution
            .evolution
            .thermal
            .volume
            .phase_states
            .iter()
            .enumerate()
        {
            poll_every(gate, cell)?;
            if let Some(phase) = phase {
                admit_state(conductivity, *phase)?;
            }
        }
        report.conduction.absolute_enthalpy_defect_j = evolution.absolute_enthalpy_defect_j;
        work.absolute_enthalpy_defect_j = evolution.absolute_enthalpy_defect_j;
        work.source_evaluations = work
            .source_evaluations
            .checked_add(darcy_source_evaluations)
            .ok_or_else(|| refuse("thermal viscosity source work overflow"))?;
        let heat_work = ImplicitResistanceHeatWork {
            heat_trials: work.trials,
            nonlinear_iterations: work.heat_nonlinear_iterations,
            linear_iterations: work.heat_linear_iterations,
            residual_evaluations: work.heat_residual_evaluations,
        };
        poll(gate)?;
        Ok((
            ImplicitThermoviscousResistiveThermalVolumeMomentumStep {
                evolution: ImplicitConjugateViscousResistiveThermalVolumeMomentumStep {
                    evolution,
                    conductive_enthalpy_increment_j: increments,
                    solid_endpoints,
                    report,
                    work: heat_work,
                },
                bulk_viscosity_pa_s: sources.bulk_pa_s,
                used_molecular_heat_j,
                feedback: work,
            },
            native_receipt,
        ))
    }
}

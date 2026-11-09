//! Bounded thermal feedback through physical molecular heat on original history.

mod rheology;
pub use rheology::ImplicitThermorheologicalThermalVolumeMomentumStep;

use super::*;
use super::super::sgs::{self, SgsThermalField};
use crate::conjugate::FvBoundary;
use crate::conjugate::flow::phase_volume::momentum::admit_viscous_thermal_input;
use fs_material::fluid_transport::viscous::NewtonianViscosity;
use fs_material::phase::transport::ResolvedPhaseTransportState;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

/// Fixed-point work and physical heat acceptance, without hidden time steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalViscosityFeedbackControls {
    /// Maximum native heat trials, all starting from the original inventory.
    /// This spans all outer Darcy trials when composed with phase resistance.
    pub maximum_trials: usize,
    /// Convex relaxation of each physical heat guess and, when selected, the
    /// frozen SGS diffusivity field, strictly in (0, 1].
    pub relaxation: f64,
    /// Sum of absolute cellwise differences between used and actual molecular
    /// heat [J], plus the non-cancelling SGS transport mismatch when selected.
    /// This also consumes the original shared thermal/work budgets.
    pub maximum_heat_defect_j: f64,
}

impl Default for ThermalViscosityFeedbackControls {
    fn default() -> Self {
        Self {
            maximum_trials: 64,
            relaxation: 0.5,
            maximum_heat_defect_j: 1e-9,
        }
    }
}

/// Whole-call work and the accepted physical coupling residual.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ThermalViscosityFeedbackReport {
    pub trials: usize,
    /// Includes the original-state predictor. Temperature-only calls reuse
    /// unchanged coefficients; endpoint rheology resolves each thermal trial.
    pub molecular_solves: usize,
    pub molecular_sweeps: usize,
    pub molecular_krylov_iterations: usize,
    pub heat_nonlinear_iterations: usize,
    pub heat_linear_iterations: usize,
    pub heat_residual_evaluations: usize,
    /// All source-resolver calls, including Darcy endpoint queries when composed.
    pub source_evaluations: usize,
    /// Signed sum Q_used - Q_actual [J], before material advection.
    pub heat_coupling_defect_j: f64,
    /// Non-cancelling molecular sum |Q_used - Q_actual| [J].
    pub absolute_heat_coupling_defect_j: f64,
    /// Twice the sum of absolute fluid-link heat changes caused by replacing
    /// the used SGS field with the accepted field at the SAME thermal endpoint.
    /// Both adjacent cell residuals are charged; opposite faces cannot cancel.
    pub absolute_sgs_heat_coupling_defect_j: f64,
    /// Explicit modeled heat transport parameter; None on molecular-only paths.
    pub turbulent_prandtl_number: Option<f64>,
    /// Cell closure actually supplied to the accepted heat graph [m^2/s].
    /// Empty unless modeled heat transport is selected; fixed solids are zero.
    pub used_cell_eddy_viscosity_m2_s: Vec<f64>,
    /// Accepted native momentum closure at that thermal endpoint [m^2/s].
    /// Its difference from the used field is accounted above, never hidden.
    pub cell_eddy_viscosity_m2_s: Vec<f64>,
    /// Accepted coupling, heat assembly, implicit endpoint and transport defects.
    pub absolute_enthalpy_defect_j: f64,
}

/// One accepted thermal constitutive endpoint and its native material advance.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitThermoviscousThermalVolumeMomentumStep {
    pub evolution: ImplicitConjugateViscousThermalVolumeMomentumStep,
    /// Exact original source receipts at the native thermal endpoint BEFORE
    /// advection; None in fixed solids. These are not transported-endpoint laws.
    pub endpoint_sources: Vec<Option<ResolvedPhaseTransportState>>,
    /// Independent bulk response evaluated with those same source states [Pa s].
    /// Fixed-solid entries are zero and never consumed by molecular stress.
    pub bulk_viscosity_pa_s: Vec<f64>,
    /// Molecular heat actually supplied to the accepted implicit heat graph.
    /// The native viscosity report retains Q_actual; their admitted difference
    /// is explicit in feedback and is never silently inserted a second time.
    pub used_molecular_heat_j: Vec<f64>,
    pub feedback: ThermalViscosityFeedbackReport,
}

pub(in super::super) struct ThermalViscositySources {
    pub(in super::super) original: Vec<Option<ResolvedPhaseTransportState>>,
    pub(in super::super) shear_pa_s: Vec<f64>,
    pub(in super::super) bulk_pa_s: Vec<f64>,
}

/// The heat driver owns the molecular state and an opaque constitutive receipt.
pub(in super::super) struct PreparedThermalFeedback<S> {
    pub(in super::super) viscosity: PreparedMaterialViscosity,
    pub(in super::super) sources: S,
}

/// Converged native heat and original-inertia molecular momentum, before any
/// pressure, Darcy or material advance. The accepted source phase belongs to
/// this thermal endpoint; it never replaces the original M/V inertia.
pub(in super::super) struct ConvergedThermalFeedback<S> {
    pub(in super::super) mechanical: PreparedThermalFeedback<S>,
    pub(in super::super) heat: PreparedImplicitHeat,
    pub(in super::super) used_molecular_heat_j: Vec<f64>,
    pub(in super::super) assembly: Dd,
    pub(in super::super) coupling_defect_j: f64,
    pub(in super::super) heat_budget_j: f64,
}

struct AcceptedThermalFeedback<S> {
    evolution: ImplicitConjugateViscousThermalVolumeMomentumStep,
    sources: S,
    used_molecular_heat_j: Vec<f64>,
    feedback: ThermalViscosityFeedbackReport,
}

impl FlowField {
    /// Resolve thermal shear and independent bulk viscosity at the same native
    /// implicit fluid/solid caloric endpoint that receives their physical heat.
    ///
    /// Each trial solves ORIGINAL fluid M/h and held-solid histories with
    /// authored heat plus a guessed molecular Q. The resolver receives the
    /// exact resulting pretransport phase state and returns its original
    /// source receipt plus independent bulk viscosity [Pa s]. The source must
    /// match that state and chart density exactly, with fully liquid mobility
    /// and no unrepresented Darcy resistance. Its actual shear viscosity and
    /// supplied bulk response drive one native kick of ORIGINAL M/P with
    /// actual M/V inertia and stress 2 mu dev(D) + zeta div(u) I.
    ///
    /// The non-cancelling mismatch between that kick's physical Q and the Q
    /// used by the graph must fit both the explicit feedback tolerance and the
    /// shared thermal/viscosity energy gates. Rejected trials relax only Q;
    /// no trial advances material or solid history again. The accepted graph
    /// keeps its used Q, its exact endpoint sources and its explicit mismatch.
    /// Force/pressure and conservative material advection run once afterward.
    ///
    /// The resolver owns fixed composition/pressure context and a thermal
    /// viscosity law, not a law of unknown endpoint velocity or shear rate.
    /// This is a bounded heat fixed point, not a simultaneous nonlinear
    /// rheology solve. Original, incoming, heated and transported fluid states
    /// must remain fully liquid; finite solids retain their own caloric owner.
    ///
    /// maximum residual evaluations and linear iterations in solve, and
    /// max_sweeps in viscosity_controls, are WHOLE-CALL ceilings here.
    /// Per-trial Newton/backtrack and native Krylov bounds remain unchanged.
    /// Exact repeated coefficient fields reuse their original molecular kick.
    ///
    /// # Errors
    /// Source/support, work exhaustion, nonconvergence, conservation or
    /// cancellation refusal returns neither fluid nor finite-solid endpoint.
    /// Inputs remain untouched; workspace is scratch for a fresh retry.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_isotropic_implicit_conjugate_thermoviscous_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    ) -> Result<ImplicitThermoviscousThermalVolumeMomentumStep, ChtError> {
        let curve = conductivity.caloric_curve();
        let accepted = self.advance_implicit_thermal_feedback(
            domain,
            reservoir,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            boundaries,
            None,
            controls,
            viscosity_controls,
            solve,
            feedback,
            workspace,
            gate,
            |faces,
             mouth,
             heat,
             previous: Option<PreparedThermalFeedback<ThermalViscositySources>>,
             work| {
                let sources = resolve_thermal_sources(
                    domain,
                    curve,
                    None,
                    |cell| feedback_phase_at(curve, input, heat, cell),
                    &mut resolve,
                    work,
                    gate,
                )?;
                let viscosity = match previous {
                    Some(previous)
                        if previous.sources.shear_pa_s == sources.shear_pa_s
                            && previous.sources.bulk_pa_s == sources.bulk_pa_s =>
                    {
                        previous.viscosity
                    }
                    _ => self.prepare_thermoviscous_momentum(
                        domain,
                        faces,
                        mouth,
                        input,
                        boundaries,
                        controls,
                        viscosity_controls,
                        &sources,
                        work,
                        gate,
                    )?,
                };
                Ok(PreparedThermalFeedback { viscosity, sources })
            },
        )?;
        Ok(ImplicitThermoviscousThermalVolumeMomentumStep {
            evolution: accepted.evolution,
            endpoint_sources: accepted.sources.original,
            bulk_viscosity_pa_s: accepted.sources.bulk_pa_s,
            used_molecular_heat_j: accepted.used_molecular_heat_j,
            feedback: accepted.feedback,
        })
    }

    /// One thermal fixed-point owner for frozen-rate and endpoint-rate laws.
    /// A constitutive callback prepares ORIGINAL momentum at each thermal trial;
    /// all accepted heat, finite-solid and material-flow work stays here.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_implicit_thermal_feedback<S>(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        turbulent_prandtl_number: Option<f64>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        mut prepare: impl FnMut(
            [FlowFace; 6],
            Option<usize>,
            Option<&PreparedImplicitHeat>,
            Option<PreparedThermalFeedback<S>>,
            &mut ThermalViscosityFeedbackReport,
        ) -> Result<PreparedThermalFeedback<S>, ChtError>,
    ) -> Result<AcceptedThermalFeedback<S>, ChtError> {
        let curve = conductivity.caloric_curve();
        for phase in input.volume.boundary_inflow.iter().flatten() {
            admit_liquid_phase(curve, *phase)?;
        }
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
        let mut work = ThermalViscosityFeedbackReport::default();
        let accepted = self.prepare_implicit_thermal_feedback_with_sgs(
            domain,
            faces,
            mouth,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            turbulent_prandtl_number,
            controls,
            viscosity_controls,
            solve,
            feedback,
            workspace,
            gate,
            &mut work,
            &mut prepare,
        )?;
        // The conservative heat representation may differ from the native
        // endpoint by its admitted residual. It too must stay liquid before
        // advection can mix or conceal that departure.
        for cell in 0..domain.cell_count() {
            poll_every(gate, cell)?;
            if domain.is_fluid(cell) {
                let original = input.volume.phase_states[cell]
                    .ok_or_else(|| refuse("missing original liquid phase"))?;
                let represented = curve
                    .advance_specific_energy(
                        original,
                        accepted.heat.combined[cell] / input.volume.mass_kg[cell],
                    )
                    .map_err(|error| refuse(error.to_string()))?;
                admit_liquid_phase(curve, represented)?;
            }
        }
        let PreparedThermalFeedback { viscosity, sources } = accepted.mechanical;
        let evolution = self.finish_implicit_viscous_thermal_volume(
            domain,
            faces,
            profile,
            conductivity,
            input,
            controls,
            accepted.heat,
            viscosity,
            accepted.assembly,
            accepted.coupling_defect_j,
            accepted.heat_budget_j,
            gate,
        )?;
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
                admit_liquid_phase(curve, *phase)?;
            }
        }
        work.absolute_enthalpy_defect_j = evolution.report.conduction.absolute_enthalpy_defect_j;
        poll(gate)?;
        Ok(AcceptedThermalFeedback {
            evolution,
            sources,
            used_molecular_heat_j: accepted.used_molecular_heat_j,
            feedback: work,
        })
    }

    /// Shared native thermal/molecular preparation, retaining original M/P/h.
    /// Callers own phase-support admission and the subsequent mechanical split.
    /// The supplied work ledger spans every nested heat and molecular trial.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn prepare_implicit_thermal_feedback<S>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        mouth: Option<usize>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        work: &mut ThermalViscosityFeedbackReport,
        prepare: impl FnMut(
            [FlowFace; 6],
            Option<usize>,
            Option<&PreparedImplicitHeat>,
            Option<PreparedThermalFeedback<S>>,
            &mut ThermalViscosityFeedbackReport,
        ) -> Result<PreparedThermalFeedback<S>, ChtError>,
    ) -> Result<ConvergedThermalFeedback<S>, ChtError> {
        self.prepare_implicit_thermal_feedback_with_sgs(
            domain, faces, mouth, conductivity, baths, solids, contacts, input,
            None, controls, viscosity_controls, solve, feedback, workspace, gate,
            work, prepare,
        )
    }

    /// Close molecular heating and optional modeled fluid heat transport on one
    /// original-history thermal/mechanical fixed point.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn prepare_implicit_thermal_feedback_with_sgs<S>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        mouth: Option<usize>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: ThermalVolumeMomentumInput<'_>,
        turbulent_prandtl_number: Option<f64>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        feedback: ThermalViscosityFeedbackControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        work: &mut ThermalViscosityFeedbackReport,
        mut prepare: impl FnMut(
            [FlowFace; 6],
            Option<usize>,
            Option<&PreparedImplicitHeat>,
            Option<PreparedThermalFeedback<S>>,
            &mut ThermalViscosityFeedbackReport,
        ) -> Result<PreparedThermalFeedback<S>, ChtError>,
    ) -> Result<ConvergedThermalFeedback<S>, ChtError> {
        let heat_budget = admit_viscous_thermal_input(domain, input, controls.evolution, gate)?;
        if feedback.maximum_trials == 0
            || !feedback.relaxation.is_finite()
            || feedback.relaxation <= 0.0
            || feedback.relaxation > 1.0
            || !feedback.maximum_heat_defect_j.is_finite()
            || feedback.maximum_heat_defect_j < 0.0
            || domain.cell_count() > viscosity_controls.max_cells
            || domain.cell_count() > controls.evolution.mechanics.reconstruction.max_cells
        {
            return Err(refuse(
                "invalid bounded thermal-viscosity feedback controls",
            ));
        }
        if work.trials >= feedback.maximum_trials {
            return Err(refuse("thermal viscosity whole-call trial work exhausted"));
        }
        if let Some(prandtl) = turbulent_prandtl_number {
            sgs::admit_prandtl_number(prandtl)?;
        }
        let mut mechanical = prepare(faces, mouth, None, None, work)?;
        let mut used_heat = mechanical.viscosity.report.cell_molecular_heat_j.clone();
        let mut used_eddy = Vec::new();
        if turbulent_prandtl_number.is_some() {
            let actual = &mechanical.viscosity.cell_eddy_viscosity_m2_s;
            if actual.len() != domain.cell_count() {
                return Err(refuse("modeled heat transport requires the accepted native SGS field"));
            }
            used_eddy.try_reserve_exact(actual.len())
                .map_err(|_| refuse("SGS thermal closure allocation failed"))?;
            for (cell, &value) in actual.iter().enumerate() {
                poll_every(gate, cell)?;
                finite("flow.thermal_sgs_eddy", value)?;
                if value < 0.0 {
                    return Err(refuse("SGS thermal diffusivity must be nonnegative"));
                }
                used_eddy.push(value);
            }
        }
        while work.trials < feedback.maximum_trials {
            poll(gate)?;
            let mut remaining = solve;
            remaining.maximum_residual_evaluations = solve
                .maximum_residual_evaluations
                .saturating_sub(work.heat_residual_evaluations);
            remaining.linear.maximum_iterations = solve
                .linear
                .maximum_iterations
                .saturating_sub(work.heat_linear_iterations);
            if remaining.maximum_residual_evaluations == 0
                || remaining.linear.maximum_iterations == 0
            {
                return Err(refuse("thermal viscosity whole-call heat work exhausted"));
            }
            let (combined, assembly) =
                combine_guessed_heat(input.volume.net_enthalpy_increment_j, &used_heat, gate)?;
            let thermal_sgs = turbulent_prandtl_number.map(|prandtl| SgsThermalField {
                eddy_viscosity_m2_s: &used_eddy,
                turbulent_prandtl_number: prandtl,
            });
            let prepared = self.prepare_implicit_conjugate_heat_with_sgs(
                domain,
                conductivity,
                baths,
                solids,
                contacts,
                CoexistenceVolumeInput {
                    net_enthalpy_increment_j: &combined,
                    ..input.volume
                },
                thermal_sgs,
                ConductiveThermalVolumeControls {
                    volume: controls.evolution.volume,
                    max_conductive_faces: controls.max_conductive_faces,
                },
                remaining,
                workspace,
                gate,
            )?;
            let native = &prepared.report.solve;
            for (count, added) in [
                (&mut work.trials, 1),
                (
                    &mut work.heat_nonlinear_iterations,
                    native.nonlinear_iterations,
                ),
                (&mut work.heat_linear_iterations, native.linear_iterations),
                (
                    &mut work.heat_residual_evaluations,
                    native.residual_evaluations,
                ),
            ] {
                add_work(count, added)?;
            }
            mechanical = prepare(faces, mouth, Some(&prepared), Some(mechanical), work)?;
            let (signed, absolute) = physical_heat_mismatch(
                &used_heat,
                &mechanical.viscosity.report.cell_molecular_heat_j,
                gate,
            )?;
            let sgs_defect = match thermal_sgs {
                Some(used) => sgs::heat_flux_mismatch_j(
                    domain,
                    conductivity,
                    input.volume.mass_kg,
                    used,
                    &mechanical.viscosity.cell_eddy_viscosity_m2_s,
                    &prepared.endpoint_samples,
                    controls.evolution.mechanics.projection.duration_s,
                    controls.max_conductive_faces,
                    gate,
                )?,
                None => 0.0,
            };
            let coupling_defect = if sgs_defect == 0.0 {
                absolute
            } else {
                (Dd::from_f64(absolute) + Dd::from_f64(sgs_defect))
                    .to_f64()
                    .next_up()
            };
            finite("flow.thermoviscous_combined_heat_defect", coupling_defect)?;
            let shared = (assembly
                + Dd::from_f64(coupling_defect)
                + Dd::from_f64(prepared.report.absolute_enthalpy_defect_j))
            .to_f64();
            let molecular = (assembly
                + Dd::from_f64(absolute)
                + Dd::from_f64(mechanical.viscosity.report.absolute_energy_defect_j))
            .to_f64();
            finite("flow.thermoviscous_shared_heat_defect", shared)?;
            finite("flow.thermoviscous_molecular_heat_defect", molecular)?;
            if coupling_defect <= feedback.maximum_heat_defect_j
                && shared <= heat_budget
                && molecular <= viscosity_controls.maximum_energy_defect_j
            {
                work.heat_coupling_defect_j = signed;
                work.absolute_heat_coupling_defect_j = absolute;
                work.absolute_sgs_heat_coupling_defect_j = sgs_defect;
                if let Some(prandtl) = turbulent_prandtl_number {
                    work.turbulent_prandtl_number = Some(prandtl);
                    work.used_cell_eddy_viscosity_m2_s = used_eddy;
                    work.cell_eddy_viscosity_m2_s =
                        core::mem::take(&mut mechanical.viscosity.cell_eddy_viscosity_m2_s);
                }
                poll(gate)?;
                return Ok(ConvergedThermalFeedback {
                    mechanical,
                    heat: prepared,
                    used_molecular_heat_j: used_heat,
                    assembly,
                    coupling_defect_j: coupling_defect,
                    heat_budget_j: heat_budget,
                });
            }
            if work.trials < feedback.maximum_trials {
                for (cell, (used, actual)) in used_heat
                    .iter_mut()
                    .zip(&mechanical.viscosity.report.cell_molecular_heat_j)
                    .enumerate()
                {
                    poll_every(gate, cell)?;
                    *used = if feedback.relaxation == 1.0 {
                        *actual
                    } else {
                        (1.0 - feedback.relaxation) * *used + feedback.relaxation * *actual
                    };
                    finite("flow.thermoviscous_relaxed_heat", *used)?;
                    if *used < 0.0 {
                        return Err(refuse(
                            "thermal viscosity requires nonnegative molecular heat",
                        ));
                    }
                }
                if turbulent_prandtl_number.is_some() {
                    let actual = &mechanical.viscosity.cell_eddy_viscosity_m2_s;
                    if actual.len() != used_eddy.len() {
                        return Err(refuse("SGS thermal closure changed field shape"));
                    }
                    for (cell, (used, actual)) in used_eddy.iter_mut().zip(actual).enumerate() {
                        poll_every(gate, cell)?;
                        *used = if feedback.relaxation == 1.0 {
                            *actual
                        } else {
                            (1.0 - feedback.relaxation) * *used + feedback.relaxation * *actual
                        };
                        finite("flow.thermoviscous_relaxed_eddy", *used)?;
                        if *used < 0.0 {
                            return Err(refuse("SGS thermal diffusivity must be nonnegative"));
                        }
                    }
                }
            }
        }
        Err(refuse(
            "thermal viscosity heat feedback did not converge within its trial budget",
        ))
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn prepare_thermoviscous_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        mouth: Option<usize>,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        sources: &ThermalViscositySources,
        work: &mut ThermalViscosityFeedbackReport,
        gate: &CancelGate,
    ) -> Result<PreparedMaterialViscosity, ChtError> {
        let viscosity_controls = remaining_molecular_controls(viscosity_controls, work)?;
        let viscous = self.prepare_material_viscosity_with_source(
            domain,
            faces,
            mouth,
            input.volume.mass_kg,
            input.momentum_kg_m_s,
            boundaries,
            Some(&sources.bulk_pa_s),
            controls.evolution.mechanics.reconstruction,
            controls.evolution.mechanics.projection.duration_s,
            viscosity_controls,
            gate,
            |_| Ok(sources.shear_pa_s.clone()),
        )?;
        add_molecular_work(work, &viscous)?;
        Ok(viscous)
    }
}

pub(in super::super) fn feedback_phase_at(
    curve: &EquilibriumEnthalpyPhaseCurve,
    input: ThermalVolumeMomentumInput<'_>,
    heat: Option<&PreparedImplicitHeat>,
    cell: usize,
) -> Result<EquilibriumPhaseState, ChtError> {
    if let Some(heat) = heat {
        let sample = heat.endpoint_samples[cell]
            .ok_or_else(|| refuse("missing native fluid thermal endpoint"))?;
        curve
            .state_at_specific_enthalpy(sample.specific_enthalpy_j_kg())
            .map_err(|error| refuse(error.to_string()))
    } else {
        input.volume.phase_states[cell]
            .ok_or_else(|| refuse("thermal viscosity needs an original fluid phase"))
    }
}

fn remaining_molecular_controls(
    mut controls: MaterialViscosityControls,
    work: &ThermalViscosityFeedbackReport,
) -> Result<MaterialViscosityControls, ChtError> {
    controls.max_sweeps = controls.max_sweeps.saturating_sub(work.molecular_sweeps);
    if controls.max_sweeps == 0 {
        return Err(refuse(
            "thermal viscosity whole-call molecular work exhausted",
        ));
    }
    Ok(controls)
}

fn add_molecular_work(
    work: &mut ThermalViscosityFeedbackReport,
    viscous: &PreparedMaterialViscosity,
) -> Result<(), ChtError> {
    for (count, added) in [
        (&mut work.molecular_solves, 1),
        (&mut work.molecular_sweeps, viscous.report.sweeps),
        (
            &mut work.molecular_krylov_iterations,
            viscous.report.krylov_iterations,
        ),
    ] {
        add_work(count, added)?;
    }
    Ok(())
}

fn admit_liquid_phase(
    curve: &EquilibriumEnthalpyPhaseCurve,
    phase: EquilibriumPhaseState,
) -> Result<(), ChtError> {
    if phase.phase_curve_identity() != curve.identity()
        || phase.material_card_identity() != curve.material_card_identity()
        || phase.liquid_mass_fraction() != 1.0
    {
        return Err(refuse(
            "thermal viscosity needs the original fully liquid chart state",
        ));
    }
    Ok(())
}

pub(in super::super) fn resolve_thermal_sources(
    domain: &VoxelDomain,
    curve: &EquilibriumEnthalpyPhaseCurve,
    mobile_law: Option<&fs_material::phase::mushy::MushyPhaseLaw<'_>>,
    phase_at: impl Fn(usize) -> Result<EquilibriumPhaseState, ChtError>,
    resolve: &mut impl FnMut(
        usize,
        &EquilibriumPhaseState,
    ) -> Result<(ResolvedPhaseTransportState, f64), ChtError>,
    work: &mut ThermalViscosityFeedbackReport,
    gate: &CancelGate,
) -> Result<ThermalViscositySources, ChtError> {
    let n = domain.cell_count();
    let mut original = Vec::with_capacity(n);
    let mut shear_pa_s = vec![0.0; n];
    let mut bulk_pa_s = vec![0.0; n];
    for cell in 0..n {
        poll_every(gate, cell)?;
        if !domain.is_fluid(cell) {
            original.push(None);
            continue;
        }
        let phase = phase_at(cell)?;
        let mobile = match mobile_law {
            Some(law) => Some(
                crate::conjugate::flow::phase_volume::momentum::admit_mobile_phase(law, &phase)?,
            ),
            None => {
                admit_liquid_phase(curve, phase)?;
                None
            }
        };
        poll(gate)?;
        add_work(&mut work.source_evaluations, 1)?;
        let resolved = resolve(cell, &phase);
        poll(gate)?;
        let (source, bulk) = resolved?;
        if *source.phase_state() != phase
            || source.source_density_kg_m3().to_bits() != phase.bulk_density_kg_m3().to_bits()
        {
            return Err(refuse(
                "thermal viscosity source must match the exact endpoint phase and density",
            ));
        }
        if let Some(mobile) = mobile {
            crate::conjugate::flow::phase_volume::momentum::source_resistance(mobile, &source)?;
        } else if source.mobility().mobile().is_none_or(|state| {
            state.inverse_permeability_per_m2 != 0.0 || state.volumes.liquid != 1.0
        }) {
            return Err(refuse(
                "thermal viscosity needs liquid without unrepresented Darcy support",
            ));
        }
        let mu = source
            .mobile_transport()
            .ok_or_else(|| refuse("liquid thermal viscosity source has no mobile transport"))?
            .dynamic_viscosity_pa_s();
        let response =
            NewtonianViscosity::new(mu, bulk).map_err(|error| refuse(error.to_string()))?;
        shear_pa_s[cell] = response.shear_viscosity_pa_s();
        bulk_pa_s[cell] = response.bulk_viscosity_pa_s();
        original.push(Some(source));
    }
    Ok(ThermalViscositySources {
        original,
        shear_pa_s,
        bulk_pa_s,
    })
}

fn combine_guessed_heat(
    authored: &[f64],
    used: &[f64],
    gate: &CancelGate,
) -> Result<(Vec<f64>, Dd), ChtError> {
    let mut combined = Vec::with_capacity(authored.len());
    let mut assembly = Dd::ZERO;
    for (cell, (&external, &molecular)) in authored.iter().zip(used).enumerate() {
        poll_every(gate, cell)?;
        let (sum, residual) = accounted_add(external, molecular)?;
        assembly = assembly + Dd::from_f64(residual.abs());
        combined.push(sum);
    }
    finite("flow.thermoviscous_heat_assembly", assembly.to_f64())?;
    Ok((combined, assembly))
}

fn physical_heat_mismatch(
    used: &[f64],
    actual: &[f64],
    gate: &CancelGate,
) -> Result<(f64, f64), ChtError> {
    let mut signed = Dd::ZERO;
    let mut absolute = Dd::ZERO;
    for (cell, (&used, &actual)) in used.iter().zip(actual).enumerate() {
        poll_every(gate, cell)?;
        let mismatch = Dd::from_f64(used) - Dd::from_f64(actual);
        signed = signed + mismatch;
        absolute = absolute + mismatch.abs();
    }
    let signed = signed.to_f64();
    let absolute = absolute.to_f64();
    finite("flow.thermoviscous_signed_heat_mismatch", signed)?;
    finite("flow.thermoviscous_absolute_heat_mismatch", absolute)?;
    Ok((signed, absolute))
}

fn add_work(count: &mut usize, added: usize) -> Result<(), ChtError> {
    *count = count
        .checked_add(added)
        .ok_or_else(|| refuse("thermal viscosity work count overflow"))?;
    Ok(())
}

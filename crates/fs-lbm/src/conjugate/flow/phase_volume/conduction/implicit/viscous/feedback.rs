//! Bounded thermal feedback through physical molecular heat on original history.

use super::*;
use crate::conjugate::FvBoundary;
use crate::conjugate::flow::phase_volume::momentum::admit_viscous_thermal_input;
use fs_material::fluid_transport::viscous::NewtonianViscosity;
use fs_material::phase::transport::ResolvedPhaseTransportState;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

/// Fixed-point work and physical heat acceptance, without hidden time steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalViscosityFeedbackControls {
    /// Maximum native heat trials, all starting from the original inventory.
    pub maximum_trials: usize,
    /// Convex relaxation of each physical heat guess, strictly in (0, 1].
    pub relaxation: f64,
    /// Sum of absolute cellwise differences between used and actual molecular
    /// heat [J]. This also consumes the original shared thermal/work budgets.
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
    /// Includes the original-state predictor; unchanged coefficients reuse it.
    pub molecular_solves: usize,
    pub molecular_sweeps: usize,
    pub molecular_krylov_iterations: usize,
    pub heat_nonlinear_iterations: usize,
    pub heat_linear_iterations: usize,
    pub heat_residual_evaluations: usize,
    pub source_evaluations: usize,
    /// Signed sum Q_used - Q_actual [J], before material advection.
    pub heat_coupling_defect_j: f64,
    /// Non-cancelling sum |Q_used - Q_actual| [J], used for acceptance.
    pub absolute_heat_coupling_defect_j: f64,
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

struct ThermalViscositySources {
    original: Vec<Option<ResolvedPhaseTransportState>>,
    shear_pa_s: Vec<f64>,
    bulk_pa_s: Vec<f64>,
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
        let mut sources = resolve_thermal_sources(
            domain,
            curve,
            |cell| {
                input.volume.phase_states[cell]
                    .ok_or_else(|| refuse("thermal viscosity needs an original fluid phase"))
            },
            &mut resolve,
            &mut work,
            gate,
        )?;
        let mut viscous = self.prepare_thermoviscous_momentum(
            domain,
            faces,
            mouth,
            input,
            boundaries,
            controls,
            viscosity_controls,
            &sources,
            &mut work,
            gate,
        )?;
        let mut used_heat = viscous.report.cell_molecular_heat_j.clone();
        for _ in 0..feedback.maximum_trials {
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
            let endpoint_sources = resolve_thermal_sources(
                domain,
                curve,
                |cell| {
                    let sample = prepared.endpoint_samples[cell]
                        .ok_or_else(|| refuse("missing native fluid thermal endpoint"))?;
                    curve
                        .state_at_specific_enthalpy(sample.specific_enthalpy_j_kg())
                        .map_err(|error| refuse(error.to_string()))
                },
                &mut resolve,
                &mut work,
                gate,
            )?;
            if sources.shear_pa_s != endpoint_sources.shear_pa_s
                || sources.bulk_pa_s != endpoint_sources.bulk_pa_s
            {
                viscous = self.prepare_thermoviscous_momentum(
                    domain,
                    faces,
                    mouth,
                    input,
                    boundaries,
                    controls,
                    viscosity_controls,
                    &endpoint_sources,
                    &mut work,
                    gate,
                )?;
            }
            sources = endpoint_sources;
            let (signed, absolute) =
                physical_heat_mismatch(&used_heat, &viscous.report.cell_molecular_heat_j, gate)?;
            let shared = (assembly
                + Dd::from_f64(absolute)
                + Dd::from_f64(prepared.report.absolute_enthalpy_defect_j))
            .to_f64();
            let molecular = (assembly
                + Dd::from_f64(absolute)
                + Dd::from_f64(viscous.report.absolute_energy_defect_j))
            .to_f64();
            finite("flow.thermoviscous_shared_heat_defect", shared)?;
            finite("flow.thermoviscous_molecular_heat_defect", molecular)?;
            if absolute <= feedback.maximum_heat_defect_j
                && shared <= heat_budget
                && molecular <= viscosity_controls.maximum_energy_defect_j
            {
                // The conservative heat representation may differ from the
                // native endpoint by its admitted residual. It too must stay
                // liquid before advection can mix or conceal that departure.
                for cell in 0..domain.cell_count() {
                    poll_every(gate, cell)?;
                    if domain.is_fluid(cell) {
                        let original = input.volume.phase_states[cell]
                            .ok_or_else(|| refuse("missing original liquid phase"))?;
                        let represented = curve
                            .advance_specific_energy(
                                original,
                                prepared.combined[cell] / input.volume.mass_kg[cell],
                            )
                            .map_err(|error| refuse(error.to_string()))?;
                        admit_liquid_phase(curve, represented)?;
                    }
                }
                let evolution = self.finish_implicit_viscous_thermal_volume(
                    domain,
                    faces,
                    profile,
                    conductivity,
                    input,
                    controls,
                    prepared,
                    viscous,
                    assembly,
                    absolute,
                    heat_budget,
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
                work.heat_coupling_defect_j = signed;
                work.absolute_heat_coupling_defect_j = absolute;
                work.absolute_enthalpy_defect_j =
                    evolution.report.conduction.absolute_enthalpy_defect_j;
                poll(gate)?;
                return Ok(ImplicitThermoviscousThermalVolumeMomentumStep {
                    evolution,
                    endpoint_sources: sources.original,
                    bulk_viscosity_pa_s: sources.bulk_pa_s,
                    used_molecular_heat_j: used_heat,
                    feedback: work,
                });
            }
            if work.trials < feedback.maximum_trials {
                for (cell, (used, actual)) in used_heat
                    .iter_mut()
                    .zip(&viscous.report.cell_molecular_heat_j)
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
            }
        }
        Err(refuse(
            "thermal viscosity heat feedback did not converge within its trial budget",
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_thermoviscous_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        mouth: Option<usize>,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ConductiveThermalVolumeMomentumControls,
        mut viscosity_controls: MaterialViscosityControls,
        sources: &ThermalViscositySources,
        work: &mut ThermalViscosityFeedbackReport,
        gate: &CancelGate,
    ) -> Result<PreparedMaterialViscosity, ChtError> {
        viscosity_controls.max_sweeps = viscosity_controls
            .max_sweeps
            .saturating_sub(work.molecular_sweeps);
        if viscosity_controls.max_sweeps == 0 {
            return Err(refuse(
                "thermal viscosity whole-call molecular work exhausted",
            ));
        }
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
        Ok(viscous)
    }
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

fn resolve_thermal_sources(
    domain: &VoxelDomain,
    curve: &EquilibriumEnthalpyPhaseCurve,
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
        admit_liquid_phase(curve, phase)?;
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
        if source.mobility().mobile().is_none_or(|state| {
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

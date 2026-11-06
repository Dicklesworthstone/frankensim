//! Implicit heat inside the original phase/drag/pressure fixed-point owner.

use super::*;
use crate::conjugate::flow::density_projection::BoundaryPressure;
use crate::conjugate::flow::{
    ResolvedResistiveThermalVolumeMomentumStep, StaggeredBoundaryPressure,
    ThermalVolumeMomentumInput,
};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::transport::ResolvedPhaseTransportState;

/// Work across ALL caloric trials, including rejected resistance guesses.
/// Per-direction Krylov cycles and per-trial Newton/backtrack limits retain
/// their native meanings; residual and Krylov ceilings cover the whole call.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ImplicitResistanceHeatWork {
    pub heat_trials: usize,
    pub nonlinear_iterations: usize,
    pub linear_iterations: usize,
    pub residual_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConductiveResistiveThermalVolumeMomentumStep {
    pub evolution: ResolvedResistiveThermalVolumeMomentumStep,
    /// Endpoint Fourier/bath heat from the ACCEPTED drag/volume trial [J].
    pub conductive_enthalpy_increment_j: Vec<f64>,
    /// Accepted heat solve; the total defect also includes caloric transport
    /// and the native physical-versus-used drag/constitutive heat residuals.
    pub report: ImplicitConductiveThermalVolumeReport,
    pub work: ImplicitResistanceHeatWork,
}

impl FlowField {
    /// Close implicit conduction/baths, phase resistance and physical drag heat
    /// through the SAME original chart, pressure and material-volume iteration.
    ///
    /// Each trial solves one backward-Euler heat interval from original M/h,
    /// using authored heat plus the current physical-drag heat guess. Accepted
    /// endpoint Fourier heat then enters the original volume/pressure trial.
    /// Its returned phase resolves mu and permeability; native drag heat and
    /// constitutive impulse/heat residuals must agree before publication.
    /// Only this accepted trial scatters pressure/drag impulse and transports
    /// momentum. Rejected trial heat never accumulates in physical inventory.
    ///
    /// Source k(T,f), bath tangents and latent storage use the existing analytic
    /// heat bridge. Conduction is implicit at the pre-advection caloric state;
    /// permeability is resolved at the transported endpoint. This retains the
    /// original first-order heat/flow split and requires time refinement.
    /// Bath heat, pressure/source work and numerical mechanical losses keep
    /// separate physical accounts. Exact solids still need an immobile owner.
    ///
    /// The supplied residual-evaluation and total Krylov-iteration ceilings
    /// cover all outer trials. Other heat controls apply to each trial; the
    /// existing resistance control bounds their number. Workspace is scratch
    /// and must have one entry per fluid cell. The optional mechanical pressure
    /// profile retains its original owner and does not alter the material EOS.
    ///
    /// # Errors
    /// Any source, phase, thermal, pressure, coupling/work budget or cancellation
    /// refusal returns no endpoint and changes none of the input inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conductive_phase_resistance_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ImplicitConductiveResistiveThermalVolumeMomentumStep, ChtError> {
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
        let (evolution, (increments, mut report)) = self
            .advance_resolved_phase_resistance_thermal_volume_with_heat(
                domain,
                faces,
                profile.map(BoundaryPressure::Profile),
                law,
                input,
                controls.resistance,
                gate,
                resolve,
                |heat| {
                    let mut remaining = solve;
                    remaining.maximum_residual_evaluations = solve
                        .maximum_residual_evaluations
                        .saturating_sub(work.residual_evaluations);
                    remaining.linear.maximum_iterations = solve
                        .linear
                        .maximum_iterations
                        .saturating_sub(work.linear_iterations);
                    if remaining.maximum_residual_evaluations == 0
                        || remaining.linear.maximum_iterations == 0
                    {
                        return Err(refuse("implicit resistance whole-call heat work exhausted"));
                    }
                    let prepared = self.prepare_implicit_conductive_heat(
                        domain,
                        conductivity,
                        baths,
                        CoexistenceVolumeInput {
                            net_enthalpy_increment_j: &heat,
                            ..input.volume
                        },
                        ConductiveThermalVolumeControls {
                            volume: controls.resistance.thermal.volume,
                            max_conductive_faces: controls.max_conductive_faces,
                        },
                        remaining,
                        workspace,
                        gate,
                    )?;
                    let native = prepared.report.solve;
                    for (count, added) in [
                        (&mut work.heat_trials, 1),
                        (&mut work.nonlinear_iterations, native.nonlinear_iterations),
                        (&mut work.linear_iterations, native.linear_iterations),
                        (&mut work.residual_evaluations, native.residual_evaluations),
                    ] {
                        *count = count
                            .checked_add(added)
                            .ok_or_else(|| refuse("implicit resistance work count overflow"))?;
                    }
                    Ok((
                        prepared.combined,
                        prepared.report.absolute_enthalpy_defect_j,
                        (prepared.increments, prepared.report),
                    ))
                },
            )?;
        let native = &evolution.evolution;
        // These are the same non-cancelling terms already admitted together by
        // the native outer iteration, including its mechanics energy gate.
        let coupling = &native.coupling;
        report.absolute_enthalpy_defect_j =
            (Dd::from_f64(coupling.absolute_heat_coupling_defect_j)
                + Dd::from_f64(coupling.absolute_constitutive_heat_defect_j)
                + Dd::from_f64(native.thermal.volume.report.absolute_enthalpy_defect_j)
                + Dd::from_f64(report.absolute_enthalpy_defect_j))
            .to_f64();
        for (cell, state) in native.thermal.volume.phase_states.iter().enumerate() {
            poll_every(gate, cell)?;
            if let Some(state) = state {
                admit_state(conductivity, *state)?;
            }
        }
        poll(gate)?;
        Ok(ImplicitConductiveResistiveThermalVolumeMomentumStep {
            evolution,
            conductive_enthalpy_increment_j: increments,
            report,
            work,
        })
    }
}

#[cfg(test)]
mod tests;

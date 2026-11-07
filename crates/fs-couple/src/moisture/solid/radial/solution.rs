//! Finite aqueous activity and its dependent native cylindrical solid endpoint.

use super::*;
use crate::moisture::solution::{
    FiniteMoistureSolution, RadialSolutionMoistureCandidate, RadialSolutionSdirk2Candidate,
    SolutionMoistureError, SolutionMoistureSolve, SolutionMoistureStep,
};

#[allow(clippy::too_many_arguments)]
impl MoistureRadialSolid {
    /// Exchange bound water with the original finite solution through the
    /// existing equilibrated vapor path, then resolve compatible native shape,
    /// stress and the wet prestressed vibration pencil from that exact endpoint.
    /// A late material or mechanical refusal publishes neither inventory.
    /// Solution osmotic pressure remains a diagnostic, not an applied traction.
    pub fn advance_finite_solution(
        &self,
        reservoir: &FiniteMoistureSolution,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<(RadialSolutionMoistureCandidate, MoistureSolidState), SolutionMoistureError> {
        let candidate =
            reservoir.advance_radial(moisture, dt_s, faces, transport, nonlinear, physical, cx)?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }

    /// Advance both finite-solution radial SDIRK stages before resolving native
    /// mechanics. The final moisture field alone selects sourced swelling,
    /// elasticity and inertia; no intermediate mechanical state is published.
    /// Diffusion stays on its original reference chart. Mechanical work does
    /// not feed chemical potential, and temperature remains prescribed.
    pub fn advance_finite_solution_sdirk2(
        &self,
        reservoir: &FiniteMoistureSolution,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<(RadialSolutionSdirk2Candidate, MoistureSolidState), SolutionMoistureError> {
        let candidate = reservoir
            .advance_radial_sdirk2(moisture, dt_s, faces, transport, nonlinear, physical, cx)?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl MoistureRadialSolid {
    /// Complete adaptive finite-solution backward-Euler transport, then resolve
    /// native shape, stress and wet inertia at the accepted final moisture.
    /// Mechanics is evaluated once after the full water duration; a late source
    /// or mechanical refusal publishes neither the field nor finite bath.
    pub fn advance_finite_solution_adaptive_to(
        &self,
        reservoir: &FiniteMoistureSolution,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            FiniteMoistureSolution,
            MoistureRadialSolidStep<MoistureAdvanceReport<SolutionMoistureStep, 2>>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        let candidate = reservoir.advance_radial_adaptive_to(
            moisture, target_s, faces, transport, nonlinear, physical, time, cx,
        )?;
        self.resolve_solution_adaptive(candidate, cx)
    }

    /// Complete adaptive finite-solution SDIRK2 transport and publish its native
    /// mechanical endpoint atomically. Only the final accepted moisture selects
    /// sourced swelling, elasticity and inertia. The original reference-chart,
    /// small-strain and prescribed-temperature assumptions remain unchanged;
    /// neither chemical-potential feedback nor osmotic traction is inferred.
    pub fn advance_finite_solution_sdirk2_adaptive_to(
        &self,
        reservoir: &FiniteMoistureSolution,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            FiniteMoistureSolution,
            MoistureRadialSolidStep<
                MoistureAdvanceReport<SolutionMoistureStep<MoistureSdirk2Report>, 2>,
            >,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        let candidate = reservoir.advance_radial_sdirk2_adaptive_to(
            moisture, target_s, faces, transport, nonlinear, physical, time, cx,
        )?;
        self.resolve_solution_adaptive(candidate, cx)
    }

    fn resolve_solution_adaptive<Report>(
        &self,
        (reservoir, moisture, transport): (FiniteMoistureSolution, MoistureRadial, Report),
        cx: &Cx<'_>,
    ) -> Result<
        (FiniteMoistureSolution, MoistureRadialSolidStep<Report>),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        let map = |error| match error {
            MoistureError::Cancelled => MoistureAdvanceError::Cancelled,
            error => MoistureAdvanceError::FinalInventory(error.into()),
        };
        let solid = self.resolve(&moisture, cx).map_err(map)?;
        poll(cx).map_err(map)?;
        Ok((
            reservoir,
            MoistureRadialSolidStep {
                moisture,
                transport,
                solid,
            },
        ))
    }
}

//! Finite aqueous activity and its dependent native cylindrical solid endpoint.

use super::*;
use crate::moisture::solution::{
    FiniteMoistureSolution, RadialSolutionMoistureCandidate, RadialSolutionSdirk2Candidate,
    SolutionMoistureError, SolutionMoistureSolve,
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

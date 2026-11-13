//! Original finite-air time integration and its dependent native solid endpoint.

use super::*;
use crate::moisture::radial::atmosphere::RadialAtmosphereSdirk2Candidate;

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl MoistureRadialSolid {
    /// Publish finite vapor, both radial SDIRK stages and native solid mechanics
    /// together. Source-dependent humidity is solved at each actual RK stage.
    /// A late source, equilibrium or geometry refusal preserves both inputs.
    pub fn advance_finite_atmosphere_sdirk2(
        &self,
        air: &atmosphere::FiniteMoistureAtmosphere,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: atmosphere::AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        (RadialAtmosphereSdirk2Candidate, MoistureSolidState),
        atmosphere::AtmosphereMoistureError,
    > {
        let candidate =
            air.advance_radial_sdirk2(moisture, dt_s, faces, transport, nonlinear, physical, cx)?;
        let solid = self.resolve(candidate.moisture(), cx)?;
        poll(cx)?;
        Ok((candidate, solid))
    }

    /// Complete an accuracy-controlled finite-air backward-Euler duration and
    /// resolve its compatible native shape, stress and wet vibration pencil.
    /// No accepted private water prefix escapes a later mechanical refusal.
    pub fn advance_finite_atmosphere_adaptive_to(
        &self,
        air: &atmosphere::FiniteMoistureAtmosphere,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: atmosphere::AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            atmosphere::FiniteMoistureAtmosphere,
            MoistureRadialSolidStep<MoistureAdvanceReport<atmosphere::AtmosphereMoistureStep, 2>>,
        ),
        MoistureAdvanceError<atmosphere::AtmosphereMoistureError>,
    > {
        let candidate = air.advance_radial_adaptive_to(
            moisture, target_s, faces, transport, nonlinear, physical, time, cx,
        )?;
        self.resolve_finite_adaptive(candidate, cx)
    }

    /// The same atomic final mechanical response with second-order finite-air
    /// time-error control. Geometry and material axes retain the original
    /// small-strain/reference-chart model. Gas pressure is not automatically
    /// applied as traction, and sorption heat is not reconstructed from loss.
    pub fn advance_finite_atmosphere_sdirk2_adaptive_to(
        &self,
        air: &atmosphere::FiniteMoistureAtmosphere,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<atmosphere::DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: atmosphere::AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            atmosphere::FiniteMoistureAtmosphere,
            MoistureRadialSolidStep<
                MoistureAdvanceReport<atmosphere::AtmosphereMoistureStep<MoistureSdirk2Report>, 2>,
            >,
        ),
        MoistureAdvanceError<atmosphere::AtmosphereMoistureError>,
    > {
        let candidate = air.advance_radial_sdirk2_adaptive_to(
            moisture, target_s, faces, transport, nonlinear, physical, time, cx,
        )?;
        self.resolve_finite_adaptive(candidate, cx)
    }

    fn resolve_finite_adaptive<Report>(
        &self,
        (air, moisture, transport): (atmosphere::FiniteMoistureAtmosphere, MoistureRadial, Report),
        cx: &Cx<'_>,
    ) -> Result<
        (
            atmosphere::FiniteMoistureAtmosphere,
            MoistureRadialSolidStep<Report>,
        ),
        MoistureAdvanceError<atmosphere::AtmosphereMoistureError>,
    > {
        let map = |error| match error {
            MoistureError::Cancelled => MoistureAdvanceError::Cancelled,
            error => MoistureAdvanceError::FinalInventory(error.into()),
        };
        let solid = self.resolve(&moisture, cx).map_err(map)?;
        poll(cx).map_err(map)?;
        Ok((
            air,
            MoistureRadialSolidStep {
                moisture,
                transport,
                solid,
            },
        ))
    }
}

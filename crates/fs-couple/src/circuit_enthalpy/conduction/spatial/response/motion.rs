//! Original nonlinear conductor-frame response in finite thermal feedback.
//!
//! The stationary endpoint owner retains all source/phase refresh, caloric
//! admission, contact routing and work budgets. This adapter supplies only the
//! actual native moving field and its independently recovered terminal work.

use super::*;
use fs_em::spatial::conduction::{
    motion::{SpatialConductionMotion, SpatialMovingConductionSnapshot},
    response::motion::SpatialMovingConductionResponseDrive,
};

/// Original constitutive closure plus native laboratory/conductor-frame fields,
/// currents, Lorentz force and signed work at the accepted thermal endpoint.
pub type SpatialMovingConductionResponseThermalStep =
    SpatialConductionResponseThermalStep<SpatialMovingConductionSnapshot>;

/// A moving nonlinear source state and its original finite thermal history.
/// The complete heat provider sample remains typed and available in the report.
pub type SpatialMovingConductionResponseEnthalpyCandidate<'t, S = MaterialNetworkSample> =
    SpatialConductionResponseEnthalpyCandidate<'t, S, SpatialMovingConductionSnapshot>;

impl ResponseHeatField for SpatialMovingConductionSnapshot {
    fn joule_power_w(&self) -> f64 {
        SpatialMovingConductionSnapshot::joule_power_w(self)
    }
    fn cell_joule_power_w(&self) -> &[f64] {
        SpatialMovingConductionSnapshot::cell_joule_power_w(self)
    }
    fn contact_response(&self) -> &[SpatialConductionContactResponse] {
        SpatialMovingConductionSnapshot::contact_response(self)
    }
    fn mechanical_power_w(&self) -> f64 {
        SpatialMovingConductionSnapshot::mechanical_power_w(self)
    }
}

impl SpatialConductionThermalNetwork<'_> {
    /// Close the original rho/sigma(T, phase, |E'|, |J|), native moving current
    /// and finite enthalpy at one backward-Euler endpoint, where E'=E+u cross B.
    /// Only native conductor-frame Joule loss heats the original material cells.
    /// Signed terminal input minus signed Lorentz output closes that heat in
    /// motoring, braking and generation; numerical closure defects add no heat.
    ///
    /// Prescribed motion/B are held over this interval on the original fixed
    /// domain. Material inventory must remain represented by its original
    /// cells; this path supplies no mesh advection, induction, magnetic storage
    /// or momentum evolution. Direct boundary and grouped-terminal controls
    /// keep their native moving-work conventions and shared solve budgets.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_response_candidate<'t, 'd>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConductionResponse,
        drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionResponseEnthalpyCandidate<'t>, Error> {
        self.step_moving_response_contact_candidate(
            thermal,
            target_s,
            problem,
            drive,
            motion,
            &[],
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Include the original finite stationary contact laws and their explicit
    /// two-body heat partitions. Both contact source temperatures/phases refresh
    /// at every thermal proposal; nonlinear trace closure and adaptive contact
    /// quadrature remain native and consume the same global field-solve budget.
    ///
    /// Both cells adjacent to an electrical contact must have zero prescribed
    /// velocity, as required by the native stationary surface-work law. Other
    /// cells may move. A failed material query, work gate or cancellation
    /// publishes neither electrical state nor thermal history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_response_contact_candidate<'t, 'd>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConductionResponse,
        drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionResponseEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        self.admit_mesh_step(
            thermal,
            target_s,
            problem.mesh(),
            problem.contacts().len(),
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_moving_response_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            motion,
            contact_bindings,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy::conduction) fn step_moving_response_resolved_candidate<
        't,
        'd,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &SpatialConductionResponse,
        drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
        motion: SpatialConductionMotion<'_>,
        contact_bindings: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'t>,
        )
            -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialMovingConductionResponseEnthalpyCandidate<'t, S>, CircuitEnthalpyError<E>>
    {
        let drive = drive.into();
        self.step_response_driver_candidate(
            thermal,
            target_s,
            problem,
            contact_bindings,
            assembly,
            response,
            controls,
            initial_energy_increments_j,
            resolve_sources,
            sample,
            |sources, contacts, bounded| {
                let candidate = problem
                    .solve_with_motion_and_materials_and_contacts_candidate(
                        sources,
                        contacts,
                        drive,
                        motion,
                        assembly,
                        field_controls,
                        bounded,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| response_solve_error(error, bounded, response))?;
                let (next, field) = candidate.into_parts();
                // Both values are native independent boundary/electrode work,
                // never reconstructed by adding heat and mechanical output.
                let supplied_power_w = match &field.terminals {
                    Some(terminals) => terminals.supplied_power_w,
                    None => field.field.terminal_power_w(),
                };
                Ok((next, field, supplied_power_w))
            },
            cancel,
        )
    }
}

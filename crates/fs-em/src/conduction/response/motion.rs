//! Original nonlinear source response in a prescribed moving planar conductor.
//!
//! The existing source/contact iteration uses the native conductor-frame field
//! E'=E+u cross B. Native laboratory E, Joule heat, Lorentz output and signed
//! terminal work remain distinct. Geometry and motion are held snapshots; this
//! owner neither moves a surface nor evolves induction or mechanical inventory.

use super::*;
use crate::conduction::motion::{ConductionMotion, ConductionMotionReport};

/// Original native field and source closure paired with mandatory motion work.
/// `field.cell_power_w` and `field.power_w` are irreversible heat; laboratory
/// electrical work and terminal supply may be negative during generation.
/// `field.power_balance_residual_w` is heat + mechanical - terminal work.
/// Source samples in `field.evaluated` use the actual conductor-frame E'.
#[derive(Clone, Debug)]
pub struct MovingConductionResponseReport {
    pub field: ConductionResponseReport,
    pub motion: ConductionMotionReport,
}

/// Accepted original source owner and its actual moving-field response.
/// Consume both together to continue. No scalable nonlinear resistance or
/// independent physical clock is produced.
pub struct MovingConductionResponseCandidate {
    problem: PlanarConductionResponse,
    report: MovingConductionResponseReport,
}
impl MovingConductionResponseCandidate {
    pub fn problem(&self) -> &PlanarConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &MovingConductionResponseReport {
        &self.report
    }
    pub fn into_parts(self) -> (PlanarConductionResponse, MovingConductionResponseReport) {
        (self.problem, self.report)
    }
}

impl PlanarConductionResponse {
    /// Close original scalar bulk/contact sources under prescribed planar
    /// motion and actual boundary voltage/current-density traces.
    ///
    /// The native field enforces charge continuity with J=Sigma(E+u cross B).
    /// Original bulk laws see E' and J; stationary contact laws retain their
    /// actual trace jump and density. Both adjacent cells of every contact must
    /// have zero prescribed velocity, as required by the original motion owner.
    /// No moving/sliding interface or axisymmetric-motion law is introduced.
    ///
    /// Native local/global work gates, original source admission, coordinate
    /// relaxation, contact refinement and their complete work budgets remain
    /// active. All discarded trials share the same accepted physical history.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_motion_candidate(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        current_traces: &[ConductionCurrentTrace],
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionResponseCandidate, ConductionResponseError> {
        self.solve_with_motion_and_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            boundary_voltage_v,
            current_traces,
            motion,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// The same closure at explicitly refreshed external thermal/phase states.
    /// Complete source pins, geometry and fixed coordinates remain owned by
    /// this response owner. Only the accepted material and native work escape.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_motion_and_materials_and_contacts_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        boundary_voltage_v: &[(usize, f64)],
        current_traces: &[ConductionCurrentTrace],
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionResponseCandidate, ConductionResponseError> {
        poll(cancelled)?;
        if !matches!(self.mesh().metric(), Metric2::Planar { .. }) {
            return Err(ConductionResponseError::Invalid(
                "planar motion response requires the original planar metric",
            ));
        }
        let (candidate, motion) = self.solve_source_frame_candidate(
            seeds,
            contacts,
            assembly,
            field_controls,
            controls,
            cancelled,
            |native, cancelled| {
                native
                    .solve_with_motion(
                        boundary_voltage_v,
                        current_traces,
                        motion,
                        field_controls,
                        &mut || cancelled(),
                    )
                    .map(|field| field.into_response_parts())
                    .map_err(ConductionResponseError::from)
            },
            |_, motion: &ConductionMotionReport, cell| motion.conductor_field_v_m[cell],
            |motion: &ConductionMotionReport, observation, cancelled| {
                // Refinement must also resolve signed work; equal heat alone
                // cannot hide unresolved opposing Lorentz transfers.
                for &power in &motion.cell_mechanical_power_w {
                    if cancelled() {
                        return Err(ConductionResponseError::Cancelled);
                    }
                    observation.power_w.push(power);
                }
                observation.power_w.extend([
                    motion.mechanical_power_w,
                    motion.electrical_power_w,
                    motion.terminal_power_w,
                ]);
                Ok(())
            },
        )?;
        let (problem, field) = candidate.into_parts();
        Ok(MovingConductionResponseCandidate {
            problem,
            report: MovingConductionResponseReport { field, motion },
        })
    }
}

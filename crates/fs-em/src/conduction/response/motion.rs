//! Original nonlinear source response in a prescribed moving planar conductor.
//!
//! The existing source/contact iteration uses the native conductor-frame field
//! E'=E+u cross B. Native laboratory E, Joule heat, Lorentz output and signed
//! terminal work remain distinct. Geometry and motion are held snapshots; this
//! owner neither moves a surface nor evolves induction or mechanical inventory.

use super::terminals::{ConductionTerminalDrive, ConductionTerminalReport};
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
    /// Original ordered direct-electrode reactions and authored source work.
    /// Absent for voltage/current-density traces.
    pub terminals: Option<ConductionTerminalReport>,
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

#[derive(Clone, Copy)]
enum MovingDrive<'a> {
    Traces {
        voltage: &'a [(usize, f64)],
        currents: &'a [ConductionCurrentTrace],
    },
    Terminals(ConductionTerminalDrive<'a>),
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
        self.solve_moving_driven_candidate(
            seeds,
            contacts,
            MovingDrive::Traces {
                voltage: boundary_voltage_v,
                currents: current_traces,
            },
            motion,
            assembly,
            field_controls,
            controls,
            usize::MAX,
            cancelled,
        )
    }

    /// Close original nonlinear sources under direct grouped moving electrodes.
    /// Each trial solves once at the actual voltage/current/floating/mixed drive.
    /// The accepted report retains the independently authored source work.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_terminals_with_motion_candidate(
        &self,
        drive: ConductionTerminalDrive<'_>,
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        maximum_field_solves: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionResponseCandidate, ConductionResponseError> {
        self.solve_terminals_with_motion_and_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            drive,
            motion,
            assembly,
            field_controls,
            controls,
            maximum_field_solves,
            cancelled,
        )
    }

    /// Direct moving electrodes with complete externally refreshed original
    /// bulk/contact sources, such as one proposed thermal/phase endpoint.
    ///
    /// The original source pins, operating axes, relaxation, quadrature and
    /// native field gates remain binding. Every discarded constitutive and
    /// contact-refinement trial consumes maximum_field_solves; zero refuses
    /// before a native solve. Original sources and geometry are never rescaled.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_terminals_with_motion_and_materials_and_contacts_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        drive: ConductionTerminalDrive<'_>,
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        maximum_field_solves: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionResponseCandidate, ConductionResponseError> {
        self.solve_moving_driven_candidate(
            seeds,
            contacts,
            MovingDrive::Terminals(drive),
            motion,
            assembly,
            field_controls,
            controls,
            maximum_field_solves,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_moving_driven_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        drive: MovingDrive<'_>,
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        maximum_field_solves: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionResponseCandidate, ConductionResponseError> {
        poll(cancelled)?;
        controls.validate()?;
        if maximum_field_solves == 0 {
            return Err(ConductionResponseError::FieldBudget {
                completed: 0,
                required: 1,
                maximum: maximum_field_solves,
            });
        }
        if !matches!(self.mesh().metric(), Metric2::Planar { .. }) {
            return Err(ConductionResponseError::Invalid(
                "planar motion response requires the original planar metric",
            ));
        }
        let mut completed = 0;
        let (candidate, (motion, terminals)) = self.solve_source_frame_candidate(
            seeds,
            contacts,
            assembly,
            field_controls,
            controls,
            cancelled,
            |native, cancelled| {
                if completed >= maximum_field_solves {
                    return Err(ConductionResponseError::FieldBudget {
                        completed,
                        required: 1,
                        maximum: maximum_field_solves,
                    });
                }
                let (field, terminals) = match drive {
                    MovingDrive::Traces { voltage, currents } => (
                        native.solve_with_motion(
                            voltage,
                            currents,
                            motion,
                            field_controls,
                            &mut || cancelled(),
                        )?,
                        None,
                    ),
                    MovingDrive::Terminals(drive) => {
                        let (field, report) = native.solve_terminals_with_motion(
                            drive.terminals,
                            drive.control,
                            drive.limits,
                            motion,
                            field_controls,
                            &mut || cancelled(),
                        )?;
                        (field, Some(report))
                    }
                };
                completed += 1;
                let (field, motion) = field.into_response_parts();
                Ok((field, (motion, terminals)))
            },
            |_, driver: &(ConductionMotionReport, Option<ConductionTerminalReport>), cell| {
                driver.0.conductor_field_v_m[cell]
            },
            |driver: &(ConductionMotionReport, Option<ConductionTerminalReport>),
             observation,
             cancelled| {
                let motion = &driver.0;
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
                if let Some(terminals) = &driver.1 {
                    // Preserve all direct current constraints, including the
                    // net reaction of each floating electrode, at refinement.
                    observation
                        .voltage_differences_v
                        .extend(&terminals.voltage_v);
                    observation.current_a.extend(&terminals.current_a);
                    observation.power_w.push(terminals.supplied_power_w);
                }
                Ok(())
            },
        )?;
        let (problem, field) = candidate.into_parts();
        Ok(MovingConductionResponseCandidate {
            problem,
            report: MovingConductionResponseReport {
                field,
                motion,
                terminals,
            },
        })
    }
}

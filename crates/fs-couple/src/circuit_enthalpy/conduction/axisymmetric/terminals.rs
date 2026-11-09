//! Source-refreshed multi-electrode fields and the original finite heat solve.
//!
//! Native conductance bases are rebuilt at every thermal trial. Simultaneous
//! controls recover one actual field, retaining local Joule cross terms and
//! floating-electrode current transfer. No terminal response is frozen across
//! a temperature or phase change and no independent heat integrator is added.

use super::*;
use crate::circuit_enthalpy::CircuitEnthalpyError;
use fs_em::conduction::axisymmetric::{
    AxisymmetricConductanceControls, AxisymmetricMixedConductanceControls,
};
use fs_math::dd::Dd;

/// Held ideal controls in the authored terminal order. The native electrical
/// owner admits complete, finite controls and independently checks reactions
/// and source work at every recovered field.
#[derive(Clone, Copy, Debug)]
pub enum AxisymmetricTerminalControl<'a> {
    /// Potential [V] of every terminal; only differences affect the field.
    Voltages(&'a [f64]),
    /// Total inward current [A] of every terminal. Their sum must vanish within
    /// the native tolerance. Zero current floats an equipotential electrode;
    /// it does not insulate its individual faces.
    Currents(&'a [f64]),
    /// A nonempty proper subset prescribes total inward currents. Current
    /// values follow `current_terminals`; voltage values follow its complement
    /// in original terminal order. Voltage-source reactions close the circuit.
    Mixed {
        current_terminals: &'a [usize],
        currents_a: &'a [f64],
        voltages_v: &'a [f64],
        controls: AxisymmetricMixedConductanceControls,
    },
}

/// Complete ideal terminal drive for one finite electrothermal interval.
/// Each terminal is a disjoint union of complete positive-area swept surfaces.
/// Shared electrodes may join disconnected mesh pieces; the native owner
/// requires one connected terminal response through bulk, contacts or electrodes.
#[derive(Clone, Copy, Debug)]
pub struct AxisymmetricTerminalDrive<'a> {
    pub terminals: &'a [&'a [usize]],
    pub control: AxisymmetricTerminalControl<'a>,
    /// Native limits for N-1 basis solves, retained values and dense work.
    /// The physical probe voltage must remain inside each refreshed source.
    pub response: AxisymmetricConductanceControls,
}

#[derive(Clone, Copy)]
pub(super) enum Drive<'a> {
    Field(ConductionDrive<'a>),
    Terminals(AxisymmetricTerminalDrive<'a>),
}

impl Drive<'_> {
    // Charge every native basis solve before starting a trial. Superposition
    // and mixed-factor recovery add no Krylov solve; native storage/dense-work
    // bounds and cancellation still apply to those operations.
    pub(super) fn field_trials(self, remaining: usize) -> Result<usize, Error> {
        let count = match self {
            Self::Field(_) => 1,
            Self::Terminals(drive) => drive
                .terminals
                .len()
                .checked_sub(1)
                .filter(|&n| n > 0)
                .ok_or_else(|| invalid("at least two physical terminals required"))?,
        };
        if count > remaining {
            return Err(CircuitEnthalpyError::CircuitBudget);
        }
        Ok(count)
    }

    pub(super) fn solve<'p, C: PlanarConductorMaterial>(
        self,
        problem: &'p AxisymmetricConduction<C>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'p, C>, Error> {
        match self {
            Self::Field(drive) => drive.solve_axisymmetric(problem, controls, cancelled),
            Self::Terminals(drive) => {
                let response = problem
                    .conductance_matrix(drive.terminals, drive.response, controls, cancelled)
                    .map_err(field_error)?;
                match drive.control {
                    AxisymmetricTerminalControl::Voltages(voltages) => response.field_at_voltages(
                        voltages,
                        controls.current_tolerance_a,
                        controls.power_tolerance_w,
                        cancelled,
                    ),
                    AxisymmetricTerminalControl::Currents(currents) => response.field_at_currents(
                        currents,
                        controls.current_tolerance_a,
                        controls.power_tolerance_w,
                        cancelled,
                    ),
                    AxisymmetricTerminalControl::Mixed {
                        current_terminals,
                        currents_a,
                        voltages_v,
                        controls: mixed_controls,
                    } => {
                        let mixed = response
                            .prepare_current_control(current_terminals, mixed_controls, cancelled)
                            .map_err(field_error)?;
                        mixed.field_at_controls(
                            voltages_v,
                            currents_a,
                            controls.current_tolerance_a,
                            controls.power_tolerance_w,
                            cancelled,
                        )
                    }
                }
                .map_err(field_error)
            }
        }
    }

    pub(super) fn supplied_power_w<C: PlanarConductorMaterial>(
        self,
        field: &AxisymmetricConductionField<'_, C>,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<f64, Error> {
        match self {
            Self::Field(drive) => drive.supplied_power_w(
                field.voltage_v(),
                field.power_w(),
                field.power_balance_residual_w(),
            ),
            Self::Terminals(drive) => drive.supplied_power_w(field, cancelled),
        }
    }
}

impl AxisymmetricTerminalDrive<'_> {
    // Current-source work uses AUTHORED totals, independently of recovered
    // reactions. Voltage-source work uses the actual native reactions. The
    // native field's retained relative/reference potentials supply physical
    // voltage differences; common bias never multiplies a current residual.
    fn supplied_power_w<C: PlanarConductorMaterial>(
        self,
        field: &AxisymmetricConductionField<'_, C>,
        cancelled: &RefCell<impl FnMut() -> bool>,
    ) -> Result<f64, Error> {
        let gauge = self
            .terminals
            .first()
            .and_then(|vertices| vertices.first())
            .copied()
            .ok_or_else(|| invalid("terminal voltage reference is absent"))?;
        let relative = field.relative_voltage_v();
        let references = field.reference_voltage_v();
        let mut work = Dd::from_f64(0.0);
        for (terminal, vertices) in self.terminals.iter().enumerate() {
            check(cancelled)?;
            let vertex = *vertices
                .first()
                .ok_or_else(|| invalid("terminal lost its actual field potential"))?;
            // Each electrode was admitted by the native response. Subtract
            // references before combining with the retained small potentials.
            let voltage = (Dd::from_f64(relative[vertex]) - Dd::from_f64(relative[gauge]))
                + (Dd::from_f64(references[vertex]) - Dd::from_f64(references[gauge]));
            let authored =
                match self.control {
                    AxisymmetricTerminalControl::Voltages(_) => None,
                    AxisymmetricTerminalControl::Currents(currents) => Some(
                        *currents
                            .get(terminal)
                            .ok_or_else(|| invalid("terminal lost its authored current"))?,
                    ),
                    AxisymmetricTerminalControl::Mixed {
                        current_terminals,
                        currents_a,
                        ..
                    } => {
                        if let Some(index) = current_terminals
                            .iter()
                            .position(|&index| index == terminal)
                        {
                            Some(*currents_a.get(index).ok_or_else(|| {
                                invalid("mixed terminal lost its authored current")
                            })?)
                        } else {
                            None
                        }
                    }
                };
            let current = if let Some(current) = authored {
                Dd::from_f64(current)
            } else {
                let mut total = Dd::from_f64(0.0);
                for &vertex in *vertices {
                    check(cancelled)?;
                    let value = field
                        .nodal_boundary_current_a()
                        .get(vertex)
                        .copied()
                        .flatten()
                        .ok_or_else(|| invalid("terminal lost its native current reaction"))?;
                    total = total + Dd::from_f64(value);
                }
                total
            };
            work = work + voltage * current;
        }
        Ok(finite(work.to_f64())?)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Advance a multi-electrode axisymmetric field and original finite thermal
    /// cells to one common endpoint. Ideal total-current and mixed controls,
    /// including floating electrodes, are resolved against freshly prepared
    /// material/contact conductance at EVERY residual evaluation.
    ///
    /// Each interface requires its existing explicit source/heat partition;
    /// use an empty contact slice only when the field has no interfaces.
    /// Complete swept-volume routes, fixed source selections, actual thermal
    /// and declared phase coordinates, Fourier links and atomic publication
    /// retain the original axisymmetric coupling semantics.
    ///
    /// Native response controls bound each preparation. In this entry,
    /// `maximum_circuit_trials` and the report count every native basis solve:
    /// N-1 per thermal evaluation, including stencils and fresh final acceptance.
    /// Insufficient remaining budget refuses before beginning those solves.
    /// Recovered fields use terminal zero as the reporting voltage gauge.
    /// This is a held ideal-source, quasistatic, first-order thermal interval;
    /// it introduces no electrical storage or terminal contact model.
    #[allow(clippy::too_many_arguments)]
    pub fn step_axisymmetric_terminal_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConduction<C>,
        drive: AxisymmetricTerminalDrive<'_>,
        contacts: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        self.step_axisymmetric_with_contacts_candidate(
            thermal,
            target_s,
            problem,
            Drive::Terminals(drive),
            contacts,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            &cancel,
        )
    }
}

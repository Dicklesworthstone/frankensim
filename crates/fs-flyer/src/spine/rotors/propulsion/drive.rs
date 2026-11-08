//! Configured shaft sources and finite mechanical supplies in native flight.
//!
//! The existing aerodynamic and geared-motion owners still do every physical
//! force and motion solve. This adapter evaluates the configured source at
//! their two actual stages, transforms its torque through an explicit gear
//! coordinate, and debits only the source work of the accepted motion.

use super::{GearedPropulsion, GearedState, PropulsionInputs, PropulsionStep};
use crate::Refusal;
use fs_airscrew::drive::{DriveOperatingPoint, DriveWork};
use fs_math::dd::Dd;

/// Runtime engineering controls and environment for one real motion stage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DriveInputs {
    /// Positive air density [kg/m³].
    pub density_kg_m3: f64,
    /// Prescribed canard setting [rad], in the aerodynamic model's domain.
    pub canard_rad: f64,
    /// Requested source command in `[0,1]`, replacing the design's trim command.
    /// This does not infer a historical pilot throttle-control topology.
    pub command: f64,
    /// World-NED gravitational acceleration [m/s²], included exactly once.
    pub gravity_world_mps2: [f64; 3],
}

/// Accepted aircraft/shaft state and its explicitly selected mechanical supply.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenState {
    /// Full native gear/carrier state, including continuing shaft momentum.
    pub motion: GearedState,
    /// Remaining usable ENGINE-shaft energy before drivetrain losses [J].
    /// `None` is an externally sustained source; `Some(0.0)` removes source
    /// torque while retaining bearing resistance. This is not chemical fuel
    /// energy, battery charge, or an inferred conversion efficiency.
    pub remaining_engine_energy_j: Option<f64>,
}

/// Accepted native motion with source, loss and supply accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct DrivenStep<const N: usize> {
    /// Complete accepted state; the input remains unchanged on any refusal.
    pub state: DrivenState,
    /// Original two aerodynamic stages and full physical mechanical accounting.
    pub propulsion: PropulsionStep<N>,
    /// Source used by the start-of-step predictor.
    pub initial_drive: DriveOperatingPoint,
    /// Source used by the accepted midpoint motion step.
    pub midpoint_drive: DriveOperatingPoint,
    /// Accepted relative travel of the drive's own output coordinate [rad].
    pub drive_output_travel_rad: f64,
    /// Work of the accepted midpoint source and its two distinct physical losses.
    pub work: DriveWork,
    /// Finite supply decrease minus engine work [J], or zero for an external
    /// source. This is floating-point accounting error, not physical loss.
    pub supply_residual_j: f64,
    /// Drive output work minus the native motor-port work [J]. A gear-coordinate
    /// transform can introduce roundoff; it must not become another loss term.
    pub drive_work_residual_j: f64,
}

/// A configured design source bound to a native propeller/airframe model.
pub struct DrivenPropulsion<'model, 'design, const N: usize> {
    propulsion: &'model GearedPropulsion<'design, N>,
    output_per_gear_ratio: f64,
}

fn refuse(code: &'static str, message: &str, repair: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![repair.into()],
    }
}

fn source_error(error: fs_airscrew::Refusal) -> Refusal {
    Refusal {
        code: error.code,
        message: error.message,
        ranked_repairs: error.ranked_repairs,
    }
}

impl<'model, 'design, const N: usize> DrivenPropulsion<'model, 'design, N> {
    /// Bind `OpenLoopDesign::shaft_drive` to a declared coordinate transform:
    /// `drive_output_angle = output_per_gear_ratio * common_gear_angle`.
    /// Thus output speed is multiplied by this ratio, and the generalized motor
    /// torque is multiplied by the SAME ratio. Their work is invariant.
    /// Use `1.0` only when the assembly's common angle is the drive output angle.
    /// No rotor inertia, engine inertia or propeller ratio is inferred.
    ///
    /// # Errors
    /// `driven-propulsion-invalid` for a nonpositive or nonfinite ratio.
    pub fn new(
        propulsion: &'model GearedPropulsion<'design, N>,
        output_per_gear_ratio: f64,
    ) -> Result<Self, Refusal> {
        if !output_per_gear_ratio.is_finite() || output_per_gear_ratio <= 0.0 {
            return Err(refuse(
                "driven-propulsion-invalid",
                "drive output/common gear ratio must be finite and positive",
                "declare the mechanical source coordinate relative to the admitted gear train",
            ));
        }
        Ok(Self {
            propulsion,
            output_per_gear_ratio,
        })
    }

    /// Advance using runtime commands and the configured source at the original
    /// start/midpoint motion stages. Each callback receives the actual stage
    /// time and full gear/carrier state; no already accepted motion is rescaled.
    ///
    /// Finite supplies debit ENGINE work, including energy subsequently lost in
    /// the transmission. An empty supply commands zero source at both stages;
    /// bearing drag and aerodynamic shaft loads remain physical. Requested
    /// commands must still be valid even when the supply is empty.
    ///
    /// The source/loss ledger uses the midpoint torques times the native step's
    /// accepted relative travel, not instantaneous power times `dt`. This keeps
    /// its work consistent with the mechanical owner. Smooth source dynamics
    /// inherit the original second-order staging; the instantaneous source
    /// rating is not an exact nonlinear power bound over an arbitrary large
    /// frozen-torque interval. Stopped/reversed propeller aerodynamics retain
    /// their original domain refusals; this does not add a windmilling model.
    ///
    /// # Errors
    /// Existing native/source refusals propagate. `shaft-energy-exhausted`
    /// refuses a depletion crossing: it publishes no state or consumed supply.
    /// Reduce/localize the interval before continuing; this method does not
    /// invent an exhaustion time or silently lower the source command.
    /// `shaft-energy-unrepresentable` refuses an erased positive supply debit.
    /// Reverse endpoint motion refuses rather than claiming positive loss over
    /// an unresolved reversal. Callback side effects remain caller-owned.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        state: &DrivenState,
        t_s: f64,
        dt_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<[f64; 2]>,
        mut inputs: C,
        wind_world: W,
    ) -> Result<DrivenStep<N>, Refusal>
    where
        C: FnMut(f64, &GearedState) -> Result<DriveInputs, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        if state
            .remaining_engine_energy_j
            .is_some_and(|value| !value.is_finite() || value < 0.0)
        {
            return Err(refuse(
                "driven-propulsion-invalid",
                "finite source requires nonnegative finite remaining engine energy",
                "supply usable engine-shaft joules or explicitly select an external source",
            ));
        }
        let empty = state.remaining_engine_energy_j == Some(0.0);
        let drive = self.propulsion.design.shaft_drive;
        let ratio = self.output_per_gear_ratio;
        let mut stages = [None; 2];
        let mut stage_index = 0;
        let propulsion = self.propulsion.step(
            &state.motion,
            t_s,
            dt_s,
            angular_origin_world,
            warm,
            |time, stage| {
                let input = inputs(time, stage)?;
                let output_speed = ratio * self.propulsion.assembly.speed(stage)?;
                let requested = drive
                    .evaluate(output_speed, input.command)
                    .map_err(source_error)?;
                let (command, point) = if empty {
                    (
                        0.0,
                        drive.evaluate(output_speed, 0.0).map_err(source_error)?,
                    )
                } else {
                    (input.command, requested)
                };
                let motor_nm = ratio * point.prop_torque_nm;
                if !motor_nm.is_finite() || (point.prop_torque_nm != 0.0 && motor_nm == 0.0) {
                    return Err(refuse(
                        "driven-propulsion-invalid",
                        "drive torque cannot be represented in the common gear coordinate",
                        "use a representable output/common gear ratio and mechanical source",
                    ));
                }
                stages[stage_index] = Some((output_speed, command, point));
                stage_index += 1;
                Ok(PropulsionInputs {
                    density_kg_m3: input.density_kg_m3,
                    canard_rad: input.canard_rad,
                    motor_nm,
                    gravity_world_mps2: input.gravity_world_mps2,
                })
            },
            wind_world,
        )?;
        let speed_before = self.propulsion.assembly.speed(&state.motion)?;
        let speed_after = propulsion.motion.speed_rad_s;
        if speed_after < 0.0 {
            return Err(refuse(
                "driven-propulsion-invalid",
                "accepted motion would cross the unidirectional source's zero shaft speed",
                "localize the stop before continuing with an admitted stopped/reverse model",
            ));
        }
        // This is the native mechanical owner's increment, before addition to
        // the unwrapped phase. Subtracting large accumulated phases loses work.
        let gear_travel = dt_s * (0.5 * speed_before + 0.5 * speed_after);
        let drive_output_travel_rad = ratio * gear_travel;
        let initial_drive = stages[0].expect("native propulsion evaluates two stages").2;
        let (midpoint_speed, midpoint_command, midpoint_drive) =
            stages[1].expect("native propulsion evaluates two stages");
        let work = drive
            .work_over_angle(midpoint_speed, midpoint_command, drive_output_travel_rad)
            .map_err(source_error)?;
        let (remaining_engine_energy_j, supply_residual_j) = if let Some(available) =
            state.remaining_engine_energy_j
        {
            if work.engine_work_j > available {
                return Err(refuse(
                    "shaft-energy-exhausted",
                    &format!(
                        "step requires {} J at the engine shaft; only {available} J remains",
                        work.engine_work_j
                    ),
                    "shorten/localize this step at source exhaustion; no endpoint was committed",
                ));
            }
            let remaining = available - work.engine_work_j;
            if work.engine_work_j > 0.0 && remaining == available {
                return Err(refuse(
                    "shaft-energy-unrepresentable",
                    "positive engine work is erased by the finite supply's energy scale",
                    "use a representable energy scale or accumulate a larger physical interval",
                ));
            }
            let residual = (Dd::from_f64(available)
                - Dd::from_f64(remaining)
                - Dd::from_f64(work.engine_work_j))
            .to_f64();
            (Some(remaining), residual)
        } else {
            (None, 0.0)
        };
        let drive_work_residual_j = work.prop_work_j - propulsion.motion.motor_work_j;
        if !supply_residual_j.is_finite() || !drive_work_residual_j.is_finite() {
            return Err(refuse(
                "driven-propulsion-invalid",
                "source/gear work accounting overflow",
                "reduce the mechanical source or interval to a representable work scale",
            ));
        }
        Ok(DrivenStep {
            state: DrivenState {
                motion: propulsion.motion.state,
                remaining_engine_energy_j,
            },
            propulsion,
            initial_drive,
            midpoint_drive,
            drive_output_travel_rad,
            work,
            supply_residual_j,
            drive_work_residual_j,
        })
    }
}

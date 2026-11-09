//! Continuing winding current and full native body/control acceleration.
//!
//! For full duration dt, i_mid=(2 L i0+dt V-dt K omega_mid)/(2 L+dt R).
//! Hence K*i_mid is a constant joint effort minus a positive impedance.
//! Eliminating that impedance in the ORIGINAL articulated mass solve retains
//! dry/fluid cross terms, back-EMF feedback, and equal-and-opposite body recoil.
//! No scalar substitute inertia, extra aero evaluation, or clipped current is
//! used. Actual endpoint current, source power and winding energy admit the
//! entire result before either mechanical or electrical state is returned.

use super::{
    ArticulatedLoads, ArticulatedModel, ArticulatedState, ArticulatedStep, ArticulatedStepConfig,
    GeneralizedKinematics, JointImpedance, MAX_CONTROL_COORDS, Refusal, cancelled, refuse,
};
use fs_math::dd::Dd;
use fs_mbd::articulated::JointKind;
use fs_mbd::rotordynamics::transient::drive::RotorCurrentDrive;
use fs_mbd::rotordynamics::transient::drive::voltage::{RotorVoltageControl, RotorWindingMotor};
use std::cell::Cell;

/// One held terminal-voltage command and its original motor/source owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointVoltageDrive {
    /// Native reciprocal Kt=Ke, winding storage and hardware capabilities.
    /// `shaft_coordinate` selects a compact revolute/helical joint, not a base
    /// coordinate. Constants are referred to that joint's angle and rate.
    pub motor: RotorWindingMotor,
    /// Requested terminal voltage [V], saturated to the motor's voltage cap.
    pub command_voltage_v: f64,
    /// Original source limits and numerical current/voltage/power allowances.
    pub control: RotorVoltageControl,
}

/// Continuing physical state; currents follow the supplied drive order.
#[derive(Clone, Debug, PartialEq)]
pub struct VoltageDrivenState {
    /// Original body, joint and aerodynamic-memory state.
    pub mechanical: ArticulatedState,
    /// One winding current [A] per drive; magnetic energy is derived from L i²/2.
    pub currents_a: Vec<f64>,
}

/// One actuator's accepted mechanical, magnetic and source energy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct JointVoltageReceipt {
    /// Compact native rotational coordinate.
    pub joint: usize,
    /// Held voltage after the hardware voltage limit [V].
    pub applied_voltage_v: f64,
    /// Original winding current [A].
    pub current_before_a: f64,
    /// Continuing winding current [A].
    pub current_after_a: f64,
    /// Mean of the actual endpoint currents [A].
    pub midpoint_current_a: f64,
    /// The accepted native load stage's relative joint rate [rad/s].
    pub midpoint_rate_rad_s: f64,
    /// Actual affine motor effort in the shared mechanical equation [N m].
    pub torque_nm: f64,
    /// Motor work at that same native midpoint [J].
    pub joint_work_j: f64,
    /// Held-voltage work V i_mid dt [J].
    pub terminal_work_j: f64,
    /// Polarized magnetic change L i_mid (i_after-i_before) [J].
    pub magnetic_energy_change_j: f64,
    /// Midpoint Joule quadrature R i_mid² dt [J], not exact RMS-current heat.
    pub winding_heat_j: f64,
    /// Native bridge conversion/braking heat [J].
    pub bridge_heat_j: f64,
    /// Positive source draw, retaining any reconstructed current reversal [J].
    pub electrical_draw_j: f64,
    /// Positive energy returned to the source [J].
    pub electrical_return_j: f64,
    /// V-L delta(i)/dt-R i_mid-K omega_mid from actual returned states [V].
    pub winding_residual_v: f64,
    /// Draw-return-magnetic change-joint work-winding heat-bridge heat [J].
    pub electrical_balance_residual_j: f64,
}

/// One atomic native mechanical/electrical result.
#[derive(Clone, Debug, PartialEq)]
pub struct VoltageDrivenStep {
    /// Original articulated step; its nonlinear mechanical accuracy is unchanged.
    pub mechanical: ArticulatedStep,
    /// Accepted drive accounts, in original input order.
    pub drives: Vec<JointVoltageReceipt>,
}

impl VoltageDrivenStep {
    /// Carry the actual native state and winding currents into the next tick.
    #[must_use]
    pub fn into_state(self) -> VoltageDrivenState {
        VoltageDrivenState {
            mechanical: self.mechanical.state,
            currents_a: self
                .drives
                .into_iter()
                .map(|drive| drive.current_after_a)
                .collect(),
        }
    }
}

fn invalid(message: impl ToString) -> Refusal {
    refuse(
        "articulated-drive-invalid",
        message.to_string(),
        "use original motor/source parameters and one current per native rotational drive",
    )
}

fn finite(value: Dd) -> Result<f64, Refusal> {
    let rounded = value.to_f64();
    if !value.hi.is_finite() || !value.lo.is_finite() || !rounded.is_finite() {
        return Err(invalid("unrepresentable winding coupling or work"));
    }
    Ok(rounded)
}

#[derive(Clone, Copy)]
struct PreparedDrive {
    input: JointVoltageDrive,
    initial_current: f64,
    voltage: f64,
    free_current: f64,
    current_per_speed: f64,
    torque_bias: f64,
    damping: f64,
}

fn admit_current(drive: JointVoltageDrive, current: f64) -> Result<(), Refusal> {
    let excess = current.abs() - drive.motor.drive.maximum_current_a;
    if !current.is_finite() || excess > drive.control.absolute_current_tolerance_a {
        return Err(refuse(
            "articulated-drive-current-limit",
            format!("winding current {current:e} A exceeds its admitted capability"),
            "reduce commanded voltage or timestep; a current limit never clips magnetic storage",
        ));
    }
    Ok(())
}

/// Advance finite-inductance drives and the original articulated body together.
///
/// One drive may bind each native revolute/helical coordinate. The native tree
/// supplies all link/rotor inertia; no motor inertia is invented or added twice.
/// Aerodynamic hinge loads and other applied loads come from `evaluate` at the
/// SAME midpoint as the electrical reaction. Its `joint_effort` excludes the
/// motors below. The original mechanical joint and iteration limits still apply.
///
/// Voltage commands saturate at the authored terminal capability. Initial/final
/// current and source-power bounds admit the declared linear-current segment;
/// unresolved continuous-time peaks are not bounded. An infeasible current or
/// source demand refuses the complete borrowed-state step. Independent drives
/// have their own held-voltage sources; shared battery/bus dynamics, magnetic
/// saturation, thermal parameter changes and switching transients are separate.
/// Winding and joint work use midpoint quadrature. The existing nonlinear body
/// integrator retains its refinement limits; no new exact global-energy claim
/// is inferred from the electrical balance.
#[allow(clippy::too_many_arguments)]
pub fn step_with_voltage_drives<F>(
    model: &ArticulatedModel,
    impedance: &[JointImpedance],
    state: &VoltageDrivenState,
    memory_inputs: &[f64],
    drives: &[JointVoltageDrive],
    t_s: f64,
    dt_s: f64,
    config: ArticulatedStepConfig<'_>,
    evaluate: F,
    mut is_cancelled: impl FnMut() -> bool,
) -> Result<VoltageDrivenStep, Refusal>
where
    F: Fn(f64, &ArticulatedState, &GeneralizedKinematics) -> Result<ArticulatedLoads, Refusal>,
{
    if is_cancelled() {
        return Err(cancelled());
    }
    let nc = model.dof_count();
    if nc > MAX_CONTROL_COORDS
        || impedance.len() != nc
        || drives.len() > nc
        || state.currents_a.len() != drives.len()
        || state.mechanical.motion.joint_velocity.len() != nc
        || !dt_s.is_finite()
        || dt_s <= 0.0
    {
        return Err(invalid(
            "invalid drive count, native joint shape or duration",
        ));
    }
    // Adding motor damping must not hide an invalid authored passive element.
    if !impedance.iter().all(|p| {
        p.damping.is_finite()
            && p.damping >= 0.0
            && p.stiffness.is_finite()
            && p.stiffness >= 0.0
            && p.rest_position.is_finite()
    }) {
        return Err(invalid("invalid original passive joint impedance"));
    }
    let mut rotational = [false; MAX_CONTROL_COORDS];
    for (link, data) in model.links().iter().enumerate() {
        if let Some(joint) = model.dof_index(link) {
            rotational[joint] = matches!(
                data.joint().kind(),
                JointKind::Revolute | JointKind::Helical
            );
        }
    }
    let mut occupied = [false; MAX_CONTROL_COORDS];
    let mut effective = impedance.to_vec();
    let mut prepared = Vec::with_capacity(drives.len());
    let half = Dd::from_f64(dt_s) * Dd::from_f64(0.5);
    for (&input, &initial_current) in drives.iter().zip(&state.currents_a) {
        if is_cancelled() {
            return Err(cancelled());
        }
        let motor = input.motor;
        let drive = motor.drive;
        drive
            .admit(nc, input.control.source, input.command_voltage_v)
            .map_err(invalid)?;
        let joint = drive.shaft_coordinate;
        if occupied[joint]
            || !rotational[joint]
            || !motor.inductance_h.is_finite()
            || motor.inductance_h <= 0.0
            || !input.control.absolute_current_tolerance_a.is_finite()
            || input.control.absolute_current_tolerance_a < 0.0
        {
            return Err(invalid(
                "duplicate/nonrotational drive or invalid winding storage/current allowance",
            ));
        }
        occupied[joint] = true;
        admit_current(input, initial_current)?;
        drive
            .admit_speed(state.mechanical.motion.joint_velocity[joint])
            .map_err(invalid)?;
        let voltage = input.command_voltage_v.clamp(
            -drive.maximum_terminal_voltage_v,
            drive.maximum_terminal_voltage_v,
        );
        let l = Dd::from_f64(motor.inductance_h);
        let k = Dd::from_f64(drive.torque_constant_nm_per_a);
        let denominator = l + half * Dd::from_f64(drive.winding_resistance_ohm);
        if finite(denominator)? <= 0.0 {
            return Err(invalid("unresolved winding midpoint impedance"));
        }
        let free_current = finite(
            (l * Dd::from_f64(initial_current) + half * Dd::from_f64(voltage)) / denominator,
        )?;
        let current_per_speed = finite(half * k / denominator)?;
        let torque_bias = finite(k * Dd::from_f64(free_current))?;
        let damping = finite(k * Dd::from_f64(current_per_speed))?;
        effective[joint].damping =
            finite(Dd::from_f64(effective[joint].damping) + Dd::from_f64(damping))?;
        prepared.push(PreparedDrive {
            input,
            initial_current,
            voltage,
            free_current,
            current_per_speed,
            torque_bias,
            damping,
        });
    }
    // The existing owner accepts the candidate whose midpoint was evaluated.
    // Retain those actual rates; do not reconstruct them by subtracting large
    // joint positions or invoke the load callback again after acceptance.
    let midpoint_rates = Cell::new([0.0; MAX_CONTROL_COORDS]);
    let mechanical = super::step(
        model,
        &effective,
        &state.mechanical,
        memory_inputs,
        t_s,
        dt_s,
        config,
        |time, stage, map| {
            let mut loads = evaluate(time, stage, map)?;
            if loads.joint_effort.len() != nc {
                return Err(invalid("load joint shape differs from drive owner"));
            }
            for drive in &prepared {
                loads.joint_effort[drive.input.motor.drive.shaft_coordinate] += drive.torque_bias;
            }
            let mut rates = [0.0; MAX_CONTROL_COORDS];
            rates[..nc].copy_from_slice(&stage.motion.joint_velocity);
            midpoint_rates.set(rates);
            Ok(loads)
        },
        &mut is_cancelled,
    )?;
    let rates = midpoint_rates.get();
    let mut receipts = Vec::with_capacity(prepared.len());
    for drive in prepared {
        if is_cancelled() {
            return Err(cancelled());
        }
        let joint = drive.input.motor.drive.shaft_coordinate;
        drive
            .input
            .motor
            .drive
            .admit_speed(mechanical.state.motion.joint_velocity[joint])
            .map_err(invalid)?;
        receipts.push(drive.receipt(rates[joint], dt_s)?);
    }
    if is_cancelled() {
        return Err(cancelled());
    }
    Ok(VoltageDrivenStep {
        mechanical,
        drives: receipts,
    })
}

impl PreparedDrive {
    fn receipt(self, speed: f64, dt_s: f64) -> Result<JointVoltageReceipt, Refusal> {
        let motor = self.input.motor;
        let drive = motor.drive;
        let dd = Dd::from_f64;
        let half_current = dd(self.free_current) - dd(self.current_per_speed) * dd(speed);
        let current_after_a = finite(half_current * dd(2.0) - dd(self.initial_current))?;
        admit_current(self.input, current_after_a)?;
        let midpoint_current_a =
            finite((dd(self.initial_current) + dd(current_after_a)) * dd(0.5))?;
        let current = dd(midpoint_current_a);
        let torque_nm = finite(dd(self.torque_bias) - dd(self.damping) * dd(speed))?;
        let joint_work_j = finite(dd(dt_s) * dd(torque_nm) * dd(speed))?;
        let terminal_work_j = finite(dd(dt_s) * dd(self.voltage) * current)?;
        let magnetic_energy_change_j = finite(
            dd(motor.inductance_h) * current * (dd(current_after_a) - dd(self.initial_current)),
        )?;
        let winding_heat_j =
            finite(dd(dt_s) * dd(drive.winding_resistance_ohm) * current * current)?;
        let winding_residual_v = finite(
            dd(self.voltage)
                - dd(motor.inductance_h) * (dd(current_after_a) - dd(self.initial_current))
                    / dd(dt_s)
                - dd(drive.winding_resistance_ohm) * current
                - dd(drive.torque_constant_nm_per_a) * dd(speed),
        )?;
        let initial_power = finite(dd(self.voltage) * dd(self.initial_current))?;
        let final_power = finite(dd(self.voltage) * dd(current_after_a))?;
        let (positive, negative) =
            RotorCurrentDrive::directional_average(initial_power, final_power).map_err(invalid)?;
        let draw_w = finite(dd(positive) / dd(drive.motoring_efficiency))?;
        let return_w = finite(dd(negative) * dd(drive.regenerative_efficiency))?;
        let electrical_draw_j = finite(dd(dt_s) * dd(draw_w))?;
        let electrical_return_j = finite(dd(dt_s) * dd(return_w))?;
        let bridge_heat_j = finite(
            dd(dt_s)
                * (dd(draw_w) * dd(1.0 - drive.motoring_efficiency)
                    + dd(negative) * dd(1.0 - drive.regenerative_efficiency)),
        )?;
        let electrical_balance_residual_j = finite(
            dd(electrical_draw_j)
                - dd(electrical_return_j)
                - dd(magnetic_energy_change_j)
                - dd(joint_work_j)
                - dd(winding_heat_j)
                - dd(bridge_heat_j),
        )?;
        let control = self.input.control.source;
        let mut excess = 0.0_f64;
        for terminal in [initial_power, final_power] {
            let source = drive.source_power(terminal).map_err(invalid)?;
            excess = excess
                .max(source - control.maximum_source_input_w)
                .max(-source - control.maximum_source_return_w);
        }
        if excess > control.absolute_power_tolerance_w {
            return Err(refuse(
                "articulated-drive-power-limit",
                format!("source excess {excess:e} W"),
                "reduce voltage or timestep, or supply the actual source capability",
            ));
        }
        if winding_residual_v.abs() > control.absolute_voltage_tolerance_v
            || (electrical_balance_residual_j / dt_s).abs() > control.absolute_power_tolerance_w
        {
            return Err(refuse(
                "articulated-drive-winding-balance",
                format!(
                    "winding residual {winding_residual_v:e} V, energy residual {electrical_balance_residual_j:e} J"
                ),
                "use representable winding scales and explicit numerical allowances",
            ));
        }
        Ok(JointVoltageReceipt {
            joint: drive.shaft_coordinate,
            applied_voltage_v: self.voltage,
            current_before_a: self.initial_current,
            current_after_a,
            midpoint_current_a,
            midpoint_rate_rad_s: speed,
            torque_nm,
            joint_work_j,
            terminal_work_j,
            magnetic_energy_change_j,
            winding_heat_j,
            bridge_heat_j,
            electrical_draw_j,
            electrical_return_j,
            winding_residual_v,
            electrical_balance_residual_j,
        })
    }
}

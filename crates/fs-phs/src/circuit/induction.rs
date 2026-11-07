//! Balanced squirrel-cage induction motor/generator on original flux states.
//!
//! State is `[lambda_s_alpha, lambda_s_beta, lambda_r_d, lambda_r_q,
//! theta, p_rotor, p_stator]`. Rotor quantities are referred to the stator and
//! expressed in rotor-fixed axes; stator quantities use stator-fixed axes.
//! `theta` is the relative MECHANICAL angle and `phi = pole_pairs * theta`.
//! With the positive rotation `R(phi)`, the original winding law is
//! `lambda_s = (Lls+Lm) i_s + Lm R i_r`,
//! `lambda_r = (Llr+Lm) i_r + Lm R^T i_s`.
//! Its single magnetic energy supplies both currents and equal/opposite shaft
//! torque. The shorted rotor obeys `lambda_r_dot = -Rr i_r`; rotor current and
//! its magnetic lag are retained even at standstill. No prescribed slip/torque
//! curve, efficiency map, or synchronous-rotor assumption enters the solve.
//!
//! All electrical axes use the power-invariant Clarke normalization, including
//! referred rotor parameters. Torque is `p Lm cross(R i_r, i_s)` without an
//! extra 3/2 factor. See Krause, Wasynczuk, Sudhoff & Pekarek, *Analysis of
//! Electric Machinery and Drive Systems*, 3rd ed. (2013), induction-machine
//! equations; here the power-invariant transform absorbs their phase factor.
//!
//! A source-owned symmetric discrete gradient and its analytic Newton action
//! reuse the native seven-state solver. Positive physical leakage keeps the
//! inductance invertible, and separate stator/rotor copper heat and internal
//! bearing heat retain loss ownership for thermal coupling. `into_system`
//! exposes the original states for native finite-bus/converter composition.
//!
//! This is a numerical reduced model with sinusoidal balanced windings,
//! constant referred inductances/resistances and finite rotor/stator inertia.
//! It does not resolve cage-bar skin effect, slots, core saturation/hysteresis,
//! temperature dependence, harmonics, neutral current, or a 3D field solve.

use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, ResistanceWitness,
    SolveBudget, StepWorkspace,
};
use fs_math::det;

mod storage;
use storage::{InductionStorage, rotate};

/// Reuse the same physical three-phase/power-invariant two-axis convention.
pub use super::synchronous::phase_voltage_to_alpha_beta;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InductionMotorParameters {
    /// Uncoupled stator leakage [H], strictly positive.
    pub stator_leakage_inductance_h: f64,
    /// Uncoupled referred rotor leakage [H], strictly positive.
    pub rotor_leakage_inductance_h: f64,
    /// Referred magnetizing inductance [H], nonnegative. Zero decouples windings.
    pub magnetizing_inductance_h: f64,
    pub stator_resistance_ohm: f64,
    /// Shorted cage resistance referred to the stator [ohm].
    pub rotor_resistance_ohm: f64,
    pub pole_pairs: u32,
    pub rotor_inertia_kg_m2: f64,
    pub stator_inertia_kg_m2: f64,
    /// Internal bearing drag on relative shaft speed [N m s/rad].
    pub bearing_drag_nm_s_per_rad: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InductionOperatingPoint {
    pub stator_current_alpha_beta_a: [f64; 2],
    /// Referred cage current in rotor-fixed axes.
    pub rotor_current_dq_a: [f64; 2],
    pub rotor_current_alpha_beta_a: [f64; 2],
    /// Physical balanced stator phase currents [A].
    pub stator_phase_current_a: [f64; 3],
    pub shaft_speed_rad_s: [f64; 2],
    /// Torque acting on rotor; the stator receives its negative [N m].
    pub electromagnetic_torque_nm: f64,
    pub stator_copper_power_w: f64,
    pub rotor_copper_power_w: f64,
}

pub struct InductionMotor {
    parameters: InductionMotorParameters,
    system: PortHamiltonian,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InductionStep {
    pub native: PreparedStepRecord,
    pub electrical_work_j: [f64; 2],
    pub shaft_work_j: [f64; 2],
    pub stator_copper_heat_j: f64,
    pub rotor_copper_heat_j: f64,
    pub bearing_heat_j: f64,
    /// Signed impulse on rotor; stator receives the opposite impulse.
    pub electromagnetic_angular_impulse_nm_s: f64,
    pub bearing_angular_impulse_nm_s: f64,
    pub shaft_angle_increment_rad: [f64; 2],
}

fn invalid(what: &'static str) -> PhsError {
    PhsError::Dimension { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PreparedStepError> {
    if cancelled() {
        Err(PreparedStepError::Cancelled)
    } else {
        Ok(())
    }
}

impl InductionMotor {
    pub fn new(parameters: InductionMotorParameters) -> Result<Self, PhsError> {
        let storage = InductionStorage::new(parameters)?;
        let p = parameters;
        let losses = [
            p.stator_resistance_ohm,
            p.rotor_resistance_ohm,
            p.bearing_drag_nm_s_per_rad,
        ];
        if losses.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(invalid(
                "induction winding and bearing losses must be finite nonnegative",
            ));
        }
        let mut skew = vec![0.0; 49];
        skew[4 * 7 + 5] = 1.0;
        skew[4 * 7 + 6] = -1.0;
        skew[5 * 7 + 4] = -1.0;
        skew[6 * 7 + 4] = 1.0;
        let mut resistance = vec![0.0; 49];
        let mut source_rows = vec![0.0; 35];
        for i in 0..4 {
            let resistance_value = losses[i / 2];
            resistance[i * 7 + i] = resistance_value;
            source_rows[i * 7 + i] = det::sqrt(resistance_value);
        }
        let bearing = losses[2];
        resistance[5 * 7 + 5] = bearing;
        resistance[5 * 7 + 6] = -bearing;
        resistance[6 * 7 + 5] = -bearing;
        resistance[6 * 7 + 6] = bearing;
        source_rows[4 * 7 + 5] = det::sqrt(bearing);
        source_rows[4 * 7 + 6] = -det::sqrt(bearing);
        let witness = ResistanceWitness::try_from_source_rows(7, resistance.clone(), source_rows)?;
        // Voltage drives only the stator. The cage is electrically shorted.
        let mut ports = vec![0.0; 28];
        ports[0] = 1.0;
        ports[5] = 1.0;
        ports[5 * 4 + 2] = 1.0;
        ports[6 * 4 + 3] = 1.0;
        let system = PortHamiltonian::new(7, 4, skew, resistance, ports, Box::new(storage))?
            .with_resistance_witness(witness)?;
        Ok(Self { parameters, system })
    }

    pub const fn parameters(&self) -> InductionMotorParameters {
        self.parameters
    }
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Compose the four work-conjugate ports `[v_s_alpha,v_s_beta,tau_r,tau_s]`
    /// with native circuit/mechanical models while retaining all seven states.
    /// For an averaged bridge, `interconnect_mapped` with ports `[0,1]` supplies
    /// voltage from finite bus charge and returns its reciprocal load current.
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }

    /// Initialize stator current in stator axes and referred rotor current in
    /// rotor axes, relative mechanical angle, and absolute shaft speeds.
    pub fn state(
        &self,
        stator_current_alpha_beta_a: [f64; 2],
        rotor_current_dq_a: [f64; 2],
        relative_angle_rad: f64,
        shaft_speed_rad_s: [f64; 2],
    ) -> Result<[f64; 7], PhsError> {
        let p = self.parameters;
        let (s, c) = det::sin_cos(f64::from(p.pole_pairs) * relative_angle_rad);
        let rotor_ab = rotate(s, c, rotor_current_dq_a);
        let stator_dq = rotate(-s, c, stator_current_alpha_beta_a);
        let mut state = [0.0; 7];
        for i in 0..2 {
            state[i] = p.stator_leakage_inductance_h * stator_current_alpha_beta_a[i]
                + p.magnetizing_inductance_h * (stator_current_alpha_beta_a[i] + rotor_ab[i]);
            state[2 + i] = p.rotor_leakage_inductance_h * rotor_current_dq_a[i]
                + p.magnetizing_inductance_h * (rotor_current_dq_a[i] + stator_dq[i]);
        }
        state[4] = relative_angle_rad;
        state[5] = p.rotor_inertia_kg_m2 * shaft_speed_rad_s[0];
        state[6] = p.stator_inertia_kg_m2 * shaft_speed_rad_s[1];
        if stator_current_alpha_beta_a
            .iter()
            .chain(&rotor_current_dq_a)
            .chain(&shaft_speed_rad_s)
            .chain(&state)
            .any(|v| !v.is_finite())
            || !self.system.hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "induction initial current, angle, momentum or energy is unrepresentable",
            ));
        }
        Ok(state)
    }

    pub fn operating_point(&self, state: &[f64; 7]) -> Result<InductionOperatingPoint, PhsError> {
        let p = self.parameters;
        let mut effort = [0.0; 7];
        self.system.storage().gradient(state, &mut effort);
        let stator_current_alpha_beta_a = [effort[0], effort[1]];
        let rotor_current_dq_a = [effort[2], effort[3]];
        let (s, c) = det::sin_cos(f64::from(p.pole_pairs) * state[4]);
        let rotor_current_alpha_beta_a = rotate(s, c, rotor_current_dq_a);
        let a = det::sqrt(2.0 / 3.0) * effort[0];
        let q = det::sqrt(0.5) * effort[1];
        let stator_phase_current_a = [a, -0.5 * a + q, -0.5 * a - q];
        let stator_copper_power_w =
            p.stator_resistance_ohm * effort[0].mul_add(effort[0], effort[1] * effort[1]);
        let rotor_copper_power_w =
            p.rotor_resistance_ohm * effort[2].mul_add(effort[2], effort[3] * effort[3]);
        if effort
            .iter()
            .chain(&rotor_current_alpha_beta_a)
            .chain(&stator_phase_current_a)
            .chain(&[stator_copper_power_w, rotor_copper_power_w])
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "induction current, torque, speed or copper loss is unrepresentable",
            ));
        }
        Ok(InductionOperatingPoint {
            stator_current_alpha_beta_a,
            rotor_current_dq_a,
            rotor_current_alpha_beta_a,
            stator_phase_current_a,
            shaft_speed_rad_s: [effort[5], effort[6]],
            electromagnetic_torque_nm: -effort[4],
            stator_copper_power_w,
            rotor_copper_power_w,
        })
    }

    /// Allocate reusable native implicit scratch. Source queries, analytic
    /// tangent actions and subsequent accepted/trial steps allocate nothing.
    pub fn prepare_controlled(
        self,
        budget: SolveBudget,
        backtracking_halvings: u8,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedInductionMotor, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(budget)?;
        workspace.set_backtracking_limit(backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedInductionMotor {
            model: self,
            workspace,
        })
    }
}

pub struct PreparedInductionMotor {
    model: InductionMotor,
    workspace: StepWorkspace,
}

impl PreparedInductionMotor {
    pub fn system(&self) -> &PortHamiltonian {
        &self.model.system
    }
    pub const fn parameters(&self) -> InductionMotorParameters {
        self.model.parameters
    }
    pub fn operating_point(&self, state: &[f64; 7]) -> Result<InductionOperatingPoint, PhsError> {
        self.model.operating_point(state)
    }

    /// Hold `[v_s_alpha,v_s_beta,tau_r,tau_s]` over the nonnegative duration.
    /// Outputs are work-conjugate discrete stator currents and absolute shaft
    /// speeds. Rotor heat uses the retained discrete rotor current, including
    /// standstill, motoring and generating. Returned heat is the dissipation
    /// of the native discrete step, not an exact-flow integration claim.
    /// A refused/cancelled step publishes neither caller state nor outputs.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 7],
        inputs: &[f64; 4],
        duration_s: f64,
        next: &mut [f64; 7],
        outputs: &mut [f64; 4],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<InductionStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let (mut candidate, mut effort_ports) = ([0.0; 7], [0.0; 4]);
        let native = self.workspace.step_into_controlled(
            &self.model.system,
            initial,
            inputs,
            duration_s,
            &mut candidate,
            &mut effort_ports,
            &mut cancelled,
        )?;
        let effort = self
            .workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("induction accepted constitutive effort is unavailable"))?;
        let p = self.model.parameters;
        let electrical_work_j = [
            duration_s * effort_ports[0] * inputs[0],
            duration_s * effort_ports[1] * inputs[1],
        ];
        let shaft_angle_increment_rad =
            [duration_s * effort_ports[2], duration_s * effort_ports[3]];
        let shaft_work_j = [
            shaft_angle_increment_rad[0] * inputs[2],
            shaft_angle_increment_rad[1] * inputs[3],
        ];
        let stator_copper_heat_j = duration_s
            * (effort[0] * p.stator_resistance_ohm * effort[0]
                + effort[1] * p.stator_resistance_ohm * effort[1]);
        let rotor_copper_heat_j = duration_s
            * (effort[2] * p.rotor_resistance_ohm * effort[2]
                + effort[3] * p.rotor_resistance_ohm * effort[3]);
        let relative_speed = effort[5] - effort[6];
        let bearing_heat_j =
            duration_s * relative_speed * p.bearing_drag_nm_s_per_rad * relative_speed;
        let electromagnetic_angular_impulse_nm_s = -duration_s * effort[4];
        let bearing_angular_impulse_nm_s =
            -duration_s * p.bearing_drag_nm_s_per_rad * relative_speed;
        if electrical_work_j
            .iter()
            .chain(&shaft_work_j)
            .chain(&shaft_angle_increment_rad)
            .chain(&[
                stator_copper_heat_j,
                rotor_copper_heat_j,
                bearing_heat_j,
                electromagnetic_angular_impulse_nm_s,
                bearing_angular_impulse_nm_s,
            ])
            .any(|v| !v.is_finite())
        {
            return Err(
                invalid("induction work, heat or angular impulse is unrepresentable").into(),
            );
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        outputs.copy_from_slice(&effort_ports);
        Ok(InductionStep {
            native,
            electrical_work_j,
            shaft_work_j,
            stator_copper_heat_j,
            rotor_copper_heat_j,
            bearing_heat_j,
            electromagnetic_angular_impulse_nm_s,
            bearing_angular_impulse_nm_s,
            shaft_angle_increment_rad,
        })
    }
}

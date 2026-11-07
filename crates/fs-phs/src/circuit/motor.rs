//! Reciprocal commutated motor/generator with a finite rotor AND stator.
//!
//! Original states are `[lambda, p_rotor, p_stator]`, with flux linkage
//! `lambda = L*i` and absolute axial angular momenta `p = J*omega`. One signed
//! SI coefficient k owns both `emf = k*(omega_rotor-omega_stator)` and the
//! electromagnetic rotor torque `k*i`. The stator receives its exact opposite.
//! With voltage and external rotor/stator torques entering their ports,
//! `lambda' = V - R*i - emf`, `p_rotor' = k*i - b*omega_relative + tau_rotor`,
//! and `p_stator' = -k*i + b*omega_relative + tau_stator`.
//!
//! Thus electrical-to-mechanical conversion is a skew power exchange, bearing
//! drag is internal, and only external torques change total angular momentum.
//! Copper and bearing heat are the separate positive trajectory integrals of
//! `R*i^2` and `b*omega_relative^2`. Regeneration needs no sign switch, clipping,
//! efficiency multiplier or separate torque/current integration.
//!
//! This is a coaxial, frozen-parameter commutated DC machine, with linear
//! inductance and viscous rotor-to-stator friction. It resolves electrical and
//! both inertial storage, not commutation ripple, saturation, iron loss,
//! brushes, a thermal law or 3D moving axes. Data and applicability are supplied
//! by the caller. A supported stator needs an actual external reaction or a
//! mechanical constraint; it is never silently pinned or assigned infinite J.
//!
//! SI reciprocity and the armature law: MIT MAS.865, Motors and Actuators,
//! <https://fab.cba.mit.edu/classes/MAS.865/topics/motors_and_actuators/index.html>.

use crate::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactZoh, PreparedSourceLoss,
};
use crate::{
    PhsError, PortHamiltonian, PreparedStepError, PreparedStepRecord, QuadraticStorage,
    ResistanceWitness,
};
use fs_math::det;

/// Bounded nonlinear winding storage on the same reciprocal three-state motor.
pub mod nonlinear;

/// Caller-admitted, frozen constitutive data in SI units. No fitted defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DcMotorParameters {
    pub inductance_h: f64,
    pub winding_resistance_ohm: f64,
    pub rotor_inertia_kg_m2: f64,
    pub stator_inertia_kg_m2: f64,
    /// Torque/current [N m/A], also back EMF/relative speed [V s/rad].
    /// Its sign selects winding orientation; zero explicitly decouples it.
    pub coupling_nm_per_a: f64,
    /// Internal viscous bearing coefficient [N m s/rad], not drag to ground.
    pub bearing_drag_nm_s_per_rad: f64,
}

/// Native three-state motor/generator with one electrical and two torque ports.
/// Inputs are `[V, tau_rotor, tau_stator]`; outputs are `[i, omega_r, omega_s]`.
/// All efforts and flows are signed; both shaft speeds are absolute.
pub struct DcMotor {
    parameters: DcMotorParameters,
    system: PortHamiltonian,
}

/// One interval on the SAME electrical/mechanical trajectory. Heat entries
/// describe the native dissipation, never additional debits to stored energy.
#[derive(Clone, Debug)]
pub struct MotorStep {
    pub native: PreparedStepRecord,
    /// Signed electrical work INTO the armature [J]; negative is regeneration.
    pub electrical_work_j: f64,
    /// Signed external work INTO rotor and stator, respectively [J].
    pub shaft_work_j: [f64; 2],
    pub copper_heat_j: f64,
    pub bearing_heat_j: f64,
    /// Mean electromagnetic torque ON the rotor [N m]; stator gets its negative.
    pub mean_electromagnetic_torque_nm: f64,
    /// Mean viscous bearing torque ON the rotor [N m]; stator gets its negative.
    pub mean_bearing_torque_nm: f64,
    /// Mean armature back EMF [V], using relative rotor/stator speed.
    pub mean_back_emf_v: f64,
    /// Integrated absolute rotor and stator speeds [rad], not endpoint speeds.
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

impl DcMotor {
    /// Assemble the original energy `lambda^2/(2L) + p_r^2/(2J_r) + p_s^2/(2J_s)`.
    /// The same signed k enters both sides of the electrical/mechanical exchange.
    /// Bearing source squares retain their exact common-rotation loss nullspace.
    ///
    /// # Errors
    /// Refuses nonpositive/unrepresentable storage, negative/nonfinite losses,
    /// nonfinite coupling, or the original skew/PSD admission failure.
    pub fn new(parameters: DcMotorParameters) -> Result<Self, PhsError> {
        let storage = [
            parameters.inductance_h,
            parameters.rotor_inertia_kg_m2,
            parameters.stator_inertia_kg_m2,
        ];
        if storage.iter().any(|&value| {
            !value.is_finite() || value <= 0.0 || !(1.0 / value).is_finite() || 1.0 / value <= 0.0
        }) {
            return Err(invalid(
                "motor requires finite positive representable L and inertias",
            ));
        }
        let r = parameters.winding_resistance_ohm;
        let b = parameters.bearing_drag_nm_s_per_rad;
        let k = parameters.coupling_nm_per_a;
        if [r, b]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "motor requires nonnegative finite winding and bearing losses",
            ));
        }
        if !k.is_finite() {
            return Err(invalid(
                "motor requires one finite signed reciprocal coefficient",
            ));
        }
        let resistance = vec![r, 0.0, 0.0, 0.0, b, -b, 0.0, -b, b];
        let witness = ResistanceWitness::try_from_source_rows(
            3,
            resistance.clone(),
            vec![det::sqrt(r), 0.0, 0.0, 0.0, det::sqrt(b), -det::sqrt(b)],
        )?;
        let system = PortHamiltonian::new(
            3,
            3,
            vec![0.0, -k, k, k, 0.0, 0.0, -k, 0.0, 0.0],
            resistance,
            vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            Box::new(QuadraticStorage::new(
                vec![
                    1.0 / storage[0],
                    0.0,
                    0.0,
                    0.0,
                    1.0 / storage[1],
                    0.0,
                    0.0,
                    0.0,
                    1.0 / storage[2],
                ],
                3,
            )?),
        )?
        .with_resistance_witness(witness)?;
        Ok(Self { parameters, system })
    }

    pub const fn parameters(&self) -> DcMotorParameters {
        self.parameters
    }

    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Move this SAME owner into native circuit/mechanical interconnection.
    /// For a converter bus (current input, voltage output), connect its bus to
    /// motor port 0 with `interconnect_mapped` and map `[1]`: the motor receives
    /// bus voltage and its current leaves the bus, including regenerative flow.
    /// The composite exact solve retains all flux, charge and momentum states.
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }

    /// Initial current [A] and absolute rotor/stator speeds [rad/s] to original
    /// flux/momenta. No relative-speed substitution or sign projection occurs.
    pub fn state(&self, current_a: f64, shaft_speed_rad_s: [f64; 2]) -> Result<[f64; 3], PhsError> {
        if [current_a, shaft_speed_rad_s[0], shaft_speed_rad_s[1]]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "motor initial current and shaft speeds must be finite",
            ));
        }
        let state = [
            self.parameters.inductance_h * current_a,
            self.parameters.rotor_inertia_kg_m2 * shaft_speed_rad_s[0],
            self.parameters.stator_inertia_kg_m2 * shaft_speed_rad_s[1],
        ];
        if state.iter().any(|value| !value.is_finite())
            || !self.system.hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "motor initial flux, momenta or energy is unrepresentable",
            ));
        }
        Ok(state)
    }

    /// Cache one positive-duration exact held-voltage/held-torque interval.
    /// Copper and bearing source squares are integrated independently on this
    /// unchanged generator, not inferred from energy loss or fixed fractions.
    /// `heat_budget` bounds each existing six-coordinate source preparation.
    /// Repeated successful calls allocate nothing and compute no exponential.
    pub fn prepare_controlled(
        self,
        duration_s: f64,
        budget: LinearSplitBudget,
        heat_budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDcMotor, PreparedStepError> {
        poll(&mut cancelled)?;
        let p = self.parameters;
        let inverse_storage = [
            1.0 / p.inductance_h,
            1.0 / p.rotor_inertia_kg_m2,
            1.0 / p.stator_inertia_kg_m2,
        ];
        let (j, r, g) = self.system.structure();
        // Q is the original diagonal storage Hessian. These observers own no
        // state and advance no clock: A = (J-R)Q and B = G of this same system.
        let flow: [f64; 9] = core::array::from_fn(|i| (j[i] - r[i]) * inverse_storage[i % 3]);
        let copper = PreparedSourceLoss::new_controlled(
            3,
            3,
            &flow,
            g,
            &[
                det::sqrt(p.winding_resistance_ohm) * inverse_storage[0],
                0.0,
                0.0,
            ],
            duration_s,
            heat_budget,
            &mut cancelled,
        )?;
        let bearing = PreparedSourceLoss::new_controlled(
            3,
            3,
            &flow,
            g,
            &[
                0.0,
                det::sqrt(p.bearing_drag_nm_s_per_rad) * inverse_storage[1],
                -det::sqrt(p.bearing_drag_nm_s_per_rad) * inverse_storage[2],
            ],
            duration_s,
            heat_budget,
            &mut cancelled,
        )?;
        let native = PreparedExactZoh::new_owned(self.system, duration_s, budget)?;
        poll(&mut cancelled)?;
        Ok(PreparedDcMotor {
            parameters: p,
            native,
            copper,
            bearing,
            heat_scratch: [0.0; 6],
        })
    }
}

/// Cached simultaneous electrical and two-shaft evolution. Accepted state and
/// angle history belong to the caller. Refused or cancelled trials are retryable.
pub struct PreparedDcMotor {
    parameters: DcMotorParameters,
    native: PreparedExactZoh<'static>,
    copper: PreparedSourceLoss,
    bearing: PreparedSourceLoss,
    heat_scratch: [f64; 6],
}

impl PreparedDcMotor {
    pub fn system(&self) -> &PortHamiltonian {
        self.native.system()
    }

    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }

    /// Hold `[voltage, external rotor torque, external stator torque]` over
    /// the interval. Outputs are actual means `[current, rotor speed, stator
    /// speed]`; their products with held inputs give physical port work.
    /// Native energy admission, both positive heat channels, finite reports and
    /// final cancellation all precede publication of either caller buffer.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        inputs: &[f64; 3],
        next: &mut [f64; 3],
        average_outputs: &mut [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MotorStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let (mut candidate, mut mean) = ([0.0; 3], [0.0; 3]);
        let native = self
            .native
            .step_into(initial, inputs, &mut candidate, &mut mean)?;
        let copper_heat_j = self.copper.evaluate_into_controlled(
            initial,
            inputs,
            &mut self.heat_scratch,
            &mut cancelled,
        )?;
        let bearing_heat_j = self.bearing.evaluate_into_controlled(
            initial,
            inputs,
            &mut self.heat_scratch,
            &mut cancelled,
        )?;
        let dt = self.native.duration();
        let electrical_work_j = (dt * mean[0]) * inputs[0];
        let shaft_angle_increment_rad = [dt * mean[1], dt * mean[2]];
        let shaft_work_j = core::array::from_fn(|i| shaft_angle_increment_rad[i] * inputs[i + 1]);
        let relative_speed = mean[1] - mean[2];
        let mean_electromagnetic_torque_nm = self.parameters.coupling_nm_per_a * mean[0];
        let mean_bearing_torque_nm = -self.parameters.bearing_drag_nm_s_per_rad * relative_speed;
        let mean_back_emf_v = self.parameters.coupling_nm_per_a * relative_speed;
        if [
            electrical_work_j,
            shaft_work_j[0],
            shaft_work_j[1],
            shaft_angle_increment_rad[0],
            shaft_angle_increment_rad[1],
            copper_heat_j,
            bearing_heat_j,
            mean_electromagnetic_torque_nm,
            mean_bearing_torque_nm,
            mean_back_emf_v,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(
                invalid("motor physical work, heat or torque report is unrepresentable").into(),
            );
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        average_outputs.copy_from_slice(&mean);
        Ok(MotorStep {
            native,
            electrical_work_j,
            shaft_work_j,
            copper_heat_j,
            bearing_heat_j,
            mean_electromagnetic_torque_nm,
            mean_bearing_torque_nm,
            mean_back_emf_v,
            shaft_angle_increment_rad,
        })
    }
}

//! Finite-inductance voltage drive coupled to original rotor motion.
//!
//! The winding carries its current between steps, storing `L I² / 2`. A held
//! terminal voltage advances `L I' + R I + Kt omega = V` together with the
//! original full mechanical midpoint equation, using `tau_mid = Kt I_mid`.
//! Torque and back-EMF share the SAME SI constant. A scalar Schur complement
//! caches the mechanical response to one ampere; warm steps need no new matrix
//! factorization, fixed-point loop or nonlinear iteration.
//!
//! Mechanical work and `R I_mid² dt` are MIDPOINT QUADRATURE. Winding heat is
//! not an exact RMS-current integral, including on the linearly reconstructed
//! current segment. The discrete electrical, magnetic and mechanical energy
//! balance is checked independently on the returned states. Endpoint current
//! and source-power limits bound the declared linear-current reconstruction,
//! not unresolved continuous-time peaks. Directional source work and bridge
//! heat split at reconstructed current reversals; zero net terminal work does
//! not erase bridge losses. Temporal accuracy still requires refinement.
//!
//! Motor parameters, nominal rotor spin and M/C/K/G are held at their declared
//! operating conditions. Saturation, thermal parameter changes, bridge switch
//! transients, source storage and quiescent losses are not inferred.

use super::{
    Dd, MAX_ROTOR_DOFS, RotorCurrentDrive, RotorDriveError, RotorDrivePowerControl, RotorError,
    RotorMidpoint, RotorPowerControl, RotorPowerError, RotorTransientState, RotorTransientStep,
    directional_average, finite, midpoint, pairing, poll, power_fraction, voltage_interval,
};

/// Constant reciprocal motor parameters plus strictly positive winding storage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorWindingMotor {
    /// Static motor/bridge parameters, current/voltage capabilities and declared
    /// speed domain. This owner replaces that model's held-current assumption.
    pub drive: RotorCurrentDrive,
    /// Positive winding inductance [H]; magnetic storage is L I² / 2.
    pub inductance_h: f64,
}

/// Continuing physical state: original rotor coordinates and winding current.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorWindingState {
    /// Original full-coordinate rotor state.
    pub rotor: RotorTransientState,
    /// Signed winding current [A]. Magnetic energy derives from this state
    /// and the declared inductance, rather than an independent mutable balance.
    pub current_a: f64,
}

/// Source capability and explicit allowances for a coupled voltage step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorVoltageControl {
    /// Source power limits and numerical power allowance. Its voltage allowance
    /// also gates the independently evaluated discrete winding equation [V].
    pub source: RotorDrivePowerControl,
    /// Absolute allowance for initial/final current-capability excess [A].
    pub absolute_current_tolerance_a: f64,
}

/// Cached electrical coupling borrowing an already prepared native rotor step.
#[derive(Debug, Clone)]
pub struct RotorVoltageMotor<'step, 'model> {
    stepper: &'step RotorMidpoint<'model>,
    motor: RotorWindingMotor,
    back_emf_impulse_per_speed: Dd,
    feedback_impedance_time: Dd,
    current_increment_denominator: Dd,
}

/// Accepted coupled motion with independent magnetic, source and heat accounts.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorVoltageStep {
    /// Native motion under uncurtailed external loads and actual midpoint torque.
    pub motion: RotorTransientStep,
    /// Continuing winding current [A], obtained from the coupled finite step.
    pub winding_current_after_a: f64,
    /// Fraction in [0,1] of the original signed voltage command; one for zero.
    pub voltage_fraction: f64,
    /// Held winding terminal voltage after voltage/current/source limiting [V].
    pub applied_voltage_v: f64,
    /// Mean of the original and returned winding currents [A].
    pub midpoint_current_a: f64,
    /// Actual Kt I_mid used by the original mechanical midpoint equation [N m].
    pub midpoint_torque_nm: f64,
    /// Initial signed electrical source power [W].
    pub initial_source_power_w: f64,
    /// Final signed electrical source power [W].
    pub final_source_power_w: f64,
    /// Nonnegative source draw integrated over reconstructed positive power [J].
    pub electrical_draw_j: f64,
    /// Nonnegative source return integrated over reconstructed negative power [J].
    pub electrical_return_j: f64,
    /// Initial winding magnetic energy L I_before² / 2 [J].
    pub magnetic_energy_before_j: f64,
    /// Final winding magnetic energy L I_after² / 2 [J].
    pub magnetic_energy_after_j: f64,
    /// Polarized change L I_mid (I_after - I_before) from actual states [J].
    pub magnetic_energy_change_j: f64,
    /// Midpoint torque/velocity work on the original shaft coordinate [J].
    pub shaft_work_j: f64,
    /// Independently paired midpoint work of external generalized loads [J].
    pub external_work_j: f64,
    /// Nonnegative midpoint Joule quadrature R I_mid² dt [J]; not exact RMS heat.
    pub winding_heat_j: f64,
    /// Nonnegative bridge conversion/braking heat, separate from winding heat [J].
    pub bridge_heat_j: f64,
    /// Held voltage minus inductive, resistive and back-EMF midpoint terms [V].
    pub winding_voltage_residual_v: f64,
    /// Native load work minus independently paired external and shaft work [J].
    pub load_partition_residual_j: f64,
    /// Draw minus return, magnetic change, shaft work and both heat accounts [J].
    pub electrical_balance_residual_j: f64,
}

impl RotorVoltageStep {
    /// Move the accepted physical state into the next step without cloning it.
    #[must_use]
    pub fn into_state(self) -> RotorWindingState {
        RotorWindingState {
            rotor: self.motion.state_after,
            current_a: self.winding_current_after_a,
        }
    }
}

/// A failed coupled step publishes no partial current, motion or source debit.
#[derive(Debug, Clone, PartialEq)]
pub enum RotorVoltageError {
    /// Existing motor, source, dynamics, energy or cancellation refusal.
    Drive(RotorDriveError),
    /// Invalid inductance, winding state, control or scalar representation.
    Invalid(&'static str),
    /// The original coupled midpoint equations have a zero scalar Schur pivot.
    SingularCoupling,
    /// No uniformly curtailed voltage command meets the final current limit.
    NoCurrentFeasibility,
    /// Actual original/returned current exceeds its declared allowance.
    CurrentCompliance {
        /// Original winding current [A].
        initial_current_a: f64,
        /// Returned or initially proposed winding current [A].
        final_current_a: f64,
        /// Positive current-magnitude excess above the declared capability [A].
        excess_a: f64,
        /// Authored numerical allowance [A].
        tolerance_a: f64,
    },
    /// Actual returned current/motion does not satisfy the winding equation.
    WindingBalance {
        /// Independent voltage residual [V].
        residual_v: f64,
        /// Authored numerical allowance [V].
        tolerance_v: f64,
    },
}

impl From<RotorDriveError> for RotorVoltageError {
    fn from(error: RotorDriveError) -> Self {
        Self::Drive(error)
    }
}

impl From<RotorError> for RotorVoltageError {
    fn from(error: RotorError) -> Self {
        Self::Drive(error.into())
    }
}

impl From<RotorPowerError> for RotorVoltageError {
    fn from(error: RotorPowerError) -> Self {
        Self::Drive(error.into())
    }
}

impl std::fmt::Display for RotorVoltageError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Drive(error) => write!(f, "{error}"),
            Self::Invalid(reason) => write!(f, "invalid rotor winding drive: {reason}"),
            Self::SingularCoupling => {
                write!(f, "singular original rotor/winding midpoint coupling")
            }
            Self::NoCurrentFeasibility => write!(f, "no current-compliant rotor voltage command"),
            Self::CurrentCompliance {
                excess_a,
                tolerance_a,
                ..
            } => write!(
                f,
                "rotor winding current excess {excess_a} A exceeds {tolerance_a} A"
            ),
            Self::WindingBalance {
                residual_v,
                tolerance_v,
            } => write!(
                f,
                "rotor winding voltage residual {residual_v} V exceeds {tolerance_v} V"
            ),
        }
    }
}

impl std::error::Error for RotorVoltageError {}

impl<'model> RotorMidpoint<'model> {
    /// Cache a finite-inductance motor's scalar coupling to this original rotor.
    ///
    /// One solve through the existing mechanical factors obtains the shaft
    /// velocity increment per ampere. No matrix is refactored or inverted.
    /// Negative active-bearing response is retained; a zero Schur pivot refuses
    /// rather than regularizing the original coupled physical equations.
    ///
    /// # Errors
    /// Invalid motor/inductance, singular or unrepresentable coupling, original
    /// native refusal, or cancellation. The borrowed rotor remains unchanged.
    pub fn voltage_motor<'step>(
        &'step self,
        motor: RotorWindingMotor,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorVoltageMotor<'step, 'model>, RotorVoltageError> {
        poll(&mut cancelled)?;
        motor.drive.admit(
            self.model.n,
            RotorDrivePowerControl {
                maximum_source_input_w: 0.0,
                maximum_source_return_w: 0.0,
                absolute_power_tolerance_w: 0.0,
                absolute_voltage_tolerance_v: 0.0,
            },
            0.0,
        )?;
        if !motor.inductance_h.is_finite() || motor.inductance_h <= 0.0 {
            return Err(RotorVoltageError::Invalid(
                "finite positive winding inductance required",
            ));
        }
        let shaft = motor.drive.shaft_coordinate;
        let h = self.duration_seconds;
        let mut current_response = [0.0; MAX_ROTOR_DOFS];
        current_response[shaft] = ((h * motor.drive.torque_constant_nm_per_a)
            / self.mass_scale[shaft])
            / self.step_row_scale[shaft];
        self.solve_velocity_increment(&mut current_response[..self.model.n], &mut cancelled)?;
        let back_emf_impulse_per_speed =
            Dd::from_f64(h) * Dd::from_f64(motor.drive.torque_constant_nm_per_a);
        let feedback_impedance_time = Dd::from_f64(h)
            * Dd::from_f64(motor.drive.winding_resistance_ohm)
            + Dd::from_f64(0.5)
                * back_emf_impulse_per_speed
                * Dd::from_f64(current_response[shaft]);
        let current_increment_denominator =
            Dd::from_f64(motor.inductance_h) + Dd::from_f64(0.5) * feedback_impedance_time;
        for value in [
            back_emf_impulse_per_speed,
            feedback_impedance_time,
            current_increment_denominator,
        ] {
            finite(value.hi)?;
            finite(value.lo)?;
        }
        if current_increment_denominator == Dd::ZERO {
            return Err(RotorVoltageError::SingularCoupling);
        }
        poll(&mut cancelled)?;
        Ok(RotorVoltageMotor {
            stepper: self,
            motor,
            back_emf_impulse_per_speed,
            feedback_impedance_time,
            current_increment_denominator,
        })
    }
}

impl RotorVoltageMotor<'_, '_> {
    /// Advance current and rotor motion under a held, uniformly curtailed voltage.
    ///
    /// The original command first meets the terminal-voltage capability. The
    /// coupled response is affine in the remaining voltage fraction: current
    /// bounds give one interval, initial terminal power is linear, and final
    /// terminal power is quadratic. Stable existing power roots select the
    /// largest source-feasible fraction under that interval's upper bound.
    /// A source cutback that makes winding current infeasible REFUSES; it never
    /// clips stored current or discards its magnetic energy. External loads
    /// remain uncurtailed. The returned actual current and original native
    /// motion independently admit voltage, source power and energy balance.
    ///
    /// # Errors
    /// Invalid inputs, cancellation, native refusal, no feasible current/source
    /// operating point, unsupported speed, or failed winding/energy admission.
    /// The borrowed mechanical and electrical state is unchanged on all paths.
    pub fn step(
        &self,
        state: &RotorWindingState,
        external_force: &[f64],
        commanded_voltage_v: f64,
        control: RotorVoltageControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorVoltageStep, RotorVoltageError> {
        poll(&mut cancelled)?;
        let stepper = self.stepper;
        let drive = self.motor.drive;
        stepper.admit_state(&state.rotor)?;
        stepper.admit_vector(external_force)?;
        drive.admit(stepper.model.n, control.source, commanded_voltage_v)?;
        if !state.current_a.is_finite()
            || !control.absolute_current_tolerance_a.is_finite()
            || control.absolute_current_tolerance_a < 0.0
        {
            return Err(RotorVoltageError::Invalid(
                "finite winding state and nonnegative current allowance required",
            ));
        }
        self.admit_current(
            state.current_a,
            state.current_a,
            control.absolute_current_tolerance_a,
        )?;
        let shaft = drive.shaft_coordinate;
        drive.admit_speed(state.rotor.velocity[shaft])?;
        let h = stepper.duration_seconds;
        let bounded_voltage = commanded_voltage_v.clamp(
            -drive.maximum_terminal_voltage_v,
            drive.maximum_terminal_voltage_v,
        );
        let free_increment =
            stepper.velocity_increment(&state.rotor, external_force, &mut cancelled)?;
        let free_mid_speed = Dd::from_f64(state.rotor.velocity[shaft])
            + Dd::from_f64(0.5) * Dd::from_f64(free_increment[shaft]);
        // Eliminate the mechanical response without subtracting a large
        // unchanged L I term: solve directly for the winding-current increment.
        let free_rhs = -self.back_emf_impulse_per_speed * free_mid_speed
            - self.feedback_impedance_time * Dd::from_f64(state.current_a);
        let free_current_increment = free_rhs / self.current_increment_denominator;
        let voltage_current_increment =
            (Dd::from_f64(h) * Dd::from_f64(bounded_voltage)) / self.current_increment_denominator;
        let final_current_constant = Dd::from_f64(state.current_a) + free_current_increment;
        let initial_power = Dd::from_f64(bounded_voltage) * Dd::from_f64(state.current_a);
        let final_power_linear = Dd::from_f64(bounded_voltage) * final_current_constant;
        let final_power_quadratic = Dd::from_f64(bounded_voltage) * voltage_current_increment;
        for value in [
            free_rhs,
            free_current_increment,
            voltage_current_increment,
            final_current_constant,
            initial_power,
            final_power_linear,
            final_power_quadratic,
        ] {
            finite(value.hi)?;
            finite(value.lo)?;
        }
        let mut minimum_fraction = 0.0;
        let mut fraction = 1.0;
        // The shared affine limiter is dimension independent; expose its
        // feasibility failure as CURRENT rather than voltage in this owner.
        voltage_interval(
            final_current_constant,
            voltage_current_increment,
            drive.maximum_current_a,
            &mut minimum_fraction,
            &mut fraction,
        )
        .map_err(|error| match error {
            RotorDriveError::NoVoltageFeasibility => RotorVoltageError::NoCurrentFeasibility,
            other => RotorVoltageError::Drive(other),
        })?;
        let terminal_limits = RotorPowerControl {
            maximum_motoring_power_w: control.source.maximum_source_input_w
                * drive.motoring_efficiency,
            maximum_regenerative_power_w: if drive.regenerative_efficiency == 0.0 {
                f64::MAX
            } else {
                (control.source.maximum_source_return_w / drive.regenerative_efficiency)
                    .min(f64::MAX)
            },
            absolute_power_tolerance_w: control.source.absolute_power_tolerance_w,
        };
        let cap = Dd::from_f64(fraction);
        fraction *= power_fraction(Dd::ZERO, initial_power * cap, terminal_limits)?;
        let cap = Dd::from_f64(fraction);
        fraction *= power_fraction(
            final_power_quadratic * cap * cap,
            final_power_linear * cap,
            terminal_limits,
        )?;
        let applied_voltage_v = finite(bounded_voltage * fraction)?;
        let current_increment = (free_rhs + Dd::from_f64(h) * Dd::from_f64(applied_voltage_v))
            / self.current_increment_denominator;
        let winding_current_after_a =
            finite((Dd::from_f64(state.current_a) + current_increment).to_f64())?;
        self.admit_current(
            state.current_a,
            winding_current_after_a,
            control.absolute_current_tolerance_a,
        )?;
        let midpoint_current_a = finite(0.5 * state.current_a + 0.5 * winding_current_after_a)?;
        let midpoint_torque_nm = finite(drive.torque_constant_nm_per_a * midpoint_current_a)?;
        let n = stepper.model.n;
        let mut applied_force = [0.0; MAX_ROTOR_DOFS];
        applied_force[..n].copy_from_slice(external_force);
        applied_force[shaft] = finite(applied_force[shaft] + midpoint_torque_nm)?;
        let motion = stepper.step(&state.rotor, &applied_force[..n], &mut cancelled)?;
        drive.admit_speed(motion.state_after.velocity[shaft])?;
        let mid_velocity = midpoint(&state.rotor.velocity, &motion.state_after.velocity);
        let actual_current_increment =
            Dd::from_f64(winding_current_after_a) - Dd::from_f64(state.current_a);
        let winding_voltage_residual_v = finite(
            ((Dd::from_f64(h) * Dd::from_f64(applied_voltage_v)
                - Dd::from_f64(h)
                    * Dd::from_f64(drive.winding_resistance_ohm)
                    * Dd::from_f64(midpoint_current_a)
                - self.back_emf_impulse_per_speed * Dd::from_f64(mid_velocity[shaft])
                - Dd::from_f64(self.motor.inductance_h) * actual_current_increment)
                / Dd::from_f64(h))
            .to_f64(),
        )?;
        if winding_voltage_residual_v.abs() > control.source.absolute_voltage_tolerance_v {
            return Err(RotorVoltageError::WindingBalance {
                residual_v: winding_voltage_residual_v,
                tolerance_v: control.source.absolute_voltage_tolerance_v,
            });
        }
        let initial_terminal_power = finite(applied_voltage_v * state.current_a)?;
        let final_terminal_power = finite(applied_voltage_v * winding_current_after_a)?;
        let initial_source_power_w = drive.source_power(initial_terminal_power)?;
        let final_source_power_w = drive.source_power(final_terminal_power)?;
        let (positive_terminal, negative_terminal) =
            directional_average(initial_terminal_power, final_terminal_power);
        let draw_w = finite(positive_terminal / drive.motoring_efficiency)?;
        let return_w = finite(negative_terminal * drive.regenerative_efficiency)?;
        let electrical_draw_j = finite(h * draw_w)?;
        let electrical_return_j = finite(h * return_w)?;
        let winding_heat_j = finite(
            h * finite((drive.winding_resistance_ohm * midpoint_current_a) * midpoint_current_a)?,
        )?;
        let bridge_heat_j = finite(
            h * finite(
                (Dd::from_f64(draw_w) * Dd::from_f64(1.0 - drive.motoring_efficiency)
                    + Dd::from_f64(negative_terminal)
                        * Dd::from_f64(1.0 - drive.regenerative_efficiency))
                .to_f64(),
            )?,
        )?;
        let magnetic_energy_before_j = self.magnetic_energy(state.current_a)?;
        let magnetic_energy_after_j = self.magnetic_energy(winding_current_after_a)?;
        let magnetic_energy_change_j = finite(
            (Dd::from_f64(self.motor.inductance_h)
                * Dd::from_f64(midpoint_current_a)
                * actual_current_increment)
                .to_f64(),
        )?;
        let shaft_work_j = finite(
            (Dd::from_f64(h)
                * Dd::from_f64(midpoint_torque_nm)
                * Dd::from_f64(mid_velocity[shaft]))
            .to_f64(),
        )?;
        let mut external_impulse = [0.0; MAX_ROTOR_DOFS];
        for i in 0..n {
            external_impulse[i] = finite(h * external_force[i])?;
        }
        let external_work_j = pairing(&external_impulse[..n], &mid_velocity[..n])?;
        let load_partition_residual_j = finite(
            (Dd::from_f64(motion.input_work_j)
                - Dd::from_f64(external_work_j)
                - Dd::from_f64(shaft_work_j))
            .to_f64(),
        )?;
        let electrical_balance_residual_j = finite(
            (Dd::from_f64(electrical_draw_j)
                - Dd::from_f64(electrical_return_j)
                - Dd::from_f64(magnetic_energy_change_j)
                - Dd::from_f64(shaft_work_j)
                - Dd::from_f64(winding_heat_j)
                - Dd::from_f64(bridge_heat_j))
            .to_f64(),
        )?;
        let source_excess_w = [initial_source_power_w, final_source_power_w]
            .into_iter()
            .map(|p| {
                (p - control.source.maximum_source_input_w)
                    .max(-p - control.source.maximum_source_return_w)
            })
            .fold(0.0, f64::max);
        let native_residual_w = finite(motion.work_energy_residual_j / h)?;
        let load_partition_residual_w = finite(load_partition_residual_j / h)?;
        let electrical_residual_w = finite(electrical_balance_residual_j / h)?;
        if [
            source_excess_w,
            native_residual_w.abs(),
            load_partition_residual_w.abs(),
            electrical_residual_w.abs(),
        ]
        .into_iter()
        .any(|x| x > control.source.absolute_power_tolerance_w)
        {
            return Err(RotorDriveError::PowerBalance {
                source_excess_w,
                native_residual_w,
                load_partition_residual_w,
                electrical_residual_w,
                tolerance_w: control.source.absolute_power_tolerance_w,
            }
            .into());
        }
        poll(&mut cancelled)?;
        Ok(RotorVoltageStep {
            motion,
            winding_current_after_a,
            voltage_fraction: if commanded_voltage_v == 0.0 {
                1.0
            } else {
                (bounded_voltage / commanded_voltage_v) * fraction
            },
            applied_voltage_v,
            midpoint_current_a,
            midpoint_torque_nm,
            initial_source_power_w,
            final_source_power_w,
            electrical_draw_j,
            electrical_return_j,
            magnetic_energy_before_j,
            magnetic_energy_after_j,
            magnetic_energy_change_j,
            shaft_work_j,
            external_work_j,
            winding_heat_j,
            bridge_heat_j,
            winding_voltage_residual_v,
            load_partition_residual_j,
            electrical_balance_residual_j,
        })
    }

    fn admit_current(
        &self,
        initial: f64,
        final_value: f64,
        tolerance_a: f64,
    ) -> Result<(), RotorVoltageError> {
        let excess_a = initial.abs().max(final_value.abs()) - self.motor.drive.maximum_current_a;
        if excess_a > tolerance_a {
            return Err(RotorVoltageError::CurrentCompliance {
                initial_current_a: initial,
                final_current_a: final_value,
                excess_a,
                tolerance_a,
            });
        }
        Ok(())
    }

    fn magnetic_energy(&self, current: f64) -> Result<f64, RotorError> {
        finite(
            (Dd::from_f64(self.motor.inductance_h)
                * Dd::from_f64(0.5 * current)
                * Dd::from_f64(current))
            .to_f64(),
        )
    }
}

//! Current-regulated electrical actuation of one original rotor coordinate.
//!
//! The caller identifies a rotational coordinate in radians and supplies an
//! explicit constant-parameter motor/bridge model. SI reciprocity uses the SAME
//! constant for torque `tau = Kt I` and back-EMF `e = Kt omega`. Winding voltage
//! is `V = e + R I`; winding heat is independently evaluated as `R I²`. Bridge
//! efficiency follows the sign of TERMINAL power `V I`, which can be positive
//! even while the shaft is braking. A holding torque therefore draws power and
//! heats the winding at stall when R is nonzero.
//!
//! Current is held over one midpoint step: winding inductance, current ramps,
//! magnetic saturation, quiescent loss, source storage and thermal feedback
//! are not inferred. The explicit speed domain, current capability, terminal
//! voltage compliance and source power limits apply to this declared rung.
//! The prepared rotor retains its nominal spin and M/C/K/G coefficients. An
//! evolving nominal spin requires the caller to prepare the appropriate rotor
//! model; this actuator does not infer a gyroscopic or thermal material update.
//! The source/voltage peak guarantee concerns the native linear-velocity
//! reconstruction; unresolved continuous-time peaks still require refinement.

use super::{
    Dd, MAX_ROTOR_DOFS, RotorError, RotorMidpoint, RotorPowerControl, RotorPowerError,
    RotorTransientState, RotorTransientStep, finite, midpoint, pairing, poll, power_fraction,
};

pub mod voltage;

/// Authored, reciprocal, zero-inductance motor and bidirectional bridge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorCurrentDrive {
    /// Original coordinate explicitly declared to be a shaft angle [rad].
    pub shaft_coordinate: usize,
    /// Positive torque constant [N m/A], also back-EMF constant [V s/rad].
    pub torque_constant_nm_per_a: f64,
    /// Nonnegative winding resistance at the declared operating conditions [ohm].
    pub winding_resistance_ohm: f64,
    /// Nonnegative magnitude of permitted held winding current [A].
    pub maximum_current_a: f64,
    /// Nonnegative available magnitude of winding terminal voltage [V].
    pub maximum_terminal_voltage_v: f64,
    /// Inclusive lower speed where these constant parameters are valid [rad/s].
    pub minimum_speed_rad_s: f64,
    /// Inclusive upper speed where these constant parameters are valid [rad/s].
    pub maximum_speed_rad_s: f64,
    /// Terminal input divided by source input when terminal power is positive,
    /// in (0,1]. Winding loss is already included in terminal input.
    pub motoring_efficiency: f64,
    /// Source return divided by recovered terminal power, in [0,1]. Zero
    /// explicitly declares a dissipative brake path with no source recovery.
    pub regenerative_efficiency: f64,
}

/// Electrical supply capability and explicit numerical acceptance allowances.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorDrivePowerControl {
    /// Nonnegative peak electrical draw from the external source [W].
    pub maximum_source_input_w: f64,
    /// Nonnegative peak electrical power returned to the source [W].
    pub maximum_source_return_w: f64,
    /// Absolute allowance [W] for source-limit excess, original mechanical
    /// work residual, load partition and independent electrical energy balance.
    pub absolute_power_tolerance_w: f64,
    /// Absolute allowance for final terminal-voltage limit excess [V].
    pub absolute_voltage_tolerance_v: f64,
}

/// Accepted original motion with separate electrical, shaft and heat accounts.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorCurrentDriveStep {
    /// Original full-coordinate motion under external loads plus actual torque.
    pub motion: RotorTransientStep,
    /// Fraction in [0,1] of the original current command; one for a zero command.
    pub current_fraction: f64,
    /// Actual held signed winding current [A].
    pub current_a: f64,
    /// Actual held torque applied to the declared shaft coordinate [N m].
    pub torque_nm: f64,
    /// Initial winding terminal voltage, evaluated as Kt omega + R I [V].
    pub initial_terminal_voltage_v: f64,
    /// Final winding terminal voltage from the returned native velocity [V].
    pub final_terminal_voltage_v: f64,
    /// Initial signed power drawn from the electrical source [W].
    pub initial_source_power_w: f64,
    /// Final signed power drawn from the electrical source [W].
    pub final_source_power_w: f64,
    /// Signed torque work on the original shaft coordinate [J].
    pub shaft_work_j: f64,
    /// Independently paired work of the supplied external generalized loads [J].
    pub external_work_j: f64,
    /// Nonnegative electrical energy drawn, integrating only positive power [J].
    pub electrical_draw_j: f64,
    /// Nonnegative energy returned, integrating only negative source power [J].
    pub electrical_return_j: f64,
    /// Nonnegative winding heat R I² dt, including holding loss at stall [J].
    pub winding_heat_j: f64,
    /// Nonnegative conversion/braking heat, separate from winding heat [J].
    pub bridge_heat_j: f64,
    /// Native total load work minus independently paired external/shaft work [J].
    pub load_partition_residual_j: f64,
    /// Draw minus return, shaft work, winding heat and bridge heat [J].
    pub electrical_balance_residual_j: f64,
}

/// A refused drive step publishes neither motion nor an electrical debit.
#[derive(Debug, Clone, PartialEq)]
pub enum RotorDriveError {
    /// Original dynamics, power-root or cancellation refusal.
    Power(RotorPowerError),
    /// Invalid declared motor, source, state or scalar representation.
    Invalid(&'static str),
    /// An endpoint is outside the authored constant-parameter speed domain.
    SpeedDomain {
        /// Refused shaft speed [rad/s].
        speed_rad_s: f64,
        /// Inclusive declared lower speed [rad/s].
        minimum_speed_rad_s: f64,
        /// Inclusive declared upper speed [rad/s].
        maximum_speed_rad_s: f64,
    },
    /// No uniformly curtailed command satisfies both endpoint voltage bounds.
    NoVoltageFeasibility,
    /// Actual accepted-step voltages exceed the declared voltage allowance.
    VoltageCompliance {
        /// Actual initial voltage [V].
        initial_voltage_v: f64,
        /// Actual final voltage [V].
        final_voltage_v: f64,
        /// Excess above the available voltage magnitude [V].
        excess_v: f64,
        /// Authored numerical allowance [V].
        tolerance_v: f64,
    },
    /// Independent final power or energy accounting exceeds its allowance.
    PowerBalance {
        /// Positive source draw or return limit excess [W].
        source_excess_w: f64,
        /// Native full-coordinate mechanical work residual per second [W].
        native_residual_w: f64,
        /// Load partition residual per second [W].
        load_partition_residual_w: f64,
        /// Electrical balance residual per second [W].
        electrical_residual_w: f64,
        /// Authored numerical allowance [W].
        tolerance_w: f64,
    },
}

impl From<RotorError> for RotorDriveError {
    fn from(error: RotorError) -> Self {
        Self::Power(RotorPowerError::Dynamics(error))
    }
}

impl From<RotorPowerError> for RotorDriveError {
    fn from(error: RotorPowerError) -> Self {
        Self::Power(error)
    }
}

impl std::fmt::Display for RotorDriveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Power(error) => write!(f, "{error}"),
            Self::Invalid(reason) => write!(f, "invalid rotor electrical drive: {reason}"),
            Self::SpeedDomain { speed_rad_s, .. } => {
                write!(
                    f,
                    "rotor electrical drive speed {speed_rad_s} is outside its declared domain"
                )
            }
            Self::NoVoltageFeasibility => write!(f, "no voltage-compliant rotor drive effort"),
            Self::VoltageCompliance {
                excess_v,
                tolerance_v,
                ..
            } => write!(
                f,
                "rotor drive voltage excess {excess_v} V exceeds {tolerance_v} V"
            ),
            Self::PowerBalance { tolerance_w, .. } => write!(
                f,
                "rotor electrical drive power/energy balance exceeds {tolerance_w} W"
            ),
        }
    }
}

impl std::error::Error for RotorDriveError {}

impl RotorMidpoint<'_> {
    /// Drive one declared shaft with a held, current-limited electrical command.
    ///
    /// External loads retain their complete original coordinates and are never
    /// curtailed. The motor current is first bounded by its authored capability.
    /// Two cached native response solves then give affine endpoint voltages and
    /// quadratic terminal powers versus the remaining current fraction. Voltage
    /// bounds intersect to one interval. Four bounded analytic power projections
    /// find the largest source-feasible fraction within its upper bound, even
    /// when active bearings create disconnected feasible sets. The original
    /// native step is solved once with the selected physical torque: no new
    /// factorization or dynamic nonlinear iteration is introduced; cost is O(n²).
    ///
    /// Actual voltage/current and torque/velocity pairings independently supply
    /// the returned source and shaft accounts. Positive and negative terminal
    /// work are integrated separately across a reversal, so net-zero mechanical
    /// work cannot conceal winding or bridge heat. Speed-domain, voltage, source
    /// and independent energy admission precede the final cancellation poll.
    ///
    /// # Errors
    /// Invalid model/control, native refusal, cancellation, an unavailable
    /// voltage/current operating point, or a speed/power/energy refusal. Borrowed
    /// state, external loads and model remain unchanged on every path.
    pub fn step_current_drive(
        &self,
        state: &RotorTransientState,
        external_force: &[f64],
        commanded_current_a: f64,
        drive: RotorCurrentDrive,
        control: RotorDrivePowerControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorCurrentDriveStep, RotorDriveError> {
        poll(&mut cancelled)?;
        self.admit_state(state)?;
        self.admit_vector(external_force)?;
        drive.admit(self.model.n, control, commanded_current_a)?;
        let shaft = drive.shaft_coordinate;
        drive.admit_speed(state.velocity[shaft])?;
        let bounded_current =
            commanded_current_a.clamp(-drive.maximum_current_a, drive.maximum_current_a);
        let bounded_torque = finite(drive.torque_constant_nm_per_a * bounded_current)?;
        let n = self.model.n;
        let free_increment = self.velocity_increment(state, external_force, &mut cancelled)?;
        let mut force_increment = [0.0; MAX_ROTOR_DOFS];
        force_increment[shaft] = ((self.duration_seconds * bounded_torque)
            / self.mass_scale[shaft])
            / self.step_row_scale[shaft];
        self.solve_velocity_increment(&mut force_increment[..n], &mut cancelled)?;

        let kt = Dd::from_f64(drive.torque_constant_nm_per_a);
        let current = Dd::from_f64(bounded_current);
        let resistive_voltage = Dd::from_f64(drive.winding_resistance_ohm) * current;
        let voltages = [
            (kt * Dd::from_f64(state.velocity[shaft]), resistive_voltage),
            (
                kt * (Dd::from_f64(state.velocity[shaft]) + Dd::from_f64(free_increment[shaft])),
                kt * Dd::from_f64(force_increment[shaft]) + resistive_voltage,
            ),
        ];
        let mut minimum_fraction = 0.0;
        let mut fraction = 1.0;
        for &(constant, slope) in &voltages {
            for value in [constant, slope, current * constant, current * slope] {
                finite(value.hi)?;
                finite(value.lo)?;
            }
            voltage_interval(
                constant,
                slope,
                drive.maximum_terminal_voltage_v,
                &mut minimum_fraction,
                &mut fraction,
            )?;
        }
        let terminal_limits = RotorPowerControl {
            maximum_motoring_power_w: control.maximum_source_input_w * drive.motoring_efficiency,
            maximum_regenerative_power_w: if drive.regenerative_efficiency == 0.0 {
                f64::MAX
            } else {
                (control.maximum_source_return_w / drive.regenerative_efficiency).min(f64::MAX)
            },
            absolute_power_tolerance_w: control.absolute_power_tolerance_w,
        };
        // Each quadratic's feasible set contains zero and at most two intervals.
        // A second pass suffices: after leaving an upper interval, every further
        // reduction remains in that constraint's lower interval. A single pass
        // could land in the first endpoint's forbidden regeneration interval.
        for _ in 0..2 {
            for &(constant, slope) in &voltages {
                let cap = Dd::from_f64(fraction);
                fraction *= power_fraction(
                    current * slope * cap * cap,
                    current * constant * cap,
                    terminal_limits,
                )?;
            }
        }
        let current_a = finite(bounded_current * fraction)?;
        let torque_nm = finite(drive.torque_constant_nm_per_a * current_a)?;
        let mut applied = [0.0; MAX_ROTOR_DOFS];
        applied[..n].copy_from_slice(external_force);
        applied[shaft] = finite(applied[shaft] + torque_nm)?;
        let motion = self.step(state, &applied[..n], &mut cancelled)?;
        drive.admit_speed(motion.state_after.velocity[shaft])?;
        let initial_terminal_voltage_v = drive.voltage(state.velocity[shaft], current_a)?;
        let final_terminal_voltage_v =
            drive.voltage(motion.state_after.velocity[shaft], current_a)?;
        let excess_v = initial_terminal_voltage_v
            .abs()
            .max(final_terminal_voltage_v.abs())
            - drive.maximum_terminal_voltage_v;
        if excess_v > control.absolute_voltage_tolerance_v {
            return Err(RotorDriveError::VoltageCompliance {
                initial_voltage_v: initial_terminal_voltage_v,
                final_voltage_v: final_terminal_voltage_v,
                excess_v,
                tolerance_v: control.absolute_voltage_tolerance_v,
            });
        }
        let initial_terminal_power = finite(initial_terminal_voltage_v * current_a)?;
        let final_terminal_power = finite(final_terminal_voltage_v * current_a)?;
        let initial_source_power_w = drive.source_power(initial_terminal_power)?;
        let final_source_power_w = drive.source_power(final_terminal_power)?;
        let (positive_terminal, negative_terminal) =
            directional_average(initial_terminal_power, final_terminal_power);
        let draw_w = finite(positive_terminal / drive.motoring_efficiency)?;
        let return_w = finite(negative_terminal * drive.regenerative_efficiency)?;
        let h = self.duration_seconds;
        let electrical_draw_j = finite(h * draw_w)?;
        let electrical_return_j = finite(h * return_w)?;
        let winding_heat_j =
            finite(h * finite((drive.winding_resistance_ohm * current_a) * current_a)?)?;
        let bridge_heat_j = finite(
            h * finite(
                (Dd::from_f64(draw_w) * Dd::from_f64(1.0 - drive.motoring_efficiency)
                    + Dd::from_f64(negative_terminal)
                        * Dd::from_f64(1.0 - drive.regenerative_efficiency))
                .to_f64(),
            )?,
        )?;
        let mid_velocity = midpoint(&state.velocity, &motion.state_after.velocity);
        let shaft_work_j = finite(
            (Dd::from_f64(h) * Dd::from_f64(torque_nm) * Dd::from_f64(mid_velocity[shaft]))
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
                - Dd::from_f64(shaft_work_j)
                - Dd::from_f64(winding_heat_j)
                - Dd::from_f64(bridge_heat_j))
            .to_f64(),
        )?;
        let source_excess_w = [initial_source_power_w, final_source_power_w]
            .into_iter()
            .map(|p| (p - control.maximum_source_input_w).max(-p - control.maximum_source_return_w))
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
        .any(|x| x > control.absolute_power_tolerance_w)
        {
            return Err(RotorDriveError::PowerBalance {
                source_excess_w,
                native_residual_w,
                load_partition_residual_w,
                electrical_residual_w,
                tolerance_w: control.absolute_power_tolerance_w,
            });
        }
        poll(&mut cancelled)?;
        Ok(RotorCurrentDriveStep {
            motion,
            current_fraction: if commanded_current_a == 0.0 {
                1.0
            } else {
                (bounded_current / commanded_current_a) * fraction
            },
            current_a,
            torque_nm,
            initial_terminal_voltage_v,
            final_terminal_voltage_v,
            initial_source_power_w,
            final_source_power_w,
            shaft_work_j,
            external_work_j,
            electrical_draw_j,
            electrical_return_j,
            winding_heat_j,
            bridge_heat_j,
            load_partition_residual_j,
            electrical_balance_residual_j,
        })
    }
}

impl RotorCurrentDrive {
    fn admit(
        self,
        n: usize,
        control: RotorDrivePowerControl,
        command: f64,
    ) -> Result<(), RotorDriveError> {
        if self.shaft_coordinate >= n
            || !command.is_finite()
            || !self.torque_constant_nm_per_a.is_finite()
            || self.torque_constant_nm_per_a <= 0.0
            || !self.minimum_speed_rad_s.is_finite()
            || !self.maximum_speed_rad_s.is_finite()
            || self.minimum_speed_rad_s > self.maximum_speed_rad_s
            || !self.motoring_efficiency.is_finite()
            || !(0.0 < self.motoring_efficiency && self.motoring_efficiency <= 1.0)
            || !self.regenerative_efficiency.is_finite()
            || !(0.0..=1.0).contains(&self.regenerative_efficiency)
            || [
                self.winding_resistance_ohm,
                self.maximum_current_a,
                self.maximum_terminal_voltage_v,
                control.maximum_source_input_w,
                control.maximum_source_return_w,
                control.absolute_power_tolerance_w,
                control.absolute_voltage_tolerance_v,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
        {
            return Err(RotorDriveError::Invalid(
                "finite reciprocal motor parameters, declared speed domain and nonnegative capabilities/allowances required",
            ));
        }
        Ok(())
    }

    fn admit_speed(self, speed_rad_s: f64) -> Result<(), RotorDriveError> {
        if speed_rad_s < self.minimum_speed_rad_s || speed_rad_s > self.maximum_speed_rad_s {
            return Err(RotorDriveError::SpeedDomain {
                speed_rad_s,
                minimum_speed_rad_s: self.minimum_speed_rad_s,
                maximum_speed_rad_s: self.maximum_speed_rad_s,
            });
        }
        Ok(())
    }

    fn voltage(self, speed: f64, current: f64) -> Result<f64, RotorError> {
        finite(
            (Dd::from_f64(self.torque_constant_nm_per_a) * Dd::from_f64(speed)
                + Dd::from_f64(self.winding_resistance_ohm) * Dd::from_f64(current))
            .to_f64(),
        )
    }

    fn source_power(self, terminal_power: f64) -> Result<f64, RotorError> {
        finite(if terminal_power >= 0.0 {
            terminal_power / self.motoring_efficiency
        } else {
            terminal_power * self.regenerative_efficiency
        })
    }
}

// Intersect an affine endpoint voltage bound with the already admitted interval.
// Scale before taking a root; only divide for a crossing known to lie inside
// the current interval, avoiding irrelevant huge/overflowing roots.
fn voltage_interval(
    constant: Dd,
    slope: Dd,
    limit: f64,
    minimum: &mut f64,
    maximum: &mut f64,
) -> Result<(), RotorDriveError> {
    let scale = constant.abs().hi.max(slope.abs().hi).max(limit);
    if scale == 0.0 {
        return Ok(());
    }
    let scale = Dd::from_f64(scale);
    let c = constant / scale;
    let s = slope / scale;
    let l = Dd::from_f64(limit) / scale;
    let at_min = c + s * Dd::from_f64(*minimum);
    let at_max = c + s * Dd::from_f64(*maximum);
    let increasing = !s.lt(Dd::ZERO);
    let (low, high) = if increasing {
        (at_min, at_max)
    } else {
        (at_max, at_min)
    };
    if high.lt(-l) || l.lt(low) {
        return Err(RotorDriveError::NoVoltageFeasibility);
    }
    if increasing {
        if at_min.lt(-l) {
            *minimum = finite(((-l - c) / s).to_f64())?.clamp(*minimum, *maximum);
        }
        if l.lt(at_max) {
            *maximum = finite(((l - c) / s).to_f64())?.clamp(*minimum, *maximum);
        }
    } else {
        if l.lt(at_min) {
            *minimum = finite(((l - c) / s).to_f64())?.clamp(*minimum, *maximum);
        }
        if at_max.lt(-l) {
            *maximum = finite(((-l - c) / s).to_f64())?.clamp(*minimum, *maximum);
        }
    }
    Ok(())
}

// Exact positive/negative areas for the declared affine power segment, divided
// by its duration. Normalize the crossing time without squaring a power or
// overflowing the sum of opposite endpoint magnitudes.
fn directional_average(initial: f64, final_value: f64) -> (f64, f64) {
    if initial >= 0.0 && final_value >= 0.0 {
        return (0.5 * initial + 0.5 * final_value, 0.0);
    }
    if initial <= 0.0 && final_value <= 0.0 {
        return (0.0, -0.5 * initial - 0.5 * final_value);
    }
    let positive = initial.max(final_value);
    let negative = -initial.min(final_value);
    let scale = positive.max(negative);
    let sum = positive / scale + negative / scale;
    (
        0.5 * positive * ((positive / scale) / sum),
        0.5 * negative * ((negative / scale) / sum),
    )
}

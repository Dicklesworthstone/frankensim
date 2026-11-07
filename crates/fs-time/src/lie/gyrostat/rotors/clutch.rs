//! Finite-capacity dry shaft friction with resolved slip/hold transitions.
//!
//! A reciprocal port has slip v = b.q' and shaft torque -b^T*effort.
//! With attitude fixed, v' = a - H*effort, where
//! H = b (diag(1/J) + A^T I^-1 A) b^T > 0. Held drives therefore give
//! piecewise constant acceleration: one slip-to-zero event suffices before
//! either static hold or opposite-sense sliding. The original rotor bank
//! owns every simultaneous impulse, its actual work and carrier recoil.

use super::{
    HeldTorque, LieStepError, RotorBank, RotorState, RotorStep, components, dd, dot_dd, increment,
    invalid, product, scalar, subtract3, vector, work_residual,
};
use fs_ga::{So3, Vec3};
use fs_math::dd::Dd;

const ROUNDING: f64 = 8192.0 * f64::EPSILON;

/// Original relative-shaft coordinates on which dry friction acts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DryClutchPort {
    /// Shaft relative to the carrier; positive effort brakes positive spin.
    Brake { rotor: usize },
    /// Slip is q'_first - ratio*q'_second. Torques are (-effort, ratio*effort).
    /// The shafts must be distinct; the finite ratio may have either sign.
    Shafts {
        first_rotor: usize,
        second_rotor: usize,
        ratio: f64,
    },
}

impl DryClutchPort {
    fn coefficient(self, rotor: usize) -> f64 {
        match self {
            Self::Brake { rotor: selected } => {
                if rotor == selected {
                    1.0
                } else {
                    0.0
                }
            }
            Self::Shafts {
                first_rotor,
                second_rotor,
                ratio,
            } => {
                if rotor == first_rotor {
                    1.0
                } else if rotor == second_rotor {
                    -ratio
                } else {
                    0.0
                }
            }
        }
    }

    fn mapped<const N: usize>(self, values: &[f64; N]) -> Dd {
        match self {
            Self::Brake { rotor } => dd(values[rotor]),
            Self::Shafts {
                first_rotor,
                second_rotor,
                ratio,
            } => dd(values[first_rotor]) - dd(ratio) * dd(values[second_rotor]),
        }
    }

    // Integrate the ONE physical angle of an already-admitted sticking pair.
    // This never projects canonical momenta or captures a nonzero initial slip.
    fn sticking_phase_rates<const N: usize>(self, rates: &mut [f64; N]) {
        match self {
            Self::Brake { rotor } => rates[rotor] = 0.0,
            Self::Shafts {
                first_rotor,
                second_rotor,
                ratio,
            } => {
                if ratio.abs() <= 1.0 {
                    rates[first_rotor] = ratio * rates[second_rotor];
                } else {
                    rates[second_rotor] = rates[first_rotor] / ratio;
                }
            }
        }
    }
}

/// Caller-authored isothermal dry friction, with capacities conjugate to slip.
///
/// Static capacity must be at least kinetic capacity, both finite and
/// nonnegative. Sliding effort is kinetic_capacity*sign(slip); static effort
/// is solved within the static interval. Both zero explicitly release the
/// clutch. No coefficient, contact area, normal force, temperature dependence
/// or thermal destination is inferred from these torque capacities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryClutchLaw {
    pub port: DryClutchPort,
    pub static_capacity_nm: f64,
    pub kinetic_capacity_nm: f64,
}

/// Retained physical branch, without a user-selected velocity capture window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DryClutchMode {
    /// Both authored capacities are zero; no reaction is applied.
    Open,
    /// A zero-slip state or an accepted synchronization event permits hold.
    Sticking,
    /// Sliding, or zero slip with load above the static capacity.
    Slipping,
}

/// Native mechanical state plus its explicit dry-friction branch.
///
/// Use `DryClutchStepper::start` after an externally applied mechanical change.
/// Retaining this value across steps preserves a proved synchronization despite
/// the native floating-point endpoint remainder; it does not add a deadband.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryClutchState<const N: usize> {
    rotors: RotorState<N>,
    port: DryClutchPort,
    mode: DryClutchMode,
}

impl<const N: usize> DryClutchState<N> {
    pub const fn rotors(&self) -> &RotorState<N> {
        &self.rotors
    }
    pub const fn mode(&self) -> DryClutchMode {
        self.mode
    }
    /// Release or transfer this mechanical state without an impulse or reset.
    pub const fn into_rotors(self) -> RotorState<N> {
        self.rotors
    }
}

/// Event data for one frozen-attitude shaft half-stage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryClutchStage {
    /// Slip-to-zero time measured from THIS half-stage's start, in seconds.
    /// None means no such event; this is not a certified full-trajectory event.
    pub synchronization_time_s: Option<f64>,
    pub sticking_duration_s: f64,
    pub mode_after: DryClutchMode,
}

/// Native motion, signed source work and physical dry-friction dissipation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryClutchStep<const N: usize> {
    pub motion: RotorStep<N>,
    pub mode_after: DryClutchMode,
    pub drive_work_j: [f64; N],
    /// Signed friction work on each original shaft; its sum is minus heat.
    pub friction_work_j: [f64; N],
    /// Nonnegative integrated kinetic-friction heat; static reactions add none.
    pub dissipated_energy_j: f64,
    /// Integral of signed effort. Port torques are -b^T times this value.
    pub friction_impulse_nms: f64,
    /// Integral of absolute slip speed during the two shaft half-stages.
    pub sliding_angle_rad: f64,
    /// Exact scalar events of the two subflows, in execution order.
    pub half_stages: [DryClutchStage; 2],
    /// Delta kinetic + heat - drive - external work - native drift defect.
    pub work_energy_residual_j: f64,
    port: DryClutchPort,
}

impl<const N: usize> DryClutchStep<N> {
    pub const fn state_after(&self) -> DryClutchState<N> {
        DryClutchState {
            rotors: self.motion.state_after,
            port: self.port,
            mode: self.mode_after,
        }
    }
}

/// One prepared finite-capacity clutch/brake on a native full-inertia bank.
///
/// A step has at most four constant-torque pieces and one original gyrostat
/// drift. All work is O(N), bounded by the bank's 64-rotor limit, with fixed
/// arrays, no allocation, no iteration and no independent pairwise kicks.
/// Held-load scalar subflows resolve stopping and breakaway analytically.
/// Smooth sliding is a second-order split approximation; general precessing
/// stick/slip transitions require refinement and have no asserted order or
/// exact event-time bound. Precession-created slip is re-evaluated explicitly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryClutchStepper<const N: usize> {
    bank: RotorBank<N>,
    law: DryClutchLaw,
    duration_s: f64,
    mobility: Dd,
    inverse_axis: Vec3,
}

impl<const N: usize> RotorBank<N> {
    pub fn dry_clutch_stepper(
        &self,
        law: DryClutchLaw,
        duration_s: f64,
    ) -> Result<DryClutchStepper<N>, LieStepError> {
        scalar(duration_s, "dry clutch duration")?;
        if duration_s <= 0.0 || 0.5 * duration_s == 0.0 {
            return Err(invalid("positive resolved dry clutch duration"));
        }
        match law.port {
            DryClutchPort::Brake { rotor } if rotor < N => {}
            DryClutchPort::Shafts {
                first_rotor,
                second_rotor,
                ratio,
            } if first_rotor < N && second_rotor < N && first_rotor != second_rotor => {
                scalar(ratio, "dry clutch transmission ratio")?;
            }
            _ => return Err(invalid("original dry clutch or brake shafts")),
        }
        let static_capacity = scalar(law.static_capacity_nm, "dry clutch static capacity")?;
        let kinetic_capacity = scalar(law.kinetic_capacity_nm, "dry clutch kinetic capacity")?;
        if kinetic_capacity < 0.0 || static_capacity < kinetic_capacity {
            return Err(invalid("static capacity >= kinetic capacity >= zero"));
        }
        let mut axis = [dd(0.0); 3];
        let mut mobility = dd(0.0);
        for (i, rotor) in self.rotors.iter().enumerate() {
            let b = dd(law.port.coefficient(i));
            mobility = mobility + b * b / dd(rotor.axial_inertia_kg_m2);
            for (component, a) in axis.iter_mut().zip(components(rotor.axis_body)) {
                *component = *component + b * dd(a);
            }
        }
        let axis = vector(
            Vec3::new(axis[0].to_f64(), axis[1].to_f64(), axis[2].to_f64()),
            "dry clutch carrier reaction axis",
        )?;
        let inverse_axis = self.inertia.velocity(axis)?;
        let mobility = mobility + dot_dd(axis, inverse_axis);
        if scalar(mobility.to_f64(), "dry clutch effective mobility")? <= 0.0 {
            return Err(invalid("positive dry clutch mobility"));
        }
        Ok(DryClutchStepper {
            bank: *self,
            law,
            duration_s,
            mobility,
            inverse_axis,
        })
    }
}

impl<const N: usize> DryClutchStepper<N> {
    pub const fn law(&self) -> DryClutchLaw {
        self.law
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Change in slip rate per unit reciprocal port impulse, [1/(kg m^2)].
    pub const fn mobility(&self) -> f64 {
        self.mobility.to_f64()
    }

    /// Bind an unmodified native state. Every nonzero initial slip slides;
    /// initialization never uses the sticking admission tolerance to capture it.
    pub fn start(&self, rotors: RotorState<N>) -> Result<DryClutchState<N>, LieStepError> {
        let slip = self.slip(&rotors)?.0;
        Ok(DryClutchState {
            rotors,
            port: self.law.port,
            mode: if self.law.static_capacity_nm == 0.0 {
                DryClutchMode::Open
            } else if slip == 0.0 {
                DryClutchMode::Sticking
            } else {
                DryClutchMode::Slipping
            },
        })
    }

    pub fn step(
        &self,
        state: &DryClutchState<N>,
        drive_torque_nm: [f64; N],
        external_torque: HeldTorque,
    ) -> Result<DryClutchStep<N>, LieStepError> {
        if state.port != self.law.port {
            return Err(invalid("retained dry clutch port changed"));
        }
        if self.law.static_capacity_nm == 0.0 {
            let motion = self.bank.step(
                &state.rotors,
                drive_torque_nm,
                external_torque,
                self.duration_s,
            )?;
            return Ok(DryClutchStep {
                drive_work_j: motion.shaft_work_j,
                friction_work_j: [0.0; N],
                dissipated_energy_j: 0.0,
                friction_impulse_nms: 0.0,
                sliding_angle_rad: 0.0,
                half_stages: [DryClutchStage {
                    synchronization_time_s: None,
                    sticking_duration_s: 0.0,
                    mode_after: DryClutchMode::Open,
                }; 2],
                work_energy_residual_j: motion.work_energy_residual_j,
                motion,
                mode_after: DryClutchMode::Open,
                port: self.law.port,
            });
        }
        let first = self.stage(
            &state.rotors,
            state.mode,
            drive_torque_nm,
            self.body_torque(state.rotors.rotation, external_torque)?,
        )?;
        let mut drifted = first.state;
        (drifted.rotation, drifted.total_momentum_body) = self.bank.inertia.momentum_step(
            drifted.rotation,
            drifted.total_momentum_body,
            self.bank.rotor_momentum(&drifted.axial_momentum_nms)?,
            self.duration_s,
        )?;
        let drift_defect = self.bank.energy_change(&first.state, &drifted)?;
        // A prior static reaction cannot silently constrain the free native
        // drift. Preserve its rounding remainder only when actual slip has not
        // changed; otherwise the next subflow begins from the returned slip.
        let mode = if first.event.mode_after == DryClutchMode::Sticking
            && self.slip(&first.state)?.0 == self.slip(&drifted)?.0
        {
            DryClutchMode::Sticking
        } else {
            DryClutchMode::Slipping
        };
        let second = self.stage(
            &drifted,
            mode,
            drive_torque_nm,
            self.body_torque(drifted.rotation, external_torque)?,
        )?;
        let mut shaft_work = [0.0; N];
        let mut drive_work = [0.0; N];
        let mut friction_work = [0.0; N];
        for i in 0..N {
            shaft_work[i] = scalar(
                (first.shaft_work[i] + second.shaft_work[i]).to_f64(),
                "dry clutch net shaft work",
            )?;
            drive_work[i] = scalar(
                (first.drive_work[i] + second.drive_work[i]).to_f64(),
                "dry clutch authored drive work",
            )?;
            friction_work[i] = scalar(
                (first.friction_work[i] + second.friction_work[i]).to_f64(),
                "dry clutch signed friction work",
            )?;
        }
        let external_work = scalar(
            (first.external_work + second.external_work).to_f64(),
            "dry clutch external work",
        )?;
        let external_impulse = state.rotors.rotation.rotate(first.external_impulse()?)?
            + drifted.rotation.rotate(second.external_impulse()?)?;
        let change = self.bank.energy_change(&state.rotors, &second.state)?;
        let motion = RotorStep {
            state_after: second.state,
            axial_momentum_change_nms: core::array::from_fn(|i| {
                second.state.axial_momentum_nms[i] - state.rotors.axial_momentum_nms[i]
            }),
            external_angular_impulse_world: vector(external_impulse, "clutch spatial impulse")?,
            shaft_work_j: shaft_work,
            external_work_j: external_work,
            kinetic_energy_change_j: change,
            drift_energy_defect_j: drift_defect,
            work_energy_residual_j: work_residual(
                change,
                &shaft_work,
                external_work,
                drift_defect,
            )?,
            angular_momentum_residual_world: vector(
                subtract3(
                    second
                        .state
                        .rotation
                        .rotate(second.state.total_momentum_body)?,
                    state
                        .rotors
                        .rotation
                        .rotate(state.rotors.total_momentum_body)?,
                    external_impulse,
                ),
                "dry clutch spatial momentum balance",
            )?,
        };
        let heat = scalar((first.heat + second.heat).to_f64(), "dry clutch total heat")?;
        let mut balance = dd(change) + dd(heat) - dd(external_work) - dd(drift_defect);
        for work in drive_work {
            balance = balance - dd(work);
        }
        Ok(DryClutchStep {
            motion,
            mode_after: second.event.mode_after,
            drive_work_j: drive_work,
            friction_work_j: friction_work,
            dissipated_energy_j: heat,
            friction_impulse_nms: scalar(
                (first.effort_impulse + second.effort_impulse).to_f64(),
                "dry clutch integrated effort",
            )?,
            sliding_angle_rad: scalar(
                (first.sliding_angle + second.sliding_angle).to_f64(),
                "dry clutch sliding angle",
            )?,
            half_stages: [first.event, second.event],
            work_energy_residual_j: scalar(balance.to_f64(), "dry clutch complete work balance")?,
            port: self.law.port,
        })
    }

    fn body_torque(&self, rotation: So3, torque: HeldTorque) -> Result<Vec3, LieStepError> {
        match torque {
            HeldTorque::Body(value) => vector(value, "dry clutch external body torque"),
            HeldTorque::Space(value) => Ok(rotation.inverse().rotate(value)?),
        }
    }

    fn slip(&self, state: &RotorState<N>) -> Result<(f64, f64), LieStepError> {
        let rates = self.bank.relative_rates(state)?;
        let omega = self.bank.angular_velocity(state)?;
        let mut scale = dd(0.0);
        for (i, rotor) in self.bank.rotors.iter().enumerate() {
            scale = scale
                + dd(self.law.port.coefficient(i).abs())
                    * (dd((state.axial_momentum_nms[i] / rotor.axial_inertia_kg_m2).abs())
                        + dd(rotor.axis_body.dot(omega).abs()));
        }
        Ok((
            scalar(self.law.port.mapped(&rates).to_f64(), "dry clutch slip")?,
            scalar(scale.to_f64(), "dry clutch slip resolution")?,
        ))
    }

    fn rest_effort(&self, acceleration: Dd) -> Result<(f64, DryClutchMode), LieStepError> {
        let quotient = acceleration.to_f64() / self.mobility.to_f64();
        let required = if quotient.is_finite() {
            acceleration / self.mobility
        } else {
            dd(quotient)
        };
        if required.to_f64().abs() <= self.law.static_capacity_nm {
            let required = scalar(required.to_f64(), "dry clutch static reaction")?;
            if acceleration.to_f64() != 0.0 && required == 0.0 {
                return Err(invalid("unresolved dry clutch static reaction"));
            }
            Ok((required, DryClutchMode::Sticking))
        } else {
            Ok((
                acceleration.to_f64().signum() * self.law.kinetic_capacity_nm,
                DryClutchMode::Slipping,
            ))
        }
    }

    fn stage(
        &self,
        state: &RotorState<N>,
        mode: DryClutchMode,
        drive: [f64; N],
        external: Vec3,
    ) -> Result<Stage<N>, LieStepError> {
        let (actual_slip, slip_scale) = self.slip(state)?;
        let mut slip = actual_slip;
        if mode == DryClutchMode::Sticking {
            admit(slip, slip_scale, "retained dry clutch synchronization")?;
            slip = 0.0;
        }
        let mut carrier = components(external).map(dd);
        let mut acceleration = dd(0.0);
        for (i, rotor) in self.bank.rotors.iter().enumerate() {
            scalar(drive[i], "dry clutch shaft drive")?;
            acceleration = acceleration
                + dd(self.law.port.coefficient(i)) * dd(drive[i]) / dd(rotor.axial_inertia_kg_m2);
            for (component, a) in carrier.iter_mut().zip(components(rotor.axis_body)) {
                *component = *component - dd(a) * dd(drive[i]);
            }
        }
        for (value, inverse) in carrier.into_iter().zip(components(self.inverse_axis)) {
            acceleration = acceleration - value * dd(inverse);
        }
        scalar(acceleration.to_f64(), "dry clutch free acceleration")?;
        let duration = 0.5 * self.duration_s;
        let mut result = Stage::new(*state);
        if slip == 0.0 {
            let (effort, mode) = self.rest_effort(acceleration)?;
            self.piece(
                &mut result,
                drive,
                external,
                duration,
                effort,
                mode,
                0.0,
                acceleration,
            )?;
            result.event.mode_after = mode;
            if mode == DryClutchMode::Sticking {
                result.event.sticking_duration_s = duration;
            }
            return Ok(result);
        }
        let effort = slip.signum() * self.law.kinetic_capacity_nm;
        let deceleration = scalar(
            (acceleration - self.mobility * dd(effort)).to_f64(),
            "dry clutch sliding acceleration",
        )?;
        let stopping = if deceleration != 0.0 && deceleration.signum() != slip.signum() {
            let time = -slip / deceleration;
            if time.is_finite() {
                scalar(
                    (-dd(slip) / dd(deceleration)).to_f64(),
                    "dry clutch stopping time",
                )?
            } else {
                f64::INFINITY
            }
        } else {
            f64::INFINITY
        };
        if stopping <= duration {
            if stopping <= 0.0 {
                return Err(invalid("resolved positive dry clutch stopping time"));
            }
            self.piece(
                &mut result,
                drive,
                external,
                stopping,
                effort,
                DryClutchMode::Slipping,
                slip,
                acceleration,
            )?;
            let (at_stop, scale) = self.slip(&result.state)?;
            admit(
                at_stop,
                scale.max(slip.abs()),
                "resolved dry clutch synchronization",
            )?;
            result.event.synchronization_time_s = Some(stopping);
            let (hold_effort, mode) = self.rest_effort(acceleration)?;
            let remaining = scalar(duration - stopping, "dry clutch post-event duration")?;
            if remaining > 0.0 {
                self.piece(
                    &mut result,
                    drive,
                    external,
                    remaining,
                    hold_effort,
                    mode,
                    0.0,
                    acceleration,
                )?;
            }
            result.event.mode_after = mode;
            if mode == DryClutchMode::Sticking {
                result.event.sticking_duration_s = remaining;
            }
        } else {
            self.piece(
                &mut result,
                drive,
                external,
                duration,
                effort,
                DryClutchMode::Slipping,
                slip,
                acceleration,
            )?;
        }
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn piece(
        &self,
        output: &mut Stage<N>,
        drive: [f64; N],
        external: Vec3,
        duration: f64,
        effort: f64,
        mode: DryClutchMode,
        initial_slip: f64,
        free_acceleration: Dd,
    ) -> Result<(), LieStepError> {
        let mut requested = [0.0; N];
        let mut drive_impulse = [dd(0.0); N];
        let mut friction_impulse = [dd(0.0); N];
        for i in 0..N {
            product(duration, drive[i], "dry clutch drive impulse")?;
            drive_impulse[i] = dd(duration) * dd(drive[i]);
            friction_impulse[i] = -dd(duration) * dd(self.law.port.coefficient(i)) * dd(effort);
            requested[i] = scalar(
                (drive_impulse[i] + friction_impulse[i]).to_f64(),
                "dry clutch net shaft impulse",
            )?;
        }
        let external_impulse = Vec3::new(
            product(duration, external.x, "dry clutch external impulse")?,
            product(duration, external.y, "dry clutch external impulse")?,
            product(duration, external.z, "dry clutch external impulse")?,
        );
        let impulse = self
            .bank
            .apply_impulse(&output.state, requested, external_impulse)?;
        let (before_slip, before_scale) = self.slip(&output.state)?;
        let (after_slip, after_scale) = self.slip(&impulse.state_after)?;
        let expected_slip = if mode == DryClutchMode::Sticking {
            admit(before_slip, before_scale, "dry clutch sticking entry")?;
            0.0
        } else {
            scalar(
                (dd(initial_slip)
                    + dd(duration) * (free_acceleration - self.mobility * dd(effort)))
                .to_f64(),
                "dry clutch scalar flow endpoint",
            )?
        };
        admit(
            after_slip - expected_slip,
            before_scale.max(after_scale).max(expected_slip.abs()),
            "dry clutch scalar/native velocity agreement",
        )?;
        let mut rates = impulse.mean_relative_rate_rad_s;
        if mode == DryClutchMode::Sticking {
            self.law.port.sticking_phase_rates(&mut rates);
        }
        let initial_rates = self.bank.relative_rates(&output.state)?;
        let final_rates = self.bank.relative_rates(&impulse.state_after)?;
        let mut after = impulse.state_after;
        let mut balance = dd(impulse.kinetic_energy_change_j) - dd(impulse.external_work_j);
        let mut scale = impulse
            .kinetic_energy_change_j
            .abs()
            .max(impulse.external_work_j.abs());
        for i in 0..N {
            admit(
                impulse.axial_momentum_change_nms[i] - requested[i],
                drive_impulse[i]
                    .to_f64()
                    .abs()
                    .max(friction_impulse[i].to_f64().abs())
                    .max(impulse.axial_momentum_change_nms[i].abs()),
                "dry clutch original shaft torque",
            )?;
            let angle = product(duration, rates[i], "dry clutch phase increment")?;
            let reverses = (initial_rates[i] > 0.0 && final_rates[i] < 0.0)
                || (initial_rates[i] < 0.0 && final_rates[i] > 0.0);
            // A shaft can reverse DURING a sliding clutch segment while the
            // clutch's own slip keeps its sign. Its net phase may then vanish
            // by cancellation. Bound that storage residual by the two local
            // motion terms, never by the ill-conditioned net angle or |q|.
            // Monotone motion keeps the original erased-increment refusal.
            let phase_scale = if reverses {
                scalar(
                    (dd(duration)
                        * dd(0.5)
                        * (dd(initial_rates[i].abs()) + dd(final_rates[i].abs())))
                    .to_f64(),
                    "dry clutch reversing phase scale",
                )?
            } else {
                angle.abs()
            };
            after.relative_angle_rad[i] = if reverses {
                scalar(
                    output.state.relative_angle_rad[i] + angle,
                    "dry clutch relative phase",
                )?
            } else {
                increment(
                    output.state.relative_angle_rad[i],
                    angle,
                    "dry clutch relative phase",
                )?
            };
            let actual_angle =
                (dd(after.relative_angle_rad[i]) - dd(output.state.relative_angle_rad[i])).to_f64();
            admit(
                actual_angle - angle,
                phase_scale,
                "dry clutch phase resolution",
            )?;
            let drive_work = drive_impulse[i] * dd(rates[i]);
            let friction_work = friction_impulse[i] * dd(rates[i]);
            output.drive_work[i] = output.drive_work[i] + drive_work;
            output.friction_work[i] = output.friction_work[i] + friction_work;
            output.shaft_work[i] = output.shaft_work[i] + dd(impulse.shaft_work_j[i]);
            balance = balance - drive_work;
            scale = scale
                .max(drive_work.to_f64().abs())
                .max(friction_work.to_f64().abs());
        }
        let slip_angle = dd(duration) * self.law.port.mapped(&rates);
        let heat = if mode == DryClutchMode::Sticking {
            dd(0.0)
        } else {
            dd(effort) * slip_angle
        };
        let heat_value = scalar(heat.to_f64(), "dry clutch physical heat")?;
        if heat_value < 0.0 {
            return Err(invalid("nonnegative dry clutch physical heat"));
        }
        admit(
            (balance + heat).to_f64(),
            scale.max(heat_value),
            "dry clutch original work balance",
        )?;
        output.heat = output.heat + heat;
        if mode != DryClutchMode::Sticking {
            output.sliding_angle = output.sliding_angle + dd(slip_angle.to_f64().abs());
        }
        output.effort_impulse = output.effort_impulse + dd(duration) * dd(effort);
        output.external_work = output.external_work + dd(impulse.external_work_j);
        for (sum, value) in output
            .external
            .iter_mut()
            .zip(components(impulse.external_angular_impulse_body))
        {
            *sum = *sum + dd(value);
        }
        output.state = after;
        Ok(())
    }
}

struct Stage<const N: usize> {
    state: RotorState<N>,
    event: DryClutchStage,
    shaft_work: [Dd; N],
    drive_work: [Dd; N],
    friction_work: [Dd; N],
    external: [Dd; 3],
    external_work: Dd,
    heat: Dd,
    effort_impulse: Dd,
    sliding_angle: Dd,
}

impl<const N: usize> Stage<N> {
    fn new(state: RotorState<N>) -> Self {
        Self {
            state,
            event: DryClutchStage {
                synchronization_time_s: None,
                sticking_duration_s: 0.0,
                mode_after: DryClutchMode::Slipping,
            },
            shaft_work: [dd(0.0); N],
            drive_work: [dd(0.0); N],
            friction_work: [dd(0.0); N],
            external: [dd(0.0); 3],
            external_work: dd(0.0),
            heat: dd(0.0),
            effort_impulse: dd(0.0),
            sliding_angle: dd(0.0),
        }
    }
    fn external_impulse(&self) -> Result<Vec3, LieStepError> {
        vector(
            Vec3::new(
                self.external[0].to_f64(),
                self.external[1].to_f64(),
                self.external[2].to_f64(),
            ),
            "dry clutch accumulated external impulse",
        )
    }
}

fn admit(residual: f64, scale: f64, context: &'static str) -> Result<(), LieStepError> {
    scalar(residual, context)?;
    scalar(scale, context)?;
    if residual.abs() > ROUNDING * scale {
        Err(invalid(context))
    } else {
        Ok(())
    }
}

//! Finite winding storage and voltage drive on the original damped eccentric body.

use super::{EccentricRotorDriveError, ElasticEccentricCurrentMotor};
use crate::rotordynamics::RotorError;
use crate::rotordynamics::rigid::{
    DampedEccentricRotor, EccentricBearingDampingStep, EccentricRotorDampingControl,
    EccentricRotorDampingError, ElasticEccentricRotorError, ElasticEccentricRotorStep,
    RotorSpinDriveAxis,
};
use crate::rotordynamics::transient::RotorTransientState;
use crate::rotordynamics::transient::drive::{
    RotorDriveError,
    voltage::{
        RotorVoltageControl, RotorVoltageError, RotorVoltageMotor, RotorVoltageStep,
        RotorWindingMotor, RotorWindingState,
    },
};
use crate::{
    DynamicsError, MomentumImpulseReceipt, RigidBodyIntegrator, RigidBodyState,
    RigidKineticEnergyChange, Vec3,
};
use fs_math::dd::Dd;

/// Continuing physical state; magnetic storage derives from the carried current.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EccentricRotorWindingState {
    /// Sole original full-inertia rigid body.
    pub body: RigidBodyState,
    /// Signed actual winding current [A], unchanged by purely mechanical sources.
    pub current_a: f64,
}

/// Held voltage request and an independent mechanical shaft load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EccentricRotorVoltageCommand {
    /// Requested signed terminal voltage [V]; original voltage/current/power limits apply.
    pub voltage_v: f64,
    /// External torque following actual body +z [N m], separate from the motor.
    pub external_shaft_torque_nm: f64,
}

/// Original source budgets and full-body transfer/whole-trajectory allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EccentricRotorVoltageControl {
    /// Sole physical trajectory interval [s].
    pub duration_seconds: f64,
    /// Original scalar full-mass midpoint relative tolerance in (0,1).
    pub relative_tolerance: f64,
    /// Original finite-inductance current, source-power and winding-equation limits.
    pub electrical: RotorVoltageControl,
    /// Original passive-source gates; total_energy_tolerance_j gates this entire
    /// mechanical, magnetic, bearing-heat and electrical-source composition.
    pub damping: EccentricRotorDampingControl,
    /// Maximum actual-versus-commanded angular impulse component defect [N m s].
    pub angular_impulse_tolerance_nms: f64,
    /// Maximum actual-versus-scalar source endpoint shaft-rate defect [rad/s].
    pub shaft_rate_tolerance_rad_s: f64,
}

/// One original winding source applied once to the original body momentum.
#[derive(Clone, Debug, PartialEq)]
pub struct EccentricRotorVoltageSource {
    /// Native coupled current, voltage, magnetic storage, work and heat receipt.
    /// Its displacement is only a source quadrature, never a physical shaft phase.
    pub electrical: RotorVoltageStep,
    /// Original full-body impulse; pose stays fixed and no clock advances here.
    pub impulse: MomentumImpulseReceipt,
    /// Actual full-body shaft speed after this source [rad/s].
    pub shaft_rate_after_rad_s: f64,
    /// Actual rate minus the original scalar source prediction [rad/s].
    pub shaft_rate_residual_rad_s: f64,
    /// Actual native kinetic change minus independently paired motor shaft work [J].
    pub source_work_residual_j: f64,
    /// Original winding equation re-evaluated using ACTUAL body midpoint speed [V].
    pub winding_voltage_residual_v: f64,
}

/// One admitted voltage-driven damped trajectory and continuing winding current.
#[derive(Clone, Debug, PartialEq)]
pub struct EccentricRotorVoltageStep {
    /// Sole physical endpoint, including the actual carried winding current.
    pub state_after: EccentricRotorWindingState,
    /// Native inductive half sources, with the first current passed into the second.
    pub sources: [EccentricRotorVoltageSource; 2],
    /// Original passive half sources at the actual initial and final poses.
    pub damping_sources: [EccentricBearingDampingStep; 2],
    /// Original full-body elastic/gravity trajectory under only the external torque.
    pub mechanical: ElasticEccentricRotorStep,
    /// Complete original-body endpoint kinetic increment [J].
    pub kinetic_energy_change_j: f64,
    /// Polarized L*(I_before+I_after)/2*(I_after-I_before), using actual endpoints [J].
    pub magnetic_energy_change_j: f64,
    /// Independent original bearing dissipation [J].
    pub bearing_heat_j: f64,
    /// Positive electrical source draw, separate from return [J].
    pub electrical_draw_j: f64,
    /// Positive electrical energy returned [J].
    pub electrical_return_j: f64,
    /// Original R*I_mid^2 source quadrature [J], not an exact RMS-current integral.
    pub winding_heat_j: f64,
    /// Original independent bridge conversion/braking loss [J].
    pub bridge_heat_j: f64,
    /// Actual accepted terminal voltages at the initial/final endpoints [V].
    /// These include inductive voltage; they are NOT Kt*omega + R*I.
    pub endpoint_terminal_voltage_v: [f64; 2],
    /// Actual signed source power at initial/final V,I [W].
    pub endpoint_source_power_w: [f64; 2],
    /// delta(K+V_gravity+V_bearing+magnetic) + all physical heat
    /// + electrical return - draw - external shaft work [J].
    /// Signed numerical defect, never converted to heat.
    pub total_energy_residual_j: f64,
}

/// A refused operation publishes no partial body, current or source debit.
#[derive(Clone, Debug, PartialEq)]
pub enum EccentricRotorVoltageError {
    /// Invalid native interval/budget or unrepresentable response.
    Invalid(&'static str),
    /// Original winding, current, terminal-source or winding-balance refusal.
    Electrical(RotorVoltageError),
    /// Original rigid-body/elastic or electrical-source transfer refusal.
    Mechanical(EccentricRotorDriveError),
    /// Original passive bearing source refusal.
    Damping(EccentricRotorDampingError),
    /// Complete original mechanical/magnetic balance exceeds its allowance.
    EnergyBalance {
        /// Signed complete balance defect [J].
        residual_j: f64,
        /// Authored total-energy allowance [J].
        tolerance_j: f64,
    },
    /// Cancellation observed before complete publication.
    Cancelled,
}
impl From<RotorVoltageError> for EccentricRotorVoltageError {
    fn from(e: RotorVoltageError) -> Self {
        Self::Electrical(e)
    }
}
impl From<RotorDriveError> for EccentricRotorVoltageError {
    fn from(e: RotorDriveError) -> Self {
        Self::Electrical(e.into())
    }
}
impl From<RotorError> for EccentricRotorVoltageError {
    fn from(e: RotorError) -> Self {
        Self::Electrical(e.into())
    }
}
impl From<EccentricRotorDriveError> for EccentricRotorVoltageError {
    fn from(e: EccentricRotorDriveError) -> Self {
        Self::Mechanical(e)
    }
}
impl From<ElasticEccentricRotorError> for EccentricRotorVoltageError {
    fn from(e: ElasticEccentricRotorError) -> Self {
        Self::Mechanical(e.into())
    }
}
impl From<DynamicsError> for EccentricRotorVoltageError {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanical(e.into())
    }
}
impl From<EccentricRotorDampingError> for EccentricRotorVoltageError {
    fn from(e: EccentricRotorDampingError) -> Self {
        Self::Damping(e)
    }
}
impl core::fmt::Display for EccentricRotorVoltageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "eccentric rotor winding drive refused: {self:?}")
    }
}
impl std::error::Error for EccentricRotorVoltageError {}

/// Original finite-inductance source on one actual passive eccentric body.
///
/// The parent's prepared shaft-mobility model is reused unchanged. Only the
/// native electrical owner changes from held current to continuing winding
/// dynamics. The source inertia is 1/(I_body^-1)_zz, never a trajectory inertia.
#[derive(Clone, Debug)]
pub struct DampedEccentricVoltageMotor<'a> {
    rotor: &'a DampedEccentricRotor,
    source: ElasticEccentricCurrentMotor<'a>,
    motor: RotorWindingMotor,
}
impl DampedEccentricRotor {
    /// Bind explicit shaft-axis actuation to the original full-body source.
    ///
    /// Original winding/static parameter admission occurs when preparing the
    /// interval's native voltage motor. No equilibrium current is inferred.
    pub fn voltage_spin_motor(
        &self,
        motor: RotorWindingMotor,
        axis: RotorSpinDriveAxis,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DampedEccentricVoltageMotor<'_>, EccentricRotorVoltageError> {
        poll(&mut cancelled)?;
        let source =
            self.elastic_component()
                .current_spin_motor(motor.drive, axis, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(DampedEccentricVoltageMotor {
            rotor: self,
            source,
            motor,
        })
    }
}
impl DampedEccentricVoltageMotor<'_> {
    /// Original full-tensor effective inertia for the fixed-pose shaft source [kg m²].
    #[must_use]
    pub const fn source_inertia_kg_m2(&self) -> f64 {
        self.source.source_inertia_kg_m2()
    }

    /// Advance D(h/2), inductive E(h/2), elastic(h), E(h/2), D(h/2).
    ///
    /// The exact original current state passes through both E sources and into
    /// the next step; purely mechanical flows leave current/magnetic storage
    /// unchanged. The native voltage motor selects a feasible held terminal
    /// voltage and solves its original midpoint winding/mechanical equations.
    /// Actual full-body speed, impulse work and the winding equation are checked
    /// independently after applying each source once to native momentum.
    ///
    /// Terminal voltage and power use accepted V and actual I. Back-EMF may
    /// exceed terminal voltage while L*dI/dt supplies the difference: no
    /// zero-inductance compliance condition is imposed. Both outer physical
    /// endpoints admit the actual speed/current/terminal-source limits again.
    ///
    /// Smooth fixed coefficients and uncurtailed held voltage give second-order
    /// splitting when electrical, bearing and free-rotation timescales are
    /// resolved. Midpoint is not L-stable; stiff winding transients and active
    /// limiting require refinement. Endpoint gates are not a continuous peak
    /// certificate. Battery storage, switching and thermal feedback are not inferred.
    pub fn step(
        &self,
        integrator: RigidBodyIntegrator,
        state: EccentricRotorWindingState,
        command: EccentricRotorVoltageCommand,
        control: EccentricRotorVoltageControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EccentricRotorVoltageStep, EccentricRotorVoltageError> {
        poll(&mut cancelled)?;
        let h = control.duration_seconds;
        let half = 0.5 * h;
        if !half.is_finite() || half <= 0.0 || 2.0 * half != h {
            return Err(EccentricRotorVoltageError::Invalid(
                "representable positive source half interval required",
            ));
        }
        for allowance in [
            control.angular_impulse_tolerance_nms,
            control.shaft_rate_tolerance_rad_s,
        ] {
            if !allowance.is_finite() || allowance < 0.0 {
                return Err(EccentricRotorVoltageError::Invalid(
                    "finite nonnegative body transfer allowances required",
                ));
            }
        }
        finite(command.external_shaft_torque_nm)?;
        self.motor.drive.admit_speed(self.speed(state.body)?)?;
        let operating = self.source.source.at_speed(0.0, &mut cancelled)?;
        let midpoint =
            operating.midpoint_stepper(half, control.relative_tolerance, &mut cancelled)?;
        let winding = midpoint.voltage_motor(self.motor, &mut cancelled)?;
        let first_damping =
            self.rotor
                .relax_damping(state.body, half, control.damping, &mut cancelled)?;
        let first = self.apply_source(
            &winding,
            EccentricRotorWindingState {
                body: first_damping.impulse.state_after,
                ..state
            },
            command.voltage_v,
            control,
            &mut cancelled,
        )?;
        let mechanical = self.rotor.elastic_component().step(
            integrator,
            first.impulse.state_after,
            command.external_shaft_torque_nm,
            h,
            &mut cancelled,
        )?;
        let second = self.apply_source(
            &winding,
            EccentricRotorWindingState {
                body: mechanical.native.state_after,
                current_a: first.electrical.winding_current_after_a,
            },
            command.voltage_v,
            control,
            &mut cancelled,
        )?;
        let last_damping = self.rotor.relax_damping(
            second.impulse.state_after,
            half,
            control.damping,
            &mut cancelled,
        )?;
        let state_after = EccentricRotorWindingState {
            body: last_damping.impulse.state_after,
            current_a: second.electrical.winding_current_after_a,
        };
        let endpoint_voltage = [
            first.electrical.applied_voltage_v,
            second.electrical.applied_voltage_v,
        ];
        let endpoints = [state, state_after];
        let mut endpoint_power = [0.0; 2];
        for i in 0..2 {
            self.motor
                .drive
                .admit_speed(self.speed(endpoints[i].body)?)?;
            endpoint_power[i] = self.motor.drive.source_power(finite(
                (Dd::from_f64(endpoint_voltage[i]) * Dd::from_f64(endpoints[i].current_a)).to_f64(),
            )?)?;
        }
        let limits = control.electrical.source;
        let excess_v = endpoint_voltage[0].abs().max(endpoint_voltage[1].abs())
            - self.motor.drive.maximum_terminal_voltage_v;
        if excess_v > limits.absolute_voltage_tolerance_v {
            return Err(RotorDriveError::VoltageCompliance {
                initial_voltage_v: endpoint_voltage[0],
                final_voltage_v: endpoint_voltage[1],
                excess_v,
                tolerance_v: limits.absolute_voltage_tolerance_v,
            }
            .into());
        }
        let excess_a = state.current_a.abs().max(state_after.current_a.abs())
            - self.motor.drive.maximum_current_a;
        if excess_a > control.electrical.absolute_current_tolerance_a {
            return Err(RotorVoltageError::CurrentCompliance {
                initial_current_a: state.current_a,
                final_current_a: state_after.current_a,
                excess_a,
                tolerance_a: control.electrical.absolute_current_tolerance_a,
            }
            .into());
        }
        let source_excess_w = endpoint_power.into_iter().fold(0.0_f64, |e, p| {
            e.max(p - limits.maximum_source_input_w)
                .max(-p - limits.maximum_source_return_w)
        });
        if source_excess_w > limits.absolute_power_tolerance_w {
            return Err(RotorDriveError::PowerBalance {
                source_excess_w,
                native_residual_w: 0.0,
                load_partition_residual_w: 0.0,
                electrical_residual_w: 0.0,
                tolerance_w: limits.absolute_power_tolerance_w,
            }
            .into());
        }
        let sources = [first, second];
        let totals: [Dd; 4] = core::array::from_fn(|i| {
            sources.iter().fold(Dd::ZERO, |sum, s| {
                let e = &s.electrical;
                sum + Dd::from_f64(
                    [
                        e.electrical_draw_j,
                        e.electrical_return_j,
                        e.winding_heat_j,
                        e.bridge_heat_j,
                    ][i],
                )
            })
        });
        let mut kinetic = RigidKineticEnergyChange::new();
        kinetic.add_body(
            self.rotor.elastic_component().mass_properties(),
            state.body,
            state_after.body,
        )?;
        let kinetic = kinetic.finish_retained()?;
        let magnetic = Dd::from_f64(self.motor.inductance_h)
            * (Dd::from_f64(0.5 * state.current_a) + Dd::from_f64(0.5 * state_after.current_a))
            * (Dd::from_f64(state_after.current_a) - Dd::from_f64(state.current_a));
        let heat = Dd::from_f64(first_damping.dissipated_energy_j)
            + Dd::from_f64(last_damping.dissipated_energy_j);
        let residual = finite(
            (kinetic
                + magnetic
                + heat
                + Dd::from_f64(
                    mechanical
                        .native
                        .work
                        .gravitational_potential_energy_change_j,
                )
                + Dd::from_f64(mechanical.bearing_potential_change_j)
                + totals[2]
                + totals[3]
                + totals[1]
                - totals[0]
                - Dd::from_f64(mechanical.shaft_work_j))
            .to_f64(),
        )?;
        if residual.abs() > control.damping.total_energy_tolerance_j {
            return Err(EccentricRotorVoltageError::EnergyBalance {
                residual_j: residual,
                tolerance_j: control.damping.total_energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(EccentricRotorVoltageStep {
            state_after,
            sources,
            damping_sources: [first_damping, last_damping],
            mechanical,
            kinetic_energy_change_j: finite(kinetic.to_f64())?,
            magnetic_energy_change_j: finite(magnetic.to_f64())?,
            bearing_heat_j: finite(heat.to_f64())?,
            electrical_draw_j: finite(totals[0].to_f64())?,
            electrical_return_j: finite(totals[1].to_f64())?,
            winding_heat_j: finite(totals[2].to_f64())?,
            bridge_heat_j: finite(totals[3].to_f64())?,
            endpoint_terminal_voltage_v: endpoint_voltage,
            endpoint_source_power_w: endpoint_power,
            total_energy_residual_j: residual,
        })
    }

    fn speed(&self, body: RigidBodyState) -> Result<f64, EccentricRotorVoltageError> {
        Ok(self
            .rotor
            .elastic_component()
            .mass_properties()
            .angular_velocity_body_checked(body.angular_momentum_body())?
            .z)
    }
    fn apply_source(
        &self,
        winding: &RotorVoltageMotor<'_, '_>,
        state: EccentricRotorWindingState,
        voltage: f64,
        control: EccentricRotorVoltageControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EccentricRotorVoltageSource, EccentricRotorVoltageError> {
        poll(cancelled)?;
        let speed = self.speed(state.body)?;
        let electrical = winding.step(
            &RotorWindingState {
                rotor: RotorTransientState {
                    displacement: vec![0.0],
                    velocity: vec![speed],
                },
                current_a: state.current_a,
            },
            &[0.0],
            voltage,
            control.electrical,
            &mut *cancelled,
        )?;
        let duration = winding.duration_s();
        let impulse = state.body.apply_momentum_impulse(
            self.rotor.elastic_component().mass_properties(),
            Vec3::ZERO,
            Vec3::new(
                0.0,
                0.0,
                finite(
                    (Dd::from_f64(duration) * Dd::from_f64(electrical.midpoint_torque_nm)).to_f64(),
                )?,
            ),
        )?;
        let actual_speed = self.speed(impulse.state_after)?;
        self.motor.drive.admit_speed(actual_speed)?;
        let speed_error = finite(actual_speed - electrical.motion.state_after.velocity[0])?;
        let work_error = finite(
            (Dd::from_f64(impulse.work.kinetic_energy_change_j)
                - Dd::from_f64(electrical.shaft_work_j))
            .to_f64(),
        )?;
        let error = impulse.angular_residual_body;
        let impulse_error = error.x.abs().max(error.y.abs()).max(error.z.abs());
        if impulse_error > control.angular_impulse_tolerance_nms
            || speed_error.abs() > control.shaft_rate_tolerance_rad_s
            || work_error.abs()
                > finite(duration * control.electrical.source.absolute_power_tolerance_w)?
        {
            return Err(EccentricRotorDriveError::SourceMismatch {
                angular_impulse_residual_nms: impulse_error,
                shaft_rate_residual_rad_s: speed_error,
                source_work_residual_j: work_error,
            }
            .into());
        }
        let mid_current = Dd::from_f64(0.5 * state.current_a)
            + Dd::from_f64(0.5 * electrical.winding_current_after_a);
        let mid_speed = Dd::from_f64(0.5 * speed) + Dd::from_f64(0.5 * actual_speed);
        let h = Dd::from_f64(duration);
        let winding_residual = finite(
            ((h * Dd::from_f64(electrical.applied_voltage_v)
                - h * Dd::from_f64(self.motor.drive.winding_resistance_ohm) * mid_current
                - h * Dd::from_f64(self.motor.drive.torque_constant_nm_per_a) * mid_speed
                - Dd::from_f64(self.motor.inductance_h)
                    * (Dd::from_f64(electrical.winding_current_after_a)
                        - Dd::from_f64(state.current_a)))
                / h)
                .to_f64(),
        )?;
        if winding_residual.abs() > control.electrical.source.absolute_voltage_tolerance_v {
            return Err(RotorVoltageError::WindingBalance {
                residual_v: winding_residual,
                tolerance_v: control.electrical.source.absolute_voltage_tolerance_v,
            }
            .into());
        }
        poll(cancelled)?;
        Ok(EccentricRotorVoltageSource {
            electrical,
            impulse,
            shaft_rate_after_rad_s: actual_speed,
            shaft_rate_residual_rad_s: speed_error,
            source_work_residual_j: work_error,
            winding_voltage_residual_v: winding_residual,
        })
    }
}
fn finite(x: f64) -> Result<f64, EccentricRotorVoltageError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(EccentricRotorVoltageError::Invalid(
            "nonfinite winding/body response",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EccentricRotorVoltageError> {
    if cancelled() {
        Err(EccentricRotorVoltageError::Cancelled)
    } else {
        Ok(())
    }
}

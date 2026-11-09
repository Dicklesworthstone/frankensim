//! Electrical momentum sources on the original eccentric full-inertia body.

use super::{
    ElasticEccentricRotor, ElasticEccentricRotorError, ElasticEccentricRotorStep,
    RotorSpinDriveAxis,
};
use crate::rotordynamics::transient::drive::{
    RotorCurrentDrive, RotorCurrentDriveStep, RotorDriveError, RotorDrivePowerControl,
};
use crate::rotordynamics::transient::{RotorMidpoint, RotorTransientState};
use crate::rotordynamics::{RotorBearingModel, RotorError};
use crate::{
    DynamicsError, MomentumImpulseReceipt, RigidBodyIntegrator, RigidBodyState,
    RigidKineticEnergyChange, Vec3,
};
use fs_math::dd::Dd;

/// Held current request and an independently supplied mechanical shaft load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EccentricRotorCurrentCommand {
    /// Requested signed winding current [A]; native electrical limits may curtail it.
    pub current_a: f64,
    /// External torque following body +z [N m], separate from the motor and bearings.
    pub external_shaft_torque_nm: f64,
}

/// Original electrical limits plus explicit native-body transfer and energy gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EccentricRotorDriveControl {
    /// Complete physical trajectory duration [s].
    pub duration_seconds: f64,
    /// Original scalar midpoint solve relative tolerance.
    pub relative_tolerance: f64,
    /// Current-source power, voltage and independent electrical balance limits.
    pub electrical: RotorDrivePowerControl,
    /// Maximum actual-versus-commanded source angular impulse discrepancy [N m s].
    pub angular_impulse_tolerance_nms: f64,
    /// Maximum native-versus-scalar source endpoint shaft-rate discrepancy [rad/s].
    pub shaft_rate_tolerance_rad_s: f64,
    /// Maximum complete mechanical/electrical energy defect [J]. This
    /// includes elastic/native-motion truncation error, not physical heat.
    pub total_energy_tolerance_j: f64,
}

/// One electrical source interval, with the actual full-body impulse retained.
#[derive(Clone, Debug, PartialEq)]
pub struct EccentricRotorElectricalSource {
    /// Original current, voltage, back-EMF, source energy and physical heat owner.
    /// Its scalar displacement is a source-work quadrature, NOT shaft phase.
    /// The native impulse below leaves pose fixed during this source interval.
    pub electrical: RotorCurrentDriveStep,
    /// One original body momentum update, including the complete tensor storage.
    pub impulse: MomentumImpulseReceipt,
    /// Actual full-body shaft speed after the source [rad/s].
    pub shaft_rate_after_rad_s: f64,
    /// Actual full-body speed minus the original scalar source prediction [rad/s].
    pub shaft_rate_residual_rad_s: f64,
    /// Original terminal law evaluated at the actual native endpoint speed [V].
    pub terminal_voltage_after_v: f64,
    /// Original signed source power at the actual native endpoint [W].
    pub source_power_after_w: f64,
    /// Actual full-body kinetic change minus electrical shaft work [J].
    pub source_work_residual_j: f64,
}

/// Complete eccentric-body motion and separately owned electrical transfers.
#[derive(Clone, Debug, PartialEq)]
pub struct EccentricRotorCurrentStep {
    /// Sole physical endpoint. The two electrical source coordinates add no pose.
    pub state_after: RigidBodyState,
    /// Source half intervals before and after the native elastic trajectory.
    pub sources: [EccentricRotorElectricalSource; 2],
    /// Original elastic/gravity trajectory under the external mechanical torque.
    /// Its shaft work belongs to that external load, not to the electrical motor.
    pub mechanical: ElasticEccentricRotorStep,
    /// Complete original-body endpoint kinetic change, including transverse terms [J].
    pub kinetic_energy_change_j: f64,
    /// Sum of the independent electrical shaft-work accounts [J].
    pub electrical_shaft_work_j: f64,
    /// Positive source draw, never netted against regeneration [J].
    pub electrical_draw_j: f64,
    /// Positive energy returned to the source [J].
    pub electrical_return_j: f64,
    /// Original winding Joule heat over both half sources [J].
    pub winding_heat_j: f64,
    /// Original conversion/braking heat over both half sources [J].
    pub bridge_heat_j: f64,
    /// `delta(K+V_gravity+V_bearing) + heat + return - draw - external_work` [J].
    /// Signed numerical defect; never added to winding or bridge heat.
    pub total_energy_residual_j: f64,
}

/// A refused composition returns neither body motion nor an electrical debit.
#[derive(Clone, Debug, PartialEq)]
pub enum EccentricRotorDriveError {
    /// Invalid actuation axis, scalar source coordinate or numerical allowance.
    Invalid(&'static str),
    /// Original electrical source refusal.
    Electrical(RotorDriveError),
    /// Original full-body mechanical refusal.
    Mechanical(ElasticEccentricRotorError),
    /// The accepted scalar source was not resolved by the original native body.
    SourceMismatch {
        /// Maximum native angular impulse component discrepancy [N m s].
        angular_impulse_residual_nms: f64,
        /// Actual shaft speed minus source prediction [rad/s].
        shaft_rate_residual_rad_s: f64,
        /// Actual full-body kinetic change minus source shaft work [J].
        source_work_residual_j: f64,
    },
    /// The complete original storage/source balance exceeds its explicit allowance.
    EnergyBalance {
        /// Signed numerical balance defect [J].
        residual_j: f64,
        /// Authored absolute allowance [J].
        tolerance_j: f64,
    },
    /// Cancellation observed before publishing a complete candidate.
    Cancelled,
}
impl From<RotorDriveError> for EccentricRotorDriveError {
    fn from(e: RotorDriveError) -> Self {
        Self::Electrical(e)
    }
}
impl From<RotorError> for EccentricRotorDriveError {
    fn from(e: RotorError) -> Self {
        Self::Electrical(e.into())
    }
}
impl From<ElasticEccentricRotorError> for EccentricRotorDriveError {
    fn from(e: ElasticEccentricRotorError) -> Self {
        Self::Mechanical(e)
    }
}
impl From<DynamicsError> for EccentricRotorDriveError {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanical(e.into())
    }
}
impl core::fmt::Display for EccentricRotorDriveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "eccentric rotor electrical drive refused: {self:?}")
    }
}
impl std::error::Error for EccentricRotorDriveError {}

/// Electrical drive bound to the original eccentric rotor and full body inertia.
///
/// At fixed pose, a shaft impulse `j e_z` changes `omega_z` by
/// `j (I_body^-1)_zz`. The one-dimensional source therefore has effective
/// inertia `1/(I_body^-1)_zz`, which need not equal `I_zz`. Other angular
/// momentum components remain fixed in that source; transverse angular
/// velocity and its kinetic coupling remain in the actual full-body update.
/// This scalar is only a source response, never a substitute trajectory inertia.
#[derive(Clone, Debug)]
pub struct ElasticEccentricCurrentMotor<'a> {
    rotor: &'a ElasticEccentricRotor,
    source: RotorBearingModel,
    motor: RotorCurrentDrive,
    source_inertia_kg_m2: f64,
}

impl ElasticEccentricRotor {
    /// Bind the original reciprocal current/voltage owner to actual body +z.
    ///
    /// `Shaft` must be explicit and the source's scalar coordinate must be zero.
    /// No axisymmetry or balanced-COM approximation is made. Existing conservative
    /// bearing restrictions remain authoritative; damping is not inferred here.
    pub fn current_spin_motor(
        &self,
        motor: RotorCurrentDrive,
        drive_axis: RotorSpinDriveAxis,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ElasticEccentricCurrentMotor<'_>, EccentricRotorDriveError> {
        poll(&mut cancelled)?;
        if drive_axis != RotorSpinDriveAxis::Shaft || motor.shaft_coordinate != 0 {
            return Err(EccentricRotorDriveError::Invalid(
                "shaft-axis actuation and source coordinate zero required",
            ));
        }
        let mobility = self
            .mass_properties()
            .angular_velocity_body_checked(Vec3::new(0.0, 0.0, 1.0))?
            .z;
        let inertia = finite(1.0 / mobility)?;
        if mobility <= 0.0 || inertia <= 0.0 {
            return Err(EccentricRotorDriveError::Invalid(
                "unresolved original shaft mobility",
            ));
        }
        let source = RotorBearingModel::new(1, &[inertia], &[0.0], &[0.0], &[0.0], &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(ElasticEccentricCurrentMotor {
            rotor: self,
            source,
            motor,
            source_inertia_kg_m2: inertia,
        })
    }
}

impl ElasticEccentricCurrentMotor<'_> {
    /// Original tensor's effective inertia for a fixed-pose shaft impulse [kg m²].
    #[must_use]
    pub const fn source_inertia_kg_m2(&self) -> f64 {
        self.source_inertia_kg_m2
    }

    /// Advance native body motion with reciprocal electrical half sources.
    ///
    /// `E(h/2), mechanical(h), E(h/2)` prepares one scalar native midpoint map,
    /// applies each accepted shaft impulse once to the actual body, and retains
    /// both original electrical receipts. The middle step owns all pose motion,
    /// eccentric bearing reaction, full-inertia gyroscopic dynamics and gravity.
    /// Both sources use the actual incoming body-frame `omega_z` as back-EMF.
    ///
    /// The symmetric composition is second order for smooth held-current
    /// operation. Active curtailment, controller switches and continuous voltage
    /// or power peaks require refinement; source-stage endpoint admission is
    /// not a continuous bound on the intervening nonlinear trajectory. There is
    /// no finite winding inductance, battery storage or thermal update here.
    /// The caller's state and source remain unchanged on every refusal.
    pub fn step(
        &self,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        command: EccentricRotorCurrentCommand,
        control: EccentricRotorDriveControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EccentricRotorCurrentStep, EccentricRotorDriveError> {
        poll(&mut cancelled)?;
        let h = control.duration_seconds;
        let half = 0.5 * h;
        if !h.is_finite() || h <= 0.0 || half == 0.0 || half * 2.0 != h {
            return Err(EccentricRotorDriveError::Invalid(
                "representable positive source half interval required",
            ));
        }
        finite(command.external_shaft_torque_nm)?;
        for tolerance in [
            control.angular_impulse_tolerance_nms,
            control.shaft_rate_tolerance_rad_s,
            control.total_energy_tolerance_j,
        ] {
            if !tolerance.is_finite() || tolerance < 0.0 {
                return Err(EccentricRotorDriveError::Invalid(
                    "finite nonnegative native-body allowances required",
                ));
            }
        }
        let operating = self.source.at_speed(0.0, &mut cancelled)?;
        let midpoint =
            operating.midpoint_stepper(half, control.relative_tolerance, &mut cancelled)?;
        let first =
            self.apply_source(&midpoint, state, command.current_a, control, &mut cancelled)?;
        let mechanical = self.rotor.step(
            integrator,
            first.impulse.state_after,
            command.external_shaft_torque_nm,
            h,
            &mut cancelled,
        )?;
        let second = self.apply_source(
            &midpoint,
            mechanical.native.state_after,
            command.current_a,
            control,
            &mut cancelled,
        )?;
        let state_after = second.impulse.state_after;
        let sources = [first, second];
        let totals: [Dd; 5] = core::array::from_fn(|i| {
            sources.iter().fold(Dd::ZERO, |sum, source| {
                let d = &source.electrical;
                sum + Dd::from_f64(
                    [
                        d.shaft_work_j,
                        d.electrical_draw_j,
                        d.electrical_return_j,
                        d.winding_heat_j,
                        d.bridge_heat_j,
                    ][i],
                )
            })
        });
        let mut storage = RigidKineticEnergyChange::new();
        storage.add_body(self.rotor.mass_properties(), state, state_after)?;
        let kinetic = storage.finish_retained()?;
        let residual = kinetic
            + Dd::from_f64(
                mechanical
                    .native
                    .work
                    .gravitational_potential_energy_change_j,
            )
            + Dd::from_f64(mechanical.bearing_potential_change_j)
            + totals[3]
            + totals[4]
            + totals[2]
            - totals[1]
            - Dd::from_f64(mechanical.shaft_work_j);
        let residual_j = finite(residual.to_f64())?;
        if residual_j.abs() > control.total_energy_tolerance_j {
            return Err(EccentricRotorDriveError::EnergyBalance {
                residual_j,
                tolerance_j: control.total_energy_tolerance_j,
            });
        }
        let result = EccentricRotorCurrentStep {
            state_after,
            sources,
            mechanical,
            kinetic_energy_change_j: finite(kinetic.to_f64())?,
            electrical_shaft_work_j: finite(totals[0].to_f64())?,
            electrical_draw_j: finite(totals[1].to_f64())?,
            electrical_return_j: finite(totals[2].to_f64())?,
            winding_heat_j: finite(totals[3].to_f64())?,
            bridge_heat_j: finite(totals[4].to_f64())?,
            total_energy_residual_j: residual_j,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn apply_source(
        &self,
        midpoint: &RotorMidpoint<'_>,
        state: RigidBodyState,
        current_a: f64,
        control: EccentricRotorDriveControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EccentricRotorElectricalSource, EccentricRotorDriveError> {
        poll(cancelled)?;
        let properties = self.rotor.mass_properties();
        let speed = properties
            .angular_velocity_body_checked(state.angular_momentum_body())?
            .z;
        let scalar = RotorTransientState {
            displacement: vec![0.0],
            velocity: vec![speed],
        };
        let electrical = midpoint.step_current_drive(
            &scalar,
            &[0.0],
            current_a,
            self.motor,
            control.electrical,
            &mut *cancelled,
        )?;
        let half = 0.5 * control.duration_seconds;
        let angular_impulse =
            finite((Dd::from_f64(electrical.torque_nm) * Dd::from_f64(half)).to_f64())?;
        let impulse = state.apply_momentum_impulse(
            properties,
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, angular_impulse),
        )?;
        let speed_after = properties
            .angular_velocity_body_checked(impulse.state_after.angular_momentum_body())?
            .z;
        // Source projection tolerances cannot waive the original electrical
        // endpoint domain/capability gates on the actual represented body.
        self.motor.admit_speed(speed_after)?;
        let terminal_voltage_after_v = self.motor.voltage(speed_after, electrical.current_a)?;
        let excess_v = terminal_voltage_after_v.abs() - self.motor.maximum_terminal_voltage_v;
        if excess_v > control.electrical.absolute_voltage_tolerance_v {
            return Err(RotorDriveError::VoltageCompliance {
                initial_voltage_v: electrical.initial_terminal_voltage_v,
                final_voltage_v: terminal_voltage_after_v,
                excess_v,
                tolerance_v: control.electrical.absolute_voltage_tolerance_v,
            }
            .into());
        }
        let source_power_after_w = self
            .motor
            .source_power(finite(terminal_voltage_after_v * electrical.current_a)?)?;
        let source_excess_w = (source_power_after_w - control.electrical.maximum_source_input_w)
            .max(-source_power_after_w - control.electrical.maximum_source_return_w)
            .max(0.0);
        if source_excess_w > control.electrical.absolute_power_tolerance_w {
            return Err(RotorDriveError::PowerBalance {
                source_excess_w,
                native_residual_w: 0.0,
                load_partition_residual_w: 0.0,
                electrical_residual_w: 0.0,
                tolerance_w: control.electrical.absolute_power_tolerance_w,
            }
            .into());
        }
        let speed_residual = finite(speed_after - electrical.motion.state_after.velocity[0])?;
        let work_residual = finite(
            (Dd::from_f64(impulse.work.kinetic_energy_change_j)
                - Dd::from_f64(electrical.shaft_work_j))
            .to_f64(),
        )?;
        let error = impulse.angular_residual_body;
        let angular_residual = error.x.abs().max(error.y.abs()).max(error.z.abs());
        let work_tolerance = finite(half * control.electrical.absolute_power_tolerance_w)?;
        if angular_residual > control.angular_impulse_tolerance_nms
            || speed_residual.abs() > control.shaft_rate_tolerance_rad_s
            || work_residual.abs() > work_tolerance
        {
            return Err(EccentricRotorDriveError::SourceMismatch {
                angular_impulse_residual_nms: angular_residual,
                shaft_rate_residual_rad_s: speed_residual,
                source_work_residual_j: work_residual,
            });
        }
        poll(cancelled)?;
        Ok(EccentricRotorElectricalSource {
            electrical,
            impulse,
            shaft_rate_after_rad_s: speed_after,
            shaft_rate_residual_rad_s: speed_residual,
            terminal_voltage_after_v,
            source_power_after_w,
            source_work_residual_j: work_residual,
        })
    }
}

fn finite(value: f64) -> Result<f64, EccentricRotorDriveError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EccentricRotorDriveError::Invalid(
            "nonfinite electrical/body source response",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EccentricRotorDriveError> {
    if cancelled() {
        Err(EccentricRotorDriveError::Cancelled)
    } else {
        Ok(())
    }
}

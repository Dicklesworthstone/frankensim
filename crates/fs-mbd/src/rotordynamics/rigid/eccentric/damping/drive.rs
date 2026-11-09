//! One admitted electrical and passive-bearing trajectory of the original body.

use super::{
    DampedEccentricRotor, EccentricBearingDampingStep, EccentricRotorDampingControl,
    EccentricRotorDampingError,
};
use crate::rotordynamics::RotorError;
use crate::rotordynamics::rigid::{
    EccentricRotorCurrentCommand, EccentricRotorCurrentStep, EccentricRotorDriveControl,
    EccentricRotorDriveError, ElasticEccentricCurrentMotor, RotorSpinDriveAxis,
};
use crate::rotordynamics::transient::drive::{RotorCurrentDrive, RotorDriveError};
use crate::{DynamicsError, RigidBodyIntegrator, RigidBodyState, RigidKineticEnergyChange};
use fs_math::dd::Dd;

/// Original source/transfer controls and the complete composition's energy gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DampedEccentricDriveControl {
    /// Original current-driven elastic trajectory controls and middle-step gate.
    /// This owner's duration is the sole complete physical interval.
    pub drive: EccentricRotorDriveControl,
    /// Original passive-source controls. Its total-energy allowance gates the
    /// COMPLETE electrical/damping trajectory, including every native loss.
    pub damping: EccentricRotorDampingControl,
}

/// Original motor law at an actual endpoint, without an extra source debit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DampedEccentricMotorEndpoint {
    /// Current actually accepted by the adjacent first/last electrical source [A].
    pub current_a: f64,
    /// Actual full-body body-z angular velocity [rad/s].
    pub shaft_rate_rad_s: f64,
    /// Original back-EMF plus winding resistive drop at that physical endpoint [V].
    pub terminal_voltage_v: f64,
    /// Original signed source power: positive draw, negative return [W].
    pub source_power_w: f64,
}

/// One physical endpoint with all original source receipts and an outer ledger.
#[derive(Clone, Debug, PartialEq)]
pub struct DampedEccentricCurrentStep {
    /// Sole physical endpoint, after the final passive source.
    pub state_after: RigidBodyState,
    /// Original passive half sources at the initial/final physical poses.
    pub damping_sources: [EccentricBearingDampingStep; 2],
    /// Original electrical/elastic trajectory and its winding/bridge accounts.
    /// Its endpoint precedes the final passive source; use state_after above.
    /// Electrical shaft work is internal exchange, not an extra external debit.
    pub current: EccentricRotorCurrentStep,
    /// Original-law admission at the complete initial and final native states.
    /// These are diagnostics, not extra source intervals or energy transfers.
    pub electrical_endpoints: [DampedEccentricMotorEndpoint; 2],
    /// Complete initial-to-final original full-tensor kinetic change [J].
    pub kinetic_energy_change_j: f64,
    /// Independent physical heat from both passive bearing sources [J].
    pub bearing_heat_j: f64,
    /// delta(K+V_gravity+V_bearing) + bearing/winding/bridge heat
    /// + electrical return - draw - actual external shaft work [J].
    /// Signed numerical defect; never assigned to heat or efficiency.
    pub total_energy_residual_j: f64,
}

/// A refused coupled candidate publishes neither body motion nor source debit.
#[derive(Clone, Debug, PartialEq)]
pub enum DampedEccentricDriveError {
    /// Invalid numerical interval/budget or unrepresentable response.
    Invalid(&'static str),
    /// Original passive-source admission or transfer refused.
    Damping(EccentricRotorDampingError),
    /// Original electrical/elastic owner or actual endpoint admission refused.
    Drive(EccentricRotorDriveError),
    /// Complete endpoint storage/source balance exceeds its authored allowance.
    EnergyBalance {
        /// Signed complete balance defect [J].
        residual_j: f64,
        /// Original damping control's whole-composition energy allowance [J].
        tolerance_j: f64,
    },
    /// Cancellation observed before publication.
    Cancelled,
}
impl From<EccentricRotorDampingError> for DampedEccentricDriveError {
    fn from(e: EccentricRotorDampingError) -> Self {
        Self::Damping(e)
    }
}
impl From<EccentricRotorDriveError> for DampedEccentricDriveError {
    fn from(e: EccentricRotorDriveError) -> Self {
        Self::Drive(e)
    }
}
impl From<RotorDriveError> for DampedEccentricDriveError {
    fn from(e: RotorDriveError) -> Self {
        Self::Drive(e.into())
    }
}
impl From<RotorError> for DampedEccentricDriveError {
    fn from(e: RotorError) -> Self {
        Self::Drive(e.into())
    }
}
impl From<DynamicsError> for DampedEccentricDriveError {
    fn from(e: DynamicsError) -> Self {
        Self::Drive(e.into())
    }
}
impl core::fmt::Display for DampedEccentricDriveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "damped eccentric rotor electrical drive refused: {self:?}"
        )
    }
}
impl std::error::Error for DampedEccentricDriveError {}

/// Current-regulated shaft actuation on the SAME original passive eccentric body.
///
/// Both prepared owners borrow one elastic geometry/inertia. Electrical and
/// passive sources retain their original laws and numerical gates; no additional
/// rotor, reduced trajectory inertia, integrator or iterative controller is added.
#[derive(Clone, Debug)]
pub struct DampedEccentricCurrentMotor<'a> {
    rotor: &'a DampedEccentricRotor,
    current: ElasticEccentricCurrentMotor<'a>,
    motor: RotorCurrentDrive,
}
impl DampedEccentricRotor {
    /// Bind the original current/voltage owner to explicit actual body-z actuation.
    ///
    /// The original electrical owner admits axis and scalar-coordinate support.
    /// Static parameter, actual-state and source limits are admitted at step time.
    pub fn current_spin_motor(
        &self,
        motor: RotorCurrentDrive,
        drive_axis: RotorSpinDriveAxis,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DampedEccentricCurrentMotor<'_>, DampedEccentricDriveError> {
        poll(&mut cancelled)?;
        let current =
            self.elastic_component()
                .current_spin_motor(motor, drive_axis, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(DampedEccentricCurrentMotor {
            rotor: self,
            current,
            motor,
        })
    }
}
impl DampedEccentricCurrentMotor<'_> {
    /// Original full tensor's effective inertia for the electrical source only [kg m²].
    #[must_use]
    pub const fn source_inertia_kg_m2(&self) -> f64 {
        self.current.source_inertia_kg_m2()
    }

    /// Advance D(h/2), original current-driven elastic trajectory(h), D(h/2).
    ///
    /// Only the original middle elastic trajectory advances pose. Its two motor
    /// sources and both outer bearing sources use the actual incoming momenta.
    /// Physical bearing heat and original winding/bridge heat are kept separate;
    /// original electrical shaft work is exchanged internally with native storage.
    ///
    /// Actual initial/final speed, voltage and signed source power are admitted
    /// using the adjacent actually accepted current. Damping can increase shaft
    /// speed through reciprocal translation/rotation coupling even while losing
    /// total kinetic energy, so the middle drive's endpoint is not authoritative.
    /// An incompatible outer endpoint refuses atomically; current is not silently
    /// reselected after work has been integrated. No new source energy is charged
    /// for these endpoint diagnostics.
    ///
    /// The composition is second order for smooth fixed coefficients and held
    /// uncurtailed current with resolved mechanical/source timescales. Active
    /// current curtailment, controller transitions and nonlinear interstage power
    /// peaks need refinement. Endpoint admission is not a continuous-time peak
    /// guarantee. Inductance, battery storage and thermal feedback are not inferred.
    pub fn step(
        &self,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        command: EccentricRotorCurrentCommand,
        control: DampedEccentricDriveControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DampedEccentricCurrentStep, DampedEccentricDriveError> {
        poll(&mut cancelled)?;
        super::validate(control.damping)?;
        self.motor
            .admit(1, control.drive.electrical, command.current_a)?;
        self.motor.admit_speed(self.speed(state)?)?;
        let h = control.drive.duration_seconds;
        let half = 0.5 * h;
        if !half.is_finite() || half <= 0.0 || 2.0 * half != h {
            return Err(DampedEccentricDriveError::Invalid(
                "representable positive source half interval required",
            ));
        }
        let first = self
            .rotor
            .relax_damping(state, half, control.damping, &mut cancelled)?;
        let current = self.current.step(
            integrator,
            first.impulse.state_after,
            command,
            control.drive,
            &mut cancelled,
        )?;
        let last =
            self.rotor
                .relax_damping(current.state_after, half, control.damping, &mut cancelled)?;
        let state_after = last.impulse.state_after;
        let endpoints = [
            self.endpoint(state, current.sources[0].electrical.current_a)?,
            self.endpoint(state_after, current.sources[1].electrical.current_a)?,
        ];
        let excess_v = endpoints[0]
            .terminal_voltage_v
            .abs()
            .max(endpoints[1].terminal_voltage_v.abs())
            - self.motor.maximum_terminal_voltage_v;
        let electrical = control.drive.electrical;
        if excess_v > electrical.absolute_voltage_tolerance_v {
            return Err(RotorDriveError::VoltageCompliance {
                initial_voltage_v: endpoints[0].terminal_voltage_v,
                final_voltage_v: endpoints[1].terminal_voltage_v,
                excess_v,
                tolerance_v: electrical.absolute_voltage_tolerance_v,
            }
            .into());
        }
        let source_excess_w = endpoints.iter().fold(0.0_f64, |excess, endpoint| {
            excess
                .max(endpoint.source_power_w - electrical.maximum_source_input_w)
                .max(-endpoint.source_power_w - electrical.maximum_source_return_w)
        });
        if source_excess_w > electrical.absolute_power_tolerance_w {
            return Err(RotorDriveError::PowerBalance {
                source_excess_w,
                native_residual_w: 0.0,
                load_partition_residual_w: 0.0,
                electrical_residual_w: 0.0,
                tolerance_w: electrical.absolute_power_tolerance_w,
            }
            .into());
        }
        let mut storage = RigidKineticEnergyChange::new();
        storage.add_body(
            self.rotor.elastic_component().mass_properties(),
            state,
            state_after,
        )?;
        let kinetic = storage.finish_retained()?;
        let bearing_heat =
            Dd::from_f64(first.dissipated_energy_j) + Dd::from_f64(last.dissipated_energy_j);
        // Reduce original half-source accounts directly. Motor shaft work is
        // already internal to full-body storage and must not be counted twice.
        let electrical_exchange = current.sources.iter().fold(Dd::ZERO, |sum, source| {
            let s = &source.electrical;
            sum + Dd::from_f64(s.winding_heat_j)
                + Dd::from_f64(s.bridge_heat_j)
                + Dd::from_f64(s.electrical_return_j)
                - Dd::from_f64(s.electrical_draw_j)
        });
        let residual = finite(
            (kinetic
                + Dd::from_f64(
                    current
                        .mechanical
                        .native
                        .work
                        .gravitational_potential_energy_change_j,
                )
                + Dd::from_f64(current.mechanical.bearing_potential_change_j)
                + bearing_heat
                + electrical_exchange
                - Dd::from_f64(current.mechanical.shaft_work_j))
            .to_f64(),
        )?;
        if residual.abs() > control.damping.total_energy_tolerance_j {
            return Err(DampedEccentricDriveError::EnergyBalance {
                residual_j: residual,
                tolerance_j: control.damping.total_energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(DampedEccentricCurrentStep {
            state_after,
            damping_sources: [first, last],
            current,
            electrical_endpoints: endpoints,
            kinetic_energy_change_j: finite(kinetic.to_f64())?,
            bearing_heat_j: finite(bearing_heat.to_f64())?,
            total_energy_residual_j: residual,
        })
    }

    fn speed(&self, state: RigidBodyState) -> Result<f64, DampedEccentricDriveError> {
        Ok(self
            .rotor
            .elastic_component()
            .mass_properties()
            .angular_velocity_body_checked(state.angular_momentum_body())?
            .z)
    }
    fn endpoint(
        &self,
        state: RigidBodyState,
        current_a: f64,
    ) -> Result<DampedEccentricMotorEndpoint, DampedEccentricDriveError> {
        let speed = self.speed(state)?;
        self.motor.admit_speed(speed)?;
        let voltage = self.motor.voltage(speed, current_a)?;
        let power = self.motor.source_power(finite(
            (Dd::from_f64(voltage) * Dd::from_f64(current_a)).to_f64(),
        )?)?;
        Ok(DampedEccentricMotorEndpoint {
            current_a,
            shaft_rate_rad_s: speed,
            terminal_voltage_v: voltage,
            source_power_w: power,
        })
    }
}
fn finite(value: f64) -> Result<f64, DampedEccentricDriveError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DampedEccentricDriveError::Invalid(
            "nonfinite coupled electrical/damping response",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DampedEccentricDriveError> {
    if cancelled() {
        Err(DampedEccentricDriveError::Cancelled)
    } else {
        Ok(())
    }
}

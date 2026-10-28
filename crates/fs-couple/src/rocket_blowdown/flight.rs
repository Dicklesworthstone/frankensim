//! Source atmosphere and nonlinear aerodynamic loading on a finite-gas rocket.
//!
//! A tentative half step predicts the environment midpoint. The accepted full
//! step starts from the ORIGINAL gas/rigid state under those midpoint loads;
//! the trial never consumes accepted propellant. The native gas midpoint solve
//! and native attitude substeps retain their original controls and reports.
//! Endpoint atmosphere and the caller coefficient law are admitted before either
//! physical owner is published. The smooth composition is second order and needs
//! temporal refinement; it does not locate impact, cutoff or source-domain events.

use super::environment::spherical::SphericalStandardAtmosphere;
use super::environment::{AerodynamicLoad, FlightEnvironmentError, RelativeAirflow};
use super::rigid::{
    GasRigidRocketCoastEnvironment, GasRigidRocketEnvironment, GasRigidRocketError,
    GasRigidRocketSolve, GasRigidRocketState, GasRigidRocketStep,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::RocketNozzle;
use fs_mbd::{
    Vec3,
    rocket::rigid::{RigidRocketError, RigidRocketStep},
};

/// Source sampling and aerodynamic reference geometry, in the body frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketFlightGeometry {
    /// Aerodynamic resultant application point relative to the common COM [m].
    pub aerodynamic_application_arm_body_m: Vec3,
    /// Explicit characteristic length used by the supplied coefficient law [m].
    pub aerodynamic_reference_length_m: f64,
    /// Explicit area used by that coefficient law [m²].
    pub aerodynamic_reference_area_m2: f64,
    /// Point where ambient nozzle pressure is sampled, relative to the COM [m].
    pub nozzle_ambient_arm_body_m: Vec3,
}

/// Explicit propulsion control mode; failed engine/source solves never select coast.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RocketFlightPropulsion {
    /// Closed valve: the identical insulated gas state remains aboard the vehicle.
    ClosedValve,
    /// The actual source-bearing reservoir/nozzle supplies accepted mass and impulse.
    Powered {
        /// Original throat/exit geometry.
        nozzle: RocketNozzle,
        /// Unit body-frame vector along the total nozzle thrust.
        thrust_direction_body: Vec3,
        /// Effective force application point relative to the common COM [m].
        thrust_application_arm_body_m: Vec3,
        /// Explicit diagonal outflow co-rotation moment per mass [m²].
        co_rotating_exhaust_specific_inertia_m2: Vec3,
        /// Additional exported specific body spin [m²/s], excluding axial jet torque.
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3,
    },
}

/// Propulsion command and additional loads held fixed in their declared frames.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketFlightForcing {
    /// Explicit powered or closed-valve mode.
    pub propulsion: RocketFlightPropulsion,
    /// World force [N] excluding atmosphere, gravity and engine thrust.
    pub additional_force_world_n: Vec3,
    /// Body torque [N m] excluding the computed aerodynamic and engine moments.
    pub additional_torque_body_n_m: Vec3,
}

/// Actual local source and aerodynamic queries at one complete candidate state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketFlightEnvironmentSample {
    /// Source geometric altitude of the vehicle COM [m].
    pub center_altitude_m: f64,
    /// Source geometric altitude at the drag point [m].
    pub aerodynamic_altitude_m: f64,
    /// Source geometric altitude at the nozzle ambient point [m].
    pub nozzle_altitude_m: f64,
    /// Ambient pressure at that actual nozzle point [Pa].
    pub ambient_pressure_pa: f64,
    /// Actual relative point flow, including rotation and source-bound wind.
    pub airflow: RelativeAirflow,
    /// Aerodynamic force, matching body torque and relative mechanical power.
    pub aerodynamic_load: AerodynamicLoad,
    /// The actual caller coefficient evaluated at this flow state.
    pub drag_coefficient: f64,
    /// Source gravity at the actual vehicle COM [m/s²].
    pub gravity_world_m_s2: Vec3,
    /// Source specific gravitational potential relative to infinity [J/kg].
    pub specific_gravitational_potential_j_kg: f64,
}

/// Original accepted powered or closed-valve physics report.
#[derive(Clone, Debug, PartialEq)]
pub enum RocketFlightStep {
    /// Full accepted finite-gas blowdown and rigid motion report.
    Powered(GasRigidRocketStep),
    /// Full accepted native coast; the reservoir is unchanged.
    Coast(RigidRocketStep),
}
impl RocketFlightStep {
    /// Original native rigid motion, scalar/world impulse and spin accounting.
    #[must_use]
    pub const fn mechanics(&self) -> &RigidRocketStep {
        match self {
            Self::Powered(step) => &step.mechanics,
            Self::Coast(step) => step,
        }
    }
}

/// Failure before publication of the joint atmosphere-admitted gas/rigid endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum AtmosphericFlightError {
    /// Invalid geometry, fixed load or interval.
    Invalid(&'static str),
    /// Original source/airflow/coefficient owner refused.
    Environment(FlightEnvironmentError),
    /// Original gas/rigid coupling refused.
    Coupling(GasRigidRocketError),
    /// Original native mass or body-point kinematics refused.
    Mechanics(RigidRocketError),
    /// Cancellation at the atmosphere composition boundary.
    Cancelled,
}
impl core::fmt::Display for AtmosphericFlightError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "atmospheric finite-gas rocket: {self:?}")
    }
}
impl std::error::Error for AtmosphericFlightError {}
impl From<FlightEnvironmentError> for AtmosphericFlightError {
    fn from(value: FlightEnvironmentError) -> Self {
        Self::Environment(value)
    }
}
impl From<GasRigidRocketError> for AtmosphericFlightError {
    fn from(value: GasRigidRocketError) -> Self {
        Self::Coupling(value)
    }
}
impl From<RigidRocketError> for AtmosphericFlightError {
    fn from(value: RigidRocketError) -> Self {
        Self::Mechanics(value)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AtmosphericFlightError> {
    if cancelled() {
        Err(AtmosphericFlightError::Cancelled)
    } else {
        Ok(())
    }
}

/// Accepted endpoint, native report and the actual midpoint/endpoint environment.
#[derive(Clone, Debug, PartialEq)]
pub struct AtmosphericFlightCandidate {
    endpoint: GasRigidRocketState,
    step: RocketFlightStep,
    midpoint_environment: RocketFlightEnvironmentSample,
    endpoint_environment: RocketFlightEnvironmentSample,
}
impl AtmosphericFlightCandidate {
    /// Joint endpoint after final source/coefficient-domain admission.
    #[must_use]
    pub const fn endpoint(&self) -> &GasRigidRocketState {
        &self.endpoint
    }
    /// Original accepted gas/rigid owner report, excluding the discarded predictor.
    #[must_use]
    pub const fn step(&self) -> &RocketFlightStep {
        &self.step
    }
    /// Environment used for the accepted full interval.
    #[must_use]
    pub const fn midpoint_environment(&self) -> &RocketFlightEnvironmentSample {
        &self.midpoint_environment
    }
    /// Independently evaluated/admitted endpoint environment; not used to resample thrust.
    #[must_use]
    pub const fn endpoint_environment(&self) -> &RocketFlightEnvironmentSample {
        &self.endpoint_environment
    }
    /// Accept the complete joint state and original reports without another solve.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasRigidRocketState,
        RocketFlightStep,
        RocketFlightEnvironmentSample,
        RocketFlightEnvironmentSample,
    ) {
        (
            self.endpoint,
            self.step,
            self.midpoint_environment,
            self.endpoint_environment,
        )
    }
}

fn sample_environment(
    state: &GasRigidRocketState,
    atmosphere: &SphericalStandardAtmosphere,
    geometry: RocketFlightGeometry,
    drag_coefficient: &mut impl FnMut(&RelativeAirflow) -> Result<f64, FlightEnvironmentError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RocketFlightEnvironmentSample, AtmosphericFlightError> {
    poll(cancelled)?;
    let frame = state.inertial_frame();
    let vehicle = state.vehicle();
    let native_body = vehicle.rigid_body()?;
    let pose = native_body.pose();
    let center = atmosphere.sample_at(frame, pose.position_world())?;
    let aerodynamic_point = pose
        .point_world_from_body(geometry.aerodynamic_application_arm_body_m)
        .map_err(RigidRocketError::from)?;
    let nozzle_point = pose
        .point_world_from_body(geometry.nozzle_ambient_arm_body_m)
        .map_err(RigidRocketError::from)?;
    let aerodynamic = atmosphere.sample_at(frame, aerodynamic_point)?;
    let nozzle = atmosphere.sample_at(frame, nozzle_point)?;
    let properties = vehicle.mass_properties()?;
    let omega = properties
        .angular_velocity_body_checked(vehicle.angular_momentum_body_kg_m2_s())
        .map_err(RigidRocketError::from)?;
    let airflow = aerodynamic.air.at_body_point(
        frame,
        vehicle.translation().velocity_m_s,
        vehicle.orientation(),
        omega,
        geometry.aerodynamic_application_arm_body_m,
        geometry.aerodynamic_reference_length_m,
    )?;
    let coefficient = drag_coefficient(&airflow)?;
    let load = airflow.drag(geometry.aerodynamic_reference_area_m2, coefficient)?;
    let (gravity, potential) = atmosphere.gravity_at(frame, pose.position_world())?;
    poll(cancelled)?;
    Ok(RocketFlightEnvironmentSample {
        center_altitude_m: center.geometric_altitude_m,
        aerodynamic_altitude_m: aerodynamic.geometric_altitude_m,
        nozzle_altitude_m: nozzle.geometric_altitude_m,
        ambient_pressure_pa: nozzle.air.pressure_pa(),
        airflow,
        aerodynamic_load: load,
        drag_coefficient: coefficient,
        gravity_world_m_s2: gravity,
        specific_gravitational_potential_j_kg: potential,
    })
}

fn advance_under_sample(
    state: &GasRigidRocketState,
    sample: RocketFlightEnvironmentSample,
    forcing: RocketFlightForcing,
    duration_s: f64,
    controls: GasRigidRocketSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(GasRigidRocketState, RocketFlightStep), AtmosphericFlightError> {
    let force = sample
        .aerodynamic_load
        .force_world_n
        .add(forcing.additional_force_world_n);
    let torque = sample
        .aerodynamic_load
        .torque_body_n_m
        .add(forcing.additional_torque_body_n_m);
    match forcing.propulsion {
        RocketFlightPropulsion::ClosedValve => {
            let candidate = state.coast(
                GasRigidRocketCoastEnvironment {
                    inertial_frame: state.inertial_frame(),
                    gravity_world_m_s2: sample.gravity_world_m_s2,
                    external_force_world_n: force,
                    external_torque_body_n_m: torque,
                },
                duration_s,
                controls.mechanics,
                cancelled,
            )?;
            let (endpoint, mechanics) = candidate.into_parts();
            Ok((endpoint, RocketFlightStep::Coast(mechanics)))
        }
        RocketFlightPropulsion::Powered {
            nozzle,
            thrust_direction_body,
            thrust_application_arm_body_m,
            co_rotating_exhaust_specific_inertia_m2,
            exhaust_swirl_specific_angular_momentum_body_m2_s,
        } => {
            let candidate = state.advance(
                GasRigidRocketEnvironment {
                    inertial_frame: state.inertial_frame(),
                    ambient_pressure_pa: sample.ambient_pressure_pa,
                    thrust_direction_body,
                    thrust_application_arm_body_m,
                    co_rotating_exhaust_specific_inertia_m2,
                    exhaust_swirl_specific_angular_momentum_body_m2_s,
                    gravity_world_m_s2: sample.gravity_world_m_s2,
                    external_force_world_n: force,
                    external_torque_body_n_m: torque,
                },
                nozzle,
                duration_s,
                controls,
                cancelled,
            )?;
            let (endpoint, report) = candidate.into_parts();
            Ok((endpoint, RocketFlightStep::Powered(report)))
        }
    }
}

/// Advance one bounded interval with midpoint atmosphere, drag and central gravity.
///
/// Atmosphere is sampled at actual COM, drag and nozzle points. The supplied
/// coefficient callback must be pure/repeatable: it may be called on a discarded
/// predictor and final admission state, and its external side effects cannot be
/// rolled back. Mach, Reynolds and body-relative flow are the actual local values;
/// no empirical coefficient, lift law or extrapolation is invented here.
///
/// The unaccepted half-step predictor is discarded. The full gas/rigid candidate
/// restarts from the original state, and final atmosphere/coefficient admission
/// precedes publication. Refinement is required for smoothly changing fields;
/// source boundaries, ground contact, automatic cutoff and staging events are
/// not localized by this operation.
pub fn advance_atmospheric_rocket(
    state: &GasRigidRocketState,
    atmosphere: &SphericalStandardAtmosphere,
    geometry: RocketFlightGeometry,
    forcing: RocketFlightForcing,
    duration_s: f64,
    controls: GasRigidRocketSolve,
    mut drag_coefficient: impl FnMut(&RelativeAirflow) -> Result<f64, FlightEnvironmentError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AtmosphericFlightCandidate, AtmosphericFlightError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || 0.5 * duration_s == 0.0
        || !geometry.aerodynamic_application_arm_body_m.is_finite()
        || !geometry.nozzle_ambient_arm_body_m.is_finite()
        || !geometry.aerodynamic_reference_length_m.is_finite()
        || geometry.aerodynamic_reference_length_m <= 0.0
        || !geometry.aerodynamic_reference_area_m2.is_finite()
        || geometry.aerodynamic_reference_area_m2 <= 0.0
        || !forcing.additional_force_world_n.is_finite()
        || !forcing.additional_torque_body_n_m.is_finite()
    {
        return Err(AtmosphericFlightError::Invalid(
            "finite body geometry and loads, positive area/length and representable positive interval required",
        ));
    }
    let initial = sample_environment(
        state,
        atmosphere,
        geometry,
        &mut drag_coefficient,
        &mut cancelled,
    )?;
    let (predictor, _) = advance_under_sample(
        state,
        initial,
        forcing,
        0.5 * duration_s,
        controls,
        &mut cancelled,
    )?;
    let midpoint = sample_environment(
        &predictor,
        atmosphere,
        geometry,
        &mut drag_coefficient,
        &mut cancelled,
    )?;
    let (endpoint, step) = advance_under_sample(
        state,
        midpoint,
        forcing,
        duration_s,
        controls,
        &mut cancelled,
    )?;
    let final_sample = sample_environment(
        &endpoint,
        atmosphere,
        geometry,
        &mut drag_coefficient,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(AtmosphericFlightCandidate {
        endpoint,
        step,
        midpoint_environment: midpoint,
        endpoint_environment: final_sample,
    })
}

#[cfg(test)]
#[path = "flight/tests.rs"]
mod tests;

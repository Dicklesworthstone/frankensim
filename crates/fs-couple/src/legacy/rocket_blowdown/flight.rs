//! Source atmosphere and nonlinear aerodynamic loading on a finite-gas rocket.
//!
//! A tentative half step predicts the environment midpoint. The accepted full
//! step starts from the ORIGINAL gas/rigid state under those midpoint loads;
//! the trial never consumes accepted propellant. The native gas midpoint solve
//! and native attitude substeps retain their original controls and reports.
//! Endpoint atmosphere and the caller coefficient law are admitted before either
//! physical owner is published. The smooth composition is second order and needs
//! temporal refinement; it does not locate impact, cutoff or source-domain events.
//! The outer interval must also resolve local quadratic drag at the initial and
//! predicted midpoint states. Native attitude substeps do not resample drag.

use super::environment::spherical::SphericalStandardAtmosphere;
use super::environment::{AerodynamicLoad, FlightEnvironmentError, RelativeAirflow};
use super::rigid::{
    GasRigidRocketCoastEnvironment, GasRigidRocketEnvironment, GasRigidRocketError,
    GasRigidRocketSolve, GasRigidRocketState, GasRigidRocketStep,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::RocketNozzle;
use fs_mbd::{
    Vec3, Wrench,
    rocket::rigid::{RigidRocketError, RigidRocketInertia, RigidRocketStep, RocketInertiaModel},
    spatial_loads::{PreparedSpatialLoads, SpatialLoadError, WorldPointForce},
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
pub enum RocketFlightPropulsion<E = Vec3> {
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
        /// Explicit native outflow co-rotation inertia per mass [m²].
        /// The type matches the vehicle's diagonal or full-tensor inertia model.
        co_rotating_exhaust_specific_inertia_m2: E,
        /// Additional exported specific body spin [m²/s], excluding axial jet torque.
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3,
    },
}

/// Propulsion command and additional loads held fixed in their declared frames.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketFlightForcing<E = Vec3> {
    /// Explicit powered or closed-valve mode.
    pub propulsion: RocketFlightPropulsion<E>,
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
pub enum RocketFlightStep<I = RigidRocketInertia> {
    /// Full accepted finite-gas blowdown and rigid motion report.
    Powered(GasRigidRocketStep<I>),
    /// Full accepted native coast; the reservoir is unchanged.
    Coast(RigidRocketStep<I>),
}
impl<I> RocketFlightStep<I> {
    /// Original native rigid motion, scalar/world impulse and spin accounting.
    #[must_use]
    pub const fn mechanics(&self) -> &RigidRocketStep<I> {
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
    /// The original prepared spatial-load owner refused force/arm arithmetic.
    SpatialLoads(SpatialLoadError),
    /// The outer explicit interval does not resolve the sampled quadratic drag.
    /// Reduce the outer interval; smaller frozen-load mechanics substeps do not
    /// resolve this atmosphere/drag update.
    DragResolution {
        /// Sample at which the local drag timescale was exceeded.
        sample: &'static str,
        /// Requested complete outer interval [s].
        duration_s: f64,
        /// Maximum interval from this sample's drag and point mobility [s].
        maximum_duration_s: f64,
    },
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
impl From<SpatialLoadError> for AtmosphericFlightError {
    fn from(value: SpatialLoadError) -> Self {
        match value {
            SpatialLoadError::Cancelled => Self::Cancelled,
            other => Self::SpatialLoads(other),
        }
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
pub struct AtmosphericFlightCandidate<I = RigidRocketInertia> {
    endpoint: GasRigidRocketState<I>,
    step: RocketFlightStep<I>,
    midpoint_environment: RocketFlightEnvironmentSample,
    endpoint_environment: RocketFlightEnvironmentSample,
}
impl<I> AtmosphericFlightCandidate<I> {
    /// Joint endpoint after final source/coefficient-domain admission.
    #[must_use]
    pub const fn endpoint(&self) -> &GasRigidRocketState<I> {
        &self.endpoint
    }
    /// Original accepted gas/rigid owner report, excluding the discarded predictor.
    #[must_use]
    pub const fn step(&self) -> &RocketFlightStep<I> {
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
        GasRigidRocketState<I>,
        RocketFlightStep<I>,
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

fn sample_environment<I: RocketInertiaModel>(
    state: &GasRigidRocketState<I>,
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

// For frozen rho and Cd, F(u)=-k|u|u has derivative norm 2|F|/|u|.
// A force at body arm r changes point velocity through
// G = I/m - [r]x I_body^-1 [r]x. The rotational term is positive
// semidefinite, so its trace bounds its largest eigenvalue. Evaluating that
// trace in the ORIGINAL principal frame keeps this bound basis independent.
// h*||dF/du||*||G|| <= 1 is a conservative local explicit-resolution rule,
// not an error certificate for arbitrary Cd callbacks or evolving wind.
fn admit_drag_interval<I: RocketInertiaModel>(
    state: &GasRigidRocketState<I>,
    sample: RocketFlightEnvironmentSample,
    geometry: RocketFlightGeometry,
    duration_s: f64,
    sample_name: &'static str,
) -> Result<(), AtmosphericFlightError> {
    let speed = sample.airflow.speed_m_s();
    let magnitude = sample.airflow.dynamic_pressure_pa()
        * geometry.aerodynamic_reference_area_m2
        * sample.drag_coefficient;
    if speed == 0.0 || magnitude == 0.0 {
        return Ok(());
    }
    let properties = state.vehicle().mass_properties()?;
    let inertia = properties.principal_inertia_body();
    let arm = properties
        .principal_axes_to_body()
        .rotate_world_to_body(geometry.aerodynamic_application_arm_body_m);
    let rotational_mobility = (arm.y * arm.y + arm.z * arm.z) / inertia.x
        + (arm.x * arm.x + arm.z * arm.z) / inertia.y
        + (arm.x * arm.x + arm.y * arm.y) / inertia.z;
    let mobility_bound = 1.0 / properties.mass() + rotational_mobility;
    let rate_bound = 2.0 * (magnitude / speed) * mobility_bound;
    if !mobility_bound.is_finite() || !rate_bound.is_finite() {
        return Err(AtmosphericFlightError::Invalid(
            "quadratic drag relaxation rate is unrepresentable",
        ));
    }
    if rate_bound > 0.0 {
        let maximum_duration_s = 1.0 / rate_bound;
        if duration_s > maximum_duration_s {
            return Err(AtmosphericFlightError::DragResolution {
                sample: sample_name,
                duration_s,
                maximum_duration_s,
            });
        }
    }
    Ok(())
}

fn advance_under_sample<I: RocketInertiaModel>(
    state: &GasRigidRocketState<I>,
    sample: RocketFlightEnvironmentSample,
    geometry: RocketFlightGeometry,
    forcing: RocketFlightForcing<I::ExhaustInertia>,
    duration_s: f64,
    controls: GasRigidRocketSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(GasRigidRocketState<I>, RocketFlightStep<I>), AtmosphericFlightError> {
    // The sampled force is held in WORLD axes during this outer stage. Its
    // moment therefore follows the actual rotating attachment at each native
    // kick, rather than becoming a frozen body-frame follower torque.
    let loads = PreparedSpatialLoads::new(
        Wrench::ZERO,
        Vec3::ZERO,
        &[WorldPointForce {
            arm_body: geometry.aerodynamic_application_arm_body_m,
            force_world: sample.aerodynamic_load.force_world_n,
        }],
        &mut *cancelled,
    )?;
    match forcing.propulsion {
        RocketFlightPropulsion::ClosedValve => {
            let candidate = state.coast_spatial(
                GasRigidRocketCoastEnvironment {
                    inertial_frame: state.inertial_frame(),
                    gravity_world_m_s2: sample.gravity_world_m_s2,
                    external_force_world_n: forcing.additional_force_world_n,
                    external_torque_body_n_m: forcing.additional_torque_body_n_m,
                },
                &loads,
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
            let candidate = state.advance_spatial(
                GasRigidRocketEnvironment {
                    inertial_frame: state.inertial_frame(),
                    ambient_pressure_pa: sample.ambient_pressure_pa,
                    thrust_direction_body,
                    thrust_application_arm_body_m,
                    co_rotating_exhaust_specific_inertia_m2,
                    exhaust_swirl_specific_angular_momentum_body_m2_s,
                    gravity_world_m_s2: sample.gravity_world_m_s2,
                    external_force_world_n: forcing.additional_force_world_n,
                    external_torque_body_n_m: forcing.additional_torque_body_n_m,
                },
                nozzle,
                &loads,
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
/// Each stage holds the sampled aerodynamic force in world axes while its
/// material attachment rotates under the original native dynamics. Every native
/// torque kick reevaluates the arm's moment; the sampled torque remains an
/// observation at that sampled pose, not a body-fixed torque for the interval.
///
/// The unaccepted half-step predictor is discarded. The full gas/rigid candidate
/// restarts from the original state, and final atmosphere/coefficient admission
/// precedes publication. Refinement is required for smoothly changing fields;
/// source boundaries, ground contact, automatic cutoff and staging events are
/// not localized by this operation.
///
/// The outer interval must resolve the initial AND predicted midpoint drag:
/// `dt * 2|F|/|u| * (1/m + tr(-[r]x I_body^-1 [r]x)) <= 1`. This accounts
/// for translational and rotational response at the actual application point.
/// The rate is for locally frozen density and drag coefficient; it does not
/// bound the caller coefficient's derivatives, changing wind, or other loads.
/// [`AtmosphericFlightError::DragResolution`] reports an actionable outer-step
/// limit without publishing either owner. `controls.mechanics` only resolves
/// motion under held loads and cannot satisfy this condition by itself.
pub fn advance_atmospheric_rocket<I: RocketInertiaModel>(
    state: &GasRigidRocketState<I>,
    atmosphere: &SphericalStandardAtmosphere,
    geometry: RocketFlightGeometry,
    forcing: RocketFlightForcing<I::ExhaustInertia>,
    duration_s: f64,
    controls: GasRigidRocketSolve,
    mut drag_coefficient: impl FnMut(&RelativeAirflow) -> Result<f64, FlightEnvironmentError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AtmosphericFlightCandidate<I>, AtmosphericFlightError> {
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
    admit_drag_interval(state, initial, geometry, duration_s, "initial")?;
    let (predictor, _) = advance_under_sample(
        state,
        initial,
        geometry,
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
    admit_drag_interval(&predictor, midpoint, geometry, duration_s, "midpoint")?;
    let (endpoint, step) = advance_under_sample(
        state,
        midpoint,
        geometry,
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

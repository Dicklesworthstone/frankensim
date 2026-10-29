//! Automatic sustained sphere/plane contact for balanced articulated carriers.
//!
//! The original force/contact/drift owner advances every base and joint. The
//! co-located carrier model makes its unforced numerical COM drift affine in
//! real arithmetic even while its shell and internal rotors rotate. Every plane
//! receives that complete affine-path check and an actual rounded endpoint
//! check; material witnesses are regenerated at each accepted step.

use super::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectorySystem, CarrierPlane, CarrierPlaneConfig,
    CarrierPlaneError, Work, admit, finite_vector, observe, plane_normal, plane_position,
};
use crate::articulated::motion::FreeFloatingMotionState;
use crate::articulated_contacts::step::{
    ArticulatedContactPath, ArticulatedContactStepConfig, ArticulatedContactStepError,
    step_articulated_frictional_contacts,
};
use crate::articulated_contacts::trajectory::ArticulatedTrajectoryMomentum;
use fs_ga::Vec3;

/// Bounds and physical gates for a complete, atomically published window.
#[derive(Debug, Clone, Copy)]
pub struct CarrierSupportConfig {
    /// `duration` is the complete window. Its momentum/work gates apply both
    /// to each original native step and to the accumulated complete window.
    pub step: ArticulatedContactStepConfig,
    /// Positive maximum native force/contact/drift step, in seconds.
    pub maximum_step: f64,
    /// Total attempted native steps, including refused geometry trials, 1..=4096.
    pub maximum_attempts: usize,
    /// Positive sphere/plane activation resolution, in metres. No negative gap
    /// is accepted; positive gaps can support impulses only within this band.
    pub activation_tolerance: f64,
    /// Total bounded observations of declared planes, 1..=1_000_000.
    pub maximum_geometry_queries: usize,
}

/// Complete native state and independently accumulated physical accounting.
#[derive(Debug, Clone)]
pub struct CarrierSupportResult {
    pub state_after: FreeFloatingMotionState,
    pub momentum_before: ArticulatedTrajectoryMomentum,
    pub momentum_after: ArticulatedTrajectoryMomentum,
    /// All original gravity, torque and support impulses, in N s.
    pub external_linear_impulse_world: Vec3,
    /// External angular impulse about `step.angular_origin_world`, in N m s.
    pub external_angular_impulse_about_origin_world: Vec3,
    /// Joint effort paired with native average pre/post-kick rates, in J.
    pub actuator_impulse_work: f64,
    /// Uniform gravity paired with native average pre/post-kick velocities, in J.
    pub gravity_impulse_work: f64,
    /// Held original link-body torque impulse work, in J.
    pub body_torque_impulse_work: f64,
    /// Work supplied by the prescribed translating supports, in J.
    pub support_work: f64,
    /// Native contact loss, including removal of unsupported free load kicks.
    /// This numerical splitting loss is not automatically frictional heat.
    pub contact_loss: f64,
    /// Held joint effort times actual joint displacement, in J.
    pub actuator_displacement_work: f64,
    /// Gravity times actual link-COM displacement, in J.
    pub gravity_displacement_work: f64,
    /// Sum of actual zero-load native drift energy defects, in J.
    pub drift_energy_defect: f64,
    /// Complete kinetic change minus impulse/support work plus contact loss, in J.
    pub kinetic_work_residual: f64,
    pub linear_momentum_residual: Vec3,
    pub angular_momentum_residual: Vec3,
    /// Accepted advances and all attempted native advances, respectively.
    pub accepted_steps: usize,
    pub attempted_steps: usize,
    pub geometry_queries: usize,
    /// Minimum observed affine-path, native endpoint and material-point gap.
    /// This floating observation does not enclose all fractional-time rounding.
    pub minimum_gap: f64,
}

/// The automatic geometry owner's complete-path refusal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CarrierSupportPathError {
    Invalid,
    Penetration { plane_id: u64, gap: f64 },
    GeometryBudget,
}

/// Refusal publishes no state, including after several successful private steps.
#[derive(Debug)]
pub enum CarrierSupportError {
    Carrier(CarrierPlaneError),
    Step(ArticulatedContactStepError<CarrierSupportPathError>),
    /// Only zero-restitution contact is admitted by this sustained lane.
    Restitution {
        plane_id: u64,
    },
    InitialGeometry,
    AttemptBudget,
    NoTimeProgress,
    WindowBalance,
}
impl core::fmt::Display for CarrierSupportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "supported articulated carrier: {self:?}")
    }
}
impl std::error::Error for CarrierSupportError {}
impl From<CarrierPlaneError> for CarrierSupportError {
    fn from(error: CarrierPlaneError) -> Self {
        Self::Carrier(error)
    }
}

/// Advance an articulated spherical carrier through inelastic landing and support.
///
/// The root-attached sphere encloses all fixed/revolute links, whose origins and
/// COMs coincide. The original carrier admission owns these restrictions. Gravity
/// supplies translation; finite body torques and internal actuator efforts retain
/// original coupled inertia and base recoil. Every declared plane must have zero
/// restitution. Bouncing impacts use the separate original event owner.
///
/// All nearby planes supply simultaneous circular Coulomb contact. The native
/// unforced drift then transports every body and joint. Its complete
/// real-arithmetic affine COM path and actual rounded endpoint are checked
/// against all planes, including inactive planes. A refused crossing halves
/// the attempted step, without publishing its
/// state; newly near planes receive regenerated radial material witnesses on the
/// next accepted step. No configuration, contact arm or velocity is projected.
///
/// The force/contact/drift scheme is first order. Activation is an explicitly
/// finite spatial resolution, not exact event time. This is floating analytic
/// sphere/plane coverage, not interval-certified CCD, unknown-shape admission,
/// general articulated collision detection, a no-slip projection, or a complete
/// Moreau--Jean solver. Fractional-time floating evaluations have no interval
/// enclosure here; small observed negative roundoff gaps still refuse.
/// Plane positions belong to this window's start and translate exactly once.
pub fn advance_supported_carrier_planes(
    system: ArticulatedTrajectorySystem<'_>,
    radius: f64,
    planes: &[CarrierPlane],
    config: CarrierSupportConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CarrierSupportResult, CarrierSupportError> {
    if cancelled() {
        return Err(CarrierSupportError::Step(
            ArticulatedContactStepError::Cancelled,
        ));
    }
    let original = CarrierPlaneConfig {
        trajectory: ArticulatedTrajectoryConfig {
            duration: config.step.duration,
            maximum_smooth_step: config.maximum_step,
            maximum_smooth_intervals: config.maximum_attempts,
            angular_origin_world: config.step.angular_origin_world,
            contact: config.step.contact,
        },
        activation_tolerance: config.activation_tolerance,
        maximum_events: 0,
        maximum_geometry_queries: config.maximum_geometry_queries,
    };
    let planes = admit(system, radius, planes, original)?;
    for plane in &planes {
        if plane.restitution != 0.0 {
            return Err(CarrierSupportError::Restitution { plane_id: plane.id });
        }
    }
    let mass = system
        .model
        .links()
        .iter()
        .fold(0.0, |total, link| total + link.inertia().mass());
    let mut state = system.state.clone();
    let mut work = Work {
        smooth: 0,
        geometry: 0,
        accepted: 0,
        minimum_gap: f64::INFINITY,
    };
    let mut time = 0.0;
    let mut result: Option<CarrierSupportResult> = None;
    while time < config.step.duration {
        let current = ArticulatedTrajectorySystem {
            state: &state,
            ..system
        };
        let (valid, gap, contacts) = observe(
            current,
            &state,
            &state,
            &planes,
            radius,
            time,
            0.0,
            true,
            original,
            &mut work,
            &mut cancelled,
        )?;
        if !valid {
            return Err(CarrierSupportError::InitialGeometry);
        }
        work.minimum_gap = work.minimum_gap.min(gap);
        let mut duration = (config.step.duration - time).min(config.maximum_step);
        let step = loop {
            if work.smooth == config.maximum_attempts {
                return Err(CarrierSupportError::AttemptBudget);
            }
            if !duration.is_finite() || duration <= 0.0 || time + duration <= time {
                return Err(CarrierSupportError::NoTimeProgress);
            }
            work.smooth += 1;
            let attempted = step_articulated_frictional_contacts(
                &[current],
                &contacts,
                ArticulatedContactStepConfig {
                    duration,
                    ..config.step
                },
                |path| admit_path(path, &planes, radius, mass, time, original, &mut work),
                &mut cancelled,
            );
            match attempted {
                Ok(step) => break step,
                Err(ArticulatedContactStepError::Path(CarrierSupportPathError::Penetration {
                    ..
                }))
                | Err(ArticulatedContactStepError::Penetration { endpoint: true, .. }) => {
                    duration *= 0.5;
                }
                Err(error) => return Err(CarrierSupportError::Step(error)),
            }
        };
        work.accepted += 1;
        work.minimum_gap = work.minimum_gap.min(step.path_admission);
        let total = result.get_or_insert_with(|| CarrierSupportResult {
            state_after: state.clone(),
            momentum_before: step.momentum_before,
            momentum_after: step.momentum_after,
            external_linear_impulse_world: Vec3::default(),
            external_angular_impulse_about_origin_world: Vec3::default(),
            actuator_impulse_work: 0.0,
            gravity_impulse_work: 0.0,
            body_torque_impulse_work: 0.0,
            support_work: 0.0,
            contact_loss: 0.0,
            actuator_displacement_work: 0.0,
            gravity_displacement_work: 0.0,
            drift_energy_defect: 0.0,
            kinetic_work_residual: 0.0,
            linear_momentum_residual: Vec3::default(),
            angular_momentum_residual: Vec3::default(),
            accepted_steps: 0,
            attempted_steps: 0,
            geometry_queries: 0,
            minimum_gap: f64::INFINITY,
        });
        total.momentum_after = step.momentum_after;
        total.external_linear_impulse_world =
            total.external_linear_impulse_world + step.external_linear_impulse_world;
        total.external_angular_impulse_about_origin_world = total
            .external_angular_impulse_about_origin_world
            + step.external_angular_impulse_about_origin_world;
        total.actuator_impulse_work += step.actuator_impulse_work;
        total.gravity_impulse_work += step.gravity_impulse_work;
        total.body_torque_impulse_work += step.body_wrench_impulse_work;
        total.support_work += step.contact.plane_work;
        total.contact_loss += step.contact.contact.dissipated_energy;
        total.actuator_displacement_work += step.actuator_displacement_work;
        total.gravity_displacement_work += step.gravity_displacement_work;
        total.drift_energy_defect += step.drift_energy_defect;
        state = step
            .states_after
            .into_iter()
            .next()
            .expect("one admitted carrier");
        time += duration;
    }
    let mut result = result.ok_or(CarrierSupportError::NoTimeProgress)?;
    result.kinetic_work_residual = result.momentum_after.kinetic_energy
        - result.momentum_before.kinetic_energy
        - result.actuator_impulse_work
        - result.gravity_impulse_work
        - result.body_torque_impulse_work
        - result.support_work
        + result.contact_loss;
    result.linear_momentum_residual = result.momentum_after.linear_world
        - result.momentum_before.linear_world
        - result.external_linear_impulse_world;
    result.angular_momentum_residual = result.momentum_after.angular_about_origin_world
        - result.momentum_before.angular_about_origin_world
        - result.external_angular_impulse_about_origin_world;
    if ![
        result.actuator_impulse_work,
        result.gravity_impulse_work,
        result.body_torque_impulse_work,
        result.support_work,
        result.contact_loss,
        result.actuator_displacement_work,
        result.gravity_displacement_work,
        result.drift_energy_defect,
        result.kinetic_work_residual,
    ]
    .into_iter()
    .all(f64::is_finite)
        || !finite_vector(result.linear_momentum_residual)
        || !finite_vector(result.angular_momentum_residual)
        || result.kinetic_work_residual.abs() > config.step.energy_tolerance
        || result.drift_energy_defect.abs() > config.step.drift_energy_tolerance
        || result.linear_momentum_residual.norm() > config.step.linear_momentum_tolerance
        || result.angular_momentum_residual.norm() > config.step.angular_momentum_tolerance
    {
        return Err(CarrierSupportError::WindowBalance);
    }
    if cancelled() {
        return Err(CarrierSupportError::Step(
            ArticulatedContactStepError::Cancelled,
        ));
    }
    result.state_after = state;
    result.accepted_steps = work.accepted;
    result.attempted_steps = work.smooth;
    result.geometry_queries = work.geometry;
    result.minimum_gap = work.minimum_gap;
    Ok(result)
}

fn admit_path(
    path: &ArticulatedContactPath<'_>,
    planes: &[CarrierPlane],
    radius: f64,
    mass: f64,
    time: f64,
    config: CarrierPlaneConfig,
    work: &mut Work,
) -> Result<f64, CarrierSupportPathError> {
    let first = &path.post_contact[0];
    let last = &path.endpoint[0];
    // Reproduce MomentumMap::coordinates and step_free_floating's native COM
    // arithmetic for this admitted zero-COM, zero-shape-linear-momentum model:
    // X(h) = X0 + R(M*v_body)*(h/M). There is no joint/attitude curvature in
    // these coordinates, though fractional floating evaluations can round away
    // from the real-arithmetic affine path. Do not substitute (R*v_body)*h.
    let momentum = first
        .base
        .world_from_base
        .rotation()
        .rotate(first.base.twist_body.linear.scale(mass))
        .map_err(|_| CarrierSupportPathError::Invalid)?;
    let affine_endpoint =
        first.base.world_from_base.translation() + momentum.scale(path.duration / mass);
    let mut minimum = f64::INFINITY;
    for plane in planes {
        work.query(config)
            .map_err(|_| CarrierSupportPathError::GeometryBudget)?;
        let normal = plane_normal(plane).map_err(|_| CarrierSupportPathError::Invalid)?;
        let first_plane =
            plane_position(plane, time).map_err(|_| CarrierSupportPathError::Invalid)?;
        let last_plane = plane_position(plane, time + path.duration)
            .map_err(|_| CarrierSupportPathError::Invalid)?;
        let start = normal.dot(first.base.world_from_base.translation() - first_plane) - radius;
        let end = normal.dot(affine_endpoint - last_plane) - radius;
        let actual = normal.dot(last.base.world_from_base.translation() - last_plane) - radius;
        if ![start, end, actual].into_iter().all(f64::is_finite) {
            return Err(CarrierSupportPathError::Invalid);
        }
        let gap = start.min(end).min(actual);
        if gap < 0.0 {
            return Err(CarrierSupportPathError::Penetration {
                plane_id: plane.id,
                gap,
            });
        }
        minimum = minimum.min(gap);
    }
    Ok(minimum)
}

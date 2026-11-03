//! Simultaneous finite ejection at specified retained-body-relative velocities.
//!
//! The original body is the rigidly locked mass distribution immediately before
//! release. Each point parcel's final velocity relative to the retained rigid
//! velocity field is prescribed. All jets therefore share the same translation
//! and rotation recoil; independent reduced-mass impulses are generally wrong.

use super::{
    PointMassRelease, PointMassReleaseBatchResult, relative_momentum, release_point_masses,
};
use crate::assembly::{
    AssemblyControls, AssemblyError, Sum, VectorSum, finite, maximum, poll, vector,
};
use crate::{MassProperties, RigidBodyState, Vec3};

/// Material already included in the original body's mass and full inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassEjection {
    /// Nonnegative point-parcel mass [kg]. Zero mass produces no parcel.
    pub mass_kg: f64,
    /// Original COM-relative location, in original body axes [m].
    pub attachment_arm_body_m: Vec3,
    /// Final parcel velocity MINUS the retained rigid velocity at that same
    /// point [m/s], in original body axes. This includes the retained body's
    /// simultaneous translation and rotation recoil, not its old velocity.
    /// Finite values on a zero-mass descriptor are validated but inactive.
    pub relative_velocity_body_m_s: Vec3,
}

/// Native event gates plus the velocity boundary condition and energy source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassEjectionControls {
    /// Original mass/impulse budgets and actual momentum/energy balance gates.
    pub assembly: AssemblyControls,
    /// Maximum component of actual relative-velocity error in body axes [m/s].
    pub relative_velocity_tolerance_m_s: f64,
    /// Available mechanical separation energy [J], finite and nonnegative.
    /// The composing source must debit the accepted native supplied work.
    pub maximum_supplied_work_j: f64,
}

/// One simultaneous velocity-controlled event, with the original native result.
#[derive(Clone, Debug, PartialEq)]
pub struct PointMassEjectionBatchResult {
    /// Native retained state, parcels, COM shift, work and actual balances.
    pub release: PointMassReleaseBatchResult,
    /// Solved native impulse descriptors in original input order, including
    /// zero-mass entries. This reuses the staging storage needed by release.
    pub impulses: Vec<PointMassRelease>,
    /// Largest component of the actual returned relative-velocity error [m/s].
    pub maximum_relative_velocity_residual_m_s: f64,
}

/// No failed ejection publishes changed bodies, parcels or source work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PointMassEjectionError {
    /// Original native admission, budget, balance or cancellation failure.
    Assembly(AssemblyError),
    /// The actual returned canonical states do not realize a requested velocity.
    RelativeVelocity {
        /// Original descriptor index.
        input_index: usize,
        /// Actual minus requested relative velocity, in body axes [m/s].
        residual_body_m_s: Vec3,
        /// Declared maximum absolute component [m/s].
        tolerance_m_s: f64,
    },
    /// The physically required separation work exceeds its available source.
    SourceWork {
        /// Native jointly computed required work [J].
        required_j: f64,
        /// Declared available work [J].
        maximum_j: f64,
    },
}

impl From<AssemblyError> for PointMassEjectionError {
    fn from(error: AssemblyError) -> Self {
        Self::Assembly(error)
    }
}
impl core::fmt::Display for PointMassEjectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid point ejection: {self:?}")
    }
}
impl std::error::Error for PointMassEjectionError {}

/// Release all point parcels with prescribed final relative velocities.
///
/// Let T and I be the ORIGINAL total mass and full COM inertia, r_i the original
/// arms and c_i the prescribed final relative velocities. Conservation gives
/// a common rigid-field change dv=-sum(m_i*c_i)/T and
/// dw=-I^-1 sum(r_i cross (m_i*c_i)). Thus the required impulses are
/// J_i=m_i*(c_i+dv+dw cross r_i). Using the original composite inertia includes
/// every shared-body cross term without a dense 3N constraint matrix.
///
/// The velocity projection uses fixed-order compensated O(inputs) reductions
/// and one original angular-mobility action. The native release owner still
/// admits the retained COM/inertia, applies all impulses and computes their
/// joint nonnegative work. At most O(inputs) staging/output storage is needed.
/// Every descriptor, including zero mass, consumes one native part slot; every
/// NONZERO SOLVED impulse consumes an impulse slot. A parcel with c_i=0 may
/// still need impulse to follow the recoiling body's local velocity.
///
/// The final velocity gate uses actual stored momenta in a common translating
/// chart and actual returned parcel/body positions. A large common boost cannot
/// hide an unrepresented ejection jump. Native conservation, this velocity gate,
/// the source-work limit and cancellation must all pass before publication.
///
/// This is an instantaneous point-detachment boundary condition. It neither
/// integrates a continuous rocket/nozzle law nor supplies pressure thrust,
/// propellant thermodynamics, fracture energy or a finite parcel inertia.
pub fn eject_point_masses(
    body: &RigidBodyState,
    properties: MassProperties,
    outgoing: &[PointMassEjection],
    controls: PointMassEjectionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PointMassEjectionBatchResult, PointMassEjectionError> {
    poll(&mut cancelled)?;
    let parts = outgoing.len().checked_add(1).ok_or(AssemblyError::Invalid(
        "point-ejection input count overflow",
    ))?;
    controls.assembly.validate(parts, 0)?;
    if [
        controls.relative_velocity_tolerance_m_s,
        controls.maximum_supplied_work_j,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value < 0.0)
    {
        return Err(AssemblyError::Invalid(
            "finite nonnegative ejection velocity tolerance and source work required",
        )
        .into());
    }
    properties.validate().map_err(AssemblyError::from)?;
    let total = properties.mass();
    let mut removed = Sum::default();
    let mut anchor_mass = 0.0;
    let mut anchor = Vec3::ZERO;
    for input in outgoing {
        poll(&mut cancelled)?;
        if !input.mass_kg.is_finite() || input.mass_kg < 0.0 {
            return Err(AssemblyError::Invalid("finite nonnegative ejected mass required").into());
        }
        vector(input.attachment_arm_body_m)?;
        vector(input.relative_velocity_body_m_s)?;
        removed.add(input.mass_kg);
        if input.mass_kg > anchor_mass {
            anchor_mass = input.mass_kg;
            anchor = input.relative_velocity_body_m_s;
        }
    }
    let removed = finite(removed.value())?;
    let retained = finite(total - removed)?;
    if anchor_mass > 0.0 && (retained <= 0.0 || retained >= total) {
        return Err(AssemblyError::Invalid(
            "positive represented retained and ejected masses required",
        )
        .into());
    }

    // c_i - mean(c) = (c_i-anchor) + anchor*M/T
    //                  - sum(m_j*(c_j-anchor))/T.
    // This retains the small reduced-mass response near total depletion.
    let mut translation = VectorSum::default();
    translation.add(vector(anchor.scale(retained / total))?);
    let mut angular = VectorSum::default();
    for input in outgoing {
        poll(&mut cancelled)?;
        if input.mass_kg == 0.0 {
            continue;
        }
        let relative = vector(input.relative_velocity_body_m_s.sub(anchor))?;
        translation.add(vector(relative.scale(-input.mass_kg / total))?);
        if input.attachment_arm_body_m != Vec3::ZERO {
            let carried = vector(input.relative_velocity_body_m_s.scale(input.mass_kg))?;
            angular.add(vector(input.attachment_arm_body_m.cross(carried))?);
        }
    }
    let translation = vector(translation.value())?;
    let rotation = vector(
        properties
            .angular_velocity_body_checked(vector(angular.value())?)
            .map_err(AssemblyError::from)?
            .scale(-1.0),
    )?;
    let mut impulses = Vec::with_capacity(outgoing.len());
    for input in outgoing {
        poll(&mut cancelled)?;
        let impulse = if input.mass_kg == 0.0 {
            Vec3::ZERO
        } else {
            let mut change = VectorSum::default();
            change.add(vector(input.relative_velocity_body_m_s.sub(anchor))?);
            change.add(translation);
            change.add(vector(rotation.cross(input.attachment_arm_body_m))?);
            vector(vector(change.value())?.scale(input.mass_kg))?
        };
        impulses.push(PointMassRelease {
            mass_kg: input.mass_kg,
            attachment_arm_body_m: input.attachment_arm_body_m,
            impulse_on_parcel_body_ns: impulse,
        });
    }
    let release = release_point_masses(
        body,
        properties,
        &impulses,
        controls.assembly,
        &mut cancelled,
    )?;
    if release.supplied_work_j > controls.maximum_supplied_work_j {
        return Err(PointMassEjectionError::SourceWork {
            required_j: release.supplied_work_j,
            maximum_j: controls.maximum_supplied_work_j,
        });
    }

    // Subtract the same boost before division: the audit concerns actual
    // canonical momenta, not two independently rounded laboratory velocities.
    let boost = body
        .center_of_mass_velocity_world(properties)
        .map_err(AssemblyError::from)?;
    let retained_velocity = velocity_in_chart(
        release.state.linear_momentum_world(),
        release.properties.mass(),
        boost,
    )?;
    let orientation = release.state.pose().orientation();
    let omega = release
        .properties
        .angular_velocity_body_checked(release.state.angular_momentum_body())
        .map_err(AssemblyError::from)?;
    let mut maximum_error: f64 = 0.0;
    let active = outgoing
        .iter()
        .enumerate()
        .filter(|(_, input)| input.mass_kg > 0.0);
    for ((input_index, input), parcel) in active.zip(&release.parcels) {
        poll(&mut cancelled)?;
        let arm_world = vector(
            parcel
                .position_world_m
                .sub(release.state.pose().position_world()),
        )?;
        let arm = vector(orientation.rotate_world_to_body(arm_world))?;
        let relative = vector(
            orientation.rotate_world_to_body(vector(
                velocity_in_chart(parcel.linear_momentum_world_ns, parcel.mass_kg, boost)?
                    .sub(retained_velocity),
            )?),
        )?;
        let residual = vector(
            relative
                .sub(omega.cross(arm))
                .sub(input.relative_velocity_body_m_s),
        )?;
        maximum_error = maximum_error.max(maximum(residual));
        if maximum(residual) > controls.relative_velocity_tolerance_m_s {
            return Err(PointMassEjectionError::RelativeVelocity {
                input_index,
                residual_body_m_s: residual,
                tolerance_m_s: controls.relative_velocity_tolerance_m_s,
            });
        }
    }
    poll(&mut cancelled)?;
    Ok(PointMassEjectionBatchResult {
        release,
        impulses,
        maximum_relative_velocity_residual_m_s: maximum_error,
    })
}

fn velocity_in_chart(momentum: Vec3, mass: f64, boost: Vec3) -> Result<Vec3, AssemblyError> {
    let relative = relative_momentum(momentum, mass, boost)?;
    vector(Vec3::new(
        relative.x / mass,
        relative.y / mass,
        relative.z / mass,
    ))
}

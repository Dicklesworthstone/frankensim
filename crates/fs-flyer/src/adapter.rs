//! fs-mbd ↔ fs-geom type-adapter seam (bead wf-root-guzez.4.2.3,
//! E3.2-iii). The spine speaks raw arrays ([f64; 3]/[f64; 4]); fs-mbd
//! speaks `Vec3`/`UnitQuaternion` (canonical double-cover representative,
//! bit-preserving transport constructor); fs-geom speaks `Point3`/`Vec3`.
//! This module is the ONLY sanctioned crossing.
//!
//! The load-bearing subtlety: fs-mbd's `from_canonical_components` REFUSES
//! the negative representative of the quaternion double cover, while both
//! representatives describe the same physical rotation. On REPLAY paths
//! bit-preservation is mandatory, so the adapter offers explicit transports:
//!
//! - [`quat_to_mbd_canonical`] — bit-preserving; the negative
//!   representative is a TYPED refusal naming the repair.
//! - [`quat_to_mbd_normalizing`] — explicit renormalize+resign via
//!   `UnitQuaternion::new` for presentation/geometry consumers, where the
//!   double-cover sign is meaningless.
//!
//! The physical observer applies the spine's group admission first, then
//! explicitly normalizes only its native observation copy. This accepts the
//! spine's admitted unit-norm roundoff without changing source replay bits
//! or admitting a non-unit source quaternion.

use crate::Refusal;
use crate::spine::{Loads, RigidBody, SixDofState, rotation_from_boundary};
use fs_geom::Point3;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, Pose, PreparedRigidPointMotion, RigidBodyIntegrator,
    RigidBodyState, UnitQuaternion, Vec3 as MbdVec3, Wrench,
};

/// Raw world-position triple → fs-geom `Point3` (bit-preserving).
#[must_use]
pub const fn pos_to_geom(pos_m: [f64; 3]) -> Point3 {
    Point3::new(pos_m[0], pos_m[1], pos_m[2])
}

/// fs-geom `Point3` → raw triple (bit-preserving).
#[must_use]
pub const fn geom_to_pos(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

/// Raw triple → fs-mbd `Vec3` (bit-preserving).
#[must_use]
pub const fn vec_to_mbd(v: [f64; 3]) -> MbdVec3 {
    MbdVec3 {
        x: v[0],
        y: v[1],
        z: v[2],
    }
}

/// fs-mbd `Vec3` → raw triple (bit-preserving).
#[must_use]
pub const fn mbd_to_vec(v: MbdVec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// Spine quaternion → fs-mbd `UnitQuaternion`, BIT-PRESERVING (replay
/// path). Uses the transport constructor; norm drift beyond fs-mbd's
/// admission or a negative double-cover representative is a typed refusal
/// — never a silent renormalization on this path.
///
/// # Errors
/// `quaternion-not-canonical` with both repairs ranked.
pub fn quat_to_mbd_canonical(quat: [f64; 4]) -> Result<UnitQuaternion, Refusal> {
    UnitQuaternion::from_canonical_components(quat).map_err(|e| Refusal {
        code: "quaternion-not-canonical",
        message: format!(
            "spine quaternion {quat:?} refused by the canonical transport constructor: {e:?}"
        ),
        ranked_repairs: vec![
            "for presentation/geometry use quat_to_mbd_normalizing (sign is meaningless there)"
                .into(),
            "on replay paths a non-canonical quaternion means state corruption — verify the \
             digest trace"
                .into(),
        ],
    })
}

/// Spine quaternion → fs-mbd `UnitQuaternion` via explicit renormalize +
/// canonical resign (presentation/geometry path; NOT bit-preserving).
///
/// # Errors
/// `quaternion-invalid` (non-finite or zero).
pub fn quat_to_mbd_normalizing(quat: [f64; 4]) -> Result<UnitQuaternion, Refusal> {
    UnitQuaternion::new(quat[0], quat[1], quat[2], quat[3]).map_err(|e| Refusal {
        code: "quaternion-invalid",
        message: format!("spine quaternion {quat:?} is not normalizable: {e:?}"),
        ranked_repairs: vec!["a NaN/zero quaternion upstream means the integrator diverged".into()],
    })
}

/// fs-mbd `UnitQuaternion` → raw components (bit-preserving).
#[must_use]
pub const fn mbd_to_quat(q: UnitQuaternion) -> [f64; 4] {
    q.components()
}

/// One-stop view of a spine state for fs-mbd consumers (replay-strict
/// quaternion path).
///
/// # Errors
/// As [`quat_to_mbd_canonical`].
pub fn state_to_mbd(
    state: &SixDofState,
) -> Result<(MbdVec3, MbdVec3, UnitQuaternion, MbdVec3), Refusal> {
    Ok((
        vec_to_mbd(state.pos_m),
        vec_to_mbd(state.vel_mps),
        quat_to_mbd_canonical(state.quat)?,
        vec_to_mbd(state.omega_body),
    ))
}

/// Rebuild a spine state from fs-mbd pieces (bit-preserving).
#[must_use]
pub const fn mbd_to_state(
    pos: MbdVec3,
    vel: MbdVec3,
    quat: UnitQuaternion,
    omega: MbdVec3,
) -> SixDofState {
    SixDofState {
        pos_m: mbd_to_vec(pos),
        vel_mps: mbd_to_vec(vel),
        quat: mbd_to_quat(quat),
        omega_body: mbd_to_vec(omega),
    }
}

impl RigidBody {
    /// Observe smooth six-DOF motion through the native rigid-body owner.
    ///
    /// `state` and `loads` must describe the same accepted instant. As in the
    /// spine, force is WORLD/NED and moment is BODY/FRD about the COM; the
    /// supplied force already includes any gravity and support reactions.
    /// Native gravity is therefore zero here, so weight is never added twice.
    /// The returned preparation shares full-inertia Euler acceleration across
    /// any number of body-fixed or prescribed moving point queries. Its point
    /// acceleration includes gravity; subtract the declared physical gravity
    /// separately to obtain an ideal accelerometer's specific force.
    ///
    /// The existing principal frame is transported without diagonalizing the
    /// tensor again. This read-only query allocates nothing on success and
    /// changes neither motion nor the spine's two-call load schedule. Do not
    /// feed it a predictor's loads paired with a corrected endpoint state.
    /// Quaternion admission is shared with the spine: both double-cover
    /// representatives and admitted unit-norm roundoff describe valid poses.
    /// Only the observation's native copy is normalized to match fs-mbd's
    /// tighter unit representation; the supplied state remains bit-for-bit
    /// unchanged and strict replay transport keeps its original semantics.
    ///
    /// # Errors
    /// Native state, wrench and derived-motion refusals are mapped to
    /// `rigid-motion-refused`; invalid source quaternions use the spine's
    /// `orientation-outside-group` refusal. This describes a rigid body with
    /// physically consistent principal inertia. Reduced rotor-carrier inertia,
    /// added fluid inertia, moving joints and impacts need their coupled owner.
    pub fn prepare_point_motion(
        &self,
        state: &SixDofState,
        loads: Loads,
    ) -> Result<PreparedRigidPointMotion, Refusal> {
        self.admit()?;
        let tensor = self.inertia_tensor();
        let principal = tensor.principal_moments();
        let axes = tensor.principal_to_body().as_quat();
        let properties = MassProperties::from_principal_frame(
            self.mass_kg,
            MbdVec3::ZERO,
            MbdVec3::new(principal.x, principal.y, principal.z),
            quat_to_mbd_canonical([axes.w, axes.x, axes.y, axes.z])?,
        )
        .map_err(motion_refusal)?;
        let rotation = rotation_from_boundary(state.quat)?;
        let quaternion = rotation.as_quat();
        let orientation =
            quat_to_mbd_normalizing([quaternion.w, quaternion.x, quaternion.y, quaternion.z])?;
        let position = vec_to_mbd(state.pos_m);
        let velocity = vec_to_mbd(state.vel_mps);
        let omega = vec_to_mbd(state.omega_body);
        let native = RigidBodyState::new(
            Pose::new(position, orientation).map_err(motion_refusal)?,
            velocity.scale(self.mass_kg),
            properties
                .angular_momentum_body_checked(omega)
                .map_err(motion_refusal)?,
        )
        .map_err(motion_refusal)?;
        RigidBodyIntegrator::new(Gravity::ZERO)
            .prepare_point_motion(
                native,
                properties,
                Wrench {
                    force_world: vec_to_mbd(loads.force_n),
                    torque_body: vec_to_mbd(loads.moment_nm),
                },
            )
            .map_err(motion_refusal)
    }
}

fn motion_refusal(error: DynamicsError) -> Refusal {
    Refusal {
        code: "rigid-motion-refused",
        message: format!("native rigid-point observation refused: {error}"),
        ranked_repairs: vec![
            "supply matching finite state and total loads for an ordinary rigid body".into(),
        ],
    }
}

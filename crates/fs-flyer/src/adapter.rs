//! fs-mbd ↔ fs-geom type-adapter seam (bead wf-root-guzez.4.2.3,
//! E3.2-iii). The spine speaks raw arrays ([f64; 3]/[f64; 4]); fs-mbd
//! speaks `Vec3`/`UnitQuaternion` (canonical double-cover representative,
//! bit-preserving transport constructor); fs-geom speaks `Point3`/`Vec3`.
//! This module is the ONLY sanctioned crossing.
//!
//! The load-bearing subtlety: fs-mbd's `from_canonical_components` REFUSES
//! the negative representative of the quaternion double cover, while the
//! spine's Lie step can walk a trajectory into that half. On REPLAY paths
//! bit-preservation is mandatory, so the adapter offers exactly two
//! spellings and no silent middle ground:
//!
//! - [`quat_to_mbd_canonical`] — bit-preserving; the negative
//!   representative is a TYPED refusal naming the repair.
//! - [`quat_to_mbd_normalizing`] — explicit renormalize+resign via
//!   `UnitQuaternion::new` for presentation/geometry consumers, where the
//!   double-cover sign is meaningless.

use crate::Refusal;
use crate::spine::SixDofState;
use fs_geom::Point3;
use fs_mbd::{UnitQuaternion, Vec3 as MbdVec3};

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

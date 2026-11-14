//! Exact native rigid-body words and complete numerical-law identities shared
//! by the Component adapters. Component owners additionally admit their energy
//! and memory accounts before publishing or accepting these words.

use crate::spatial_loads::internal_impulse::viscous::{
    PassiveBushingMatrix, ViscousConnectorControls,
};
use crate::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

/// Position, canonical quaternion, linear momentum, angular momentum.
pub(crate) const BODY_WORDS: usize = 13;

fn body_values(state: RigidBodyState) -> [f64; BODY_WORDS] {
    let p = state.pose().position_world();
    let [w, x, y, z] = state.pose().orientation().components();
    let linear = state.linear_momentum_world();
    let angular = state.angular_momentum_body();
    [
        p.x, p.y, p.z, w, x, y, z, linear.x, linear.y, linear.z, angular.x, angular.y, angular.z,
    ]
}

/// Encode no count prefix. A refusal leaves the caller's word vector intact.
pub(crate) fn encode_bodies(states: &[RigidBodyState], words: &mut Vec<u64>) -> bool {
    if states.iter().any(|state| {
        body_values(*state).iter().any(|value| !value.is_finite())
            || UnitQuaternion::from_canonical_components(state.pose().orientation().components())
                .is_err()
    }) {
        return false;
    }
    for state in states {
        words.extend(body_values(*state).map(f64::to_bits));
    }
    true
}

/// Require exactly the admitted body's count; never normalize a saved rotation.
pub(crate) fn decode_bodies(words: &[u64], count: usize) -> Option<Vec<RigidBodyState>> {
    if count.checked_mul(BODY_WORDS)? != words.len() {
        return None;
    }
    words
        .chunks_exact(BODY_WORDS)
        .map(|row| {
            let value = |index| f64::from_bits(row[index]);
            let vector = |start| Vec3::new(value(start), value(start + 1), value(start + 2));
            let orientation =
                UnitQuaternion::from_canonical_components([value(3), value(4), value(5), value(6)])
                    .ok()?;
            let pose = Pose::new(vector(0), orientation).ok()?;
            RigidBodyState::new(pose, vector(7), vector(10)).ok()
        })
        .collect()
}

pub(crate) fn push_vec3(value: Vec3, words: &mut Vec<u64>) {
    words.extend([value.x, value.y, value.z].map(f64::to_bits));
}

/// Original principal representation, including its fixed body-frame rotation.
pub(crate) fn push_properties(properties: &[MassProperties], words: &mut Vec<u64>) {
    words.push(properties.len() as u64);
    for properties in properties {
        words.push(properties.mass().to_bits());
        push_vec3(properties.center_of_mass_body(), words);
        push_vec3(properties.principal_inertia_body(), words);
        words.extend(
            properties
                .principal_axes_to_body()
                .components()
                .map(f64::to_bits),
        );
    }
}

/// Retain supplied factors and evaluation mode as well as the native SI form.
pub(crate) fn push_matrix(matrix: &PassiveBushingMatrix, words: &mut Vec<u64>) {
    matrix.append_checkpoint_identity(words);
}

pub(crate) fn push_source_controls(controls: ViscousConnectorControls, words: &mut Vec<u64>) {
    push_vec3(controls.impulse.angular_origin_world_m, words);
    words.extend(
        [
            controls.impulse.linear_tolerance_ns,
            controls.impulse.angular_tolerance_nms,
            controls.impulse.work_tolerance_j,
            controls.linear_speed_tolerance_m_s,
            controls.angular_speed_tolerance_rad_s,
        ]
        .map(f64::to_bits),
    );
}

pub(crate) fn push_torques(torques: &[super::TorquePort], words: &mut Vec<u64>) {
    words.push(torques.len() as u64);
    for port in torques {
        words.extend([port.body as u64, port.input as u64, port.output as u64]);
    }
}

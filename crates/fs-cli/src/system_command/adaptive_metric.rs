//! Complete native network state accounts for the generic step-doubling driver.

use fs_math::det::{atan2, hypot};
use fs_mbd::component::relaxation::RelaxationState;
use fs_mbd::{MassProperties, RigidBodyState};
use fs_system::{Checkpoint, ErrorNorm, ErrorTolerance, NodeId};

#[derive(Clone, Copy)]
pub(super) enum BodyQuantity {
    Position,
    Velocity,
    AngularVelocity,
    Orientation,
}

impl BodyQuantity {
    pub(super) const ACCOUNTS: [(&'static str, Self); 4] = [
        ("position_m", Self::Position),
        ("velocity_m_s", Self::Velocity),
        ("angular_velocity_rad_s", Self::AngularVelocity),
        ("orientation_rad", Self::Orientation),
    ];

    fn vector(self, body: &RigidBodyState, properties: MassProperties) -> Option<[f64; 3]> {
        let vector = match self {
            Self::Position => body.pose().position_world(),
            Self::Velocity => body.center_of_mass_velocity_world(properties).ok()?,
            Self::AngularVelocity => properties
                .angular_velocity_body_checked(body.angular_momentum_body())
                .ok()?,
            Self::Orientation => return None,
        };
        Some([vector.x, vector.y, vector.z])
    }
}

#[derive(Clone, Copy)]
pub(super) enum NetworkCoordinate {
    Body {
        body: usize,
        properties: MassProperties,
        quantity: BodyQuantity,
    },
    Branch {
        connector: usize,
        branch: usize,
    },
}

impl NetworkCoordinate {
    pub(super) fn include(
        self,
        norm: &mut ErrorNorm,
        node: NodeId,
        initial: &Checkpoint,
        coarse: &Checkpoint,
        refined: &Checkpoint,
        tolerance: ErrorTolerance,
    ) {
        let (error, scale) = self
            .estimate(node, initial, coarse, refined)
            .unwrap_or((f64::NAN, f64::NAN));
        norm.include_error(error, scale, tolerance);
    }

    fn estimate(
        self,
        node: NodeId,
        initial: &Checkpoint,
        coarse: &Checkpoint,
        refined: &Checkpoint,
    ) -> Option<(f64, f64)> {
        match self {
            Self::Body {
                body,
                properties,
                quantity,
            } => {
                let initial = read_body(initial, node, body)?;
                let coarse = read_body(coarse, node, body)?;
                let refined = read_body(refined, node, body)?;
                if matches!(quantity, BodyQuantity::Orientation) {
                    Some((
                        rotation_distance(
                            coarse.pose().orientation().components(),
                            refined.pose().orientation().components(),
                        ),
                        // Relative angular tolerance follows this interval's
                        // motion, not the arbitrary orientation of the world.
                        rotation_distance(
                            initial.pose().orientation().components(),
                            refined.pose().orientation().components(),
                        ),
                    ))
                } else {
                    Some(vector_error(
                        quantity.vector(initial, properties)?,
                        quantity.vector(coarse, properties)?,
                        quantity.vector(refined, properties)?,
                    ))
                }
            }
            Self::Branch { connector, branch } => {
                let read = |checkpoint: &Checkpoint| {
                    checkpoint
                        .state::<RelaxationState>(node)?
                        .memory
                        .get(connector)?
                        .get(branch)
                        .copied()
                };
                // Original signed factor amplitudes include all six entries,
                // even unloaded rows. Equal branch energy is not equal state.
                Some(vector_error(read(initial)?, read(coarse)?, read(refined)?))
            }
        }
    }
}

fn read_body(checkpoint: &Checkpoint, node: NodeId, body: usize) -> Option<&RigidBodyState> {
    checkpoint
        .state::<Vec<RigidBodyState>>(node)
        .map(Vec::as_slice)
        .or_else(|| {
            checkpoint
                .state::<RelaxationState>(node)
                .map(|state| state.bodies.as_slice())
        })?
        .get(body)
}

fn vector_error<const N: usize>(
    initial: [f64; N],
    coarse: [f64; N],
    refined: [f64; N],
) -> (f64, f64) {
    if initial
        .iter()
        .chain(&coarse)
        .chain(&refined)
        .any(|x| !x.is_finite())
    {
        return (f64::NAN, f64::NAN);
    }
    let difference: [f64; N] = core::array::from_fn(|i| coarse[i] - refined[i]);
    let (error, start, end) = (norm(&difference), norm(&initial), norm(&refined));
    if !(error.is_finite() && start.is_finite() && end.is_finite()) {
        return (f64::NAN, f64::NAN);
    }
    (error, start.max(end))
}

fn norm(values: &[f64]) -> f64 {
    if values.iter().any(|x| !x.is_finite()) {
        return f64::NAN;
    }
    values.iter().fold(0.0, |length, &x| hypot(length, x))
}

fn normalized(q: [f64; 4]) -> Option<[f64; 4]> {
    let length = norm(&q);
    (length.is_finite() && length > 0.0).then(|| q.map(|x| x / length))
}

/// Shortest SO(3) angle, including the quaternion double cover and permitted
/// norm drift. Only metric-local copies are normalized; accepted bits stay put.
fn rotation_distance(first: [f64; 4], second: [f64; 4]) -> f64 {
    let (Some(first), Some(second)) = (normalized(first), normalized(second)) else {
        return f64::NAN;
    };
    let minus: [f64; 4] = core::array::from_fn(|i| first[i] - second[i]);
    let plus: [f64; 4] = core::array::from_fn(|i| first[i] + second[i]);
    let (minus, plus) = (norm(&minus), norm(&plus));
    if !(minus.is_finite() && plus.is_finite()) {
        return f64::NAN;
    }
    // Swapping either quaternion's sign exchanges these distances exactly.
    // The chord formula resolves tiny rotations without subtracting near-one
    // dot products, and stays continuous across canonical half-turn signs.
    4.0 * atan2(minus.min(plus), minus.max(plus))
}

#[cfg(test)]
mod tests;

//! Spatial loads reduced about the ORIGINAL rod material nodes.
//!
//! A point `x + R r` carries world force `F` and nodal moment `(R r) x F`.
//! The node is not its inertia centroid. The native rod kick already performs
//! that distinct moment transport, once, together with supported reactions.
//! This adapter reuses the mass-independent `fs-mbd` resultant/first-moment
//! reducer with the material node as its reference point; it never invokes a
//! rigid-body integrator or substitutes a new inertia model.

use fs_exec::CancelGate;
use fs_mbd::spatial_loads::{
    MAX_WORLD_POINT_FORCES, PreparedSpatialLoads, SpatialLoadError, WorldPointForce,
};
use fs_mbd::{UnitQuaternion, Vec3, Wrench};
use fs_solid::rod::MAX_MATERIAL_ROD_SEGMENTS;

use super::{RodMotionError, RodNodalWrench, So3, SpatialRodError, checkpoint, finite};

/// A force fixed in world axes and attached to a rotating material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodPointForce {
    /// Node whose ORIGINAL material frame owns the attachment.
    pub node: usize,
    /// Offset from the material node, not the mass centroid, in body axes [m].
    pub offset_body_m: [f64; 3],
    /// Force held in world axes [N].
    pub force_world_n: [f64; 3],
}

/// Immutable spatial loads with constant storage and evaluation cost per node.
///
/// Cold preparation is linear in nodes plus points, retains the supplied order
/// within each node, and uses the existing compensated first-moment reduction.
/// Both native momentum kicks evaluate the rotating arms and body couples at
/// their actual frames. Attachment count does not enter the step's hot loop.
///
/// Point forces are held in the world; body-held and world-held free couples
/// keep their distinct meanings. These prescribed loads are not a drag,
/// friction, spring, contact, motor circuit or other velocity/state feedback
/// law. Replacing such a law by a frozen load does not preserve its accuracy.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedRodLoads {
    nodes: Vec<PreparedSpatialLoads>,
}

impl PreparedRodLoads {
    /// Prepare nodal world forces/couples, extra body-following free couples,
    /// and point forces. The two nodal slices must match the rod's node count;
    /// all point offsets refer to those same original material-node origins.
    /// Gravity and buoyancy remain with the rod's existing environment owner.
    ///
    /// At most `MAX_WORLD_POINT_FORCES` points are admitted over the whole rod.
    /// Bad input or cancellation publishes no partially prepared load.
    pub fn new(
        held_world: &[RodNodalWrench],
        couples_body_nm: &[[f64; 3]],
        points: &[RodPointForce],
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        checkpoint(gate)?;
        let n = held_world.len();
        if !(2..=MAX_MATERIAL_ROD_SEGMENTS + 1).contains(&n)
            || couples_body_nm.len() != n
            || points.len() > MAX_WORLD_POINT_FORCES
        {
            return Err(RodMotionError::Invalid(
                "spatial rod loads need matching bounded nodes and point forces",
            ));
        }
        let mut grouped = vec![Vec::new(); n];
        for point in points {
            checkpoint(gate)?;
            let bucket = grouped.get_mut(point.node).ok_or(RodMotionError::Invalid(
                "spatial point force refers to an absent rod node",
            ))?;
            bucket.push(WorldPointForce {
                // The reducer is independent of inertia. Its reference point
                // is explicitly this material node, not the displaced COM.
                arm_body: vector(point.offset_body_m),
                force_world: vector(point.force_world_n),
            });
        }
        let mut nodes = Vec::with_capacity(n);
        for ((held, couple), points) in held_world.iter().zip(couples_body_nm).zip(grouped) {
            checkpoint(gate)?;
            nodes.push(
                PreparedSpatialLoads::new(
                    Wrench {
                        force_world: vector(held.force_n),
                        torque_body: vector(*couple),
                    },
                    vector(held.torque_nm),
                    &points,
                    || gate.is_requested(),
                )
                .map_err(|error| match error {
                    SpatialLoadError::Cancelled => SpatialRodError::Cancelled.into(),
                    other => RodMotionError::SpatialLoads(other),
                })?,
            );
        }
        checkpoint(gate)?;
        Ok(Self { nodes })
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Resultant force and moment ABOUT the material node, both in world axes.
    /// This is an instantaneous inspection query; freezing the returned moment
    /// over a rotating interval changes the point force's physical action.
    /// Use `RodMotion::step_spatial_candidate` for the prepared motion path.
    pub fn wrench_at(&self, node: usize, frame: So3) -> Result<RodNodalWrench, RodMotionError> {
        self.stage_at(node, frame).map(|(world, _)| world)
    }

    fn stage_at(
        &self,
        node: usize,
        frame: So3,
    ) -> Result<(RodNodalWrench, [f64; 3]), RodMotionError> {
        let loads = self.nodes.get(node).ok_or(RodMotionError::Invalid(
            "spatial load query refers to an absent rod node",
        ))?;
        let q = frame.as_quat();
        let orientation = UnitQuaternion::new(q.w, q.x, q.y, q.z)
            .map_err(|e| RodMotionError::SpatialLoads(e.into()))?;
        let wrench = loads
            .wrench_at(orientation)
            .map_err(|e| RodMotionError::SpatialLoads(e.into()))?;
        let force_n = components(wrench.force_world);
        let torque_nm = components(orientation.rotate_body_to_world(wrench.torque_body));
        finite(&force_n)?;
        finite(&torque_nm)?;
        Ok((
            RodNodalWrench { force_n, torque_nm },
            components(wrench.torque_body),
        ))
    }
}

#[derive(Clone, Copy)]
pub(super) enum StageLoads<'a> {
    World(&'a [RodNodalWrench]),
    Spatial(&'a PreparedRodLoads),
}

impl StageLoads<'_> {
    pub(super) fn len(self) -> usize {
        match self {
            Self::World(loads) => loads.len(),
            Self::Spatial(loads) => loads.node_count(),
        }
    }

    pub(super) fn validate(self) -> Result<(), RodMotionError> {
        if let Self::World(loads) = self {
            for load in loads {
                finite(&load.force_n)?;
                finite(&load.torque_nm)?;
            }
        }
        Ok(())
    }

    pub(super) fn stage_at(
        self,
        node: usize,
        frame: So3,
    ) -> Result<(RodNodalWrench, [f64; 3]), RodMotionError> {
        match self {
            Self::World(loads) => Ok((
                loads[node],
                super::rotate(frame.inverse(), loads[node].torque_nm)?,
            )),
            // Reuse the reducer's body moment in the native kick, avoiding
            // a body -> world -> body round trip at every node and stage.
            Self::Spatial(loads) => loads.stage_at(node, frame),
        }
    }
}

fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

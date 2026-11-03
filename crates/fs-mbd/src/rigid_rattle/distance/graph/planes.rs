//! Finite moving planes in the original coupled native-body graph.

use super::super::super::finite_plane::PreparedPlane;
use super::{DistanceStepError, Geometry, MassProperties, RigidBodyState, Row, RowKind, Vec3};

/// A material point of one body slides in a plane carried by a finite body.
///
/// The first body owns the point and the second carries the plane. Reaction
/// forces act at the actual projected sliding point on the carrier, so its
/// changing angular Jacobian, recoil and angular momentum are retained. Both
/// bodies may belong to other joints in the same simultaneous graph.
///
/// This contributes one smooth bilateral scalar row. It permits tensile
/// reactions and unrestricted in-plane motion, with free relative attitude.
/// Collision discovery, finite rail extent and friction require their own
/// laws; this joint does not infer those from plane geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointPlaneJoint {
    bodies: [usize; 2],
    prepared: PreparedPlane,
}

impl PointPlaneJoint {
    /// Prepare a finite point/carrier pair. The nonzero body-local plane
    /// normal is normalized once and only its anchor's normal offset is kept.
    /// Actual body-slice indices and local tolerances are admitted by the graph.
    pub fn new(
        bodies: [usize; 2],
        point_body_m: Vec3,
        plane_point_body_m: Vec3,
        plane_normal_body: Vec3,
    ) -> Result<Self, DistanceStepError> {
        if bodies[0] == bodies[1] {
            return Err(DistanceStepError::Invalid("point/plane bodies must differ"));
        }
        Ok(Self {
            bodies,
            prepared: PreparedPlane::new(point_body_m, plane_point_body_m, plane_normal_body)?,
        })
    }

    /// Point body, then carrier body, in the original graph body order.
    pub const fn bodies(self) -> [usize; 2] {
        self.bodies
    }

    // A prismatic joint groups two orthogonal plane rows for vector admission.
    // Standalone point/plane joints supply component=0 and width=1.
    pub(super) fn row(self, component: usize, width: usize) -> Row {
        Row {
            bodies: self.bodies,
            points_body_m: [self.prepared.point_body_m, Vec3::ZERO],
            kind: RowKind::Plane(PlaneRow {
                prepared: self.prepared,
                component,
                width,
            }),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct PlaneRow {
    prepared: PreparedPlane,
    pub component: usize,
    pub width: usize,
}

pub(super) fn geometry(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    row: PlaneRow,
    position_tolerance_m: f64,
) -> Result<Geometry, DistanceStepError> {
    row.prepared
        .geometry(states, properties, position_tolerance_m)
}

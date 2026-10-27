//! Native external impulses at points carried by the original material regions.
//!
//! A load is an explicit event at the accepted thermal/mechanical time. It does
//! not invent a force law, integration duration, heat destination or collision
//! response. The caller owns the equal reaction and external work source.

use super::*;
use fs_mbd::{ImpulseReceipt, PointKinematics};

/// A point in an original material region, before source thermal expansion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAttachment {
    pub cell: CellId,
    /// Reference body coordinates [m], expanded about this region's fixed point.
    pub reference_point_m: fs_ga::Vec3,
}

/// Externally computed momentum input in N s, at one declared material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialImpulse {
    pub attachment: MaterialAttachment,
    pub impulse_world: fs_mbd::Vec3,
}

/// Explicit event work and allocation bounds; no physical loss coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialImpulseControls {
    pub maximum_impulses: usize,
    /// Absolute allowance [J] for every native and combined work residual.
    pub work_residual_tolerance_j: f64,
}

/// Native event receipts and a separate external-work account at unchanged time.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialImpulseAudit {
    receipts: Vec<ImpulseReceipt>,
    pub time_s: f64,
    /// Sum of J dot midpoint point velocity, with both signs retained [J].
    pub external_work_j: f64,
    /// Independently evaluated total rigid-body kinetic storage change [J].
    pub kinetic_energy_change_j: f64,
    pub work_energy_residual_j: f64,
    pub maximum_impulse_residual_j: f64,
}
impl MaterialImpulseAudit {
    /// Ordered native receipts for the exact events applied to the candidate.
    pub fn receipts(&self) -> &[ImpulseReceipt] {
        &self.receipts
    }
}

impl HeatedRigidBody<'_> {
    /// Actual expanded point and velocity for an external load/geometry owner.
    /// This uses the same principal-frame map as damping and impulse application.
    pub fn material_point_kinematics(
        &self,
        attachment: MaterialAttachment,
    ) -> Result<PointKinematics, HeatedMotionError> {
        Ok(self.state().point_kinematics(
            self.principal_frame().properties(),
            attachment_arm(self, attachment)?,
        )?)
    }

    /// Apply declared instantaneous loads without rebuilding the heated body.
    ///
    /// Events use slice order at the SAME accepted clock. Each point follows
    /// its original region's expansion, centre of mass and full principal frame.
    /// Pose, thermal cells, boundaries and event cursor remain those of the
    /// accepted state. The returned body continues through `advance_to`.
    /// Physical loss/impact heat is not inferred from signed external work.
    ///
    /// Every input, native event and work balance must pass before publication.
    /// A late refusal or cancellation leaves the original entire history intact.
    pub fn apply_material_impulses(
        &self,
        loads: &[MaterialImpulse],
        controls: MaterialImpulseControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, MaterialImpulseAudit), HeatedMotionError> {
        let mut checkpoint = || {
            if cancelled() {
                Err(HeatedMotionError::from(DepositionError::Cancelled))
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        let tolerance = controls.work_residual_tolerance_j;
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(HeatedMotionError::Invalid(
                "impulse work tolerance must be finite nonnegative joules",
            ));
        }
        if loads.len() > controls.maximum_impulses {
            return Err(HeatedMotionError::Budget("material impulse count"));
        }
        let before = self.diagnostics()?;
        let mut next = self.clone();
        let mut receipts = Vec::with_capacity(loads.len());
        let mut external_work_j = 0.0;
        let mut maximum_impulse_residual_j: f64 = 0.0;
        for load in loads {
            checkpoint()?;
            let arm = attachment_arm(&next, load.attachment)?;
            let receipt = next.state().apply_impulse_at_body_point(
                next.principal_frame().properties(),
                arm,
                load.impulse_world,
            )?;
            let a = receipt.point_kinematics_before.point_velocity_world;
            let b = receipt.point_kinematics_after.point_velocity_world;
            let work = finite(load.impulse_world.dot(a.scale(0.5).add(b.scale(0.5))))?;
            let residual = finite(receipt.work.kinetic_energy_change_j - work)?;
            let native_difference = finite(receipt.work.impulse_work_j - work)?;
            maximum_impulse_residual_j = maximum_impulse_residual_j
                .max(residual.abs())
                .max(native_difference.abs())
                .max(receipt.work.work_energy_residual_j.abs());
            if maximum_impulse_residual_j > tolerance {
                return Err(HeatedMotionError::Budget("material impulse work residual"));
            }
            external_work_j = finite(external_work_j + work)?;
            next.motion.state = receipt.state_after;
            receipts.push(receipt);
        }
        let after = next.diagnostics()?;
        let kinetic_energy_change_j = finite(
            (after.translational_kinetic_energy - before.translational_kinetic_energy)
                + (after.rotational_kinetic_energy - before.rotational_kinetic_energy),
        )?;
        let work_energy_residual_j = finite(kinetic_energy_change_j - external_work_j)?;
        if work_energy_residual_j.abs() > tolerance {
            return Err(HeatedMotionError::Budget(
                "combined material impulse work residual",
            ));
        }
        checkpoint()?;
        Ok((
            next,
            MaterialImpulseAudit {
                receipts,
                time_s: self.thermal.time_s(),
                external_work_j,
                kinetic_energy_change_j,
                work_energy_residual_j,
                maximum_impulse_residual_j,
            },
        ))
    }
}

pub(super) fn attachment_region<'b, 'a>(
    body: &'b HeatedRigidBody<'a>,
    attachment: MaterialAttachment,
) -> Result<&'b ThermalRegion<'a>, HeatedMotionError> {
    if ![
        attachment.reference_point_m.x,
        attachment.reference_point_m.y,
        attachment.reference_point_m.z,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(HeatedMotionError::Invalid("nonfinite material attachment"));
    }
    body.motion
        .expansion
        .regions
        .iter()
        .find(|r| r.cell_id == attachment.cell)
        .ok_or(HeatedMotionError::Invalid(
            "attachment must name an original material/thermal region",
        ))
}

pub(super) fn attachment_arm(
    body: &HeatedRigidBody<'_>,
    attachment: MaterialAttachment,
) -> Result<fs_mbd::Vec3, HeatedMotionError> {
    let binding = attachment_region(body, attachment)?;
    let index = body
        .reference()
        .regions()
        .iter()
        .position(|r| r.key() == binding.key)
        .ok_or(HeatedMotionError::Invalid(
            "missing original attachment region",
        ))?;
    let region = &body.configuration().regions()[index];
    let point = region.fixed_point_m()
        + (attachment.reference_point_m - region.fixed_point_m())
            .scale(1.0 + region.expansion().free_linear_strain());
    Ok(body.principal_frame().arm_for_body_point(point)?)
}

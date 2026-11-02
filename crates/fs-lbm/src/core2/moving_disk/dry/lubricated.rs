//! Implicit normal cylinder/plane films on the physical dry-disk candidate.

use super::*;
use crate::core2::{UnderResolvedCylinderPair2, UnderResolvedLubricationExchange2};

/// A translating prescribed plane, with its unit normal directed into fluid.
///
/// Coordinates describe the start of a lattice tick. This is a momentum/work
/// reservoir, not a finite-inertia wall or a resolved-fluid boundary. The caller
/// supplies the corresponding resolved geometry and advances this descriptor
/// only after an accepted tick. Tangential drag, contact and plane rotation are
/// outside the original normal-cylinder film model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LubricationPlane2 {
    point: [f64; 2],
    normal: [f64; 2],
    velocity: [f64; 2],
    cutoff_gap: f64,
    minimum_gap: f64,
}

impl LubricationPlane2 {
    /// Normalize a nonzero normal and admit the original film geometry/speed.
    pub fn new(
        point: [f64; 2],
        normal: [f64; 2],
        velocity: [f64; 2],
        cutoff_gap: f64,
        minimum_gap: f64,
    ) -> Result<Self, DryDiskError2> {
        let scale = normal[0].abs().max(normal[1].abs());
        let speed_sq = velocity[0].mul_add(velocity[0], velocity[1] * velocity[1]);
        if !point
            .into_iter()
            .chain(normal)
            .chain(velocity)
            .all(f64::is_finite)
            || !(scale > 0.0 && speed_sq < 0.03)
            || !(cutoff_gap.is_finite()
                && cutoff_gap > 0.0
                && minimum_gap.is_finite()
                && minimum_gap > 0.0
                && minimum_gap < cutoff_gap)
        {
            return Err(DryDiskError2::Invalid("invalid dry-disk lubrication plane"));
        }
        let normal = normal.map(|x| x / scale);
        let length = normal[0].mul_add(normal[0], normal[1] * normal[1]).sqrt();
        Ok(Self {
            point,
            normal: normal.map(|x| x / length),
            velocity,
            cutoff_gap,
            minimum_gap,
        })
    }

    /// One point on the plane at the start of the next supplied tick.
    pub const fn point(self) -> [f64; 2] {
        self.point
    }
    /// Unit normal into the fluid half-space.
    pub const fn normal(self) -> [f64; 2] {
        self.normal
    }
    /// Prescribed constant translational velocity in lattice units.
    pub const fn velocity(self) -> [f64; 2] {
        self.velocity
    }
    /// Distance above which the missing resistance is zero.
    pub const fn cutoff_gap(self) -> f64 {
        self.cutoff_gap
    }
    /// Positive floor for the original normal-resistance singularity.
    pub const fn minimum_gap(self) -> f64 {
        self.minimum_gap
    }

    /// Advance only this descriptor through a nonnegative lattice duration.
    /// No body, fluid state or accepted clock is modified.
    pub fn translated(self, duration: f64) -> Result<Self, DryDiskError2> {
        if !(duration.is_finite() && duration >= 0.0) {
            return Err(DryDiskError2::Invalid("invalid plane translation duration"));
        }
        let point = core::array::from_fn(|a| duration.mul_add(self.velocity[a], self.point[a]));
        if !point.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid("plane translation overflowed"));
        }
        Ok(Self { point, ..self })
    }
}

/// One original normal-film exchange evaluated at the common solved midpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct PlaneLubricationExchange2 {
    /// The disk is the first solid; the plane reservoir is the second. Both
    /// angular impulses use the caller's fixed origin and are exactly opposite.
    pub native: UnderResolvedLubricationExchange2,
    /// Impulse work delivered to the dry disk at its solved midpoint velocity.
    pub body_work: f64,
    /// Impulse work received by the prescribed plane reservoir.
    pub plane_work: f64,
}

/// One physical dry-disk/fluid tick with all normal plane films solved together.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct DryLubricatedDiskExchange2 {
    /// Native fluid and original dry inertia. Its interior correction excludes
    /// these films; its residuals include their separately reported transfer.
    pub motion: DryDiskExchange2,
    /// Plane exchanges in the same order as the supplied descriptors.
    pub planes: Vec<PlaneLubricationExchange2>,
    /// Sum of film impulses received by the disk, about the requested origin.
    pub body_impulse: [f64; 3],
    /// Film work received by the disk, separate from authored/source work.
    pub body_work: f64,
    /// Film work received by the prescribed plane reservoirs.
    pub plane_work: f64,
    /// Nonnegative missing-film loss: `-(body_work + plane_work)` up to roundoff.
    /// No thermal destination is invented or advanced by this owner.
    pub dissipated_energy: f64,
}

impl DryDisk2 {
    /// Advance physical dry inertia and all prescribed normal plane films in
    /// one coupled fluid candidate, using the existing cylinder resistance.
    ///
    /// The three unknowns are the total supplemental impulse/couple. Their
    /// residual includes the actual fictitious-fluid change and every film
    /// evaluated at the SAME native midpoint velocity. The native dry mass and
    /// polar inertia remain positive and unchanged, including light solids.
    /// Film force acts along each center-to-plane normal, so it adds no spin
    /// couple. Its fixed-origin torque is the center lever arm of that force.
    ///
    /// `dynamic_viscosity` is the positive fluid dynamic viscosity in lattice
    /// units; all film loads/work are per unit out-of-plane depth per tick.
    /// At most 64 distinct planes are admitted. Start, midpoint and accepted
    /// endpoint gaps must remain strictly positive through both free drifts.
    /// Resistance uses the midpoint gap and the supplied floor/cutoff. Neither
    /// contacts nor a destination for the dissipated heat are supplied.
    ///
    /// Planes and accepted disk/fluid state stay unchanged on refusal. Advance
    /// supplied plane coordinates with `translated(1.0)` only after success.
    /// Empty or entirely inactive films preserve the ordinary dry-step result.
    /// The original native collision and cylinder-law assertion APIs still
    /// apply to unrepresentable trial states/resistance before publication.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubrication_planes(
        &mut self,
        grid: &mut Grid,
        external_impulse: [f64; 2],
        external_couple: f64,
        moment_origin: [f64; 2],
        dynamic_viscosity: f64,
        planes: &[LubricationPlane2],
        config: DryDiskCouplingConfig2,
        cancel: &CancelGate,
    ) -> Result<DryLubricatedDiskExchange2, DryDiskError2> {
        checkpoint(cancel)?;
        config.admit()?;
        if !(dynamic_viscosity.is_finite() && dynamic_viscosity > 0.0)
            || !moment_origin.into_iter().all(f64::is_finite)
            || planes.len() > 64
        {
            return Err(DryDiskError2::Invalid(
                "invalid dry-disk film viscosity, origin or plane count",
            ));
        }
        let prepared = PreparedPlanes2::new(self, dynamic_viscosity, planes, cancel)?;
        self.step_inner(
            grid,
            external_impulse,
            external_couple,
            moment_origin,
            config,
            cancel,
            Some(&prepared),
        )
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct FilmTotals2 {
    pub(super) impulse: [f64; 3],
    pub(super) body_work: f64,
    pub(super) plane_work: f64,
    pub(super) dissipation: f64,
    pub(super) active: usize,
}
struct PreparedPlane2 {
    pair: UnderResolvedCylinderPair2,
    endpoint: LubricationPlane2,
}
pub(super) struct PreparedPlanes2(Vec<PreparedPlane2>);
impl PreparedPlanes2 {
    fn new(
        disk: &DryDisk2,
        viscosity: f64,
        planes: &[LubricationPlane2],
        cancel: &CancelGate,
    ) -> Result<Self, DryDiskError2> {
        let body = disk.body();
        let center = core::array::from_fn(|a| body.center[a] + 0.5 * body.velocity[a]);
        let radius = disk.radius();
        let mut prepared = Vec::with_capacity(planes.len());
        for (i, &plane) in planes.iter().enumerate() {
            checkpoint(cancel)?;
            if planes[..i].iter().any(|other| {
                other.normal == plane.normal
                    && dot(subtract(other.point, plane.point), plane.normal) == 0.0
            }) {
                return Err(DryDiskError2::Invalid(
                    "duplicate dry-disk lubrication plane",
                ));
            }
            positive_gap(body.center, plane, radius)?;
            let midpoint = plane.translated(0.5)?;
            let gap = positive_gap(center, midpoint, radius)?;
            let pair = UnderResolvedCylinderPair2::new(
                gap,
                plane.cutoff_gap,
                plane.minimum_gap,
                radius,
                viscosity,
                plane.normal.map(|x| -x),
                center,
                [[0.0; 2], plane.velocity],
            );
            prepared.push(PreparedPlane2 {
                pair,
                endpoint: plane.translated(1.0)?,
            });
        }
        Ok(Self(prepared))
    }

    pub(super) fn evaluate(
        &self,
        velocity: [f64; 2],
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<FilmTotals2, DryDiskError2> {
        let mut sums = [Dd::ZERO; 6];
        let mut active = 0;
        for plane in &self.0 {
            checkpoint(cancel)?;
            let receipt = plane.exchange(velocity, origin)?;
            let native = receipt.native;
            if native.resistance_increment != 0.0 {
                active += 1;
                let values = [
                    native.first_solid_impulse[0],
                    native.first_solid_impulse[1],
                    native.first_solid_angular_impulse,
                    receipt.body_work,
                    receipt.plane_work,
                    -native.solid_pair_work,
                ];
                for (sum, value) in sums.iter_mut().zip(values) {
                    *sum = *sum + Dd::from_f64(value);
                }
            }
        }
        let values = sums.map(Dd::to_f64);
        if !values.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid("dry-disk film sum overflowed"));
        }
        Ok(FilmTotals2 {
            impulse: [values[0], values[1], values[2]],
            body_work: values[3],
            plane_work: values[4],
            dissipation: values[5],
            active,
        })
    }

    pub(super) fn admit_endpoint(
        &self,
        center: [f64; 2],
        radius: f64,
        cancel: &CancelGate,
    ) -> Result<(), DryDiskError2> {
        for plane in &self.0 {
            checkpoint(cancel)?;
            positive_gap(center, plane.endpoint, radius)?;
        }
        Ok(())
    }

    pub(super) fn receipts(
        &self,
        velocity: [f64; 2],
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<Vec<PlaneLubricationExchange2>, DryDiskError2> {
        self.0
            .iter()
            .map(|plane| {
                checkpoint(cancel)?;
                plane.exchange(velocity, origin)
            })
            .collect()
    }
}
impl PreparedPlane2 {
    fn exchange(
        &self,
        velocity: [f64; 2],
        origin: [f64; 2],
    ) -> Result<PlaneLubricationExchange2, DryDiskError2> {
        let mut pair = self.pair;
        pair.surface_velocities[0] = velocity;
        let native = pair.correction(origin);
        let body_work = dot(native.first_solid_impulse, velocity);
        let plane_work = dot(native.second_solid_impulse, pair.surface_velocities[1]);
        let scale = body_work
            .abs()
            .max(plane_work.abs())
            .max(native.solid_pair_work.abs());
        if !body_work.is_finite()
            || !plane_work.is_finite()
            || (body_work + plane_work - native.solid_pair_work).abs()
                > 4096.0 * f64::EPSILON * scale
        {
            return Err(DryDiskError2::Invalid("dry-disk film work is unresolved"));
        }
        Ok(PlaneLubricationExchange2 {
            native,
            body_work,
            plane_work,
        })
    }
}
fn positive_gap(
    center: [f64; 2],
    plane: LubricationPlane2,
    radius: f64,
) -> Result<f64, DryDiskError2> {
    let gap = dot(subtract(center, plane.point), plane.normal) - radius;
    if !(gap.is_finite() && gap > 64.0 * f64::EPSILON * radius) {
        Err(DryDiskError2::Invalid(
            "dry-disk plane contact requires a separate contact owner",
        ))
    } else {
        Ok(gap)
    }
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[0], a[1] * b[1])
}

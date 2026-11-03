//! Implicit cylinder/plane squeeze, sliding and spin films on the dry disk.

use super::*;
use crate::core2::{UnderResolvedCylinderPair2, UnderResolvedLubricationExchange2};

/// A translating prescribed plane, with its unit normal directed into fluid.
///
/// Coordinates describe the start of a lattice tick. This is a momentum/work
/// reservoir, not a finite-inertia wall or a resolved-fluid boundary. The caller
/// supplies the corresponding resolved geometry and advances this descriptor
/// only after an accepted tick. Normal resistance is always available;
/// [Self::with_tangential] opts into the cylinder-specific sliding/spin law.
/// Contact, plane rotation and finite wall recoil remain separate owners.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LubricationPlane2 {
    point: [f64; 2],
    normal: [f64; 2],
    velocity: [f64; 2],
    cutoff_gap: f64,
    minimum_gap: f64,
    tangential: Option<(f64, f64)>,
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
            tangential: None,
        })
    }

    /// Add the leading small-gap sliding AND spinning cylinder resistance.
    ///
    /// Per unit depth, the missing resistance is
    /// zeta = 2*pi*mu*sqrt(2*R)*(1/sqrt(h_eff) - 1/sqrt(cutoff)).
    /// The disk receives tangent force -zeta*(V-Vwall).t*t and center couple
    /// -zeta*R^2*omega. Translation and spin DECOUPLE for an infinite cylinder:
    /// a dashpot on contact-point slip, or the sphere's mixed resistance,
    /// would incorrectly remove rolling dissipation or introduce cross terms.
    ///
    /// These independent tangent/spin scales are lattice lengths. The positive
    /// floor gives h_eff = max(gap, minimum_gap); the correction vanishes at and
    /// above cutoff. At solve time cutoff must be <= 0.02*R, a small-gap policy
    /// for this leading-order Stokes model, not a finite-gap accuracy bound.
    /// Normal resistance retains the constructor's own cutoff and floor.
    ///
    /// Jeffrey & Onishi, QJMAM 34 (1981), 129-137,
    /// <https://doi.org/10.1093/qjmam/34.2.129>; the leading coefficients are
    /// reproduced in Reboux, Richardson & Jensen, equation (2.1),
    /// <https://personalpages.manchester.ac.uk/staff/oliver.jensen/OpenArchive/RebouxRichardsonJensen.pdf>.
    /// No contact, cavitation, finite-span end effects or heat state is added.
    pub fn with_tangential(
        mut self,
        cutoff_gap: f64,
        minimum_gap: f64,
    ) -> Result<Self, DryDiskError2> {
        if !(cutoff_gap.is_finite() && cutoff_gap > 0.0)
            || !(minimum_gap.is_finite() && minimum_gap > 0.0)
            || minimum_gap > cutoff_gap
        {
            return Err(DryDiskError2::Invalid(
                "invalid cylinder tangential film scales",
            ));
        }
        self.tangential = Some((cutoff_gap, minimum_gap));
        Ok(self)
    }

    /// Optional sliding/spin (cutoff, positive floor), separate from squeeze.
    pub const fn tangential_scales(self) -> Option<(f64, f64)> {
        self.tangential
    }

    fn tangential_resistance(
        self,
        radius: f64,
        gap: f64,
        viscosity: f64,
    ) -> Result<f64, DryDiskError2> {
        let Some((cutoff, floor)) = self.tangential else {
            return Ok(0.0);
        };
        if cutoff > 0.02 * radius {
            return Err(DryDiskError2::Invalid(
                "tangential lubrication cutoff exceeds 0.02 cylinder radius",
            ));
        }
        let effective = gap.max(floor);
        if effective >= cutoff {
            return Ok(0.0);
        }
        // Rationalize 1-sqrt(h/cutoff) so close-to-cutoff films do not lose
        // their small positive resistance by subtracting two square roots.
        let deficit = ((cutoff - effective) / cutoff) / (1.0 + (effective / cutoff).sqrt());
        let resistance =
            (2.0 * core::f64::consts::PI * core::f64::consts::SQRT_2) * viscosity * radius.sqrt()
                / effective.sqrt()
                * deficit;
        if !(resistance.is_finite() && resistance > 0.0) {
            return Err(DryDiskError2::Invalid(
                "cylinder tangential resistance is unrepresentable",
            ));
        }
        Ok(resistance)
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

/// One complete film exchange evaluated at the common solved midpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct PlaneLubricationExchange2 {
    /// Original NORMAL-only cylinder law. The disk is the first solid and the
    /// plane the second; its moments use the same fixed origin as the totals.
    pub native: UnderResolvedLubricationExchange2,
    /// Complete [linear_x, linear_y, angular] impulse received by the disk.
    /// The angular component includes its center couple and force lever arm
    /// about the requested fixed origin.
    pub body_impulse: [f64; 3],
    /// Exactly opposite complete impulse received by the prescribed plane,
    /// about that SAME fixed origin.
    pub plane_impulse: [f64; 3],
    /// Intrinsic disk-center angular impulse; normal and pure sliding films
    /// add zero, while the spinning cylinder receives -zeta*R^2*omega.
    pub body_center_couple: f64,
    /// Projection of the collision center onto the midpoint plane.
    pub plane_footpoint: [f64; 2],
    /// Free reaction couple at that footpoint. Together with the plane force
    /// lever arm it equals plane_impulse[2], retaining total angular impulse
    /// for the distributed tangential traction, including pure sliding.
    pub plane_couple: f64,
    /// Sliding resistance increment; spin resistance is R^2 times this value.
    /// Zero when disabled, at/above cutoff, or floor equals cutoff.
    pub tangential_resistance_increment: f64,
    /// Impulse work delivered to the dry disk at its solved midpoint velocity.
    pub body_work: f64,
    /// Impulse work received by the prescribed plane reservoir.
    pub plane_work: f64,
    /// Nonnegative normal plus sliding/spinning loss per lattice tick.
    pub dissipated_energy: f64,
}

/// One physical dry-disk/fluid tick with all selected plane films solved together.
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
    /// Advance physical dry inertia and all prescribed plane films in one
    /// coupled fluid candidate. Sliding/spin films require with_tangential.
    ///
    /// The three unknowns are the total supplemental impulse/couple. Their
    /// residual includes the actual fictitious-fluid change and every film
    /// evaluated at the SAME native midpoint translation and spin. The dry mass and
    /// polar inertia remain positive and unchanged, including light solids.
    /// Normal films add no center couple. The optional cylinder spin resistance
    /// enters the same implicit angular solve; each plane receives the full
    /// opposing force and moment, including its distributed-traction couple.
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
    pub(super) center_couple: f64,
    pub(super) body_work: f64,
    pub(super) plane_work: f64,
    pub(super) dissipation: f64,
    pub(super) active: usize,
}
struct PreparedPlane2 {
    pair: UnderResolvedCylinderPair2,
    tangent: [f64; 2],
    tangential_resistance: f64,
    footpoint: [f64; 2],
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
            prepared.push(PreparedPlane2::new(plane, center, radius, viscosity)?);
        }
        Ok(Self(prepared))
    }

    pub(super) fn evaluate(
        &self,
        velocity: [f64; 2],
        angular_velocity: f64,
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<FilmTotals2, DryDiskError2> {
        let mut sums = [Dd::ZERO; 7];
        let mut active = 0;
        for plane in &self.0 {
            checkpoint(cancel)?;
            let receipt = plane.exchange(velocity, angular_velocity, origin)?;
            let native = receipt.native;
            if native.resistance_increment != 0.0 || plane.tangential_resistance != 0.0 {
                active += 1;
                let values = [
                    receipt.body_impulse[0],
                    receipt.body_impulse[1],
                    receipt.body_impulse[2],
                    receipt.body_work,
                    receipt.plane_work,
                    receipt.dissipated_energy,
                    receipt.body_center_couple,
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
            center_couple: values[6],
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
        angular_velocity: f64,
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<Vec<PlaneLubricationExchange2>, DryDiskError2> {
        self.0
            .iter()
            .map(|plane| {
                checkpoint(cancel)?;
                plane.exchange(velocity, angular_velocity, origin)
            })
            .collect()
    }
}
impl PreparedPlane2 {
    fn new(
        plane: LubricationPlane2,
        center: [f64; 2],
        radius: f64,
        viscosity: f64,
    ) -> Result<Self, DryDiskError2> {
        let gap = positive_gap(center, plane.translated(0.5)?, radius)?;
        let tangential_resistance = plane.tangential_resistance(radius, gap, viscosity)?;
        let footpoint = core::array::from_fn(|a| center[a] - (radius + gap) * plane.normal[a]);
        if !footpoint.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid("film footpoint is unrepresentable"));
        }
        Ok(Self {
            pair: UnderResolvedCylinderPair2::new(
                gap,
                plane.cutoff_gap,
                plane.minimum_gap,
                radius,
                viscosity,
                plane.normal.map(|x| -x),
                center,
                [[0.0; 2], plane.velocity],
            ),
            tangent: [-plane.normal[1], plane.normal[0]],
            tangential_resistance,
            footpoint,
            endpoint: plane.translated(1.0)?,
        })
    }

    fn exchange(
        &self,
        velocity: [f64; 2],
        angular_velocity: f64,
        origin: [f64; 2],
    ) -> Result<PlaneLubricationExchange2, DryDiskError2> {
        if !velocity
            .into_iter()
            .chain([angular_velocity])
            .all(f64::is_finite)
        {
            return Err(DryDiskError2::Invalid("nonfinite film midpoint velocity"));
        }
        let mut pair = self.pair;
        pair.surface_velocities[0] = velocity;
        let native = pair.correction(origin);
        let mut body_work = dot(native.first_solid_impulse, velocity);
        let mut plane_work = dot(native.second_solid_impulse, pair.surface_velocities[1]);
        let mut body_impulse = [
            native.first_solid_impulse[0],
            native.first_solid_impulse[1],
            native.first_solid_angular_impulse,
        ];
        let mut body_center_couple = 0.0;
        let mut plane_couple = 0.0;
        let mut dissipated_energy = -native.solid_pair_work;
        if self.tangential_resistance != 0.0 {
            let k = self.tangential_resistance;
            let slip = dot(subtract(velocity, pair.surface_velocities[1]), self.tangent);
            let spin_speed = pair.reduced_radius * angular_velocity;
            let tangent_impulse = self.tangent.map(|t| -k * slip * t);
            body_center_couple = if spin_speed == 0.0 {
                0.0
            } else {
                -(k * spin_speed) * pair.reduced_radius
            };
            body_impulse[0] += tangent_impulse[0];
            body_impulse[1] += tangent_impulse[1];
            body_impulse[2] += cross(subtract(pair.line_of_action_point, origin), tangent_impulse)
                + body_center_couple;
            plane_couple = -body_center_couple
                + cross(
                    subtract(self.footpoint, pair.line_of_action_point),
                    tangent_impulse,
                );
            body_work += dot(tangent_impulse, velocity) + body_center_couple * angular_velocity;
            plane_work -= dot(tangent_impulse, pair.surface_velocities[1]);
            dissipated_energy += k * (slip * slip + spin_speed * spin_speed);
        }
        let plane_impulse = body_impulse.map(|j| if j == 0.0 { 0.0 } else { -j });
        let scale = body_work.abs().max(plane_work.abs()).max(dissipated_energy);
        if !body_impulse
            .into_iter()
            .chain([
                body_center_couple,
                plane_couple,
                body_work,
                plane_work,
                dissipated_energy,
            ])
            .all(f64::is_finite)
            || dissipated_energy < 0.0
            || (body_work + plane_work + dissipated_energy).abs() > 4096.0 * f64::EPSILON * scale
        {
            return Err(DryDiskError2::Invalid("dry-disk film work is unresolved"));
        }
        Ok(PlaneLubricationExchange2 {
            native,
            body_impulse,
            plane_impulse,
            body_center_couple,
            plane_footpoint: self.footpoint,
            plane_couple,
            tangential_resistance_increment: self.tangential_resistance,
            body_work,
            plane_work,
            dissipated_energy,
        })
    }
}

#[cfg(test)]
mod tests;
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

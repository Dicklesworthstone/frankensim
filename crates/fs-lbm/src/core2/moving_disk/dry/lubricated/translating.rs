//! Finite normal and tangential plane recoil inside the native dry-disk solve.

use super::slider::CondensedResistance2;
use super::*;

/// A finite plane carrier free to translate in both directions, with fixed orientation.
///
/// The descriptor carries its accepted point and actual velocity. Its normal
/// stays fixed, while the center of mass and plane point undergo the same two
/// free half-drifts. Mass is supplied per out-of-plane depth; an infinite local
/// plane does not define a material volume or mass by itself.
/// The orientation guide supplies a zero-work couple only, never linear force.
/// This film-only owner does not create a resolved moving-fluid boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteLubricationPlane2 {
    plane: LubricationPlane2,
    mass: f64,
    center: [f64; 2],
}

impl FiniteLubricationPlane2 {
    /// Bind original plane/film geometry, positive mass and an independent
    /// center of mass. The plane descriptor's velocity is the carrier velocity.
    pub fn new(
        plane: LubricationPlane2,
        mass: f64,
        center: [f64; 2],
    ) -> Result<Self, DryDiskError2> {
        if !(mass.is_finite() && mass > 0.0)
            || !center.into_iter().chain(plane.point).all(f64::is_finite)
            || !admitted_velocity(plane.velocity)
        {
            return Err(DryDiskError2::Invalid(
                "finite plane needs positive mass and finite low-Mach translation",
            ));
        }
        Ok(Self {
            plane,
            mass,
            center,
        })
    }

    /// Authoritative accepted plane point, fixed normal and actual velocity.
    pub const fn plane(self) -> LubricationPlane2 {
        self.plane
    }
    /// Positive physical mass per unit out-of-plane depth.
    pub const fn mass(self) -> f64 {
        self.mass
    }
    /// Accepted center of mass in lattice coordinates.
    pub const fn center(self) -> [f64; 2] {
        self.center
    }
    /// Actual accepted world velocity, including normal recoil.
    pub const fn velocity(self) -> [f64; 2] {
        self.plane.velocity
    }
}

/// One finite translating carrier's actual trajectory and storage balance.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct TranslatingPlaneExchange2 {
    /// Complete endpoint assigned to the supplied carrier after success.
    pub state_after: FiniteLubricationPlane2,
    /// Actual center of mass at the common disk/fluid/film impulse.
    pub collision_center: [f64; 2],
    /// Solved common midpoint translation, in world lattice coordinates.
    pub velocity_midpoint: [f64; 2],
    /// Actual [Px, Py, Lz] increment about the requested fixed origin.
    pub momentum_change: [f64; 3],
    /// Authored center impulse and its fixed-origin moment.
    pub external_impulse: [f64; 3],
    /// Pure couple DELIVERED by the orientation guide. It has zero work;
    /// no angular velocity, rotary inertia or linear support reaction is inferred.
    pub guide_couple: f64,
    /// Work of the authored center impulse at the solved midpoint.
    pub external_work: f64,
    /// Actual change of carrier kinetic energy.
    pub kinetic_energy_change: f64,
    /// Actual momentum change minus film, authored and guide impulse.
    pub momentum_residual: [f64; 3],
    /// Kinetic increment minus authored work minus received film work.
    pub work_residual: f64,
}

/// Native disk/fluid motion and simultaneous finite translating plane carriers.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct DryTranslatingPlanesExchange2 {
    /// Original native physical-dry result and per-plane film reactions.
    /// Plane work enters finite kinetic storage, not a prescribed reservoir.
    pub motion: DryLubricatedDiskExchange2,
    /// Finite carrier responses in exactly the supplied plane order.
    pub planes: Vec<TranslatingPlaneExchange2>,
    /// Exterior fluid plus dry disk plus all finite carriers, after subtracting
    /// authored inputs, guide couples and the original fluid/boundary sources.
    pub combined_momentum_residual: [f64; 3],
}

impl DryDisk2 {
    /// Advance native fluid, a dry disk and up to 64 translating finite planes.
    ///
    /// Each plane receives its own authored WORLD center impulse, separate
    /// from film loads. Original normal and tangent resistance are each
    /// condensed with the same positive mass using k/(1+k/(2m)); all planes
    /// participate in every existing three-unknown dry-fluid candidate.
    /// The original cylinder spin resistance is retained. The guide opposes
    /// the complete center couple and supplies no linear impulse or work.
    ///
    /// The first half-drift uses the accepted initial velocities to determine
    /// the actual midpoint gap. The second uses both solved endpoint velocities.
    /// Start, midpoint and endpoint gaps must be positive: since orientation
    /// is fixed, each half-drift's relative gap is affine and cannot cross
    /// between admitted endpoints. A held-wall endpoint is never substituted.
    ///
    /// Every carrier state and balance is prepared before the sole native
    /// disk/fluid publication; subsequent carrier assignments are infallible.
    /// Refusal and cancellation preserve every accepted state. No face extent,
    /// contact, resolved-wall motion or film-heat destination is supplied.
    /// The original implicit midpoint method is not L-stable; stiff unresolved
    /// resistance may reverse relative velocity and still needs time resolution.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_translating_lubrication_planes(
        &mut self,
        grid: &mut Grid,
        external_impulse: [f64; 2],
        external_couple: f64,
        moment_origin: [f64; 2],
        dynamic_viscosity: f64,
        planes: &mut [FiniteLubricationPlane2],
        plane_external_impulses: &[[f64; 2]],
        config: DryDiskCouplingConfig2,
        cancel: &CancelGate,
    ) -> Result<DryTranslatingPlanesExchange2, DryDiskError2> {
        checkpoint(cancel)?;
        config.admit()?;
        if !(dynamic_viscosity.is_finite() && dynamic_viscosity > 0.0)
            || !moment_origin.into_iter().all(f64::is_finite)
            || planes.len() > 64
            || plane_external_impulses.len() != planes.len()
        {
            return Err(DryDiskError2::Invalid(
                "invalid finite-plane viscosity, origin, count or load count",
            ));
        }
        let descriptors: Vec<_> = planes.iter().map(|plane| plane.plane).collect();
        let mut prepared = PreparedPlanes2::new(self, dynamic_viscosity, &descriptors, cancel)?;
        for ((film, &plane), &authored) in prepared
            .planes
            .iter_mut()
            .zip(planes.iter())
            .zip(plane_external_impulses)
        {
            checkpoint(cancel)?;
            film.translation = Some(PreparedTranslation2::new(
                plane,
                authored,
                film.pair,
                film.tangential_resistance,
            )?);
        }
        let receipt = self.step_inner(
            grid,
            external_impulse,
            external_couple,
            moment_origin,
            config,
            cancel,
            Some(&prepared),
            |receipt| prepared.finish_translations(receipt, moment_origin, cancel),
        )?;
        for (plane, response) in planes.iter_mut().zip(&receipt.planes) {
            *plane = response.state_after;
        }
        Ok(receipt)
    }
}

pub(super) struct PreparedTranslation2 {
    initial: FiniteLubricationPlane2,
    center: [f64; 2],
    midpoint_plane: LubricationPlane2,
    disk_center: [f64; 2],
    external: [f64; 2],
    free_midpoint: [f64; 2],
    normal_reference: UnderResolvedLubricationExchange2,
    normal: CondensedResistance2,
    tangent: CondensedResistance2,
}

#[derive(Clone, Copy)]
pub(super) struct TranslationResponse2 {
    pub(super) velocity_midpoint: [f64; 2],
    pub(super) normal_impulse: f64,
    pub(super) normal_rate: f64,
    pub(super) tangent_impulse: f64,
    pub(super) tangent_rate: f64,
    pub(super) tangent_midpoint: f64,
}

impl PreparedTranslation2 {
    fn new(
        initial: FiniteLubricationPlane2,
        external: [f64; 2],
        mut pair: UnderResolvedCylinderPair2,
        tangent_resistance: f64,
    ) -> Result<Self, DryDiskError2> {
        // Obtain the original normal coefficient/floor/cutoff from its ONE
        // native owner. No duplicate normal resistance law is introduced.
        pair.surface_velocities = [[0.0; 2]; 2];
        let normal_reference = pair.correction([0.0; 2]);
        let normal =
            CondensedResistance2::new(initial.mass, normal_reference.resistance_increment)?;
        let tangent = CondensedResistance2::new(initial.mass, tangent_resistance)?;
        let free_midpoint =
            core::array::from_fn(|a| initial.velocity()[a] + (0.5 * external[a]) / initial.mass);
        let center = core::array::from_fn(|a| initial.center[a] + 0.5 * initial.velocity()[a]);
        if !external
            .into_iter()
            .chain(free_midpoint)
            .chain(center)
            .all(f64::is_finite)
        {
            return Err(DryDiskError2::Invalid("finite-plane candidate overflowed"));
        }
        Ok(Self {
            initial,
            center,
            midpoint_plane: initial.plane.translated(0.5)?,
            disk_center: pair.line_of_action_point,
            external,
            free_midpoint,
            normal_reference,
            normal,
            tangent,
        })
    }

    pub(super) fn response(&self, velocity: [f64; 2]) -> TranslationResponse2 {
        let n = self.initial.plane.normal;
        let t = [-n[1], n[0]];
        let free_n = dot(self.free_midpoint, n);
        let free_t = dot(self.free_midpoint, t);
        let disk_n = dot(velocity, n);
        let disk_t = dot(velocity, t);
        let delta_n = disk_n - free_n;
        let delta_t = disk_t - free_t;
        let normal_midpoint = self.normal.midpoint(free_n, disk_n);
        let tangent_midpoint = self.tangent.midpoint(free_t, disk_t);
        // Preserve the exact free/common-motion path instead of projecting
        // it into a numerically rounded orthogonal basis and back.
        let velocity_midpoint = if (self.normal.effective == 0.0 && self.tangent.effective == 0.0)
            || (delta_n == 0.0 && delta_t == 0.0)
        {
            self.free_midpoint
        } else {
            core::array::from_fn(|a| normal_midpoint * n[a] + tangent_midpoint * t[a])
        };
        TranslationResponse2 {
            velocity_midpoint,
            normal_impulse: -self.normal.effective * delta_n,
            normal_rate: self.normal.retained * delta_n,
            tangent_impulse: -self.tangent.effective * delta_t,
            tangent_rate: self.tangent.retained * delta_t,
            tangent_midpoint,
        }
    }

    pub(super) fn normal_exchange(
        &self,
        response: TranslationResponse2,
        origin: [f64; 2],
    ) -> Result<UnderResolvedLubricationExchange2, DryDiskError2> {
        normal_exchange(
            self.initial.plane.normal,
            self.disk_center,
            self.normal_reference,
            response,
            origin,
        )
    }

    pub(super) fn endpoint(
        &self,
        velocity: [f64; 2],
    ) -> Result<(FiniteLubricationPlane2, [f64; 2]), DryDiskError2> {
        let midpoint = self.response(velocity).velocity_midpoint;
        let after = core::array::from_fn(|a| 2.0 * midpoint[a] - self.initial.velocity()[a]);
        if !admitted_velocity(midpoint) {
            return Err(DryDiskError2::Invalid(
                "finite-plane midpoint exceeds low-Mach limit",
            ));
        }
        let center = core::array::from_fn(|a| self.center[a] + 0.5 * after[a]);
        let point = core::array::from_fn(|a| self.midpoint_plane.point[a] + 0.5 * after[a]);
        let plane = LubricationPlane2 {
            point,
            velocity: after,
            ..self.initial.plane
        };
        Ok((
            FiniteLubricationPlane2::new(plane, self.initial.mass, center)?,
            midpoint,
        ))
    }

    fn finish(
        &self,
        film: PlaneLubricationExchange2,
        disk_velocity: [f64; 2],
        origin: [f64; 2],
    ) -> Result<TranslatingPlaneExchange2, DryDiskError2> {
        let initial = self.initial;
        let (after, midpoint) = self.endpoint(disk_velocity)?;
        let delta_v = subtract(after.velocity(), initial.velocity());
        let delta = delta_v.map(|v| initial.mass * v);
        let old_momentum = initial.velocity().map(|v| initial.mass * v);
        let momentum_change = [
            delta[0],
            delta[1],
            cross(subtract(after.center, origin), delta)
                + cross(subtract(after.center, initial.center), old_momentum),
        ];
        let arm = subtract(self.center, origin);
        let external_impulse = [
            self.external[0],
            self.external[1],
            cross(arm, self.external),
        ];
        let plane_force = [film.plane_impulse[0], film.plane_impulse[1]];
        let guide_couple =
            -film.plane_couple - cross(subtract(film.plane_footpoint, self.center), plane_force);
        let momentum_residual = core::array::from_fn(|a| {
            momentum_change[a]
                - film.plane_impulse[a]
                - external_impulse[a]
                - if a == 2 { guide_couple } else { 0.0 }
        });
        let external_work = dot(self.external, midpoint);
        let sum_v = core::array::from_fn(|a| after.velocity()[a] + initial.velocity()[a]);
        let kinetic_energy_change = 0.5 * initial.mass * dot(delta_v, sum_v);
        let work_residual = kinetic_energy_change - external_work - film.plane_work;
        let scale = norm(momentum_change)
            .max(norm(film.plane_impulse))
            .max(norm(external_impulse))
            .max(guide_couple.abs())
            .max(norm([
                old_momentum[0],
                old_momentum[1],
                cross(arm, old_momentum),
            ]));
        let energy_scale = kinetic_energy_change
            .abs()
            .max(external_work.abs())
            .max(film.plane_work.abs())
            .max(initial.mass * dot(initial.velocity(), initial.velocity()))
            .max(initial.mass * dot(after.velocity(), after.velocity()));
        if !momentum_change
            .into_iter()
            .chain(external_impulse)
            .chain(momentum_residual)
            .chain([
                guide_couple,
                external_work,
                kinetic_energy_change,
                work_residual,
                scale,
                energy_scale,
            ])
            .all(f64::is_finite)
            || norm(momentum_residual) > 16384.0 * f64::EPSILON * scale
            || work_residual.abs() > 16384.0 * f64::EPSILON * energy_scale
        {
            return Err(DryDiskError2::Invalid(
                "finite-plane momentum/work balance is unresolved",
            ));
        }
        Ok(TranslatingPlaneExchange2 {
            state_after: after,
            collision_center: self.center,
            velocity_midpoint: midpoint,
            momentum_change,
            external_impulse,
            guide_couple,
            external_work,
            kinetic_energy_change,
            momentum_residual,
            work_residual,
        })
    }
}

impl PreparedPlanes2 {
    fn finish_translations(
        &self,
        motion: DryLubricatedDiskExchange2,
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<DryTranslatingPlanesExchange2, DryDiskError2> {
        let mut planes = Vec::with_capacity(self.planes.len());
        for (plane, &film) in self.planes.iter().zip(&motion.planes) {
            checkpoint(cancel)?;
            let Some(translation) = &plane.translation else {
                return Err(DryDiskError2::Invalid("missing finite-plane response"));
            };
            planes.push(translation.finish(
                film,
                motion.motion.native.velocity_midpoint,
                origin,
            )?);
        }
        let mut combined = motion.motion.combined_momentum_residual;
        for a in 0..3 {
            let mut sum = Dd::from_f64(combined[a]) + Dd::from_f64(motion.body_impulse[a]);
            for plane in &planes {
                sum = sum + Dd::from_f64(plane.momentum_change[a])
                    - Dd::from_f64(plane.external_impulse[a])
                    - Dd::from_f64(if a == 2 { plane.guide_couple } else { 0.0 });
            }
            combined[a] = sum.to_f64();
        }
        if !combined.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid(
                "finite-plane aggregate balance overflowed",
            ));
        }
        Ok(DryTranslatingPlanesExchange2 {
            motion,
            planes,
            combined_momentum_residual: combined,
        })
    }
}

pub(super) fn admitted_velocity(velocity: [f64; 2]) -> bool {
    velocity.into_iter().all(f64::is_finite) && dot(velocity, velocity) < 0.03
}

pub(super) fn normal_exchange(
    normal: [f64; 2],
    disk_center: [f64; 2],
    reference: UnderResolvedLubricationExchange2,
    response: TranslationResponse2,
    origin: [f64; 2],
) -> Result<UnderResolvedLubricationExchange2, DryDiskError2> {
    let impulse = if response.normal_impulse == 0.0 {
        [0.0; 2]
    } else {
        normal.map(|n| response.normal_impulse * n)
    };
    let angular = cross(subtract(disk_center, origin), impulse);
    let work = response.normal_impulse * response.normal_rate;
    if !impulse
        .into_iter()
        .chain([angular, work, response.normal_rate])
        .all(f64::is_finite)
        || work > 0.0
    {
        return Err(DryDiskError2::Invalid(
            "finite-plane normal response overflowed",
        ));
    }
    Ok(UnderResolvedLubricationExchange2 {
        first_solid_impulse: impulse,
        second_solid_impulse: impulse.map(|j| if j == 0.0 { 0.0 } else { -j }),
        first_solid_angular_impulse: angular,
        second_solid_angular_impulse: if angular == 0.0 { 0.0 } else { -angular },
        normal_gap_velocity: response.normal_rate,
        solid_pair_work: work,
        ..reference
    })
}

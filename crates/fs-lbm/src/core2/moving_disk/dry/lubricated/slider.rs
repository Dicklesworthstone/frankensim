//! Finite tangential wall recoil condensed into the original dry-disk solve.

use super::*;

/// A finite mass sliding parallel to its attached infinite-plane approximation.
///
/// The guide fixes normal position and orientation. It receives the normal
/// film reaction and the full reaction couple; those constraints do no work.
/// Only tangential translation has finite inertia. The supplied mass is per
/// out-of-plane depth, not inferred from the infinite geometric surface.
/// No resolved-fluid wall, rotary inertia, contact or thermal state is added.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteLubricationSlider2 {
    plane: LubricationPlane2,
    mass: f64,
    center: [f64; 2],
    tangential_velocity: f64,
}

impl FiniteLubricationSlider2 {
    /// Bind a stationary plane descriptor, positive mass, center of mass and
    /// signed velocity along tangent [-normal_y, normal_x], in lattice units.
    /// The center may be offset from the plane. Its offset is retained during
    /// travel; translation along an infinite plane leaves its geometry exact.
    pub fn new(
        plane: LubricationPlane2,
        mass: f64,
        center: [f64; 2],
        tangential_velocity: f64,
    ) -> Result<Self, DryDiskError2> {
        if plane.velocity != [0.0; 2]
            || !(mass.is_finite() && mass > 0.0)
            || !center.into_iter().all(f64::is_finite)
            || !admitted_speed(tangential_velocity)
        {
            return Err(DryDiskError2::Invalid(
                "finite slider needs stationary geometry, positive mass and finite low-Mach motion",
            ));
        }
        Ok(Self {
            plane,
            mass,
            center,
            tangential_velocity,
        })
    }

    /// Fixed plane geometry and original normal/tangent/spin film scales.
    pub const fn plane(self) -> LubricationPlane2 {
        self.plane
    }
    /// Physical translating mass per unit out-of-plane depth.
    pub const fn mass(self) -> f64 {
        self.mass
    }
    /// Accepted center of mass in lattice coordinates.
    pub const fn center(self) -> [f64; 2] {
        self.center
    }
    /// Accepted signed speed along [-normal_y, normal_x].
    pub const fn tangential_velocity(self) -> f64 {
        self.tangential_velocity
    }
    /// Actual accepted world velocity; the guide permits no normal motion.
    pub fn velocity(self) -> [f64; 2] {
        self.tangent().map(|t| self.tangential_velocity * t)
    }
    fn tangent(self) -> [f64; 2] {
        [-self.plane.normal[1], self.plane.normal[0]]
    }
}

/// One finite slider's actual trajectory, storage and external guide account.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct SliderLubricationExchange2 {
    /// Complete accepted endpoint, also assigned to the supplied slider.
    pub state_after: FiniteLubricationSlider2,
    /// Center at the common disk/fluid/film impulse, between two free half-drifts.
    pub collision_center: [f64; 2],
    /// Signed solved tangential midpoint velocity.
    pub velocity_midpoint: f64,
    /// Actual [Px, Py, Lz] increment about the same requested fixed origin.
    pub momentum_change: [f64; 3],
    /// Authored tangential impulse, including its collision-center moment.
    pub external_impulse: [f64; 3],
    /// Impulse DELIVERED by the guide to this slider, about the fixed origin.
    /// This is an external port for the disk/fluid/slider system.
    pub guide_impulse: [f64; 3],
    /// Guide's center couple, separately retaining the distributed-traction moment.
    pub guide_center_couple: f64,
    /// Work of the authored tangential impulse at the solved midpoint.
    pub external_work: f64,
    /// Actual change of this finite slider's kinetic storage.
    pub kinetic_energy_change: f64,
    /// Actual change minus film, authored and guide impulse.
    pub momentum_residual: [f64; 3],
    /// Kinetic change minus authored work minus received plane film work.
    /// The ideal guide's normal/couple work is identically zero.
    pub work_residual: f64,
}

/// One native dry-disk/fluid step and simultaneous finite-slider motion.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct DryLubricatedSlidersExchange2 {
    /// Original disk/fluid and per-plane films. Here plane work enters finite
    /// slider storage; it is not an infinite reservoir or a second film impulse.
    pub motion: DryLubricatedDiskExchange2,
    /// Finite responses in exactly the supplied slider/plane order.
    pub sliders: Vec<SliderLubricationExchange2>,
    /// Exterior fluid plus dry disk plus every finite slider, after subtracting
    /// all authored loads, guide reactions and original fluid/boundary sources.
    pub combined_momentum_residual: [f64; 3],
}

impl DryDisk2 {
    /// Advance disk, actual populations and up to 64 finite guided plane sliders.
    ///
    /// Each external slider impulse acts tangentially at its own center and
    /// excludes film/guide reactions. For the original sliding coefficient k,
    /// condense the slider midpoint exactly:
    /// k_eff = k/(1+k/(2m)), u_free = u_before+J_external/(2m).
    /// The disk receives -k_eff*(v_disk,t,mid-u_free) in EVERY original
    /// three-unknown dry-fluid candidate. True slider midpoint work uses the
    /// original resistance; the cylinder spin coefficient remains unchanged.
    /// No later reaction kick or extra fluid collision is performed.
    ///
    /// Normal position/orientation stay fixed, so the original positive-gap
    /// checks remain exact under slider travel. Guide normal impulses and
    /// couples are explicit, zero-work external ports. This does not supply
    /// unconstrained wall recoil, finite face extents, resolved wall geometry,
    /// contact or a heat destination. Implicit midpoint is passive, but very
    /// stiff unresolved films may reverse relative velocity within a tick;
    /// time resolution is still required to resolve physical relaxation.
    ///
    /// Every response and endpoint is admitted before publication. Failure
    /// or cancellation preserves all accepted populations, disk and sliders.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubrication_sliders(
        &mut self,
        grid: &mut Grid,
        external_impulse: [f64; 2],
        external_couple: f64,
        moment_origin: [f64; 2],
        dynamic_viscosity: f64,
        sliders: &mut [FiniteLubricationSlider2],
        slider_external_impulses: &[f64],
        config: DryDiskCouplingConfig2,
        cancel: &CancelGate,
    ) -> Result<DryLubricatedSlidersExchange2, DryDiskError2> {
        checkpoint(cancel)?;
        config.admit()?;
        if !(dynamic_viscosity.is_finite() && dynamic_viscosity > 0.0)
            || !moment_origin.into_iter().all(f64::is_finite)
            || sliders.len() > 64
            || slider_external_impulses.len() != sliders.len()
        {
            return Err(DryDiskError2::Invalid(
                "invalid slider viscosity, origin, count or load count",
            ));
        }
        let descriptors: Vec<_> = sliders.iter().map(|slider| slider.plane).collect();
        let mut prepared = PreparedPlanes2::new(self, dynamic_viscosity, &descriptors, cancel)?;
        for ((plane, &slider), &authored) in prepared
            .planes
            .iter_mut()
            .zip(sliders.iter())
            .zip(slider_external_impulses)
        {
            checkpoint(cancel)?;
            plane.slider = Some(PreparedSlider2::new(
                slider,
                authored,
                plane.tangential_resistance,
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
            |receipt| prepared.finish_sliders(receipt, moment_origin, cancel),
        )?;
        // All fallible work precedes the existing native publication. These
        // copies publish exactly the already-admitted finite endpoints.
        for (slider, response) in sliders.iter_mut().zip(&receipt.sliders) {
            *slider = response.state_after;
        }
        Ok(receipt)
    }
}

/// Shared finite-mass midpoint elimination for each orthogonal film mode.
pub(super) struct CondensedResistance2 {
    pub(super) effective: f64,
    pub(super) retained: f64,
    pub(super) recoil: f64,
}

impl CondensedResistance2 {
    pub(super) fn new(m: f64, k: f64) -> Result<Self, DryDiskError2> {
        if !(m.is_finite() && m > 0.0 && k.is_finite() && k >= 0.0) {
            return Err(DryDiskError2::Invalid(
                "invalid finite film mass or resistance",
            ));
        }
        // Never form k/(2m) on the stiff/light branch: it can overflow even
        // though the condensed resistance has the finite limit 2m.
        let (effective, retained, recoil) = if k <= m {
            let ratio = 0.5 * (k / m);
            (
                k / (1.0 + ratio),
                1.0 / (1.0 + ratio),
                ratio / (1.0 + ratio),
            )
        } else {
            let ratio = m / k;
            (
                m / (0.5 + ratio),
                ratio / (0.5 + ratio),
                0.5 / (0.5 + ratio),
            )
        };
        Ok(Self {
            effective,
            retained,
            recoil,
        })
    }

    pub(super) fn midpoint(&self, free: f64, speed: f64) -> f64 {
        self.retained * free + self.recoil * speed
    }
}

pub(super) struct PreparedSlider2 {
    initial: FiniteLubricationSlider2,
    center: [f64; 2],
    external: f64,
    pub(super) free_midpoint: f64,
    pub(super) effective_resistance: f64,
    pub(super) retained_fraction: f64,
    recoil_fraction: f64,
}

impl PreparedSlider2 {
    fn new(
        initial: FiniteLubricationSlider2,
        external: f64,
        k: f64,
    ) -> Result<Self, DryDiskError2> {
        let m = initial.mass;
        let CondensedResistance2 {
            effective: effective_resistance,
            retained: retained_fraction,
            recoil: recoil_fraction,
        } = CondensedResistance2::new(m, k)?;
        let free_midpoint = initial.tangential_velocity + (0.5 * external) / m;
        let velocity = initial.velocity();
        let center = core::array::from_fn(|a| (0.5 * velocity[a]) + initial.center[a]);
        if !external.is_finite()
            || !free_midpoint.is_finite()
            || !center.into_iter().all(f64::is_finite)
            || ![effective_resistance, retained_fraction, recoil_fraction]
                .into_iter()
                .all(f64::is_finite)
        {
            return Err(DryDiskError2::Invalid("finite slider candidate overflowed"));
        }
        Ok(Self {
            initial,
            center,
            external,
            free_midpoint,
            effective_resistance,
            retained_fraction,
            recoil_fraction,
        })
    }

    pub(super) fn midpoint(&self, disk_tangent_speed: f64) -> f64 {
        self.retained_fraction * self.free_midpoint + self.recoil_fraction * disk_tangent_speed
    }

    fn finish(
        &self,
        film: PlaneLubricationExchange2,
        disk_velocity: [f64; 2],
        origin: [f64; 2],
    ) -> Result<SliderLubricationExchange2, DryDiskError2> {
        let initial = self.initial;
        let tangent = initial.tangent();
        let midpoint = self.midpoint(dot(disk_velocity, tangent));
        let endpoint = 2.0 * midpoint - initial.tangential_velocity;
        let center = core::array::from_fn(|a| (0.5 * endpoint * tangent[a]) + self.center[a]);
        if !admitted_speed(midpoint) {
            return Err(DryDiskError2::Invalid(
                "finite slider midpoint exceeds low-Mach limit",
            ));
        }
        let after = FiniteLubricationSlider2::new(initial.plane, initial.mass, center, endpoint)?;
        let delta = tangent.map(|t| initial.mass * (endpoint - initial.tangential_velocity) * t);
        let old_momentum = initial.velocity().map(|v| initial.mass * v);
        let momentum_change = [
            delta[0],
            delta[1],
            cross(subtract(center, origin), delta)
                + cross(subtract(center, initial.center), old_momentum),
        ];
        let arm = subtract(self.center, origin);
        let external = tangent.map(|t| self.external * t);
        let external_impulse = [external[0], external[1], cross(arm, external)];
        let plane_force = [film.plane_impulse[0], film.plane_impulse[1]];
        let guide_normal = -dot(plane_force, initial.plane.normal);
        let guide_force = initial.plane.normal.map(|n| guide_normal * n);
        let guide_center_couple =
            -film.plane_couple - cross(subtract(film.plane_footpoint, self.center), plane_force);
        let guide_impulse = [
            guide_force[0],
            guide_force[1],
            cross(arm, guide_force) + guide_center_couple,
        ];
        let momentum_residual = core::array::from_fn(|a| {
            momentum_change[a] - film.plane_impulse[a] - external_impulse[a] - guide_impulse[a]
        });
        let external_work = self.external * midpoint;
        let kinetic_energy_change = 0.5
            * initial.mass
            * (endpoint - initial.tangential_velocity)
            * (endpoint + initial.tangential_velocity);
        let work_residual = kinetic_energy_change - external_work - film.plane_work;
        let scale = norm(momentum_change)
            .max(norm(film.plane_impulse))
            .max(norm(external_impulse))
            .max(norm(guide_impulse))
            .max(norm([
                old_momentum[0],
                old_momentum[1],
                cross(arm, old_momentum),
            ]));
        let energy_scale = kinetic_energy_change
            .abs()
            .max(external_work.abs())
            .max(film.plane_work.abs())
            .max(initial.mass * initial.tangential_velocity.powi(2))
            .max(initial.mass * endpoint.powi(2));
        if !momentum_change
            .into_iter()
            .chain(external_impulse)
            .chain(guide_impulse)
            .chain(momentum_residual)
            .chain([
                guide_center_couple,
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
                "finite slider momentum/work balance is unresolved",
            ));
        }
        Ok(SliderLubricationExchange2 {
            state_after: after,
            collision_center: self.center,
            velocity_midpoint: midpoint,
            momentum_change,
            external_impulse,
            guide_impulse,
            guide_center_couple,
            external_work,
            kinetic_energy_change,
            momentum_residual,
            work_residual,
        })
    }
}

impl PreparedPlanes2 {
    pub(in super::super) fn finish_sliders(
        &self,
        motion: DryLubricatedDiskExchange2,
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<DryLubricatedSlidersExchange2, DryDiskError2> {
        let mut sliders = Vec::new();
        for (plane, &film) in self.planes.iter().zip(&motion.planes) {
            if let Some(slider) = &plane.slider {
                checkpoint(cancel)?;
                sliders.push(slider.finish(
                    film,
                    motion.motion.native.velocity_midpoint,
                    origin,
                )?);
            }
        }
        let mut combined = motion.motion.combined_momentum_residual;
        if !sliders.is_empty() {
            for a in 0..3 {
                let mut sum = Dd::from_f64(combined[a]) + Dd::from_f64(motion.body_impulse[a]);
                for slider in &sliders {
                    sum = sum + Dd::from_f64(slider.momentum_change[a])
                        - Dd::from_f64(slider.external_impulse[a])
                        - Dd::from_f64(slider.guide_impulse[a]);
                }
                combined[a] = sum.to_f64();
            }
        }
        if !combined.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid(
                "finite slider aggregate balance overflowed",
            ));
        }
        Ok(DryLubricatedSlidersExchange2 {
            motion,
            sliders,
            combined_momentum_residual: combined,
        })
    }
}

fn admitted_speed(speed: f64) -> bool {
    speed.is_finite() && speed * speed < 0.03
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_slider_condensation_retains_a_bounded_nonzero_stiff_limit() {
        let plane = LubricationPlane2::new([0.0; 2], [0.0, 1.0], [0.0; 2], 0.2, 0.01).unwrap();
        for (m, k) in [(2.0, 0.0), (2.0, 0.1), (0.01, 10.0), (1e-300, 1e308)] {
            let slider = FiniteLubricationSlider2::new(plane, m, [0.0; 2], 0.004).unwrap();
            let prepared = PreparedSlider2::new(slider, 0.0, k).unwrap();
            assert!(prepared.effective_resistance.is_finite());
            assert!(prepared.effective_resistance <= k);
            assert!(prepared.effective_resistance <= 2.0 * m);
            if k > 0.0 {
                assert!(prepared.effective_resistance > 0.0);
            }
            let midpoint = prepared.midpoint(-0.003);
            assert!((-0.003..=0.004).contains(&midpoint));
            if m == 1e-300 {
                assert_eq!(prepared.effective_resistance, 2.0 * m);
                assert_eq!(midpoint, -0.003);
            }
        }
    }
}

//! Prescribed planar surfaces for the shared normal-film resistance solve.

use super::*;

/// One infinite planar surface with constant prescribed translation.
///
/// The normal points into the fluid half-space occupied by the spheres. The
/// surface supplies a missing-film correction only: it does not create voxel
/// geometry, impose a resolved no-slip condition, or move the fluid boundary.
/// A complete wall model must supply that resolved boundary separately. This
/// is particularly useful beside an existing stationary halfway voxel wall.
///
/// Wall and sphere-pair resolution cutoffs need independent calibration. The
/// wall scales therefore belong to this surface; the dynamic viscosity and
/// combined active-pair budget come from [`SphereLubrication3`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LubricationPlane3 {
    pub(super) point: [f64; 3],
    pub(super) normal: [f64; 3],
    pub(super) velocity: [f64; 3],
    pub(super) cutoff_gap: f64,
    pub(super) roughness_gap: f64,
    pub(super) tangential: Option<(f64, f64)>,
}

impl LubricationPlane3 {
    /// Admit a finite point, nonzero normal, low-Mach translation and positive
    /// `roughness_gap <= cutoff_gap`, all in lattice units. The normal is
    /// normalized once. For a moving-sphere step, `point` is at the tick START;
    /// for a frozen-geometry exchange it is at the exchange pose.
    ///
    /// # Errors
    /// Refuses invalid geometry, unsupported speeds and unrepresentable scales.
    pub fn new(
        point: [f64; 3],
        normal: [f64; 3],
        velocity: [f64; 3],
        cutoff_gap: f64,
        roughness_gap: f64,
    ) -> Result<Self> {
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        if !point.into_iter().all(f64::is_finite)
            || !(length.is_finite() && length > 0.0)
            || !low_mach(velocity)
            || !(cutoff_gap.is_finite() && cutoff_gap > 0.0)
            || !(roughness_gap.is_finite() && roughness_gap > 0.0)
            || roughness_gap > cutoff_gap
        {
            return Err(invalid(
                "invalid lubrication plane geometry, speed or scales",
            ));
        }
        Ok(Self {
            point,
            normal: normal.map(|x| x / length),
            velocity,
            cutoff_gap,
            roughness_gap,
            tangential: None,
        })
    }

    /// Enable the leading logarithmic translation/parallel-spin wall correction.
    ///
    /// Both tangent directions retain the reciprocal force/torque block. With
    /// `s = n cross t`, `u = (v-Vwall).t`, `w = R*omega.s`, its resistance is
    /// `k * [[4,-1],[-1,4]]`, where `k = 4*pi*mu*R/5 * ln(cutoff/h_eff)`.
    /// The positive eigenvalues 3k and 5k retain finite rolling-mode resistance;
    /// this is not a dashpot acting only on contact-point slip. Rotation about
    /// the normal has no leading singular resistance and gets no correction.
    ///
    /// These separate tangential cutoff/roughness inputs use lattice lengths.
    /// At solve time the cutoff MUST be <= 0.02 times every supplied sphere's
    /// radius, retaining the small-gap regime of the Goldman-Cox-Brenner
    /// asymptotics (coefficients reproduced in Ryu et al., Phys. Rev. E 99,
    /// 022802, equations 7-10, https://doi.org/10.1103/PhysRevE.99.022802).
    /// Actual gaps remain strictly positive. No finite-gap complete
    /// resistance, cavitation, contact friction or experimental calibration is
    /// implied. Each active tangent block adds four scalar resistance modes.
    ///
    /// # Errors
    /// Refuses nonfinite/nonpositive scales or roughness above the cutoff.
    pub fn with_tangential(mut self, cutoff_gap: f64, roughness_gap: f64) -> Result<Self> {
        if !(cutoff_gap.is_finite() && cutoff_gap > 0.0)
            || !(roughness_gap.is_finite() && roughness_gap > 0.0)
            || roughness_gap > cutoff_gap
        {
            return Err(invalid("invalid tangential lubrication scales"));
        }
        self.tangential = Some((cutoff_gap, roughness_gap));
        Ok(self)
    }

    /// A point on the plane at the supplied window's start.
    pub fn point(self) -> [f64; 3] {
        self.point
    }

    /// Unit normal into the permitted half-space.
    pub fn normal(self) -> [f64; 3] {
        self.normal
    }

    /// Constant prescribed world translation in lattice velocity units.
    pub fn velocity(self) -> [f64; 3] {
        self.velocity
    }

    /// Advance this prescribed geometry for a subsequent window.
    ///
    /// # Errors
    /// Refuses negative/nonfinite time or overflow, without changing this value.
    pub fn translated(self, duration: f64) -> Result<Self> {
        let point = core::array::from_fn(|i| self.velocity[i].mul_add(duration, self.point[i]));
        if !(duration.is_finite() && duration >= 0.0) || !point.into_iter().all(f64::is_finite) {
            return Err(invalid("lubrication plane translation is unrepresentable"));
        }
        Ok(Self { point, ..self })
    }

    pub(super) fn distance(self, center: [f64; 3]) -> f64 {
        dot(
            self.normal,
            core::array::from_fn(|i| center[i] - self.point[i]),
        )
    }

    pub(super) fn tangential_modes(
        self,
        radius: f64,
        gap: f64,
        viscosity: f64,
    ) -> Result<Option<[TangentialMode; 4]>> {
        let Some((cutoff, roughness)) = self.tangential else {
            return Ok(None);
        };
        if cutoff > 0.02 * radius {
            return Err(invalid(
                "tangential lubrication cutoff exceeds 0.02 sphere radius",
            ));
        }
        let effective = gap.max(roughness);
        if effective >= cutoff {
            return Ok(None);
        }
        let k = (4.0 * core::f64::consts::PI / 5.0)
            * viscosity
            * radius
            * fs_math::det::ln(cutoff / effective);
        if !(k.is_finite() && k > 0.0) {
            return Err(invalid(
                "tangential lubrication resistance is unrepresentable",
            ));
        }
        // Choose a well-conditioned auxiliary axis. The two tangent directions
        // are isotropic: their sum is independent of this coordinate gauge.
        let axis = (0..3)
            .min_by(|&a, &b| self.normal[a].abs().total_cmp(&self.normal[b].abs()))
            .expect("three axes");
        let mut hint = [0.0; 3];
        hint[axis] = 1.0;
        let tangent = cross(self.normal, hint);
        let length = tangent[0].hypot(tangent[1]).hypot(tangent[2]);
        let tangent = tangent.map(|x| x / length);
        let basis = [tangent, cross(self.normal, tangent)];
        let modes = core::array::from_fn(|i| {
            let tangent = basis[i / 2];
            let spin = cross(self.normal, tangent);
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            TangentialMode {
                linear: tangent.map(|x| x * core::f64::consts::FRAC_1_SQRT_2),
                angular: spin.map(|x| sign * radius * x * core::f64::consts::FRAC_1_SQRT_2),
                resistance: if sign > 0.0 { 3.0 * k } else { 5.0 * k },
            }
        });
        if modes.iter().any(|mode| {
            !mode.resistance.is_finite() || !mode.angular.into_iter().all(f64::is_finite)
        }) {
            return Err(invalid("tangential lubrication modes overflowed"));
        }
        Ok(Some(modes))
    }
}

pub(super) struct TangentialMode {
    pub linear: [f64; 3],
    pub angular: [f64; 3],
    pub resistance: f64,
}

/// One prescribed surface's reaction, in the caller's plane order.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PlaneLubricationExchange3 {
    /// Missing-film impulse RECEIVED by this plane; opposite to the spheres.
    pub impulse: [f64; 3],
    /// Angular impulse about the caller's world origin, including the force at
    /// each actual plane footpoint and its balancing reaction couple.
    pub angular_impulse: [f64; 3],
    /// Sum of the free reaction couples at the individual plane footpoints.
    /// These retain total angular momentum for noncentral tangential traction.
    pub couple_impulse: [f64; 3],
    /// Work received by the prescribed surface, `impulse dot velocity`.
    /// Its driving source supplies the negative of this value.
    pub work: f64,
    /// Number of this plane's positive-resistance sphere films.
    pub active_spheres: usize,
    /// Number of this plane's active tangential/parallel-spin blocks.
    pub tangential_spheres: usize,
}

/// Native fluid/body exchange and separately attributable prescribed reactions.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct LubricatedPlanesExchange3<const N: usize> {
    /// One simultaneous native body/fluid solve, including all gap impulses.
    pub coupling: PartialSaturationBodiesExchange3<N>,
    /// Plane reactions in exactly the input order, including inactive surfaces.
    pub planes: Vec<PlaneLubricationExchange3>,
}

impl BoundaryGrid3 {
    /// Solve sphere-pair normal and sphere-plane missing resistance together.
    ///
    /// Body fractions and resolved wall geometry remain caller-owned. Each plane
    /// uses its own positive-gap/roughness/cutoff scales and the model's dynamic
    /// viscosity. All sphere pairs and every sphere/plane gap must be positive,
    /// even outside the cutoff. At most 64 planes and 64 TOTAL active films are
    /// admitted. Opt-in tangent blocks retain translation/rotation reciprocity
    /// and require their independent scalar-mode budget. There is one shared
    /// response solve, not sequential wall kicks.
    ///
    /// # Errors
    /// Inherits original immersed assembly admission; additionally refuses
    /// overlap, duplicate planes, exhausted film budgets or invalid solves.
    /// No body or population is published on any refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_partial_saturation_planes<const N: usize>(
        &mut self,
        fractions: [&[f64]; N],
        radii: [f64; N],
        bodies: &mut [PartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        planes: &[LubricationPlane3],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<LubricatedPlanesExchange3<N>> {
        self.step_with_partial_saturation_resistance_checked(
            fractions,
            bodies,
            external,
            origin,
            gate,
            Some((model, radii, planes)),
            |_, _, _| Ok(()),
        )
        .map(|(coupling, (), planes)| LubricatedPlanesExchange3 { coupling, planes })
    }
}

pub(in crate::d3q19::boundary) fn admit_planes(planes: &[LubricationPlane3]) -> Result<()> {
    if planes.len() > 64 {
        return Err(invalid("lubrication supports at most 64 prescribed planes"));
    }
    for (i, plane) in planes.iter().enumerate() {
        if planes[..i].iter().any(|other| {
            plane.normal == other.normal
                && dot(
                    plane.normal,
                    core::array::from_fn(|j| plane.point[j] - other.point[j]),
                ) == 0.0
        }) {
            return Err(invalid("duplicate lubrication plane"));
        }
    }
    Ok(())
}

// Relative center/plane distance is affine on each original rigid half-drift.
// Both endpoints positive therefore admit the WHOLE numerical segment; this
// does not enclose the unknown continuous physical trajectory.
pub(in crate::d3q19::boundary) fn admit_plane_positions<const N: usize>(
    planes: &[LubricationPlane3],
    centers: [[f64; 3]; N],
    radii: [f64; N],
    time: f64,
    gate: &CancelGate,
) -> Result<()> {
    for plane in planes {
        checkpoint(gate)?;
        let plane = plane.translated(time)?;
        for (center, radius) in centers.into_iter().zip(radii) {
            let gap = plane.distance(center) - radius;
            if !(gap.is_finite() && gap > 64.0 * f64::EPSILON * radius) {
                return Err(invalid("moving sphere meets prescribed lubrication plane"));
            }
        }
    }
    Ok(())
}

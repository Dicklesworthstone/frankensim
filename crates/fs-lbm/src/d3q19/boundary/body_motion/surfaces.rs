//! Moving finite film carriers on the original rigid and resistance owners.

use super::super::partial_saturation::lubrication::planes::admit_planes;
use super::super::{
    FiniteLubricationPlane3, LubricationPlane3, PlaneLubricationExchange3, SphereLubrication3,
};
use super::*;

mod shared;
pub use shared::{
    DryLubricatedCarriersExchange3, MovingCarrierExchange3, MovingCarrierSurface3,
    MovingLubricatedCarriersExchange3,
};

/// One body-fixed planar film patch carried by finite rigid inertia.
///
/// This is only the existing missing-film port. The carrier has no sampled
/// fluid interior or inferred shape, displaced volume, buoyancy or resolved
/// wall. Each descriptor is an independent carrier; several faces of one body
/// require a shared-carrier model and must not be represented by copied slots.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingFiniteLubricationPlane3 {
    body: MovingPartialSaturationBody3,
    point_body: [f64; 3],
    normal_body: [f64; 3],
    geometry: LubricationPlane3,
    cutoff_gap: f64,
    roughness_gap: f64,
    tangential: Option<(f64, f64)>,
}

impl MovingFiniteLubricationPlane3 {
    /// Attach a body-coordinate point and normal to original cached inertia.
    /// Coordinates use lattice units. The normal points toward the spheres.
    ///
    /// # Errors
    /// Refuses invalid plane scales, nonfinite pose/momenta, nonpositive mass,
    /// inadmissible center speed or unrepresentable full world inertia.
    pub fn new(
        body: MovingPartialSaturationBody3,
        point_body: [f64; 3],
        normal_body: [f64; 3],
        cutoff_gap: f64,
        roughness_gap: f64,
    ) -> Result<Self> {
        let local =
            LubricationPlane3::new(point_body, normal_body, [0.0; 3], cutoff_gap, roughness_gap)?;
        let mut value = Self {
            body,
            point_body,
            normal_body: local.normal(),
            geometry: local,
            cutoff_gap,
            roughness_gap,
            tangential: None,
        };
        value = value.at_body(body)?;
        // Original full-mobility admission also checks mass, velocity and spin.
        value.frozen([0.0; 6])?;
        Ok(value)
    }

    /// Enable the same native tangent/parallel-spin missing resistance.
    ///
    /// # Errors
    /// Refuses invalid tangent cutoff/roughness or a nonpositive missing law.
    pub fn with_tangential(mut self, cutoff_gap: f64, roughness_gap: f64) -> Result<Self> {
        self.geometry = self.geometry.with_tangential(cutoff_gap, roughness_gap)?;
        self.tangential = Some((cutoff_gap, roughness_gap));
        Ok(self)
    }
    /// Actual accepted pose, linear momentum and body-frame spin momentum.
    pub const fn body(self) -> MovingPartialSaturationBody3 {
        self.body
    }
    /// Accepted world point on the body-fixed plane.
    pub fn point(self) -> [f64; 3] {
        self.geometry.point()
    }
    /// Accepted unit world normal.
    pub fn normal(self) -> [f64; 3] {
        self.geometry.normal()
    }
    /// Local rigid velocity at a world point; original low-Mach admission applies.
    ///
    /// # Errors
    /// Refuses nonfinite points or an unrepresentable/inadmissible local speed.
    pub fn velocity_at(self, point: [f64; 3]) -> Result<[f64; 3]> {
        let velocity = point_velocity(
            self.body.velocity,
            self.body.angular_velocity_world()?,
            core::array::from_fn(|i| point[i] - self.body.center[i]),
        );
        if !low_mach(velocity) {
            return Err(invalid("moving finite surface point speed is inadmissible"));
        }
        Ok(velocity)
    }

    fn at_body(mut self, body: MovingPartialSaturationBody3) -> Result<Self> {
        let arm = array(native(body.orientation.rotate(vector(self.point_body)))?);
        let normal = array(native(body.orientation.rotate(vector(self.normal_body)))?);
        self.geometry = LubricationPlane3::new(
            core::array::from_fn(|i| body.center[i] + arm[i]),
            normal,
            body.velocity,
            self.cutoff_gap,
            self.roughness_gap,
        )?;
        if let Some((cutoff, roughness)) = self.tangential {
            self.geometry = self.geometry.with_tangential(cutoff, roughness)?;
        }
        self.body = body;
        Ok(self)
    }

    fn frozen(self, external: [f64; 6]) -> Result<LubricationPlane3> {
        let mut plane = FiniteLubricationPlane3::new(
            self.point(),
            self.normal(),
            world_body(self.body)?,
            self.cutoff_gap,
            self.roughness_gap,
        )?;
        if let Some((cutoff, roughness)) = self.tangential {
            plane = plane.with_tangential(cutoff, roughness)?;
        }
        plane.driven(external)
    }
}

/// Stable surface order shared by moving bodies, dry correction and film solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MovingLubricationSurface3 {
    /// Prescribed translation and optional axial spin; actuator owns received work.
    Prescribed(LubricationPlane3),
    /// Finite body-fixed plane; original free rotation advances its normal/inertia.
    Finite(MovingFiniteLubricationPlane3),
}

/// A finite carrier's complete two-half-drift trajectory account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingFinitePlaneExchange3 {
    /// Carrier center at the common resistance solve.
    pub midpoint_center: [f64; 3],
    /// Native free-drift orientation at that same midpoint.
    pub midpoint_orientation: So3,
    /// Authored impulse and its moment about the requested fixed world origin.
    pub external_impulse: [f64; 6],
    /// Actual full-trajectory change of linear and fixed-origin angular momentum.
    pub momentum_change: [f64; 6],
    /// Actual change minus authored and received film wrench.
    pub momentum_residual: [f64; 6],
    /// Actual translational and rotational kinetic-energy change over the tick.
    pub kinetic_energy_change: f64,
    /// Numerical free-rotation energy change, not viscous loss or heat.
    pub free_rotation_energy_change: f64,
    /// Kinetic change minus film/external work and free-rotation defect.
    pub work_residual: f64,
}

/// Original film exchange and, for a finite carrier, its actual moving account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingSurfaceExchange3 {
    /// Unchanged native film reaction and its finite or prescribed work account.
    pub film: PlaneLubricationExchange3,
    /// Actual evolving carrier trajectory; absent for a prescribed reservoir.
    pub finite: Option<MovingFinitePlaneExchange3>,
}

/// One shared moving-sphere/fluid tick and every stable-order surface response.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct MovingLubricatedSurfacesExchange3<const N: usize> {
    /// Original native sphere/fluid motion, including all common film impulses.
    pub motion: MovingPartialSaturationBodiesExchange3<N>,
    /// Surface reactions and actual finite trajectories in the input order.
    pub surfaces: Vec<MovingSurfaceExchange3>,
}

/// One physical dry-sphere/fluid tick and every stable-order surface response.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct DryLubricatedSurfacesExchange3<const N: usize> {
    /// The existing sphere/exterior-fluid subsystem account retains ALL received
    /// surface reactions. Finite carriers' momentum/work is reported below;
    /// their film work must not be relabeled as prescribed reservoir work.
    pub motion: DryBodiesExchange3<N>,
    /// Surface reactions and actual finite trajectories in the input order.
    pub surfaces: Vec<MovingSurfaceExchange3>,
}

impl BoundaryGrid3 {
    /// Advance physical dry spheres and mixed moving film carriers together.
    ///
    /// The original coupled 6N retained-fluid correction evaluates every film
    /// carrier from the SAME accepted surface state on each probe. Finite
    /// mobility therefore affects that correction, all other film reactions,
    /// and the one streamed fluid candidate. Only the converged, independently
    /// admitted result publishes. Original dry limits and the moving-surface
    /// bounds of [`Self::step_with_lubricated_spheres_and_surfaces`] apply.
    ///
    /// # Errors
    /// Inherits moving-surface admission and the original dry-body correction's
    /// convergence, condition, evaluation and sampled-volume bounds. Every
    /// refusal or cancellation preserves accepted bodies, surfaces and fluid.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dry_lubricated_spheres_and_surfaces<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        surfaces: &mut [MovingLubricationSurface3],
        surface_external: &[[f64; 6]],
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryLubricatedSurfacesExchange3<N>> {
        model.admit(radii)?;
        let prepared = PreparedSurfaces::new(surfaces, surface_external, gate)?;
        let shapes = radii.map(|r| BodyFrameEllipsoid3::sphere(r).expect("admitted radius"));
        let (report, candidate) = self.step_with_dry_partial_saturation_resistance(
            samplers,
            core::array::from_fn(|a| &shapes[a] as &dyn BodyFrameShape3),
            bodies,
            external,
            origin,
            config,
            gate,
            Some((model, &prepared.initial_planes)),
            Some(&prepared),
        )?;
        let candidate = candidate.expect("mixed surfaces prepare their accepted state");
        surfaces.copy_from_slice(&candidate.states);
        Ok(DryLubricatedSurfacesExchange3 {
            motion: report.motion,
            surfaces: candidate.reports,
        })
    }

    /// Advance spheres, mixed film carriers and fluid through one native solve.
    ///
    /// Every finite plane uses the same cached full-inertia free half-drifts as
    /// the spheres. Its actual midpoint geometry/mobility enters the original
    /// simultaneous resistance solve. Prescribed planes also advance to the
    /// accepted endpoint; callers must not translate them a second time.
    ///
    /// At most 64 surfaces are admitted. Loads are center impulse/couple pairs
    /// in world axes; prescribed slots require zero. Conservative angular-travel
    /// bounds admit both complete numerical half-drifts, and can refuse a valid
    /// near-gap trajectory. This is not exact CCD, carrier-carrier contact,
    /// resolved moving-wall geometry, a GCL or a heat deposition model.
    /// Any refusal preserves all accepted surfaces, sphere states and populations.
    ///
    /// # Errors
    /// Refuses mismatched loads, excessive surfaces, nonzero prescribed loads,
    /// unrepresentable rigid motion, unexcluded contact/duplicates, inadmissible
    /// native film/fluid solves, allocation failure or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_spheres_and_surfaces<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        surfaces: &mut [MovingLubricationSurface3],
        surface_external: &[[f64; 6]],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MovingLubricatedSurfacesExchange3<N>> {
        model.admit(radii)?;
        let prepared = PreparedSurfaces::new(surfaces, surface_external, gate)?;
        let shapes = radii.map(|r| BodyFrameEllipsoid3::sphere(r).expect("admitted radius"));
        let (report, candidate) = self.moving_partial_saturation_resistance_transition(
            samplers,
            core::array::from_fn(|a| &shapes[a] as &dyn BodyFrameShape3),
            bodies,
            external,
            origin,
            gate,
            Some((model, &prepared.initial_planes)),
            true,
            Some(&prepared),
        )?;
        let candidate = candidate.expect("mixed surfaces prepare their accepted state");
        surfaces.copy_from_slice(&candidate.states);
        Ok(MovingLubricatedSurfacesExchange3 {
            motion: report.motion,
            surfaces: candidate.reports,
        })
    }
}

pub(super) struct PreparedSurfaces {
    initial: Vec<MovingLubricationSurface3>,
    midpoint: Vec<MovingLubricationSurface3>,
    external: Vec<[f64; 6]>,
    pub initial_planes: Vec<LubricationPlane3>,
    pub midpoint_planes: Vec<LubricationPlane3>,
    shared: Option<shared::PreparedSharedSurfaces>,
}
pub(super) struct SurfaceCandidate {
    pub states: Vec<MovingLubricationSurface3>,
    pub reports: Vec<MovingSurfaceExchange3>,
    shared: Option<shared::SharedSurfaceCandidate>,
}

impl PreparedSurfaces {
    pub fn new(
        surfaces: &[MovingLubricationSurface3],
        external: &[[f64; 6]],
        gate: &CancelGate,
    ) -> Result<Self> {
        checkpoint(gate)?;
        if surfaces.len() > 64 || surfaces.len() != external.len() {
            return Err(invalid(
                "moving surface count or load count is inadmissible",
            ));
        }
        let mut initial = reserved(surfaces.len())?;
        let mut midpoint = reserved(surfaces.len())?;
        let mut loads = reserved(surfaces.len())?;
        let mut initial_planes = reserved(surfaces.len())?;
        let mut midpoint_planes = reserved(surfaces.len())?;
        for (&surface, &load) in surfaces.iter().zip(external) {
            checkpoint(gate)?;
            let (start, middle, mid) = match surface {
                MovingLubricationSurface3::Prescribed(p) => {
                    if load.into_iter().any(|x| x != 0.0) {
                        return Err(invalid(
                            "prescribed moving surface cannot receive an authored impulse",
                        ));
                    }
                    let middle = p.translated(0.5)?;
                    (p, middle, MovingLubricationSurface3::Prescribed(middle))
                }
                MovingLubricationSurface3::Finite(p) => {
                    let mid = p.at_body(drift(p.body, 0.5)?)?;
                    (
                        p.frozen(load)?,
                        mid.frozen(load)?,
                        MovingLubricationSurface3::Finite(mid),
                    )
                }
            };
            initial.push(surface);
            midpoint.push(mid);
            loads.push(load);
            initial_planes.push(start);
            midpoint_planes.push(middle);
        }
        admit_planes(&initial_planes)?;
        admit_planes(&midpoint_planes)?;
        admit_distinct_sweep(&initial, &midpoint, gate)?;
        Ok(Self {
            initial,
            midpoint,
            external: loads,
            initial_planes,
            midpoint_planes,
            shared: None,
        })
    }

    pub fn empty_candidate(&self) -> SurfaceCandidate {
        debug_assert!(self.initial_planes.is_empty());
        SurfaceCandidate {
            states: Vec::new(),
            reports: Vec::new(),
            shared: self
                .shared
                .as_ref()
                .map(|_| shared::SharedSurfaceCandidate::empty()),
        }
    }

    pub fn admit_first_half<const N: usize>(
        &self,
        before: [[f64; 3]; N],
        after: [[f64; 3]; N],
        radii: [f64; N],
        gate: &CancelGate,
    ) -> Result<()> {
        admit_spheres(&self.initial, &self.midpoint, before, after, radii, gate)
    }

    #[allow(clippy::too_many_lines)]
    pub fn finish<const N: usize>(
        &self,
        exchanges: &[PlaneLubricationExchange3],
        before: [[f64; 3]; N],
        after: [[f64; 3]; N],
        radii: [f64; N],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<SurfaceCandidate> {
        if let Some(shared) = &self.shared {
            return shared.finish(exchanges, before, after, radii, origin, gate);
        }
        if exchanges.len() != self.initial.len() {
            return Err(invalid("missing moving surface response"));
        }
        let mut states = reserved(exchanges.len())?;
        let mut reports = reserved(exchanges.len())?;
        // A second free half-drift starts with the EXCHANGED momentum, not the
        // incoming midpoint rate. Retain that state separately for sweep bounds.
        let mut exchanged = reserved(exchanges.len())?;
        for (index, &film) in exchanges.iter().enumerate() {
            checkpoint(gate)?;
            let (start, end, finite) = match (self.initial[index], self.midpoint[index]) {
                (
                    MovingLubricationSurface3::Prescribed(p),
                    MovingLubricationSurface3::Prescribed(mid),
                ) => {
                    if film.finite_body.is_some() {
                        return Err(invalid("unexpected finite surface response"));
                    }
                    (
                        MovingLubricationSurface3::Prescribed(mid),
                        MovingLubricationSurface3::Prescribed(p.translated(1.0)?),
                        None,
                    )
                }
                (MovingLubricationSurface3::Finite(p), MovingLubricationSurface3::Finite(mid)) => {
                    let (kicked, accepted, account) =
                        finish_carrier(p.body, mid.body, film, self.external[index], origin)?;
                    let end = mid.at_body(accepted)?;
                    let start = mid.at_body(kicked)?;
                    (
                        MovingLubricationSurface3::Finite(start),
                        MovingLubricationSurface3::Finite(end),
                        Some(account),
                    )
                }
                _ => return Err(invalid("moving surface kind changed during a trial")),
            };
            exchanged.push(start);
            states.push(end);
            reports.push(MovingSurfaceExchange3 { film, finite });
        }
        admit_spheres(&exchanged, &states, before, after, radii, gate)?;
        admit_distinct_sweep(&exchanged, &states, gate)?;
        checkpoint(gate)?;
        Ok(SurfaceCandidate {
            states,
            reports,
            shared: None,
        })
    }
}

// One rigid momentum update and trajectory account, whether the carrier owns
// one face or receives an aggregate wrench from several attached patches.
fn finish_carrier(
    initial: MovingPartialSaturationBody3,
    midpoint: MovingPartialSaturationBody3,
    film: PlaneLubricationExchange3,
    load: [f64; 6],
    origin: [f64; 3],
) -> Result<(
    MovingPartialSaturationBody3,
    MovingPartialSaturationBody3,
    MovingFinitePlaneExchange3,
)> {
    let response = film
        .finite_body
        .ok_or_else(|| invalid("missing finite surface response"))?;
    // Use the original directly accumulated COM spin impulse, avoiding
    // subtraction of close endpoint angular velocities.
    let impulse = core::array::from_fn(|i| response.spin_impulse[i] + load[i + 3]);
    let mut kicked = midpoint;
    kicked.velocity = response.body_after.velocity;
    kicked.angular_momentum_body = midpoint.angular_momentum_body
        + native(midpoint.orientation.inverse().rotate(vector(impulse)))?;
    let accepted = drift(kicked, 0.5)?;
    let external = at_origin(
        load,
        core::array::from_fn(|i| midpoint.center[i] - origin[i]),
    );
    let old_l = total_angular_momentum(initial, origin)?;
    let new_l = total_angular_momentum(accepted, origin)?;
    let change = core::array::from_fn(|i| {
        if i < 3 {
            initial.mass * (accepted.velocity[i] - initial.velocity[i])
        } else {
            new_l[i - 3] - old_l[i - 3]
        }
    });
    let residual = core::array::from_fn(|i| {
        change[i]
            - external[i]
            - if i < 3 {
                film.impulse[i]
            } else {
                film.angular_impulse[i - 3]
            }
    });
    let e0 = rotational_energy(initial.inertia, initial.angular_momentum_body)?;
    let em = rotational_energy(midpoint.inertia, midpoint.angular_momentum_body)?;
    let ek = rotational_energy(kicked.inertia, kicked.angular_momentum_body)?;
    let e1 = rotational_energy(accepted.inertia, accepted.angular_momentum_body)?;
    let free = (em - e0) + (e1 - ek);
    let kinetic = 0.5
        * initial.mass
        * (dot(accepted.velocity, accepted.velocity) - dot(initial.velocity, initial.velocity))
        + e1
        - e0;
    let work = kinetic - response.external_work - film.work - free;
    if !external
        .into_iter()
        .chain(change)
        .chain(residual)
        .chain([free, kinetic, work])
        .all(f64::is_finite)
    {
        return Err(invalid("moving finite surface momentum/work overflowed"));
    }
    Ok((
        kicked,
        accepted,
        MovingFinitePlaneExchange3 {
            midpoint_center: midpoint.center,
            midpoint_orientation: midpoint.orientation,
            external_impulse: external,
            momentum_change: change,
            momentum_residual: residual,
            kinetic_energy_change: kinetic,
            free_rotation_energy_change: free,
            work_residual: work,
        },
    ))
}

fn reserved<T>(count: usize) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("moving surface allocation refused"))?;
    Ok(values)
}

fn drift(
    mut body: MovingPartialSaturationBody3,
    duration: f64,
) -> Result<MovingPartialSaturationBody3> {
    body.center = core::array::from_fn(|i| body.center[i] + duration * body.velocity[i]);
    (body.orientation, body.angular_momentum_body) = native(body.inertia.momentum_step(
        body.orientation,
        body.angular_momentum_body,
        Vec3::new(0.0, 0.0, 0.0),
        duration,
    ))?;
    if !body.center.into_iter().all(f64::is_finite) || !low_mach(body.velocity) {
        return Err(invalid("moving surface drift is unrepresentable"));
    }
    Ok(body)
}

fn world_body(body: MovingPartialSaturationBody3) -> Result<PartialSaturationBody3> {
    let r = body.orientation.matrix();
    let tensor = r.compose(body.inertia.matrix()).compose(r.transpose());
    Ok(PartialSaturationBody3 {
        center: body.center,
        mass: body.mass,
        velocity: body.velocity,
        angular_velocity: body.angular_velocity_world()?,
        inertia_world: core::array::from_fn(|i| {
            core::array::from_fn(|j| tensor.m[3 * i.max(j) + i.min(j)])
        }),
    })
}

fn at_origin(impulse: [f64; 6], arm: [f64; 3]) -> [f64; 6] {
    let orbital = cross(arm, impulse[..3].try_into().unwrap());
    core::array::from_fn(|i| {
        if i < 3 {
            impulse[i]
        } else {
            impulse[i] + orbital[i - 3]
        }
    })
}
fn length(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

impl MovingLubricationSurface3 {
    fn geometry(self) -> LubricationPlane3 {
        match self {
            Self::Prescribed(p) => p,
            Self::Finite(p) => p.geometry,
        }
    }
    fn center(self) -> [f64; 3] {
        match self {
            Self::Prescribed(p) => p.point(),
            Self::Finite(p) => p.body.center,
        }
    }
    fn offset(self) -> f64 {
        match self {
            Self::Prescribed(_) => 0.0,
            Self::Finite(p) => dot(p.normal_body, p.point_body),
        }
    }
    // Each axis in the original symmetric split rotates for a total half tick.
    // |L_i| <= |L| gives this bound without an eigensolve or sampled CCD path.
    // It also bounds the exact spherical/axisymmetric special-case subflows.
    fn normal_travel(self) -> Result<f64> {
        let value = match self {
            Self::Prescribed(_) => 0.0,
            Self::Finite(p) => {
                let inertia = array(p.body.inertia.principal_moments());
                0.5 * length(array(p.body.angular_momentum_body))
                    * inertia.into_iter().map(f64::recip).sum::<f64>()
            }
        };
        if !value.is_finite() {
            return Err(invalid("moving plane angular-travel bound overflowed"));
        }
        Ok(value.min(2.0))
    }
}

fn admit_spheres<const N: usize>(
    surfaces: &[MovingLubricationSurface3],
    endpoints: &[MovingLubricationSurface3],
    before: [[f64; 3]; N],
    after: [[f64; 3]; N],
    radii: [f64; N],
    gate: &CancelGate,
) -> Result<()> {
    for (&surface, &endpoint) in surfaces.iter().zip(endpoints) {
        checkpoint(gate)?;
        let plane = surface.geometry();
        let travel = surface.normal_travel()?;
        for ((c0, c1), radius) in before.into_iter().zip(after).zip(radii) {
            let r0 = core::array::from_fn(|i| c0[i] - surface.center()[i]);
            let r1 = core::array::from_fn(|i| c1[i] - endpoint.center()[i]);
            let scale = length(r0).max(length(r1)) + surface.offset().abs() + radius;
            let floor = dot(plane.normal(), r0).min(dot(plane.normal(), r1))
                - surface.offset()
                - radius
                - travel * length(r0).max(length(r1));
            let end_gap = dot(
                endpoint.geometry().normal(),
                core::array::from_fn(|i| c1[i] - endpoint.geometry().point()[i]),
            ) - radius;
            if ![floor, end_gap, scale].into_iter().all(f64::is_finite)
                || floor <= 256.0 * f64::EPSILON * scale
                || end_gap <= 256.0 * f64::EPSILON * scale
            {
                return Err(invalid(
                    "moving sphere/film half-drift cannot exclude contact",
                ));
            }
        }
    }
    Ok(())
}

// A repeated coincident patch would duplicate resistance. For normal-changing
// pairs, prove distinct normals or separated offsets throughout the half-drift.
// This is intentionally conservative; no carrier-carrier collision is inferred.
fn admit_distinct_sweep(
    before: &[MovingLubricationSurface3],
    after: &[MovingLubricationSurface3],
    gate: &CancelGate,
) -> Result<()> {
    for a in 0..before.len() {
        checkpoint(gate)?;
        for b in 0..a {
            let na = before[a].geometry().normal();
            let nb = before[b].geometry().normal();
            let travel_a = before[a].normal_travel()?;
            let travel_b = before[b].normal_travel()?;
            if length(core::array::from_fn(|i| na[i] - nb[i]))
                > travel_a + travel_b + 256.0 * f64::EPSILON
            {
                continue;
            }
            let reference = before[a].center();
            let offset = |s: MovingLubricationSurface3, n| {
                dot(n, core::array::from_fn(|i| s.center()[i] - reference[i])) + s.offset()
            };
            let start = offset(before[a], na) - offset(before[b], nb);
            let finish = offset(after[a], na) - offset(after[b], nb);
            let rotation = |i: usize, travel: f64| {
                let start_arm = length(core::array::from_fn(|d| {
                    before[i].center()[d] - reference[d]
                }));
                let end_arm = length(core::array::from_fn(|d| {
                    after[i].center()[d] - reference[d]
                }));
                travel * start_arm.max(end_arm)
            };
            let bound = rotation(a, travel_a) + rotation(b, travel_b);
            let scale = start.abs() + finish.abs() + bound;
            if ![start, finish, bound, scale]
                .into_iter()
                .all(f64::is_finite)
                || start.signum() != finish.signum()
                || start.abs().min(finish.abs()) <= bound + 256.0 * f64::EPSILON * scale
            {
                return Err(invalid("moving film surfaces cannot be proven distinct"));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;

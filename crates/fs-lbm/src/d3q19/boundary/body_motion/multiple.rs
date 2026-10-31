//! Joint publication of several moving immersed bodies and one fluid tick.

use super::super::partial_saturation::lubrication::planes::{admit_plane_positions, admit_planes};
use super::super::{
    LubricationPlane3, PartialSaturationBodiesExchange3, PlaneLubricationExchange3,
    SphereLubrication3,
};
use super::*;

/// Accepted assembly trajectory with per-body pose, impulse and work accounts.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct MovingPartialSaturationBodiesExchange3<const N: usize> {
    /// Simultaneous cell-law response, optional normal lubrication and stationary
    /// walls, each counted once. Gap work/impulses belong to `coupling.lubrication`.
    pub coupling: PartialSaturationBodiesExchange3<N>,
    /// Per-body trajectory accounts in stable input order. The nested
    /// `coupling.walls` and `coupling.forcing` fields are zero; use the assembly
    /// wall and imposed-fluid-force accounts above.
    pub bodies: [MovingPartialSaturationExchange3; N],
}

/// Moving native bodies/fluid and prescribed planar-film reactions.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct MovingLubricatedPlanesExchange3<const N: usize> {
    /// The original rigid half-drifts and simultaneous fluid/body exchange.
    pub motion: MovingPartialSaturationBodiesExchange3<N>,
    /// Plane reactions in the input order, evaluated at their midpoint poses.
    pub planes: Vec<PlaneLubricationExchange3>,
}

impl BoundaryGrid3 {
    /// Advance all body positions, orientations and momenta with one fluid tick.
    ///
    /// Each body uses the native free half-drift / implicit midpoint exchange /
    /// free half-drift split, with its own cached full body-frame inertia and
    /// reusable sampler. Different shapes may share a lattice cell, but their
    /// enclosing spheres must stay strictly disjoint throughout each affine
    /// center half-step. This conservative test can reject separated elongated
    /// shapes. Wall/halo admission is inherited from the single-body owner.
    ///
    /// The shared native collision sees every midpoint pose simultaneously.
    /// All final poses, geometry, balance diagnostics and cancellation checks
    /// complete before any fluid or body is published. There is no grid clone
    /// and no sequential per-body stream. Sampling remains `O(N * grid cells)`
    /// including clearing the caller-owned fraction fields; pair admission is
    /// `O(N^2)`. This targets small rigid assemblies, not a sparse particle index.
    ///
    /// # Errors
    /// Refuses any invalid body, sampler/shape, overlapping swept enclosing
    /// spheres, inadmissible wall/halo path, native coupling failure or
    /// cancellation. Only sampling and fluid scratch may change on failure.
    /// This is retained-interior-fluid coupling with additional solid inertia,
    /// without solid contact, lubrication correction, displaced-fluid closure,
    /// fresh-cell remap, GCL or a full-system temporal-order certificate.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_moving_partial_saturation_bodies<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        shapes: [&dyn BodyFrameShape3; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MovingPartialSaturationBodiesExchange3<N>> {
        self.step_with_moving_partial_saturation_resistance(
            samplers, shapes, bodies, external, origin, gate, None,
        )
        .map(|report| report.motion)
    }

    /// Advance actual separated spheres with simultaneous missing-film
    /// resistance, native fluid exchange, and the original rigid half-drifts.
    ///
    /// The supplied radii define the sampled spheres; they are not independent
    /// curvature hints for arbitrary shapes. Resistance uses actual midpoint
    /// gaps. Opt-in unequal-sphere tangent/rotation modes use both original
    /// full-inertia spin responses. Original swept nonintersection and wall/halo
    /// admission remain.
    /// Missing-film impulses and viscous work are separate from resolved-fluid
    /// exchange in `coupling.lubrication`; all per-body balance residuals include
    /// both. The model is retained-interior-fluid plus additional solid inertia.
    ///
    /// # Errors
    /// Inherits original moving-assembly and [`SphereLubrication3`] admission.
    /// All body/fluid publication remains atomic, including late endpoint or
    /// cancellation failures. Sampler and fluid scratch may change on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_spheres<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MovingPartialSaturationBodiesExchange3<N>> {
        model.admit(radii)?;
        let spheres = radii
            .map(|r| BodyFrameEllipsoid3::sphere(r).expect("lubrication admitted every radius"));
        self.step_with_moving_partial_saturation_resistance(
            samplers,
            core::array::from_fn(|a| &spheres[a] as &dyn BodyFrameShape3),
            bodies,
            external,
            origin,
            gate,
            Some((model, &[])),
        )
        .map(|report| report.motion)
    }

    /// Advance actual spheres with simultaneous sphere-pair and prescribed-plane
    /// films, including each plane's optional reciprocal tangential/spin law.
    /// Plane points are at this tick's START and translate at their
    /// declared velocity. The exchange uses all midpoint gaps; actual start,
    /// midpoint and endpoint gaps admit both original affine relative drifts.
    ///
    /// These planes supply subgrid resistance and its external reaction/work;
    /// they do not move or create the resolved fluid boundaries. Existing voxel
    /// walls remain independently authoritative. Callers of moving prescribed
    /// surfaces must supply any required resolved-fluid wall model separately.
    ///
    /// # Errors
    /// Inherits [`Self::step_with_lubricated_spheres`] and plane admission. Any
    /// late plane overlap or cancellation leaves all populations/bodies intact.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_lubricated_spheres_and_planes<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        radii: [f64; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        model: SphereLubrication3,
        planes: &[LubricationPlane3],
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MovingLubricatedPlanesExchange3<N>> {
        model.admit(radii)?;
        admit_planes(planes)?;
        let spheres = radii
            .map(|r| BodyFrameEllipsoid3::sphere(r).expect("lubrication admitted every radius"));
        self.step_with_moving_partial_saturation_resistance(
            samplers,
            core::array::from_fn(|a| &spheres[a] as &dyn BodyFrameShape3),
            bodies,
            external,
            origin,
            gate,
            Some((model, planes)),
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn step_with_moving_partial_saturation_resistance<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        shapes: [&dyn BodyFrameShape3; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        origin: [f64; 3],
        gate: &CancelGate,
        lubrication: Option<(SphereLubrication3, &[LubricationPlane3])>,
    ) -> Result<MovingLubricatedPlanesExchange3<N>> {
        checkpoint(gate)?;
        if N == 0 || !origin.into_iter().all(f64::is_finite) {
            return Err(invalid(
                "moving assembly requires bodies and a finite origin",
            ));
        }
        let initial = *bodies;
        let planes = lubrication.map_or(&[][..], |(_, planes)| planes);
        let mut midpoint_planes = Vec::new();
        midpoint_planes
            .try_reserve_exact(planes.len())
            .map_err(|_| invalid("moving lubrication plane allocation refused"))?;
        for plane in planes {
            checkpoint(gate)?;
            midpoint_planes.push(plane.translated(0.5)?);
        }
        let mut midpoint = initial;
        let mut radii = [0.0; N];
        let mut volumes = [[0.0; 3]; N];
        let mut frozen: [PartialSaturationBody3; N] =
            core::array::from_fn(|a| PartialSaturationBody3 {
                center: initial[a].center,
                mass: initial[a].mass,
                inertia_world: [[0.0; 3]; 3],
                velocity: initial[a].velocity,
                angular_velocity: [0.0; 3],
            });
        for a in 0..N {
            checkpoint(gate)?;
            let body = initial[a];
            if samplers[a].dims != self.dimensions()
                || !(body.mass.is_finite() && body.mass > 0.0)
                || !body
                    .center
                    .into_iter()
                    .chain(external[a])
                    .all(f64::is_finite)
                || !low_mach(body.velocity)
            {
                return Err(invalid(
                    "moving assembly body, sampler or impulse is inadmissible",
                ));
            }
            radii[a] = shape_radius(shapes[a])?;
            midpoint[a].center = core::array::from_fn(|i| body.center[i] + 0.5 * body.velocity[i]);
            admit_sweep(self, body.center, midpoint[a].center, radii[a], gate)?;
            (midpoint[a].orientation, midpoint[a].angular_momentum_body) =
                native(body.inertia.momentum_step(
                    body.orientation,
                    body.angular_momentum_body,
                    Vec3::new(0.0, 0.0, 0.0),
                    0.5,
                ))?;
            let omega = midpoint[a].angular_velocity_world()?;
            volumes[a][0] = sample_geometry(
                samplers[a].dims,
                samplers[a].samples_per_axis,
                body.center,
                body.orientation,
                shapes[a],
                Some((body.velocity, body.angular_velocity_world()?)),
                None,
                gate,
            )?;
            volumes[a][1] = sample_geometry(
                samplers[a].dims,
                samplers[a].samples_per_axis,
                midpoint[a].center,
                midpoint[a].orientation,
                shapes[a],
                Some((body.velocity, omega)),
                Some(&mut samplers[a].fractions),
                gate,
            )?;
            let rotation = midpoint[a].orientation.matrix();
            let tensor = rotation
                .compose(body.inertia.matrix())
                .compose(rotation.transpose());
            frozen[a] = PartialSaturationBody3 {
                center: midpoint[a].center,
                mass: body.mass,
                inertia_world: core::array::from_fn(|i| {
                    core::array::from_fn(|j| tensor.m[3 * i.max(j) + i.min(j)])
                }),
                velocity: body.velocity,
                angular_velocity: omega,
            };
        }
        admit_pairs(
            initial.map(|b| b.center),
            midpoint.map(|b| b.center),
            radii,
            gate,
        )?;
        admit_plane_positions(planes, initial.map(|b| b.center), radii, 0.0, gate)?;
        admit_plane_positions(planes, midpoint.map(|b| b.center), radii, 0.5, gate)?;
        let frozen_before = frozen;
        let fields = core::array::from_fn(|a| samplers[a].fractions.as_slice());
        let (_, (accepted, report), planes) = self
            .step_with_partial_saturation_resistance_checked(
                fields,
                &mut frozen,
                external,
                origin,
                gate,
                lubrication.map(|(model, _)| (model, radii, midpoint_planes.as_slice())),
                |grid, exchanged, coupling| {
                    let mut accepted = initial;
                    let mut reports = [None; N];
                    for a in 0..N {
                        checkpoint(gate)?;
                        let body = initial[a];
                        let mid = midpoint[a];
                        let response = coupling.bodies[a];
                        let dw = core::array::from_fn(|i| {
                            exchanged[a].angular_velocity[i] - frozen_before[a].angular_velocity[i]
                        });
                        let spin_impulse = frozen_before[a].inertia_world.map(|row| dot(row, dw));
                        let after_exchange_momentum = mid.angular_momentum_body
                            + native(mid.orientation.inverse().rotate(vector(spin_impulse)))?;
                        let (orientation, angular_momentum_body) =
                            native(body.inertia.momentum_step(
                                mid.orientation,
                                after_exchange_momentum,
                                Vec3::new(0.0, 0.0, 0.0),
                                0.5,
                            ))?;
                        accepted[a] = MovingPartialSaturationBody3 {
                            center: core::array::from_fn(|i| {
                                mid.center[i] + 0.5 * exchanged[a].velocity[i]
                            }),
                            orientation,
                            velocity: exchanged[a].velocity,
                            angular_momentum_body,
                            ..body
                        };
                        admit_sweep(grid, mid.center, accepted[a].center, radii[a], gate)?;
                        volumes[a][2] = sample_geometry(
                            samplers[a].dims,
                            samplers[a].samples_per_axis,
                            accepted[a].center,
                            orientation,
                            shapes[a],
                            Some((accepted[a].velocity, accepted[a].angular_velocity_world()?)),
                            None,
                            gate,
                        )?;
                        let arm = core::array::from_fn(|i| mid.center[i] - origin[i]);
                        let orbital = cross(arm, external[a][..3].try_into().unwrap());
                        let external_angular_impulse =
                            core::array::from_fn(|i| external[a][i + 3] + orbital[i]);
                        let before_h = total_angular_momentum(body, origin)?;
                        let after_h = total_angular_momentum(accepted[a], origin)?;
                        let momentum_residual = core::array::from_fn(|i| {
                            if i < 3 {
                                body.mass * (accepted[a].velocity[i] - body.velocity[i])
                                    - external[a][i]
                                    - response.hydrodynamic.body_impulse[i]
                                    - coupling.lubrication.body_impulses[a][i]
                            } else {
                                after_h[i - 3]
                                    - before_h[i - 3]
                                    - external_angular_impulse[i - 3]
                                    - response.hydrodynamic.body_angular_impulse[i - 3]
                                    - coupling.lubrication.body_angular_impulses[a][i - 3]
                            }
                        });
                        let spin_before =
                            rotational_energy(body.inertia, body.angular_momentum_body)?;
                        let spin_midpoint =
                            rotational_energy(body.inertia, mid.angular_momentum_body)?;
                        let spin_exchanged =
                            rotational_energy(body.inertia, after_exchange_momentum)?;
                        let spin_after = rotational_energy(body.inertia, angular_momentum_body)?;
                        let kinetic_energy_change = 0.5
                            * body.mass
                            * (dot(accepted[a].velocity, accepted[a].velocity)
                                - dot(body.velocity, body.velocity))
                            + spin_after
                            - spin_before;
                        let free_rotation_energy_change =
                            (spin_midpoint - spin_before) + (spin_after - spin_exchanged);
                        let work_residual = kinetic_energy_change
                            - response.external_work
                            - response.hydrodynamic.body_work
                            - coupling.lubrication.body_work[a]
                            - free_rotation_energy_change;
                        if !volumes[a]
                            .into_iter()
                            .chain(external_angular_impulse)
                            .chain(momentum_residual)
                            .chain([
                                kinetic_energy_change,
                                free_rotation_energy_change,
                                work_residual,
                            ])
                            .all(f64::is_finite)
                        {
                            return Err(invalid("moving assembly endpoint account overflowed"));
                        }
                        reports[a] = Some(MovingPartialSaturationExchange3 {
                            coupling: PartialSaturationBodyExchange3 {
                                hydrodynamic: response.hydrodynamic,
                                walls: super::super::WallMomentumExchange3::default(),
                                periodic: super::super::PeriodicMomentumTransport3::default(),
                                forcing: super::super::PartialSaturationForcing3::default(),
                                velocity_midpoint: response.velocity_midpoint,
                                angular_velocity_midpoint: response.angular_velocity_midpoint,
                                kinetic_energy_change: response.kinetic_energy_change,
                                external_work: response.external_work,
                                momentum_residual: response.momentum_residual,
                            },
                            midpoint_center: mid.center,
                            midpoint_orientation: mid.orientation,
                            sampled_solid_volume: volumes[a],
                            external_angular_impulse,
                            momentum_residual,
                            kinetic_energy_change,
                            free_rotation_energy_change,
                            work_residual,
                        });
                    }
                    admit_pairs(
                        midpoint.map(|b| b.center),
                        accepted.map(|b| b.center),
                        radii,
                        gate,
                    )?;
                    admit_plane_positions(planes, accepted.map(|b| b.center), radii, 1.0, gate)?;
                    let report = MovingPartialSaturationBodiesExchange3 {
                        coupling: *coupling,
                        bodies: reports
                            .map(|report| report.expect("every assembly body was finalized")),
                    };
                    Ok((accepted, report))
                },
            )?;
        *bodies = accepted;
        Ok(MovingLubricatedPlanesExchange3 {
            motion: report,
            planes,
        })
    }
}

// Each half-step has affine centers. Minimize the squared relative center
// distance on that SAME time interval, rather than comparing unrelated times
// on two swept segments. A small floating-point margin rejects grazing cases;
// this is a numerical bounding-sphere admission, not an interval CCD theorem.
fn admit_pairs<const N: usize>(
    before: [[f64; 3]; N],
    after: [[f64; 3]; N],
    radii: [f64; N],
    gate: &CancelGate,
) -> Result<()> {
    for a in 0..N {
        checkpoint(gate)?;
        for b in 0..a {
            if b.is_multiple_of(TILE_CELLS) {
                checkpoint(gate)?;
            }
            let r = core::array::from_fn(|i| before[a][i] - before[b][i]);
            let d = core::array::from_fn(|i| {
                (after[a][i] - before[a][i]) - (after[b][i] - before[b][i])
            });
            let speed2 = dot(d, d);
            let t = if speed2 > 0.0 {
                (-dot(r, d) / speed2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let closest = core::array::from_fn(|i| r[i] + t * d[i]);
            let radius = radii[a] + radii[b];
            let distance2 = dot(closest, closest);
            let threshold = radius * radius * (1.0 + 64.0 * f64::EPSILON);
            if !speed2.is_finite()
                || !distance2.is_finite()
                || !threshold.is_finite()
                || distance2 <= threshold
            {
                return Err(invalid("moving assembly swept enclosing spheres meet"));
            }
        }
    }
    checkpoint(gate)
}

//! Simultaneous physical dry bodies with coupled retained-fluid corrections.

use super::super::surfaces::{PreparedSurfaces, SurfaceCandidate};
use super::*;
use crate::d3q19::{
    LubricationPlane3, PartialSaturationBodiesExchange3, PlaneLubricationExchange3,
    SphereLubrication3, SphereLubricationExchange3,
};

mod lubricated;
pub use lubricated::DryLubricatedPlanesExchange3;

/// One accepted fluid tick and all physical dry bodies, in caller order.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct DryBodiesExchange3<const N: usize> {
    /// Native simultaneous collision and motion. Its body loads include the
    /// solved interior corrections; its full-grid force/wall accounts occur once.
    pub native: MovingPartialSaturationBodiesExchange3<N>,
    /// Full-grid physical Guo inventories before/after streaming, counted once.
    pub fluid: [ImmersedMomentum3; 2],
    /// Per-body endpoint-weighted fictitious-fluid inventories before/after.
    pub interior: [[ImmersedMomentum3; 2]; N],
    /// Per-body RMS fluid-minus-rigid-point speed, before/after.
    pub interior_slip_rms: [[f64; 2]; N],
    /// Solved corrections about the requested fixed world origin.
    pub interior_correction_impulses: [[f64; 6]; N],
    /// Imposed force/torque on each distinct collision-midpoint interior.
    pub interior_forcing_impulses: [[f64; 6]; N],
    /// Full fluid forcing minus the SUM of all fictitious-interior sources.
    pub represented_fluid_forcing_impulse: [f64; 6],
    /// Authored body loads expressed about the requested fixed world origin.
    pub physical_external_impulses: [[f64; 6]; N],
    /// Native reactions plus measured interior changes minus interior forcing.
    pub corrected_hydrodynamic_impulses: [[f64; 6]; N],
    /// Actual dry-body linear/fixed-origin angular momentum changes.
    pub body_momentum_changes: [[f64; 6]; N],
    /// Per-body changes minus authored, corrected hydrodynamic and film impulses.
    pub body_momentum_residuals: [[f64; 6]; N],
    /// Exterior fluid plus all dry bodies minus all physical sources and fluxes,
    /// including the opposite reaction received by every film surface. Finite
    /// carriers are outside this sphere/exterior-fluid subsystem account.
    pub combined_momentum_residual: [f64; 6],
    /// Authored load work at each native midpoint rigid velocity.
    pub physical_external_work: [f64; N],
    /// Interior correction work; neither authored external power nor heat.
    pub interior_correction_work: [f64; N],
    /// All dry bodies plus exterior fluid MACROSCOPIC kinetic change.
    pub represented_kinetic_energy_change: f64,
    /// Accepted simultaneous Newton corrections.
    pub iterations: usize,
    /// Full native assembly evaluations, including probes and backtracks.
    pub candidate_evaluations: usize,
    /// Last estimated condition number of the full scaled 6N-by-6N Jacobian.
    pub jacobian_condition: Option<f64>,
}

// Adapts a borrowed heterogeneous shape to the original sized single-body API.
struct ShapeRef<'a>(&'a dyn BodyFrameShape3);
impl BodyFrameShape3 for ShapeRef<'_> {
    fn bounding_radius(&self) -> f64 {
        self.0.bounding_radius()
    }
    fn contains(&self, point: Vec3) -> bool {
        self.0.contains(point)
    }
}

struct AssemblyStep<'a, const N: usize> {
    shapes: [&'a dyn BodyFrameShape3; N],
    initial: [MovingPartialSaturationBody3; N],
    external: [[f64; 6]; N],
    origin: [f64; 3],
    arms: [[f64; 3]; N],
    radii: [f64; N],
    inside_before: [ImmersedMomentum3; N],
    interior_forcing: Option<[[f64; 6]; N]>,
    lubrication: Option<(SphereLubrication3, &'a [LubricationPlane3])>,
    surfaces: Option<&'a PreparedSurfaces>,
    maximum_relative_volume_change: f64,
    gate: &'a CancelGate,
}

struct AssemblyCandidate<const N: usize> {
    bodies: [MovingPartialSaturationBody3; N],
    native: MovingPartialSaturationBodiesExchange3<N>,
    planes: Vec<PlaneLubricationExchange3>,
    surfaces: Option<SurfaceCandidate>,
    interior: [ImmersedMomentum3; N],
    slip: [f64; N],
    residual: [[f64; 6]; N],
    target: [[f64; 6]; N],
}

impl BoundaryGrid3 {
    /// Advance 1 through 8 physical dry bodies and the resolved fluid together.
    ///
    /// Every original positive dry mass and full tensor is retained. The 6N
    /// corrections solve `J_a = Delta Pi_inside_a - J_force_inside_a` using
    /// physical Guo inventories from ONE simultaneous native fluid candidate.
    /// All cross-body Jacobian entries are retained: independently or
    /// sequentially solving bodies would miss momentum streamed between them.
    /// Each external WORLD impulse acts at that body's midpoint COM, followed
    /// by its center couple. Fluid forces are held through the tick; pressure
    /// and source-frame conventions match the single-body dry API.
    ///
    /// One body delegates exactly to that owner. Larger assemblies use bounded
    /// safeguarded Newton: each iteration needs 6N full candidate probes plus
    /// its backtracks; `maximum_evaluations` caps ALL candidates. The dense
    /// Jacobian has at most 48 squared entries. Native fluid/sampler scratch is
    /// reused; no grid backup, extra stream per body, or reduced inertia is used.
    /// Each candidate samples O(N * cells) geometry; this targets small assemblies.
    ///
    /// Swept enclosing spheres must stay disjoint. Endpoint sampled occupancy
    /// may share distinct subvolumes of a cell, but its sum must not exceed one.
    /// Collision-midpoint forcing is cached once per body, then subtracted once
    /// from the common full-fluid source. Sampled volume/slip limitations match
    /// the single-body owner; this does not add contact, lubrication or a GCL.
    ///
    /// # Errors
    /// Refuses empty/oversized assemblies, invalid geometry or original body
    /// data, overlapping sweeps/occupancy, unacceptable quadrature volume change,
    /// singular/ill-conditioned or unconverged solves, work exhaustion and
    /// cancellation. Every fallible check precedes the single atomic publication
    /// of all fluid populations and bodies. Only reusable scratch may change.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_dry_partial_saturation_bodies<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        shapes: [&dyn BodyFrameShape3; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryBodiesExchange3<N>> {
        self.step_with_dry_partial_saturation_resistance(
            samplers, shapes, bodies, external, origin, config, gate, None, None,
        )
        .map(|(report, _)| report.motion)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::d3q19::boundary::body_motion) fn step_with_dry_partial_saturation_resistance<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        shapes: [&dyn BodyFrameShape3; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external: [[f64; 6]; N],
        origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
        lubrication: Option<(SphereLubrication3, &[LubricationPlane3])>,
        surfaces: Option<&PreparedSurfaces>,
    ) -> Result<(DryLubricatedPlanesExchange3<N>, Option<SurfaceCandidate>)> {
        checkpoint(gate)?;
        config.admit()?;
        if !(1..=8).contains(&N) || !origin.into_iter().all(f64::is_finite) {
            return Err(invalid(
                "dry assembly requires 1 through 8 bodies and a finite origin",
            ));
        }
        // A lone sphere without planes has no film partner. Keep the original
        // single-dry arithmetic exactly, including its convergence policy.
        if N == 1 && lubrication.is_none_or(|(_, planes)| planes.is_empty()) {
            let single = self.step_with_dry_partial_saturation_body(
                &mut samplers[0],
                &ShapeRef(shapes[0]),
                &mut bodies[0],
                external[0],
                origin,
                config,
                gate,
            )?;
            return Ok((
                DryLubricatedPlanesExchange3 {
                    motion: singleton(single),
                    planes: Vec::new(),
                },
                surfaces.map(|_| SurfaceCandidate {
                    states: Vec::new(),
                    reports: Vec::new(),
                }),
            ));
        }
        let initial = *bodies;
        let mut radii = [0.0; N];
        let mut inside_before = [ImmersedMomentum3::default(); N];
        let mut slip_before = [0.0; N];
        let mut probe_scale = [0.0; N];
        let mut fluid_before = ImmersedMomentum3::default();
        let arms = core::array::from_fn(|a| {
            core::array::from_fn(|i| {
                initial[a].center[i] + 0.5 * initial[a].velocity[i] - origin[i]
            })
        });
        for a in 0..N {
            if samplers[a].dims != self.dimensions()
                || !(initial[a].mass.is_finite() && initial[a].mass > 0.0)
                || !external[a].into_iter().all(f64::is_finite)
            {
                return Err(invalid(
                    "dry assembly requires matching samplers, positive masses and finite loads",
                ));
            }
            radii[a] = shape_radius(shapes[a])?;
            samplers[a].sample(
                initial[a].center,
                initial[a].orientation,
                &ShapeRef(shapes[a]),
                gate,
            )?;
            let (all, inside, slip) =
                self.dry_inventory(false, a == 0, &samplers[a], initial[a], origin, gate)?;
            if a == 0 {
                fluid_before = all;
            }
            inside_before[a] = inside;
            slip_before[a] = slip;
            probe_scale[a] = (initial[a].mass + inside.mass) * 0.01;
            if !(inside.mass > 0.0 && probe_scale[a].is_finite() && probe_scale[a] > 0.0) {
                return Err(invalid(
                    "dry assembly needs resolved interiors and finite Jacobian scales",
                ));
            }
        }
        admit_inventory_occupancy(samplers, gate)?;
        let mut step = AssemblyStep {
            shapes,
            initial,
            external,
            origin,
            arms,
            radii,
            inside_before,
            interior_forcing: None,
            lubrication,
            surfaces,
            maximum_relative_volume_change: config.maximum_relative_volume_change,
            gate,
        };
        let mut correction = [[0.0; 6]; N];
        let mut evaluations = 0;
        let mut iterations = 0;
        let mut condition = None;
        let mut current = self.dry_assembly_candidate(
            samplers,
            &mut step,
            correction,
            &mut evaluations,
            config.maximum_evaluations,
        )?;
        let size = 6 * N;
        let mut jacobian = vec![0.0; size * size];
        let mut update = vec![0.0; size];
        loop {
            checkpoint(gate)?;
            let residual = residual_ratio(&current, correction, &step, config)?;
            if residual <= 1.0 {
                break;
            }
            if iterations == config.maximum_iterations {
                return Err(invalid("dry assembly correction did not converge"));
            }
            for column in 0..size {
                let (a, i) = (column / 6, column % 6);
                let delta = config.jacobian_step * (probe_scale[a] + correction[a][i].abs());
                if !(delta.is_finite() && delta > 0.0) {
                    return Err(invalid("dry assembly Jacobian perturbation overflowed"));
                }
                let mut probe = correction;
                probe[a][i] += delta;
                let candidate = self.dry_assembly_candidate(
                    samplers,
                    &mut step,
                    probe,
                    &mut evaluations,
                    config.maximum_evaluations,
                )?;
                // Every body row responds to this column's shared fluid step.
                for row in 0..size {
                    jacobian[size * row + column] = (candidate.residual[row / 6][row % 6]
                        - current.residual[row / 6][row % 6])
                        / delta;
                }
            }
            if !jacobian.iter().all(|x| x.is_finite()) {
                return Err(invalid("dry assembly Jacobian overflowed"));
            }
            let factor = fs_la::factor::lu(&jacobian, size)
                .map_err(|_| invalid("dry assembly correction Jacobian is singular"))?;
            let estimate = factor.condition_1(&jacobian);
            if !(estimate.is_finite() && estimate <= config.maximum_condition_number) {
                return Err(invalid(
                    "dry assembly correction Jacobian is ill-conditioned",
                ));
            }
            condition = Some(estimate);
            for (row, value) in update.iter_mut().enumerate() {
                *value = -current.residual[row / 6][row % 6];
            }
            factor.solve(&mut update);
            if !update.iter().all(|x| x.is_finite()) {
                return Err(invalid("dry assembly Newton update overflowed"));
            }
            let mut accepted = None;
            let mut alpha = 1.0;
            for _ in 0..=config.maximum_backtracks {
                let trial = core::array::from_fn(|a| {
                    core::array::from_fn(|i| correction[a][i] + alpha * update[6 * a + i])
                });
                if evaluations == config.maximum_evaluations {
                    return Err(invalid(
                        "dry assembly candidate evaluation budget exhausted",
                    ));
                }
                match self.dry_assembly_candidate(
                    samplers,
                    &mut step,
                    trial,
                    &mut evaluations,
                    config.maximum_evaluations,
                ) {
                    Ok(candidate)
                        if residual_ratio(&candidate, trial, &step, config)? < residual =>
                    {
                        accepted = Some((trial, candidate));
                        break;
                    }
                    Err(PartialSaturationError3::Cancelled) => {
                        return Err(PartialSaturationError3::Cancelled);
                    }
                    Ok(_) | Err(_) => {}
                }
                alpha *= 0.5;
            }
            let Some((trial, candidate)) = accepted else {
                return Err(invalid(
                    "dry assembly Newton backtracking did not decrease residual",
                ));
            };
            correction = trial;
            current = candidate;
            iterations += 1;
        }
        // No probes follow acceptance: next and endpoint samplers belong to
        // current. Independent physical accounts still precede the only swap.
        let receipt = self.dry_assembly_receipt(
            samplers,
            &step,
            &current,
            correction,
            fluid_before,
            slip_before,
            config,
            iterations,
            evaluations,
            condition,
        )?;
        checkpoint(gate)?;
        let mut next = self
            .next
            .take()
            .ok_or_else(|| invalid("dry assembly streamed candidate is unavailable"))?;
        core::mem::swap(&mut self.f, &mut next);
        self.next = Some(next);
        self.topology_locked = true;
        *bodies = current.bodies;
        Ok((
            DryLubricatedPlanesExchange3 {
                motion: receipt,
                planes: current.planes,
            },
            current.surfaces,
        ))
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn dry_assembly_receipt<const N: usize>(
        &self,
        samplers: &[BodyFrameSampler3; N],
        step: &AssemblyStep<'_, N>,
        current: &AssemblyCandidate<N>,
        correction: [[f64; 6]; N],
        fluid_before: ImmersedMomentum3,
        slip_before: [f64; N],
        config: DryBodyCouplingConfig3,
        iterations: usize,
        evaluations: usize,
        condition: Option<f64>,
    ) -> Result<DryBodiesExchange3<N>> {
        let forcing = step
            .interior_forcing
            .ok_or_else(|| invalid("missing dry assembly midpoint interior sources"))?;
        let native = current.native.coupling;
        let film = native.lubrication;
        let plane_reaction: [f64; 6] = core::array::from_fn(|i| {
            if i < 3 {
                film.plane_impulse[i]
            } else {
                film.plane_angular_impulse[i - 3]
            }
        });
        let boundary = core::array::from_fn(|i| {
            if i < 3 {
                native.walls.fluid_population_impulse[i] + native.open.impulse[i]
            } else {
                native.walls.fluid_population_angular_impulse[i - 3]
                    + native.periodic.angular_momentum_change[i - 3]
                    + native.open.angular_impulse[i - 3]
            }
        });
        let mut represented_forcing = core::array::from_fn(|i| {
            if i < 3 {
                native.forcing.impulse[i]
            } else {
                native.forcing.angular_impulse[i - 3]
            }
        });
        let mut receipt = DryBodiesExchange3 {
            native: current.native,
            fluid: [fluid_before, ImmersedMomentum3::default()],
            interior: core::array::from_fn(|a| [step.inside_before[a], current.interior[a]]),
            interior_slip_rms: core::array::from_fn(|a| [slip_before[a], current.slip[a]]),
            interior_correction_impulses: [[0.0; 6]; N],
            interior_forcing_impulses: forcing,
            represented_fluid_forcing_impulse: [0.0; 6],
            physical_external_impulses: [[0.0; 6]; N],
            corrected_hydrodynamic_impulses: [[0.0; 6]; N],
            body_momentum_changes: [[0.0; 6]; N],
            body_momentum_residuals: [[0.0; 6]; N],
            combined_momentum_residual: [0.0; 6],
            physical_external_work: [0.0; N],
            interior_correction_work: [0.0; N],
            represented_kinetic_energy_change: 0.0,
            iterations,
            candidate_evaluations: evaluations,
            jacobian_condition: condition,
        };
        let account_radius = step.radii.into_iter().fold(0.0, f64::max);
        let mut account_scale = norm(normalize(boundary, account_radius))
            + norm(normalize(plane_reaction, account_radius));
        for a in 0..N {
            let (all, inside, slip) = self.dry_inventory(
                true,
                a == 0,
                &samplers[a],
                current.bodies[a],
                step.origin,
                step.gate,
            )?;
            if inside != current.interior[a] || slip != current.slip[a] {
                return Err(invalid(
                    "dry assembly final inventory changed after convergence",
                ));
            }
            if a == 0 {
                receipt.fluid[1] = all;
                receipt.represented_kinetic_energy_change =
                    all.kinetic_energy - fluid_before.kinetic_energy;
            }
            let inside_delta = difference(inside.momentum, step.inside_before[a].momentum);
            let correction_center = denormalize(correction[a], step.radii[a]);
            let correction_origin = at_origin(correction_center, step.arms[a]);
            let physical_external = at_origin(step.external[a], step.arms[a]);
            let before_angular = total_angular_momentum(step.initial[a], step.origin)?;
            let after_angular = total_angular_momentum(current.bodies[a], step.origin)?;
            let body_delta = core::array::from_fn(|i| {
                if i < 3 {
                    step.initial[a].mass
                        * (current.bodies[a].velocity[i] - step.initial[a].velocity[i])
                } else {
                    after_angular[i - 3] - before_angular[i - 3]
                }
            });
            let response = native.bodies[a];
            let hydro = core::array::from_fn(|i| {
                let reaction = if i < 3 {
                    response.hydrodynamic.body_impulse[i]
                } else {
                    response.hydrodynamic.body_angular_impulse[i - 3]
                };
                reaction + inside_delta[i] - forcing[a][i]
            });
            let film_impulse: [f64; 6] = core::array::from_fn(|i| {
                if i < 3 {
                    film.body_impulses[a][i]
                } else {
                    film.body_angular_impulses[a][i - 3]
                }
            });
            let body_residual = core::array::from_fn(|i| {
                body_delta[i] - physical_external[i] - hydro[i] - film_impulse[i]
            });
            let body_scale = [
                body_delta,
                physical_external,
                hydro,
                inside_delta,
                forcing[a],
                film_impulse,
            ]
            .into_iter()
            .map(|x| norm(normalize(x, step.radii[a])))
            .fold(0.0, f64::max);
            if norm(normalize(body_residual, step.radii[a])) > config.tolerance(body_scale)? {
                return Err(invalid(
                    "independent dry assembly body momentum account did not close",
                ));
            }
            let work = |impulse: [f64; 6]| {
                dot(impulse[..3].try_into().unwrap(), response.velocity_midpoint)
                    + dot(
                        impulse[3..].try_into().unwrap(),
                        response.angular_velocity_midpoint,
                    )
            };
            receipt.interior_correction_impulses[a] = correction_origin;
            receipt.physical_external_impulses[a] = physical_external;
            receipt.corrected_hydrodynamic_impulses[a] = hydro;
            receipt.body_momentum_changes[a] = body_delta;
            receipt.body_momentum_residuals[a] = body_residual;
            receipt.physical_external_work[a] = work(step.external[a]);
            receipt.interior_correction_work[a] = work(correction_center);
            receipt.represented_kinetic_energy_change += current.native.bodies[a]
                .kinetic_energy_change
                - (inside.kinetic_energy - step.inside_before[a].kinetic_energy);
            for i in 0..6 {
                receipt.combined_momentum_residual[i] +=
                    body_delta[i] - inside_delta[i] - physical_external[i];
                represented_forcing[i] -= forcing[a][i];
            }
            account_scale += [body_delta, inside_delta, physical_external, forcing[a]]
                .into_iter()
                .map(|x| norm(normalize(x, account_radius)))
                .sum::<f64>();
            if ![
                receipt.physical_external_work[a],
                receipt.interior_correction_work[a],
            ]
            .into_iter()
            .chain(correction_origin)
            .chain(physical_external)
            .chain(hydro)
            .chain(body_delta)
            .chain(body_residual)
            .all(f64::is_finite)
            {
                return Err(invalid("dry assembly body momentum/work overflowed"));
            }
        }
        let fluid_delta = difference(receipt.fluid[1].momentum, fluid_before.momentum);
        for i in 0..6 {
            // The plane receipt is the wrench RECEIVED by every film
            // surface. Its opposite is a source for this open subsystem;
            // internal sphere-pair film impulses must cancel independently.
            receipt.combined_momentum_residual[i] +=
                fluid_delta[i] - boundary[i] - represented_forcing[i] + plane_reaction[i];
        }
        receipt.represented_fluid_forcing_impulse = represented_forcing;
        account_scale += norm(normalize(fluid_delta, account_radius))
            + norm(normalize(represented_forcing, account_radius));
        // Sum the per-body absolute allowances; physical sources and the one
        // fluid inventory still enter the independent common account once.
        let tolerance =
            config.tolerance(account_scale)? + (N - 1) as f64 * config.absolute_impulse_tolerance;
        if ![receipt.represented_kinetic_energy_change, tolerance]
            .into_iter()
            .chain(represented_forcing)
            .chain(receipt.combined_momentum_residual)
            .all(f64::is_finite)
        {
            return Err(invalid("dry assembly combined momentum/work overflowed"));
        }
        if norm(normalize(
            receipt.combined_momentum_residual,
            account_radius,
        )) > tolerance
        {
            return Err(invalid(
                "independent dry assembly/fluid momentum account did not close",
            ));
        }
        Ok(receipt)
    }

    fn dry_assembly_candidate<const N: usize>(
        &mut self,
        samplers: &mut [BodyFrameSampler3; N],
        step: &mut AssemblyStep<'_, N>,
        correction: [[f64; 6]; N],
        evaluations: &mut usize,
        maximum_evaluations: usize,
    ) -> Result<AssemblyCandidate<N>> {
        checkpoint(step.gate)?;
        if *evaluations == maximum_evaluations {
            return Err(invalid(
                "dry assembly candidate evaluation budget exhausted",
            ));
        }
        *evaluations += 1;
        let external = core::array::from_fn(|a| {
            let impulse = denormalize(correction[a], step.radii[a]);
            core::array::from_fn(|i| step.external[a][i] + impulse[i])
        });
        let mut bodies = step.initial;
        let (coupled, surfaces) = self.moving_partial_saturation_resistance_transition(
            samplers,
            step.shapes,
            &mut bodies,
            external,
            step.origin,
            step.gate,
            step.lubrication,
            false,
            step.surfaces,
        )?;
        let native = coupled.motion;
        for body in &native.bodies {
            let volumes = body.sampled_solid_volume;
            if !(volumes[0] > 0.0)
                || volumes.into_iter().any(|v| {
                    !(v > 0.0)
                        || ((v - volumes[0]) / volumes[0]).abs()
                            > step.maximum_relative_volume_change
                })
            {
                return Err(invalid(
                    "dry assembly sampled volume change exceeds geometric admission",
                ));
            }
        }
        // Native sampling is still at the collision midpoint. It depends only
        // on incoming bodies, so this source is constant across ALL probes.
        let forcing = match step.interior_forcing {
            Some(forcing) => forcing,
            None => {
                let mut forcing = [[0.0; 6]; N];
                for a in 0..N {
                    forcing[a] = self.dry_interior_forcing(&samplers[a], step.origin, step.gate)?;
                }
                step.interior_forcing = Some(forcing);
                forcing
            }
        };
        let mut interior = [ImmersedMomentum3::default(); N];
        let mut slip = [0.0; N];
        let mut target = [[0.0; 6]; N];
        let mut residual = [[0.0; 6]; N];
        for a in 0..N {
            samplers[a].sample(
                bodies[a].center,
                bodies[a].orientation,
                &ShapeRef(step.shapes[a]),
                step.gate,
            )?;
            let (_, inside, speed) =
                self.dry_inventory(true, false, &samplers[a], bodies[a], step.origin, step.gate)?;
            interior[a] = inside;
            slip[a] = speed;
            let delta = difference(
                difference(inside.momentum, step.inside_before[a].momentum),
                forcing[a],
            );
            let orbital = cross(step.arms[a], delta[..3].try_into().unwrap());
            target[a] = normalize(
                core::array::from_fn(|i| {
                    if i < 3 {
                        delta[i]
                    } else {
                        delta[i] - orbital[i - 3]
                    }
                }),
                step.radii[a],
            );
            residual[a] = difference(correction[a], target[a]);
            if !target[a].into_iter().chain(residual[a]).all(f64::is_finite) {
                return Err(invalid("dry assembly interior correction overflowed"));
            }
        }
        admit_inventory_occupancy(samplers, step.gate)?;
        Ok(AssemblyCandidate {
            bodies,
            native,
            planes: coupled.planes,
            surfaces,
            interior,
            slip,
            residual,
            target,
        })
    }
}

fn residual_ratio<const N: usize>(
    candidate: &AssemblyCandidate<N>,
    correction: [[f64; 6]; N],
    step: &AssemblyStep<'_, N>,
    config: DryBodyCouplingConfig3,
) -> Result<f64> {
    let mut ratio: f64 = 0.0;
    for a in 0..N {
        let tolerance = config.tolerance(norm(correction[a]).max(norm(candidate.target[a])))?;
        // A center-linear defect also carries angular impulse about the fixed
        // reporting origin. Require that actual moment account to converge;
        // accepting only COM residuals can hide moment-arm amplification.
        let world = normalize(
            at_origin(
                denormalize(candidate.residual[a], step.radii[a]),
                step.arms[a],
            ),
            step.radii[a],
        );
        if !world.into_iter().all(f64::is_finite) {
            return Err(invalid(
                "dry assembly fixed-origin correction residual overflowed",
            ));
        }
        ratio = ratio.max(norm(candidate.residual[a]).max(norm(world)) / tolerance);
    }
    Ok(ratio)
}

fn admit_inventory_occupancy<const N: usize>(
    samplers: &[BodyFrameSampler3; N],
    gate: &CancelGate,
) -> Result<()> {
    for i in 0..samplers[0].fractions.len() {
        if i.is_multiple_of(TILE_CELLS) {
            checkpoint(gate)?;
        }
        let total: f64 = samplers.iter().map(|s| s.fractions[i]).sum();
        if !(total.is_finite() && total <= 1.0) {
            return Err(invalid(
                "dry assembly sampled interiors exceed cell ownership",
            ));
        }
    }
    Ok(())
}

// Called only for N == 1 after exact single-body publication. Preserve its
// physical account bit for bit and give shared walls/periodic/forcing one owner.
fn singleton<const N: usize>(single: DryBodyExchange3) -> DryBodiesExchange3<N> {
    let coupling = single.native.coupling;
    let mut body = single.native;
    body.coupling.walls = Default::default();
    body.coupling.periodic = Default::default();
    body.coupling.open = Default::default();
    body.coupling.forcing = Default::default();
    DryBodiesExchange3 {
        native: MovingPartialSaturationBodiesExchange3 {
            coupling: PartialSaturationBodiesExchange3 {
                bodies: [coupling.into(); N],
                walls: coupling.walls,
                periodic: coupling.periodic,
                open: coupling.open,
                forcing: coupling.forcing,
                lubrication: SphereLubricationExchange3::default(),
                population_mass_residual: 0.0,
                population_impulse_residual: [0.0; 3],
                population_angular_impulse_residual: [0.0; 3],
            },
            bodies: [body; N],
        },
        fluid: single.fluid,
        interior: [single.interior; N],
        interior_slip_rms: [single.interior_slip_rms; N],
        interior_correction_impulses: [single.interior_correction_impulse; N],
        interior_forcing_impulses: [single.interior_forcing_impulse; N],
        represented_fluid_forcing_impulse: single.represented_fluid_forcing_impulse,
        physical_external_impulses: [single.physical_external_impulse; N],
        corrected_hydrodynamic_impulses: [single.corrected_hydrodynamic_impulse; N],
        body_momentum_changes: [single.body_momentum_change; N],
        body_momentum_residuals: [single.body_momentum_residual; N],
        combined_momentum_residual: single.combined_momentum_residual,
        physical_external_work: [single.physical_external_work; N],
        interior_correction_work: [single.interior_correction_work; N],
        represented_kinetic_energy_change: single.represented_kinetic_energy_change,
        iterations: single.iterations,
        candidate_evaluations: single.candidate_evaluations,
        jacobian_condition: single.jacobian_condition,
    }
}

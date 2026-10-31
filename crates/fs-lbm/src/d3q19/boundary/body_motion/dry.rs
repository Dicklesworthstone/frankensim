//! Positive dry inertia with an implicit fictitious-interior momentum account.
//!
//! The native NT collision still retains populations inside the body. For an
//! ordinary immersed solid, its measured reaction must therefore be corrected
//! by the CHANGE of that retained momentum, not by an assumed rigid-fluid mass.
//! We solve J_c = Pi_inside(new(J_c)) - Pi_inside(old) using the SAME native
//! drift/collision/stream/drift candidate, with the original positive dry mass
//! and tensor. Angular inventories use one fixed world origin.
//!
//! The fictitious-fluid correction follows the momentum accounting discussed
//! by Breugem, JCP 231 (2012), doi:10.1016/j.jcp.2012.02.026, and Lacis, Taira
//! & Bagheri, JCP 305 (2016), doi:10.1016/j.jcp.2015.10.041. This NT realization
//! is not their projection method or a claim of their density-ratio stability.
//! Sampled geometry is piecewise constant; rejected/ill-conditioned solves are
//! possible. Interior slip and sampled volume changes remain explicit. Neither
//! exact no-slip, a GCL, nor a moving-free-surface model is asserted.

use super::*;
use crate::d3q19::{E3, Q3};

/// Numerical admission for one implicit dry-body/retained-fluid correction.
/// All impulses are lattice units. Angular residuals are divided by the shape's
/// enclosing radius before comparison with the linear-impulse tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DryBodyCouplingConfig3 {
    /// Maximum accepted Newton corrections, in 1..=32.
    pub maximum_iterations: usize,
    /// Total full native candidate evaluations, including probes/backtracking,
    /// in 1..=512. This is a hard work bound, not merely a convergence count.
    pub maximum_evaluations: usize,
    /// Maximum halvings after a failed Newton trial, in 0..=16.
    pub maximum_backtracks: usize,
    /// Positive absolute normalized impulse tolerance.
    pub absolute_impulse_tolerance: f64,
    /// Nonnegative relative tolerance, at most 1e-3.
    pub relative_tolerance: f64,
    /// Maximum relative change of raw sampled volume from the initial pose,
    /// including midpoint and endpoint, in [0, 0.1]. This bounds geometric
    /// leakage per step; it is not a geometric-conservation certificate.
    pub maximum_relative_volume_change: f64,
    /// Relative forward-difference step, in [1e-8, 1e-2].
    pub jacobian_step: f64,
    /// Maximum estimated 1-norm condition number of the scaled 6x6 Jacobian.
    pub maximum_condition_number: f64,
}

impl Default for DryBodyCouplingConfig3 {
    fn default() -> Self {
        Self {
            maximum_iterations: 12,
            maximum_evaluations: 128,
            maximum_backtracks: 8,
            absolute_impulse_tolerance: 1e-11,
            relative_tolerance: 1e-9,
            maximum_relative_volume_change: 0.01,
            jacobian_step: 1e-5,
            maximum_condition_number: 1e10,
        }
    }
}

impl DryBodyCouplingConfig3 {
    fn admit(self) -> Result<()> {
        if !(1..=32).contains(&self.maximum_iterations)
            || !(1..=512).contains(&self.maximum_evaluations)
            || self.maximum_backtracks > 16
            || !(self.absolute_impulse_tolerance.is_finite()
                && self.absolute_impulse_tolerance > 0.0)
            || !(self.relative_tolerance.is_finite()
                && (0.0..=1e-3).contains(&self.relative_tolerance))
            || !(self.maximum_relative_volume_change.is_finite()
                && (0.0..=0.1).contains(&self.maximum_relative_volume_change))
            || !(self.jacobian_step.is_finite() && (1e-8..=1e-2).contains(&self.jacobian_step))
            || !(self.maximum_condition_number.is_finite() && self.maximum_condition_number >= 1.0)
        {
            return Err(invalid("invalid bounded dry-body correction controls"));
        }
        Ok(())
    }

    fn tolerance(self, scale: f64) -> Result<f64> {
        let value = self.absolute_impulse_tolerance + self.relative_tolerance * scale;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("dry-body correction tolerance overflowed"))
        }
    }
}

/// Direct sum over actual population moments, with no rigid-fluid substitution.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ImmersedMomentum3 {
    /// Population mass, weighted by occupancy for an interior inventory.
    pub mass: f64,
    /// Linear momentum, followed by angular momentum about the requested fixed
    /// world origin. The interior inventory uses the endpoint occupancy.
    pub momentum: [f64; 6],
    /// Macroscopic kinetic energy sum rho |u|^2/2. This excludes thermal,
    /// acoustic/compressional and non-equilibrium population energy.
    pub kinetic_energy: f64,
}

/// Accepted dry-body motion and independently measured population inventories.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct DryBodyExchange3 {
    /// Unchanged native candidate account. Its external input includes the
    /// solved interior correction; use physical_external_impulse for authored
    /// loads. Native inertia is the ORIGINAL dry body in this API.
    pub native: MovingPartialSaturationExchange3,
    /// Full-grid inventories before and after actual streaming.
    pub fluid: [ImmersedMomentum3; 2],
    /// Actual endpoint-fraction-weighted fictitious-fluid inventories.
    pub interior: [ImmersedMomentum3; 2],
    /// Interior RMS fluid-minus-rigid-point speed before/after. Observational:
    /// a converged momentum solve does not make this slip vanish.
    pub interior_slip_rms: [f64; 2],
    /// Solved correction impulse, linear then angular about the fixed origin.
    pub interior_correction_impulse: [f64; 6],
    /// Authored external impulse, linear then angular about the fixed origin.
    pub physical_external_impulse: [f64; 6],
    /// Native hydrodynamic reaction plus the MEASURED interior momentum change.
    /// This is the represented corrected force, not an exact no-slip traction.
    pub corrected_hydrodynamic_impulse: [f64; 6],
    /// Actual endpoint dry-body linear and fixed-origin angular momentum change.
    pub body_momentum_change: [f64; 6],
    /// Dry-body change minus physical external and corrected hydrodynamic input.
    pub body_momentum_residual: [f64; 6],
    /// Independently summed all-fluid minus interior plus dry-body change,
    /// minus physical external, fixed-wall and periodic transport terms.
    pub combined_momentum_residual: [f64; 6],
    /// Work of the authored load at the native midpoint velocity.
    pub physical_external_work: f64,
    /// Work of the solved interior impulse at that same midpoint. It is not
    /// heat or a prescribed external power source.
    pub interior_correction_work: f64,
    /// Dry-body plus all-fluid minus interior MACROSCOPIC kinetic change.
    /// No dissipation or complete thermodynamic-energy closure is inferred.
    pub represented_kinetic_energy_change: f64,
    /// Accepted Newton corrections.
    pub iterations: usize,
    /// Full candidates, including every Jacobian probe and backtrack.
    pub candidate_evaluations: usize,
    /// Last estimated scaled Jacobian condition number; none for initial closure.
    pub jacobian_condition: Option<f64>,
}

struct DryStep<'a, S> {
    shape: &'a S,
    initial: MovingPartialSaturationBody3,
    external: [f64; 6],
    origin: [f64; 3],
    arm: [f64; 3],
    radius: f64,
    inside_before: ImmersedMomentum3,
    maximum_relative_volume_change: f64,
    gate: &'a CancelGate,
}

struct Candidate {
    body: MovingPartialSaturationBody3,
    native: MovingPartialSaturationExchange3,
    interior: ImmersedMomentum3,
    slip: f64,
    residual: [f64; 6],
    target: [f64; 6],
}

impl BoundaryGrid3 {
    /// Advance ONE ordinary dry body and resolved fluid without subtracting any
    /// inertia. Positive dry mass and full positive tensor are retained for
    /// light, neutral and heavy bodies; there is no density-ratio shortcut.
    ///
    /// External WORLD impulse contains linear components at the midpoint COM,
    /// then a center couple. Use REDUCED-pressure loads: e.g. (M-rho*V)g for
    /// gravity once, with no separately forced hydrostatic fluid column.
    /// This bounded lane admits no uniform or per-cell fluid forcing, only
    /// stationary walls/periodic faces and the existing BGK/low-Mach envelope.
    ///
    /// Six unknown correction impulses are solved by safeguarded Newton using
    /// exact streamed candidate inventories and bounded finite-difference
    /// probes. The native populations, shape sampling, full tensor rotation
    /// and split operators are reused. Only the existing post/next/sampler
    /// scratch changes during iteration; no full-grid backup is made.
    ///
    /// # Errors
    /// Ill-conditioned, nonconvergent, over-budget or cancelled steps publish
    /// neither populations, pose/momenta nor topology admission. Independent
    /// linear AND fixed-origin angular momentum accounts must close before the
    /// final publication. Finite interior slip and quadrature volume changes
    /// remain diagnostics, not a hidden exact-rigidity/GCL guarantee.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_dry_partial_saturation_body(
        &mut self,
        sampler: &mut BodyFrameSampler3,
        shape: &impl BodyFrameShape3,
        body: &mut MovingPartialSaturationBody3,
        external_impulse: [f64; 6],
        moment_origin: [f64; 3],
        config: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryBodyExchange3> {
        checkpoint(gate)?;
        config.admit()?;
        if sampler.dims != self.dimensions()
            || self.force != [0.0; 3]
            || self.force_field.is_some()
            || !(body.mass.is_finite() && body.mass > 0.0)
            || !external_impulse
                .into_iter()
                .chain(moment_origin)
                .all(f64::is_finite)
        {
            return Err(invalid(
                "dry-body correction requires matching sampler, positive dry mass and unforced fluid",
            ));
        }
        let initial = *body;
        let radius = shape_radius(shape)?;
        let midpoint: [f64; 3] =
            core::array::from_fn(|i| initial.center[i] + 0.5 * initial.velocity[i]);
        let arm = core::array::from_fn(|i| midpoint[i] - moment_origin[i]);
        sampler.sample(initial.center, initial.orientation, shape, gate)?;
        let (fluid_before, inside_before, slip_before) = self.dry_inventory(
            false,
            true,
            sampler.fractions(),
            initial,
            moment_origin,
            gate,
        )?;
        let step = DryStep {
            shape,
            initial,
            external: external_impulse,
            origin: moment_origin,
            arm,
            radius,
            inside_before,
            maximum_relative_volume_change: config.maximum_relative_volume_change,
            gate,
        };
        if !(inside_before.mass > 0.0) {
            return Err(invalid("dry body has no resolved interior volume"));
        }
        // Scale rotational unknowns by the enclosing radius. Jacobian columns
        // therefore compare impulses with the same units, including off-axis
        // full-tensor response. The floor corresponds to a 0.01 speed probe.
        let probe_scale = (initial.mass + inside_before.mass) * 0.01;
        if !(probe_scale.is_finite() && probe_scale > 0.0) {
            return Err(invalid("dry-body Jacobian scale overflowed"));
        }
        let mut correction = [0.0; 6];
        let mut evaluations = 0;
        let mut iterations = 0;
        let mut condition = None;
        let mut current = self.dry_candidate(
            sampler,
            &step,
            correction,
            &mut evaluations,
            config.maximum_evaluations,
        )?;
        loop {
            checkpoint(gate)?;
            let residual = norm(current.residual);
            let tolerance = config.tolerance(norm(correction).max(norm(current.target)))?;
            if residual <= tolerance {
                break;
            }
            if iterations == config.maximum_iterations {
                return Err(invalid("dry-body correction did not converge"));
            }
            let mut jacobian = [0.0; 36];
            for column in 0..6 {
                let delta = config.jacobian_step * (probe_scale + correction[column].abs());
                if !(delta.is_finite() && delta > 0.0) {
                    return Err(invalid("dry-body Jacobian perturbation overflowed"));
                }
                let mut probe = correction;
                probe[column] += delta;
                let candidate = self.dry_candidate(
                    sampler,
                    &step,
                    probe,
                    &mut evaluations,
                    config.maximum_evaluations,
                )?;
                for row in 0..6 {
                    jacobian[6 * row + column] =
                        (candidate.residual[row] - current.residual[row]) / delta;
                }
            }
            if !jacobian.into_iter().all(f64::is_finite) {
                return Err(invalid("dry-body Jacobian overflowed"));
            }
            let factor = fs_la::factor::lu(&jacobian, 6)
                .map_err(|_| invalid("dry-body correction Jacobian is singular"))?;
            let estimate = factor.condition_1(&jacobian);
            if !(estimate.is_finite() && estimate <= config.maximum_condition_number) {
                return Err(invalid("dry-body correction Jacobian is ill-conditioned"));
            }
            condition = Some(estimate);
            let mut update = current.residual.map(|x| -x);
            factor.solve(&mut update);
            if !update.into_iter().all(f64::is_finite) {
                return Err(invalid("dry-body Newton update overflowed"));
            }
            let mut accepted = None;
            let mut alpha = 1.0;
            for _ in 0..=config.maximum_backtracks {
                let trial = core::array::from_fn(|i| correction[i] + alpha * update[i]);
                // Candidate admission failures can be resolved by a shorter
                // step. Cancellation and work-budget exhaustion cannot.
                if evaluations == config.maximum_evaluations {
                    return Err(invalid("dry-body candidate evaluation budget exhausted"));
                }
                match self.dry_candidate(
                    sampler,
                    &step,
                    trial,
                    &mut evaluations,
                    config.maximum_evaluations,
                ) {
                    Ok(candidate) if norm(candidate.residual) < residual => {
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
                    "dry-body Newton backtracking did not decrease residual",
                ));
            };
            correction = trial;
            current = candidate;
            iterations += 1;
        }

        // The last successful evaluation owns next. No Jacobian/backtrack is
        // performed after this acceptance; every fallible account precedes swap.
        let (fluid_after, inside_after, slip_after) = self.dry_inventory(
            true,
            true,
            sampler.fractions(),
            current.body,
            moment_origin,
            gate,
        )?;
        if inside_after != current.interior || slip_after != current.slip {
            return Err(invalid(
                "dry-body final streamed inventory changed after convergence",
            ));
        }
        let inside_delta = difference(inside_after.momentum, inside_before.momentum);
        let fluid_delta = difference(fluid_after.momentum, fluid_before.momentum);
        let correction_center = denormalize(correction, radius);
        let correction_origin = at_origin(correction_center, arm);
        let physical_external = at_origin(external_impulse, arm);
        let before_angular = total_angular_momentum(initial, moment_origin)?;
        let after_angular = total_angular_momentum(current.body, moment_origin)?;
        let body_delta = core::array::from_fn(|i| {
            if i < 3 {
                initial.mass * (current.body.velocity[i] - initial.velocity[i])
            } else {
                after_angular[i - 3] - before_angular[i - 3]
            }
        });
        let hydrodynamic: [f64; 6] = core::array::from_fn(|i| {
            if i < 3 {
                current.native.coupling.hydrodynamic.body_impulse[i] + inside_delta[i]
            } else {
                current.native.coupling.hydrodynamic.body_angular_impulse[i - 3] + inside_delta[i]
            }
        });
        let body_residual =
            core::array::from_fn(|i| body_delta[i] - physical_external[i] - hydrodynamic[i]);
        let boundary: [f64; 6] = core::array::from_fn(|i| {
            if i < 3 {
                current.native.coupling.walls.fluid_population_impulse[i]
            } else {
                current
                    .native
                    .coupling
                    .walls
                    .fluid_population_angular_impulse[i - 3]
                    + current.native.coupling.periodic.angular_momentum_change[i - 3]
            }
        });
        let combined_residual = core::array::from_fn(|i| {
            fluid_delta[i] - inside_delta[i] + body_delta[i] - physical_external[i] - boundary[i]
        });
        let account_scale = [
            body_delta,
            physical_external,
            hydrodynamic,
            inside_delta,
            fluid_delta,
            boundary,
        ]
        .into_iter()
        .map(|x| norm(normalize(x, radius)))
        .fold(0.0, f64::max);
        let account_tolerance = config.tolerance(account_scale)?;
        if norm(normalize(body_residual, radius)) > account_tolerance
            || norm(normalize(combined_residual, radius)) > account_tolerance
        {
            return Err(invalid(
                "independent dry-body/fluid momentum account did not close",
            ));
        }
        let velocity = current.native.coupling.velocity_midpoint;
        let omega = current.native.coupling.angular_velocity_midpoint;
        let work = |impulse: [f64; 6]| {
            dot(impulse[..3].try_into().unwrap(), velocity)
                + dot(impulse[3..].try_into().unwrap(), omega)
        };
        let physical_work = work(external_impulse);
        let correction_work = work(correction_center);
        let kinetic_change = current.native.kinetic_energy_change + fluid_after.kinetic_energy
            - fluid_before.kinetic_energy
            - (inside_after.kinetic_energy - inside_before.kinetic_energy);
        if ![physical_work, correction_work, kinetic_change]
            .into_iter()
            .chain(correction_origin)
            .chain(physical_external)
            .chain(hydrodynamic)
            .chain(body_delta)
            .chain(body_residual)
            .chain(combined_residual)
            .all(f64::is_finite)
        {
            return Err(invalid(
                "dry-body physical momentum/work account overflowed",
            ));
        }
        let receipt = DryBodyExchange3 {
            native: current.native,
            fluid: [fluid_before, fluid_after],
            interior: [inside_before, inside_after],
            interior_slip_rms: [slip_before, slip_after],
            interior_correction_impulse: correction_origin,
            physical_external_impulse: physical_external,
            corrected_hydrodynamic_impulse: hydrodynamic,
            body_momentum_change: body_delta,
            body_momentum_residual: body_residual,
            combined_momentum_residual: combined_residual,
            physical_external_work: physical_work,
            interior_correction_work: correction_work,
            represented_kinetic_energy_change: kinetic_change,
            iterations,
            candidate_evaluations: evaluations,
            jacobian_condition: condition,
        };
        checkpoint(gate)?;
        let mut next = self
            .next
            .take()
            .ok_or_else(|| invalid("dry-body streamed candidate is unavailable"))?;
        core::mem::swap(&mut self.f, &mut next);
        self.next = Some(next);
        self.topology_locked = true;
        *body = current.body;
        Ok(receipt)
    }

    fn dry_candidate(
        &mut self,
        sampler: &mut BodyFrameSampler3,
        step: &DryStep<'_, impl BodyFrameShape3>,
        correction: [f64; 6],
        evaluations: &mut usize,
        maximum_evaluations: usize,
    ) -> Result<Candidate> {
        checkpoint(step.gate)?;
        if *evaluations == maximum_evaluations {
            return Err(invalid("dry-body candidate evaluation budget exhausted"));
        }
        *evaluations += 1;
        let center_impulse = denormalize(correction, step.radius);
        let combined = core::array::from_fn(|i| step.external[i] + center_impulse[i]);
        let mut body = step.initial;
        let native = self.moving_partial_saturation_transition(
            sampler,
            step.shape,
            &mut body,
            combined,
            step.origin,
            step.gate,
            false,
        )?;
        let volumes = native.sampled_solid_volume;
        if !(volumes[0] > 0.0)
            || volumes.into_iter().any(|v| {
                !(v > 0.0)
                    || ((v - volumes[0]) / volumes[0]).abs() > step.maximum_relative_volume_change
            })
        {
            return Err(invalid(
                "dry-body sampled volume change exceeds geometric admission",
            ));
        }
        sampler.sample(body.center, body.orientation, step.shape, step.gate)?;
        let (_, interior, slip) = self.dry_inventory(
            true,
            false,
            sampler.fractions(),
            body,
            step.origin,
            step.gate,
        )?;
        let delta = difference(interior.momentum, step.inside_before.momentum);
        let orbital = cross(step.arm, delta[..3].try_into().unwrap());
        let target = normalize(
            core::array::from_fn(|i| {
                if i < 3 {
                    delta[i]
                } else {
                    delta[i] - orbital[i - 3]
                }
            }),
            step.radius,
        );
        let residual = difference(correction, target);
        if !target.into_iter().chain(residual).all(f64::is_finite) {
            return Err(invalid("dry-body interior correction overflowed"));
        }
        Ok(Candidate {
            body,
            native,
            interior,
            slip,
            residual,
            target,
        })
    }

    // Raw population inventories are independent of the native impulse receipt.
    // Probes read populations ONLY at occupied cells. Full-field inventories
    // run once initially and once after convergence, without an all-grid copy.
    #[allow(clippy::too_many_arguments)]
    fn dry_inventory(
        &self,
        candidate: bool,
        all_fluid: bool,
        fractions: &[f64],
        body: MovingPartialSaturationBody3,
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(ImmersedMomentum3, ImmersedMomentum3, f64)> {
        let fields = if candidate {
            self.next
                .as_ref()
                .ok_or_else(|| invalid("missing dry-body candidate populations"))?
        } else {
            &self.f
        };
        let omega = body.angular_velocity_world()?;
        let mut all = ImmersedMomentum3::default();
        let mut interior = ImmersedMomentum3::default();
        let mut slip_squared_mass = 0.0;
        for tile in 0..fields[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let cell = self.coords(tile, lane);
                let epsilon = fractions[(cell[2] * self.ny + cell[1]) * self.nx + cell[0]];
                if !all_fluid && epsilon == 0.0 {
                    continue;
                }
                let mut rho = 0.0;
                let mut momentum = [0.0; 3];
                for q in 0..Q3 {
                    let population = fields[q][tile].0[lane];
                    rho += population;
                    momentum[0] += f64::from(E3[q].0) * population;
                    momentum[1] += f64::from(E3[q].1) * population;
                    momentum[2] += f64::from(E3[q].2) * population;
                }
                if !(rho.is_finite() && rho > 0.0) || !momentum.into_iter().all(f64::is_finite) {
                    return Err(invalid("invalid dry-body population inventory"));
                }
                let center: [f64; 3] = cell.map(|x| x as f64 + 0.5);
                let offset = core::array::from_fn(|i| center[i] - origin[i]);
                let angular = cross(offset, momentum);
                let kinetic = 0.5 * dot(momentum, momentum) / rho;
                let quantities: [f64; 6] =
                    core::array::from_fn(|i| if i < 3 { momentum[i] } else { angular[i - 3] });
                if all_fluid {
                    all.mass += rho;
                    all.kinetic_energy += kinetic;
                }
                interior.mass += epsilon * rho;
                interior.kinetic_energy += epsilon * kinetic;
                for i in 0..6 {
                    if all_fluid {
                        all.momentum[i] += quantities[i];
                    }
                    interior.momentum[i] += epsilon * quantities[i];
                }
                if epsilon > 0.0 {
                    let r = core::array::from_fn(|i| center[i] - body.center[i]);
                    let rigid = point_velocity(body.velocity, omega, r);
                    let slip = core::array::from_fn(|i| momentum[i] / rho - rigid[i]);
                    slip_squared_mass += epsilon * rho * dot(slip, slip);
                }
            }
        }
        let slip = if interior.mass > 0.0 {
            (slip_squared_mass / interior.mass).sqrt()
        } else {
            0.0
        };
        if ![
            all.mass,
            all.kinetic_energy,
            interior.mass,
            interior.kinetic_energy,
            slip,
        ]
        .into_iter()
        .chain(all.momentum)
        .chain(interior.momentum)
        .all(f64::is_finite)
        {
            return Err(invalid(
                "dry-body population momentum/energy inventory overflowed",
            ));
        }
        Ok((all, interior, slip))
    }
}

fn norm(values: [f64; 6]) -> f64 {
    values.into_iter().map(f64::abs).fold(0.0, f64::max)
}
fn normalize(mut impulse: [f64; 6], radius: f64) -> [f64; 6] {
    for x in &mut impulse[3..] {
        *x /= radius;
    }
    impulse
}
fn denormalize(mut impulse: [f64; 6], radius: f64) -> [f64; 6] {
    for x in &mut impulse[3..] {
        *x *= radius;
    }
    impulse
}
fn at_origin(mut impulse: [f64; 6], arm: [f64; 3]) -> [f64; 6] {
    let orbital = cross(arm, impulse[..3].try_into().unwrap());
    for i in 0..3 {
        impulse[i + 3] += orbital[i];
    }
    impulse
}
fn difference(after: [f64; 6], before: [f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| after[i] - before[i])
}

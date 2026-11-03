//! Physical dry disk inertia with the actual fictitious-fluid momentum removed.
//!
//! The correction is `J_c = Pi_inside(new) - Pi_inside(old) - J_inside,F`.
//! Inventories use the analytic disk fractions and physical Guo momentum, not
//! a rigid-interior mass approximation. Since D2Q9 holds `rho*g + fext` during
//! collision, its physical endpoint source is `(F_old + F_new)/2`; the density
//! after streaming participates in every correction candidate.
//! Fractional angular inventories integrate the actual circle/cell first
//! moment against each cell's physical momentum. Uniform oblique translation
//! therefore transports the disk interior's geometric moment without creating
//! a fictitious spin. Full-cell native collision/streaming torque is unchanged;
//! exterior angular momentum is its exact cut-interior complement.
//!
//! This follows the fictitious-fluid momentum account of Breugem, JCP 231
//! (2012), doi:10.1016/j.jcp.2012.02.026, and Lacis, Taira & Bagheri, JCP 305
//! (2016), doi:10.1016/j.jcp.2015.10.041. It is an NT realization, not those
//! papers' projection method or a claim of their density-ratio stability.

use super::super::{
    Cell, CollisionModel2, Grid, PartialSaturationBody2, PartialSaturationBodyExchange2,
    raw_population_moments,
};
use super::{Bounds, MovingDisk2, cell_moments};
use crate::{E, Q};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

mod lubricated;
pub use lubricated::{DryLubricatedDiskExchange2, LubricationPlane2, PlaneLubricationExchange2};
use lubricated::{FilmTotals2, PreparedPlanes2};

/// Bounded numerical admission for the three coupled correction impulses.
/// Angular impulses are divided by disk radius before tolerance comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DryDiskCouplingConfig2 {
    /// Maximum accepted Newton corrections, in 1..=32.
    pub maximum_iterations: usize,
    /// Hard limit on all native collision evaluations, in 1..=256.
    pub maximum_evaluations: usize,
    /// Maximum trial halvings per Newton correction, in 0..=16.
    pub maximum_backtracks: usize,
    /// Positive absolute normalized impulse tolerance.
    pub absolute_impulse_tolerance: f64,
    /// Relative impulse tolerance, in [0, 1e-3].
    pub relative_tolerance: f64,
    /// Relative forward difference step, in [1e-8, 1e-2].
    pub jacobian_step: f64,
    /// Maximum estimated 1-norm condition number of the scaled Jacobian.
    pub maximum_condition_number: f64,
}

impl Default for DryDiskCouplingConfig2 {
    fn default() -> Self {
        Self {
            maximum_iterations: 12,
            maximum_evaluations: 96,
            maximum_backtracks: 8,
            absolute_impulse_tolerance: 1e-11,
            relative_tolerance: 1e-9,
            jacobian_step: 1e-5,
            maximum_condition_number: 1e10,
        }
    }
}

impl DryDiskCouplingConfig2 {
    fn admit(self) -> Result<(), DryDiskError2> {
        if !(1..=32).contains(&self.maximum_iterations)
            || !(1..=256).contains(&self.maximum_evaluations)
            || self.maximum_backtracks > 16
            || !(self.absolute_impulse_tolerance.is_finite()
                && self.absolute_impulse_tolerance > 0.0)
            || !(self.relative_tolerance.is_finite()
                && (0.0..=1e-3).contains(&self.relative_tolerance))
            || !(self.jacobian_step.is_finite() && (1e-8..=1e-2).contains(&self.jacobian_step))
            || !(self.maximum_condition_number.is_finite() && self.maximum_condition_number >= 1.0)
        {
            return Err(DryDiskError2::Invalid("invalid dry-disk coupling controls"));
        }
        Ok(())
    }

    fn tolerance(self, scale: f64) -> Result<f64, DryDiskError2> {
        let value = self.absolute_impulse_tolerance + self.relative_tolerance * scale;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(DryDiskError2::Invalid("impulse tolerance overflowed"))
        }
    }
}

/// A numerical or geometric refusal leaves accepted disk and fluid unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DryDiskError2 {
    /// Unsupported state, geometry, or unresolved numerical balance.
    Invalid(&'static str),
    /// The hard native-candidate evaluation budget was consumed.
    EvaluationBudget,
    /// Cancellation was requested before publication.
    Cancelled,
}

impl core::fmt::Display for DryDiskError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => f.write_str(reason),
            Self::EvaluationBudget => f.write_str("dry-disk candidate evaluation budget exhausted"),
            Self::Cancelled => f.write_str("dry-disk coupling cancelled"),
        }
    }
}
impl std::error::Error for DryDiskError2 {}

/// Direct population inventory, with endpoint occupancy for an interior sum.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ImmersedMomentum2 {
    /// Population mass; fractional cells carry their actual fluid density.
    pub mass: f64,
    /// Physical `(Px, Py, Lz)` about the requested fixed lattice origin.
    /// Interior `Lz` integrates the analytic cut-cell first moment against the
    /// piecewise-constant physical momentum density in that lattice cell.
    pub momentum: [f64; 3],
    /// Macroscopic `sum rho |u|^2/2`; excludes thermal and compressional energy.
    pub kinetic_energy: f64,
}

/// One accepted physical dry disk and actual native fluid tick.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct DryDiskExchange2 {
    /// Original collision, with authored PLUS correction/film loads and unchanged
    /// physical dry inertia. Native forcing is the held incoming Guo source.
    pub native: PartialSaturationBodyExchange2,
    /// Location of the coupled impulse between the two free half-drifts.
    pub collision_center: [f64; 2],
    /// Full physical Guo inventories before and after native streaming.
    pub fluid: [ImmersedMomentum2; 2],
    /// Actual retained interior inventories with their own endpoint geometry.
    pub interior: [ImmersedMomentum2; 2],
    /// Solved retained-fluid correction, expressed about the requested origin.
    pub interior_correction_impulse: [f64; 3],
    /// Collision-midpoint occupancy times `(F_old + F_new)/2`, with torque
    /// integrated using the same analytic cut-cell first moment.
    pub interior_forcing_impulse: [f64; 3],
    /// Full physical endpoint source minus the fictitious interior source.
    pub represented_fluid_forcing_impulse: [f64; 3],
    /// Authored body load, including its collision-center lever arm.
    pub physical_external_impulse: [f64; 3],
    /// Native reaction plus measured interior change minus interior forcing.
    pub corrected_hydrodynamic_impulse: [f64; 3],
    /// Actual dry-body linear, orbital and spin angular momentum changes.
    pub body_momentum_change: [f64; 3],
    /// Actual body change minus authored load, corrected reaction and any
    /// separately reported lubrication impulse.
    pub body_momentum_residual: [f64; 3],
    /// Stationary wall/exterior bounce impulse received by the fluid.
    pub boundary_fluid_impulse: [f64; 3],
    /// Angular momentum transported by periodic fluid seams.
    pub periodic_angular_impulse: f64,
    /// Exterior fluid plus dry body minus physical sources, boundary fluxes and
    /// any separately reported lubrication impulse received from plane reservoirs.
    pub combined_momentum_residual: [f64; 3],
    /// Authored impulse/couple work at the native midpoint rigid velocity.
    pub physical_external_work: f64,
    /// Correction work, kept separate from authored power and heat.
    pub interior_correction_work: f64,
    /// Exterior fluid plus dry body macroscopic kinetic-energy change.
    pub represented_kinetic_energy_change: f64,
    /// Accepted Newton corrections.
    pub iterations: usize,
    /// All native collisions, including Jacobian probes and rejected trials.
    pub candidate_evaluations: usize,
    /// Last condition estimate of the scaled three-by-three Jacobian.
    pub jacobian_condition: Option<f64>,
}

/// A physical dry disk coupled to the original forced D2Q9/NT operator.
///
/// Mass and polar inertia are the original positive dry values, including for
/// solids lighter than displaced fluid. The native grid retains populations
/// inside the disk; measured interior changes remove their inertia implicitly.
/// Fixed radius, one disk, BGK/Guo, lattice `dx = dt = 1`, and a fluid-only swept
/// box with the original domain halo are required. Contact, periodic body seam
/// crossing, multiple disks and free surfaces need separate owners.
#[derive(Debug, Clone)]
pub struct DryDisk2 {
    disk: MovingDisk2,
    // Only midpoint cut cells need a correction to area * cell_center.
    // Sorted lattice indices permit one linear merge during each inventory.
    midpoint_cut_moments: Vec<(usize, [f64; 2])>,
}

impl DryDisk2 {
    /// Bind physical dry inertia and analytic disk geometry without changing fluid.
    pub fn new(
        grid: &Grid,
        radius: f64,
        body: PartialSaturationBody2,
        angle: f64,
    ) -> Result<Self, DryDiskError2> {
        admit_grid(grid, [grid.nx, grid.ny])?;
        if !(radius.is_finite()
            && radius > 0.0
            && (radius * radius).is_finite()
            && radius * radius > 0.0
            && body.mass.is_finite()
            && body.mass > 0.0
            && body.inertia.is_finite()
            && body.inertia > 0.0)
            || !body
                .center
                .into_iter()
                .chain(body.velocity)
                .chain([body.angular_velocity, angle])
                .all(f64::is_finite)
        {
            return Err(DryDiskError2::Invalid(
                "dry disk needs finite state and positive radius, mass and inertia",
            ));
        }
        admit_drift(grid, radius, body.center, body.center)?;
        Ok(Self {
            disk: MovingDisk2::new(grid, radius, body, angle),
            midpoint_cut_moments: Vec::new(),
        })
    }

    /// Accepted physical dry state; original mass/inertia never change.
    #[must_use]
    pub fn body(&self) -> PartialSaturationBody2 {
        self.disk.body()
    }
    /// Accepted material orientation, counterclockwise radians.
    #[must_use]
    pub fn angle(&self) -> f64 {
        self.disk.angle()
    }
    /// Fixed radius in lattice lengths.
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.disk.radius()
    }
    /// Original analytic circle/unit-cell intersection fraction at accepted pose.
    #[must_use]
    pub fn solid_fraction(&self, x: usize, y: usize) -> f64 {
        self.disk.solid_fraction(x, y)
    }

    /// Advance the original populations and physical dry body atomically.
    ///
    /// Loads are held impulses over one lattice tick, applied at the collision
    /// center plus an independent couple. Both center and fixed-origin angular
    /// correction residuals must converge. Each Jacobian uses three actual
    /// native collision/stream probes; no second fluid grid is allocated.
    /// Streaming is observed through the native pull and publishes only after
    /// all endpoint, cancellation and independent momentum checks succeed.
    ///
    /// Cancellation is polled at inventory rows and candidate boundaries; one
    /// existing native collision is indivisible. A refusal may replace private
    /// geometry/collision scratch, which the next attempt rebuilds.
    ///
    /// # Panics
    /// The existing D2Q9 collision retains its assertion API for low-Mach,
    /// numerical inertia and population admissibility, including trial states.
    /// Such a panic also precedes accepted disk/fluid publication. No global
    /// temporal-order, exact no-slip, calibrated drag or total-energy claim is made.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        grid: &mut Grid,
        external_impulse: [f64; 2],
        external_couple: f64,
        moment_origin: [f64; 2],
        config: DryDiskCouplingConfig2,
        cancel: &CancelGate,
    ) -> Result<DryDiskExchange2, DryDiskError2> {
        self.step_inner(
            grid,
            external_impulse,
            external_couple,
            moment_origin,
            config,
            cancel,
            None,
        )
        .map(|receipt| receipt.motion)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_inner(
        &mut self,
        grid: &mut Grid,
        external_impulse: [f64; 2],
        external_couple: f64,
        moment_origin: [f64; 2],
        config: DryDiskCouplingConfig2,
        cancel: &CancelGate,
        planes: Option<&PreparedPlanes2>,
    ) -> Result<DryLubricatedDiskExchange2, DryDiskError2> {
        checkpoint(cancel)?;
        config.admit()?;
        admit_grid(grid, self.disk.dims)?;
        if !external_impulse
            .into_iter()
            .chain([external_couple])
            .chain(moment_origin)
            .all(f64::is_finite)
        {
            return Err(DryDiskError2::Invalid(
                "dry-disk load and moment origin must be finite",
            ));
        }
        let initial = self.disk.body;
        let center = core::array::from_fn(|a| initial.center[a] + 0.5 * initial.velocity[a]);
        admit_drift(grid, self.disk.radius, initial.center, center)?;
        self.disk.rasterize(center);
        self.cache_midpoint_moments(center, cancel)?;
        let before = self.measure(grid, None, initial.center, moment_origin, cancel)?;
        let context = Context {
            initial,
            center,
            origin: moment_origin,
            external: [external_impulse[0], external_impulse[1], external_couple],
            before,
            planes,
        };
        let arm = subtract(center, moment_origin);
        let radius = self.disk.radius;
        let mut evaluations = 0;
        let mut correction = [0.0; 3];
        let mut candidate =
            self.evaluate(grid, correction, &context, config, &mut evaluations, cancel)?;
        let mut iterations = 0;
        let mut condition = None;
        loop {
            let scale = norm(correction)
                .max(norm(candidate.target))
                .max(norm(normalize(context.external, radius)));
            let tolerance = config.tolerance(scale)?;
            let global_residual = about_origin(denormalize(candidate.residual, radius), arm);
            if norm(candidate.residual) <= tolerance
                && norm(normalize(global_residual, radius)) <= tolerance
            {
                break;
            }
            if iterations == config.maximum_iterations {
                return Err(DryDiskError2::Invalid(
                    "dry-disk correction did not converge",
                ));
            }
            let mut jacobian = [0.0; 9];
            for column in 0..3 {
                let h = config.jacobian_step * correction[column].abs().max(1.0);
                let mut probe = correction;
                probe[column] += h;
                let actual_h = probe[column] - correction[column];
                let value =
                    self.evaluate(grid, probe, &context, config, &mut evaluations, cancel)?;
                for row in 0..3 {
                    jacobian[3 * row + column] =
                        (value.residual[row] - candidate.residual[row]) / actual_h;
                }
            }
            if !jacobian.into_iter().all(f64::is_finite) {
                return Err(DryDiskError2::Invalid(
                    "dry-disk correction Jacobian overflowed",
                ));
            }
            let factor = fs_la::factor::lu(&jacobian, 3)
                .map_err(|_| DryDiskError2::Invalid("singular dry-disk correction"))?;
            let estimate = factor.condition_1(&jacobian);
            if !estimate.is_finite() || estimate > config.maximum_condition_number {
                return Err(DryDiskError2::Invalid(
                    "ill-conditioned dry-disk correction",
                ));
            }
            condition = Some(estimate);
            let mut update = candidate.residual.map(|x| -x);
            factor.solve(&mut update);
            if !update.into_iter().all(f64::is_finite) {
                return Err(DryDiskError2::Invalid(
                    "dry-disk correction update overflowed",
                ));
            }
            let merit = norm(candidate.residual).max(norm(normalize(global_residual, radius)));
            let mut accepted = None;
            let mut fraction = 1.0;
            for _ in 0..=config.maximum_backtracks {
                let trial = core::array::from_fn(|k| correction[k] + fraction * update[k]);
                match self.evaluate(grid, trial, &context, config, &mut evaluations, cancel) {
                    Ok(value) => {
                        let global = about_origin(denormalize(value.residual, radius), arm);
                        if norm(value.residual).max(norm(normalize(global, radius))) < merit {
                            accepted = Some((trial, value));
                            break;
                        }
                    }
                    Err(DryDiskError2::Invalid(_)) => {}
                    Err(error) => return Err(error),
                }
                fraction *= 0.5;
            }
            let Some((trial, value)) = accepted else {
                return Err(DryDiskError2::Invalid(
                    "dry-disk correction backtracking failed",
                ));
            };
            correction = trial;
            candidate = value;
            iterations += 1;
        }
        // The last successful evaluation is precisely this accepted residual:
        // no Jacobian probe or rejected trial may leave the published post.
        let after = &candidate.measure;
        let physical_external = about_origin(context.external, arm);
        let mut correction_center = denormalize(correction, radius);
        if candidate.films.active != 0 {
            correction_center[0] -= candidate.films.impulse[0];
            correction_center[1] -= candidate.films.impulse[1];
            if candidate.films.center_couple != 0.0 {
                correction_center[2] -= candidate.films.center_couple;
            }
        }
        let correction_origin = about_origin(correction_center, arm);
        let measured_correction: [f64; 3] = core::array::from_fn(|k| {
            after.interior.momentum[k] - before.interior.momentum[k] - after.interior_force[k]
        });
        let hydro = candidate.exchange.hydrodynamic;
        let native_reaction = about_origin(
            [
                hydro.solid_impulse[0],
                hydro.solid_impulse[1],
                hydro.solid_angular_impulse,
            ],
            arm,
        );
        let corrected_hydrodynamic =
            core::array::from_fn(|k| native_reaction[k] + measured_correction[k]);
        let body_change = body_increment(initial, candidate.body, moment_origin);
        let mut body_residual = core::array::from_fn(|k| {
            body_change[k] - physical_external[k] - corrected_hydrodynamic[k]
        });
        let represented_force = core::array::from_fn(|k| after.force[k] - after.interior_force[k]);
        let mut combined = core::array::from_fn(|k| {
            after.fluid.momentum[k]
                - before.fluid.momentum[k]
                - (after.interior.momentum[k] - before.interior.momentum[k])
                + body_change[k]
                - physical_external[k]
                - represented_force[k]
                - after.boundary[k]
                - if k == 2 { after.periodic } else { 0.0 }
        });
        if candidate.films.active != 0 {
            for k in 0..3 {
                body_residual[k] -= candidate.films.impulse[k];
                combined[k] -= candidate.films.impulse[k];
            }
        }
        let scale = norm(normalize(body_change, radius))
            .max(norm(normalize(physical_external, radius)))
            .max(norm(normalize(corrected_hydrodynamic, radius)))
            .max(norm(normalize(candidate.films.impulse, radius)));
        let body_allowance = config.tolerance(scale)? + 4096.0 * f64::EPSILON * scale.max(1.0);
        let fluid_scale = norm(normalize(before.fluid.momentum, radius))
            .max(norm(normalize(after.fluid.momentum, radius)))
            .max(norm(normalize(before.interior.momentum, radius)))
            .max(norm(normalize(after.interior.momentum, radius)))
            .max(norm(normalize(after.force, radius)))
            .max(norm(normalize(after.boundary, radius)))
            .max((after.periodic / radius).abs());
        let combined_allowance = body_allowance + 4096.0 * f64::EPSILON * fluid_scale.max(1.0);
        if !body_residual
            .into_iter()
            .chain(combined)
            .all(f64::is_finite)
            || norm(normalize(body_residual, radius)) > body_allowance
            || norm(normalize(combined, radius)) > combined_allowance
        {
            return Err(DryDiskError2::Invalid(
                "dry-disk physical momentum balance is unresolved",
            ));
        }
        let work = |impulse: [f64; 3]| {
            impulse[0].mul_add(
                candidate.exchange.velocity_midpoint[0],
                impulse[1] * candidate.exchange.velocity_midpoint[1],
            ) + impulse[2] * candidate.exchange.angular_velocity_midpoint
        };
        let physical_work = work(context.external);
        let correction_work = work(correction_center);
        let kinetic = candidate.exchange.kinetic_energy_change + after.fluid.kinetic_energy
            - before.fluid.kinetic_energy
            - (after.interior.kinetic_energy - before.interior.kinetic_energy);
        if ![physical_work, correction_work, kinetic]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(DryDiskError2::Invalid("dry-disk work account overflowed"));
        }
        let receipt = DryDiskExchange2 {
            native: candidate.exchange,
            collision_center: center,
            fluid: [before.fluid, after.fluid],
            interior: [before.interior, after.interior],
            interior_correction_impulse: correction_origin,
            interior_forcing_impulse: after.interior_force,
            represented_fluid_forcing_impulse: represented_force,
            physical_external_impulse: physical_external,
            corrected_hydrodynamic_impulse: corrected_hydrodynamic,
            body_momentum_change: body_change,
            body_momentum_residual: body_residual,
            boundary_fluid_impulse: after.boundary,
            periodic_angular_impulse: after.periodic,
            combined_momentum_residual: combined,
            physical_external_work: physical_work,
            interior_correction_work: correction_work,
            represented_kinetic_energy_change: kinetic,
            iterations,
            candidate_evaluations: evaluations,
            jacobian_condition: condition,
        };
        let plane_receipts = if let Some(planes) = planes {
            planes.admit_endpoint(candidate.body.center, radius, cancel)?;
            planes.receipts(
                candidate.exchange.velocity_midpoint,
                candidate.exchange.angular_velocity_midpoint,
                moment_origin,
                cancel,
            )?
        } else {
            Vec::new()
        };
        let receipt = DryLubricatedDiskExchange2 {
            motion: receipt,
            planes: plane_receipts,
            body_impulse: candidate.films.impulse,
            body_work: candidate.films.body_work,
            plane_work: candidate.films.plane_work,
            dissipated_energy: candidate.films.dissipation,
        };
        checkpoint(cancel)?;
        grid.stream_from(&self.disk.post);
        self.disk.body = candidate.body;
        self.disk.angle = candidate.angle;
        Ok(receipt)
    }

    fn evaluate(
        &mut self,
        grid: &Grid,
        correction: [f64; 3],
        context: &Context<'_>,
        config: DryDiskCouplingConfig2,
        evaluations: &mut usize,
        cancel: &CancelGate,
    ) -> Result<Candidate, DryDiskError2> {
        checkpoint(cancel)?;
        if *evaluations == config.maximum_evaluations {
            return Err(DryDiskError2::EvaluationBudget);
        }
        let load = denormalize(correction, self.disk.radius);
        let load = core::array::from_fn::<_, 3, _>(|k| load[k] + context.external[k]);
        if !load.into_iter().all(f64::is_finite) {
            return Err(DryDiskError2::Invalid("dry-disk trial load overflowed"));
        }
        let mut body = context.initial;
        body.center = context.center;
        *evaluations += 1;
        let exchange = grid.collide_into_with_partial_saturation_body(
            &mut self.disk.post,
            &self.disk.fractions,
            &mut body,
            [load[0], load[1]],
            load[2],
        );
        body.center = core::array::from_fn(|k| context.center[k] + 0.5 * body.velocity[k]);
        let angle =
            self.disk.angle + 0.5 * context.initial.angular_velocity + 0.5 * body.angular_velocity;
        if !angle.is_finite() {
            return Err(DryDiskError2::Invalid("dry-disk angle overflowed"));
        }
        admit_drift(grid, self.disk.radius, context.center, body.center)?;
        let measure = self.measure(
            grid,
            Some(&self.disk.post),
            body.center,
            context.origin,
            cancel,
        )?;
        let target_origin = core::array::from_fn::<_, 3, _>(|k| {
            (Dd::from_f64(measure.interior.momentum[k])
                - Dd::from_f64(context.before.interior.momentum[k])
                - Dd::from_f64(measure.interior_force[k]))
            .to_f64()
        });
        let mut target = normalize(
            about_origin(target_origin, subtract(context.origin, context.center)),
            self.disk.radius,
        );
        let films = if let Some(planes) = context.planes {
            planes.evaluate(
                exchange.velocity_midpoint,
                exchange.angular_velocity_midpoint,
                context.origin,
                cancel,
            )?
        } else {
            FilmTotals2::default()
        };
        if films.active != 0 {
            // The solve uses a center couple scaled by radius, while the
            // separate film receipt keeps its full fixed-origin moment.
            target[0] += films.impulse[0];
            target[1] += films.impulse[1];
            if films.center_couple != 0.0 {
                target[2] += films.center_couple / self.disk.radius;
            }
        }
        let residual = core::array::from_fn(|k| correction[k] - target[k]);
        if !target.into_iter().chain(residual).all(f64::is_finite) {
            return Err(DryDiskError2::Invalid(
                "dry-disk interior correction overflowed",
            ));
        }
        Ok(Candidate {
            body,
            angle,
            exchange,
            measure,
            target,
            residual,
            films,
        })
    }

    fn cache_midpoint_moments(
        &mut self,
        center: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<(), DryDiskError2> {
        self.midpoint_cut_moments.clear();
        let bounds = self
            .disk
            .raster_bounds
            .expect("midpoint raster precedes its moments");
        for y in bounds.lower[1]..bounds.upper[1] {
            checkpoint(cancel)?;
            for x in bounds.lower[0]..bounds.upper[0] {
                let i = y * self.disk.dims[0] + x;
                let fraction = self.disk.fractions[i];
                if fraction > 0.0 && fraction < 1.0 {
                    let moments = cell_moments(
                        self.disk.radius,
                        [x as f64 - center[0], y as f64 - center[1]],
                    );
                    self.midpoint_cut_moments.push((i, moments.first_moment));
                }
            }
        }
        Ok(())
    }

    fn measure(
        &self,
        grid: &Grid,
        post: Option<&[[f64; Q]]>,
        center: [f64; 2],
        origin: [f64; 2],
        cancel: &CancelGate,
    ) -> Result<Measurement, DryDiskError2> {
        let mut fluid = Inventory::default();
        let mut interior = Inventory::default();
        let mut force = [Dd::ZERO; 3];
        let mut interior_force = [Dd::ZERO; 3];
        let mut boundary = [Dd::ZERO; 3];
        let mut periodic = Dd::ZERO;
        let mut source_cut = 0;
        for y in 0..grid.ny {
            checkpoint(cancel)?;
            for x in 0..grid.nx {
                let i = grid.idx(x, y);
                if grid.flags[i] != Cell::Fluid {
                    continue;
                }
                let arm = [x as f64 - origin[0], y as f64 - origin[1]];
                let populations = if let Some(post) = post {
                    core::array::from_fn(|q| {
                        let source = grid.source(x, y, q);
                        let value = grid.stationary_pull_population(post, i, q, source);
                        if source.is_none_or(|s| grid.flags[s] == Cell::Wall) {
                            add_impulse(
                                &mut boundary,
                                [
                                    2.0 * value * f64::from(E[q].0),
                                    2.0 * value * f64::from(E[q].1),
                                ],
                                arm,
                                1.0,
                            );
                        } else if let Some(s) = source {
                            let shift = [
                                x as f64 - (s % grid.nx) as f64 - f64::from(E[q].0),
                                y as f64 - (s / grid.nx) as f64 - f64::from(E[q].1),
                            ];
                            periodic = periodic
                                + Dd::from_f64(cross(
                                    shift,
                                    [value * f64::from(E[q].0), value * f64::from(E[q].1)],
                                ));
                        }
                        value
                    })
                } else {
                    grid.f[i]
                };
                let (rho, raw) = raw_population_moments(&populations);
                let source =
                    core::array::from_fn::<_, 2, _>(|k| grid.g[k].mul_add(rho, grid.fext[i][k]));
                let momentum = core::array::from_fn::<_, 2, _>(|k| raw[k] + 0.5 * source[k]);
                let kinetic = 0.5 * (momentum[0] * momentum[0] + momentum[1] * momentum[1]) / rho;
                if !(rho.is_finite() && rho > 0.0 && kinetic.is_finite())
                    || !populations
                        .into_iter()
                        .chain(momentum)
                        .chain(source)
                        .all(f64::is_finite)
                {
                    return Err(DryDiskError2::Invalid(
                        "dry-disk candidate has invalid physical fluid moments",
                    ));
                }
                let moments = cell_moments(
                    self.disk.radius,
                    [x as f64 - center[0], y as f64 - center[1]],
                );
                fluid.add(rho, momentum, kinetic, arm, 1.0);
                interior.add(rho, momentum, kinetic, arm, moments.area);
                interior.add_first_moment(moments.first_moment, momentum);
                if post.is_some() {
                    let old_rho = raw_population_moments(&grid.f[i]).0;
                    let average = core::array::from_fn(|k| {
                        0.5 * grid.g[k].mul_add(old_rho, grid.fext[i][k]) + 0.5 * source[k]
                    });
                    add_impulse(&mut force, average, arm, 1.0);
                    add_impulse(&mut interior_force, average, arm, self.disk.fractions[i]);
                    if let Some(&(index, first_moment)) = self.midpoint_cut_moments.get(source_cut)
                        && index == i
                    {
                        interior_force[2] =
                            interior_force[2] + Dd::from_f64(cross(first_moment, average));
                        source_cut += 1;
                    }
                }
            }
        }
        let result = Measurement {
            fluid: fluid.finish(),
            interior: interior.finish(),
            force: force.map(Dd::to_f64),
            interior_force: interior_force.map(Dd::to_f64),
            boundary: boundary.map(Dd::to_f64),
            periodic: periodic.to_f64(),
        };
        if !result
            .fluid
            .momentum
            .into_iter()
            .chain(result.interior.momentum)
            .chain(result.force)
            .chain(result.interior_force)
            .chain(result.boundary)
            .chain([
                result.fluid.mass,
                result.interior.mass,
                result.fluid.kinetic_energy,
                result.interior.kinetic_energy,
                result.periodic,
            ])
            .all(f64::is_finite)
        {
            return Err(DryDiskError2::Invalid(
                "dry-disk physical inventory overflowed",
            ));
        }
        Ok(result)
    }
}

#[derive(Clone, Copy)]
struct Context<'a> {
    initial: PartialSaturationBody2,
    center: [f64; 2],
    origin: [f64; 2],
    external: [f64; 3],
    before: Measurement,
    planes: Option<&'a PreparedPlanes2>,
}
struct Candidate {
    body: PartialSaturationBody2,
    angle: f64,
    exchange: PartialSaturationBodyExchange2,
    measure: Measurement,
    target: [f64; 3],
    residual: [f64; 3],
    films: FilmTotals2,
}
#[derive(Clone, Copy)]
struct Measurement {
    fluid: ImmersedMomentum2,
    interior: ImmersedMomentum2,
    force: [f64; 3],
    interior_force: [f64; 3],
    boundary: [f64; 3],
    periodic: f64,
}
#[derive(Default)]
struct Inventory([Dd; 5]);
impl Inventory {
    fn add(&mut self, rho: f64, momentum: [f64; 2], kinetic: f64, arm: [f64; 2], fraction: f64) {
        let values = [rho, momentum[0], momentum[1], cross(arm, momentum), kinetic];
        for (sum, value) in self.0.iter_mut().zip(values) {
            *sum = *sum + Dd::from_f64(fraction * value);
        }
    }
    fn finish(self) -> ImmersedMomentum2 {
        let x = self.0.map(Dd::to_f64);
        ImmersedMomentum2 {
            mass: x[0],
            momentum: [x[1], x[2], x[3]],
            kinetic_energy: x[4],
        }
    }
    fn add_first_moment(&mut self, first_moment: [f64; 2], momentum: [f64; 2]) {
        if first_moment != [0.0; 2] {
            self.0[3] = self.0[3] + Dd::from_f64(cross(first_moment, momentum));
        }
    }
}
fn add_impulse(sum: &mut [Dd; 3], impulse: [f64; 2], arm: [f64; 2], fraction: f64) {
    for (sum, value) in sum
        .iter_mut()
        .zip([impulse[0], impulse[1], cross(arm, impulse)])
    {
        *sum = *sum + Dd::from_f64(fraction * value);
    }
}
fn checkpoint(cancel: &CancelGate) -> Result<(), DryDiskError2> {
    if cancel.is_requested() {
        Err(DryDiskError2::Cancelled)
    } else {
        Ok(())
    }
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}
fn subtract(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn about_origin(impulse: [f64; 3], arm: [f64; 2]) -> [f64; 3] {
    [
        impulse[0],
        impulse[1],
        impulse[2] + cross(arm, [impulse[0], impulse[1]]),
    ]
}
fn normalize(x: [f64; 3], radius: f64) -> [f64; 3] {
    [x[0], x[1], x[2] / radius]
}
fn denormalize(x: [f64; 3], radius: f64) -> [f64; 3] {
    [x[0], x[1], x[2] * radius]
}
fn norm(x: [f64; 3]) -> f64 {
    x.into_iter().map(f64::abs).fold(0.0, f64::max)
}
fn body_increment(
    before: PartialSaturationBody2,
    after: PartialSaturationBody2,
    origin: [f64; 2],
) -> [f64; 3] {
    let delta = core::array::from_fn(|k| before.mass * (after.velocity[k] - before.velocity[k]));
    let old_momentum = before.velocity.map(|v| before.mass * v);
    [
        delta[0],
        delta[1],
        before.inertia * (after.angular_velocity - before.angular_velocity)
            + cross(subtract(after.center, origin), delta)
            + cross(subtract(after.center, before.center), old_momentum),
    ]
}
fn admit_grid(grid: &Grid, dims: [usize; 2]) -> Result<(), DryDiskError2> {
    let Some(count) = grid.nx.checked_mul(grid.ny) else {
        return Err(DryDiskError2::Invalid("dry-disk grid size overflowed"));
    };
    if grid.nx < 3
        || grid.ny < 3
        || [grid.nx, grid.ny] != dims
        || grid.flags.len() != count
        || grid.f.len() != count
        || grid.tau.len() != count
        || grid.fext.len() != count
        || grid.collision != CollisionModel2::Bgk
        || !grid.g.into_iter().all(f64::is_finite)
        || !grid
            .flags
            .iter()
            .all(|c| matches!(c, Cell::Fluid | Cell::Wall))
    {
        return Err(DryDiskError2::Invalid(
            "dry disk requires unchanged BGK fluid/wall grid dimensions and fields",
        ));
    }
    Ok(())
}
fn admit_drift(
    grid: &Grid,
    radius: f64,
    from: [f64; 2],
    to: [f64; 2],
) -> Result<(), DryDiskError2> {
    let dims = [grid.nx, grid.ny];
    for axis in 0..2 {
        if !(from[axis].is_finite()
            && to[axis].is_finite()
            && from[axis].min(to[axis]) - radius >= 0.5
            && from[axis].max(to[axis]) + radius <= dims[axis] as f64 - 1.5)
        {
            return Err(DryDiskError2::Invalid(
                "physical dry disk left the interior lattice domain",
            ));
        }
    }
    let bounds = Bounds {
        lower: core::array::from_fn(|k| (from[k].min(to[k]) - radius + 0.5).floor() as usize),
        upper: core::array::from_fn(|k| (from[k].max(to[k]) + radius + 0.5).ceil() as usize),
    };
    for y in bounds.lower[1]..bounds.upper[1] {
        for x in bounds.lower[0]..bounds.upper[0] {
            if grid.flags[grid.idx(x, y)] != Cell::Fluid {
                return Err(DryDiskError2::Invalid(
                    "physical dry-disk sweep requires fluid cells; contact is not modeled",
                ));
            }
        }
    }
    Ok(())
}

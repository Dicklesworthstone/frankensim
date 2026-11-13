//! Work-conjugate finite-body motion through the original directional contact.
//!
//! Midpoint compliance is `h²/(2 m_reduced)`. The native ellipse solves both
//! tangent components with its fixed-size scalar projection. Only sources that
//! own speed/path require another scalar consistency root; no vector motion
//! Newton solve is introduced. Static reaction is used only at accepted rest.
//! Moving intervals use the kinetic ellipse even below the larger static limit.
//!
//! Arrest probes close both impulse components and physical heat. Turning is
//! subdivided before accepting mirror velocities or crossing a force-free
//! material-axis velocity through zero. These are midpoint paths, not exact
//! arbitrary anisotropic trajectories. Refine the maximum step independently
//! of residual tolerances; nonlinear weakening need not have a unique root.

use super::*;
use crate::tribothermal::dynamics::{
    MotionDrive,
    planar::{PlanarMotionSolve, root, root_side},
};
use fs_mbd_material::ResolvedBody;
use fs_thermal::enthalpy::DepositedCell;
use fs_tribo::FrictionRegime;

#[derive(Clone)]
pub struct DirectionalPair<'a> {
    contact: DirectionalSlidingContact<'a>,
    inventory: [DepositedCell; 2],
    material_bodies: Option<Arc<[ResolvedBody; 2]>>,
    velocities_mps: [[f64; 2]; 2],
}

pub struct DirectionalMotionSegment<S> {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub relative_velocity_before_mps: [f64; 2],
    pub relative_velocity_after_mps: [f64; 2],
    pub reached_stop: bool,
    pub sliding: DirectionalSlidingReport<S>,
    pub body_displacements_m: [[f64; 2]; 2],
    pub body_momentum_residual_n_s: [[f64; 2]; 2],
    /// Reduced-mass impulse equivalent of source versus actual path increment.
    pub source_motion_residual_n_s: f64,
    pub kinetic_energy_change_j: f64,
    pub external_work_j: f64,
    pub total_energy_residual_j: f64,
}

pub struct DirectionalMotionReport<S> {
    pub segments: Vec<DirectionalMotionSegment<S>>,
    pub motion_trials: usize,
    pub friction_trials: usize,
    pub halvings: usize,
    pub body_displacements_m: [[f64; 2]; 2],
    pub body_momentum_residual_n_s: [[f64; 2]; 2],
    pub kinetic_energy_change_j: f64,
    pub thermal_energy_change_j: f64,
    pub deposited_heat_j: f64,
    pub external_work_j: f64,
    pub total_energy_residual_j: f64,
}

pub struct DirectionalMotionCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    pair: DirectionalPair<'a>,
    report: DirectionalMotionReport<S>,
}
impl<'a, 't, S> DirectionalMotionCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn pair(&self) -> &DirectionalPair<'a> {
        &self.pair
    }
    pub const fn report(&self) -> &DirectionalMotionReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        DirectionalPair<'a>,
        DirectionalMotionReport<S>,
    ) {
        (self.thermal, self.pair, self.report)
    }
}

fn norm<E>(v: [f64; 2]) -> Result<f64, SlidingError<E>> {
    number(det::hypot(v[0], v[1]))
}
fn dot<E>(a: [f64; 2], b: [f64; 2]) -> Result<f64, SlidingError<E>> {
    number(number(a[0] * b[0])? + number(a[1] * b[1])?)
}
fn principal<E>(axis: [f64; 2], v: [f64; 2]) -> Result<[f64; 2], SlidingError<E>> {
    Ok([dot(axis, v)?, dot([-axis[1], axis[0]], v)?])
}

fn projection_roundoff<E>(axis: [f64; 2], v: [f64; 2]) -> Result<[f64; 2], SlidingError<E>> {
    let a = axis.map(f64::abs);
    let v = v.map(f64::abs);
    Ok([dot(a, v)?, dot([a[1], a[0]], v)?].map(|x| 64.0 * f64::EPSILON * x))
}

/// A small transverse impulse cannot select a discontinuously different static
/// branch above kinetic capacity, or erase an unforced frictionless velocity.
fn may_arrest<E>(
    before: [f64; 2],
    force: [f64; 2],
    normal: f64,
    friction: &TangentialStep<SourceThermalEllipticState>,
) -> Result<bool, SlidingError<E>> {
    let axis = friction.state.constitutive_state().principal_axis();
    let v = principal(axis, before)?;
    let f = principal(axis, force)?;
    let caps = friction
        .state
        .constitutive_law()
        .kinetic_coefficients()
        .map(|mu| mu * normal);
    let mut scaled = [0.0; 2];
    let mut inside = true;
    for i in 0..2 {
        number::<E>(caps[i])?;
        if caps[i] == 0.0 {
            if f[i] == 0.0 && v[i] != 0.0 {
                return Ok(false);
            }
            if f[i] != 0.0 {
                inside = false;
            }
        } else {
            scaled[i] = f[i] / caps[i];
        }
    }
    inside &= det::hypot(scaled[0], scaled[1]) < 1.0;
    if inside {
        return Ok(true);
    }
    let net = [
        number(force[0] - friction.traction_n[0])?,
        number(force[1] - friction.traction_n[1])?,
    ];
    let speed = norm(before)?;
    let effort = norm(net)?;
    if speed == 0.0 || effort == 0.0 {
        return Ok(false);
    }
    let cross = (before[0] / speed) * (net[1] / effort) - (before[1] / speed) * (net[0] / effort);
    Ok(cross.abs() <= 64.0 * f64::EPSILON)
}

fn admissible_turn<E>(
    before: [f64; 2],
    after: [f64; 2],
    force: [f64; 2],
    axis: [f64; 2],
) -> Result<bool, SlidingError<E>> {
    if dot(before, after)? < 0.0 {
        return Ok(false);
    }
    let speed_before = norm(before)?;
    let speed_after = norm(after)?;
    if speed_before == 0.0 || speed_after == 0.0 {
        return Ok(true);
    }
    let b = principal(axis, before)?;
    let a = principal(axis, after)?;
    let rb = projection_roundoff(axis, before)?;
    let ra = projection_roundoff(axis, after)?;
    let f = principal(axis, force)?;
    for i in 0..2 {
        let crossed = b[i] > 0.0 && a[i] < 0.0 || b[i] < 0.0 && a[i] > 0.0;
        let against_force = f[i] == 0.0 || b[i] > 0.0 && f[i] > 0.0 || b[i] < 0.0 && f[i] < 0.0;
        if crossed && against_force && b[i].abs() > rb[i] && a[i].abs() > ra[i] {
            return Ok(false);
        }
    }
    Ok(true)
}

impl<'a> DirectionalPair<'a> {
    /// Retain complete finite caloric masses as tangent inertias.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        contact: DirectionalSlidingContact<'a>,
        velocities_mps: [[f64; 2]; 2],
    ) -> Result<Self, SlidingError> {
        if contact.time_s != thermal.time_s()
            || velocities_mps.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(SlidingError::Invalid(
                "directional bodies require synchronized clocks and finite velocities",
            ));
        }
        let cell = |i: usize| {
            thermal
                .cells()
                .get(&contact.data.cells[i])
                .cloned()
                .ok_or(SlidingError::Invalid("inactive directional body"))
        };
        let inventory = [cell(0)?, cell(1)?];
        let pair = Self {
            contact,
            inventory,
            material_bodies: None,
            velocities_mps,
        };
        for i in 0..2 {
            number::<Infallible>(velocities_mps[0][i] - velocities_mps[1][i])?;
        }
        pair.kinetic_energy_j()?;
        Ok(pair)
    }

    /// Preserve original source mass, geometry and density receipts. Complete
    /// body mass/volume must match the ordered caloric inventory to roundoff.
    pub fn from_material_bodies(
        thermal: &DepositionSimulation<'_>,
        contact: DirectionalSlidingContact<'a>,
        bodies: [ResolvedBody; 2],
        velocities_mps: [[f64; 2]; 2],
    ) -> Result<Self, SlidingError> {
        let mut pair = Self::new(thermal, contact, velocities_mps)?;
        for (i, body) in bodies.iter().enumerate() {
            let identity = pair.contact.data.point.materials[i].content_hash();
            if body
                .regions()
                .iter()
                .any(|r| r.material().card_identity() != identity)
            {
                return Err(SlidingError::Invalid(
                    "directional body must retain its ordered original material",
                ));
            }
            let volume = body.regions().iter().try_fold(0.0, |sum, region| {
                number::<Infallible>(sum + region.geometry().volume_m3())
            })?;
            let cell = &pair.inventory[i];
            for (actual, expected) in [
                (body.inertia().mass(), cell.mass_kg),
                (volume, cell.geometry.volume_m3),
            ] {
                if (actual - expected).abs()
                    > 128.0 * f64::EPSILON * actual.abs().max(expected.abs())
                {
                    return Err(SlidingError::Invalid(
                        "directional source mass and volume must match its complete caloric cell",
                    ));
                }
            }
        }
        pair.material_bodies = Some(Arc::new(bodies));
        pair.kinetic_energy_j()?;
        Ok(pair)
    }
    pub const fn contact(&self) -> &DirectionalSlidingContact<'a> {
        &self.contact
    }
    pub const fn velocities_mps(&self) -> [[f64; 2]; 2] {
        self.velocities_mps
    }
    pub fn material_bodies(&self) -> Option<&[ResolvedBody; 2]> {
        self.material_bodies.as_deref()
    }
    pub fn masses_kg(&self) -> [f64; 2] {
        self.material_bodies.as_deref().map_or_else(
            || self.inventory.each_ref().map(|c| c.mass_kg),
            |bodies| bodies.each_ref().map(|body| body.inertia().mass()),
        )
    }
    pub fn kinetic_energy_j(&self) -> Result<f64, SlidingError> {
        let mut energy = 0.0;
        for (mass, velocity) in self.masses_kg().into_iter().zip(self.velocities_mps) {
            for component in velocity {
                energy = number(energy + number(0.5 * mass * component * component)?)?;
            }
        }
        Ok(energy)
    }

    /// Advance a held-force interval with physical midpoint work and heat.
    /// Normal load, patch geometry and tangent/material frame stay fixed.
    /// Turning and temperature/speed/path-dependent data need time refinement.
    /// A previously stuck segment does not locate weakening inside that segment.
    /// Every source, arithmetic, cancellation and work-budget refusal propagates.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_planar_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: MotionDrive,
        controls: PlanarMotionSolve,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DirectionalMotionCandidate<'a, 't, S>, SlidingError<E>> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() {
                Err(SlidingError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        if controls.maximum_segments == 0
            || controls.maximum_motion_trials == 0
            || !(1..=128).contains(&controls.maximum_bisections)
            || [
                controls.maximum_step_s,
                controls.impulse_tolerance_n_s,
                controls.event_time_tolerance_s,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || !target_s.is_finite()
            || target_s <= thermal.time_s()
            || thermal.time_s() != self.contact.time_s
            || energy_guess_j.len() != thermal.cells().len()
            || energy_guess_j.iter().any(|v| !v.is_finite())
            || drive.body_forces_n.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(SlidingError::Invalid(
                "invalid directional clocks, controls, forces or thermal guess",
            ));
        }
        for previous in &self.inventory {
            if thermal.cells().get(&previous.id) != Some(previous) {
                return Err(SlidingError::Invalid(
                    "directional motion must retain original finite thermal inventory",
                ));
            }
        }
        let masses = self.masses_kg();
        let total_mass = number(masses[0] + masses[1])?;
        let reduced_mass = number((masses[0] / total_mass) * masses[1])?;
        if reduced_mass <= 0.0 {
            return Err(SlidingError::Invalid(
                "positive finite reduced inertia required",
            ));
        }
        let relative_force = [
            number(
                reduced_mass
                    * number(
                        drive.body_forces_n[0][0] / masses[0]
                            - drive.body_forces_n[1][0] / masses[1],
                    )?,
            )?,
            number(
                reduced_mass
                    * number(
                        drive.body_forces_n[0][1] / masses[0]
                            - drive.body_forces_n[1][1] / masses[1],
                    )?,
            )?,
        ];
        let mut pair = self.clone();
        let mut heat = thermal.clone();
        let mut guess = energy_guess_j.to_vec();
        let trials = Cell::new(0usize);
        let friction_trials = Cell::new(0usize);
        let mut report = DirectionalMotionReport {
            segments: Vec::new(),
            motion_trials: 0,
            friction_trials: 0,
            halvings: 0,
            body_displacements_m: [[0.0; 2]; 2],
            body_momentum_residual_n_s: [[0.0; 2]; 2],
            kinetic_energy_change_j: 0.0,
            thermal_energy_change_j: 0.0,
            deposited_heat_j: 0.0,
            external_work_j: 0.0,
            total_energy_residual_j: 0.0,
        };
        while heat.time_s() < target_s {
            check()?;
            if report.segments.len() >= controls.maximum_segments {
                return Err(SlidingError::MotionBudget);
            }
            let start = heat.time_s();
            let remaining = number(target_s - start)?;
            let before = [
                number(pair.velocities_mps[0][0] - pair.velocities_mps[1][0])?,
                number(pair.velocities_mps[0][1] - pair.velocities_mps[1][1])?,
            ];
            let initial_speed = norm(before)?;
            let source_motion = pair.contact.motion_dependent();
            let solve = |duration: f64,
                         increment: [f64; 2],
                         compliance: f64,
                         source_distance: f64,
                         mode| {
                check()?;
                if trials.get() >= controls.maximum_motion_trials {
                    return Err(SlidingError::MotionBudget);
                }
                trials.set(trials.get() + 1);
                let end = if duration == remaining {
                    target_s
                } else {
                    number(start + duration)?
                };
                if end <= start {
                    return Err(SlidingError::Invalid(
                        "directional segment has no representable duration",
                    ));
                }
                let candidate = pair.contact.step_candidate_with_motion(
                    &heat,
                    end,
                    SlidingDrive {
                        normal_load_n: drive.normal_load_n,
                        relative_increment_m: increment,
                    },
                    compliance,
                    source_distance,
                    mode,
                    controls.sliding,
                    &guess,
                    |state| resolve(state),
                    || (*cancel.borrow_mut())(),
                )?;
                friction_trials.set(
                    friction_trials
                        .get()
                        .checked_add(candidate.report.friction_trials)
                        .ok_or(SlidingError::MotionBudget)?,
                );
                Ok(candidate)
            };
            let actual_duration = |h: f64| {
                if h == remaining {
                    Ok(remaining)
                } else {
                    number(number(start + h)? - start)
                }
            };
            let proposed = remaining.min(controls.maximum_step_s);
            let end = number(start + proposed)?;
            let mut duration = actual_duration(if end > start && end.next_up() == target_s {
                remaining
            } else {
                proposed
            })?;
            let (candidate, after, reached_stop) = loop {
                check()?;
                if duration <= 0.0 {
                    return Err(SlidingError::Invalid(
                        "directional timestep is unrepresentable",
                    ));
                }
                if initial_speed > 0.0 {
                    let direction = before.map(|v| v / initial_speed);
                    let at_stop = |h: f64| {
                        let h = actual_duration(h)?;
                        let increment =
                            [number(0.5 * h * before[0])?, number(0.5 * h * before[1])?];
                        let candidate = solve(
                            h,
                            increment,
                            0.0,
                            norm(increment)?,
                            EllipticFrictionMode::Kinetic,
                        )?;
                        let residual = number(
                            -reduced_mass * initial_speed
                                - h * number(
                                    dot(relative_force, direction)?
                                        - dot(candidate.report.friction.traction_n, direction)?,
                                )?,
                        )?;
                        Ok((candidate, residual))
                    };
                    let (stop, stop_residual) = at_stop(duration)?;
                    if stop_residual >= 0.0 {
                        let (stop, _) = root(
                            duration,
                            -reduced_mass * initial_speed,
                            stop,
                            stop_residual,
                            controls,
                            true,
                            at_stop,
                        )?;
                        let h = number(stop.thermal.time_s() - start)?;
                        let traction = stop.report.friction.traction_n;
                        let residual = [
                            number(
                                -reduced_mass * before[0] - h * (relative_force[0] - traction[0]),
                            )?,
                            number(
                                -reduced_mass * before[1] - h * (relative_force[1] - traction[1]),
                            )?,
                        ];
                        if may_arrest(
                            before,
                            relative_force,
                            drive.normal_load_n,
                            &stop.report.friction,
                        )? && norm(residual)? <= 0.5 * controls.impulse_tolerance_n_s
                        {
                            break (stop, [0.0; 2], true);
                        }
                    }
                }
                // A rank-one ellipse can arrest its dissipative component
                // while the orthogonal frictionless component keeps moving.
                // Resolve that event rather than freezing the whole vector or
                // repeatedly bisecting through a component mirror velocity.
                let mu = pair.contact.state.constitutive_law().kinetic_coefficients();
                let active = if mu[0] == 0.0 && mu[1] > 0.0 {
                    Some(1)
                } else if mu[1] == 0.0 && mu[0] > 0.0 {
                    Some(0)
                } else {
                    None
                };
                if let Some(active) =
                    active.filter(|_| initial_speed > 0.0 && drive.normal_load_n > 0.0)
                {
                    let axis = pair.contact.state.constitutive_state().principal_axis();
                    let b = principal(axis, before)?;
                    let f = principal(axis, relative_force)?;
                    // Bound this projection's own roundoff, not unrelated free-axis speed.
                    if b[active].abs() > projection_roundoff(axis, before)?[active] {
                        let free = 1 - active;
                        let free_axis = if free == 0 { axis } else { [-axis[1], axis[0]] };
                        let active_axis = if active == 0 {
                            axis
                        } else {
                            [-axis[1], axis[0]]
                        };
                        let direction = active_axis.map(|x| x * b[active].signum());
                        let endpoint = |h: f64| {
                            let v = number(b[free] + h * f[free] / reduced_mass)?;
                            Ok::<_, SlidingError<E>>([
                                number(free_axis[0] * v)?,
                                number(free_axis[1] * v)?,
                            ])
                        };
                        let at_component_stop = |h: f64| {
                            let h = actual_duration(h)?;
                            let after = endpoint(h)?;
                            let increment = [
                                number(0.5 * h * (before[0] + after[0]))?,
                                number(0.5 * h * (before[1] + after[1]))?,
                            ];
                            let candidate = solve(
                                h,
                                increment,
                                0.0,
                                norm(increment)?,
                                EllipticFrictionMode::Kinetic,
                            )?;
                            let residual = number(
                                -reduced_mass * b[active].abs()
                                    - h * number(
                                        dot(relative_force, direction)?
                                            - dot(candidate.report.friction.traction_n, direction)?,
                                    )?,
                            )?;
                            Ok((candidate, residual))
                        };
                        let (stop, residual) = at_component_stop(duration)?;
                        if residual >= 0.0 {
                            let (stop, _) = root(
                                duration,
                                -reduced_mass * b[active].abs(),
                                stop,
                                residual,
                                controls,
                                true,
                                at_component_stop,
                            )?;
                            let h = number(stop.thermal.time_s() - start)?;
                            let after = endpoint(h)?;
                            let t = stop.report.friction.traction_n;
                            let defect = [
                                number(
                                    reduced_mass * (after[0] - before[0])
                                        - h * (relative_force[0] - t[0]),
                                )?,
                                number(
                                    reduced_mass * (after[1] - before[1])
                                        - h * (relative_force[1] - t[1]),
                                )?,
                            ];
                            let updated = stop
                                .report
                                .friction
                                .state
                                .constitutive_law()
                                .kinetic_coefficients();
                            if updated[free] == 0.0
                                && updated[active] > 0.0
                                && norm(defect)? <= 0.5 * controls.impulse_tolerance_n_s
                            {
                                break (stop, after, after == [0.0; 2]);
                            }
                        }
                    }
                }
                let compliance = number(0.5 * duration * duration / reduced_mass)?;
                if compliance <= 0.0 {
                    return Err(SlidingError::Invalid(
                        "directional midpoint compliance is unrepresentable",
                    ));
                }
                let released = [
                    number(duration * before[0] + compliance * relative_force[0])?,
                    number(duration * before[1] + compliance * relative_force[1])?,
                ];
                if initial_speed == 0.0 {
                    let static_probe = solve(
                        duration,
                        released,
                        compliance,
                        0.0,
                        EllipticFrictionMode::StaticOrSliding,
                    )?;
                    if static_probe.report.friction.regime == FrictionRegime::Sticking
                        || drive.normal_load_n == 0.0
                            && static_probe.report.friction.slip_increment_m == [0.0; 2]
                    {
                        if static_probe.report.friction.slip_increment_m != [0.0; 2] {
                            return Err(SlidingError::Invalid(
                                "directional static probe retained nonzero slip",
                            ));
                        }
                        break (static_probe, [0.0; 2], false);
                    }
                }
                let free_distance = norm(released)?;
                let mut proposal = None;
                if free_distance > 0.0 {
                    let lower = solve(
                        duration,
                        released,
                        compliance,
                        0.0,
                        EllipticFrictionMode::Kinetic,
                    )?;
                    if lower.report.friction.regime != FrictionRegime::Sticking {
                        if source_motion {
                            let scale = number(2.0 * reduced_mass / duration)?;
                            let low =
                                number(-scale * norm(lower.report.friction.slip_increment_m)?)?;
                            drop(lower);
                            let at_distance = |distance: f64| {
                                let candidate = solve(
                                    duration,
                                    released,
                                    compliance,
                                    distance,
                                    EllipticFrictionMode::Kinetic,
                                )?;
                                let residual = number(
                                    scale
                                        * number(
                                            distance
                                                - norm(candidate.report.friction.slip_increment_m)?,
                                        )?,
                                )?;
                                Ok((candidate, residual))
                            };
                            let (upper, high) = at_distance(free_distance)?;
                            proposal = Some(
                                root_side(
                                    free_distance,
                                    low,
                                    upper,
                                    high,
                                    controls,
                                    false,
                                    true,
                                    at_distance,
                                )?
                                .0,
                            );
                        } else {
                            proposal = Some(lower);
                        }
                    }
                }
                if let Some(candidate) = proposal {
                    let slip = candidate.report.friction.slip_increment_m;
                    let after = [
                        number(2.0 * slip[0] / duration - before[0])?,
                        number(2.0 * slip[1] / duration - before[1])?,
                    ];
                    let axis = candidate
                        .report
                        .friction
                        .state
                        .constitutive_state()
                        .principal_axis();
                    if initial_speed == 0.0 || admissible_turn(before, after, relative_force, axis)?
                    {
                        let stopped = after == [0.0; 2];
                        if !stopped
                            || may_arrest(
                                before,
                                relative_force,
                                drive.normal_load_n,
                                &candidate.report.friction,
                            )?
                        {
                            break (candidate, after, stopped);
                        }
                    }
                }
                if report.halvings >= controls.maximum_halvings {
                    return Err(SlidingError::MotionBudget);
                }
                report.halvings += 1;
                let next = actual_duration(0.5 * duration)?;
                if next <= 0.0 || next >= duration {
                    return Err(SlidingError::Invalid(
                        "directional turning subdivision is unrepresentable",
                    ));
                }
                duration = next;
            };
            let duration = number(candidate.thermal.time_s() - start)?;
            let actual_distance = norm(candidate.report.friction.slip_increment_m)?;
            let source_residual = if source_motion {
                number(
                    (2.0 * reduced_mass / duration)
                        * number(candidate.report.source_distance_m - actual_distance)?,
                )?
            } else {
                0.0
            };
            if source_residual.abs() > controls.impulse_tolerance_n_s {
                return Err(SlidingError::Invalid(
                    "directional source speed/path does not match actual solved motion",
                ));
            }
            let mut velocities = [[0.0; 2]; 2];
            let mut displacement = [[0.0; 2]; 2];
            let mut momentum = [[0.0; 2]; 2];
            let mut kinetic = 0.0;
            let mut work = 0.0;
            for axis in 0..2 {
                let initial_momentum = number(
                    number(masses[0] * pair.velocities_mps[0][axis])?
                        + number(masses[1] * pair.velocities_mps[1][axis])?,
                )?;
                let centre = number(
                    number(
                        initial_momentum
                            + number(
                                duration
                                    * number(
                                        drive.body_forces_n[0][axis] + drive.body_forces_n[1][axis],
                                    )?,
                            )?,
                    )? / total_mass,
                )?;
                velocities[0][axis] = number(centre + (masses[1] / total_mass) * after[axis])?;
                velocities[1][axis] = number(centre - (masses[0] / total_mass) * after[axis])?;
                if after[axis] != 0.0 && velocities[0][axis] == velocities[1][axis] {
                    return Err(SlidingError::Invalid(
                        "native directional relative motion is unrepresentable",
                    ));
                }
                for body in 0..2 {
                    let old = pair.velocities_mps[body][axis];
                    let next = velocities[body][axis];
                    let mean = number(0.5 * old + 0.5 * next)?;
                    let dv = number(next - old)?;
                    displacement[body][axis] = number(duration * mean)?;
                    kinetic = number(kinetic + number(masses[body] * dv * mean)?)?;
                    work = number(
                        work + number(drive.body_forces_n[body][axis] * displacement[body][axis])?,
                    )?;
                    let sign = if body == 0 { -1.0 } else { 1.0 };
                    momentum[body][axis] = number(
                        number(masses[body] * dv)?
                            - number(
                                duration
                                    * number(
                                        drive.body_forces_n[body][axis]
                                            + sign * candidate.report.friction.traction_n[axis],
                                    )?,
                            )?,
                    )?;
                }
            }
            for residual in momentum {
                if norm(residual)? > controls.impulse_tolerance_n_s {
                    return Err(SlidingError::Invalid(
                        "directional midpoint body impulse exceeds its independent gate",
                    ));
                }
            }
            let (next_heat, contact, sliding) = candidate.into_parts();
            let residual = number(
                kinetic
                    + sliding.thermal.stored_energy_change_j
                    + sliding.friction.stored_change_j
                    + sliding.partition.other_work_j()
                    - work
                    - sliding.external_heat_j,
            )?;
            if residual.abs() > controls.total_energy_tolerance_j {
                return Err(SlidingError::Invalid(
                    "directional kinetic, finite heat and actual force work do not close",
                ));
            }
            report.kinetic_energy_change_j = number(report.kinetic_energy_change_j + kinetic)?;
            report.thermal_energy_change_j =
                number(report.thermal_energy_change_j + sliding.thermal.stored_energy_change_j)?;
            report.deposited_heat_j = number(report.deposited_heat_j + sliding.deposited_heat_j)?;
            report.external_work_j = number(report.external_work_j + work)?;
            report.total_energy_residual_j = number(report.total_energy_residual_j + residual)?;
            for body in 0..2 {
                for axis in 0..2 {
                    report.body_displacements_m[body][axis] =
                        number(report.body_displacements_m[body][axis] + displacement[body][axis])?;
                    report.body_momentum_residual_n_s[body][axis] = number(
                        report.body_momentum_residual_n_s[body][axis] + momentum[body][axis],
                    )?;
                }
            }
            report.segments.push(DirectionalMotionSegment {
                start_time_s: start,
                end_time_s: next_heat.time_s(),
                relative_velocity_before_mps: before,
                relative_velocity_after_mps: after,
                reached_stop,
                sliding,
                body_displacements_m: displacement,
                body_momentum_residual_n_s: momentum,
                source_motion_residual_n_s: source_residual,
                kinetic_energy_change_j: kinetic,
                external_work_j: work,
                total_energy_residual_j: residual,
            });
            pair.inventory = pair
                .inventory
                .each_ref()
                .map(|old| next_heat.cells()[&old.id].clone());
            pair.contact = contact;
            pair.velocities_mps = velocities;
            guess.fill(0.0);
            heat = next_heat;
        }
        if report.total_energy_residual_j.abs() > controls.total_energy_tolerance_j
            || report
                .body_momentum_residual_n_s
                .iter()
                .any(|r| det::hypot(r[0], r[1]) > controls.impulse_tolerance_n_s)
        {
            return Err(SlidingError::Invalid(
                "directional cumulative physical work or vector impulse exceeds its gate",
            ));
        }
        report.motion_trials = trials.get();
        report.friction_trials = friction_trials.get();
        check()?;
        Ok(DirectionalMotionCandidate {
            thermal: heat,
            pair,
            report,
        })
    }
}

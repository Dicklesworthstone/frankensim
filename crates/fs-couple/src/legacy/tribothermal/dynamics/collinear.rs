//! Physical stopping/reversal through the original material and heat owners.
//!
//! Relative motion occupies one declared axis of the unchanged tangent frame;
//! common transverse motion is retained. Each sliding segment uses
//! `du = h*(v_relative_before + v_relative_after)/2`. The original prescribed-slip
//! contact solves traction and finite endpoint temperatures for that SAME du.
//! A bounded scalar impulse solve determines the endpoint velocity. A segment
//! that would reverse is cut at zero relative velocity, then the original static
//! return map decides whether the remaining interval sticks or starts sliding.
//! Thus stopping heat is traction times actual stopping distance, never a
//! kinetic defect reassigned to heat. Held Coulomb coefficients/forces recover
//! the exact constant-acceleration solution, including the stopping instant.
//!
//! Material coefficients are still queried at endpoint temperature and segment
//! mean slip speed. Thermal transport remains the original backward-Euler model.
//! Source-dependent coefficients require timestep refinement; a sign bracket
//! is not a uniqueness proof or a continuous first-event certificate for an
//! arbitrary nonmonotone material curve. Thermal weakening inside a previously
//! stuck interval is not independently event-localized. Source/topology events
//! must already be cut by the caller. No normal work, rotation or 2-D turning
//! trajectory is inferred. The general planar endpoint path remains available.

use super::*;

/// Budgets and independent physical gates for the collinear event path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CollinearMotionSolve {
    pub sliding: SlidingSolve,
    /// Zero or one, selecting an axis of the contact's existing tangent frame.
    pub tangent_axis: usize,
    /// Accepted sliding/sticking segments in this one atomic call.
    pub maximum_segments: usize,
    /// Complete original contact/thermal solves, including discarded probes.
    pub maximum_motion_trials: usize,
    /// Safeguarded secant/bisection iterations per scalar root, in 1..=128.
    pub maximum_bisections: usize,
    pub impulse_tolerance_n_s: f64,
    /// Final stopping-time bracket width, unless an exact residual root is found.
    pub event_time_tolerance_s: f64,
    pub total_energy_tolerance_j: f64,
}

/// One actual held-effort sliding or sticking interval.
pub struct CollinearMotionSegment<S> {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub relative_velocity_before_mps: f64,
    pub relative_velocity_after_mps: f64,
    /// The interval ended at an explicitly solved stop, including a final stop.
    pub reached_stop: bool,
    /// Original source, traction, slip, heat partition and thermal accounts.
    pub sliding: SlidingReport<S>,
    pub body_displacements_m: [[f64; 2]; 2],
    pub body_momentum_residual_n_s: [[f64; 2]; 2],
    pub kinetic_energy_change_j: f64,
    pub external_work_j: f64,
    /// Kinetic + caloric + interface storage + outside share - external supply.
    pub total_energy_residual_j: f64,
}

pub struct CollinearMotionReport<S> {
    pub segments: Vec<CollinearMotionSegment<S>>,
    pub motion_trials: usize,
    pub friction_trials: usize,
    pub body_displacements_m: [[f64; 2]; 2],
    pub kinetic_energy_change_j: f64,
    pub thermal_energy_change_j: f64,
    /// Actual physical heat deposited in the two finite bodies.
    pub deposited_heat_j: f64,
    pub external_work_j: f64,
    pub total_energy_residual_j: f64,
}

/// No independently publishable mechanical, contact-history or thermal half-step.
pub struct CollinearMotionCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    pair: SlidingPair<'a>,
    report: CollinearMotionReport<S>,
}
impl<'a, 't, S> CollinearMotionCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub const fn pair(&self) -> &SlidingPair<'a> { &self.pair }
    pub const fn report(&self) -> &CollinearMotionReport<S> { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>, SlidingPair<'a>, CollinearMotionReport<S>) {
        (self.thermal, self.pair, self.report)
    }
}

/// A sign-bracketed scalar solve; every trial is a complete original candidate.
/// No state is published on a source error, failed bracket or exhausted budget.
fn bisect<T, E>(
    mut lower: f64, mut upper: f64, mut lower_residual: f64,
    upper_candidate: T, mut upper_residual: f64,
    controls: CollinearMotionSolve, time_root: bool,
    mut evaluate: impl FnMut(f64) -> Result<(T, f64), SlidingError<E>>,
) -> Result<(T, f64), SlidingError<E>> {
    if !lower_residual.is_finite() || !upper_residual.is_finite()
        || lower_residual > 0.0 || upper_residual < 0.0 || lower >= upper {
        return Err(SlidingError::Invalid("collinear impulse root is not sign bracketed"));
    }
    if upper_residual == 0.0 {
        return Ok((upper_candidate, upper));
    }
    drop(upper_candidate);
    for _ in 0..controls.maximum_bisections {
        let width = number(upper - lower)?;
        // Constant-Coulomb impulse is affine: its secant gives the physical
        // root in one trial. Keep bisection progress for nonlinear source laws.
        let secant = lower - lower_residual * width / (upper_residual - lower_residual);
        let middle = if secant.is_finite()
            && secant > lower + width / 16.0 && secant < upper - width / 16.0 {
            secant
        } else { number(lower + 0.5 * width)? };
        if middle <= lower || middle >= upper {
            return Err(SlidingError::Invalid("collinear root interval is unrepresentable"));
        }
        let (candidate, residual) = evaluate(middle)?;
        number::<E>(residual)?;
        if residual == 0.0 || (residual.abs() <= controls.impulse_tolerance_n_s
            && (!time_root || upper - lower <= controls.event_time_tolerance_s)) {
            return Ok((candidate, middle));
        }
        if residual < 0.0 { lower = middle; lower_residual = residual; }
        else { upper = middle; upper_residual = residual; }
    }
    Err(SlidingError::MotionBudget)
}

impl<'a> SlidingPair<'a> {
    /// Advance collinear relative motion, resolving a stop before reversal.
    ///
    /// External forces and normal load are held. Transverse relative velocity
    /// and relative external acceleration must be exactly zero; they are never
    /// projected away. Material curves must cover zero speed and the trial
    /// bracket. Every segment keeps the original caloric laws, body inventory,
    /// source selection and friction history. `resolve` supplies only OTHER
    /// heat/conduction and must not mutate accepted histories.
    ///
    /// `energy_guess_j` contains energy INCREMENTS in sorted thermal-cell order.
    /// Subsequent event segments start from zero increments against their actual
    /// accepted predecessor, preserving narrow material temperature domains.
    /// All segments publish together; even a refusal after the stop preserves
    /// the original caller states. No numerical kinetic-loss heat is introduced.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_collinear_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: MotionDrive,
        controls: CollinearMotionSolve,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CollinearMotionCandidate<'a, 't, S>, SlidingError<E>> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() { Err(SlidingError::Cancelled) } else { Ok(()) }
        };
        check()?;
        if controls.tangent_axis > 1 || controls.maximum_segments == 0
            || controls.maximum_motion_trials == 0
            || !(1..=128).contains(&controls.maximum_bisections)
            || [controls.impulse_tolerance_n_s, controls.event_time_tolerance_s,
                controls.total_energy_tolerance_j].iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !target_s.is_finite() || target_s <= thermal.time_s()
            || thermal.time_s() != self.contact.time_s()
            || energy_guess_j.len() != thermal.cells().len()
            || drive.body_forces_n.iter().flatten().any(|v| !v.is_finite()) {
            return Err(SlidingError::Invalid("invalid collinear motion clocks, forces or solve controls"));
        }
        for old in &self.inventory {
            if thermal.cells().get(&old.id) != Some(old) {
                return Err(SlidingError::Invalid("collinear motion must retain the accepted finite thermal inventory"));
            }
        }
        let axis = controls.tangent_axis;
        let transverse = 1 - axis;
        let masses = self.masses_kg();
        let total_mass = number(masses[0] + masses[1])?;
        let reduced_mass = number((masses[0] / total_mass) * masses[1])?;
        if reduced_mass <= 0.0
            || self.velocities_mps[0][transverse] != self.velocities_mps[1][transverse]
            || drive.body_forces_n[0][transverse] / masses[0]
                != drive.body_forces_n[1][transverse] / masses[1] {
            return Err(SlidingError::Invalid("collinear motion requires positive inertia and zero transverse relative motion/acceleration"));
        }
        let relative_force = number(reduced_mass * number(
            drive.body_forces_n[0][axis] / masses[0] - drive.body_forces_n[1][axis] / masses[1])?)?;
        let mut pair = self.clone();
        let mut heat = thermal.clone();
        let mut guess = energy_guess_j.to_vec();
        let trials = Cell::new(0usize);
        let friction_trials = Cell::new(0usize);
        let mut report = CollinearMotionReport { segments: Vec::new(), motion_trials: 0,
            friction_trials: 0, body_displacements_m: [[0.0; 2]; 2],
            kinetic_energy_change_j: 0.0, thermal_energy_change_j: 0.0,
            deposited_heat_j: 0.0, external_work_j: 0.0, total_energy_residual_j: 0.0 };
        while heat.time_s() < target_s {
            check()?;
            if report.segments.len() >= controls.maximum_segments { return Err(SlidingError::MotionBudget); }
            let start = heat.time_s();
            let remaining = number(target_s - start)?;
            let before = number(pair.velocities_mps[0][axis] - pair.velocities_mps[1][axis])?;
            let solve = |duration: f64, increment: f64, compliance: f64, material_increment: f64| {
                check()?;
                if trials.get() >= controls.maximum_motion_trials { return Err(SlidingError::MotionBudget); }
                trials.set(trials.get() + 1);
                let end = if duration == remaining { target_s } else { number(start + duration)? };
                if end <= start { return Err(SlidingError::Invalid("collinear segment has no representable duration")); }
                let mut released = [0.0; 2]; released[axis] = increment;
                let mut source = [0.0; 2]; source[axis] = material_increment;
                let candidate = pair.contact.step_candidate_with_motion(&heat, end,
                    SlidingDrive { normal_load_n: drive.normal_load_n, relative_increment_m: released },
                    compliance, source, controls.sliding, &guess, |state, _| resolve(state),
                    || (*cancel.borrow_mut())())?;
                friction_trials.set(friction_trials.get().checked_add(candidate.report().friction_trials)
                    .ok_or(SlidingError::MotionBudget)?);
                Ok(candidate)
            };
            let at_velocity = |duration: f64, after: f64| {
                let actual_duration = if duration == remaining { remaining }
                    else { number(number(start + duration)? - start)? };
                let increment = number(actual_duration * number(0.5 * before + 0.5 * after)?)?;
                let candidate = solve(duration, increment, 0.0, increment)?;
                let traction = candidate.report().friction.traction_n[axis];
                let residual = number(number(reduced_mass * number(after - before)?)?
                    - number(actual_duration * number(relative_force - traction)?)?)?;
                Ok((candidate, residual))
            };
            // The static probe is admissible only when the body already rests.
            // Applying it to nonzero before-velocity would erase stopping heat.
            let static_candidate = if before == 0.0 {
                let compliance = number(0.5 * remaining * remaining / reduced_mass)?;
                if compliance <= 0.0 { return Err(SlidingError::Invalid("collinear static compliance is unrepresentable")); }
                let candidate = solve(remaining, number(compliance * relative_force)?, compliance, 0.0)?;
                if candidate.report().friction.regime == FrictionRegime::Sticking {
                    if candidate.report().friction.slip_increment_m != [0.0; 2] {
                        return Err(SlidingError::Invalid("native static probe retained nonzero slip"));
                    }
                    Some(candidate)
                } else { None }
            } else { None };
            let (candidate, after, reached_stop) = if let Some(candidate) = static_candidate {
                (candidate, 0.0, false)
            } else {
                let direction = if before != 0.0 { before.signum() } else { relative_force.signum() };
                if direction == 0.0 {
                    return Err(SlidingError::Invalid("unforced resting contact did not admit native sticking"));
                }
                let (zero, zero_residual) = at_velocity(remaining, 0.0)?;
                let lower_residual = number(direction * zero_residual)?;
                if before != 0.0 && lower_residual >= 0.0 {
                    let (stopping, _) = bisect(0.0, remaining, -reduced_mass * before.abs(),
                        zero, lower_residual, controls, true,
                        |duration| at_velocity(duration, 0.0).and_then(|(candidate, residual)|
                            Ok((candidate, number(direction * residual)?))))?;
                    (stopping, 0.0, true)
                } else {
                    drop(zero);
                    let free = number(direction * number(before + remaining * relative_force / reduced_mass)?)?;
                    if free <= 0.0 { return Err(SlidingError::Invalid("sliding velocity bracket does not retain its direction")); }
                    let (upper, upper_residual) = at_velocity(remaining, direction * free)?;
                    let (sliding, speed) = bisect(0.0, free, lower_residual,
                        upper, number(direction * upper_residual)?, controls, false,
                        |speed| at_velocity(remaining, direction * speed).and_then(|(candidate, residual)|
                            Ok((candidate, number(direction * residual)?))))?;
                    (sliding, direction * speed, false)
                }
            };
            let duration = number(candidate.thermal().time_s() - start)?;
            let mut velocities = [[0.0; 2]; 2];
            let mut displacement = [[0.0; 2]; 2];
            let mut momentum = [[0.0; 2]; 2];
            let mut kinetic = 0.0;
            let mut work = 0.0;
            for component in 0..2 {
                let initial_momentum = number(number(masses[0] * pair.velocities_mps[0][component])?
                    + number(masses[1] * pair.velocities_mps[1][component])?)?;
                let centre = number(number(initial_momentum + number(duration * number(
                    drive.body_forces_n[0][component] + drive.body_forces_n[1][component])?)?)? / total_mass)?;
                let relative = if component == axis { after } else { 0.0 };
                velocities[0][component] = number(centre + (masses[1] / total_mass) * relative)?;
                velocities[1][component] = number(centre - (masses[0] / total_mass) * relative)?;
                for body in 0..2 {
                    let old = pair.velocities_mps[body][component];
                    let next = velocities[body][component];
                    let mean = number(0.5 * old + 0.5 * next)?;
                    let dv = number(next - old)?;
                    displacement[body][component] = number(duration * mean)?;
                    kinetic = number(kinetic + number(masses[body] * dv * mean)?)?;
                    work = number(work + number(drive.body_forces_n[body][component] * displacement[body][component])?)?;
                    let sign = if body == 0 { -1.0 } else { 1.0 };
                    momentum[body][component] = number(number(masses[body] * dv)? - number(duration * number(
                        drive.body_forces_n[body][component] + sign * candidate.report().friction.traction_n[component])?)?)?;
                }
            }
            if momentum.iter().any(|r| det::hypot(r[0], r[1]) > controls.impulse_tolerance_n_s) {
                return Err(SlidingError::Invalid("collinear midpoint body impulse exceeds its independent gate"));
            }
            let (next_heat, contact, sliding) = candidate.into_parts();
            let residual = number(kinetic + sliding.thermal.stored_energy_change_j
                + sliding.friction.stored_change_j + sliding.partition.other_work_j()
                - work - sliding.external_heat_j)?;
            if residual.abs() > controls.total_energy_tolerance_j {
                return Err(SlidingError::Invalid("collinear kinetic, finite heat and actual force work do not close"));
            }
            report.kinetic_energy_change_j = number(report.kinetic_energy_change_j + kinetic)?;
            report.thermal_energy_change_j = number(report.thermal_energy_change_j + sliding.thermal.stored_energy_change_j)?;
            report.deposited_heat_j = number(report.deposited_heat_j + sliding.deposited_heat_j)?;
            report.external_work_j = number(report.external_work_j + work)?;
            report.total_energy_residual_j = number(report.total_energy_residual_j + residual)?;
            for body in 0..2 { for component in 0..2 {
                report.body_displacements_m[body][component] = number(
                    report.body_displacements_m[body][component] + displacement[body][component])?;
            }}
            report.segments.push(CollinearMotionSegment { start_time_s: start,
                end_time_s: next_heat.time_s(), relative_velocity_before_mps: before,
                relative_velocity_after_mps: after, reached_stop, sliding,
                body_displacements_m: displacement, body_momentum_residual_n_s: momentum,
                kinetic_energy_change_j: kinetic, external_work_j: work, total_energy_residual_j: residual });
            pair.inventory = pair.inventory.each_ref().map(|old| next_heat.cells()[&old.id].clone());
            pair.contact = contact;
            pair.velocities_mps = velocities;
            guess.fill(0.0);
            heat = next_heat;
        }
        if report.total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(SlidingError::Invalid("collinear cumulative physical energy residual exceeds its gate"));
        }
        report.motion_trials = trials.get();
        report.friction_trials = friction_trials.get();
        check()?;
        Ok(CollinearMotionCandidate { thermal: heat, pair, report })
    }
}

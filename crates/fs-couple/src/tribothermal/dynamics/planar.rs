//! Planar dry sliding with physical work, turning and slip-arrest events.
//!
//! Each segment uses `du = h*(v_before + v_after)/2` and the ORIGINAL
//! material/contact/enthalpy candidate. Its traction heats the finite bodies
//! through that same displacement, so there is no numerical kinetic loss to
//! hide missing stopping heat. For the current isotropic material law the
//! midpoint slip is parallel to `v_before + h*F_relative/(2*m_reduced)`; only
//! its scalar magnitude needs a root solve, even while the trajectory turns.
//!
//! A proposed arrest is solved at zero final relative velocity before static
//! friction is considered. Both components of its impulse must pass the
//! independent gate. Oblique motion is subdivided until its turning path or
//! limiting arrest is resolved; no velocity threshold selects sticking.
//! Constant collinear Coulomb motion recovers its physical stopping time and
//! distance, at every orientation of the tangent frame. General turning uses
//! convergent midpoint segments, with a caller maximum duration and explicit
//! subdivision/work budgets. The summed vector impulse and energy defects
//! are checked over the complete call before any history can be published.
//!
//! Material data are still queried at segment mean speed, accumulated path
//! length and endpoint temperatures. Thermal transport remains backward Euler;
//! temperature/history-dependent laws require time refinement. Scalar brackets
//! do not prove uniqueness for arbitrary weakening data. A previously stuck
//! segment does not localize thermal weakening within that segment. Fixed area,
//! fixed normal load and the original isotropic dry law remain assumptions;
//! no evolving normal geometry, anisotropic law or experimental claim is added.

use super::*;

/// Temporal resolution and bounded work, separate from physical acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarMotionSolve {
    /// Original material, constitutive and finite-enthalpy controls.
    pub sliding: SlidingSolve,
    /// Largest segment [s], subject to representable endpoint rounding.
    /// Refine this independently of root accuracy.
    pub maximum_step_s: f64,
    /// Accepted sliding/sticking segments in this atomic call.
    pub maximum_segments: usize,
    /// Complete original contact/thermal solves, including discarded probes.
    pub maximum_motion_trials: usize,
    /// Total rejected segment halvings within this call.
    pub maximum_halvings: usize,
    /// Safeguarded secant/bisection iterations per root, in 1..=128.
    pub maximum_bisections: usize,
    /// Euclidean impulse defect of each body, both per segment and cumulatively.
    pub impulse_tolerance_n_s: f64,
    /// Stopping-time bracket width, unless an exact scalar root is found.
    pub event_time_tolerance_s: f64,
    /// Absolute physical energy defect per segment and over the complete call [J].
    pub total_energy_tolerance_j: f64,
}

/// One accepted path segment, retaining the original source and thermal result.
pub struct PlanarMotionSegment<S> {
    /// Original synchronized mechanical/thermal interval start [s].
    pub start_time_s: f64,
    /// Accepted interval end, before any pending source/topology event [s].
    pub end_time_s: f64,
    /// A-minus-B tangent velocity at the start [m/s].
    pub relative_velocity_before_mps: [f64; 2],
    /// A-minus-B tangent velocity at the accepted endpoint [m/s].
    pub relative_velocity_after_mps: [f64; 2],
    /// A solved zero-velocity endpoint whose full vector impulse passed its gate.
    pub reached_stop: bool,
    /// Original source, slip, heat partition and finite-thermal accounts.
    pub sliding: SlidingReport<S>,
    /// Each body's midpoint displacement in the existing tangent frame [m].
    pub body_displacements_m: [[f64; 2]; 2],
    /// Original body momentum change minus applied/contact impulses [N s].
    pub body_momentum_residual_n_s: [[f64; 2]; 2],
    /// Both original inertial masses' kinetic energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Applied forces dotted with actual body displacement [J].
    pub external_work_j: f64,
    /// Kinetic + caloric + interface + outside loss - force work - outside heat [J].
    pub total_energy_residual_j: f64,
}

/// Complete interval accounting, without numerical kinetic loss or extra stores.
pub struct PlanarMotionReport<S> {
    /// Actual accepted path, in chronological order.
    pub segments: Vec<PlanarMotionSegment<S>>,
    /// Complete material/contact/thermal solves, including discarded candidates.
    pub motion_trials: usize,
    /// Original law evaluations, including all discarded thermal candidates.
    pub friction_trials: usize,
    /// Rejected turning/event segment subdivisions.
    pub halvings: usize,
    /// Sum of the accepted physical midpoint displacements [m].
    pub body_displacements_m: [[f64; 2]; 2],
    /// Sum of each body's independently gated vector impulse defects [N s].
    pub body_momentum_residual_n_s: [[f64; 2]; 2],
    /// Complete change in the original bodies' kinetic energy [J].
    pub kinetic_energy_change_j: f64,
    /// Complete change in the original finite enthalpy network [J].
    pub thermal_energy_change_j: f64,
    /// Physical friction heat entering the two original bodies [J].
    pub deposited_heat_j: f64,
    /// Applied force work over the accepted piecewise path [J].
    pub external_work_j: f64,
    /// Sum of physical energy defects, independently gated before publication [J].
    pub total_energy_residual_j: f64,
}

/// All segments' mechanical, material-history and finite-thermal state together.
pub struct PlanarMotionCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    pair: SlidingPair<'a>,
    report: PlanarMotionReport<S>,
}
impl<'a, 't, S> PlanarMotionCandidate<'a, 't, S> {
    /// Immutable accepted finite-thermal proposal.
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    /// Immutable accepted mechanical and material-history proposal.
    pub const fn pair(&self) -> &SlidingPair<'a> {
        &self.pair
    }
    /// Physical accounting and actual bounded solver work.
    pub const fn report(&self) -> &PlanarMotionReport<S> {
        &self.report
    }
    /// Consume the joint proposal; no independently committed half-step exists.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SlidingPair<'a>,
        PlanarMotionReport<S>,
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

/// The residual always has impulse units, for both speed and event-time roots.
fn root<T, E>(
    mut upper: f64,
    mut lower_residual: f64,
    upper_candidate: T,
    mut upper_residual: f64,
    controls: PlanarMotionSolve,
    time_root: bool,
    mut evaluate: impl FnMut(f64) -> Result<(T, f64), SlidingError<E>>,
) -> Result<(T, f64), SlidingError<E>> {
    if !upper.is_finite()
        || upper <= 0.0
        || !lower_residual.is_finite()
        || !upper_residual.is_finite()
        || lower_residual >= 0.0
        || upper_residual < 0.0
    {
        return Err(SlidingError::Invalid(
            "planar impulse root is not sign bracketed",
        ));
    }
    if upper_residual == 0.0 {
        return Ok((upper_candidate, upper));
    }
    drop(upper_candidate);
    let mut lower = 0.0;
    // A midpoint residual contributes twice its magnitude to the body impulse.
    // Leave room for accumulation and independently reconstructed roundoff.
    let tolerance = controls.impulse_tolerance_n_s / (8.0 * controls.maximum_segments as f64);
    for _ in 0..controls.maximum_bisections {
        let width = number(upper - lower)?;
        let secant = lower - lower_residual * width / (upper_residual - lower_residual);
        let middle =
            if secant.is_finite() && secant > lower + width / 16.0 && secant < upper - width / 16.0
            {
                secant
            } else {
                lower + 0.5 * width
            };
        if middle <= lower || middle >= upper {
            return Err(SlidingError::Invalid(
                "planar root interval is unrepresentable",
            ));
        }
        let (candidate, residual) = evaluate(middle)?;
        number::<E>(residual)?;
        if residual == 0.0
            || (residual.abs() <= tolerance
                && (!time_root || width <= controls.event_time_tolerance_s))
        {
            return Ok((candidate, middle));
        }
        if residual < 0.0 {
            lower = middle;
            lower_residual = residual;
        } else {
            upper = middle;
            upper_residual = residual;
        }
    }
    Err(SlidingError::MotionBudget)
}

impl<'a> SlidingPair<'a> {
    /// Advance general planar relative motion through its original physical owners.
    ///
    /// Forces and normal load are held in the existing fixed tangent frame.
    /// Source/topology events must already be cut by the caller. Each accepted
    /// segment retains physical midpoint work and path length; static reaction
    /// is probed only at an actual accepted zero relative velocity. Segments
    /// turning by more than a right angle are subdivided instead of accepting
    /// the implicit midpoint method's spurious mirror-velocity sticking root.
    ///
    /// An oblique limiting stop is an approximate event with its complete vector
    /// impulse explicitly bounded, not an exact analytic stopping trajectory.
    /// Decrease `maximum_step_s` for path accuracy and the impulse/time gates for
    /// event accuracy. An unresolved event, material/source refusal, cancellation
    /// or exhausted work budget leaves every caller-owned state unchanged.
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
    ) -> Result<PlanarMotionCandidate<'a, 't, S>, SlidingError<E>> {
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
            || thermal.time_s() != self.contact.time_s()
            || energy_guess_j.len() != thermal.cells().len()
            || energy_guess_j.iter().any(|v| !v.is_finite())
            || drive.body_forces_n.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(SlidingError::Invalid(
                "invalid planar motion clocks, forces or solve controls",
            ));
        }
        for old in &self.inventory {
            if thermal.cells().get(&old.id) != Some(old) {
                return Err(SlidingError::Invalid(
                    "planar motion must retain the original finite thermal inventory",
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
        let mut relative_force = [0.0; 2];
        for axis in 0..2 {
            relative_force[axis] = number(
                reduced_mass
                    * number(
                        drive.body_forces_n[0][axis] / masses[0]
                            - drive.body_forces_n[1][axis] / masses[1],
                    )?,
            )?;
        }
        let force_magnitude = norm(relative_force)?;
        let force_direction = if force_magnitude == 0.0 {
            [0.0; 2]
        } else {
            [
                relative_force[0] / force_magnitude,
                relative_force[1] / force_magnitude,
            ]
        };
        let mut pair = self.clone();
        let mut heat = thermal.clone();
        let mut guess = energy_guess_j.to_vec();
        let trials = Cell::new(0usize);
        let friction_trials = Cell::new(0usize);
        let mut report = PlanarMotionReport {
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
            let solve = |duration: f64,
                         increment: [f64; 2],
                         compliance: f64,
                         source_increment: [f64; 2]| {
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
                        "planar segment has no representable duration",
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
                    source_increment,
                    controls.sliding,
                    &guess,
                    |state, _| resolve(state),
                    || (*cancel.borrow_mut())(),
                )?;
                friction_trials.set(
                    friction_trials
                        .get()
                        .checked_add(candidate.report().friction_trials)
                        .ok_or(SlidingError::MotionBudget)?,
                );
                Ok(candidate)
            };
            let actual_duration = |duration: f64| {
                if duration == remaining {
                    Ok(remaining)
                } else {
                    number(number(start + duration)? - start)
                }
            };
            let proposed = remaining.min(controls.maximum_step_s);
            let proposed_end = number(start + proposed)?;
            // Repeated floating-point steps can land one clock ULP short of
            // the target. Resolve that final endpoint in this real segment,
            // instead of inventing a tiny extra segment whose slip cannot
            // advance the original accumulated history. A step that cannot
            // advance the clock still refuses below; no source/history guard
            // is relaxed and no representable multi-ULP interval is skipped.
            let proposed = if proposed_end > start && proposed_end.next_up() == target_s {
                remaining
            } else {
                proposed
            };
            let mut duration = actual_duration(proposed)?;
            let (candidate, after, reached_stop) = loop {
                check()?;
                if duration <= 0.0 {
                    return Err(SlidingError::Invalid("planar timestep is unrepresentable"));
                }
                if initial_speed == 0.0 {
                    let compliance = number(0.5 * duration * duration / reduced_mass)?;
                    if compliance <= 0.0 {
                        return Err(SlidingError::Invalid(
                            "planar static compliance is unrepresentable",
                        ));
                    }
                    let released = [
                        number(compliance * relative_force[0])?,
                        number(compliance * relative_force[1])?,
                    ];
                    let candidate = solve(duration, released, compliance, [0.0; 2])?;
                    if candidate.report().friction.regime == FrictionRegime::Sticking {
                        if candidate.report().friction.slip_increment_m != [0.0; 2] {
                            return Err(SlidingError::Invalid(
                                "native static probe retained nonzero planar slip",
                            ));
                        }
                        break (candidate, [0.0; 2], false);
                    }
                } else {
                    let direction = [before[0] / initial_speed, before[1] / initial_speed];
                    let at_stop = |h: f64| {
                        let h = actual_duration(h)?;
                        let increment =
                            [number(0.5 * h * before[0])?, number(0.5 * h * before[1])?];
                        let candidate = solve(h, increment, 0.0, increment)?;
                        let traction = candidate.report().friction.traction_n;
                        let residual = number(
                            -reduced_mass * initial_speed
                                - h * number(
                                    dot(relative_force, direction)? - dot(traction, direction)?,
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
                        let h = number(stop.thermal().time_s() - start)?;
                        let traction = stop.report().friction.traction_n;
                        let residual = [
                            number(
                                -reduced_mass * before[0] - h * (relative_force[0] - traction[0]),
                            )?,
                            number(
                                -reduced_mass * before[1] - h * (relative_force[1] - traction[1]),
                            )?,
                        ];
                        let cross =
                            direction[0] * force_direction[1] - direction[1] * force_direction[0];
                        // Above kinetic capacity an oblique Coulomb path turns
                        // through a nonzero speed minimum instead of stopping.
                        // A small transverse impulse alone cannot select the
                        // discontinuously different static branch. Only a
                        // collinear event (to direction roundoff), or a drive
                        // below kinetic capacity, admits this arrest limit.
                        let arrest_branch =
                            force_magnitude < norm(traction)? || cross.abs() <= 64.0 * f64::EPSILON;
                        if arrest_branch && norm(residual)? <= 0.5 * controls.impulse_tolerance_n_s
                        {
                            break (stop, [0.0; 2], true);
                        }
                    }
                }
                let free = [
                    number(before[0] + 0.5 * duration * relative_force[0] / reduced_mass)?,
                    number(before[1] + 0.5 * duration * relative_force[1] / reduced_mass)?,
                ];
                let free_speed = norm(free)?;
                let mut accepted = None;
                if free_speed > 0.0 {
                    let direction = [free[0] / free_speed, free[1] / free_speed];
                    let zero = solve(duration, [0.0; 2], 0.0, [0.0; 2])?;
                    let kinetic_limit = number(
                        drive.normal_load_n
                            * zero
                                .report()
                                .material
                                .law
                                .zero_speed_kinetic_coefficient()?,
                    )?;
                    let lower_residual =
                        number(0.5 * duration * kinetic_limit - reduced_mass * free_speed)?;
                    drop(zero);
                    // A nonpositive sliding root is the mirror-velocity artifact
                    // of midpoint capture, not a physical arrest of a moving body.
                    if lower_residual < 0.0 {
                        let at_speed = |speed: f64| {
                            let increment = [
                                number(duration * speed * direction[0])?,
                                number(duration * speed * direction[1])?,
                            ];
                            let candidate = solve(duration, increment, 0.0, increment)?;
                            let residual = number(
                                reduced_mass * number(speed - free_speed)?
                                    + 0.5
                                        * duration
                                        * dot(candidate.report().friction.traction_n, direction)?,
                            )?;
                            Ok((candidate, residual))
                        };
                        let (upper, upper_residual) = at_speed(free_speed)?;
                        let (candidate, speed) = root(
                            free_speed,
                            lower_residual,
                            upper,
                            upper_residual,
                            controls,
                            false,
                            at_speed,
                        )?;
                        let after = [
                            number(2.0 * speed * direction[0] - before[0])?,
                            number(2.0 * speed * direction[1] - before[1])?,
                        ];
                        let after_speed = norm(after)?;
                        let turn_before = if initial_speed == 0.0 {
                            0.0
                        } else {
                            (before[0] / initial_speed) * force_direction[1]
                                - (before[1] / initial_speed) * force_direction[0]
                        };
                        let turn_after = if after_speed == 0.0 {
                            0.0
                        } else {
                            (after[0] / after_speed) * force_direction[1]
                                - (after[1] / after_speed) * force_direction[0]
                        };
                        // With a held force and isotropic drag the velocity
                        // direction approaches the force direction without
                        // crossing it. Midpoint can overshoot this angular
                        // relaxation when its step is too long near arrest.
                        let monotone_turn = turn_before * turn_after >= 0.0
                            || turn_before.abs() <= 64.0 * f64::EPSILON
                            || turn_after.abs() <= 64.0 * f64::EPSILON;
                        if initial_speed == 0.0 || dot(before, after)? >= 0.0 && monotone_turn {
                            accepted = Some((candidate, after, false));
                        }
                    }
                }
                if let Some(accepted) = accepted {
                    break accepted;
                }
                if report.halvings >= controls.maximum_halvings {
                    return Err(SlidingError::MotionBudget);
                }
                report.halvings += 1;
                let next = actual_duration(0.5 * duration)?;
                if next <= 0.0 || next >= duration {
                    return Err(SlidingError::Invalid(
                        "planar turning/event subdivision is unrepresentable",
                    ));
                }
                duration = next;
            };
            let duration = number(candidate.thermal().time_s() - start)?;
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
                        "native planar body motion is unrepresentable",
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
                                            + sign * candidate.report().friction.traction_n[axis],
                                    )?,
                            )?,
                    )?;
                }
            }
            for residual in momentum {
                if norm(residual)? > controls.impulse_tolerance_n_s {
                    return Err(SlidingError::Invalid(
                        "planar midpoint body impulse exceeds its independent gate",
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
                    "planar kinetic, finite heat and actual force work do not close",
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
            report.segments.push(PlanarMotionSegment {
                start_time_s: start,
                end_time_s: next_heat.time_s(),
                relative_velocity_before_mps: before,
                relative_velocity_after_mps: after,
                reached_stop,
                sliding,
                body_displacements_m: displacement,
                body_momentum_residual_n_s: momentum,
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
        if report.total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(SlidingError::Invalid(
                "planar cumulative physical energy residual exceeds its gate",
            ));
        }
        for residual in report.body_momentum_residual_n_s {
            if norm(residual)? > controls.impulse_tolerance_n_s {
                return Err(SlidingError::Invalid(
                    "planar cumulative vector impulse exceeds its gate",
                ));
            }
        }
        report.motion_trials = trials.get();
        report.friction_trials = friction_trials.get();
        check()?;
        Ok(PlanarMotionCandidate {
            thermal: heat,
            pair,
            report,
        })
    }
}

//! Actual piecewise affine contact motion and the original force-law solves.

use super::*;

pub(super) struct Resolved {
    pub phases: Vec<SolidCoulombPhase>,
    pub stop_time: Option<f64>,
    pub stop_residual: Option<f64>,
}

pub(super) struct Engine<'a, 'contact, 'cx> {
    pub mobility: [[f64; 3]; 3],
    pub acceleration: [f64; 3],
    pub axis: [f64; 2],
    pub contact: &'a SolidSustainedContact<'contact>,
    pub kinetic: f64,
    pub duration: f64,
    pub controls: SolidSustainedControls,
    pub native_terms: u128,
    pub evaluations: usize,
    pub static_trials: usize,
    pub cx: &'a Cx<'cx>,
}

pub(super) fn collinear<E, N>(
    id: u64,
    axis: [f64; 2],
    velocity: [f64; 3],
    allowance: f64,
) -> Result<(), SolidSustainedError<E, N>> {
    let perpendicular = Dd::from_f64(axis[0]) * Dd::from_f64(velocity[2])
        - Dd::from_f64(axis[1]) * Dd::from_f64(velocity[1]);
    if finite(perpendicular)?.abs() > allowance {
        return Err(mode(id, "contact motion leaves the authored sliding axis"));
    }
    Ok(())
}

impl Engine<'_, '_, '_> {
    fn acceleration_for(&self, force: [f64; 3]) -> Result<[f64; 3], &'static str> {
        let mut result = [0.0; 3];
        for (a, value) in result.iter_mut().enumerate() {
            *value = finite(Dd::from_f64(self.acceleration[a]) + dot(&self.mobility[a], &force))?;
        }
        Ok(result)
    }

    fn along(&self, vector: [f64; 3]) -> Result<f64, &'static str> {
        finite(dot(&self.axis, &vector[1..]))
    }

    // Normal dwell plus the kinetic boundary determines a force without any
    // integration or regularization. The shared solver still verifies this
    // warm start against the actual duration and full coupled mobility.
    fn kinetic_force<E, N>(&self, sign: f64) -> Result<[f64; 3], SolidSustainedError<E, N>> {
        let denominator = Dd::from_f64(self.mobility[0][0])
            - Dd::from_f64(sign * self.kinetic) * dot(&self.mobility[0][1..], &self.axis);
        let denominator = finite(denominator)?;
        if denominator <= 0.0 {
            return Err(mode(
                self.contact.id,
                "positive kinetic normal mobility required",
            ));
        }
        let normal = -self.acceleration[0] / denominator;
        let force = [
            normal,
            -sign * self.kinetic * normal * self.axis[0],
            -sign * self.kinetic * normal * self.axis[1],
        ];
        if force.iter().any(|x| !x.is_finite()) || normal <= 0.0 {
            return Err(mode(
                self.contact.id,
                "strictly compressive finite contact force required",
            ));
        }
        Ok(force)
    }

    pub(super) fn resolve<B: Native, E>(
        &mut self,
        mut velocity: [f64; 3],
        mut initial_force: [f64; 3],
    ) -> Result<Resolved, SolidSustainedError<E, B::Error>> {
        let mut result = Resolved {
            phases: Vec::with_capacity(2),
            stop_time: None,
            stop_residual: None,
        };
        let mut offset = 0.0;
        let mut gap = self.contact.gap_m;
        // Every nonzero incoming slip first receives kinetic friction, even
        // below the numerical audit tolerance. Static capture is attempted only
        // from actual initial rest or after this owner's located kinetic stop.
        if velocity[1] != 0.0 || velocity[2] != 0.0 {
            let speed = self.along(velocity)?;
            if speed == 0.0 {
                return Err(mode(
                    self.contact.id,
                    "nonzero slip must resolve along the sliding axis",
                ));
            }
            let force = self.kinetic_force(speed.signum())?;
            let acceleration = self.acceleration_for(force)?;
            let along_acceleration = self.along(acceleration)?;
            let stopping_time =
                if speed.signum() != along_acceleration.signum() && along_acceleration != 0.0 {
                    let time = -speed / along_acceleration;
                    (time.is_finite() && time <= self.duration).then_some(time)
                } else {
                    None
                };
            let duration = stopping_time.unwrap_or(self.duration);
            let exact_duration = Dd::from_f64(duration);
            let solve = self.solve::<B, E>(self.kinetic, force, gap, velocity, exact_duration)?;
            let phase = self.phase::<E, B::Error>(
                CoulombPhaseMode::Sliding,
                0.0,
                exact_duration,
                self.kinetic,
                gap,
                velocity,
                solve,
            )?;
            velocity = phase.relative_velocity_after_mps;
            gap = phase.normal_gap_after_m;
            initial_force = phase.force_n;
            offset = duration;
            if let Some(time) = stopping_time {
                let residual = velocity[1].hypot(velocity[2]);
                if residual > self.controls.velocity_tolerance_mps {
                    return Err(mode(
                        self.contact.id,
                        "located stop exceeds velocity residual allowance",
                    ));
                }
                result.stop_time = Some(time);
                result.stop_residual = Some(residual);
            }
            result.phases.push(phase);
        }
        if offset < self.duration {
            let duration = Dd::from_f64(self.duration) - Dd::from_f64(offset);
            self.static_trials += 1;
            let trial = self.solve::<B, E>(
                self.contact.friction_coefficient,
                initial_force,
                gap,
                velocity,
                duration,
            )?;
            let acceleration = self.acceleration_for(trial.forces_n[0])?;
            let endpoint = advance_velocity(velocity, acceleration, duration)?;
            let slip = trial.paths_m[0][1].hypot(trial.paths_m[0][2]);
            let captures = slip <= self.controls.displacement_tolerance_m
                && endpoint[1].hypot(endpoint[2]) <= self.controls.velocity_tolerance_mps;
            let phase = if captures {
                self.phase::<E, B::Error>(
                    CoulombPhaseMode::Sticking,
                    offset,
                    duration,
                    self.contact.friction_coefficient,
                    gap,
                    velocity,
                    trial,
                )?
            } else {
                collinear::<E, B::Error>(
                    self.contact.id,
                    self.axis,
                    endpoint,
                    self.controls.velocity_tolerance_mps,
                )?;
                let direction = self.along([0.0, trial.paths_m[0][1], trial.paths_m[0][2]])?;
                if direction == 0.0 {
                    return Err(mode(self.contact.id, "breakaway direction is unresolved"));
                }
                let force = self.kinetic_force(direction.signum())?;
                let solve = self.solve::<B, E>(self.kinetic, force, gap, velocity, duration)?;
                self.phase::<E, B::Error>(
                    CoulombPhaseMode::Sliding,
                    offset,
                    duration,
                    self.kinetic,
                    gap,
                    velocity,
                    solve,
                )?
            };
            result.phases.push(phase);
        }
        Ok(result)
    }

    fn solve<B: Native, E>(
        &mut self,
        mu: f64,
        force: [f64; 3],
        gap: f64,
        velocity: [f64; 3],
        duration: Dd,
    ) -> Result<SolveReport, SolidSustainedError<E, B::Error>> {
        poll(self.cx)?;
        let rounded_duration = finite(duration)?;
        if rounded_duration <= 0.0 {
            return Err("positive representable Coulomb segment duration required".into());
        }
        let ratio = self.duration / rounded_duration;
        let scales = NaturalMapScales {
            normal_n_per_m: self.contact.scales.normal_n_per_m * ratio * ratio,
            tangent_n_per_m: self.contact.scales.tangent_n_per_m * ratio * ratio,
        };
        let mut controls = self.controls.solve;
        controls.maximum_response_evaluations = controls
            .maximum_response_evaluations
            .saturating_sub(self.evaluations);
        if controls.maximum_response_evaluations < 2 {
            return Err(SolidSustainedError::Budget {
                resource: "Coulomb event response evaluations",
                required: self.evaluations as u128 + 2,
                maximum: self.controls.solve.maximum_response_evaluations as u128,
            });
        }
        let previous = self.evaluations;
        let mut exhausted = None;
        let mut replay = None;
        let mut evaluations = self.evaluations;
        let mut cancelled = || self.cx.checkpoint().is_err();
        let report = fs_contact::sustained::solve::solve(
            &[mu],
            &[scales],
            &[force],
            controls,
            &mut |forces, stop| {
                if stop() {
                    return Err(B::cancelled());
                }
                let required = self.native_terms + 128 * (evaluations as u128 + 2);
                if required > self.controls.maximum_response_multiply_terms as u128 {
                    exhausted = Some(required);
                    return Err(B::invalid("Coulomb event response work budget exhausted"));
                }
                evaluations += 1;
                let acceleration = self.acceleration_for(forces[0]).map_err(B::invalid)?;
                let mut path =
                    displacement(velocity, acceleration, duration).map_err(B::invalid)?;
                path[0] = finite(Dd::from_f64(path[0]) + Dd::from_f64(gap)).map_err(B::invalid)?;
                replay = Some((forces[0], path));
                Ok(vec![path])
            },
            &mut cancelled,
        );
        self.evaluations = evaluations;
        if let Some(required) = exhausted {
            return Err(SolidSustainedError::Budget {
                resource: "Coulomb event multiply terms",
                required,
                maximum: self.controls.maximum_response_multiply_terms as u128,
            });
        }
        let report = report.map_err(|error| match error {
            SolveError::Cancelled => SolidSustainedError::Cancelled,
            SolveError::Response(error) if B::is_cancelled(&error) => {
                SolidSustainedError::Cancelled
            }
            error => SolidSustainedError::Solve(Box::new(error)),
        })?;
        if evaluations - previous != report.response_evaluations
            || replay != Some((report.forces_n[0], report.paths_m[0]))
        {
            return Err("missing final segment force replay".into());
        }
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    fn phase<E, N>(
        &self,
        mode: CoulombPhaseMode,
        start: f64,
        duration: Dd,
        mu: f64,
        gap: f64,
        velocity: [f64; 3],
        solve: SolveReport,
    ) -> Result<SolidCoulombPhase, SolidSustainedError<E, N>> {
        let force = solve.forces_n[0];
        let acceleration = self.acceleration_for(force)?;
        let after = advance_velocity(velocity, acceleration, duration)?;
        let travel = displacement(velocity, acceleration, duration)?;
        let phase = SolidCoulombPhase {
            mode,
            start_offset_s: start,
            duration_s: finite(duration)?,
            friction_coefficient: mu,
            force_n: force,
            relative_velocity_before_mps: velocity,
            relative_velocity_after_mps: after,
            relative_displacement_m: travel,
            normal_gap_before_m: gap,
            normal_gap_after_m: solve.paths_m[0][0],
            frictional_work_j: finite(dot(&force[1..], &travel[1..]))?,
            solve,
            duration,
        };
        audit::phase(&phase, self.axis, self.contact.id, self.controls)?;
        Ok(phase)
    }
}

fn advance_velocity(v: [f64; 3], a: [f64; 3], h: Dd) -> Result<[f64; 3], &'static str> {
    let mut result = [0.0; 3];
    for i in 0..3 {
        result[i] = finite(Dd::from_f64(v[i]) + Dd::from_f64(a[i]) * h)?;
    }
    Ok(result)
}

fn displacement(v: [f64; 3], a: [f64; 3], h: Dd) -> Result<[f64; 3], &'static str> {
    let mut result = [0.0; 3];
    for i in 0..3 {
        result[i] = finite((Dd::from_f64(v[i]) + Dd::from_f64(0.5) * Dd::from_f64(a[i]) * h) * h)?;
    }
    Ok(result)
}

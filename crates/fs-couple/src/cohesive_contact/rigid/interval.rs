//! Bounded continuation of the original finite-rotation rigid cohesive owner.
//!
//! Only its explicit unresolved Coulomb reversal triggers subdivision. Every
//! attempt uses the same native body, objective material and friction equations.
//! No velocity projection, extra impulse or numerical heating resolves a stop.
//! This brackets an admissible boundary of the discrete path model; it does not
//! certify an earliest continuous event or a unique softening branch.

use super::*;

/// Subdivision controls and cumulative limits for all accepted AND rejected
/// native solves. Physical acceptance tolerances remain those of the owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCohesiveIntervalConfig {
    /// Smallest refined trial; a shorter final requested remainder is allowed.
    pub minimum_step_s: f64,
    /// Admitted/refused bracket width, not a certified stopping-time error.
    pub reversal_time_tolerance_s: f64,
    pub maximum_trials: usize,
    pub maximum_substeps: usize,
    /// Attempts reserve remaining native allowances. Successful trials refund
    /// unused work; failed trials conservatively charge the full reservation.
    pub maximum_total_evaluations: usize,
    pub maximum_total_material_evaluations: usize,
    pub maximum_total_iterations: usize,
}

/// Sums along accepted original material paths, in the original site order.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RigidCohesiveIntervalSite {
    pub fracture_dissipation_j: f64,
    pub viscous_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub cohesive_spatial_work_j: f64,
}

/// Complete requested interval. The native owner's observations retain the
/// final accepted substep; this report contains the entire interval's accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidCohesiveInterval {
    pub initial_time_s: f64,
    pub time_s: f64,
    pub accepted_substeps: usize,
    pub attempted_substeps: usize,
    pub charged_evaluations: usize,
    pub charged_material_evaluations: usize,
    pub charged_iterations: usize,
    pub energy_before_j: f64,
    pub energy_after_j: f64,
    pub external_work_j: f64,
    pub fracture_dissipation_j: f64,
    pub viscous_dissipation_j: f64,
    pub friction_dissipation_j: f64,
    pub sites: Vec<RigidCohesiveIntervalSite>,
    /// Independently reconstructed from endpoint storage and accepted work.
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub external_impulse_n_s: Vec3,
    pub external_angular_impulse_n_m_s: Vec3,
    pub linear_momentum_residual_n_s: Vec3,
    pub linear_momentum_tolerance_n_s: f64,
    pub angular_momentum_residual_n_m_s: Vec3,
    pub angular_momentum_tolerance_n_m_s: f64,
}

impl<'a> RigidCohesiveSystem<'a> {
    /// Advance a complete interval under CONSTANT original body loads. Each
    /// native substep is at most the owner's nominal step_duration_s(). Sources
    /// remain at their currently accepted temperature, and every original
    /// damage, viscosity and static/kinetic history continues from that state.
    ///
    /// Any refusal or cancellation restores the entire incoming accepted state,
    /// including observations, friction history and the nominal clock segment.
    /// Only unresolved friction reversals refine; all other refusals propagate.
    pub fn advance_interval(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        duration_s: f64,
        controls: RigidCohesiveIntervalConfig,
    ) -> Result<RigidCohesiveInterval, RigidCohesiveError> {
        self.interval_inner(loads, duration_s, controls, None)
    }

    pub fn advance_interval_under_gate(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        duration_s: f64,
        controls: RigidCohesiveIntervalConfig,
        gate: &CancelGate,
    ) -> Result<RigidCohesiveInterval, RigidCohesiveError> {
        self.interval_inner(loads, duration_s, controls, Some(gate))
    }

    fn interval_inner(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        duration_s: f64,
        controls: RigidCohesiveIntervalConfig,
        gate: Option<&CancelGate>,
    ) -> Result<RigidCohesiveInterval, RigidCohesiveError> {
        poll(gate)?;
        positive(duration_s)?;
        positive(controls.minimum_step_s)?;
        positive(controls.reversal_time_tolerance_s)?;
        if controls.maximum_trials == 0
            || controls.maximum_substeps == 0
            || controls.maximum_total_evaluations == 0
            || controls.maximum_total_material_evaluations == 0
            || controls.maximum_total_iterations == 0
        {
            return Err(RigidCohesiveError::Budget(
                "nonzero rigid interval work limits required",
            ));
        }
        let end_time = finite(self.time_s + duration_s)?;
        if end_time <= self.time_s || self.thermal_active {
            return Err(RigidCohesiveError::Invalid(
                "rigid interval clock or active temperature stage",
            ));
        }
        let original = Checkpoint::new(self)?;
        let result = self.run_interval(loads, end_time, controls, gate);
        if result.is_err() {
            original.restore(self);
        }
        result
    }

    fn run_interval(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        end_time: f64,
        controls: RigidCohesiveIntervalConfig,
        gate: Option<&CancelGate>,
    ) -> Result<RigidCohesiveInterval, RigidCohesiveError> {
        // Allocate three small accepted-state buffers once, never native solver
        // matrices or a second mechanics owner and never inside Newton trials.
        let mut predecessor = Checkpoint::new(self)?;
        let mut best = Checkpoint::new(self)?;
        let initial_momentum = momentum(self)?;
        let mut result = RigidCohesiveInterval::new(self)?;
        let mut work = Work {
            controls,
            trials: 0,
            evaluations: 0,
            material: 0,
            iterations: 0,
        };
        while self.time_s < end_time {
            poll(gate)?;
            if result.accepted_substeps >= controls.maximum_substeps {
                return Err(RigidCohesiveError::Budget("rigid interval substeps"));
            }
            predecessor.save(self);
            let remaining = end_time - self.time_s;
            let h = remaining.min(self.dt_s);
            // Publish the requested endpoint exactly. Adding a rounded
            // (end-current) back to current can otherwise overshoot by one ULP.
            let at = if h == remaining {
                end_time
            } else {
                finite(self.time_s + h)?
            };
            let (frame, accepted_h) = match work.attempt(self, loads, h, at, gate) {
                Ok(frame) => (frame, h),
                Err(RigidCohesiveError::UnresolvedFrictionReversal { mut site }) => {
                    let mut low = 0.0;
                    let mut high = h;
                    loop {
                        if low > 0.0 && high - low <= controls.reversal_time_tolerance_s {
                            best.restore(self);
                            break (*self.last_frame.as_ref().expect("admitted bracket"), low);
                        }
                        let mid = low + (high - low) * 0.5;
                        if mid < controls.minimum_step_s || mid <= low || mid >= high {
                            return Err(RigidCohesiveError::Budget(
                                "rigid reversal time resolution",
                            ));
                        }
                        predecessor.restore(self);
                        let at = finite(self.time_s + mid)?;
                        match work.attempt(self, loads, mid, at, gate) {
                            Ok(frame) => {
                                if self
                                    .friction
                                    .as_ref()
                                    .expect("friction reversal owner")
                                    .resolved_stop(self, site)?
                                {
                                    break (frame, mid);
                                }
                                low = mid;
                                best.save(self);
                            }
                            Err(RigidCohesiveError::UnresolvedFrictionReversal { site: next }) => {
                                high = mid;
                                site = next;
                            }
                            Err(error) => return Err(error),
                        }
                    }
                }
                Err(error) => return Err(error),
            };
            result.append(self, &predecessor, loads, frame, accepted_h)?;
        }
        result.attempted_substeps = work.trials;
        result.charged_evaluations = work.evaluations;
        result.charged_material_evaluations = work.material;
        result.charged_iterations = work.iterations;
        result.energy_after_j = self.total_energy_j()?;
        result.energy_residual_j = finite(
            result.energy_after_j - result.energy_before_j
                + result.fracture_dissipation_j
                + result.viscous_dissipation_j
                + result.friction_dissipation_j
                - result.external_work_j,
        )?;
        balance(
            "rigid interval energy",
            result.energy_residual_j,
            result.energy_tolerance_j,
        )?;
        let final_momentum = momentum(self)?;
        result.linear_momentum_residual_n_s = final_momentum
            .0
            .sub(initial_momentum.0)
            .sub(result.external_impulse_n_s);
        result.angular_momentum_residual_n_m_s = final_momentum
            .1
            .sub(initial_momentum.1)
            .sub(result.external_angular_impulse_n_m_s);
        balance(
            "rigid interval linear momentum",
            norm(result.linear_momentum_residual_n_s),
            result.linear_momentum_tolerance_n_s,
        )?;
        balance(
            "rigid interval angular momentum",
            norm(result.angular_momentum_residual_n_m_s),
            result.angular_momentum_tolerance_n_m_s,
        )?;
        poll(gate)?;
        self.clock_anchor = Some((self.steps, self.time_s));
        Ok(result)
    }
}

struct Work {
    controls: RigidCohesiveIntervalConfig,
    trials: usize,
    evaluations: usize,
    material: usize,
    iterations: usize,
}
impl Work {
    fn attempt(
        &mut self,
        system: &mut RigidCohesiveSystem<'_>,
        loads: [RigidCohesiveLoad; 2],
        h: f64,
        at: f64,
        gate: Option<&CancelGate>,
    ) -> Result<RigidCohesiveFrame, RigidCohesiveError> {
        poll(gate)?;
        if self.trials >= self.controls.maximum_trials {
            return Err(RigidCohesiveError::Budget("rigid interval trials"));
        }
        let config = system.config;
        let nominal_h = system.dt_s;
        let iteration_limit = system.friction.as_ref().map_or(
            config.maximum_iterations,
            friction::Friction::iteration_limit,
        );
        let evaluations = config
            .maximum_evaluations
            .min(self.controls.maximum_total_evaluations - self.evaluations);
        let material = config
            .maximum_material_evaluations
            .min(self.controls.maximum_total_material_evaluations - self.material);
        let iterations =
            iteration_limit.min(self.controls.maximum_total_iterations - self.iterations);
        if evaluations == 0 || material == 0 || iterations == 0 {
            return Err(RigidCohesiveError::Budget("rigid interval cumulative work"));
        }
        self.trials += 1;
        self.evaluations += evaluations;
        self.material += material;
        self.iterations += iterations;
        system.dt_s = h;
        system.config.maximum_evaluations = evaluations;
        system.config.maximum_material_evaluations = material;
        system.config.maximum_iterations = iterations;
        if let Some(friction) = &mut system.friction {
            friction.set_iteration_limit(iterations);
        }
        let trial = system.step_at_with_acceptance(loads, at, gate, |frame, _| {
            Ok::<_, RigidCohesiveError>(*frame)
        });
        // Restore all numerical controls before propagating ANY native refusal.
        system.dt_s = nominal_h;
        system.config = config;
        if let Some(friction) = &mut system.friction {
            friction.set_iteration_limit(iteration_limit);
        }
        let frame = trial?;
        self.evaluations -= evaluations - frame.evaluations;
        self.material -= material - frame.material_evaluations;
        self.iterations -= iterations - frame.iterations;
        Ok(frame)
    }
}

/// Fixed-temperature accepted state only. Original source pointers survive in
/// the actual points and friction laws; no bond or branch is reconstructed.
struct Checkpoint<'a> {
    bodies: [RigidBodyState; 2],
    points: Vec<SourceBoundObjectiveCohesivePoint<'a>>,
    frames: Vec<Option<RigidCohesiveSiteFrame>>,
    friction: Option<friction::Checkpoint<'a>>,
    steps: u64,
    time: f64,
    anchor: Option<(u64, f64)>,
    wrenches: [Wrench; 2],
    frame: Option<RigidCohesiveFrame>,
}
impl<'a> Checkpoint<'a> {
    fn new(system: &RigidCohesiveSystem<'a>) -> Result<Self, RigidCohesiveError> {
        let mut points = Vec::new();
        let mut frames = Vec::new();
        points
            .try_reserve_exact(system.sites.len())
            .map_err(|_| RigidCohesiveError::Budget("interval material checkpoint"))?;
        frames
            .try_reserve_exact(system.sites.len())
            .map_err(|_| RigidCohesiveError::Budget("interval observation checkpoint"))?;
        points.extend(system.sites.iter().map(|s| s.point));
        frames.extend_from_slice(&system.site_frames);
        Ok(Self {
            bodies: system.bodies.each_ref().map(|b| b.state),
            points,
            frames,
            friction: system
                .friction
                .as_ref()
                .map(friction::Friction::checkpoint)
                .transpose()?,
            steps: system.steps,
            time: system.time_s,
            anchor: system.clock_anchor,
            wrenches: system.last_wrenches,
            frame: system.last_frame,
        })
    }
    fn save(&mut self, system: &RigidCohesiveSystem<'a>) {
        self.bodies = system.bodies.each_ref().map(|b| b.state);
        for (point, site) in self.points.iter_mut().zip(&system.sites) {
            *point = site.point;
        }
        self.frames.copy_from_slice(&system.site_frames);
        if let (Some(owner), Some(saved)) = (&system.friction, &mut self.friction) {
            owner.save(saved);
        }
        self.steps = system.steps;
        self.time = system.time_s;
        self.anchor = system.clock_anchor;
        self.wrenches = system.last_wrenches;
        self.frame = system.last_frame;
    }
    fn restore(&self, system: &mut RigidCohesiveSystem<'a>) {
        for (body, state) in system.bodies.iter_mut().zip(self.bodies) {
            body.state = state;
        }
        for (site, point) in system.sites.iter_mut().zip(&self.points) {
            site.point = *point;
        }
        system.site_frames.copy_from_slice(&self.frames);
        if let (Some(owner), Some(saved)) = (&mut system.friction, &self.friction) {
            owner.restore(saved);
        }
        system.steps = self.steps;
        system.time_s = self.time;
        system.clock_anchor = self.anchor;
        system.last_wrenches = self.wrenches;
        system.last_frame = self.frame;
    }
}

impl RigidCohesiveInterval {
    fn new(system: &RigidCohesiveSystem<'_>) -> Result<Self, RigidCohesiveError> {
        let mut sites = Vec::new();
        sites
            .try_reserve_exact(system.sites.len())
            .map_err(|_| RigidCohesiveError::Budget("interval site accounts"))?;
        sites.resize(system.sites.len(), RigidCohesiveIntervalSite::default());
        Ok(Self {
            initial_time_s: system.time_s,
            time_s: system.time_s,
            accepted_substeps: 0,
            attempted_substeps: 0,
            charged_evaluations: 0,
            charged_material_evaluations: 0,
            charged_iterations: 0,
            energy_before_j: system.total_energy_j()?,
            energy_after_j: 0.0,
            external_work_j: 0.0,
            fracture_dissipation_j: 0.0,
            viscous_dissipation_j: 0.0,
            friction_dissipation_j: 0.0,
            sites,
            energy_residual_j: 0.0,
            energy_tolerance_j: 0.0,
            external_impulse_n_s: Vec3::ZERO,
            external_angular_impulse_n_m_s: Vec3::ZERO,
            linear_momentum_residual_n_s: Vec3::ZERO,
            linear_momentum_tolerance_n_s: 0.0,
            angular_momentum_residual_n_m_s: Vec3::ZERO,
            angular_momentum_tolerance_n_m_s: 0.0,
        })
    }
    fn append(
        &mut self,
        system: &RigidCohesiveSystem<'_>,
        before: &Checkpoint<'_>,
        loads: [RigidCohesiveLoad; 2],
        frame: RigidCohesiveFrame,
        h: f64,
    ) -> Result<(), RigidCohesiveError> {
        self.accepted_substeps += 1;
        self.time_s = frame.time_s;
        self.external_work_j = finite(self.external_work_j + frame.external_work_j)?;
        self.fracture_dissipation_j =
            finite(self.fracture_dissipation_j + frame.fracture_dissipation_j)?;
        self.viscous_dissipation_j =
            finite(self.viscous_dissipation_j + frame.viscous_dissipation_j)?;
        self.friction_dissipation_j =
            finite(self.friction_dissipation_j + frame.friction_dissipation_j)?;
        self.energy_tolerance_j = finite(self.energy_tolerance_j + frame.energy_tolerance_j)?;
        self.linear_momentum_tolerance_n_s =
            finite(self.linear_momentum_tolerance_n_s + frame.linear_momentum_tolerance_n_s)?;
        self.angular_momentum_tolerance_n_m_s =
            finite(self.angular_momentum_tolerance_n_m_s + frame.angular_momentum_tolerance_n_m_s)?;
        for (total, site) in self.sites.iter_mut().zip(&system.site_frames) {
            let site = site
                .as_ref()
                .expect("complete accepted native observations");
            total.fracture_dissipation_j =
                finite(total.fracture_dissipation_j + site.fracture_dissipation_j)?;
            total.viscous_dissipation_j =
                finite(total.viscous_dissipation_j + site.viscous_dissipation_j)?;
            total.friction_dissipation_j = finite(
                total.friction_dissipation_j
                    + site.friction.map_or(0.0, |f| f.friction_dissipation_j),
            )?;
            total.cohesive_spatial_work_j =
                finite(total.cohesive_spatial_work_j + site.spatial_work_j)?;
        }
        for body in 0..2 {
            let force = loads[body].force_world_n.add(
                system
                    .gravity()
                    .acceleration_world()
                    .scale(system.bodies[body].properties.mass()),
            );
            let mid = before.bodies[body]
                .pose()
                .position_world()
                .scale(0.5)
                .add(system.bodies[body].state.pose().position_world().scale(0.5));
            self.external_impulse_n_s = self.external_impulse_n_s.add(force.scale(h));
            self.external_angular_impulse_n_m_s = self
                .external_angular_impulse_n_m_s
                .add(mid.cross(force).add(loads[body].moment_world_n_m).scale(h));
        }
        finite_vector(self.external_impulse_n_s)?;
        finite_vector(self.external_angular_impulse_n_m_s)
    }
}

fn momentum(system: &RigidCohesiveSystem<'_>) -> Result<(Vec3, Vec3), RigidCohesiveError> {
    let mut p = Vec3::ZERO;
    let mut l = Vec3::ZERO;
    for body in &system.bodies {
        let d = system.integrator.diagnostics(body.state, body.properties)?;
        p = p.add(d.linear_momentum_world);
        l = l.add(d.angular_momentum_world);
    }
    finite_vector(p)?;
    finite_vector(l)?;
    Ok((p, l))
}

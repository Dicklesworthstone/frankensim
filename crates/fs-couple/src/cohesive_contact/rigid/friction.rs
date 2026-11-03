//! Original-body, objective damaged-area friction with joint tangent disks.
//!
//! A site's mean capacity is mu times the ORIGINAL material integral of
//! damage times compression, over the same local path as its cohesive law.
//! Two tangent efforts per site join the twelve native body-wrench unknowns.
//! The cone scale only conditions a natural-map residual; it is not a spring.
//! Static capture, kinetic sliding and reclosure use the shared damaged-area
//! law. An already sliding site uses its kinetic capture disk for this step.
//!
//! Damped least-squares corrections permit redundant sticking sites without
//! adding physical compliance. They select a numerical representative of any
//! indeterminate local reactions. All accepted forces satisfy the ORIGINAL
//! unregularized disks and objective work equations. This is a discrete path
//! model, not an exact stopping-time or unique softening-branch claim. A
//! reversal of resolved endpoint slip refuses; subdivide that mechanical step
//! to resolve the stop. Velocity resolution is independent of timestep, so
//! resolved sliding retains its physical loss even for tiny net increments.

use super::solve::{Evaluation, Work, pack, terminal, unpack};
use super::*;
use crate::cohesive_contact::CohesiveContactConfig;
use crate::cohesive_contact::frictional_damage::local::LocalDamageFriction;
use crate::cohesive_contact::frictional_damage::{
    DamageFrictionConfig, DamageFrictionRegime, DamageFrictionState,
};
use crate::cohesive_contact::reclosure::{ReclosureActivation, ResidualCoulombSource};
use fs_blake3::DomainHasher;

const MAX_DIM: usize = 12 + 2 * 64;

/// Explicit numerical controls additional to the original body/source limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDamageFrictionConfig {
    /// All joint Newton corrections and static-to-kinetic switches, at most 128.
    pub maximum_iterations: usize,
    /// Absolute actual local sticking-slip tolerance [m].
    pub slip_absolute_tolerance_m: f64,
    pub slip_relative_tolerance: f64,
    /// Absolute objective tangential velocity resolution [m/s], independent of
    /// timestep. A small displacement alone cannot establish static capture.
    pub velocity_absolute_tolerance_m_s: f64,
    /// Dimensionless damping of scaled Newton NORMAL equations, in (0,1).
    /// Does not enter physical force, work, storage or contact compliance.
    pub newton_damping: f64,
}

/// One original residual-surface source per existing cohesive quadrature site.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDamageFrictionSite<'a> {
    pub source: &'a ResidualCoulombSource,
    /// Explicit adoption of the existing interface's tangential occupancy.
    pub activation: ReclosureActivation,
}

/// Accepted finite-step objective efforts and physical friction loss.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDamageFrictionSiteFrame {
    pub regime: DamageFrictionRegime,
    pub local_force_n: [f64; 2],
    pub force_on_a_n: Vec3,
    pub frame_couple_on_a_n_m: Vec3,
    /// Original path-mean damage times compression, multiplied by site area.
    pub damaged_compression_force_n: f64,
    pub slip_increment_m: [f64; 2],
    pub friction_dissipation_j: f64,
    pub spatial_work_j: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Trial {
    pub(super) frame: RigidDamageFrictionSiteFrame,
}
impl Trial {
    pub(super) const EMPTY: Self = Self {
        frame: RigidDamageFrictionSiteFrame {
            regime: DamageFrictionRegime::Inactive,
            local_force_n: [0.0; 2],
            force_on_a_n: Vec3::ZERO,
            frame_couple_on_a_n_m: Vec3::ZERO,
            damaged_compression_force_n: 0.0,
            slip_increment_m: [0.0; 2],
            friction_dissipation_j: 0.0,
            spatial_work_j: 0.0,
        },
    };

    pub(super) fn new(
        site: &RigidCohesiveSite<'_>,
        numerical: &super::super::objective::ObjectiveCohesiveStep<'_>,
        pressure_pa: f64,
        force_local: [f64; 2],
    ) -> Result<Self, RigidCohesiveError> {
        let old = site.point.numerical();
        let new = &numerical.point;
        let local = Vec3::new(0.0, force_local[0], force_local[1]);
        let force = old
            .frame()
            .rotate_body_to_world(local)
            .scale(0.5)
            .add(new.frame().rotate_body_to_world(local).scale(0.5));
        let middle = old
            .separation_world_m()
            .scale(0.5)
            .add(new.separation_world_m().scale(0.5));
        let couple = middle.cross(force);
        let jump0 = old.material().jump_m();
        let jump1 = new.material().jump_m();
        let slip = [jump1[1] - jump0[1], jump1[2] - jump0[2]];
        let spatial_work = force.dot(new.separation_world_m().sub(old.separation_world_m()))
            - couple.dot(numerical.cayley_increment_world);
        finite_vector(force)?;
        finite_vector(couple)?;
        finite(spatial_work)?;
        Ok(Self {
            frame: RigidDamageFrictionSiteFrame {
                local_force_n: force_local,
                force_on_a_n: force,
                frame_couple_on_a_n_m: couple,
                damaged_compression_force_n: finite(old.area_m2() * pressure_pa)?,
                slip_increment_m: slip,
                spatial_work_j: spatial_work,
                ..Self::EMPTY.frame
            },
        })
    }
}

impl<'a> RigidCohesiveSystem<'a> {
    /// Select damaged-area friction without recreating any body, bond, material
    /// state or clock. All original source assignments are checked first.
    pub fn enable_damaged_area_friction(
        &mut self,
        sites: Vec<RigidDamageFrictionSite<'a>>,
        config: RigidDamageFrictionConfig,
        gate: &CancelGate,
    ) -> Result<(), RigidCohesiveError> {
        poll(Some(gate))?;
        if self.friction.is_some() {
            return Err(RigidCohesiveError::Invalid(
                "damaged-area friction already selected",
            ));
        }
        let owner = Friction::new(self, sites, config, gate)?;
        let mut staged = Vec::new();
        staged
            .try_reserve_exact(self.sites.len())
            .map_err(|_| RigidCohesiveError::Budget("friction observation allocation"))?;
        staged.resize(self.sites.len(), Trial::EMPTY);
        poll(Some(gate))?;
        self.friction = Some(owner);
        self.staged_friction = staged;
        Ok(())
    }

    pub fn damaged_area_friction_sources(&self) -> Option<&[RigidDamageFrictionSite<'a>]> {
        self.friction.as_ref().map(|f| f.sources.as_slice())
    }
    pub fn damaged_area_friction_states(&self) -> Option<&[DamageFrictionState<'a>]> {
        self.friction.as_ref().map(|f| f.states.as_slice())
    }
}

pub(super) struct Friction<'a> {
    sources: Vec<RigidDamageFrictionSite<'a>>,
    laws: Vec<LocalDamageFriction<'a>>,
    states: Vec<DamageFrictionState<'a>>,
    staged: Vec<DamageFrictionState<'a>>,
    temperature: Vec<Option<TemperatureStage<'a>>>,
    last_forces: [[f64; 2]; 64],
    staged_forces: [[f64; 2]; 64],
    incoming_velocity: [[f64; 2]; 64],
    config: RigidDamageFrictionConfig,
    cone_scale: [f64; 64],
    coordinate_scale: [f64; MAX_DIM],
    jacobian: Vec<f64>,
    normal: Vec<f64>,
    linear: LuWorkspace,
}

#[derive(Clone, Copy)]
struct TemperatureStage<'a> {
    law: LocalDamageFriction<'a>,
    state: DamageFrictionState<'a>,
}

#[derive(Clone, Copy)]
struct Joint {
    native: Evaluation,
    residual: [f64; MAX_DIM],
    tolerance: [f64; MAX_DIM],
    slips: [f64; 64],
    pressure: [f64; 64],
    slip_tolerance: [f64; 64],
    endpoint_velocity: [[f64; 2]; 64],
}
impl Joint {
    fn worst(&self, n: usize) -> f64 {
        (0..n)
            .map(|i| self.residual[i].abs() / self.tolerance[i])
            .fold(0.0, f64::max)
    }
    fn merit(&self, n: usize) -> f64 {
        (0..n)
            .map(|i| (self.residual[i] / self.tolerance[i]).powi(2))
            .sum()
    }
}

impl<'a> Friction<'a> {
    fn new(
        system: &RigidCohesiveSystem<'a>,
        sources: Vec<RigidDamageFrictionSite<'a>>,
        config: RigidDamageFrictionConfig,
        gate: &CancelGate,
    ) -> Result<Self, RigidCohesiveError> {
        if sources.len() != system.sites.len()
            || !(1..=128).contains(&config.maximum_iterations)
            || !config.newton_damping.is_finite()
            || !(0.0..1.0).contains(&config.newton_damping)
            || config.newton_damping == 0.0
            || !config.slip_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&config.slip_relative_tolerance)
            || config.slip_relative_tolerance == 0.0
        {
            return Err(RigidCohesiveError::Invalid(
                "original friction sources and bounded numerical controls required",
            ));
        }
        positive(config.slip_absolute_tolerance_m)?;
        positive(config.velocity_absolute_tolerance_m_s)?;
        let c = system.config;
        let local_config = DamageFrictionConfig {
            contact: CohesiveContactConfig {
                max_iterations: c.maximum_iterations,
                max_material_evaluations: c.maximum_material_evaluations,
                max_backtracks: c.maximum_backtracks,
                jump_absolute_tolerance_m: config.slip_absolute_tolerance_m,
                jump_relative_tolerance: config.slip_relative_tolerance,
                maximum_jump_increment_m: c.maximum_jump_increment_m,
                maximum_force_n: c.maximum_site_force_n,
                material_energy_tolerance_j_m2: c.objective.material_energy_tolerance_j_m2,
            },
            max_friction_iterations: config.maximum_iterations,
            force_absolute_tolerance_n: c.force_absolute_tolerance_n,
            force_relative_tolerance: c.relative_tolerance,
            work_absolute_tolerance_j: c.objective.spatial_work_tolerance_j,
            work_relative_tolerance: c.energy_relative_tolerance,
        };
        let mut laws = Vec::new();
        let mut states = Vec::new();
        let count = sources.len();
        laws.try_reserve_exact(count)
            .map_err(|_| RigidCohesiveError::Budget("friction laws allocation"))?;
        states
            .try_reserve_exact(count)
            .map_err(|_| RigidCohesiveError::Budget("friction states allocation"))?;
        let mut cone_scale = [0.0; 64];
        let mut coordinate_scale = [0.0; MAX_DIM];
        for i in 0..12 {
            let absolute = if i % 6 < 3 {
                c.force_absolute_tolerance_n
            } else {
                c.torque_absolute_tolerance_n_m
            };
            let maximum = if i % 6 < 3 {
                c.maximum_total_force_n
            } else {
                c.maximum_total_torque_n_m
            };
            coordinate_scale[i] = (absolute / c.relative_tolerance).min(maximum).max(absolute);
        }
        let total_area: f64 = system
            .sites
            .iter()
            .map(|s| s.point.numerical().area_m2())
            .sum();
        let translation_response = 0.5
            * system.dt_s
            * system.dt_s
            * (1.0 / system.bodies[0].properties.mass() + 1.0 / system.bodies[1].properties.mass());
        for (index, (source, site)) in sources.iter().zip(&system.sites).enumerate() {
            poll(Some(gate))?;
            let cohesive = site.point.source();
            let area = site.point.numerical().area_m2();
            let mut hash = DomainHasher::new("org.frankensim.fs-couple.rigid-damage-friction.v1");
            hash.update(cohesive.resolved().identity().as_bytes());
            hash.update(source.source.resolved().identity().as_bytes());
            hash.update(&(index as u64).to_le_bytes());
            hash.update(&area.to_bits().to_le_bytes());
            let law = LocalDamageFriction::new(
                cohesive,
                source.source,
                area,
                local_config,
                hash.finalize(),
            )?;
            let state = law.adopt(
                site.point.numerical().material(),
                [&system.bodies[0].material, &system.bodies[1].material],
                source.activation,
            )?;
            cone_scale[index] = positive(
                1.0 / translation_response.max(
                    cohesive.law().onset_separation_m() / (area * cohesive.law().strength_pa()),
                ),
            )?;
            let scale = positive((area / total_area).sqrt() * coordinate_scale[0])?;
            coordinate_scale[12 + 2 * index] = scale;
            coordinate_scale[13 + 2 * index] = scale;
            laws.push(law);
            states.push(state);
        }
        let dim = 12 + 2 * count;
        let mut staged = Vec::new();
        staged
            .try_reserve_exact(count)
            .map_err(|_| RigidCohesiveError::Budget("friction staged states allocation"))?;
        staged.extend_from_slice(&states);
        let mut temperature = Vec::new();
        temperature
            .try_reserve_exact(count)
            .map_err(|_| RigidCohesiveError::Budget("friction temperature workspace allocation"))?;
        temperature.resize_with(count, || None);
        Ok(Self {
            sources,
            laws,
            states,
            staged,
            temperature,
            last_forces: [[0.0; 2]; 64],
            staged_forces: [[0.0; 2]; 64],
            incoming_velocity: [[0.0; 2]; 64],
            config,
            cone_scale,
            coordinate_scale,
            jacobian: zeros(dim * dim)?,
            normal: zeros(dim * dim)?,
            linear: LuWorkspace::new(dim)?,
        })
    }

    pub(super) fn stage_temperature(
        &mut self,
        system: &RigidCohesiveSystem<'a>,
        next_sources: &[&'a ResidualCoulombSource],
        temperature_axis: &str,
        gate: Option<&CancelGate>,
    ) -> Result<(), RigidCohesiveError> {
        if next_sources.len() != self.laws.len() {
            return Err(RigidCohesiveError::Invalid(
                "one temperature successor per original friction site required",
            ));
        }
        let bodies = system.bodies.each_ref().map(|body| body.state);
        for (p, source) in next_sources.iter().enumerate() {
            poll(gate)?;
            if self.states[p].point() != system.sites[p].point.numerical().material() {
                return Err(RigidCohesiveError::Invalid(
                    "original cohesive/friction histories differ before temperature transfer",
                ));
            }
            let point = system.staged_temperature[p]
                .as_ref()
                .ok_or(RigidCohesiveError::Invalid(
                    "incomplete cohesive temperature stage",
                ))?
                .point();
            let (law, state) = self.laws[p].adopt_temperature_successor(
                self.states[p],
                point.source(),
                *source,
                point.numerical().material(),
                temperature_axis,
            )?;
            // Newly damaged compressed material can already be moving. Use the
            // SAME resolved-velocity threshold as the original stop audit;
            // numerical residue from an accepted static solve cannot by itself
            // replace its retained branch. No body velocity is changed.
            let velocity = local_slip_velocity(system, p, &bodies)?;
            let state = if self.resolved_velocity(velocity) {
                state.with_endpoint_tangential_velocity(velocity)?
            } else {
                state
            };
            self.temperature[p] = Some(TemperatureStage { law, state });
        }
        Ok(())
    }

    fn history_for_step(
        &self,
        p: usize,
        thermal_active: bool,
    ) -> Result<(LocalDamageFriction<'a>, DamageFrictionState<'a>), RigidCohesiveError> {
        if thermal_active {
            self.temperature[p]
                .map(|s| (s.law, s.state))
                .ok_or(RigidCohesiveError::Invalid(
                    "incomplete friction temperature stage",
                ))
        } else {
            Ok((self.laws[p], self.states[p]))
        }
    }

    pub(super) fn temperature_complete(&self) -> bool {
        self.temperature.iter().all(Option::is_some)
    }

    fn resolved_velocity(&self, velocity: [f64; 2]) -> bool {
        velocity[0].hypot(velocity[1]) > self.config.velocity_absolute_tolerance_m_s
    }

    fn evaluate(
        &self,
        system: &mut RigidCohesiveSystem<'a>,
        x: &[f64; MAX_DIM],
        mu: &[f64; 64],
        loads: [RigidCohesiveLoad; 2],
        work: &mut Work,
        gate: Option<&CancelGate>,
    ) -> Result<Joint, RigidCohesiveError> {
        let count = self.laws.len();
        let mut forces = [[0.0; 2]; 64];
        for i in 0..count {
            forces[i] = [x[12 + 2 * i], x[13 + 2 * i]];
        }
        let native = system.evaluate(
            core::array::from_fn(|i| x[i]),
            loads,
            Some(&forces[..count]),
            work,
            gate,
        )?;
        let mut value = Joint {
            native,
            residual: [0.0; MAX_DIM],
            tolerance: [1.0; MAX_DIM],
            slips: [0.0; 64],
            pressure: [0.0; 64],
            slip_tolerance: [0.0; 64],
            endpoint_velocity: [[0.0; 2]; 64],
        };
        let endpoint_bodies = native.receipts.map(|receipt| receipt.state_after);
        for i in 0..12 {
            value.residual[i] = native.residual[i];
            value.tolerance[i] = if i % 6 < 3 {
                system.config.force_absolute_tolerance_n
            } else {
                system.config.torque_absolute_tolerance_n_m
            } + system.config.relative_tolerance
                * x[i].abs().max(native.expected[i].abs());
        }
        for p in 0..count {
            poll(gate)?;
            let (law, state) = self.history_for_step(p, system.thermal_active)?;
            let trial = system.staged_friction[p].frame;
            let delta = trial.slip_increment_m;
            let slip = delta[0].hypot(delta[1]);
            let pressure = trial.damaged_compression_force_n;
            let cap = finite(mu[p] * pressure)?;
            let k = self.cone_scale[p];
            let q = [forces[p][0] + k * delta[0], forces[p][1] + k * delta[1]];
            let length = q[0].hypot(q[1]);
            finite(length)?;
            let projected = fs_material::cohesive::friction::project_tangent_disk(q, cap).map_err(
                crate::cohesive_contact::frictional_damage::DamageFrictionError::Material,
            )?;
            let old = state.point().jump_m();
            let jump = system.staged_sites[p]
                .as_ref()
                .expect("native evaluation staged every site")
                .point()
                .numerical()
                .material()
                .jump_m();
            let slip_tol = law.jump_tolerance(old, jump, jump)?;
            value.slips[p] = slip;
            value.pressure[p] = pressure;
            value.slip_tolerance[p] = slip_tol;
            value.endpoint_velocity[p] = local_slip_velocity(system, p, &endpoint_bodies)?;
            let moving = self.resolved_velocity(self.incoming_velocity[p])
                || self.resolved_velocity(value.endpoint_velocity[p]);
            for i in 0..2 {
                let j = 12 + 2 * p + i;
                value.residual[j] = finite(forces[p][i] - projected[i])?;
                let force_tol = law.force_tolerance(forces[p][i], cap)? * 0.0625;
                value.tolerance[j] = positive(if length <= cap {
                    // A resolved moving path must reach its kinetic disk even
                    // when its displacement is below the static resolution.
                    let resolution = if moving && slip > 0.0 {
                        slip.min(slip_tol)
                    } else {
                        slip_tol
                    };
                    force_tol.min(k * resolution * 0.0625)
                } else if slip <= 0.25 * slip_tol && (!moving || slip == 0.0) {
                    force_tol.min(k * slip_tol * 0.0625)
                } else {
                    force_tol * (k * slip / (cap + k * slip))
                })?;
            }
            balance(
                "objective friction pullback",
                trial.spatial_work_j - forces[p][0] * delta[0] - forces[p][1] * delta[1],
                system.config.objective.spatial_work_tolerance_j,
            )?;
        }
        finite(value.merit(12 + 2 * count))?;
        Ok(value)
    }

    pub(super) fn solve(
        &mut self,
        system: &mut RigidCohesiveSystem<'a>,
        loads: [RigidCohesiveLoad; 2],
        gate: Option<&CancelGate>,
    ) -> Result<(Evaluation, [Wrench; 2], usize, Work), RigidCohesiveError> {
        let count = self.laws.len();
        let dim = 12 + 2 * count;
        let mut x = [0.0; MAX_DIM];
        x[..12].copy_from_slice(&pack(system.last_wrenches));
        let mut mu = [0.0; 64];
        let mut yielded = [false; 64];
        let incoming_bodies = system.bodies.each_ref().map(|body| body.state);
        for p in 0..count {
            poll(gate)?;
            let (law, state) = self.history_for_step(p, system.thermal_active)?;
            x[12 + 2 * p..14 + 2 * p].copy_from_slice(&self.last_forces[p]);
            self.incoming_velocity[p] = local_slip_velocity(system, p, &incoming_bodies)?;
            yielded[p] = state.occupancy() == DamageFrictionRegime::Sliding
                || self.resolved_velocity(self.incoming_velocity[p]);
            mu[p] = if yielded[p] {
                law.kinetic_mu
            } else {
                law.capture_mu(state)
            };
        }
        let mut work = Work {
            evaluations: 0,
            material_evaluations: if system.thermal_active { count } else { 0 },
        };
        let mut rhs = [0.0; MAX_DIM];
        let mut correction = [0.0; MAX_DIM];
        for iteration in 0..=self.config.maximum_iterations {
            poll(gate)?;
            let value = self.evaluate(system, &x, &mu, loads, &mut work, gate)?;
            let worst = value.worst(dim);
            if worst <= 1.0 {
                let mut changed = false;
                for p in 0..count {
                    if !yielded[p]
                        && value.pressure[p] > 0.0
                        && (value.slips[p] > 0.25 * value.slip_tolerance[p]
                            || self.resolved_velocity(value.endpoint_velocity[p]))
                    {
                        yielded[p] = true;
                        let (law, _) = self.history_for_step(p, system.thermal_active)?;
                        if mu[p] != law.kinetic_mu {
                            mu[p] = law.kinetic_mu;
                            changed = true;
                        }
                    }
                }
                if !changed {
                    let mut native = value.native;
                    native.friction = self.audit(system, &x, &value, &mut work, gate)?;
                    return Ok((
                        native,
                        unpack(core::array::from_fn(|i| x[i])),
                        iteration,
                        work,
                    ));
                }
                continue;
            }
            if iteration == self.config.maximum_iterations {
                return Err(RigidCohesiveError::NotConverged {
                    scaled_residual: worst,
                    iterations: iteration,
                });
            }
            for j in 0..dim {
                let increment = 6.055454452393343e-6 * x[j].abs().max(self.coordinate_scale[j]);
                let mut high_x = x;
                high_x[j] += increment;
                let mut low_x = x;
                low_x[j] -= increment;
                let high = self.evaluate(system, &high_x, &mu, loads, &mut work, gate);
                if let Err(ref e) = high {
                    terminal(e)?;
                }
                let low = self.evaluate(system, &low_x, &mu, loads, &mut work, gate);
                if let Err(ref e) = low {
                    terminal(e)?;
                }
                for i in 0..dim {
                    let derivative = match (&high, &low) {
                        (Ok(a), Ok(b)) => (a.residual[i] - b.residual[i]) / (high_x[j] - low_x[j]),
                        (Ok(a), Err(_)) => (a.residual[i] - value.residual[i]) / (high_x[j] - x[j]),
                        (Err(_), Ok(b)) => (value.residual[i] - b.residual[i]) / (x[j] - low_x[j]),
                        (Err(e), Err(_)) => return Err(e.clone()),
                    };
                    self.jacobian[i * dim + j] =
                        finite(derivative * self.coordinate_scale[j] / value.tolerance[i])?;
                }
            }
            let scale = positive(self.jacobian.iter().map(|v| v.abs()).fold(0.0, f64::max))?;
            for v in &mut self.jacobian {
                *v /= scale;
            }
            for i in 0..dim {
                rhs[i] = 0.0;
                for row in 0..dim {
                    rhs[i] -= self.jacobian[row * dim + i]
                        * (value.residual[row] / value.tolerance[row] / scale);
                }
                for j in 0..dim {
                    self.normal[i * dim + j] = finite(
                        (0..dim)
                            .map(|r| self.jacobian[r * dim + i] * self.jacobian[r * dim + j])
                            .sum::<f64>()
                            + if i == j {
                                self.config.newton_damping
                            } else {
                                0.0
                            },
                    )?;
                }
            }
            poll(gate)?;
            self.linear
                .solve_into(&self.normal, &rhs[..dim], &mut correction[..dim])?;
            let mut fraction = 1.0;
            let mut accepted = false;
            for _ in 0..system.config.maximum_backtracks {
                let mut next = x;
                for i in 0..dim {
                    next[i] += fraction * self.coordinate_scale[i] * correction[i];
                }
                match self.evaluate(system, &next, &mu, loads, &mut work, gate) {
                    Ok(next_value)
                        if next_value.merit(dim) < value.merit(dim)
                            || next_value.worst(dim) <= 1.0 =>
                    {
                        x = next;
                        accepted = true;
                        break;
                    }
                    Err(e) => terminal(&e)?,
                    _ => {}
                }
                fraction *= 0.5;
            }
            if !accepted {
                // At a short timestep the interior-disk Jacobian can see
                // almost prescribed slip and cannot reach its distant cone.
                // One joint boundary predictor changes only effort unknowns;
                // the complete original bodies and every site are reevaluated.
                let mut next = x;
                next[..12].copy_from_slice(&value.native.expected);
                // Restore the current physical trial after Jacobian/line-search
                // scratch before reading its actual local direction.
                self.evaluate(system, &x, &mu, loads, &mut work, gate)?;
                for p in 0..count {
                    let delta = system.staged_friction[p].frame.slip_increment_m;
                    let slip = delta[0].hypot(delta[1]);
                    if slip > 0.0 {
                        for i in 0..2 {
                            next[12 + 2 * p + i] = mu[p] * value.pressure[p] * delta[i] / slip;
                        }
                    }
                }
                match self.evaluate(system, &next, &mu, loads, &mut work, gate) {
                    Ok(next_value)
                        if next_value.merit(dim) < value.merit(dim)
                            || next_value.worst(dim) <= 1.0 =>
                    {
                        x = next;
                        accepted = true;
                    }
                    Err(e) => terminal(&e)?,
                    _ => {}
                }
            }
            if !accepted {
                return Err(RigidCohesiveError::NotConverged {
                    scaled_residual: worst,
                    iterations: iteration + 1,
                });
            }
        }
        Err(RigidCohesiveError::NotConverged {
            scaled_residual: f64::INFINITY,
            iterations: self.config.maximum_iterations,
        })
    }

    fn audit(
        &mut self,
        system: &mut RigidCohesiveSystem<'a>,
        x: &[f64; MAX_DIM],
        value: &Joint,
        work: &mut Work,
        gate: Option<&CancelGate>,
    ) -> Result<f64, RigidCohesiveError> {
        let mut loss = 0.0;
        for p in 0..self.laws.len() {
            poll(gate)?;
            let (law, state) = self.history_for_step(p, system.thermal_active)?;
            if work.material_evaluations >= system.config.maximum_material_evaluations {
                return Err(RigidCohesiveError::Budget(
                    "original friction audit evaluations",
                ));
            }
            work.material_evaluations += 1;
            let point = system.staged_sites[p]
                .as_ref()
                .expect("native evaluation staged every site")
                .numerical();
            let mut applied = point.mean_traction_pa.map(|f| f * law.area);
            for i in 0..2 {
                applied[i + 1] += x[12 + 2 * p + i];
            }
            let old_velocity = self.incoming_velocity[p];
            let new_velocity = value.endpoint_velocity[p];
            let incoming_moving = self.resolved_velocity(old_velocity);
            let endpoint_moving = self.resolved_velocity(new_velocity);
            let tolerance = self.config.velocity_absolute_tolerance_m_s;
            if value.pressure[p] > 0.0
                && incoming_moving
                && endpoint_moving
                && old_velocity[0] * new_velocity[0] + old_velocity[1] * new_velocity[1]
                    < -tolerance * tolerance
            {
                return Err(RigidCohesiveError::Invalid(
                    "unresolved friction stop reverses endpoint slip; subdivide the step",
                ));
            }
            let regime = if value.pressure[p] == 0.0 {
                DamageFrictionRegime::Inactive
            } else if value.slips[p] <= 0.25 * value.slip_tolerance[p]
                && !incoming_moving
                && !endpoint_moving
            {
                DamageFrictionRegime::Sticking
            } else {
                DamageFrictionRegime::Sliding
            };
            let audit = law.audit(
                state,
                point.point.material().jump_m(),
                system.dt_s,
                applied,
                regime,
            )?;
            if audit.state.point() != point.point.material() {
                return Err(RigidCohesiveError::Invalid(
                    "original cohesive/friction successors differ",
                ));
            }
            let frame = &mut system.staged_friction[p].frame;
            frame.regime = regime;
            frame.friction_dissipation_j = audit.friction_dissipation_j;
            balance(
                "actual objective friction work",
                frame.spatial_work_j - audit.friction_dissipation_j,
                system.config.objective.spatial_work_tolerance_j
                    + system.config.energy_relative_tolerance * audit.friction_dissipation_j,
            )?;
            self.staged[p] = audit.state.with_resolved_endpoint_tangential_velocity(
                new_velocity,
                self.config.velocity_absolute_tolerance_m_s,
            )?;
            self.staged_forces[p] = frame.local_force_n;
            loss += audit.friction_dissipation_j;
        }
        finite(loss)
    }

    pub(super) fn publish(&mut self, thermal_active: bool) {
        if thermal_active {
            for p in 0..self.laws.len() {
                let stage =
                    self.temperature[p].expect("complete friction temperature stage checked");
                self.laws[p] = stage.law;
                self.sources[p].source = stage.law.friction;
            }
        }
        self.states.copy_from_slice(&self.staged);
        self.last_forces = self.staged_forces;
    }
}

fn zeros(n: usize) -> Result<Vec<f64>, RigidCohesiveError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(n)
        .map_err(|_| RigidCohesiveError::Budget("joint friction workspace allocation"))?;
    result.resize(n, 0.0);
    Ok(result)
}

fn local_slip_velocity(
    system: &RigidCohesiveSystem<'_>,
    p: usize,
    states: &[RigidBodyState; 2],
) -> Result<[f64; 2], RigidCohesiveError> {
    let site = &system.sites[p];
    let a = states[0].point_kinematics(system.bodies[0].properties, site.arms_body_m[0])?;
    let b = states[1].point_kinematics(system.bodies[1].properties, site.arms_body_m[1])?;
    let velocity = b.point_velocity_world.sub(a.point_velocity_world).sub(
        a.angular_velocity_world
            .cross(separation(states, site.arms_body_m)?),
    );
    let local =
        compose(states[0].pose().orientation(), site.frame_in_a)?.rotate_world_to_body(velocity);
    finite_vector(local)?;
    Ok([local.y, local.z])
}

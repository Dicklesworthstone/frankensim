//! Simultaneous damaged-area friction over the ORIGINAL finite bonded patch.
//!
//! Unknowns are all 3N total site forces. Every jump is `g_free - S F`, with
//! the complete original cross-site/cross-axis response. Original local damage,
//! cohesive/viscous work and the analytic damaged-compression measure remain
//! with the shared local law; quadrature area is applied exactly once.
//!
//! The tangential residual is `f - project_disk(f + k delta_g, mu P)`.
//! Here k only scales a cone residual: at its zero the EXACT rigid Coulomb law
//! holds, without an elastic friction spring. Joint capture solves use each
//! site's static coefficient (kinetic if already moving); sites whose capture
//! solution slides enter their kinetic branch and the WHOLE patch is resolved.
//! All switches and Newton corrections share explicit work/evaluation ceilings.
//!
//! Rank-deficient reciprocal compliance is admitted, including coincident
//! quadrature and rigid bodies with fewer modes than contact coordinates.
//! Damped least-squares Newton corrections use the original prepared LU and
//! area-weighted force coordinates. This is numerical reaction selection for
//! potentially indeterminate static forces, not additional physical stiffness,
//! a uniqueness theorem or a stability certificate. Physical force/work gates
//! check the original unregularized equations. Refine time and quadrature.

use super::{CohesivePatchError, MAX_COHESIVE_PATCH_POINTS, filled};
use crate::cohesive_contact::frictional_damage::local::LocalDamageFriction;
use crate::cohesive_contact::frictional_damage::{
    DamageFrictionAudit, DamageFrictionConfig, DamageFrictionError, DamageFrictionRegime,
    DamageFrictionState,
};
use crate::cohesive_contact::reclosure::{ReclosureActivation, ResidualCoulombSource};
use fs_blake3::DomainHasher;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_matdb::MaterialStateId;
use fs_material::cohesive::CohesivePoint;
use fs_material::cohesive::source::ResolvedBilinearCohesive;

#[derive(Clone, Copy)]
pub struct DamageFrictionPatchSite<'a> {
    pub cohesive: &'a ResolvedBilinearCohesive,
    pub friction: &'a ResidualCoulombSource,
    pub area_m2: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DamageFrictionPatchConfig {
    pub max_points: usize,
    /// One total point-law budget, including every Jacobian and final audit.
    /// max_friction_iterations bounds ALL joint corrections/branch switches.
    pub friction: DamageFrictionConfig,
    /// Dimensionless regularization of the scaled Newton NORMAL EQUATIONS,
    /// in (0,1). Changes only numerical corrections, never physical compliance,
    /// the cone, storage or accepted force/work equations.
    pub newton_damping: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DamageFrictionPatchError {
    Invalid(&'static str),
    Setup(CohesivePatchError),
    Point {
        point: usize,
        source: DamageFrictionError,
    },
    Factor(LuWorkspaceError),
    Cancelled,
    Budget,
    NotConverged {
        residual_ratio: f64,
        iterations: usize,
    },
}
impl core::fmt::Display for DamageFrictionPatchError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint damaged cohesive patch refused: {self:?}")
    }
}
impl std::error::Error for DamageFrictionPatchError {}
impl From<CohesivePatchError> for DamageFrictionPatchError {
    fn from(e: CohesivePatchError) -> Self {
        Self::Setup(e)
    }
}
impl From<LuWorkspaceError> for DamageFrictionPatchError {
    fn from(e: LuWorkspaceError) -> Self {
        Self::Factor(e)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct DamageFrictionPatchTrial<'s, 'a> {
    pub jump_m: &'s [[f64; 3]],
    pub mean_force_n: &'s [[f64; 3]],
    pub regimes: &'s [DamageFrictionRegime],
    pub compatibility_tolerance_m: &'s [f64],
    pub iterations: usize,
    pub kinetic_switches: usize,
    pub material_evaluations: usize,
    laws: &'s [LocalDamageFriction<'a>],
}
impl<'a> DamageFrictionPatchTrial<'_, 'a> {
    /// Same source-bound local audit on actual receiving-body motion; this
    /// immutable view remains valid while its solver's scratch is borrowed.
    pub fn audit_point(
        &self,
        p: usize,
        previous: DamageFrictionState<'a>,
        jump: [f64; 3],
        dt: f64,
        applied: [f64; 3],
        regime: DamageFrictionRegime,
    ) -> Result<DamageFrictionAudit<'a>, DamageFrictionPatchError> {
        self.laws
            .get(p)
            .ok_or(DamageFrictionPatchError::Invalid("site index"))?
            .audit(previous, jump, dt, applied, regime)
            .map_err(|e| local_error(p, e))
    }
}

struct Evaluation {
    force: Vec<[f64; 3]>,
    jump: Vec<[f64; 3]>,
    delta: Vec<[f64; 3]>,
    material: Vec<[f64; 3]>,
    pressure: Vec<f64>,
    residual: Vec<[f64; 3]>,
    force_tolerance: Vec<[f64; 3]>,
    jump_tolerance: Vec<f64>,
    projected: Vec<bool>,
    regimes: Vec<DamageFrictionRegime>,
}
impl Evaluation {
    fn new(n: usize) -> Result<Self, DamageFrictionPatchError> {
        Ok(Self {
            force: filled(n, [0.0; 3])?,
            jump: filled(n, [0.0; 3])?,
            delta: filled(n, [0.0; 3])?,
            material: filled(n, [0.0; 3])?,
            pressure: filled(n, 0.0)?,
            residual: filled(n, [0.0; 3])?,
            force_tolerance: filled(n, [0.0; 3])?,
            jump_tolerance: filled(n, 0.0)?,
            projected: filled(n, false)?,
            regimes: filled(n, DamageFrictionRegime::Inactive)?,
        })
    }
    fn worst(&self) -> f64 {
        self.residual
            .iter()
            .zip(&self.force_tolerance)
            .flat_map(|(r, t)| (0..3).map(move |i| r[i].abs() / t[i]))
            .fold(0.0, f64::max)
    }
    fn merit(&self) -> f64 {
        self.residual
            .iter()
            .zip(&self.force_tolerance)
            .map(|(r, t)| (0..3).map(|i| (r[i] / t[i]).powi(2)).sum::<f64>())
            .sum()
    }
}

/// One joint solve, never a collection of disconnected site load lines.
pub struct DamageFrictionPatchLoadLine<'a> {
    laws: Vec<LocalDamageFriction<'a>>,
    compliance: Vec<f64>,
    config: DamageFrictionPatchConfig,
    cone_scale: Vec<f64>,
    coordinate_scale: Vec<f64>,
    mu: Vec<f64>,
    yielded: Vec<bool>,
    current: Evaluation,
    candidate: Evaluation,
    jacobian: Vec<f64>,
    normal: Vec<f64>,
    rhs: Vec<f64>,
    correction: Vec<f64>,
    solve: LuWorkspace,
}
impl<'a> DamageFrictionPatchLoadLine<'a> {
    pub fn new(
        sites: Vec<DamageFrictionPatchSite<'a>>,
        compliance: Vec<f64>,
        config: DamageFrictionPatchConfig,
    ) -> Result<Self, DamageFrictionPatchError> {
        let n = sites.len();
        if !(1..=MAX_COHESIVE_PATCH_POINTS).contains(&config.max_points)
            || n == 0
            || n > config.max_points
            || config.friction.contact.max_material_evaluations < 2 * n
            || !config.newton_damping.is_finite()
            || config.newton_damping <= 0.0
            || config.newton_damping >= 1.0
        {
            return Err(DamageFrictionPatchError::Invalid(
                "bounded nonempty patch and explicit numerical damping required",
            ));
        }
        let dim = 3 * n;
        if compliance.len() != dim * dim {
            return Err(DamageFrictionPatchError::Invalid(
                "complete 3N by 3N compliance required",
            ));
        }
        let mut normal = filled(dim * dim, 0.0)?;
        positive_semidefinite(&compliance, &mut normal, dim)?;
        let total_area = sites.iter().map(|s| s.area_m2).sum::<f64>();
        positive(total_area)?;
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.damage-friction-patch.v1");
        for site in &sites {
            hash.update(site.cohesive.resolved().identity().as_bytes());
            hash.update(site.friction.resolved().identity().as_bytes());
            hash.update(&site.area_m2.to_bits().to_le_bytes());
        }
        for x in &compliance {
            hash.update(&x.to_bits().to_le_bytes());
        }
        for x in [
            config.newton_damping,
            config.friction.contact.jump_absolute_tolerance_m,
            config.friction.contact.jump_relative_tolerance,
            config.friction.contact.maximum_jump_increment_m,
            config.friction.contact.maximum_force_n,
            config.friction.contact.material_energy_tolerance_j_m2,
            config.friction.force_absolute_tolerance_n,
            config.friction.force_relative_tolerance,
            config.friction.work_absolute_tolerance_j,
            config.friction.work_relative_tolerance,
        ] {
            hash.update(&x.to_bits().to_le_bytes());
        }
        for v in [
            config.max_points,
            config.friction.max_friction_iterations,
            config.friction.contact.max_iterations,
            config.friction.contact.max_backtracks,
            config.friction.contact.max_material_evaluations,
        ] {
            hash.update(&(v as u64).to_le_bytes());
        }
        let identity = hash.finalize();
        let mut laws = Vec::new();
        laws.try_reserve_exact(n)
            .map_err(|_| DamageFrictionPatchError::Invalid("patch allocation"))?;
        let mut cone_scale = filled(n, 0.0)?;
        let mut coordinate_scale = filled(n, 0.0)?;
        for (p, site) in sites.into_iter().enumerate() {
            let mut hash =
                DomainHasher::new("org.frankensim.fs-couple.damage-friction-patch.site.v1");
            hash.update(identity.as_bytes());
            hash.update(&(p as u64).to_le_bytes());
            laws.push(
                LocalDamageFriction::new(
                    site.cohesive,
                    site.friction,
                    site.area_m2,
                    config.friction,
                    hash.finalize(),
                )
                .map_err(|e| local_error(p, e))?,
            );
            let response = (1..3)
                .map(|axis| {
                    compliance[(3 * p + axis) * dim..(3 * p + axis + 1) * dim]
                        .iter()
                        .map(|x| x.abs())
                        .sum::<f64>()
                })
                .fold(0.0, f64::max);
            // Zero mechanical tangent response remains valid: the cohesive
            // onset gives a finite residual scale without altering the law.
            let material_length = site.cohesive.law().onset_separation_m();
            let material_force = site.area_m2 * site.cohesive.law().strength_pa();
            cone_scale[p] = positive(1.0 / response.max(material_length / material_force))?;
            coordinate_scale[p] = positive((site.area_m2 / total_area).sqrt())?;
        }
        Ok(Self {
            laws,
            compliance,
            config,
            cone_scale,
            coordinate_scale,
            mu: filled(n, 0.0)?,
            yielded: filled(n, false)?,
            current: Evaluation::new(n)?,
            candidate: Evaluation::new(n)?,
            jacobian: filled(dim * dim, 0.0)?,
            normal,
            rhs: filled(dim, 0.0)?,
            correction: filled(dim, 0.0)?,
            solve: LuWorkspace::new(dim)?,
        })
    }
    pub fn point_count(&self) -> usize {
        self.laws.len()
    }
    pub fn compliance_m_per_n(&self) -> &[f64] {
        &self.compliance
    }
    pub const fn config(&self) -> DamageFrictionPatchConfig {
        self.config
    }
    pub fn friction_source(&self, p: usize) -> Option<&'a ResidualCoulombSource> {
        self.laws.get(p).map(|l| l.friction)
    }
    pub fn adopt_point(
        &self,
        p: usize,
        point: CohesivePoint<'a>,
        receiving: [&MaterialStateId; 2],
        activation: ReclosureActivation,
    ) -> Result<DamageFrictionState<'a>, DamageFrictionPatchError> {
        self.laws
            .get(p)
            .ok_or(DamageFrictionPatchError::Invalid("site index"))?
            .adopt(point, receiving, activation)
            .map_err(|e| local_error(p, e))
    }
    pub fn audit_point(
        &self,
        p: usize,
        previous: DamageFrictionState<'a>,
        jump: [f64; 3],
        dt: f64,
        applied: [f64; 3],
        regime: DamageFrictionRegime,
    ) -> Result<DamageFrictionAudit<'a>, DamageFrictionPatchError> {
        self.laws
            .get(p)
            .ok_or(DamageFrictionPatchError::Invalid("site index"))?
            .audit(previous, jump, dt, applied, regime)
            .map_err(|e| local_error(p, e))
    }

    pub fn trial(
        &mut self,
        previous: &[DamageFrictionState<'a>],
        free: &[[f64; 3]],
        dt: f64,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<DamageFrictionPatchTrial<'_, 'a>, DamageFrictionPatchError> {
        poll(&mut keep_going)?;
        let n = self.laws.len();
        let dim = 3 * n;
        if previous.len() != n
            || free.len() != n
            || !dt.is_finite()
            || dt <= 0.0
            || free.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(DamageFrictionPatchError::Invalid(
                "complete original states, finite free jumps and positive timestep required",
            ));
        }
        for p in 0..n {
            self.laws[p]
                .require_state(previous[p])
                .map_err(|e| local_error(p, e))?;
            self.mu[p] = self.laws[p].capture_mu(previous[p]);
            self.yielded[p] = previous[p].occupancy() == DamageFrictionRegime::Sliding;
            self.current.force[p] = previous[p]
                .point()
                .traction_pa()
                .map(|x| x * self.laws[p].area);
        }
        let mut evaluations = 0;
        let mut switches = 0;
        let mut corrections = 0;
        for iteration in 0..=self.config.friction.max_friction_iterations {
            let ctx = Context {
                laws: &self.laws,
                compliance: &self.compliance,
                cone_scale: &self.cone_scale,
                config: self.config,
                previous,
                free,
                dt,
                mu: &self.mu,
            };
            ctx.evaluate(&mut self.current, &mut evaluations, &mut keep_going)?;
            let worst = self.current.worst();
            if worst <= 1.0 {
                let mut changed = false;
                for p in 0..n {
                    let slip = self.current.delta[p][1].hypot(self.current.delta[p][2]);
                    if !self.yielded[p]
                        && self.current.pressure[p] > 0.0
                        && slip > 0.25 * self.current.jump_tolerance[p]
                    {
                        self.yielded[p] = true;
                        if self.mu[p] != self.laws[p].kinetic_mu {
                            self.mu[p] = self.laws[p].kinetic_mu;
                            switches += 1;
                            changed = true;
                        }
                    }
                }
                if changed {
                    if iteration == self.config.friction.max_friction_iterations {
                        return Err(DamageFrictionPatchError::NotConverged {
                            residual_ratio: worst,
                            iterations: iteration,
                        });
                    }
                    continue;
                }
                for p in 0..n {
                    spend(&mut evaluations, self.config, &mut keep_going)?;
                    let regime = if self.current.pressure[p] == 0.0 {
                        DamageFrictionRegime::Inactive
                    } else if self.current.delta[p][1].hypot(self.current.delta[p][2])
                        <= 0.25 * self.current.jump_tolerance[p]
                    {
                        DamageFrictionRegime::Sticking
                    } else {
                        DamageFrictionRegime::Sliding
                    };
                    self.laws[p]
                        .audit(
                            previous[p],
                            self.current.jump[p],
                            dt,
                            self.current.force[p],
                            regime,
                        )
                        .map_err(|e| local_error(p, e))?;
                    self.current.regimes[p] = regime;
                }
                poll(&mut keep_going)?;
                return Ok(DamageFrictionPatchTrial {
                    jump_m: &self.current.jump,
                    mean_force_n: &self.current.force,
                    regimes: &self.current.regimes,
                    compatibility_tolerance_m: &self.current.jump_tolerance,
                    iterations: iteration,
                    kinetic_switches: switches,
                    material_evaluations: evaluations,
                    laws: &self.laws,
                });
            }
            if iteration == self.config.friction.max_friction_iterations
                || corrections == self.config.friction.contact.max_iterations
            {
                return Err(DamageFrictionPatchError::NotConverged {
                    residual_ratio: worst,
                    iterations: iteration,
                });
            }
            ctx.jacobian(
                &self.current,
                &mut self.jacobian,
                &mut evaluations,
                &mut keep_going,
            )?;
            // Scale rows to the actual local gates and columns to sqrt(area).
            // Minimum change of these coordinates is an area-weighted traction
            // choice when static reaction distributions are nonunique.
            let mut scale: f64 = 0.0;
            for i in 0..dim {
                for j in 0..dim {
                    let v = &mut self.jacobian[i * dim + j];
                    *v = finite(
                        *v * self.coordinate_scale[j / 3]
                            / self.current.force_tolerance[i / 3][i % 3],
                    )?;
                    scale = scale.max(v.abs());
                }
            }
            positive(scale)?;
            for v in &mut self.jacobian {
                *v /= scale;
            }
            for i in 0..dim {
                self.rhs[i] = 0.0;
                for row in 0..dim {
                    self.rhs[i] -= self.jacobian[row * dim + i]
                        * (self.current.residual[row / 3][row % 3]
                            / self.current.force_tolerance[row / 3][row % 3]
                            / scale);
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
            poll(&mut keep_going)?;
            self.solve
                .solve_into(&self.normal, &self.rhs, &mut self.correction)?;
            let merit = self.current.merit();
            let mut step = 1.0;
            let mut accepted = false;
            for _ in 0..self.config.friction.contact.max_backtracks {
                for p in 0..n {
                    for i in 0..3 {
                        self.candidate.force[p][i] = self.current.force[p][i]
                            + step * self.coordinate_scale[p] * self.correction[3 * p + i];
                    }
                }
                match ctx.evaluate(&mut self.candidate, &mut evaluations, &mut keep_going) {
                    Ok(()) if self.candidate.merit() < merit || self.candidate.worst() <= 1.0 => {
                        core::mem::swap(&mut self.current, &mut self.candidate);
                        accepted = true;
                        break;
                    }
                    Ok(()) => {}
                    Err(e) if outside(&e) => {}
                    Err(e) => return Err(e),
                }
                step *= 0.5;
            }
            if !accepted {
                // A kinematically prescribed slip can lie in a zero-compliance
                // direction: the interior-cone Jacobian then cannot reach its
                // boundary. Try the exact local cone boundary as ONE joint
                // force predictor, reevaluating the complete original S.
                // This never accepts disconnected per-site mechanics.
                for p in 0..n {
                    self.candidate.force[p] = self.current.material[p];
                    let slip = self.current.delta[p][1].hypot(self.current.delta[p][2]);
                    if slip > 0.0 {
                        let cap = self.mu[p] * self.current.pressure[p];
                        for i in 1..3 {
                            self.candidate.force[p][i] += cap * self.current.delta[p][i] / slip;
                        }
                    }
                }
                match ctx.evaluate(&mut self.candidate, &mut evaluations, &mut keep_going) {
                    Ok(()) if self.candidate.merit() < merit || self.candidate.worst() <= 1.0 => {
                        core::mem::swap(&mut self.current, &mut self.candidate);
                        accepted = true;
                    }
                    Ok(()) => {}
                    Err(e) if outside(&e) => {}
                    Err(e) => return Err(e),
                }
            }
            if !accepted {
                return Err(DamageFrictionPatchError::NotConverged {
                    residual_ratio: worst,
                    iterations: iteration + 1,
                });
            }
            corrections += 1;
        }
        unreachable!("bounded joint solve")
    }
}

struct Context<'s, 'a> {
    laws: &'s [LocalDamageFriction<'a>],
    compliance: &'s [f64],
    cone_scale: &'s [f64],
    config: DamageFrictionPatchConfig,
    previous: &'s [DamageFrictionState<'a>],
    free: &'s [[f64; 3]],
    dt: f64,
    mu: &'s [f64],
}
impl Context<'_, '_> {
    fn local(
        &self,
        p: usize,
        jump: [f64; 3],
        evaluations: &mut usize,
        keep: &mut impl FnMut() -> bool,
    ) -> Result<([f64; 3], f64), DamageFrictionPatchError> {
        let c = self.config.friction.contact;
        let old = self.previous[p].point().jump_m();
        if (0..3)
            .any(|i| !jump[i].is_finite() || (jump[i] - old[i]).abs() > c.maximum_jump_increment_m)
        {
            return Err(local_error(
                p,
                DamageFrictionError::Range("site jump increment"),
            ));
        }
        spend(evaluations, self.config, keep)?;
        let (_, mean, pressure) = self.previous[p]
            .point()
            .trial_over_time_with_damage_weighted_compression(
                jump,
                self.dt,
                c.material_energy_tolerance_j_m2,
            )
            .map_err(|e| local_error(p, DamageFrictionError::Material(e)))?;
        let force = mean.map(|x| x * self.laws[p].area);
        if force
            .iter()
            .any(|x| !x.is_finite() || x.abs() > c.maximum_force_n)
        {
            return Err(local_error(
                p,
                DamageFrictionError::Range("site cohesive force"),
            ));
        }
        Ok((force, finite(self.laws[p].area * pressure)?))
    }
    fn evaluate(
        &self,
        out: &mut Evaluation,
        evaluations: &mut usize,
        keep: &mut impl FnMut() -> bool,
    ) -> Result<(), DamageFrictionPatchError> {
        let n = self.laws.len();
        let dim = 3 * n;
        for p in 0..n {
            poll(keep)?;
            let old = self.previous[p].point().jump_m();
            for i in 0..3 {
                if !out.force[p][i].is_finite()
                    || out.force[p][i].abs() > self.config.friction.contact.maximum_force_n
                {
                    return Err(local_error(
                        p,
                        DamageFrictionError::Range("site total force"),
                    ));
                }
                let response = (0..dim)
                    .map(|j| self.compliance[(3 * p + i) * dim + j] * out.force[j / 3][j % 3])
                    .sum::<f64>();
                out.delta[p][i] = finite((self.free[p][i] - old[i]) - response)?;
                out.jump[p][i] = finite(old[i] + out.delta[p][i])?;
            }
            let (material, pressure) = self.local(p, out.jump[p], evaluations, keep)?;
            out.material[p] = material;
            out.pressure[p] = pressure;
            out.jump_tolerance[p] = self.laws[p]
                .jump_tolerance(old, self.free[p], out.jump[p])
                .map_err(|e| local_error(p, e))?;
            let f = [out.force[p][1] - material[1], out.force[p][2] - material[2]];
            let q = [
                f[0] + self.cone_scale[p] * out.delta[p][1],
                f[1] + self.cone_scale[p] * out.delta[p][2],
            ];
            let cap = finite(self.mu[p] * pressure)?;
            let length = q[0].hypot(q[1]);
            out.projected[p] = length > cap;
            let projected = if length > cap {
                q.map(|x| cap * (x / length))
            } else {
                q
            };
            out.residual[p] = [
                out.force[p][0] - material[0],
                f[0] - projected[0],
                f[1] - projected[1],
            ];
            let slip = out.delta[p][1].hypot(out.delta[p][2]);
            for i in 0..3 {
                let tol = self.laws[p]
                    .force_tolerance(out.force[p][i], material[i])
                    .map_err(|e| local_error(p, e))?
                    * 0.0625;
                out.force_tolerance[p][i] =
                    if i > 0 && (!out.projected[p] || slip <= 0.25 * out.jump_tolerance[p]) {
                        tol.min(self.cone_scale[p] * out.jump_tolerance[p] * 0.0625)
                    } else if i > 0 {
                        // Near stopping, disk residuals attenuate directional force
                        // errors by k*slip/(cap+k*slip). Scale the gate accordingly;
                        // a small cone residual alone must not mask a wrong force.
                        let motion = self.cone_scale[p] * slip;
                        tol * (motion / (cap + motion))
                    } else {
                        tol
                    };
                finite(out.residual[p][i])?;
                positive(out.force_tolerance[p][i])?;
            }
        }
        finite(out.merit())?;
        Ok(())
    }
    fn jacobian(
        &self,
        current: &Evaluation,
        out: &mut [f64],
        evaluations: &mut usize,
        keep: &mut impl FnMut() -> bool,
    ) -> Result<(), DamageFrictionPatchError> {
        let dim = 3 * self.laws.len();
        for p in 0..self.laws.len() {
            let old = self.previous[p].point().jump_m();
            let x = current.jump[p];
            let law = self.previous[p].point().law();
            let scales = law
                .mixed_mode_parameters()
                .map_or([law.onset_separation_m(); 3], |v| {
                    [
                        v.normal_strength_pa / v.normal_stiffness_pa_m,
                        v.tangential_strength_pa / v.tangential_stiffness_pa_m,
                        v.tangential_strength_pa / v.tangential_stiffness_pa_m,
                    ]
                });
            let mut dmaterial = [[0.0; 3]; 3];
            let mut dpressure = [0.0; 3];
            for axis in 0..3 {
                let h = 2.0_f64.powi(-20)
                    * scales[axis]
                        .max(old[axis].abs())
                        .max(x[axis].abs())
                        .max(self.config.friction.contact.jump_absolute_tolerance_m);
                let mut plus = x;
                plus[axis] += h;
                let mut minus = x;
                minus[axis] -= h;
                let a = self.local(p, plus, evaluations, keep);
                let b = self.local(p, minus, evaluations, keep);
                let (a, b, width) = match (a, b) {
                    (Ok(a), Ok(b)) => (a, b, plus[axis] - minus[axis]),
                    (Ok(a), Err(e)) if outside(&e) => (
                        a,
                        (current.material[p], current.pressure[p]),
                        plus[axis] - x[axis],
                    ),
                    (Err(e), Ok(b)) if outside(&e) => (
                        (current.material[p], current.pressure[p]),
                        b,
                        x[axis] - minus[axis],
                    ),
                    (Err(e), _) | (_, Err(e)) => return Err(e),
                };
                positive(width)?;
                for i in 0..3 {
                    dmaterial[i][axis] = finite((a.0[i] - b.0[i]) / width)?;
                }
                dpressure[axis] = finite((a.1 - b.1) / width)?;
            }
            let f = [
                current.force[p][1] - current.material[p][1],
                current.force[p][2] - current.material[p][2],
            ];
            let q = [
                f[0] + self.cone_scale[p] * current.delta[p][1],
                f[1] + self.cone_scale[p] * current.delta[p][2],
            ];
            let length = q[0].hypot(q[1]);
            let cap = self.mu[p] * current.pressure[p];
            for j in 0..dim {
                let dg =
                    core::array::from_fn::<_, 3, _>(|i| -self.compliance[(3 * p + i) * dim + j]);
                let df = core::array::from_fn::<_, 3, _>(|i| {
                    (if j == 3 * p + i { 1.0 } else { 0.0 }) - dot(dmaterial[i], dg)
                });
                out[(3 * p) * dim + j] = df[0];
                let dq = [
                    df[1] + self.cone_scale[p] * dg[1],
                    df[2] + self.cone_scale[p] * dg[2],
                ];
                if current.projected[p] {
                    let direction = q.map(|x| x / length);
                    let parallel = direction[0] * dq[0] + direction[1] * dq[1];
                    let dc = self.mu[p] * dot(dpressure, dg);
                    for i in 0..2 {
                        out[(3 * p + i + 1) * dim + j] = finite(
                            df[i + 1]
                                - dc * direction[i]
                                - (cap / length) * (dq[i] - direction[i] * parallel),
                        )?;
                    }
                } else {
                    for i in 1..3 {
                        out[(3 * p + i) * dim + j] = finite(-self.cone_scale[p] * dg[i])?;
                    }
                }
            }
        }
        Ok(())
    }
}

fn positive_semidefinite(
    s: &[f64],
    scratch: &mut [f64],
    n: usize,
) -> Result<(), DamageFrictionPatchError> {
    let scale = s.iter().map(|x| x.abs()).fold(0.0, f64::max);
    if !scale.is_finite() || s.iter().any(|x| !x.is_finite()) {
        return Err(DamageFrictionPatchError::Invalid(
            "finite mechanical compliance required",
        ));
    }
    if scale == 0.0 {
        return Ok(());
    }
    let tolerance = 128.0 * f64::EPSILON * n as f64;
    for i in 0..n {
        for j in 0..n {
            if (s[i * n + j] - s[j * n + i]).abs() > 64.0 * f64::EPSILON * scale {
                return Err(DamageFrictionPatchError::Invalid(
                    "reciprocal mechanical compliance required",
                ));
            }
            scratch[i * n + j] = s[i * n + j] / scale;
        }
    }
    // Pivoted Schur test, allowing a zero trailing block. This is a floating
    // admission gate, not a certified eigenvalue enclosure or regularization.
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|a, b| scratch[a * n + *a].total_cmp(&scratch[b * n + *b]))
            .unwrap();
        if scratch[pivot * n + pivot] <= tolerance {
            if (k..n).any(|i| (k..n).any(|j| scratch[i * n + j].abs() > tolerance)) {
                return Err(DamageFrictionPatchError::Invalid(
                    "nonnegative mechanical response required",
                ));
            }
            break;
        }
        for j in 0..n {
            scratch.swap(k * n + j, pivot * n + j);
        }
        for i in 0..n {
            scratch.swap(i * n + k, i * n + pivot);
        }
        for i in k + 1..n {
            for j in k + 1..n {
                scratch[i * n + j] -= scratch[i * n + k] * scratch[k * n + j] / scratch[k * n + k];
            }
        }
    }
    Ok(())
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn local_error(point: usize, source: DamageFrictionError) -> DamageFrictionPatchError {
    DamageFrictionPatchError::Point { point, source }
}
fn outside(e: &DamageFrictionPatchError) -> bool {
    matches!(
        e,
        DamageFrictionPatchError::Point {
            source: DamageFrictionError::Range(_)
                | DamageFrictionError::Material(
                    fs_material::cohesive::CohesiveError::SeparationLimit
                ),
            ..
        }
    )
}
fn finite(v: f64) -> Result<f64, DamageFrictionPatchError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(DamageFrictionPatchError::Invalid(
            "nonfinite patch arithmetic",
        ))
    }
}
fn positive(v: f64) -> Result<f64, DamageFrictionPatchError> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(DamageFrictionPatchError::Invalid(
            "positive finite patch scale required",
        ))
    }
}
fn poll(keep: &mut impl FnMut() -> bool) -> Result<(), DamageFrictionPatchError> {
    if keep() {
        Ok(())
    } else {
        Err(DamageFrictionPatchError::Cancelled)
    }
}
fn spend(
    n: &mut usize,
    c: DamageFrictionPatchConfig,
    keep: &mut impl FnMut() -> bool,
) -> Result<(), DamageFrictionPatchError> {
    poll(keep)?;
    if *n >= c.friction.contact.max_material_evaluations {
        return Err(DamageFrictionPatchError::Budget);
    }
    *n += 1;
    Ok(())
}

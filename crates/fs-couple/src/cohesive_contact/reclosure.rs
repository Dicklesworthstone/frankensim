//! Residual Coulomb contact after complete failure of the original bond.
//!
//! The original cohesive point retains its fracture history and undamaged
//! compression potential. Only fully failed points can enter this owner.
//! Friction uses separately sourced residual coefficients and the exact mean
//! compression over the declared straight jump path. A wholly open path has
//! no pressure or friction. Opening clears contact occupancy; a later closure
//! enters the fresh static/kinetic Coulomb selection. No rebonding is inferred.
//!
//! All entries of the original reciprocal positive-definite compliance enter
//! the coupled normal/two-tangent solve. A conservative normal/friction feedback
//! screen admits a monotone scalar sliding continuation. SPD alone is not a
//! Coulomb uniqueness theorem. Floating-point residual gates are not interval
//! certificates, and straight-path work closure does not resolve intrastep
//! reversals or prove trajectory accuracy. A `Sticking` step has zero net slip;
//! a receiving dynamic owner must retain kinetic occupancy when its actual
//! endpoint velocity remains nonzero. Refine the mechanical timestep.

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialStateId;
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_material::cohesive::{CohesiveError, CohesivePoint, CohesiveRateStep};
use fs_material::state_point::MaterialStatePointError;
use fs_tribo::TriboError;

pub mod source;
pub use source::{ResidualCoulombSource, resolve_residual_coulomb};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReclosureRegime {
    Open,
    Sticking,
    Sliding,
}

/// Initial occupancy cannot be inferred from a failed material point alone.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReclosureActivation {
    AtRest,
    Sliding,
}

/// Explicit mechanical, constitutive, work and nonlinear-solve admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReclosureConfig {
    pub max_iterations: usize,
    /// Includes analytical compression probes and the final original-law trial.
    pub max_evaluations: usize,
    pub jump_absolute_tolerance_m: f64,
    pub jump_relative_tolerance: f64,
    pub maximum_jump_increment_m: f64,
    pub maximum_force_n: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
    pub material_energy_tolerance_j_m2: f64,
    pub work_absolute_tolerance_j: f64,
    pub work_relative_tolerance: f64,
    /// Margin in the stated sufficient normal/friction feedback screen, (0,1).
    pub minimum_monotonicity_margin: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ReclosureError {
    Invalid(&'static str),
    Source(MaterialStatePointError),
    Material(CohesiveError),
    Friction(TriboError),
    NotFullyFailed,
    StateMismatch,
    Cancelled,
    Budget,
    Range(&'static str),
    CouplingUnestablished {
        margin: f64,
        required: f64,
    },
    NotConverged,
    ForceMismatch {
        component: usize,
        residual_n: f64,
        tolerance_n: f64,
    },
    WorkMismatch {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for ReclosureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fractured-interface reclosure refused: {self:?}")
    }
}
impl std::error::Error for ReclosureError {}
impl From<MaterialStatePointError> for ReclosureError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Source(value)
    }
}
impl From<CohesiveError> for ReclosureError {
    fn from(value: CohesiveError) -> Self {
        Self::Material(value)
    }
}
impl From<TriboError> for ReclosureError {
    fn from(value: TriboError) -> Self {
        Self::Friction(value)
    }
}

/// Accepted original material history and residual-contact occupancy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReclosureState<'a> {
    identity: ContentHash,
    point: CohesivePoint<'a>,
    occupancy: ReclosureRegime,
}
impl<'a> ReclosureState<'a> {
    #[must_use]
    pub const fn point(&self) -> CohesivePoint<'a> {
        self.point
    }
    #[must_use]
    pub const fn occupancy(&self) -> ReclosureRegime {
        self.occupancy
    }

    /// Keep a moving closed interface on its occupied kinetic branch. A zero
    /// net slip over a step need not mean zero endpoint velocity. This update
    /// changes only next-step occupancy, not the audited straight-path work or
    /// the original material history. No velocity threshold is introduced.
    pub fn with_endpoint_tangential_velocity(
        mut self,
        velocity_m_s: [f64; 2],
    ) -> Result<Self, ReclosureError> {
        for velocity in velocity_m_s {
            finite(velocity)?;
        }
        if self.point.jump_m()[0] < 0.0 && velocity_m_s.iter().any(|value| *value != 0.0) {
            self.occupancy = ReclosureRegime::Sliding;
        }
        Ok(self)
    }
}

/// Original constitutive reevaluation on actual receiving-body motion.
#[derive(Clone, Debug, PartialEq)]
pub struct ReclosureAudit<'a> {
    pub state: ReclosureState<'a>,
    pub material: CohesiveRateStep<'a>,
    pub regime: ReclosureRegime,
    pub expected_force_n: [f64; 3],
    pub force_residual_n: [f64; 3],
    pub force_tolerance_n: [f64; 3],
    /// Nonnegative residual sliding work, with no automatic heat recipient.
    pub friction_dissipation_j: f64,
    /// Applied-force path work minus original material work and friction loss.
    pub work_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReclosureStep<'a> {
    pub audit: ReclosureAudit<'a>,
    /// Resistance to positive normal opening and each positive tangent [N].
    pub mean_force_n: [f64; 3],
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    pub iterations: usize,
    pub evaluations: usize,
}

/// Small reusable solve borrowing both complete original source owners.
/// Construction and continuation allocate no scratch arrays on the heap.
pub struct ReclosureLoadLine<'a> {
    cohesive: &'a ResolvedBilinearCohesive,
    friction: &'a ResidualCoulombSource,
    area: f64,
    compliance: [[f64; 3]; 3],
    config: ReclosureConfig,
    identity: ContentHash,
    static_mu: f64,
    kinetic_mu: f64,
    tangent_scale: f64,
    monotonicity_margin: f64,
}
impl<'a> ReclosureLoadLine<'a> {
    pub fn new(
        cohesive: &'a ResolvedBilinearCohesive,
        friction: &'a ResidualCoulombSource,
        area_m2: f64,
        compliance_m_per_n: [[f64; 3]; 3],
        config: ReclosureConfig,
    ) -> Result<Self, ReclosureError> {
        validate_config(config)?;
        positive(area_m2)?;
        friction.require_cohesive_source(cohesive)?;
        let static_mu = friction.law().static_coefficient()?;
        let kinetic_mu = friction.law().zero_speed_kinetic_coefficient()?;
        let s = compliance_m_per_n;
        let scale = s.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
        positive(scale)?;
        for i in 0..3 {
            for j in 0..3 {
                finite(s[i][j])?;
                if (s[i][j] - s[j][i]).abs() > 64.0 * f64::EPSILON * scale {
                    return Err(ReclosureError::Invalid(
                        "reciprocal mechanical compliance required",
                    ));
                }
            }
        }
        positive(s[0][0])?;
        let tangent = [[s[1][1], s[1][2]], [s[2][1], s[2][2]]];
        let h = core::array::from_fn(|i| {
            core::array::from_fn(|j| tangent[i][j] - s[i + 1][0] * (s[0][j + 1] / s[0][0]))
        });
        let (t_min, t_max) = eigen_limits(tangent)?;
        let (h_min, _) = eigen_limits(h)?;
        let tangent_scale = positive((t_min * t_max).sqrt())?;
        // For sliding t=lambda*d, dt/dlambda=(I+lambda*H)^-1*d,
        // H=T-alpha*c*c^T, 0<=alpha<=1/Snn. This bound keeps
        // d(|t|-mu*p)/dlambda positive throughout the scalar branch.
        let feedback =
            finite(static_mu * (s[0][1].hypot(s[0][2]) / s[0][0]) * (t_max / h_min).sqrt())?;
        let margin = 1.0 - feedback;
        if margin < config.minimum_monotonicity_margin {
            return Err(ReclosureError::CouplingUnestablished {
                margin,
                required: config.minimum_monotonicity_margin,
            });
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.cohesive-reclosure.v1");
        hash.update(cohesive.resolved().identity().as_bytes());
        hash.update(friction.resolved().identity().as_bytes());
        hash.update(&area_m2.to_bits().to_le_bytes());
        for value in s.iter().flatten() {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for value in [
            config.jump_absolute_tolerance_m,
            config.jump_relative_tolerance,
            config.maximum_jump_increment_m,
            config.maximum_force_n,
            config.force_absolute_tolerance_n,
            config.force_relative_tolerance,
            config.material_energy_tolerance_j_m2,
            config.work_absolute_tolerance_j,
            config.work_relative_tolerance,
            config.minimum_monotonicity_margin,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&(config.max_iterations as u64).to_le_bytes());
        hash.update(&(config.max_evaluations as u64).to_le_bytes());
        Ok(Self {
            cohesive,
            friction,
            area: area_m2,
            compliance: s,
            config,
            identity: hash.finalize(),
            static_mu,
            kinetic_mu,
            tangent_scale,
            monotonicity_margin: margin,
        })
    }
    #[must_use]
    pub const fn cohesive_source(&self) -> &'a ResolvedBilinearCohesive {
        self.cohesive
    }
    #[must_use]
    pub const fn friction_source(&self) -> &'a ResidualCoulombSource {
        self.friction
    }
    #[must_use]
    pub const fn config(&self) -> ReclosureConfig {
        self.config
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area
    }
    #[must_use]
    pub const fn compliance_m_per_n(&self) -> [[f64; 3]; 3] {
        self.compliance
    }
    #[must_use]
    pub const fn monotonicity_margin(&self) -> f64 {
        self.monotonicity_margin
    }

    pub fn adopt(
        &self,
        point: CohesivePoint<'a>,
        receiving: [&MaterialStateId; 2],
        activation: ReclosureActivation,
    ) -> Result<ReclosureState<'a>, ReclosureError> {
        self.cohesive
            .require_surface_materials(receiving[0], receiving[1])?;
        self.friction.require_surface_materials(receiving)?;
        if !core::ptr::eq(point.law(), self.cohesive.law()) {
            return Err(ReclosureError::StateMismatch);
        }
        if !point.failed() {
            return Err(ReclosureError::NotFullyFailed);
        }
        let occupancy = if point.jump_m()[0] >= 0.0 {
            ReclosureRegime::Open
        } else if activation == ReclosureActivation::Sliding {
            ReclosureRegime::Sliding
        } else {
            ReclosureRegime::Sticking
        };
        Ok(ReclosureState {
            identity: self.identity,
            point,
            occupancy,
        })
    }

    pub fn trial(
        &self,
        previous: ReclosureState<'a>,
        free_jump_m: [f64; 3],
        dt_s: f64,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<ReclosureStep<'a>, ReclosureError> {
        self.require_state(previous)?;
        positive(dt_s)?;
        for x in free_jump_m {
            finite(x)?;
        }
        let old = previous.point.jump_m();
        let mut evaluations = 0;
        let free = self.candidate(old, free_jump_m, 0.0, &mut evaluations, &mut keep_going)?;
        if free.pressure == 0.0 {
            spend(&mut evaluations, self.config, &mut keep_going)?;
            let audit = self.audit(previous, free.jump, dt_s, free.force, ReclosureRegime::Open)?;
            poll(&mut keep_going)?;
            return Ok(ReclosureStep {
                audit,
                mean_force_n: free.force,
                compatibility_residual_m: [0.0; 3],
                compatibility_tolerance_m: self.jump_tolerance(old, free_jump_m, free.jump)?,
                iterations: 0,
                evaluations,
            });
        }
        let static_trial =
            self.candidate(old, free_jump_m, 1.0, &mut evaluations, &mut keep_going)?;
        let capture_mu = if previous.occupancy == ReclosureRegime::Sliding {
            self.kinetic_mu
        } else {
            self.static_mu
        };
        let static_norm = static_trial.force[1].hypot(static_trial.force[2]);
        let cap = finite(capture_mu * static_trial.pressure)?;
        let (chosen, regime, iterations) =
            if static_norm <= cap + self.force_tolerance(static_norm, cap)? {
                (static_trial, ReclosureRegime::Sticking, 0)
            } else {
                if self.kinetic_mu == 0.0 {
                    (free, ReclosureRegime::Sliding, 0)
                } else {
                    let mut lo = 0.0;
                    let mut hi = 1.0;
                    let mut found = None;
                    for iteration in 1..=self.config.max_iterations {
                        let z = f64::midpoint(lo, hi);
                        let trial =
                            self.candidate(old, free_jump_m, z, &mut evaluations, &mut keep_going)?;
                        let magnitude = trial.force[1].hypot(trial.force[2]);
                        let capacity = finite(self.kinetic_mu * trial.pressure)?;
                        let residual = magnitude - capacity;
                        if residual.abs() <= self.force_tolerance(magnitude, capacity)? {
                            found = Some((trial, ReclosureRegime::Sliding, iteration));
                            break;
                        }
                        if residual > 0.0 {
                            hi = z;
                        } else {
                            lo = z;
                        }
                        if lo.next_up() >= hi {
                            break;
                        }
                    }
                    found.ok_or(ReclosureError::NotConverged)?
                }
            };
        spend(&mut evaluations, self.config, &mut keep_going)?;
        let audit = self.audit(previous, chosen.jump, dt_s, chosen.force, regime)?;
        let mut residual = [0.0; 3];
        for i in 0..3 {
            residual[i] =
                finite(chosen.jump[i] - free_jump_m[i] + dot3(self.compliance[i], chosen.force))?;
        }
        let tolerance = self.jump_tolerance(old, free_jump_m, chosen.jump)?;
        if max_abs(residual) > tolerance {
            return Err(ReclosureError::NotConverged);
        }
        poll(&mut keep_going)?;
        Ok(ReclosureStep {
            audit,
            mean_force_n: chosen.force,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
            iterations,
            evaluations,
        })
    }

    /// Reevaluate the original material and residual friction on actual body
    /// motion. Sticking reaction is checked against its cone; sliding force is
    /// reconstructed independently from actual slip and original mean pressure.
    pub fn audit(
        &self,
        previous: ReclosureState<'a>,
        jump: [f64; 3],
        dt_s: f64,
        applied: [f64; 3],
        regime: ReclosureRegime,
    ) -> Result<ReclosureAudit<'a>, ReclosureError> {
        self.require_state(previous)?;
        positive(dt_s)?;
        let old = previous.point.jump_m();
        let delta = core::array::from_fn(|i| jump[i] - old[i]);
        for i in 0..3 {
            finite(jump[i])?;
            finite(applied[i])?;
            if delta[i].abs() > self.config.maximum_jump_increment_m {
                return Err(ReclosureError::Range("actual jump increment"));
            }
            if applied[i].abs() > self.config.maximum_force_n {
                return Err(ReclosureError::Range("actual applied force"));
            }
        }
        let (material, mean) = previous.point.trial_over_time_with_mean_traction(
            jump,
            dt_s,
            self.config.material_energy_tolerance_j_m2,
        )?;
        if !material.equilibrium.point.failed()
            || material.equilibrium.energy.dissipated != 0.0
            || material.viscous_dissipated_j_m2 != 0.0
            || mean[1] != 0.0
            || mean[2] != 0.0
        {
            return Err(ReclosureError::Invalid(
                "failed bond must retain only original compression storage",
            ));
        }
        let pressure = finite(-self.area * mean[0])?;
        if pressure < 0.0 {
            return Err(ReclosureError::Invalid(
                "nonnegative compression pressure required",
            ));
        }
        let slip = delta[1].hypot(delta[2]);
        let mut expected = [-pressure, 0.0, 0.0];
        let loss = match regime {
            ReclosureRegime::Open => {
                if pressure != 0.0 {
                    return Err(ReclosureError::Invalid("open path carries compression"));
                }
                0.0
            }
            ReclosureRegime::Sticking => {
                let tolerance = self.jump_tolerance(old, jump, jump)?;
                if slip > tolerance {
                    return Err(ReclosureError::Range("actual sticking slip"));
                }
                let mu = if previous.occupancy == ReclosureRegime::Sliding {
                    self.kinetic_mu
                } else {
                    self.static_mu
                };
                let magnitude = applied[1].hypot(applied[2]);
                let cap = finite(mu * pressure)?;
                if magnitude > cap + self.force_tolerance(magnitude, cap)? {
                    return Err(ReclosureError::Range(
                        "actual sticking reaction exceeds capture cone",
                    ));
                }
                expected[1] = applied[1];
                expected[2] = applied[2];
                0.0
            }
            ReclosureRegime::Sliding => {
                if slip == 0.0 {
                    return Err(ReclosureError::Invalid(
                        "sliding needs actual nonzero displacement",
                    ));
                }
                let cap = finite(self.kinetic_mu * pressure)?;
                for i in 1..3 {
                    expected[i] = finite(cap * (delta[i] / slip))?;
                }
                finite(cap * slip)?
            }
        };
        let mut force_residual = [0.0; 3];
        let mut force_tolerance = [0.0; 3];
        for i in 0..3 {
            force_residual[i] = finite(applied[i] - expected[i])?;
            force_tolerance[i] = self.force_tolerance(applied[i], expected[i])?;
            if force_residual[i].abs() > force_tolerance[i] {
                return Err(ReclosureError::ForceMismatch {
                    component: i,
                    residual_n: force_residual[i],
                    tolerance_n: force_tolerance[i],
                });
            }
        }
        let work = finite(dot3(applied, delta))?;
        let material_work = finite(self.area * material.total_work_j_m2())?;
        let defect = finite(work - material_work - loss)?;
        let tolerance = finite(
            self.config.work_absolute_tolerance_j
                + self.config.work_relative_tolerance
                    * work.abs().max(material_work.abs()).max(loss),
        )?;
        if defect.abs() > tolerance {
            return Err(ReclosureError::WorkMismatch {
                residual_j: defect,
                tolerance_j: tolerance,
            });
        }
        let occupancy = if jump[0] >= 0.0 {
            ReclosureRegime::Open
        } else {
            regime
        };
        Ok(ReclosureAudit {
            state: ReclosureState {
                identity: self.identity,
                point: material.equilibrium.point,
                occupancy,
            },
            material,
            regime,
            expected_force_n: expected,
            force_residual_n: force_residual,
            force_tolerance_n: force_tolerance,
            friction_dissipation_j: loss,
            work_residual_j: defect,
        })
    }

    fn require_state(&self, state: ReclosureState<'a>) -> Result<(), ReclosureError> {
        if state.identity != self.identity || !core::ptr::eq(state.point.law(), self.cohesive.law())
        {
            return Err(ReclosureError::StateMismatch);
        }
        if !state.point.failed() {
            return Err(ReclosureError::NotFullyFailed);
        }
        Ok(())
    }
    fn force_tolerance(&self, a: f64, b: f64) -> Result<f64, ReclosureError> {
        finite(
            self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance * a.abs().max(b.abs()),
        )
    }
    fn jump_tolerance(
        &self,
        old: [f64; 3],
        free: [f64; 3],
        jump: [f64; 3],
    ) -> Result<f64, ReclosureError> {
        let scale = (0..3)
            .map(|i| (free[i] - old[i]).abs().max((jump[i] - old[i]).abs()))
            .fold(0.0, f64::max);
        finite(self.config.jump_absolute_tolerance_m + self.config.jump_relative_tolerance * scale)
    }
    fn candidate(
        &self,
        old: [f64; 3],
        free: [f64; 3],
        z: f64,
        evaluations: &mut usize,
        keep_going: &mut impl FnMut() -> bool,
    ) -> Result<Candidate, ReclosureError> {
        let s = self.compliance;
        let v = [free[1] - old[1], free[2] - old[2]];
        let row = [s[0][1], s[0][2]];
        let col = [s[1][0], s[2][0]];
        let mut matrix = [[s[1][1] * z, s[1][2] * z], [s[2][1] * z, s[2][2] * z]];
        for i in 0..2 {
            matrix[i][i] += (1.0 - z) * self.tangent_scale;
        }
        let inv = inverse2(matrix)?;
        let w = mul2(inv, v).map(|x| z * x);
        let c = mul2(inv, col).map(|x| z * x);
        let free_normal = finite(free[0] - dot2(row, w))?;
        let compliance = positive(s[0][0] - dot2(row, c))?;
        let (normal, pressure) =
            self.normal(old[0], free_normal, compliance, evaluations, keep_going)?;
        let tangent = [
            finite(w[0] + c[0] * pressure)?,
            finite(w[1] + c[1] * pressure)?,
        ];
        let force = [-pressure, tangent[0], tangent[1]];
        let mut jump = [normal, 0.0, 0.0];
        for i in 1..3 {
            jump[i] = finite(free[i] - dot3(s[i], force))?;
        }
        Ok(Candidate {
            jump,
            force,
            pressure,
        })
    }
    fn normal(
        &self,
        old: f64,
        free: f64,
        compliance: f64,
        evaluations: &mut usize,
        keep_going: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ReclosureError> {
        let stiffness = positive(self.area * self.cohesive.law().compression_stiffness_pa_m())?;
        let at = |g: f64| -> Result<(f64, f64), ReclosureError> {
            let (mean, derivative) = compression_mean(old, g)?;
            Ok((finite(-stiffness * mean)?, finite(-stiffness * derivative)?))
        };
        spend(evaluations, self.config, keep_going)?;
        let (p0, _) = at(free)?;
        if p0 == 0.0 {
            return Ok((free, 0.0));
        }
        let mut lo = free;
        let mut hi = finite(free + compliance * p0)?;
        let mut gap = f64::midpoint(lo, hi);
        for _ in 0..self.config.max_iterations {
            spend(evaluations, self.config, keep_going)?;
            let (pressure, derivative) = at(gap)?;
            let residual = finite((gap - free) - compliance * pressure)?;
            let tolerance = self.config.jump_absolute_tolerance_m * 0.125
                + self.config.jump_relative_tolerance
                    * 0.125
                    * (gap - old).abs().max((free - old).abs());
            if residual.abs() <= tolerance {
                return Ok((gap, pressure));
            }
            if residual > 0.0 {
                hi = gap;
            } else {
                lo = gap;
            }
            let newton = gap - residual / (1.0 - compliance * derivative);
            let next = if newton > lo && newton < hi {
                newton
            } else {
                f64::midpoint(lo, hi)
            };
            if next == gap {
                break;
            }
            gap = next;
        }
        Err(ReclosureError::NotConverged)
    }
}

#[derive(Clone, Copy)]
struct Candidate {
    jump: [f64; 3],
    force: [f64; 3],
    pressure: f64,
}

// Mean of min(g,0) and its exact endpoint derivative on a straight path.
fn compression_mean(a: f64, b: f64) -> Result<(f64, f64), ReclosureError> {
    if a >= 0.0 && b >= 0.0 {
        return Ok((0.0, 0.0));
    }
    if a <= 0.0 && b <= 0.0 {
        return Ok((finite(0.5 * a + 0.5 * b)?, 0.5));
    }
    if a < 0.0 {
        let ratio = -a / (b - a);
        Ok((finite(0.5 * a * ratio)?, finite(0.5 * ratio * ratio)?))
    } else {
        let ratio = -b / (a - b);
        Ok((
            finite(0.5 * b * ratio)?,
            finite(ratio * (1.0 - 0.5 * ratio))?,
        ))
    }
}
fn validate_config(c: ReclosureConfig) -> Result<(), ReclosureError> {
    if !(1..=256).contains(&c.max_iterations) || c.max_evaluations < 2 {
        return Err(ReclosureError::Invalid(
            "bounded positive iteration/evaluation budgets required",
        ));
    }
    for x in [
        c.jump_absolute_tolerance_m,
        c.maximum_jump_increment_m,
        c.maximum_force_n,
        c.force_absolute_tolerance_n,
        c.material_energy_tolerance_j_m2,
        c.work_absolute_tolerance_j,
    ] {
        positive(x)?;
    }
    for x in [
        c.jump_relative_tolerance,
        c.force_relative_tolerance,
        c.work_relative_tolerance,
        c.minimum_monotonicity_margin,
    ] {
        if !x.is_finite() || x <= 0.0 || x >= 1.0 {
            return Err(ReclosureError::Invalid(
                "relative tolerances and monotonicity margin must be in (0,1)",
            ));
        }
    }
    Ok(())
}
fn eigen_limits(m: [[f64; 2]; 2]) -> Result<(f64, f64), ReclosureError> {
    positive(m[0][0])?;
    positive(m[1][1])?;
    let scale = m.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
    let a = m[0][0] / scale;
    let b = m[0][1] / scale;
    let c = m[1][0] / scale;
    let d = m[1][1] / scale;
    let determinant = positive(a * d - b * c)?;
    let discriminant = finite((a - d) * (a - d) + 4.0 * b * c)?;
    if discriminant < 0.0 {
        return Err(ReclosureError::Invalid(
            "positive reciprocal mechanical block required",
        ));
    }
    let maximum = positive(0.5 * (a + d + discriminant.sqrt()))?;
    Ok((
        positive(scale * (determinant / maximum))?,
        positive(scale * maximum)?,
    ))
}
fn inverse2(m: [[f64; 2]; 2]) -> Result<[[f64; 2]; 2], ReclosureError> {
    let scale = m.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
    positive(scale)?;
    let a = m[0][0] / scale;
    let b = m[0][1] / scale;
    let c = m[1][0] / scale;
    let d = m[1][1] / scale;
    let det = positive(a * d - b * c)?;
    let factor = positive((1.0 / scale) / det)?;
    let result = [[d * factor, -b * factor], [-c * factor, a * factor]];
    for value in result.iter().flatten() {
        finite(*value)?;
    }
    Ok(result)
}
fn mul2(m: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [dot2(m[0], x), dot2(m[1], x)]
}
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn max_abs(v: [f64; 3]) -> f64 {
    v.into_iter().map(f64::abs).fold(0.0, f64::max)
}
fn finite(x: f64) -> Result<f64, ReclosureError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(ReclosureError::Invalid(
            "nonfinite reclosure input/arithmetic",
        ))
    }
}
fn positive(x: f64) -> Result<f64, ReclosureError> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(ReclosureError::Invalid(
            "positive finite reclosure input required",
        ))
    }
}
fn poll(keep_going: &mut impl FnMut() -> bool) -> Result<(), ReclosureError> {
    if keep_going() {
        Ok(())
    } else {
        Err(ReclosureError::Cancelled)
    }
}
fn spend(
    n: &mut usize,
    config: ReclosureConfig,
    keep_going: &mut impl FnMut() -> bool,
) -> Result<(), ReclosureError> {
    poll(keep_going)?;
    if *n >= config.max_evaluations {
        return Err(ReclosureError::Budget);
    }
    *n += 1;
    Ok(())
}

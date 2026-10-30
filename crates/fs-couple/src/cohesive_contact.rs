//! Three-component cohesive work on a caller's actual held-force response.
//!
//! The material owner supplies the analytic mean traction along the complete
//! opening/shear path. This owner solves `g = g_free - S A mean_traction(g)`
//! with all entries of the original mechanical compliance S. Endpoint traction
//! is never substituted for the force whose work drives the bodies. The local
//! damped Newton solve retains the previous physical jump as its initial guess;
//! it does not certify uniqueness through softening or resolve intrastep path
//! reversals. Refine the physical timestep for trajectory convergence.

use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_material::cohesive::{CohesiveError, CohesivePoint, CohesiveRateStep};

/// Explicit physical range and bounded nonlinear work for one cohesive patch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesiveContactConfig {
    /// Maximum accepted Newton corrections, at most 128.
    pub max_iterations: usize,
    /// Complete material-path evaluations, including Jacobians and line search.
    pub max_material_evaluations: usize,
    /// Maximum backtracking trials for one correction, at most 32.
    pub max_backtracks: usize,
    /// Absolute componentwise compatibility tolerance [m].
    pub jump_absolute_tolerance_m: f64,
    /// Relative compatibility tolerance, scaled by this step's displacement.
    pub jump_relative_tolerance: f64,
    /// Maximum change of any local jump component in one step [m].
    pub maximum_jump_increment_m: f64,
    /// Maximum absolute component of the held physical patch force [N].
    pub maximum_force_n: f64,
    /// Absolute material path-work closure budget per original bonded area [J/m2].
    pub material_energy_tolerance_j_m2: f64,
}

/// A refused trial never changes the original cohesive point or caller mechanics.
#[derive(Debug, Clone, PartialEq)]
pub enum CohesiveContactError {
    Invalid(&'static str),
    Material(CohesiveError),
    Factor(LuWorkspaceError),
    Cancelled,
    Budget(&'static str),
    Range(&'static str),
    NotConverged {
        residual_m: f64,
        tolerance_m: f64,
        iterations: usize,
    },
}
impl core::fmt::Display for CohesiveContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive contact trial refused: {self:?}")
    }
}
impl std::error::Error for CohesiveContactError {}
impl From<CohesiveError> for CohesiveContactError {
    fn from(value: CohesiveError) -> Self {
        Self::Material(value)
    }
}
impl From<LuWorkspaceError> for CohesiveContactError {
    fn from(value: LuWorkspaceError) -> Self {
        Self::Factor(value)
    }
}

/// Original constitutive successor and its actual work-conjugate patch force.
#[derive(Clone, Debug, PartialEq)]
pub struct CohesiveContactStep<'a> {
    pub material: CohesiveRateStep<'a>,
    /// Force resisting positive opening and each positive shear coordinate [N].
    pub mean_force_n: [f64; 3],
    /// `g - g_free + S mean_force`, on the solved path [m].
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    pub iterations: usize,
    pub material_evaluations: usize,
}

/// Reusable 3-by-3 interface solve; allocation occurs only at construction.
///
/// The caller owns the complete mechanical map and its work audit. Compliance
/// is neither diagonalized nor symmetrized; anisotropic and cross-axis body
/// response is retained. Area is the original bonded reference area, not an
/// invented fraction of normal-load contact. Damage remains with the material.
#[derive(Debug)]
pub struct CohesiveLoadLine {
    area_m2: f64,
    compliance_m_n: [[f64; 3]; 3],
    config: CohesiveContactConfig,
    solve: LuWorkspace,
}

impl CohesiveLoadLine {
    pub fn new(
        area_m2: f64,
        compliance_m_n: [[f64; 3]; 3],
        config: CohesiveContactConfig,
    ) -> Result<Self, CohesiveContactError> {
        if !area_m2.is_finite()
            || area_m2 <= 0.0
            || compliance_m_n.iter().flatten().any(|v| !v.is_finite())
            || (0..3).any(|i| compliance_m_n[i][i] < 0.0)
            || !valid_config(config)
        {
            return Err(CohesiveContactError::Invalid(
                "finite physical area/compliance and explicit positive solve budgets required",
            ));
        }
        Ok(Self {
            area_m2,
            compliance_m_n,
            config,
            solve: LuWorkspace::new(3)?,
        })
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    #[must_use]
    pub const fn compliance_m_per_n(&self) -> [[f64; 3]; 3] {
        self.compliance_m_n
    }
    #[must_use]
    pub const fn config(&self) -> CohesiveContactConfig {
        self.config
    }

    /// Solve an immutable trial over one declared straight jump path. `free_jump`
    /// is the absolute endpoint predicted without this patch, not an increment.
    /// `keep_going` must return false when the caller's cancellation is requested;
    /// it is polled before every material evaluation and before returning success.
    /// Only scratch changes on refusal. The returned point is committed by the
    /// caller after its actual mechanics and complete energy window are checked.
    pub fn trial<'a>(
        &mut self,
        point: CohesivePoint<'a>,
        free_jump: [f64; 3],
        dt_s: f64,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<CohesiveContactStep<'a>, CohesiveContactError> {
        if free_jump.iter().any(|v| !v.is_finite()) || !positive(dt_s) {
            return Err(CohesiveContactError::Invalid(
                "finite free jump and positive timestep required",
            ));
        }
        let old = point.jump_m();
        let scales = point.law().mixed_mode_parameters().map_or(
            [point.law().onset_separation_m(); 3],
            |p| {
                [
                    p.normal_strength_pa / p.normal_stiffness_pa_m,
                    p.tangential_strength_pa / p.tangential_stiffness_pa_m,
                    p.tangential_strength_pa / p.tangential_stiffness_pa_m,
                ]
            },
        );
        let mut evaluations = 0;
        let mut current = self.evaluate(
            point,
            old,
            free_jump,
            dt_s,
            &mut evaluations,
            &mut keep_going,
        )?;
        for iteration in 0..=self.config.max_iterations {
            let residual = norm_inf(current.compatibility_residual_m);
            let tolerance = current.compatibility_tolerance_m;
            if residual <= tolerance {
                if !keep_going() {
                    return Err(CohesiveContactError::Cancelled);
                }
                current.iterations = iteration;
                current.material_evaluations = evaluations;
                return Ok(current);
            }
            if iteration == self.config.max_iterations {
                return Err(CohesiveContactError::NotConverged {
                    residual_m: residual,
                    tolerance_m: tolerance,
                    iterations: iteration,
                });
            }
            let x = current.material.equilibrium.point.jump_m();
            let mut jacobian = [0.0; 9];
            for j in 0..3 {
                // Central differences on the analytic mean, scaled to actual
                // material onset length. Endpoint tangents are different and
                // cannot be reused as though they differentiated path work.
                let h = (1.0 / 1_048_576.0)
                    * scales[j]
                        .max(old[j].abs())
                        .max(x[j].abs())
                        .max(self.config.jump_absolute_tolerance_m);
                let mut xp = x;
                xp[j] += h;
                let mut xm = x;
                xm[j] -= h;
                let plus = self.evaluate(
                    point,
                    xp,
                    free_jump,
                    dt_s,
                    &mut evaluations,
                    &mut keep_going,
                );
                let minus = self.evaluate(
                    point,
                    xm,
                    free_jump,
                    dt_s,
                    &mut evaluations,
                    &mut keep_going,
                );
                let (p, m, width) = match (plus, minus) {
                    (Ok(p), Ok(m)) => (
                        p.compatibility_residual_m,
                        m.compatibility_residual_m,
                        xp[j] - xm[j],
                    ),
                    (Ok(p), Err(e)) if outside_trial_domain(&e) => (
                        p.compatibility_residual_m,
                        current.compatibility_residual_m,
                        xp[j] - x[j],
                    ),
                    (Err(e), Ok(m)) if outside_trial_domain(&e) => (
                        current.compatibility_residual_m,
                        m.compatibility_residual_m,
                        x[j] - xm[j],
                    ),
                    (Err(e), _) | (_, Err(e)) => return Err(e),
                };
                if !(width > 0.0 && width.is_finite()) {
                    return Err(CohesiveContactError::Invalid(
                        "unrepresentable local derivative displacement",
                    ));
                }
                for i in 0..3 {
                    jacobian[i * 3 + j] = checked((p[i] - m[i]) / width)?;
                }
            }
            let rhs = current.compatibility_residual_m.map(|r| -r);
            let mut correction = [0.0; 3];
            self.solve.solve_into(&jacobian, &rhs, &mut correction)?;
            let mut scale = 1.0;
            let mut accepted = None;
            for _ in 0..self.config.max_backtracks {
                let candidate = core::array::from_fn(|i| x[i] + scale * correction[i]);
                match self.evaluate(
                    point,
                    candidate,
                    free_jump,
                    dt_s,
                    &mut evaluations,
                    &mut keep_going,
                ) {
                    Ok(next)
                        if norm_inf(next.compatibility_residual_m) < residual
                            || norm_inf(next.compatibility_residual_m)
                                <= next.compatibility_tolerance_m =>
                    {
                        accepted = Some(next);
                        break;
                    }
                    Ok(_) => {}
                    Err(e) if outside_trial_domain(&e) => {}
                    Err(e) => return Err(e),
                }
                scale *= 0.5;
            }
            current = accepted.ok_or(CohesiveContactError::NotConverged {
                residual_m: residual,
                tolerance_m: tolerance,
                iterations: iteration + 1,
            })?;
        }
        unreachable!("bounded cohesive solve returns at its final iteration")
    }

    fn evaluate<'a>(
        &self,
        point: CohesivePoint<'a>,
        jump: [f64; 3],
        free: [f64; 3],
        dt: f64,
        count: &mut usize,
        keep_going: &mut impl FnMut() -> bool,
    ) -> Result<CohesiveContactStep<'a>, CohesiveContactError> {
        if !keep_going() {
            return Err(CohesiveContactError::Cancelled);
        }
        if *count >= self.config.max_material_evaluations {
            return Err(CohesiveContactError::Budget("material path evaluations"));
        }
        *count += 1;
        let old = point.jump_m();
        let increment = core::array::from_fn(|i| jump[i] - old[i]);
        if jump.iter().chain(increment.iter()).any(|v| !v.is_finite())
            || norm_inf(increment) > self.config.maximum_jump_increment_m
        {
            return Err(CohesiveContactError::Range(
                "per-step cohesive jump increment",
            ));
        }
        let (material, mean) = point.trial_over_time_with_mean_traction(
            jump,
            dt,
            self.config.material_energy_tolerance_j_m2,
        )?;
        let force = mean.map(|v| self.area_m2 * v);
        if force.iter().any(|v| !v.is_finite()) || norm_inf(force) > self.config.maximum_force_n {
            return Err(CohesiveContactError::Range("cohesive held force"));
        }
        let response: [f64; 3] =
            core::array::from_fn(|i| (0..3).map(|j| self.compliance_m_n[i][j] * force[j]).sum());
        let residual: [f64; 3] = core::array::from_fn(|i| jump[i] - free[i] + response[i]);
        if residual.iter().any(|v| !v.is_finite()) {
            return Err(CohesiveContactError::Invalid(
                "cohesive compatibility overflow",
            ));
        }
        let free_increment = core::array::from_fn(|i| free[i] - old[i]);
        let scale = norm_inf(increment)
            .max(norm_inf(free_increment))
            .max(norm_inf(response));
        let tolerance = checked(
            self.config.jump_absolute_tolerance_m + self.config.jump_relative_tolerance * scale,
        )?;
        Ok(CohesiveContactStep {
            material,
            mean_force_n: force,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
            iterations: 0,
            material_evaluations: 0,
        })
    }
}

fn outside_trial_domain(error: &CohesiveContactError) -> bool {
    matches!(
        error,
        CohesiveContactError::Range(_)
            | CohesiveContactError::Material(
                CohesiveError::Invalid(_) | CohesiveError::SeparationLimit
            )
    )
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn valid_config(config: CohesiveContactConfig) -> bool {
    (1..=128).contains(&config.max_iterations)
        && (1..=32).contains(&config.max_backtracks)
        && config.max_material_evaluations > 0
        && positive(config.jump_absolute_tolerance_m)
        && positive(config.jump_relative_tolerance)
        && config.jump_relative_tolerance < 1.0
        && positive(config.maximum_jump_increment_m)
        && positive(config.maximum_force_n)
        && positive(config.material_energy_tolerance_j_m2)
}
fn checked(value: f64) -> Result<f64, CohesiveContactError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CohesiveContactError::Invalid(
            "derived contact value is not finite",
        ))
    }
}
fn norm_inf(v: [f64; 3]) -> f64 {
    v.into_iter().map(f64::abs).fold(0.0, f64::max)
}

#[cfg(test)]
mod tests;

/// Finite bonded patches with complete cross-site mechanical compliance.
pub mod patch;

/// Objective cohesive transport and exact finite-rotation work ports.
pub mod objective;

/// Original finite-rotation rigid bodies coupled through a cohesive patch.
pub mod rigid;

/// Residual frictional reclosure retaining completely failed bond history.
pub mod reclosure;

/// Damaged-area rigid friction during progressive cohesive degradation.
pub mod frictional_damage;

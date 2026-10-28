//! A finite bonded patch with coupled opening/shear history at every site.
//!
//! The caller supplies the original quadrature areas and complete mechanical
//! cross-compliance. This owner solves all `3N` jumps simultaneously:
//! `g_i = g_free_i - sum_j S_ij A_j mean_traction_j(g_j)`.
//! Finite lever arms therefore transmit moments and redistribute force as an
//! edge debonds. No independent-site solve or averaged damage is substituted.
//!
//! Each Jacobian column reevaluates only its local material point, then applies
//! every entry of the corresponding compliance columns. Material evaluations
//! scale linearly with site count per Jacobian; dense algebra remains quadratic
//! for assembly and cubic for the bounded prepared LU solve. All storage is
//! allocated at construction. Accepted points remain the caller's immutable
//! history until its actual mechanics/work transaction succeeds.

use super::{CohesiveContactConfig, CohesiveContactError, norm_inf, positive, valid_config};
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_material::cohesive::CohesivePoint;

/// Explicit dense-solve ceiling: at most 192 simultaneous jump coordinates.
pub const MAX_COHESIVE_PATCH_POINTS: usize = 64;

/// Spatial capacity and shared per-site physical/tolerance limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesivePatchConfig {
    /// Caller capacity in 1..=MAX_COHESIVE_PATCH_POINTS. No sites are truncated.
    pub max_points: usize,
    /// Force, increment and compatibility gates apply separately at each site.
    /// The material-evaluation budget counts TOTAL point-law calls, including
    /// finite differences and rejected line-search trials, across the patch.
    pub contact: CohesiveContactConfig,
}

/// A refusal changes scratch only, never a supplied accepted point.
#[derive(Debug, Clone, PartialEq)]
pub enum CohesivePatchError {
    Invalid(&'static str),
    Capacity,
    Point {
        point: usize,
        source: CohesiveContactError,
    },
    Factor(LuWorkspaceError),
    Cancelled,
    Budget(&'static str),
    NotConverged {
        point: usize,
        residual_m: f64,
        tolerance_m: f64,
        iterations: usize,
    },
}
impl core::fmt::Display for CohesivePatchError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive patch trial refused: {self:?}")
    }
}
impl std::error::Error for CohesivePatchError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Point { source, .. } => Some(source),
            Self::Factor(source) => Some(source),
            _ => None,
        }
    }
}
impl From<LuWorkspaceError> for CohesivePatchError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Factor(error)
    }
}

/// Borrowed complete solution, valid until the solver is borrowed mutably again.
/// The native owner must reevaluate its original points on actual staged body
/// motion before publishing damage. No material/source histories are cloned or
/// committed by this view.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesivePatchTrial<'s> {
    pub jump_m: &'s [[f64; 3]],
    /// Original site area times the analytic mean traction [N].
    pub mean_force_n: &'s [[f64; 3]],
    /// `g - g_free + S F` at every site [m].
    pub compatibility_residual_m: &'s [[f64; 3]],
    /// A separately scaled absolute-plus-relative tolerance at each site [m].
    pub compatibility_tolerance_m: &'s [f64],
    pub iterations: usize,
    pub material_evaluations: usize,
}

#[derive(Debug)]
struct Evaluation {
    jump: Vec<[f64; 3]>,
    force: Vec<[f64; 3]>,
    residual: Vec<[f64; 3]>,
    tolerance: Vec<f64>,
}
impl Evaluation {
    fn new(points: usize) -> Result<Self, CohesivePatchError> {
        Ok(Self {
            jump: filled(points, [0.0; 3])?,
            force: filled(points, [0.0; 3])?,
            residual: filled(points, [0.0; 3])?,
            tolerance: filled(points, 0.0)?,
        })
    }

    fn worst(&self) -> Result<Worst, CohesivePatchError> {
        let mut worst = Worst {
            point: 0,
            residual: 0.0,
            tolerance: self.tolerance[0],
            merit: 0.0,
        };
        for point in 0..self.jump.len() {
            let residual = norm_inf(self.residual[point]);
            let merit = checked(residual / self.tolerance[point])?;
            if merit > worst.merit {
                worst = Worst {
                    point,
                    residual,
                    tolerance: self.tolerance[point],
                    merit,
                };
            }
        }
        Ok(worst)
    }
}

#[derive(Clone, Copy)]
struct Worst {
    point: usize,
    residual: f64,
    tolerance: f64,
    merit: f64,
}
impl Worst {
    fn refusal(self, iterations: usize) -> CohesivePatchError {
        CohesivePatchError::NotConverged {
            point: self.point,
            residual_m: self.residual,
            tolerance_m: self.tolerance,
            iterations,
        }
    }
}

/// Prepared dense coupling of heterogeneous original cohesive points.
/// The compliance uses row-major `[point0 normal,t1,t2,point1 normal,t1,t2,...]`
/// order on BOTH axes. It is never diagonalized, symmetrized or regularized.
/// Rank-deficient mechanical compliance is allowed; the solved Jacobian retains
/// its identity term. Singular softening branches are explicit refusals.
#[derive(Debug)]
pub struct CohesivePatchLoadLine {
    areas: Vec<f64>,
    compliance: Vec<f64>,
    config: CohesivePatchConfig,
    current: Evaluation,
    candidate: Evaluation,
    jacobian: Vec<f64>,
    rhs: Vec<f64>,
    correction: Vec<f64>,
    solve: LuWorkspace,
}

impl CohesivePatchLoadLine {
    pub fn new(
        areas_m2: Vec<f64>,
        compliance_m_per_n: Vec<f64>,
        config: CohesivePatchConfig,
    ) -> Result<Self, CohesivePatchError> {
        let count = areas_m2.len();
        if !(1..=MAX_COHESIVE_PATCH_POINTS).contains(&config.max_points)
            || count == 0
            || count > config.max_points
            || areas_m2.iter().any(|&area| !positive(area))
            || !valid_config(config.contact)
            || config.contact.max_material_evaluations < count
        {
            return Err(CohesivePatchError::Invalid(
                "bounded nonempty sites, positive original areas and complete solve budgets required",
            ));
        }
        // count <= 64 above establishes representable fixed dense extents.
        let dimension = 3 * count;
        let extent = dimension * dimension;
        if compliance_m_per_n.len() != extent
            || compliance_m_per_n.iter().any(|value| !value.is_finite())
            || (0..dimension).any(|i| compliance_m_per_n[i * dimension + i] < 0.0)
        {
            return Err(CohesivePatchError::Invalid(
                "complete finite 3N-by-3N mechanical compliance with nonnegative diagonal required",
            ));
        }
        Ok(Self {
            areas: areas_m2,
            compliance: compliance_m_per_n,
            config,
            current: Evaluation::new(count)?,
            candidate: Evaluation::new(count)?,
            jacobian: filled(extent, 0.0)?,
            rhs: filled(dimension, 0.0)?,
            correction: filled(dimension, 0.0)?,
            solve: LuWorkspace::new(dimension)?,
        })
    }

    #[must_use]
    pub fn point_count(&self) -> usize {
        self.areas.len()
    }
    #[must_use]
    pub fn areas_m2(&self) -> &[f64] {
        &self.areas
    }
    #[must_use]
    pub fn compliance_m_per_n(&self) -> &[f64] {
        &self.compliance
    }
    #[must_use]
    pub const fn config(&self) -> CohesivePatchConfig {
        self.config
    }

    /// Solve every site together from the immutable accepted jump/history.
    /// `free_jump` contains absolute endpoint jumps before ANY patch reaction.
    /// No per-trial allocation occurs. Polls precede every material evaluation,
    /// the bounded dense solve, every backtrack, and successful return. Softening
    /// can admit multiple roots; this local continuation is not a uniqueness or
    /// timestep-convergence certificate. Split physical reversals explicitly.
    pub fn trial(
        &mut self,
        points: &[CohesivePoint<'_>],
        free_jump: &[[f64; 3]],
        dt_s: f64,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<CohesivePatchTrial<'_>, CohesivePatchError> {
        poll(&mut keep_going)?;
        if points.len() != self.areas.len()
            || free_jump.len() != self.areas.len()
            || free_jump.iter().flatten().any(|value| !value.is_finite())
            || !positive(dt_s)
        {
            return Err(CohesivePatchError::Invalid(
                "original points and finite free jumps must match every site; positive time required",
            ));
        }
        let context = Context {
            areas: &self.areas,
            compliance: &self.compliance,
            config: self.config.contact,
            points,
            free: free_jump,
            dt: dt_s,
        };
        for (jump, point) in self.current.jump.iter_mut().zip(points) {
            *jump = point.jump_m();
        }
        let mut evaluations = 0;
        context.evaluate(&mut self.current, &mut evaluations, &mut keep_going)?;
        let dimension = 3 * points.len();
        for iteration in 0..=self.config.contact.max_iterations {
            let worst = self.current.worst()?;
            if worst.merit <= 1.0 {
                poll(&mut keep_going)?;
                return Ok(CohesivePatchTrial {
                    jump_m: &self.current.jump,
                    mean_force_n: &self.current.force,
                    compatibility_residual_m: &self.current.residual,
                    compatibility_tolerance_m: &self.current.tolerance,
                    iterations: iteration,
                    material_evaluations: evaluations,
                });
            }
            if iteration == self.config.contact.max_iterations {
                return Err(worst.refusal(iteration));
            }
            for (point, old_point) in points.iter().enumerate() {
                let old = old_point.jump_m();
                let x = self.current.jump[point];
                let scales = old_point.law().mixed_mode_parameters().map_or(
                    [old_point.law().onset_separation_m(); 3],
                    |p| {
                        [
                            p.normal_strength_pa / p.normal_stiffness_pa_m,
                            p.tangential_strength_pa / p.tangential_stiffness_pa_m,
                            p.tangential_strength_pa / p.tangential_stiffness_pa_m,
                        ]
                    },
                );
                for j in 0..3 {
                    let h = (1.0 / 1_048_576.0)
                        * scales[j]
                            .max(old[j].abs())
                            .max(x[j].abs())
                            .max(context.config.jump_absolute_tolerance_m);
                    let mut xp = x;
                    xp[j] += h;
                    let mut xm = x;
                    xm[j] -= h;
                    let plus = context.force(point, xp, &mut evaluations, &mut keep_going);
                    let minus = context.force(point, xm, &mut evaluations, &mut keep_going);
                    let (a, b, width) = match (plus, minus) {
                        (Ok(a), Ok(b)) => (a, b, xp[j] - xm[j]),
                        (Ok(a), Err(e)) if outside_trial_domain(&e) => {
                            (a, self.current.force[point], xp[j] - x[j])
                        }
                        (Err(e), Ok(b)) if outside_trial_domain(&e) => {
                            (self.current.force[point], b, x[j] - xm[j])
                        }
                        (Err(e), _) | (_, Err(e)) => return Err(e),
                    };
                    if !positive(width) {
                        return Err(CohesivePatchError::Invalid(
                            "unrepresentable local derivative displacement",
                        ));
                    }
                    let mut derivative = [0.0; 3];
                    for k in 0..3 {
                        derivative[k] = checked((a[k] - b[k]) / width)?;
                    }
                    let column = 3 * point + j;
                    for row in 0..dimension {
                        let response = (0..3)
                            .map(|k| {
                                context.compliance[row * dimension + 3 * point + k] * derivative[k]
                            })
                            .sum::<f64>();
                        self.jacobian[row * dimension + column] =
                            checked(response + if row == column { 1.0 } else { 0.0 })?;
                    }
                }
            }
            for point in 0..points.len() {
                for j in 0..3 {
                    self.rhs[3 * point + j] = -self.current.residual[point][j];
                }
            }
            poll(&mut keep_going)?;
            self.solve
                .solve_into(&self.jacobian, &self.rhs, &mut self.correction)?;
            poll(&mut keep_going)?;
            let mut scale = 1.0;
            let mut accepted = false;
            for _ in 0..self.config.contact.max_backtracks {
                poll(&mut keep_going)?;
                for point in 0..points.len() {
                    for j in 0..3 {
                        self.candidate.jump[point][j] =
                            self.current.jump[point][j] + scale * self.correction[3 * point + j];
                    }
                }
                match context.evaluate(&mut self.candidate, &mut evaluations, &mut keep_going) {
                    Ok(()) => {
                        let next = self.candidate.worst()?;
                        if next.merit < worst.merit || next.merit <= 1.0 {
                            std::mem::swap(&mut self.current, &mut self.candidate);
                            accepted = true;
                            break;
                        }
                    }
                    Err(e) if outside_trial_domain(&e) => {}
                    Err(e) => return Err(e),
                }
                scale *= 0.5;
            }
            if !accepted {
                return Err(worst.refusal(iteration + 1));
            }
        }
        unreachable!("bounded patch solve returns at its final iteration")
    }
}

struct Context<'p, 'a> {
    areas: &'p [f64],
    compliance: &'p [f64],
    config: CohesiveContactConfig,
    points: &'p [CohesivePoint<'a>],
    free: &'p [[f64; 3]],
    dt: f64,
}
impl Context<'_, '_> {
    fn force(
        &self,
        point: usize,
        jump: [f64; 3],
        evaluations: &mut usize,
        keep_going: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 3], CohesivePatchError> {
        poll(keep_going)?;
        let old = self.points[point].jump_m();
        let increment = core::array::from_fn(|i| jump[i] - old[i]);
        if jump.iter().chain(increment.iter()).any(|x| !x.is_finite())
            || norm_inf(increment) > self.config.maximum_jump_increment_m
        {
            return Err(point_error(
                point,
                CohesiveContactError::Range("per-step cohesive jump increment"),
            ));
        }
        if *evaluations >= self.config.max_material_evaluations {
            return Err(CohesivePatchError::Budget("total point-law evaluations"));
        }
        *evaluations += 1;
        let (_, mean) = self.points[point]
            .trial_over_time_with_mean_traction(
                jump,
                self.dt,
                self.config.material_energy_tolerance_j_m2,
            )
            .map_err(|source| point_error(point, CohesiveContactError::Material(source)))?;
        let force = mean.map(|value| self.areas[point] * value);
        if force.iter().any(|value| !value.is_finite())
            || norm_inf(force) > self.config.maximum_force_n
        {
            return Err(point_error(
                point,
                CohesiveContactError::Range("cohesive held site force"),
            ));
        }
        Ok(force)
    }

    fn evaluate(
        &self,
        out: &mut Evaluation,
        evaluations: &mut usize,
        keep_going: &mut impl FnMut() -> bool,
    ) -> Result<(), CohesivePatchError> {
        for point in 0..self.points.len() {
            out.force[point] = self.force(point, out.jump[point], evaluations, keep_going)?;
        }
        let dimension = 3 * self.points.len();
        for point in 0..self.points.len() {
            poll(keep_going)?;
            let old = self.points[point].jump_m();
            let increment = core::array::from_fn(|j| out.jump[point][j] - old[j]);
            let free_increment = core::array::from_fn(|j| self.free[point][j] - old[j]);
            let mut response = [0.0; 3];
            for j in 0..3 {
                let row = (3 * point + j) * dimension;
                response[j] = checked(
                    (0..dimension)
                        .map(|column| {
                            self.compliance[row + column] * out.force[column / 3][column % 3]
                        })
                        .sum(),
                )?;
                out.residual[point][j] =
                    checked(out.jump[point][j] - self.free[point][j] + response[j])?;
            }
            let scale = norm_inf(increment)
                .max(norm_inf(free_increment))
                .max(norm_inf(response));
            out.tolerance[point] = checked(
                self.config.jump_absolute_tolerance_m + self.config.jump_relative_tolerance * scale,
            )?;
        }
        Ok(())
    }
}

fn point_error(point: usize, source: CohesiveContactError) -> CohesivePatchError {
    CohesivePatchError::Point { point, source }
}
fn outside_trial_domain(error: &CohesivePatchError) -> bool {
    matches!(error, CohesivePatchError::Point { source, .. } if super::outside_trial_domain(source))
}
fn checked(value: f64) -> Result<f64, CohesivePatchError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CohesivePatchError::Invalid("nonfinite patch arithmetic"))
    }
}
fn poll(keep_going: &mut impl FnMut() -> bool) -> Result<(), CohesivePatchError> {
    if keep_going() {
        Ok(())
    } else {
        Err(CohesivePatchError::Cancelled)
    }
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, CohesivePatchError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| CohesivePatchError::Capacity)?;
    result.resize(count, value);
    Ok(result)
}

#[cfg(test)]
mod tests;

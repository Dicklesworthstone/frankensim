//! Simultaneous Signorini--Coulomb impulses with the complete Delassus matrix.
//!
//! The normal law is `0 <= p_n _|_ (v_n - target) >= 0`; tangential impulses
//! maximize dissipation in the disk `|p_t| <= mu*p_n` at the FINAL velocity.
//! This is the non-associated Coulomb law, not an associated cone quadratic
//! (which would introduce normal dilation), nor successive normal/friction clips.
//! All rows participate in one Newton--Krylov natural-map solve. Existing
//! numerical owners provide matrix admission and nonlinear globalization.
//!
//! This bounded first owner admits numerically independent (SPD) contact rows.
//! A weighted contraction witness establishes determinacy within that admitted
//! numerical model; without it, distinct checked roots report nonuniqueness and
//! an otherwise checked candidate remains explicitly unresolved. Singular rigid
//! manifolds, Painleve existence decisions, poses, CCD, gap repair, persistent
//! time stepping, and classical adjoints are not supplied by this event solve.
//! Total midpoint loss is an impact dissipation receipt, NOT a friction-heat
//! allocation: prescribed motion and material energy partition belong downstream.

use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::orthant::active_set::{SpdOrthantConfig, SpdOrthantError, solve_spd_orthant};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use std::cell::{Cell, RefCell};

/// Explicit choice of normal endpoint law; these alternatives are not aliases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NormalVelocityLaw {
    /// Target `-coefficient*min(v_n_before, 0)`, with coefficient in [0, 1].
    NewtonRestitution { coefficient: f64 },
    /// Caller-owned nonnegative normal endpoint target [m/s]. This does not
    /// authorize penetration repair or an undeclared configuration step.
    Target { speed: f64 },
}

/// One candidate point contact in a caller-owned orthonormal, normal-first frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionalContact {
    pub id: u64,
    pub gap: f64,
    pub activation_tolerance: f64,
    /// Relative velocity [normal, tangent 1, tangent 2], m/s. Positive normal
    /// speed separates; positive impulse increases the corresponding speed.
    pub relative_velocity: [f64; 3],
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// Bounded work and separate numerical/physical acceptance controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionalContactConfig {
    /// At most 21 contacts (63 scalar rows), matching the numerical SPD owner.
    pub maximum_contacts: usize,
    pub maximum_evaluations: usize,
    pub maximum_newton_iterations: usize,
    /// Absolute directional-difference perturbation of impulse [N s].
    pub difference_step_impulse: f64,
    pub velocity_tolerance: f64,
    pub impulse_tolerance: f64,
    /// Normal complementarity and tangential duality-gap tolerance [J].
    pub complementarity_tolerance: f64,
    pub energy_tolerance: f64,
    pub matrix: SpdOrthantConfig,
    pub newton: NewtonKrylovConfig,
}
impl Default for FrictionalContactConfig {
    fn default() -> Self {
        Self {
            maximum_contacts: 21,
            maximum_evaluations: 32768,
            maximum_newton_iterations: 96,
            difference_step_impulse: 1e-6,
            velocity_tolerance: 1e-9,
            impulse_tolerance: 1e-9,
            complementarity_tolerance: 1e-9,
            energy_tolerance: 1e-9,
            matrix: SpdOrthantConfig::default(),
            newton: NewtonKrylovConfig {
                relative_tolerance: 1e-12,
                // Directional differences have roundoff/truncation error;
                // oversolving their inner Krylov action can falsely stall an
                // otherwise convergent outer contact solve. Physical gates and
                // the outer 1e-12 absolute target remain independent.
                forcing_minimum: 1e-6,
                ..Default::default()
            },
        }
    }
}

/// Reconstructed law residuals; no natural-map residual substitutes for these.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FrictionalContactResiduals {
    pub normal_impulse_violation: f64,
    pub normal_velocity_violation: f64,
    pub normal_complementarity: f64,
    pub cone_violation: f64,
    pub sliding_impulse_residual: f64,
    pub tangential_duality_gap: f64,
}

/// Complete accepted candidate, in the ORIGINAL caller contact order.
#[derive(Clone, Debug, PartialEq)]
pub struct FrictionalContactResult {
    pub impulses: Vec<[f64; 3]>,
    pub relative_velocities_after: Vec<[f64; 3]>,
    pub normal_targets: Vec<f64>,
    /// Original indices, ordered by stable contact ID, including separating
    /// geometrically active contacts that another contact may drive closed.
    pub active_indices: Vec<usize>,
    pub relative_work: f64,
    /// Nonnegative TOTAL impact loss; no automatic heat or friction partition.
    pub dissipated_energy: f64,
    pub passivity_residual: f64,
    pub residuals: FrictionalContactResiduals,
    /// Weighted block maximum-norm contraction bound, strictly below one.
    /// This is a sufficient floating-point admission, not an interval proof.
    pub contraction_bound: f64,
    /// Positive weights [normal magnitude, tangential Euclidean norm] for each
    /// contact in active_indices order, allowing independent witness checking.
    pub contraction_weights: Vec<f64>,
    pub minimum_normalized_pivot: f64,
    pub evaluations: usize,
    pub nonlinear: Option<NewtonReport>,
}

/// No variant publishes a partially accepted contact event.
#[derive(Debug)]
pub enum FrictionalContactError {
    Invalid(&'static str),
    PenetratingContact {
        index: usize,
        gap: f64,
    },
    /// Singular/near-dependent and indefinite systems retain the numerical
    /// owner's distinct classifications; none means proven physical nonexistence.
    MobilityNotAdmitted(SpdOrthantError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    SolverInterrupted(SolverRunProgress),
    Cancelled,
    EvaluationBudget,
    NonFinite,
    ResidualNotSatisfied(FrictionalContactResiduals),
    EnergyIncrease {
        work: f64,
        allowed: f64,
    },
    /// Both witnesses satisfy the physical law and passivity gates. Values are
    /// in original input order; no arbitrary branch is selected.
    Nonunique {
        first: Vec<[f64; 3]>,
        second: Vec<[f64; 3]>,
    },
    /// A residual-checked passive root is evidence of a candidate, not a
    /// determination. The bounded search did not certify a unique response.
    UnresolvedDeterminacy {
        candidate_impulses: Vec<[f64; 3]>,
    },
}
impl core::fmt::Display for FrictionalContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled frictional contact: {self:?}")
    }
}
impl std::error::Error for FrictionalContactError {}

/// Solve one frozen-geometry event using full row-major `3N x 3N` mobility.
/// Entries have units kg^-1; vectors are normal-first, contact-major. Every
/// mutation is private scratch. Cancellation, numerical failure, energy growth,
/// or unresolved determinacy returns no authority to update physical state.
pub fn solve_frictional_contacts(
    contacts: &[FrictionalContact],
    delassus: &[f64],
    config: FrictionalContactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FrictionalContactResult, FrictionalContactError> {
    use FrictionalContactError as E;
    poll(&mut cancelled)?;
    validate_config(config)?;
    let n = contacts.len();
    if n > config.maximum_contacts
        || delassus.len() != 9 * n * n
        || delassus.iter().any(|v| !v.is_finite())
    {
        return Err(E::Invalid(
            "contact count, mobility shape or finite entries",
        ));
    }
    let mut order: Vec<_> = (0..n).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    if order
        .windows(2)
        .any(|w| contacts[w[0]].id == contacts[w[1]].id)
    {
        return Err(E::Invalid("duplicate contact ID"));
    }
    let mut targets = vec![0.0; n];
    for (i, c) in contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        if !c.gap.is_finite()
            || !c.activation_tolerance.is_finite()
            || c.activation_tolerance < 0.0
            || c.relative_velocity.iter().any(|v| !v.is_finite())
            || !c.friction_coefficient.is_finite()
            || c.friction_coefficient < 0.0
        {
            return Err(E::Invalid("contact geometry, velocity or friction"));
        }
        if c.gap < 0.0 {
            return Err(E::PenetratingContact {
                index: i,
                gap: c.gap,
            });
        }
        targets[i] = match c.normal_law {
            NormalVelocityLaw::NewtonRestitution { coefficient }
                if coefficient.is_finite() && (0.0..=1.0).contains(&coefficient) =>
            {
                -coefficient * c.relative_velocity[0].min(0.0)
            }
            NormalVelocityLaw::Target { speed } if speed.is_finite() && speed >= 0.0 => speed,
            _ => return Err(E::Invalid("normal endpoint law")),
        };
    }
    let active: Vec<_> = order
        .iter()
        .copied()
        .filter(|&i| contacts[i].gap <= contacts[i].activation_tolerance)
        .collect();
    let d = 3 * active.len();
    let mut matrix = Vec::with_capacity(d * d);
    let mut free = Vec::with_capacity(d);
    let mut mu = Vec::with_capacity(active.len());
    for &i in &active {
        free.extend_from_slice(&contacts[i].relative_velocity);
        mu.push(contacts[i].friction_coefficient);
        for a in 0..3 {
            for &j in &active {
                for b in 0..3 {
                    matrix.push(delassus[(3 * i + a) * (3 * n) + 3 * j + b]);
                }
            }
        }
    }
    let admission = solve_spd_orthant(&matrix, &vec![0.0; d], config.matrix, &mut cancelled)
        .map_err(|e| match e {
            SpdOrthantError::Cancelled => E::Cancelled,
            e => E::MobilityNotAdmitted(e),
        })?;
    let mut steps = Vec::with_capacity(d);
    for k in 0..active.len() {
        let a = 3 * k;
        let rn = 1.0 / matrix[a * d + a];
        // One common step preserves the Euclidean Coulomb disk in both axes.
        let rt = 2.0 / (matrix[(a + 1) * d + a + 1] + matrix[(a + 2) * d + a + 2]);
        if !rn.is_finite() || !rt.is_finite() || rn <= 0.0 || rt <= 0.0 {
            return Err(E::NonFinite);
        }
        steps.extend_from_slice(&[rn, rt, rt]);
    }
    let witness = contraction_witness(&matrix, &steps, &mu, &mut cancelled)?;
    let evaluations = Cell::new(0);
    let cancel = RefCell::new(&mut cancelled);
    let solve = |seed: Vec<f64>| -> Result<(Vec<f64>, NewtonReport), E> {
        let problem =
            FallibleFiniteDifference::new(d, config.difference_step_impulse, |p: &[f64]| {
                poll(&mut **cancel.borrow_mut())?;
                if evaluations.get() == config.maximum_evaluations {
                    return Err(E::EvaluationBudget);
                }
                evaluations.set(evaluations.get() + 1);
                let mut residual = vec![0.0; d];
                for k in 0..active.len() {
                    poll(&mut **cancel.borrow_mut())?;
                    let a = 3 * k;
                    let mut v = [0.0; 3];
                    for b in 0..3 {
                        v[b] =
                            free[a + b] + fs_solver::dot(&matrix[(a + b) * d..(a + b + 1) * d], p);
                    }
                    residual[a] = p[a] - (p[a] - steps[a] * (v[0] - targets[active[k]])).max(0.0);
                    let trial = [
                        p[a + 1] - steps[a + 1] * v[1],
                        p[a + 2] - steps[a + 2] * v[2],
                    ];
                    let radius = mu[k] * p[a].max(0.0);
                    let length = trial[0].hypot(trial[1]);
                    if v.iter().chain(&trial).any(|x| !x.is_finite())
                        || !radius.is_finite()
                        || !length.is_finite()
                    {
                        return Err(E::NonFinite);
                    }
                    let fraction = if length > radius {
                        radius / length
                    } else {
                        1.0
                    };
                    residual[a + 1] = p[a + 1] - fraction * trial[0];
                    residual[a + 2] = p[a + 2] - fraction * trial[1];
                }
                if residual.iter().any(|v| !v.is_finite()) {
                    return Err(E::NonFinite);
                }
                Ok(residual)
            })
            .map_err(fd_error)?;
        let construction = NewtonKrylovState::new(&problem, seed, config.newton);
        if let Some(e) = problem.take_failure() {
            return Err(fd_error(e));
        }
        let mut state = construction.map_err(E::Newton)?;
        let report = state.run_controlled(&problem, config.maximum_newton_iterations, &mut || {
            problem.has_failed() || (**cancel.borrow_mut())()
        });
        if let Some(e) = problem.take_failure() {
            return Err(fd_error(e));
        }
        if report.progress == SolverRunProgress::Paused {
            return Err(E::Cancelled);
        }
        if report.progress != SolverRunProgress::Complete {
            return Err(E::SolverInterrupted(report.progress));
        }
        if !report.report.converged {
            return Err(E::NotConverged(report.report));
        }
        Ok((state.x, report.report))
    };
    let (point, report) = if d == 0 {
        (vec![], None)
    } else {
        let (point, report) = solve(vec![0.0; d])?;
        (point, Some(report))
    };
    let mut result = audit(
        contacts,
        delassus,
        &order,
        &active,
        &targets,
        &point,
        config,
        &mut **cancel.borrow_mut(),
    )?;
    if let Some((bound, weights)) = witness {
        result.contraction_bound = bound;
        result.contraction_weights = weights;
    } else {
        // A second, genuinely loaded branch probe can reveal the separating /
        // compressive ambiguity that a zero/frictionless pair of seeds misses.
        let mut seed = vec![0.0; d];
        for k in 0..active.len() {
            let a = 3 * k;
            let tangential_speed = free[a + 1].hypot(free[a + 2]);
            seed[a] = steps[a]
                * (free[a].abs()
                    + targets[active[k]]
                    + tangential_speed
                    + config.velocity_tolerance);
            if tangential_speed > 0.0 {
                seed[a + 1] = -mu[k] * seed[a] * (free[a + 1] / tangential_speed);
                seed[a + 2] = -mu[k] * seed[a] * (free[a + 2] / tangential_speed);
            }
        }
        match solve(seed) {
            Ok((other, _)) => {
                match audit(
                    contacts,
                    delassus,
                    &order,
                    &active,
                    &targets,
                    &other,
                    config,
                    &mut **cancel.borrow_mut(),
                ) {
                    Ok(second)
                        if point
                            .iter()
                            .zip(&other)
                            .any(|(a, b)| (a - b).abs() > 16.0 * config.impulse_tolerance) =>
                    {
                        return Err(E::Nonunique {
                            first: result.impulses,
                            second: second.impulses,
                        });
                    }
                    Err(E::Cancelled) => return Err(E::Cancelled),
                    _ => {}
                }
            }
            Err(E::Cancelled) => return Err(E::Cancelled),
            Err(E::EvaluationBudget) => return Err(E::EvaluationBudget),
            Err(_) => {}
        }
        poll(&mut **cancel.borrow_mut())?;
        return Err(E::UnresolvedDeterminacy {
            candidate_impulses: result.impulses,
        });
    }
    result.minimum_normalized_pivot = admission.min_normalized_pivot;
    result.evaluations = evaluations.get();
    result.nonlinear = report;
    poll(&mut **cancel.borrow_mut())?;
    Ok(result)
}

/// Bound the natural map in contact normal/tangential block norms. Disk
/// projection is nonexpansive in its argument and 1-Lipschitz in its radius.
/// Hence `delta T <= L delta p`; a positive w with Lw < w proves contraction
/// in max_i |delta p_i|/w_i. This is a sufficient condition, not a failed
/// uniqueness test when no witness is found. It performs no impulse iteration.
fn contraction_witness(
    w: &[f64],
    steps: &[f64],
    mu: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<(f64, Vec<f64>)>, FrictionalContactError> {
    let n = mu.len();
    let d = 3 * n;
    let blocks = 2 * n;
    if n == 0 {
        return Ok(Some((0.0, vec![])));
    }
    let mut l = vec![0.0; blocks * blocks];
    for i in 0..n {
        poll(cancelled)?;
        let a = 3 * i;
        let rn = steps[a];
        let rt = steps[a + 1];
        for j in 0..n {
            let b = 3 * j;
            let diagonal = if i == j { 1.0 } else { 0.0 };
            l[(2 * i) * blocks + 2 * j] = (diagonal - rn * w[a * d + b]).abs();
            l[(2 * i) * blocks + 2 * j + 1] = rn * w[a * d + b + 1].hypot(w[a * d + b + 2]);
            l[(2 * i + 1) * blocks + 2 * j] = rt * w[(a + 1) * d + b].hypot(w[(a + 2) * d + b])
                + if i == j { mu[i] } else { 0.0 };
            // Frobenius norm bounds the 2x2 operator norm, avoiding a second
            // eigensolver and preserving a conservative determinacy boundary.
            let x = diagonal - rt * w[(a + 1) * d + b + 1];
            let y = -rt * w[(a + 1) * d + b + 2];
            let z = -rt * w[(a + 2) * d + b + 1];
            let t = diagonal - rt * w[(a + 2) * d + b + 2];
            l[(2 * i + 1) * blocks + 2 * j + 1] = x.hypot(y).hypot(z.hypot(t));
        }
    }
    let mut weights = vec![1.0; blocks];
    for _ in 0..128 {
        poll(cancelled)?;
        let image: Vec<_> = l
            .chunks_exact(blocks)
            .map(|row| fs_solver::dot(row, &weights))
            .collect();
        if image.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
        let bound = image
            .iter()
            .zip(&weights)
            .map(|(a, b)| a / b)
            .fold(0.0, f64::max)
            * (1.0 + 256.0 * f64::EPSILON * (blocks as f64));
        if bound < 1.0 - 1e-8 {
            return Ok(Some((bound, weights)));
        }
        weights = image.into_iter().map(|v| 1.0 + v).collect();
    }
    Ok(None)
}

fn audit(
    contacts: &[FrictionalContact],
    w: &[f64],
    order: &[usize],
    active: &[usize],
    targets: &[f64],
    point: &[f64],
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FrictionalContactResult, FrictionalContactError> {
    use FrictionalContactError as E;
    let n = contacts.len();
    let d = 3 * n;
    let mut impulses = vec![[0.0; 3]; n];
    for (k, &i) in active.iter().enumerate() {
        impulses[i].copy_from_slice(&point[3 * k..3 * k + 3]);
    }
    let mut velocities: Vec<_> = contacts.iter().map(|c| c.relative_velocity).collect();
    for &i in order {
        poll(cancelled)?;
        for a in 0..3 {
            for &j in active {
                for b in 0..3 {
                    velocities[i][a] += w[(3 * i + a) * d + 3 * j + b] * impulses[j][b];
                }
            }
        }
    }
    let mut r = FrictionalContactResiduals::default();
    if impulses
        .iter()
        .chain(&velocities)
        .flatten()
        .any(|v| !v.is_finite())
    {
        return Err(E::NonFinite);
    }
    for &i in active {
        let p = impulses[i];
        let v = velocities[i];
        let g = v[0] - targets[i];
        let radius = contacts[i].friction_coefficient * p[0];
        let speed = v[1].hypot(v[2]);
        if !g.is_finite() || !radius.is_finite() || !speed.is_finite() {
            return Err(E::NonFinite);
        }
        r.normal_impulse_violation = r.normal_impulse_violation.max((-p[0]).max(0.0));
        r.normal_velocity_violation = r.normal_velocity_violation.max((-g).max(0.0));
        r.normal_complementarity = r.normal_complementarity.max((p[0] * g).abs());
        r.cone_violation = r.cone_violation.max((p[1].hypot(p[2]) - radius).max(0.0));
        let duality_gap = (radius * speed + p[1] * v[1] + p[2] * v[2]).abs();
        if !duality_gap.is_finite() {
            return Err(E::NonFinite);
        }
        r.tangential_duality_gap = r.tangential_duality_gap.max(duality_gap);
        if speed > config.velocity_tolerance {
            r.sliding_impulse_residual = r
                .sliding_impulse_residual
                .max((p[1] + radius * (v[1] / speed)).hypot(p[2] + radius * (v[2] / speed)));
        }
    }
    let mut work = 0.0;
    let mut work_scale = 0.0;
    for &i in order {
        for a in 0..3 {
            let term =
                impulses[i][a] * (0.5 * contacts[i].relative_velocity[a] + 0.5 * velocities[i][a]);
            work += term;
            work_scale += term.abs();
        }
    }
    if impulses
        .iter()
        .chain(&velocities)
        .flatten()
        .any(|v| !v.is_finite())
        || !work.is_finite()
        || !work_scale.is_finite()
        || [
            r.normal_impulse_violation,
            r.normal_velocity_violation,
            r.normal_complementarity,
            r.cone_violation,
            r.sliding_impulse_residual,
            r.tangential_duality_gap,
        ]
        .iter()
        .any(|v| !v.is_finite())
    {
        return Err(E::NonFinite);
    }
    if r.normal_impulse_violation > config.impulse_tolerance
        || r.normal_velocity_violation > config.velocity_tolerance
        || r.normal_complementarity > config.complementarity_tolerance
        || r.cone_violation > config.impulse_tolerance
        || r.sliding_impulse_residual > config.impulse_tolerance
        || r.tangential_duality_gap > config.complementarity_tolerance
    {
        return Err(E::ResidualNotSatisfied(r));
    }
    let allowed = config.energy_tolerance + 64.0 * f64::EPSILON * work_scale;
    if work > allowed {
        return Err(E::EnergyIncrease { work, allowed });
    }
    poll(cancelled)?;
    Ok(FrictionalContactResult {
        impulses,
        relative_velocities_after: velocities,
        normal_targets: targets.to_vec(),
        active_indices: active.to_vec(),
        relative_work: work,
        dissipated_energy: (-work).max(0.0),
        passivity_residual: work.max(0.0),
        residuals: r,
        contraction_bound: f64::NAN,
        contraction_weights: vec![],
        minimum_normalized_pivot: f64::NAN,
        evaluations: 0,
        nonlinear: None,
    })
}

fn validate_config(c: FrictionalContactConfig) -> Result<(), FrictionalContactError> {
    if !(1..=21).contains(&c.maximum_contacts)
        || !(1..=1_000_000).contains(&c.maximum_evaluations)
        || !(1..=256).contains(&c.maximum_newton_iterations)
        || [
            c.difference_step_impulse,
            c.velocity_tolerance,
            c.impulse_tolerance,
            c.complementarity_tolerance,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || !c.energy_tolerance.is_finite()
        || c.energy_tolerance < 0.0
        || c.newton.linear_restart > 64
        || c.newton.max_linear_cycles > 128
    {
        return Err(FrictionalContactError::Invalid("configuration"));
    }
    c.newton.validate().map_err(FrictionalContactError::Newton)
}
fn fd_error(e: FiniteDifferenceError<FrictionalContactError>) -> FrictionalContactError {
    match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => FrictionalContactError::Invalid(s),
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FrictionalContactError> {
    if cancelled() {
        Err(FrictionalContactError::Cancelled)
    } else {
        Ok(())
    }
}

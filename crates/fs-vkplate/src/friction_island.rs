//! Simultaneous compliant normal contact and unregularized Coulomb friction.
//!
//! Coordinates are `[normal forces, opposing Coulomb pairs, retained scalar
//! tangential forces]`. The caller
//! supplies the fixed Gram response `C = B^T A^-1 B`, with `B = [N, T]`, and
//! paths `[penetrations, T^T(w_next-w_previous)] = predicted - C q`. Thus the
//! tangent coordinates oppose slip; their physical traction is `-q_t`.
//!
//! Normal constitutive secants remain with the caller. Coulomb residuals and
//! fixed generalized derivatives come from `fs-contact::sustained`, while the
//! shared Newton--FGMRES owner supplies globalization and true linear residuals.
//! The pressure-dependent friction capacity makes this Jacobian nonsymmetric.
//! No positive-definite contact matrix, smoothing, or tangential spring is
//! introduced. A successful candidate establishes the stated discrete laws,
//! not existence or uniqueness of the non-associated friction problem. Retained
//! scalar laws supply their own signed force and analytic tangent from frozen
//! history. Their receiving adapter owns the final storage/work/heat check;
//! elastic energy return is not subjected to a dry-friction work-sign gate.

use std::cell::{Cell, RefCell};

use fs_contact::sustained::{
    ContactLawError, ContactLinearization, ContactPoint, NaturalMapScales, linearize,
};
use fs_solver::{
    Globalization, LineSearchConfig, NewtonError, NewtonKrylovConfig, NewtonKrylovState,
    NewtonStallDiagnosis, NonlinearProblem,
};

use crate::contact_island::RefusalSite;

/// Accuracy and bounded work for one coupled response component.
#[derive(Clone, Copy, Debug)]
pub struct Controls {
    pub relative_tolerance: f64,
    /// Absolute centered-slip tolerance [m].
    pub displacement_tolerance: f64,
    pub force_tolerance: f64,
    /// Total Newton attempts, including at most two force-scale refinements.
    /// The body adapter's rounded-pose retries spend this same budget.
    pub max_newton_iterations: usize,
    /// FGMRES columns per Newton attempt. The restart/cycle product is rounded
    /// down to this budget; no incomplete restart is added beyond it.
    pub max_linear_iterations: usize,
    pub linear_restart: usize,
    pub max_backtracks: usize,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            displacement_tolerance: 1e-14,
            force_tolerance: 1e-12,
            max_newton_iterations: 32,
            max_linear_iterations: 256,
            linear_restart: 24,
            max_backtracks: 24,
        }
    }
}

impl Controls {
    pub(crate) fn validate(self) -> Result<(), SolveError> {
        if !self.relative_tolerance.is_finite()
            || self.relative_tolerance <= 0.0
            || self.relative_tolerance >= 1.0
            || !self.displacement_tolerance.is_finite()
            || self.displacement_tolerance <= 0.0
            || !self.force_tolerance.is_finite()
            || self.force_tolerance <= 0.0
            || self.max_linear_iterations == 0
            || self.linear_restart == 0
            || self.max_backtracks > 1022
        {
            return Err(SolveError::InvalidInput(
                "invalid accuracy or work controls",
            ));
        }
        Ok(())
    }
}

/// One Coulomb pair, linked to its original compliant normal force.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Tangential {
    pub normal: usize,
    pub mu: f64,
}

/// Work counts and independently recomputed physical residuals. Residuals are
/// normalized by the declared absolute-plus-relative tolerances; both must be
/// at most one. Friction also checks capacity, duality gap and work sign.
#[derive(Clone, Copy, Debug, Default)]
pub struct SolveReport {
    pub newton_iterations: usize,
    /// FGMRES columns attempted, including a column that triggers refusal.
    pub linear_iterations: usize,
    pub compliance_applications: usize,
    pub normal_residual: f64,
    pub friction_residual: f64,
    /// Original retained-law force residual at the common final path.
    pub retained_residual: f64,
    /// Signed work on the bodies, recomputed by the adapter from the accepted
    /// rounded nodal centered slip: `-0.5 * q_t dot (w_next-w_previous)`.
    pub friction_work_j: f64,
    /// Nonnegative physical dissipation, filled by the body adapter after its
    /// final rounded-pose check. Retained-law heat comes from its original
    /// constitutive owner, separately from recoverable interface storage.
    pub dissipated_j: f64,
    /// Change in recoverable tangential interface storage [J]. A negative
    /// value can supply positive body work without violating passivity.
    pub stored_change_j: f64,
    /// Candidates the core accepted and the rounded-pose audit then refused
    /// before this one was published. Their Newton work is in the counts above.
    pub rounded_pose_retries: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SolveError {
    InvalidInput(&'static str),
    InvalidLaw {
        contact: usize,
    },
    InvalidResponse,
    Coulomb {
        friction: usize,
        error: ContactLawError,
    },
    Mindlin {
        contact: usize,
        error: fs_contact::interface_binding::mindlin::SourceMindlinError,
    },
    Newton(NewtonError),
    NoConvergence {
        site: RefusalSite,
        newton_iterations: usize,
        normal_residual: f64,
        friction_residual: f64,
        retained_residual: f64,
        diagnosis: Option<NewtonStallDiagnosis>,
    },
}

impl std::fmt::Display for SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(reason) => write!(f, "joint contact: {reason}"),
            Self::InvalidLaw { contact } => {
                write!(f, "joint contact {contact}: invalid normal force or slope")
            }
            Self::InvalidResponse => {
                write!(f, "joint contact: non-finite body response or residual")
            }
            Self::Coulomb { friction, error } => write!(f, "joint friction {friction}: {error}"),
            Self::Mindlin { contact, error } => {
                write!(f, "joint Mindlin contact {contact}: {error}")
            }
            Self::Newton(error) => write!(f, "joint contact: {error}"),
            Self::NoConvergence {
                site,
                newton_iterations,
                normal_residual,
                friction_residual,
                retained_residual,
                diagnosis,
            } => write!(
                f,
                "joint contact refused at {site:?} after {newton_iterations} Newton attempts: scaled normal {normal_residual:e}, friction {friction_residual:e}, retained {retained_residual:e}, {diagnosis:?}"
            ),
        }
    }
}

impl std::error::Error for SolveError {}

/// The same Coulomb acceptance and work convention used by the core and its
/// receiving body's final rounded-pose check.
#[derive(Clone, Copy, Debug)]
pub(crate) struct FrictionCheck {
    pub residual: f64,
    pub friction_work_j: f64,
    pub work_tolerance_j: f64,
}

fn tangent_step(diagonal_bound: [f64; 2]) -> Result<f64, ContactLawError> {
    if diagonal_bound.iter().any(|c| !c.is_finite() || *c < 0.0) {
        return Err(ContactLawError::Invalid(
            "nonnegative finite tangent response bounds required",
        ));
    }
    // Halving first avoids overflow in the sum of finite bounds.
    let average = 0.5 * diagonal_bound[0] + 0.5 * diagonal_bound[1];
    let step = if average > 0.0 { average.recip() } else { 1.0 };
    if step.is_finite() && step > 0.0 {
        Ok(step)
    } else {
        Err(ContactLawError::NonFinite("tangent natural-map scale"))
    }
}

/// What rounding the pose can add to the natural residual: the slip's rounding
/// times the tangent step, times the disk projection's slope, which is one
/// inside the disk and `capacity / |trial|` outside it. No iteration can remove
/// this, and on a lightly loaded contact far from rest it can exceed the force
/// tolerance.
fn pose_rounding_force(
    opposing_force: [f64; 2],
    centered_slip: [f64; 2],
    slip_rounding: [f64; 2],
    capacity: f64,
    step: f64,
) -> f64 {
    let error = step * slip_rounding[0].hypot(slip_rounding[1]);
    let trial = (step * centered_slip[0] + opposing_force[0])
        .hypot(step * centered_slip[1] + opposing_force[1]);
    if trial - error > capacity {
        error * capacity / (trial - error)
    } else {
        error
    }
}

#[allow(clippy::too_many_arguments)]
fn friction_audit(
    normal_force: f64,
    opposing_force: [f64; 2],
    centered_slip: [f64; 2],
    predicted_slip: [f64; 2],
    mu: f64,
    natural_residual: [f64; 2],
    pose_rounding_n: f64,
    controls: Controls,
) -> Result<FrictionCheck, ContactLawError> {
    let magnitude = opposing_force[0].hypot(opposing_force[1]);
    let slip = centered_slip[0].hypot(centered_slip[1]);
    let capacity = mu * normal_force.max(0.0);
    let force_scale =
        controls.force_tolerance + controls.relative_tolerance * magnitude.max(capacity);
    let slip_scale = controls.displacement_tolerance
        + controls.relative_tolerance * slip.max(predicted_slip[0].hypot(predicted_slip[1]));
    let work = opposing_force[0].mul_add(centered_slip[0], opposing_force[1] * centered_slip[1]);
    let duality_gap = capacity.mul_add(slip, -work);
    let work_scale = force_scale * slip
        + slip_scale * magnitude.max(capacity)
        + controls.force_tolerance * controls.displacement_tolerance;
    let residual = natural_residual[0].hypot(natural_residual[1]) / (force_scale + pose_rounding_n);
    let cone = (magnitude - capacity).max(0.0) / force_scale;
    let duality = duality_gap.abs() / work_scale;
    let wrong_work = (-work).max(0.0) / work_scale;
    if [
        force_scale,
        slip_scale,
        work_scale,
        residual,
        cone,
        duality,
        wrong_work,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || force_scale <= 0.0
        || work_scale <= 0.0
    {
        return Err(ContactLawError::NonFinite(
            "Coulomb physical residual or work budget",
        ));
    }
    Ok(FrictionCheck {
        residual: residual.max(cone).max(duality).max(wrong_work),
        friction_work_j: -0.5 * work,
        work_tolerance_j: 0.5 * work_scale,
    })
}

/// Re-evaluate the original Coulomb law at an independently formed final
/// path. A result is admitted only when `residual <= 1` and the caller also
/// admits its original compliant normal law at that same path. `slip_rounding`
/// bounds what rounding the pose put into `centered_slip`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn check_friction(
    normal_force: f64,
    opposing_force: [f64; 2],
    centered_slip: [f64; 2],
    predicted_slip: [f64; 2],
    mu: f64,
    diagonal_bound: [f64; 2],
    slip_rounding: [f64; 2],
    controls: Controls,
) -> Result<FrictionCheck, ContactLawError> {
    if predicted_slip.iter().any(|v| !v.is_finite()) {
        return Err(ContactLawError::Invalid("finite predicted slip required"));
    }
    let step = tangent_step(diagonal_bound)?;
    let block = linearize(
        ContactPoint {
            force_n: [normal_force, -opposing_force[0], -opposing_force[1]],
            gap_m: 0.0,
            slip_increment_m: centered_slip,
            friction_coefficient: mu,
        },
        NaturalMapScales {
            normal_n_per_m: 1.0,
            tangent_n_per_m: step,
        },
    )?;
    friction_audit(
        normal_force,
        opposing_force,
        centered_slip,
        predicted_slip,
        mu,
        [block.residual_n[1], block.residual_n[2]],
        pose_rounding_force(
            opposing_force,
            centered_slip,
            slip_rounding,
            mu * normal_force.max(0.0),
            step,
        ),
        controls,
    )
}

#[derive(Clone, Default)]
struct Cache {
    valid: bool,
    point: Vec<f64>,
    paths: Vec<f64>,
    residual: Vec<f64>,
    normal_force: Vec<f64>,
    slope: Vec<f64>,
    retained_force: Vec<f64>,
    retained_slope: Vec<f64>,
    tangents: Vec<ContactLinearization>,
}

/// Linear contact storage retained across samples. The shared Newton--FGMRES
/// implementation owns its Krylov workspace. Numerical continuation uses the
/// supplied previous force; no branch search or uniqueness claim is made.
#[derive(Clone, Default)]
pub(crate) struct Workspace {
    accepted: bool,
    force: Vec<f64>,
    paths: Vec<f64>,
    scales: Vec<f64>,
    tangent_steps: Vec<f64>,
    action_response: Vec<f64>,
    cache: Cache,
}

fn apply_response(
    apply: &mut impl FnMut(&[f64], &mut [f64]),
    applications: &Cell<usize>,
    x: &[f64],
    y: &mut [f64],
) -> Result<(), SolveError> {
    if x.iter().any(|v| !v.is_finite()) {
        return Err(SolveError::InvalidResponse);
    }
    if x.iter().all(|&v| v == 0.0) {
        y.fill(0.0);
        return Ok(());
    }
    applications.set(applications.get() + 1);
    y.fill(f64::NAN);
    apply(x, y);
    if y.iter().any(|v| !v.is_finite()) {
        return Err(SolveError::InvalidResponse);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    cache: &mut Cache,
    x: &[f64],
    predicted: &[f64],
    normal_count: usize,
    tangents: &[Tangential],
    steps: &[f64],
    apply: &mut impl FnMut(&[f64], &mut [f64]),
    normal_law: &mut impl FnMut(usize, f64) -> (f64, f64),
    retained_law: &mut impl FnMut(usize, f64) -> Result<(f64, f64), SolveError>,
    applications: &Cell<usize>,
) -> Result<(), SolveError> {
    if cache.valid
        && cache
            .point
            .iter()
            .zip(x)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    {
        return Ok(());
    }
    cache.valid = false;
    apply_response(apply, applications, x, &mut cache.paths)?;
    for (path, &free) in cache.paths.iter_mut().zip(predicted) {
        *path = free - *path;
        if !path.is_finite() {
            return Err(SolveError::InvalidResponse);
        }
    }
    for i in 0..normal_count {
        let (force, slope) = normal_law(i, cache.paths[i]);
        if !force.is_finite() || force < 0.0 || !slope.is_finite() || slope < 0.0 {
            return Err(SolveError::InvalidLaw { contact: i });
        }
        cache.normal_force[i] = force;
        cache.slope[i] = slope;
        cache.residual[i] = x[i] - force;
    }
    cache.tangents.clear();
    for (j, tangent) in tangents.iter().enumerate() {
        let a = normal_count + 2 * j;
        let block = linearize(
            ContactPoint {
                force_n: [x[tangent.normal], -x[a], -x[a + 1]],
                // The rigid normal row is unused: the compliant law above owns it.
                // A zero placeholder avoids irrelevant normal-trial overflow.
                gap_m: 0.0,
                slip_increment_m: [cache.paths[a], cache.paths[a + 1]],
                friction_coefficient: tangent.mu,
            },
            NaturalMapScales {
                normal_n_per_m: 1.0,
                tangent_n_per_m: steps[j],
            },
        )
        .map_err(|error| SolveError::Coulomb { friction: j, error })?;
        cache.residual[a] = -block.residual_n[1];
        cache.residual[a + 1] = -block.residual_n[2];
        cache.tangents.push(block);
    }
    let start = normal_count + 2 * tangents.len();
    for j in 0..cache.retained_force.len() {
        // The caller returns dF/d(centered path), including any midpoint
        // conversion. Unlike a unilateral reaction, stored shear can be signed.
        let (force, slope) = retained_law(j, cache.paths[start + j])?;
        if !force.is_finite() || !slope.is_finite() || slope < 0.0 {
            return Err(SolveError::InvalidResponse);
        }
        cache.retained_force[j] = force;
        cache.retained_slope[j] = slope;
        cache.residual[start + j] = x[start + j] - force;
    }
    if cache.residual.iter().any(|v| !v.is_finite()) {
        return Err(SolveError::InvalidResponse);
    }
    cache.point.copy_from_slice(x);
    cache.valid = true;
    Ok(())
}

/// Check physical force laws independently of numerical row scaling and the
/// Newton norm. The duality gap prevents a projection scale from concealing
/// finite sticking slip or a sliding force in the wrong direction.
fn audit(
    cache: &Cache,
    predicted: &[f64],
    normal_count: usize,
    tangents: &[Tangential],
    controls: Controls,
) -> Result<(f64, f64, f64, bool), SolveError> {
    let mut normal = 0.0_f64;
    let mut friction = 0.0_f64;
    let mut retained = 0.0_f64;
    let mut feasible = true;
    for i in 0..normal_count {
        let q = cache.point[i];
        feasible &= q >= 0.0;
        let scale = controls.force_tolerance
            + controls.relative_tolerance * q.abs().max(cache.normal_force[i]);
        let residual = cache.residual[i].abs() / scale;
        if !scale.is_finite() || !residual.is_finite() {
            return Err(SolveError::InvalidResponse);
        }
        normal = normal.max(residual);
    }
    for (j, tangent) in tangents.iter().enumerate() {
        let a = normal_count + 2 * j;
        let checked = friction_audit(
            cache.point[tangent.normal],
            [cache.point[a], cache.point[a + 1]],
            [cache.paths[a], cache.paths[a + 1]],
            [predicted[a], predicted[a + 1]],
            tangent.mu,
            [cache.residual[a], cache.residual[a + 1]],
            0.0,
            controls,
        )
        .map_err(|error| SolveError::Coulomb { friction: j, error })?;
        friction = friction.max(checked.residual);
    }
    let start = normal_count + 2 * tangents.len();
    for (j, &expected) in cache.retained_force.iter().enumerate() {
        let q = cache.point[start + j];
        let scale =
            controls.force_tolerance + controls.relative_tolerance * q.abs().max(expected.abs());
        let residual = cache.residual[start + j].abs() / scale;
        if !scale.is_finite() || !residual.is_finite() {
            return Err(SolveError::InvalidResponse);
        }
        retained = retained.max(residual);
    }
    Ok((normal, friction, retained, feasible))
}

struct Problem<'a, F, L, R> {
    predicted: &'a [f64],
    normal_count: usize,
    tangents: &'a [Tangential],
    steps: &'a [f64],
    scales: &'a [f64],
    diagonal_bound: &'a [f64],
    apply: RefCell<&'a mut F>,
    normal_law: RefCell<&'a mut L>,
    retained_law: RefCell<&'a mut R>,
    applications: &'a Cell<usize>,
    linear_iterations: &'a Cell<usize>,
    cache: RefCell<&'a mut Cache>,
    action_response: RefCell<&'a mut [f64]>,
    failure: RefCell<Option<SolveError>>,
}

impl<F, L, R> Problem<'_, F, L, R>
where
    F: FnMut(&[f64], &mut [f64]),
    L: FnMut(usize, f64) -> (f64, f64),
    R: FnMut(usize, f64) -> Result<(f64, f64), SolveError>,
{
    fn fail(&self, error: SolveError) {
        let mut failure = self.failure.borrow_mut();
        if failure.is_none() {
            *failure = Some(error);
        }
    }

    fn ensure(&self, x: &[f64]) -> bool {
        if self.failure.borrow().is_some() {
            return false;
        }
        if x.len() != self.predicted.len() || x.iter().any(|v| !v.is_finite()) {
            self.fail(SolveError::InvalidResponse);
            return false;
        }
        let result = evaluate(
            &mut self.cache.borrow_mut(),
            x,
            self.predicted,
            self.normal_count,
            self.tangents,
            self.steps,
            &mut **self.apply.borrow_mut(),
            &mut **self.normal_law.borrow_mut(),
            &mut **self.retained_law.borrow_mut(),
            self.applications,
        );
        if let Err(error) = result {
            self.fail(error);
            false
        } else {
            true
        }
    }

    fn physical_residuals(&self, x: &[f64], controls: Controls) -> Option<(f64, f64, f64, bool)> {
        if !self.ensure(x) {
            return None;
        }
        match audit(
            &self.cache.borrow(),
            self.predicted,
            self.normal_count,
            self.tangents,
            controls,
        ) {
            Ok(result) => Some(result),
            Err(error) => {
                self.fail(error);
                None
            }
        }
    }
}

impl<F, L, R> NonlinearProblem for Problem<'_, F, L, R>
where
    F: FnMut(&[f64], &mut [f64]),
    L: FnMut(usize, f64) -> (f64, f64),
    R: FnMut(usize, f64) -> Result<(f64, f64), SolveError>,
{
    fn dimension(&self) -> usize {
        self.predicted.len()
    }

    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.ensure(x) {
            return;
        }
        let cache = self.cache.borrow();
        for i in 0..output.len() {
            output[i] = cache.residual[i] / self.scales[i];
        }
        if output.iter().any(|v| !v.is_finite()) {
            self.fail(SolveError::InvalidResponse);
        }
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.ensure(x) {
            return;
        }
        let mut response = self.action_response.borrow_mut();
        if let Err(error) = apply_response(
            &mut **self.apply.borrow_mut(),
            self.applications,
            direction,
            &mut response,
        ) {
            self.fail(error);
            return;
        }
        let cache = self.cache.borrow();
        for i in 0..self.normal_count {
            output[i] = cache.slope[i].mul_add(response[i], direction[i]) / self.scales[i];
        }
        for (j, tangent) in self.tangents.iter().enumerate() {
            let a = self.normal_count + 2 * j;
            let action = cache.tangents[j].action(
                [direction[tangent.normal], -direction[a], -direction[a + 1]],
                [0.0, -response[a], -response[a + 1]],
            );
            match action {
                Ok(action) => {
                    output[a] = -action[1] / self.scales[a];
                    output[a + 1] = -action[2] / self.scales[a + 1];
                }
                Err(error) => {
                    self.fail(SolveError::Coulomb { friction: j, error });
                    return;
                }
            }
        }
        let start = self.normal_count + 2 * self.tangents.len();
        for (j, &slope) in cache.retained_slope.iter().enumerate() {
            let i = start + j;
            output[i] = slope.mul_add(response[i], direction[i]) / self.scales[i];
        }
        if output.iter().any(|v| !v.is_finite()) {
            self.fail(SolveError::InvalidResponse);
        }
    }

    fn preconditioner_apply(
        &self,
        x: &[f64],
        _: usize,
        _: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        self.linear_iterations.set(self.linear_iterations.get() + 1);
        output.fill(f64::NAN);
        if !self.ensure(x) {
            return;
        }
        let cache = self.cache.borrow();
        for i in 0..self.normal_count {
            let diagonal = cache.slope[i].mul_add(self.diagonal_bound[i], 1.0);
            output[i] = residual[i] * (self.scales[i] / diagonal);
        }
        for (j, block) in cache.tangents.iter().enumerate() {
            for axis in 0..2 {
                let a = self.normal_count + 2 * j + axis;
                let diagonal = block.path_derivative[axis + 1][axis + 1].mul_add(
                    self.diagonal_bound[a],
                    block.force_derivative[axis + 1][axis + 1],
                );
                // A locked sticking direction can have a zero Jacobian row.
                // Identity scaling does not alter the physical response or
                // manufacture a tangential stiffness in that force gauge.
                let inverse = if diagonal > 0.0 {
                    self.scales[a] / diagonal
                } else {
                    self.scales[a]
                };
                output[a] = residual[a] * inverse;
            }
        }
        let start = self.normal_count + 2 * self.tangents.len();
        for (j, &slope) in cache.retained_slope.iter().enumerate() {
            let i = start + j;
            let diagonal = slope.mul_add(self.diagonal_bound[i], 1.0);
            output[i] = residual[i] * (self.scales[i] / diagonal);
        }
        if output.iter().any(|v| !v.is_finite()) {
            self.fail(SolveError::InvalidResponse);
        }
    }
}

impl Workspace {
    pub(crate) fn forces(&self) -> &[f64] {
        if self.accepted { &self.force } else { &[] }
    }
    pub(crate) fn paths(&self) -> &[f64] {
        if self.accepted { &self.paths } else { &[] }
    }

    /// Solve the declared joint laws without publishing any body or material
    /// state. The callback is one fixed symmetric positive-semidefinite Gram
    /// action for this call; normal laws and histories must stay frozen.
    /// `acceptance` is the largest scaled residual, at most one, this call may
    /// return. It tightens the iteration only; `Controls` still defines the
    /// laws' tolerance, which the rounded-pose audit applies at one.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve(
        &mut self,
        predicted: &[f64],
        previous_force: &[f64],
        compliance_diagonal_bound: &[f64],
        normal_count: usize,
        tangents: &[Tangential],
        retained_count: usize,
        controls: Controls,
        acceptance: f64,
        mut apply_compliance: impl FnMut(&[f64], &mut [f64]),
        mut normal_law: impl FnMut(usize, f64) -> (f64, f64),
        mut retained_law: impl FnMut(usize, f64) -> Result<(f64, f64), SolveError>,
    ) -> Result<SolveReport, SolveError> {
        self.accepted = false;
        controls.validate()?;
        if !(acceptance > 0.0 && acceptance <= 1.0) {
            return Err(SolveError::InvalidInput(
                "acceptance must lie in (0, 1]",
            ));
        }
        let dimension = tangents
            .len()
            .checked_mul(2)
            .and_then(|n| n.checked_add(normal_count))
            .and_then(|n| n.checked_add(retained_count))
            .ok_or(SolveError::InvalidInput("contact dimensions overflow"))?;
        if predicted.len() != dimension
            || previous_force.len() != dimension
            || compliance_diagonal_bound.len() != dimension
            || predicted
                .iter()
                .chain(previous_force)
                .any(|v| !v.is_finite())
            || previous_force[..normal_count].iter().any(|q| *q < 0.0)
            || compliance_diagonal_bound
                .iter()
                .any(|c| !c.is_finite() || *c < 0.0)
            || tangents
                .iter()
                .any(|t| t.normal >= normal_count || !t.mu.is_finite() || t.mu < 0.0)
        {
            return Err(SolveError::InvalidInput(
                "invalid coordinates, friction, warm forces or diagonal bounds",
            ));
        }
        for v in [
            &mut self.force,
            &mut self.paths,
            &mut self.scales,
            &mut self.action_response,
            &mut self.cache.point,
            &mut self.cache.paths,
            &mut self.cache.residual,
        ] {
            v.resize(dimension, 0.0);
        }
        self.cache.normal_force.resize(normal_count, 0.0);
        self.cache.slope.resize(normal_count, 0.0);
        self.cache.retained_force.resize(retained_count, 0.0);
        self.cache.retained_slope.resize(retained_count, 0.0);
        self.tangent_steps.resize(tangents.len(), 0.0);
        for (j, step) in self.tangent_steps.iter_mut().enumerate() {
            let a = normal_count + 2 * j;
            *step = tangent_step([
                compliance_diagonal_bound[a],
                compliance_diagonal_bound[a + 1],
            ])
            .map_err(|error| SolveError::Coulomb { friction: j, error })?;
        }
        self.force.copy_from_slice(previous_force);
        self.cache.valid = false;
        let applications = Cell::new(0);
        let linear_iterations = Cell::new(0);
        let mut report = SolveReport::default();
        let mut diagnosis = None;
        let restart = controls
            .linear_restart
            .min(controls.max_linear_iterations)
            .min(dimension.max(1));
        let admitted = |normal: f64, friction: f64, retained: f64, feasible: bool| {
            feasible && normal <= acceptance && friction <= acceptance && retained <= acceptance
        };
        let config = NewtonKrylovConfig {
            absolute_tolerance: 0.125 * acceptance,
            relative_tolerance: 0.0,
            linear_restart: restart,
            max_linear_cycles: controls.max_linear_iterations / restart,
            forcing_minimum: 1e-8,
            globalization: Globalization::LineSearch(LineSearchConfig {
                minimum_step: 2.0_f64.powi(-(controls.max_backtracks as i32)),
                ..LineSearchConfig::default()
            }),
            ..NewtonKrylovConfig::default()
        };
        for refinement in 0..=2 {
            evaluate(
                &mut self.cache,
                &self.force,
                predicted,
                normal_count,
                tangents,
                &self.tangent_steps,
                &mut apply_compliance,
                &mut normal_law,
                &mut retained_law,
                &applications,
            )?;
            let (normal, friction, retained, feasible) =
                audit(&self.cache, predicted, normal_count, tangents, controls)?;
            report.normal_residual = normal;
            report.friction_residual = friction;
            report.retained_residual = retained;
            report.compliance_applications = applications.get();
            if admitted(normal, friction, retained, feasible) {
                self.paths.copy_from_slice(&self.cache.paths);
                self.accepted = true;
                return Ok(report);
            }
            if report.newton_iterations == controls.max_newton_iterations {
                diagnosis = Some(NewtonStallDiagnosis::BudgetExhausted);
                break;
            }
            let largest_normal = self
                .cache
                .normal_force
                .iter()
                .copied()
                .chain(self.force[..normal_count].iter().copied().map(f64::abs))
                .fold(0.0, f64::max);
            for i in 0..normal_count {
                // The first solve gives initially open contacts the same
                // load scale as their neighbours. Refinement then uses each
                // candidate's own force for the stricter physical gate.
                let reference = if refinement == 0 {
                    largest_normal
                } else {
                    self.force[i].abs().max(self.cache.normal_force[i])
                };
                self.scales[i] = controls.force_tolerance + controls.relative_tolerance * reference;
            }
            for (j, t) in tangents.iter().enumerate() {
                let a = normal_count + 2 * j;
                let normal_reference = if refinement == 0 {
                    largest_normal
                } else {
                    self.force[t.normal].max(self.cache.normal_force[t.normal])
                };
                let reference = self.force[a]
                    .hypot(self.force[a + 1])
                    .max(t.mu * normal_reference);
                let scale = controls.force_tolerance + controls.relative_tolerance * reference;
                self.scales[a] = scale;
                self.scales[a + 1] = scale;
            }
            let start = normal_count + 2 * tangents.len();
            for j in 0..retained_count {
                let reference = self.force[start + j]
                    .abs()
                    .max(self.cache.retained_force[j].abs());
                self.scales[start + j] =
                    controls.force_tolerance + controls.relative_tolerance * reference;
            }
            if self.scales.iter().any(|s| !s.is_finite() || *s <= 0.0) {
                return Err(SolveError::InvalidResponse);
            }
            let problem = Problem {
                predicted,
                normal_count,
                tangents,
                steps: &self.tangent_steps,
                scales: &self.scales,
                diagonal_bound: compliance_diagonal_bound,
                apply: RefCell::new(&mut apply_compliance),
                normal_law: RefCell::new(&mut normal_law),
                retained_law: RefCell::new(&mut retained_law),
                applications: &applications,
                linear_iterations: &linear_iterations,
                cache: RefCell::new(&mut self.cache),
                action_response: RefCell::new(&mut self.action_response),
                failure: RefCell::new(None),
            };
            let construction = NewtonKrylovState::new(&problem, self.force.clone(), config);
            if let Some(error) = problem.failure.borrow_mut().take() {
                return Err(error);
            }
            let mut state = construction.map_err(SolveError::Newton)?;
            let remaining = controls.max_newton_iterations - report.newton_iterations;
            for _ in 0..remaining {
                // Count attempted solves even when Krylov refuses before the
                // shared owner's accepted-step counter advances. Refinement
                // never resets the caller's total Newton work budget.
                report.newton_iterations += 1;
                let native = state.run(&problem, 1);
                diagnosis = native.diagnosis;
                if let Some((normal, friction, retained, feasible)) =
                    problem.physical_residuals(&state.x, controls)
                {
                    report.normal_residual = normal;
                    report.friction_residual = friction;
                    report.retained_residual = retained;
                    if admitted(normal, friction, retained, feasible) {
                        break;
                    }
                }
                if problem.failure.borrow().is_some()
                    || native.converged
                    || native
                        .diagnosis
                        .is_some_and(|d| d != NewtonStallDiagnosis::BudgetExhausted)
                {
                    break;
                }
            }
            report.linear_iterations = linear_iterations.get();
            if let Some(error) = problem.failure.borrow_mut().take() {
                return Err(error);
            }
            drop(problem);
            self.force.copy_from_slice(&state.x);
            // Restore unilateral feasibility before refinement. This is only
            // a new candidate: its complete response and every original law
            // must pass again before any result can be published.
            let mut polished = false;
            for i in 0..normal_count {
                if self.force[i] < 0.0 {
                    self.force[i] = self.cache.normal_force[i];
                    polished = true;
                }
            }
            if polished {
                self.cache.valid = false;
            }
        }
        // A final replay is needed even when the last Newton step spent the
        // remaining budget. Its physical result, not the stop label, decides.
        evaluate(
            &mut self.cache,
            &self.force,
            predicted,
            normal_count,
            tangents,
            &self.tangent_steps,
            &mut apply_compliance,
            &mut normal_law,
            &mut retained_law,
            &applications,
        )?;
        let (normal, friction, retained, feasible) =
            audit(&self.cache, predicted, normal_count, tangents, controls)?;
        report.normal_residual = normal;
        report.friction_residual = friction;
        report.retained_residual = retained;
        report.compliance_applications = applications.get();
        if admitted(normal, friction, retained, feasible) {
            self.paths.copy_from_slice(&self.cache.paths);
            self.accepted = true;
            Ok(report)
        } else {
            Err(SolveError::NoConvergence {
                site: RefusalSite::CoreAudit,
                newton_iterations: report.newton_iterations,
                normal_residual: normal,
                friction_residual: friction,
                retained_residual: retained,
                diagnosis,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense(c: &[[f64; 3]; 3], q: &[f64], y: &mut [f64]) {
        for i in 0..3 {
            y[i] = (0..3).map(|j| c[i][j] * q[j]).sum();
        }
    }

    #[test]
    fn compliant_normal_and_sliding_or_sticking_share_one_final_pose() {
        let c = [[1e-6, 1e-6, 0.0], [1e-6, 2e-6, 0.0], [0.0, 0.0, 1e-6]];
        for (mu, normal, tangent) in [
            (0.12, 1000.0 / 1.56, 1000.0 / 1.56 * 0.12),
            (0.2, 640.0, 80.0),
        ] {
            let mut workspace = Workspace::default();
            let report = workspace
                .solve(
                    &[0.001, 0.0008, 0.0],
                    &[0.0; 3],
                    &[1e-6, 2e-6, 1e-6],
                    1,
                    &[Tangential { normal: 0, mu }],
                    0,
                    Controls::default(),
                    1.0,
                    |q, y| dense(&c, q, y),
                    |_, d| (0.5e6 * (d + 0.001), 0.5e6),
                    |_, _| unreachable!(),
                )
                .unwrap();
            assert!((workspace.forces()[0] - normal).abs() < 2e-7);
            assert!((workspace.forces()[1] - tangent).abs() < 2e-7);
            assert!(report.normal_residual <= 1.0 && report.friction_residual <= 1.0);
            if mu == 0.2 {
                assert!(workspace.paths()[1].abs() < 2e-13);
            } else {
                assert!(workspace.paths()[1] > 0.0);
            }
        }
    }

    #[test]
    fn locked_tangent_direction_needs_no_invented_compliance() {
        let mut workspace = Workspace::default();
        workspace
            .solve(
                &[0.001, 0.0, 0.0],
                &[0.0; 3],
                &[1e-6, 0.0, 0.0],
                1,
                &[Tangential { normal: 0, mu: 0.4 }],
                0,
                Controls::default(),
                1.0,
                |q, y| {
                    y[0] = 1e-6 * q[0];
                    y[1] = 0.0;
                    y[2] = 0.0;
                },
                |_, d| (1e6 * d.max(0.0), if d > 0.0 { 1e6 } else { 0.0 }),
                |_, _| unreachable!(),
            )
            .unwrap();
        assert!((workspace.forces()[0] - 500.0).abs() < 1e-7);
        assert_eq!(&workspace.forces()[1..], &[0.0, 0.0]);
    }

    #[test]
    fn budget_refusal_hides_an_earlier_successful_candidate() {
        let mut workspace = Workspace::default();
        let c = [[1e-6, 0.0, 0.0], [0.0, 1e-6, 0.0], [0.0, 0.0, 1e-6]];
        let law = |_: usize, d: f64| (1e6 * d.max(0.0), if d > 0.0 { 1e6 } else { 0.0 });
        let t = [Tangential { normal: 0, mu: 0.2 }];
        workspace
            .solve(
                &[0.001, 0.001, 0.0],
                &[0.0; 3],
                &[1e-6; 3],
                1,
                &t,
                0,
                Controls::default(),
                1.0,
                |q, y| dense(&c, q, y),
                law,
                |_, _| unreachable!(),
            )
            .unwrap();
        let refused = workspace.solve(
            &[0.001, 0.001, 0.0],
            &[0.0; 3],
            &[1e-6; 3],
            1,
            &t,
            0,
            Controls {
                max_newton_iterations: 0,
                ..Controls::default()
            },
            1.0,
            |q, y| dense(&c, q, y),
            law,
            |_, _| unreachable!(),
        );
        assert!(matches!(
            refused,
            Err(SolveError::NoConvergence {
                site: RefusalSite::CoreAudit,
                ..
            })
        ));
        assert!(workspace.forces().is_empty() && workspace.paths().is_empty());
    }
}

//! Full tangential body response without replacing the constitutive law.
//!
//! The actual load line is `du = du_free - C F`, with symmetric positive
//! semidefinite `C`. A scalar-compliance call to the original law evaluates
//! `F = R_c(du_free + (c I - C) F)`. Its scalar `c` is a numerical probe only:
//! the original matrix, solved slip, retained state and work are checked before
//! a candidate is returned. No physical compliance is averaged or regularized.
//!
//! A full-matrix zero-slip reaction is tried through the original law first,
//! preserving its static priority when a distinct kinetic solution also exists.
//! Bounded inverse-Broyden iterations with backtracking accelerate the remaining
//! two-dimensional equation. Every trial starts from the same accepted history.
//! Nonmonotone laws keep their own native branch selection; this adapter does
//! not establish uniqueness or silently accept a failed constitutive solve.

use core::fmt;

use fs_math::det;

use super::{TangentialLaw, TangentialStep, TangentialStepInput};
use crate::TriboError;

/// Explicit numerical budgets for the full two-axis body/constitutive solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MatrixTangentialConfig {
    /// Maximum original-law evaluations, including refused/backtracked trials.
    /// Between 1 and 4096; each call retains the original inner solve budget.
    pub maximum_evaluations: usize,
    /// Euclidean residual tolerance for `du + C F - du_free` [m].
    pub displacement_tolerance_m: f64,
    /// Absolute native-storage, native-work and actual-body-work tolerance [J].
    pub energy_tolerance_j: f64,
    /// Relative arithmetic symmetry/PSD tolerance, between zero and 1e-8.
    /// Admitted roundoff never replaces entries in the physical matrix.
    pub matrix_relative_tolerance: f64,
}

impl MatrixTangentialConfig {
    /// Check budgets independently of any contact or candidate history.
    pub fn validate(self) -> Result<(), MatrixTangentialError> {
        if !(1..=4096).contains(&self.maximum_evaluations)
            || !self.displacement_tolerance_m.is_finite()
            || self.displacement_tolerance_m <= 0.0
            || !self.energy_tolerance_j.is_finite()
            || self.energy_tolerance_j <= 0.0
            || !self.matrix_relative_tolerance.is_finite()
            || !(0.0..=1e-8).contains(&self.matrix_relative_tolerance)
        {
            return Err(MatrixTangentialError::Invalid("matrix tangential budgets"));
        }
        Ok(())
    }
}

/// Full original-law candidate and its actual matrix/body work residual.
#[derive(Clone, Debug, PartialEq)]
pub struct MatrixTangentialStep<S> {
    /// The original constitutive state, traction, slip, storage and heat.
    /// Compatibility residuals refer to the full supplied body matrix.
    pub contact: TangentialStep<S>,
    /// Accepted outer updates; isotropic and zero responses need none.
    pub iterations: usize,
    /// Original-law evaluations, including backtracking trials.
    pub law_evaluations: usize,
    /// `Delta H + D - F dot (du_free - C F)` [J].
    pub body_work_residual_j: f64,
}

/// No failed candidate is published and the borrowed history stays unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MatrixTangentialError {
    /// Invalid controls, nonpassive matrix, or inconsistent original-law output.
    Invalid(&'static str),
    /// The original constitutive owner refused an input or candidate.
    Law(TriboError),
    /// Caller cancellation observed before or after a bounded law evaluation.
    Cancelled,
    /// The explicit solve budget did not satisfy actual matrix compatibility.
    Unresolved {
        /// Largest displacement or work-equivalent displacement residual [m].
        residual_m: f64,
        /// Requested displacement tolerance [m].
        tolerance_m: f64,
        /// Original-law calls spent before refusal.
        evaluations: usize,
    },
}

impl fmt::Display for MatrixTangentialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "matrix tangential compliance: {self:?}")
    }
}
impl std::error::Error for MatrixTangentialError {}
impl From<TriboError> for MatrixTangentialError {
    fn from(error: TriboError) -> Self {
        Self::Law(error)
    }
}

/// Couple one original two-axis constitutive law to the full body compliance.
///
/// `input.compliance_m_per_n` must be zero: the supplied matrix is the sole
/// physical body response. A fixed tangent frame, interval, load, temperature
/// and pressure retain their original constitutive meanings. Inner-law errors
/// are preserved, and cancellation is polled around every original-law call.
///
/// An exactly isotropic matrix delegates directly to the original scalar solve.
/// Otherwise a full-matrix zero-slip reaction is tried through the original law
/// first (the minimum-norm reaction for a rank-one matrix). A native refusal of
/// that trial remains a refusal, including a feasible static reaction exceeding
/// the caller's force budget; it cannot silently select sliding instead.
/// For a general matrix, the smallest positive eigenvalue selects the scalar
/// probe. Using the largest eigenvalue can incorrectly choose a static branch
/// when static and kinetic Coulomb coefficients differ. Rank-one null motion
/// bounds the probe using the explicit force budget. Neither choice changes C.
pub fn advance_with_matrix_compliance<L: TangentialLaw>(
    law: &L,
    state: &L::State,
    input: &TangentialStepInput,
    compliance_m_per_n: [[f64; 2]; 2],
    config: MatrixTangentialConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<MatrixTangentialStep<L::State>, MatrixTangentialError> {
    config.validate()?;
    input.solve.validate()?;
    if input.compliance_m_per_n != 0.0 || input.slip_increment_m.iter().any(|v| !v.is_finite()) {
        return Err(MatrixTangentialError::Invalid(
            "matrix owns the sole body compliance and requires finite released slip",
        ));
    }
    let matrix = Matrix::admit(compliance_m_per_n, config.matrix_relative_tolerance)?;
    let stored_before = number(law.stored_energy_j(state))?;
    if stored_before < 0.0 {
        return Err(MatrixTangentialError::Invalid(
            "negative original contact storage",
        ));
    }
    let probe = matrix.probe(input, config)?;
    let initial_force = matrix.zero_slip_reaction(input, config).unwrap_or([0.0; 2]);
    let initial_inverse = matrix.inverse_initial(probe);
    let mut problem = Problem {
        law,
        state,
        input,
        matrix,
        config,
        probe,
        stored_before,
        evaluations: 0,
        cancelled,
    };
    let mut current = problem.evaluate(initial_force)?;
    let mut inverse = initial_inverse;
    let mut iterations = 0;
    loop {
        if current.resolved(config) {
            problem.poll()?;
            return Ok(MatrixTangentialStep {
                contact: current.contact,
                iterations,
                law_evaluations: problem.evaluations,
                body_work_residual_j: current.body_work_residual_j,
            });
        }
        if problem.evaluations == config.maximum_evaluations {
            return Err(problem.unresolved(&current));
        }
        let mut next = None;
        // First use the secant inverse. If it does not reduce the residual,
        // restart from the bounded spectral inverse of I + C/c.
        for restart in 0..2 {
            if restart == 1 {
                inverse = initial_inverse;
            }
            let direction = mul(inverse, current.residual_n);
            let mut scale = 1.0;
            for _ in 0..16 {
                if problem.evaluations == config.maximum_evaluations {
                    break;
                }
                let trial_force = project_force(
                    [
                        current.force[0] - scale * direction[0],
                        current.force[1] - scale * direction[1],
                    ],
                    input.solve.maximum_force_n,
                );
                scale *= 0.5;
                if trial_force == current.force || trial_force.iter().any(|v| !v.is_finite()) {
                    continue;
                }
                match problem.evaluate(trial_force) {
                    Ok(trial)
                        if trial.resolved(config)
                            || trial.force_residual_n < current.force_residual_n =>
                    {
                        next = Some(trial);
                        break;
                    }
                    Ok(_)
                    | Err(MatrixTangentialError::Law(_))
                    | Err(MatrixTangentialError::Invalid(_)) => {}
                    Err(error) => return Err(error),
                }
            }
            if next.is_some() {
                break;
            }
        }
        let Some(candidate) = next else {
            return Err(problem.unresolved(&current));
        };
        inverse = secant_inverse(inverse, &current, &candidate).unwrap_or(initial_inverse);
        current = candidate;
        iterations += 1;
    }
}

struct Trial<S> {
    force: [f64; 2],
    residual_n: [f64; 2],
    force_residual_n: f64,
    displacement_residual_m: f64,
    body_work_residual_j: f64,
    contact: TangentialStep<S>,
}

impl<S> Trial<S> {
    fn resolved(&self, config: MatrixTangentialConfig) -> bool {
        self.displacement_residual_m <= config.displacement_tolerance_m
            && self.body_work_residual_j.abs() <= config.energy_tolerance_j
    }
}

struct Problem<'a, L: TangentialLaw, F: FnMut() -> bool> {
    law: &'a L,
    state: &'a L::State,
    input: &'a TangentialStepInput,
    matrix: Matrix,
    config: MatrixTangentialConfig,
    probe: f64,
    stored_before: f64,
    evaluations: usize,
    cancelled: F,
}

impl<L: TangentialLaw, F: FnMut() -> bool> Problem<'_, L, F> {
    fn poll(&mut self) -> Result<(), MatrixTangentialError> {
        if (self.cancelled)() {
            Err(MatrixTangentialError::Cancelled)
        } else {
            Ok(())
        }
    }

    fn evaluate(&mut self, force: [f64; 2]) -> Result<Trial<L::State>, MatrixTangentialError> {
        self.poll()?;
        let response = mul(self.matrix.original, force);
        let mut input = *self.input;
        input.compliance_m_per_n = self.probe;
        for i in 0..2 {
            input.slip_increment_m[i] =
                number(self.input.slip_increment_m[i] + (self.probe * force[i] - response[i]))?;
        }
        self.evaluations += 1;
        let mut contact = self.law.advance(self.state, &input)?;
        self.poll()?;
        let after = number(self.law.stored_energy_j(&contact.state))?;
        let change = number(after - self.stored_before)?;
        let heat = number(contact.dissipated_j)?;
        let force_norm = norm(contact.traction_n)?;
        if after < 0.0
            || heat < 0.0
            || force_norm > self.input.solve.maximum_force_n
            || contact.iterations > self.input.solve.max_iterations
            || number(change - contact.stored_change_j)?.abs() > self.config.energy_tolerance_j
            || number(contact.closure_residual_j)?.abs() > self.config.energy_tolerance_j
        {
            return Err(MatrixTangentialError::Invalid(
                "original constitutive storage/work/budget",
            ));
        }
        let native_work = number(dot(contact.traction_n, contact.slip_increment_m))?;
        if number(change + heat - native_work)?.abs() > self.config.energy_tolerance_j {
            return Err(MatrixTangentialError::Invalid(
                "original constitutive energy closure",
            ));
        }
        let mut residual_n = [0.0; 2];
        let mut actual = [0.0; 2];
        let actual_response = mul(self.matrix.original, contact.traction_n);
        for i in 0..2 {
            residual_n[i] = number(force[i] - contact.traction_n[i])?;
            actual[i] = number(self.input.slip_increment_m[i] - actual_response[i])?;
            contact.compatibility_residual_m[i] = number(contact.slip_increment_m[i] - actual[i])?;
        }
        Ok(Trial {
            force,
            residual_n,
            force_residual_n: norm(residual_n)?,
            displacement_residual_m: norm(contact.compatibility_residual_m)?,
            body_work_residual_j: number(change + heat - dot(contact.traction_n, actual))?,
            contact,
        })
    }

    fn unresolved(&self, trial: &Trial<L::State>) -> MatrixTangentialError {
        // Express either failed acceptance condition in the displacement units
        // of this diagnostic without hiding a stricter work-budget failure.
        let work_equivalent = self.config.displacement_tolerance_m
            * (trial.body_work_residual_j.abs() / self.config.energy_tolerance_j);
        MatrixTangentialError::Unresolved {
            residual_m: trial.displacement_residual_m.max(work_equivalent),
            tolerance_m: self.config.displacement_tolerance_m,
            evaluations: self.evaluations,
        }
    }
}

struct Matrix {
    original: [[f64; 2]; 2],
    normalized: [[f64; 2]; 2],
    scale: f64,
    minimum: f64,
    maximum: f64,
}

impl Matrix {
    fn admit(original: [[f64; 2]; 2], tolerance: f64) -> Result<Self, MatrixTangentialError> {
        let mut scale = 0.0_f64;
        for value in original.into_iter().flatten() {
            number(value)?;
            scale = scale.max(value.abs());
        }
        if scale == 0.0 {
            return Ok(Self {
                original,
                normalized: original,
                scale,
                minimum: 0.0,
                maximum: 0.0,
            });
        }
        let a = original[0][0] / scale;
        let d = original[1][1] / scale;
        let b = 0.5 * (original[0][1] / scale) + 0.5 * (original[1][0] / scale);
        if (original[0][1] / scale - original[1][0] / scale).abs() > tolerance {
            return Err(MatrixTangentialError::Invalid(
                "nonsymmetric body compliance",
            ));
        }
        let middle = 0.5 * a + 0.5 * d;
        let spread = det::hypot(0.5 * a - 0.5 * d, b);
        let maximum = middle + spread;
        let minimum = if maximum > 0.0 {
            (a * d - b * b) / maximum
        } else {
            middle - spread
        };
        if maximum <= 0.0 || minimum < -tolerance {
            return Err(MatrixTangentialError::Invalid(
                "nonpositive-semidefinite body compliance",
            ));
        }
        number(maximum * scale)?;
        Ok(Self {
            original,
            normalized: [[a, b], [b, d]],
            scale,
            minimum: minimum.max(0.0),
            maximum,
        })
    }

    fn probe(
        &self,
        input: &TangentialStepInput,
        config: MatrixTangentialConfig,
    ) -> Result<f64, MatrixTangentialError> {
        if self.scale == 0.0 {
            return Ok(0.0);
        }
        if self.original[0][1] == 0.0
            && self.original[1][0] == 0.0
            && self.original[0][0] == self.original[1][1]
        {
            return Ok(self.original[0][0]);
        }
        let mut probe = self.minimum * self.scale;
        if probe == 0.0 {
            probe = self.maximum * self.scale;
            let projected = mul(self.normalized, input.slip_increment_m);
            let null_motion = norm([
                input.slip_increment_m[0] - projected[0] / self.maximum,
                input.slip_increment_m[1] - projected[1] / self.maximum,
            ])?;
            if null_motion > config.displacement_tolerance_m {
                probe = probe.min((0.5 * null_motion) / input.solve.maximum_force_n);
            }
        }
        if !probe.is_finite() || probe <= 0.0 {
            return Err(MatrixTangentialError::Invalid(
                "unrepresentable scalar probe compliance",
            ));
        }
        Ok(probe)
    }

    // Try the full body's zero-slip reaction before a kinetic root. This uses
    // only mechanics, and acceptance still belongs to the original law. A layer
    // may return elastic motion instead; "Sticking" alone is never sufficient.
    fn zero_slip_reaction(
        &self,
        input: &TangentialStepInput,
        config: MatrixTangentialConfig,
    ) -> Option<[f64; 2]> {
        if self.scale == 0.0
            || (self.original[0][1] == 0.0
                && self.original[1][0] == 0.0
                && self.original[0][0] == self.original[1][1])
        {
            return None;
        }
        let free = input.slip_increment_m;
        let force = if self.minimum > 0.0 {
            let a = self.original[0][0] / self.scale;
            let b = self.original[0][1] / self.scale;
            let c = self.original[1][0] / self.scale;
            let d = self.original[1][1] / self.scale;
            let determinant = a * d - b * c;
            [
                ((d * free[0] - b * free[1]) / determinant) / self.scale,
                ((a * free[1] - c * free[0]) / determinant) / self.scale,
            ]
        } else {
            let projected = mul(self.normalized, free);
            [
                (projected[0] / self.maximum) / (self.maximum * self.scale),
                (projected[1] / self.maximum) / (self.maximum * self.scale),
            ]
        };
        if force.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let response = mul(self.original, force);
        let residual = det::hypot(free[0] - response[0], free[1] - response[1]);
        (residual <= config.displacement_tolerance_m).then_some(force)
    }

    fn inverse_initial(&self, probe: f64) -> [[f64; 2]; 2] {
        if self.scale == 0.0 {
            return [[1.0, 0.0], [0.0, 1.0]];
        }
        let c = probe / self.scale;
        let low = c / (c + self.minimum);
        let high = c / (c + self.maximum);
        let difference = self.maximum - self.minimum;
        if difference == 0.0 {
            return [[low, 0.0], [0.0, low]];
        }
        let factor = (high - low) / difference;
        [
            [
                low + factor * (self.normalized[0][0] - self.minimum),
                factor * self.normalized[0][1],
            ],
            [
                factor * self.normalized[1][0],
                low + factor * (self.normalized[1][1] - self.minimum),
            ],
        ]
    }
}

// The least-change inverse secant update uses vector norms/dot products, so
// rotating the fixed tangent axes does not select a preferred solver direction.
fn secant_inverse<S>(
    mut inverse: [[f64; 2]; 2],
    old: &Trial<S>,
    new: &Trial<S>,
) -> Option<[[f64; 2]; 2]> {
    let y = [
        new.residual_n[0] - old.residual_n[0],
        new.residual_n[1] - old.residual_n[1],
    ];
    let length = det::hypot(y[0], y[1]);
    if !length.is_finite() || length == 0.0 {
        return None;
    }
    let hy = mul(inverse, y);
    for i in 0..2 {
        let correction = (new.force[i] - old.force[i] - hy[i]) / length;
        for j in 0..2 {
            inverse[i][j] += correction * (y[j] / length);
        }
    }
    inverse
        .into_iter()
        .flatten()
        .all(f64::is_finite)
        .then_some(inverse)
}

fn project_force(mut force: [f64; 2], maximum: f64) -> [f64; 2] {
    let length = det::hypot(force[0], force[1]);
    if length.is_finite() && length > maximum {
        force[0] = (force[0] / length) * maximum;
        force[1] = (force[1] / length) * maximum;
    }
    force
}

fn mul(matrix: [[f64; 2]; 2], vector: [f64; 2]) -> [f64; 2] {
    [dot(matrix[0], vector), dot(matrix[1], vector)]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn norm(value: [f64; 2]) -> Result<f64, MatrixTangentialError> {
    number(det::hypot(value[0], value[1]))
}
fn number(value: f64) -> Result<f64, MatrixTangentialError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MatrixTangentialError::Invalid(
            "nonfinite matrix tangential arithmetic",
        ))
    }
}

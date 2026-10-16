//! Explicit algebraic admission, shared by reference and prepared solves.
//!
//! The arithmetic floor describes the evaluated residual, conditional on the
//! storage callbacks. It is not an error bound for a user-supplied constitutive
//! law, a discretization, or the resulting physical trajectory.

use crate::{PhsError, Storage};

/// Caller-declared absolute/relative residual and Newton-update limits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveBudget {
    /// Relative residual ceiling, multiplied by the recorded scale.
    pub relative: f64,
    /// Absolute ceiling in the residual's units; combined with relative by max.
    pub absolute: f64,
    /// Maximum number of Newton updates (initial residual is iteration zero).
    pub max_iterations: usize,
}

impl Default for SolveBudget {
    fn default() -> Self {
        Self {
            relative: crate::NEWTON_TOL,
            absolute: 0.0,
            max_iterations: crate::NEWTON_MAX,
        }
    }
}

impl SolveBudget {
    /// Validate the declaration, without changing an existing solver's policy.
    ///
    /// # Errors
    /// Negative/nonfinite tolerances or an update limit above the shared cap.
    pub fn validate(self) -> Result<Self, PhsError> {
        if !self.relative.is_finite()
            || self.relative < 0.0
            || !self.absolute.is_finite()
            || self.absolute < 0.0
            || self.max_iterations > crate::NEWTON_MAX
        {
            return Err(PhsError::Dimension {
                what: "invalid declared Newton solve budget",
            });
        }
        Ok(self)
    }

    /// Absolute algebraic ceiling at the given state/equation scale.
    #[must_use]
    pub fn threshold(self, scale: f64) -> f64 {
        self.absolute.max(self.relative * scale)
    }

    /// Record an explicit request; no rounding allowance is added here.
    #[must_use]
    pub fn certificate(
        self,
        site: NewtonSite,
        residual: f64,
        scale: f64,
        rounding_floor: f64,
        iterations: usize,
        stagnant_iterations: usize,
    ) -> SolveCertificate {
        let requested_budget = self.threshold(scale);
        SolveCertificate {
            site,
            residual,
            requested_budget,
            effective_budget: requested_budget,
            scale,
            rounding_floor,
            floor_raised: false,
            iterations,
            stagnant_iterations,
        }
    }
}

/// Physical equation whose algebraic residual was assessed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewtonSite {
    /// Descriptor differential/algebraic update.
    Descriptor,
    /// Allocating Gonzalez update.
    Step,
    /// Transactional prepared update (finite differences or analytic tangent).
    Prepared,
    /// Equal-effort Kirchhoff split.
    KirchhoffParallel,
    /// Modal aperture work closure [J].
    ModalWork,
    /// Modal aperture volume-flow closure [m^3/s].
    ModalFlow,
}

/// Actual equation, tolerance and arithmetic floor of one candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveCertificate {
    /// Assessed equation.
    pub site: NewtonSite,
    /// Absolute observed infinity-norm residual, in the equation's units.
    pub residual: f64,
    /// Absolute requested ceiling before any default-API roundoff allowance.
    pub requested_budget: f64,
    /// Actual admission ceiling; explicit declarations never raise this.
    pub effective_budget: f64,
    /// Scale used to convert the relative declaration to an absolute ceiling.
    pub scale: f64,
    /// Computed arithmetic floor, in the same units.
    pub rounding_floor: f64,
    /// True exactly when a default API raised the request to its computed floor.
    pub floor_raised: bool,
    /// Newton updates; modal balance gates record scalar residual evaluations.
    pub iterations: usize,
    /// Consecutive non-improving iterations at assessment.
    pub stagnant_iterations: usize,
}

impl SolveCertificate {
    /// Declare the default API's evaluated roundoff allowance visibly.
    #[must_use]
    pub fn with_rounding_allowance(mut self) -> Self {
        self.floor_raised = self.rounding_floor > self.requested_budget;
        self.effective_budget = self.requested_budget.max(self.rounding_floor);
        self
    }
}

/// Admit exactly at the declared cap; a next-up residual refuses.
///
/// # Errors
/// [`PhsError::BudgetBelowRoundingFloor`] if the requested ceiling cannot resolve
/// the evaluated equation. Otherwise [`PhsError::NewtonRefused`] for nonfinite
/// or negative fields, a nonpositive scale, or a residual above the ceiling.
pub fn admit_newton_iterate(c: SolveCertificate) -> Result<SolveCertificate, PhsError> {
    admit(c, false)
}

/// Admit a default-API certificate whose roundoff allowance is fully recorded.
///
/// # Errors
/// Incoherent or nonfinite certificates, or a residual above the effective cap.
/// Dropping the floor while retaining a raised cap is an incoherent certificate.
pub fn admit_default_newton_iterate(c: SolveCertificate) -> Result<SolveCertificate, PhsError> {
    admit(c, true)
}

fn admit(c: SolveCertificate, automatic_floor: bool) -> Result<SolveCertificate, PhsError> {
    if ![
        c.residual,
        c.requested_budget,
        c.effective_budget,
        c.scale,
        c.rounding_floor,
    ]
    .iter()
    .all(|v| v.is_finite())
        || c.residual < 0.0
        || c.requested_budget < 0.0
        || c.scale <= 0.0
        || c.rounding_floor < 0.0
    {
        return Err(PhsError::NewtonRefused { certificate: c });
    }
    if !automatic_floor && c.rounding_floor > c.requested_budget {
        return Err(PhsError::BudgetBelowRoundingFloor { certificate: c });
    }
    let expected = if automatic_floor {
        c.requested_budget.max(c.rounding_floor)
    } else {
        c.requested_budget
    };
    if c.effective_budget != expected
        || c.floor_raised != (automatic_floor && c.rounding_floor > c.requested_budget)
        || c.residual > c.effective_budget
    {
        return Err(PhsError::NewtonRefused { certificate: c });
    }
    Ok(c)
}

/// Standard gamma_n envelope for the supplied *evaluated* absolute operands.
/// No hardcoded noise allowance or observed residual enters this computation.
#[must_use]
pub fn arithmetic_rounding_floor(operations: usize, magnitude: f64) -> f64 {
    let ne = operations as f64 * f64::EPSILON;
    if !magnitude.is_finite() || magnitude < 0.0 || ne >= 1.0 {
        return f64::INFINITY;
    }
    ne / (1.0 - ne) * magnitude
}

/// Arithmetic floor for x1-x0-dt*((J-R)dg+Gu), or unscaled algebraic rows.
/// Gonzalez's energy-subtraction quotient contributes an additional propagated
/// envelope when the SAME cancellation guard used by the DG owner is active.
/// This does not estimate finite-difference truncation or callback model error.
#[allow(clippy::too_many_arguments)]
pub(crate) fn residual_rounding_floor(
    storage: &dyn Storage,
    x0: &[f64],
    x: &[f64],
    effort: &[f64],
    j: &[f64],
    r: &[f64],
    gu: &[f64],
    dt: f64,
    n_diff: usize,
    additional_flow: &[f64],
) -> f64 {
    let n = x.len();
    let state_scale = x0
        .iter()
        .chain(x)
        .fold(f64::MIN_POSITIVE, |a, v| a.max(v.abs()));
    let length2: f64 = x.iter().zip(x0).map(|(b, a)| (b - a) * (b - a)).sum();
    let guard = 1.0e-14 * state_scale;
    let energy_magnitude = storage.hamiltonian(x0).abs()
        + storage.hamiltonian(x).abs()
        + effort
            .iter()
            .zip(x.iter().zip(x0))
            .map(|(g, (b, a))| (g * (b - a)).abs())
            .sum::<f64>();
    let energy_floor =
        arithmetic_rounding_floor(n.saturating_mul(4).saturating_add(8), energy_magnitude);
    let mut floor = 0.0_f64;
    for row in 0..n {
        let mut flow_magnitude =
            gu[row].abs() + additional_flow.get(row).copied().unwrap_or(0.0).abs();
        let mut propagated = 0.0;
        for col in 0..n {
            let operator = j[row * n + col].abs() + r[row * n + col].abs();
            flow_magnitude += operator * effort[col].abs();
            if length2 > guard * guard && col < n_diff {
                propagated += operator * energy_floor * ((x[col] - x0[col]).abs() / length2);
            }
        }
        let factor = if row < n_diff { dt.abs() } else { 1.0 };
        let magnitude = factor * flow_magnitude
            + if row < n_diff {
                x[row].abs() + x0[row].abs()
            } else {
                0.0
            };
        let candidate = arithmetic_rounding_floor(n.saturating_mul(4).saturating_add(8), magnitude)
            + factor * propagated;
        if !candidate.is_finite() {
            return f64::INFINITY;
        }
        floor = floor.max(candidate);
    }
    floor
}

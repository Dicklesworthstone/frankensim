//! Endpoint-pressure Archard removal with the original contact equilibrium.
//!
//! The unknown is pressure at each cell with nonzero prescribed irreversible
//! slip. For each trial, the ORIGINAL wear quadrature evaluates
//! `dr_body = k(p, |ds|/dt, T) * p * |ds| / H_body`. The ORIGINAL normal owner
//! then equilibrates that recession from the accepted predecessor, retaining
//! its plastic history and its normal/parameter work accounts. The root requires
//! the quadrature pressure to equal that native endpoint pressure.
//!
//! This is a backward-Euler quadrature with mean speed and held original source
//! conditions. It is first order and still requires physical time/path refinement;
//! endpoint convergence alone does not bound finite-distance error or establish
//! uniqueness for an arbitrary pressure-dependent empirical coefficient.
//! There is no friction/slip closure, debris transport, material export or heat
//! allocation here. Constant approach, maintained force and a fixed loading
//! spring use their original native normal boundaries.

use std::cell::Cell;

use fs_exec::Cx;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};

use super::{
    SourcePlasticArchardWear, SourcePlasticContactError, SourcePlasticDriveError,
    SourcePlasticLoadError, SourcePlasticWearAdvance, SourcePlasticWearError,
    SourcePlasticWearState, checkpoint,
};
use crate::normal_patch::{
    FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt, FiniteGapPlasticLoadControls,
    FiniteGapPlasticLoadReceipt,
};

/// The accepted normal boundary is retained through every implicit trial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ImplicitWearBoundary {
    HeldApproach,
    MaintainedLoad {
        force_n: f64,
        controls: FiniteGapPlasticLoadControls,
    },
    FixedDrive {
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
    },
}

/// Original force/spring equilibrium receipt, without a duplicate work account.
#[derive(Clone, Debug, PartialEq)]
pub enum ImplicitWearBoundaryReceipt {
    HeldApproach,
    MaintainedLoad(FiniteGapPlasticLoadReceipt),
    FixedDrive(FiniteGapPlasticDriveReceipt),
}

/// Independent physical endpoint gates and bounded numerical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitWearControls {
    /// Maximum number of nonzero-slip grid cells in the nonlinear system.
    pub maximum_unknowns: usize,
    /// Includes residual, numerical Jacobian, line-search and final audit trials.
    /// Each trial also obeys the original normal owner's solve budget.
    pub maximum_evaluations: usize,
    pub maximum_newton_iterations: usize,
    /// Numerical perturbation in pressure/native-hardness-cap coordinates.
    pub difference_step: f64,
    /// Componentwise quadrature versus physical endpoint pressure gate [Pa].
    pub pressure_tolerance_pa: f64,
    /// Componentwise retained versus endpoint-law recession gate for each body [m].
    pub recession_tolerance_m: f64,
    /// Nonnegative scale using endpoint pressures and this interval's recession,
    /// never cumulative removed depth from earlier accepted history.
    pub relative_tolerance: f64,
    /// Numerical convergence acts on dimensionless pressure residuals. Physical
    /// pressure and recession gates above remain independent of this tolerance.
    pub newton: NewtonKrylovConfig,
}

/// One accepted implicit quadrature, with all original material inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitWearAdvance {
    /// Original quadrature/inventory and native normal continuation. All energy
    /// and plastic-work fields belong to this original receipt alone.
    pub wear: SourcePlasticWearAdvance,
    pub boundary: ImplicitWearBoundaryReceipt,
    /// Pressures actually used to query coefficients and integrate removal [Pa].
    /// Zero-slip entries carry the native endpoint pressure; they own no sample.
    pub quadrature_pressure_pa: Vec<f64>,
    /// Absent for a wholly zero-slip interval, which has no nonlinear unknowns.
    pub nonlinear: Option<NewtonReport>,
    pub evaluations: usize,
    /// Largest absolute quadrature/endpoint pressure difference on slipping cells.
    pub pressure_residual_pa: f64,
    /// Largest absolute per-body discrepancy from re-evaluating the original
    /// wear law at the ACTUAL returned endpoint pressure [m].
    pub recession_residual_m: f64,
    /// Maximum of the separate componentwise physical gate ratios, <= 1.
    pub maximum_normalized_residual: f64,
}

/// No refused solve publishes a partial wear or plastic history.
#[derive(Clone, Debug, PartialEq)]
pub enum ImplicitWearError {
    Invalid(&'static str),
    Step(SourcePlasticWearError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    EvaluationBudget,
    Residual {
        pressure_pa: f64,
        recession_m: f64,
        normalized: f64,
    },
}
impl From<SourcePlasticWearError> for ImplicitWearError {
    fn from(error: SourcePlasticWearError) -> Self {
        Self::Step(error)
    }
}
impl From<NewtonError> for ImplicitWearError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl core::fmt::Display for ImplicitWearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit original-source Archard wear: {self:?}")
    }
}
impl std::error::Error for ImplicitWearError {}

impl SourcePlasticArchardWear {
    /// Solve one implicit endpoint-pressure wear interval from accepted history.
    ///
    /// `initial_pressure_pa` is an explicit numerical guess on the complete
    /// original grid. The accepted predecessor's pressure is a useful starting
    /// point when its removal trial fits the existing physical/source envelope.
    /// Original material/interface temperatures and the normal boundary stay
    /// fixed. The interval supplies actual mean speed for speed-dependent k.
    ///
    /// Numerical coordinates are pressure divided by the native cap. Trial
    /// quadrature uses their projection onto [0,1], while the residual retains
    /// the UNPROJECTED coordinate: `R(x)=x-p_native(clamp(x,0,1))/H`. Roots
    /// therefore satisfy the physical pressure box. This semismooth projection
    /// prevents finite-difference probes from inventing tension or exceeding H;
    /// it does not extrapolate a coefficient outside its original source domain.
    /// Source failures are retained, including from numerical probes.
    ///
    /// Every normal trial starts from `before`, never from another trial. Final
    /// admission independently checks both pressure and each body's recession
    /// equation, so steep k(p) cannot hide a large wear error behind a small
    /// pressure defect. Return/adopt only the complete original native candidate.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_wear_implicit(
        &self,
        before: &SourcePlasticWearState,
        local_unsigned_slip_m: &[f64],
        duration_s: f64,
        boundary: ImplicitWearBoundary,
        initial_pressure_pa: &[f64],
        controls: ImplicitWearControls,
        cx: &Cx<'_>,
    ) -> Result<ImplicitWearAdvance, ImplicitWearError> {
        checkpoint(cx)?;
        self.admit(before)?;
        controls.newton.validate()?;
        let caps = &self.contact.model().config().hardness_pa;
        if local_unsigned_slip_m.len() != caps.len()
            || initial_pressure_pa.len() != caps.len()
            || local_unsigned_slip_m
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || initial_pressure_pa
                .iter()
                .zip(caps)
                .any(|(p, h)| !p.is_finite() || *p < 0.0 || p > h)
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || controls.maximum_evaluations == 0
            || controls.maximum_newton_iterations == 0
            || [
                controls.difference_step,
                controls.pressure_tolerance_pa,
                controls.recession_tolerance_m,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance < 0.0
        {
            return Err(ImplicitWearError::Invalid(
                "finite nonnegative actual slip, boxed pressure guess, positive interval and bounded controls required",
            ));
        }
        if let ImplicitWearBoundary::MaintainedLoad { controls, .. }
        | ImplicitWearBoundary::FixedDrive { controls, .. } = boundary
        {
            self.contact
                .admit_load_controls(controls)
                .map_err(SourcePlasticWearError::from)?;
        }
        let mut indices = Vec::new();
        let mut seed = Vec::new();
        for (i, slip) in local_unsigned_slip_m.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            if *slip > 0.0 {
                if indices.len() == controls.maximum_unknowns {
                    return Err(ImplicitWearError::Invalid(
                        "implicit sliding-cell unknown budget",
                    ));
                }
                indices.push(i);
                seed.push(initial_pressure_pa[i] / caps[i]);
            }
        }
        let evaluations = Cell::new(0usize);
        let evaluate_pressure = |pressure: &[f64]| {
            checkpoint(cx)?;
            if evaluations.get() >= controls.maximum_evaluations {
                return Err(ImplicitWearError::EvaluationBudget);
            }
            evaluations.set(evaluations.get() + 1);
            self.advance_wear_at_pressure_with_recession(
                before,
                pressure,
                local_unsigned_slip_m,
                Some(duration_s),
                cx,
                |increment| match boundary {
                    ImplicitWearBoundary::HeldApproach => {
                        let normal = self
                            .contact
                            .model()
                            .recede(&before.contact.contact, increment, cx)
                            .map_err(SourcePlasticContactError::from)?;
                        Ok((normal, ImplicitWearBoundaryReceipt::HeldApproach))
                    }
                    ImplicitWearBoundary::MaintainedLoad { force_n, controls } => {
                        let normal = self
                            .contact
                            .model()
                            .recede_at_load(
                                &before.contact.contact,
                                increment,
                                force_n,
                                controls,
                                cx,
                            )
                            .map_err(SourcePlasticLoadError::from)?;
                        Ok((
                            normal.continuation,
                            ImplicitWearBoundaryReceipt::MaintainedLoad(normal.load),
                        ))
                    }
                    ImplicitWearBoundary::FixedDrive { drive, controls } => {
                        let normal = self
                            .contact
                            .model()
                            .recede_with_drive(
                                &before.contact.contact,
                                increment,
                                drive,
                                controls,
                                cx,
                            )
                            .map_err(SourcePlasticDriveError::from)?;
                        Ok((
                            normal.continuation,
                            ImplicitWearBoundaryReceipt::FixedDrive(normal.drive),
                        ))
                    }
                },
            )
            .map_err(ImplicitWearError::from)
        };
        let evaluate = |x: &[f64]| {
            if x.len() != indices.len() || x.iter().any(|value| !value.is_finite()) {
                return Err(ImplicitWearError::Invalid(
                    "complete finite normalized pressure coordinates",
                ));
            }
            // Zero-slip entries are immaterial to quadrature and own no source
            // sample. Their physical pressures come only from the native solve.
            let mut pressure = vec![0.0; caps.len()];
            for (slot, &i) in indices.iter().enumerate() {
                if slot % 64 == 0 {
                    checkpoint(cx)?;
                }
                pressure[i] = x[slot].clamp(0.0, 1.0) * caps[i];
            }
            let (wear, normal) = evaluate_pressure(&pressure)?;
            let endpoint = wear.candidate.contact_state().contact_state().pressure_pa();
            let residual: Vec<_> = indices
                .iter()
                .enumerate()
                .map(|(slot, &i)| x[slot] - endpoint[i] / caps[i])
                .collect();
            Ok((residual, pressure, wear, normal))
        };
        let (solution, nonlinear) = if indices.is_empty() {
            (Vec::new(), None)
        } else {
            let problem = FallibleFiniteDifference::new(
                indices.len(),
                controls.difference_step,
                |x: &[f64]| evaluate(x).map(|(residual, _, _, _)| residual),
            )
            .map_err(numerical)?;
            let initialized = NewtonKrylovState::new(&problem, seed, controls.newton);
            if let Some(error) = problem.take_failure() {
                return Err(numerical(error));
            }
            let mut solver = initialized?;
            let result =
                solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                    cx.checkpoint().is_err() || problem.has_failed()
                });
            if let Some(error) = problem.take_failure() {
                return Err(numerical(error));
            }
            checkpoint(cx)?;
            if result.progress != SolverRunProgress::Complete {
                return Err(ImplicitWearError::Interrupted(result.progress));
            }
            if !result.report.converged {
                return Err(ImplicitWearError::NotConverged(result.report));
            }
            (solver.x, Some(result.report))
        };
        let (_, mut pressure, wear, boundary) = evaluate(&solution)?;
        let actual = wear.candidate.contact_state().contact_state().pressure_pa();
        // An independent re-evaluation at physical endpoint pressure tests the
        // actual constitutive recession, including pressure-dependent k. Its
        // candidate and work are discarded; they never advance the predecessor.
        let endpoint_wear = if indices.is_empty() {
            None
        } else {
            Some(evaluate_pressure(actual)?.0)
        };
        let mut pressure_residual_pa: f64 = 0.0;
        let mut recession_residual_m: f64 = 0.0;
        let mut normalized: f64 = 0.0;
        for (i, slip) in local_unsigned_slip_m.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            if *slip == 0.0 {
                pressure[i] = actual[i];
                continue;
            }
            let error = (pressure[i] - actual[i]).abs();
            pressure_residual_pa = pressure_residual_pa.max(error);
            normalized = normalized.max(ratio(
                error,
                controls.pressure_tolerance_pa,
                controls.relative_tolerance,
                pressure[i].abs().max(actual[i].abs()),
            )?);
            for side in 0..2 {
                let retained = wear.recession_increment_m[side][i];
                let endpoint = endpoint_wear
                    .as_ref()
                    .expect("nonzero slip has endpoint audit")
                    .recession_increment_m[side][i];
                let error = (retained - endpoint).abs();
                recession_residual_m = recession_residual_m.max(error);
                normalized = normalized.max(ratio(
                    error,
                    controls.recession_tolerance_m,
                    controls.relative_tolerance,
                    retained.abs().max(endpoint.abs()),
                )?);
            }
        }
        if normalized > 1.0 {
            return Err(ImplicitWearError::Residual {
                pressure_pa: pressure_residual_pa,
                recession_m: recession_residual_m,
                normalized,
            });
        }
        checkpoint(cx)?;
        Ok(ImplicitWearAdvance {
            wear,
            boundary,
            quadrature_pressure_pa: pressure,
            nonlinear,
            evaluations: evaluations.get(),
            pressure_residual_pa,
            recession_residual_m,
            maximum_normalized_residual: normalized,
        })
    }
}

fn ratio(error: f64, absolute: f64, relative: f64, scale: f64) -> Result<f64, ImplicitWearError> {
    let denominator = absolute + relative * scale;
    let result = error / denominator;
    if !error.is_finite() || !denominator.is_finite() || !result.is_finite() {
        Err(ImplicitWearError::Invalid(
            "representable physical residual scale",
        ))
    } else {
        Ok(result)
    }
}

fn numerical(error: FiniteDifferenceError<ImplicitWearError>) -> ImplicitWearError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => ImplicitWearError::Invalid(reason),
    }
}

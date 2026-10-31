//! Material-bearing building design under empirical peak-drift tail risk.
//!
//! Every supplied catalog candidate is evaluated against the same original
//! motions by the native distributed-plasticity frame. Selection minimizes the
//! native reference mass over the feasible finite catalog; neither response
//! monotonicity in section size nor a continuous-design optimum is assumed.
//! Empirical feasibility on this motion ensemble is not a population-risk
//! certificate or measured-building qualification.

use crate::building::BuildingFrame;
use fs_exec::Cx;
use fs_robust::{EmpiricalCvarReport, RobustError, empirical_cvar};
use fs_scenario::ScenarioError;
use fs_scenario::ensemble::{RealizationBudget, SpectrumModel, StochasticEnsemble};
use fs_scenario::signal::{Interp, TimeSignal};
use fs_solid::SolidError;
use std::borrow::Cow;

/// Original typed forcing retained by a building risk/design study.
#[derive(Debug, Clone, Copy)]
pub enum BuildingMotionSource<'a> {
    /// Seeded Kanai-Tajimi motions, synthesized once per member.
    Ensemble(&'a StochasticEnsemble),
    /// Equally weighted recorded acceleration tables, each with its own grid.
    ///
    /// Only linear, uniform SI acceleration tables beginning at zero are
    /// admitted. The first value is the native building's initial forcing;
    /// n samples yield n - 1 advancing steps. Tables carry their original
    /// units/interpolation but do not themselves certify recording provenance.
    Records(&'a [TimeSignal]),
}

/// Explicit risk target and deterministic bounds for the finite catalog study.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BuildingCvarControls {
    /// Tail confidence beta, strictly between zero and one.
    pub beta: f64,
    /// Maximum empirical CVaR of peak inter-story drift ratio.
    pub drift_limit: f64,
    /// Maximum admitted complete building candidates.
    pub max_candidates: usize,
    /// Maximum admitted generated motions or recorded tables.
    pub max_members: u32,
    /// Per-motion sample limit; max_work additionally bounds ensemble synthesis.
    pub realization_budget: RealizationBudget,
}

/// One candidate's actual mass and complete empirical drift-loss sample.
#[derive(Debug, Clone, PartialEq)]
pub struct BuildingCvarEvaluation {
    /// Full reference physical mass, including floors and supported members [kg].
    pub mass_kg: f64,
    /// Peak inter-story drift ratios in original ensemble-member order.
    pub peak_drifts: Vec<f64>,
    /// Canonical fractional-tail CVaR and its order-statistic metadata.
    pub tail: EmpiricalCvarReport,
}

/// Minimum-mass feasible candidate in the supplied finite catalog.
#[derive(Debug)]
pub struct BuildingCvarDesign<'a> {
    /// Selected candidate's original catalog index; equal masses prefer the first.
    pub selected_index: usize,
    /// Original admitted candidate, retaining material bindings and initial state.
    pub selected: &'a BuildingFrame,
    /// Original catalog in evaluation order.
    pub catalog: &'a [BuildingFrame],
    /// Original complete generation recipe or borrowed typed recording suite.
    pub motions: BuildingMotionSource<'a>,
    /// Exact controls used for this study.
    pub controls: BuildingCvarControls,
    /// Complete results aligned with the original catalog.
    pub evaluations: Vec<BuildingCvarEvaluation>,
}

/// A refused study never mutates candidates or publishes a partial optimum.
#[derive(Debug)]
pub enum BuildingCvarError {
    /// Invalid target, empty study, overflow, or exceeded candidate/member cap.
    InvalidInput(&'static str),
    /// Candidate does not bind its section laws to explicit constituent densities.
    MissingMaterialSections {
        /// Original catalog index.
        candidate: usize,
    },
    /// Candidate changes the common building envelope, floor mass, gravity or solve setup.
    IncompatibleCandidate {
        /// Original catalog index.
        candidate: usize,
    },
    /// The owning scenario generator refused a member or its work budget.
    Realization {
        /// Original member index.
        member: u32,
        /// Original scenario diagnostic.
        source: ScenarioError,
    },
    /// A recorded table or its sampling contract was refused.
    Record {
        /// Original recording index.
        member: u32,
        /// Original signal diagnostic or explicit structural sampling refusal.
        source: ScenarioError,
    },
    /// The owning nonlinear building solver refused a complete history.
    History {
        /// Original catalog index.
        candidate: usize,
        /// Original ensemble-member index.
        member: u32,
        /// Original structural diagnostic.
        source: SolidError,
    },
    /// A completed history returned an invalid peak-drift loss.
    NonFiniteLoss {
        /// Original catalog index.
        candidate: usize,
        /// Original ensemble-member index.
        member: u32,
        /// Rejected drift ratio.
        value: f64,
    },
    /// Canonical empirical tail-risk calculation refused.
    Robust(RobustError),
    /// Cancellation observed at a bounded study-work boundary.
    Cancelled,
    /// Every successfully evaluated candidate violates the declared risk target.
    NoFeasibleCandidate {
        /// Declared empirical CVaR limit.
        limit: f64,
        /// Smallest empirical CVaR observed in the catalog.
        best_cvar: f64,
    },
    /// Loss or result storage could not be reserved.
    AllocationFailed,
}

impl core::fmt::Display for BuildingCvarError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(what) => write!(f, "building CVaR study: {what}"),
            Self::MissingMaterialSections { candidate } => {
                write!(f, "candidate {candidate} requires explicit material sections")
            }
            Self::IncompatibleCandidate { candidate } => {
                write!(f, "candidate {candidate} changes the common building problem")
            }
            Self::Realization { member, source } => {
                write!(f, "ground-motion member {member}: {source}")
            }
            Self::Record { member, source } => {
                write!(f, "ground-acceleration recording {member}: {source}")
            }
            Self::History { candidate, member, source } => {
                write!(f, "building candidate {candidate}, motion {member}: {source}")
            }
            Self::NonFiniteLoss { candidate, member, value } => {
                write!(f, "invalid drift {value} for candidate {candidate}, motion {member}")
            }
            Self::Robust(source) => write!(f, "building empirical CVaR: {source:?}"),
            Self::Cancelled => write!(f, "building CVaR study cancelled"),
            Self::NoFeasibleCandidate { limit, best_cvar } => {
                write!(f, "no catalog candidate meets CVaR limit {limit}; smallest is {best_cvar}")
            }
            Self::AllocationFailed => write!(f, "building CVaR result allocation failed"),
        }
    }
}

impl std::error::Error for BuildingCvarError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Realization { source, .. } | Self::Record { source, .. } => Some(source),
            Self::History { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl From<RobustError> for BuildingCvarError {
    fn from(source: RobustError) -> Self {
        Self::Robust(source)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), BuildingCvarError> {
    cx.checkpoint().map_err(|_| BuildingCvarError::Cancelled)
}

fn same_problem(a: &BuildingFrame, b: &BuildingFrame) -> bool {
    let (a, b) = (a.spec(), b.spec());
    a.stories == b.stories
        && a.bays == b.bays
        && a.story_height == b.story_height
        && a.bay_width == b.bay_width
        && a.floor_mass == b.floor_mass
        && a.gravity_m_s2 == b.gravity_m_s2
        && a.zeta == b.zeta
        && a.sections_per_member == b.sections_per_member
}

fn record_error(member: u32, what: &str) -> BuildingCvarError {
    BuildingCvarError::Record {
        member,
        source: ScenarioError::Evaluate { what: what.into() },
    }
}

fn admit_record<'a>(
    signal: &'a TimeSignal,
    member: u32,
    max_samples: usize,
    cx: &Cx<'_>,
) -> Result<(&'a [f64], f64), BuildingCvarError> {
    checkpoint(cx)?;
    let TimeSignal::Table { times, values, dims, interp } = signal else {
        return Err(record_error(member, "a tabulated acceleration recording is required"));
    };
    // Apply cheap shape/work caps before the native signal owner scans payloads.
    if times.len() < 2 || times.len() != values.len() || times.len() > max_samples {
        return Err(record_error(member, "record needs matching sample/time counts in 2..=max_samples"));
    }
    if dims.0 != [1, 0, -2, 0, 0, 0] || *interp != Interp::Linear {
        return Err(record_error(member, "record must carry SI acceleration with linear interpolation"));
    }
    let mut violations = Vec::new();
    signal.check("building ground-acceleration recording", &mut violations);
    if let Some(first) = violations.first() {
        return Err(record_error(member, &first.what));
    }
    checkpoint(cx)?;
    if times[0] != 0.0 {
        return Err(record_error(member, "record must include its initial forcing at time zero"));
    }
    let dt = times[1];
    if !dt.is_finite() || dt <= 0.0 {
        return Err(record_error(member, "record timestep must be finite and positive"));
    }
    for (index, &time) in times.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        #[allow(clippy::cast_precision_loss)]
        let expected = index as f64 * dt;
        // Admit only roundoff around the complete nominal grid, never an
        // accumulated sequence of locally similar but drifting intervals.
        let tolerance = 64.0 * f64::EPSILON * expected.abs().max(dt);
        if !expected.is_finite() || (time - expected).abs() > tolerance {
            return Err(record_error(member, "record is not uniformly sampled on its declared grid"));
        }
    }
    Ok((values, dt))
}

/// Evaluate generated motions through the shared material-building catalog path.
///
/// This is the ensemble convenience entry point for
/// [building_cvar_from_motions]; the returned source retains the original recipe.
pub fn building_cvar_mass_min<'a>(
    ensemble: &'a StochasticEnsemble,
    catalog: &'a [BuildingFrame],
    controls: BuildingCvarControls,
    cx: &Cx<'_>,
) -> Result<BuildingCvarDesign<'a>, BuildingCvarError> {
    building_cvar_from_motions(BuildingMotionSource::Ensemble(ensemble), catalog, controls, cx)
}

/// Select minimum reference mass subject to empirical peak-drift CVaR.
///
/// Candidates must be admitted through BuildingFrame::with_sections. Section
/// dimensions, column area scale and constituent materials may vary; building
/// envelope, floor masses, gravity, damping ratio and section quadrature must agree.
/// Compare candidates using the same intended inertia formulation and initial
/// loading/history assumptions.
///
/// Generated motions are synthesized once under the native work budget.
/// Recorded tables are all admitted before any history runs, then their
/// original value slices are used without resampling or dropping the initial
/// sample. Each record has its own timestep and contributes one equally
/// weighted peak-drift loss, regardless of duration. Uniformity admits only
/// 64 machine epsilons times max(abs(i * dt), dt) around the nominal grid.
/// Each motion is applied to every candidate on an independent clone, so
/// the catalog's source bindings, gravity preload and constitutive states are
/// unchanged even on failure. Every candidate is evaluated: resonance and
/// material-dependent inertia need not give monotone risk versus section size.
///
/// Working storage is one realized motion, one trial frame/history, and one
/// scalar loss per candidate/member. At most max_candidates * max_members
/// complete histories are admitted. Cancellation is checked before/after each
/// realization and each history; latency can include one complete bounded
/// native history, not merely one Newton correction.
///
/// The winning catalog entry is already evaluated in its actual geometry and
/// materials. There is no surrogate-scale snap or unchecked final replacement.
/// Equal feasible masses deterministically select the first catalog entry.
///
/// # Errors
/// Refuses malformed inputs, missing material bindings, incompatible catalogs,
/// exceeded budgets, failed realization/history, nonfinite loss, cancellation,
/// or a catalog with no empirically feasible design. A failed candidate is not
/// silently discarded from the minimum-mass search.
pub fn building_cvar_from_motions<'a>(
    motions: BuildingMotionSource<'a>,
    catalog: &'a [BuildingFrame],
    controls: BuildingCvarControls,
    cx: &Cx<'_>,
) -> Result<BuildingCvarDesign<'a>, BuildingCvarError> {
    checkpoint(cx)?;
    let member_count = match motions {
        BuildingMotionSource::Ensemble(ensemble) => {
            if !matches!(ensemble.model, SpectrumModel::KanaiTajimi { .. }) {
                return Err(BuildingCvarError::InvalidInput(
                    "ensemble must describe Kanai-Tajimi ground acceleration",
                ));
            }
            ensemble.members
        }
        BuildingMotionSource::Records(records) => u32::try_from(records.len()).map_err(|_| {
            BuildingCvarError::InvalidInput("recording count is not representable")
        })?,
    };
    if catalog.is_empty() || member_count == 0 {
        return Err(BuildingCvarError::InvalidInput(
            "nonempty material catalog and motion suite are required",
        ));
    }
    if !controls.drift_limit.is_finite() || controls.drift_limit < 0.0 {
        return Err(BuildingCvarError::InvalidInput(
            "drift CVaR limit must be finite and nonnegative",
        ));
    }
    // Use the canonical risk owner's admission before expensive histories.
    empirical_cvar(&[0.0], controls.beta)?;
    if catalog.len() > controls.max_candidates || member_count > controls.max_members {
        return Err(BuildingCvarError::InvalidInput(
            "candidate or motion count exceeds its explicit study limit",
        ));
    }
    let members = usize::try_from(member_count).map_err(|_| {
        BuildingCvarError::InvalidInput("ensemble member count is not representable")
    })?;
    catalog.len().checked_mul(members).ok_or(
        BuildingCvarError::InvalidInput("candidate-by-member history count overflows"),
    )?;
    for (candidate, frame) in catalog.iter().enumerate() {
        checkpoint(cx)?;
        if frame.material_sections().is_none() {
            return Err(BuildingCvarError::MissingMaterialSections { candidate });
        }
        if !same_problem(&catalog[0], frame) {
            return Err(BuildingCvarError::IncompatibleCandidate { candidate });
        }
        if !frame.reference_mass_kg().is_finite() || frame.reference_mass_kg() <= 0.0 {
            return Err(BuildingCvarError::InvalidInput(
                "candidate reference physical mass must be finite and positive",
            ));
        }
    }

    let mut recorded = Vec::new();
    if let BuildingMotionSource::Records(records) = motions {
        recorded.try_reserve_exact(members)
            .map_err(|_| BuildingCvarError::AllocationFailed)?;
        let mut total_samples = 0usize;
        for (index, signal) in records.iter().enumerate() {
            let member = u32::try_from(index).map_err(|_| {
                BuildingCvarError::InvalidInput("recording index is not representable")
            })?;
            let record = admit_record(signal, member, controls.realization_budget.max_samples, cx)?;
            total_samples = total_samples.checked_add(record.0.len()).ok_or(
                BuildingCvarError::InvalidInput("recorded suite sample count overflows"),
            )?;
            recorded.push(record);
        }
        total_samples.checked_mul(catalog.len()).ok_or(
            BuildingCvarError::InvalidInput("catalog-by-record sample work overflows"),
        )?;
    }

    let mut losses = Vec::new();
    losses.try_reserve_exact(catalog.len())
        .map_err(|_| BuildingCvarError::AllocationFailed)?;
    for _ in catalog {
        checkpoint(cx)?;
        let mut sample = Vec::new();
        sample.try_reserve_exact(members)
            .map_err(|_| BuildingCvarError::AllocationFailed)?;
        losses.push(sample);
    }
    for member in 0..member_count {
        checkpoint(cx)?;
        let (values, dt): (Cow<'_, [f64]>, f64) = match motions {
            BuildingMotionSource::Ensemble(ensemble) => {
                let motion = ensemble.realize_with_budget(member, controls.realization_budget)
                    .map_err(|source| BuildingCvarError::Realization { member, source })?;
                (Cow::Owned(motion.values), ensemble.dt.value)
            }
            BuildingMotionSource::Records(_) => {
                let index = usize::try_from(member).map_err(|_| {
                    BuildingCvarError::InvalidInput("recording index is not representable")
                })?;
                (Cow::Borrowed(recorded[index].0), recorded[index].1)
            }
        };
        checkpoint(cx)?;
        for (candidate, frame) in catalog.iter().enumerate() {
            checkpoint(cx)?;
            let mut trial = frame.clone();
            let response = trial.run(values.as_ref(), dt)
                .map_err(|source| BuildingCvarError::History { candidate, member, source })?;
            checkpoint(cx)?;
            let value = response.peak_drift_ratio;
            if !value.is_finite() || value < 0.0 {
                return Err(BuildingCvarError::NonFiniteLoss { candidate, member, value });
            }
            losses[candidate].push(value);
        }
    }

    let mut evaluations: Vec<BuildingCvarEvaluation> = Vec::new();
    evaluations.try_reserve_exact(catalog.len())
        .map_err(|_| BuildingCvarError::AllocationFailed)?;
    let mut selected_index: Option<usize> = None;
    let mut best_cvar = f64::INFINITY;
    for (candidate, peak_drifts) in losses.into_iter().enumerate() {
        checkpoint(cx)?;
        let tail = empirical_cvar(&peak_drifts, controls.beta)?;
        let mass_kg = catalog[candidate].reference_mass_kg();
        best_cvar = best_cvar.min(tail.cvar());
        if tail.cvar() <= controls.drift_limit
            && selected_index.is_none_or(|previous| mass_kg < evaluations[previous].mass_kg)
        {
            selected_index = Some(candidate);
        }
        evaluations.push(BuildingCvarEvaluation { mass_kg, peak_drifts, tail });
    }
    checkpoint(cx)?;
    let selected_index = selected_index.ok_or(BuildingCvarError::NoFeasibleCandidate {
        limit: controls.drift_limit,
        best_cvar,
    })?;
    Ok(BuildingCvarDesign {
        selected_index,
        selected: &catalog[selected_index],
        catalog,
        motions,
        controls,
        evaluations,
    })
}

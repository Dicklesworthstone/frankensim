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
use fs_solid::SolidError;

/// Explicit risk target and deterministic bounds for the finite catalog study.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BuildingCvarControls {
    /// Tail confidence beta, strictly between zero and one.
    pub beta: f64,
    /// Maximum empirical CVaR of peak inter-story drift ratio.
    pub drift_limit: f64,
    /// Maximum admitted complete building candidates.
    pub max_candidates: usize,
    /// Maximum admitted ground-motion ensemble members.
    pub max_members: u32,
    /// Per-motion sample and spectral-synthesis work limits.
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
    /// Original complete motion-generation specification.
    pub ensemble: &'a StochasticEnsemble,
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
    /// Candidate changes the common building envelope, floor mass or solve setup.
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
            Self::Realization { source, .. } => Some(source),
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
        && a.zeta == b.zeta
        && a.sections_per_member == b.sections_per_member
}

/// Select minimum reference mass subject to empirical peak-drift CVaR.
///
/// Candidates must be admitted through BuildingFrame::with_sections. Section
/// dimensions, column area scale and constituent materials may vary; building
/// envelope, floor masses, damping ratio and section quadrature must agree.
/// Compare candidates using the same intended inertia formulation and initial
/// loading/history assumptions.
///
/// Each original motion is realized once under the supplied native budget and
/// applied to every candidate. Each history runs on an independent clone, so
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
pub fn building_cvar_mass_min<'a>(
    ensemble: &'a StochasticEnsemble,
    catalog: &'a [BuildingFrame],
    controls: BuildingCvarControls,
    cx: &Cx<'_>,
) -> Result<BuildingCvarDesign<'a>, BuildingCvarError> {
    checkpoint(cx)?;
    if catalog.is_empty()
        || ensemble.members == 0
        || !matches!(ensemble.model, SpectrumModel::KanaiTajimi { .. })
    {
        return Err(BuildingCvarError::InvalidInput(
            "nonempty material catalog and Kanai-Tajimi motion ensemble are required",
        ));
    }
    if !controls.drift_limit.is_finite() || controls.drift_limit < 0.0 {
        return Err(BuildingCvarError::InvalidInput(
            "drift CVaR limit must be finite and nonnegative",
        ));
    }
    // Use the canonical risk owner's admission before expensive histories.
    empirical_cvar(&[0.0], controls.beta)?;
    if catalog.len() > controls.max_candidates || ensemble.members > controls.max_members {
        return Err(BuildingCvarError::InvalidInput(
            "candidate or motion count exceeds its explicit study limit",
        ));
    }
    let members = usize::try_from(ensemble.members).map_err(|_| {
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
    for member in 0..ensemble.members {
        checkpoint(cx)?;
        let motion = ensemble.realize_with_budget(member, controls.realization_budget)
            .map_err(|source| BuildingCvarError::Realization { member, source })?;
        checkpoint(cx)?;
        for (candidate, frame) in catalog.iter().enumerate() {
            checkpoint(cx)?;
            let mut trial = frame.clone();
            let response = trial.run(&motion.values, ensemble.dt.value)
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
        ensemble,
        controls,
        evaluations,
    })
}

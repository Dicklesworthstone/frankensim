//! The `Scenario` root value: frames + BCs + load cases + combinations +
//! ensembles + environment + contact laws, with whole-scenario validity
//! checking. Violations are STRUCTURED FIXES (code, what, fix) — the
//! agent-facing refusal format — never free-text panics. Environmental
//! fields are explicit constructor arguments: nothing is defaulted
//! silently.

use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
use crate::ensemble::StochasticEnsemble;
use crate::frame::{FrameMotion, FrameTree};
use fs_exec::Cx;
use fs_qty::{Dims, QtyAny};
use std::fmt;

const ACCEL_DIMS: Dims = Dims([1, 0, -2, 0, 0, 0]);
const TEMP_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
/// Net-flux tolerance relative to the gross flux magnitude.
const FLUX_REL_TOL: f64 = 1e-9;
const FIXED_FINDING_SLOTS: usize = 12;
const FINDINGS_PER_FRAME: usize = 13;
const FINDINGS_PER_BC: usize = 8;
const FINDINGS_PER_CASE: usize = 3;
const FINDINGS_PER_COMBINATION: usize = 2;
const FINDINGS_PER_COMBINATION_TERM: usize = 4;
const FINDINGS_PER_ENSEMBLE: usize = 16;
const FINDINGS_PER_CONTACT: usize = 4;

/// Explicit deterministic limits for whole-scenario semantic validation.
///
/// Collection fields cap public `Vec` authority before validation allocates
/// indexes. `max_signal_scalars` counts table times/values and Chebyshev
/// coefficients. `max_flux_checkpoints` counts raw set-local checkpoints
/// before deterministic deduplication. `max_work` covers record visits,
/// ordered-index comparisons, signal scans, checkpoint sorting, and flux
/// evaluations in machine-independent logical units.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValidationBudget {
    /// Maximum non-world frames.
    pub max_frames: usize,
    /// Maximum always-active boundary conditions.
    pub max_base_bcs: usize,
    /// Maximum load cases.
    pub max_cases: usize,
    /// Maximum boundary conditions across all load cases.
    pub max_case_bcs: usize,
    /// Maximum load combinations.
    pub max_combinations: usize,
    /// Maximum terms across all combinations.
    pub max_combination_terms: usize,
    /// Maximum stochastic ensembles.
    pub max_ensembles: usize,
    /// Maximum contact-law rows.
    pub max_contacts: usize,
    /// Maximum dynamically sized signal scalars visited by validation.
    pub max_signal_scalars: usize,
    /// Maximum raw net-flux validation checkpoints across effective sets.
    pub max_flux_checkpoints: usize,
    /// Maximum bytes across exact string identities and references.
    pub max_identity_bytes: usize,
    /// Maximum preflighted slots in the private finding buffer.
    pub max_findings: usize,
    /// Maximum deterministic logical work units.
    pub max_work: u128,
}

/// Conservative default for [`Scenario::validate`].
pub const DEFAULT_VALIDATION_BUDGET: ValidationBudget = ValidationBudget {
    max_frames: 131_072,
    max_base_bcs: 65_536,
    max_cases: 16_384,
    max_case_bcs: 262_144,
    max_combinations: 65_536,
    max_combination_terms: 262_144,
    max_ensembles: 65_536,
    max_contacts: 65_536,
    max_signal_scalars: 1_048_576,
    max_flux_checkpoints: 1_048_576,
    max_identity_bytes: 16_777_216,
    max_findings: 8_388_608,
    max_work: 134_217_728,
};

impl Default for ValidationBudget {
    fn default() -> Self {
        DEFAULT_VALIDATION_BUDGET
    }
}

/// Preflighted semantic-validation shape retained for ledger receipts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValidationPlan {
    /// Non-world frames.
    pub frames: usize,
    /// Always-active boundary conditions.
    pub base_bcs: usize,
    /// Load cases.
    pub cases: usize,
    /// Boundary conditions across all cases.
    pub case_bcs: usize,
    /// Load combinations.
    pub combinations: usize,
    /// Terms across all combinations.
    pub combination_terms: usize,
    /// Stochastic ensembles.
    pub ensembles: usize,
    /// Contact-law rows.
    pub contacts: usize,
    /// Dynamically sized signal scalars.
    pub signal_scalars: usize,
    /// Raw net-flux checkpoints across effective sets.
    pub flux_checkpoints: usize,
    /// Exact identity/reference bytes.
    pub identity_bytes: usize,
    /// Proved worst-case slots for validation findings.
    pub finding_capacity: usize,
    /// Deterministic logical work units.
    pub planned_work: u128,
}

/// Resource-admission or cancellation refusal from semantic validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationError {
    /// One declared collection/resource exceeded its explicit cap.
    LimitExceeded {
        /// Stable resource name.
        resource: &'static str,
        /// Requested units.
        requested: usize,
        /// Admitted units.
        limit: usize,
    },
    /// Checked work-plan arithmetic overflowed.
    WorkPlanOverflow {
        /// Phase whose arithmetic overflowed.
        phase: &'static str,
    },
    /// The complete deterministic work plan exceeded the caller's budget.
    WorkExceeded {
        /// Requested logical units.
        requested: u128,
        /// Admitted logical units.
        limit: u128,
    },
    /// A preflighted scratch allocation was refused.
    AllocationRefused {
        /// Stable scratch-resource name.
        resource: &'static str,
        /// Requested elements.
        requested: usize,
    },
    /// Cancellation was observed before any result was published.
    Cancelled {
        /// Deterministic checkpoint phase.
        phase: &'static str,
        /// Conservatively completed planned work units. In-phase checkpoints
        /// report zero until the phase's plan is fully reconciled.
        completed: u128,
        /// Preflighted logical work units, or zero before preflight.
        planned: u128,
    },
}

impl fmt::Display for ValidationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded {
                resource,
                requested,
                limit,
            } => write!(
                f,
                "scenario validation {resource} request {requested} exceeds limit {limit}"
            ),
            Self::WorkPlanOverflow { phase } => {
                write!(f, "scenario validation work plan overflowed during {phase}")
            }
            Self::WorkExceeded { requested, limit } => write!(
                f,
                "scenario validation work request {requested} exceeds limit {limit}"
            ),
            Self::AllocationRefused {
                resource,
                requested,
            } => write!(
                f,
                "scenario validation allocation for {requested} {resource} elements was refused"
            ),
            Self::Cancelled {
                phase,
                completed,
                planned,
            } => write!(
                f,
                "scenario validation cancelled during {phase} after {completed}/{planned} work units"
            ),
        }
    }
}

impl core::error::Error for ValidationError {}

impl ValidationError {
    pub(crate) fn into_violation(self) -> Violation {
        Violation {
            code: "validation-resource-refused",
            what: self.to_string(),
            fix: "reduce the scenario or call validate_with_budget using an explicitly reviewed larger budget"
                .to_string(),
        }
    }
}

/// One validity finding with its structured fix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Violation {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// What is wrong, with context.
    pub what: String,
    /// How to fix it.
    pub fix: String,
}

/// Explicit environmental fields — required, never silently defaulted.
#[derive(Debug, Clone, PartialEq)]
pub struct Environment {
    /// Gravity vector (m/s², world frame).
    pub gravity: [QtyAny; 3],
    /// Ambient temperature (K).
    pub ambient_temperature: QtyAny,
    /// Ambient pressure (Pa).
    pub ambient_pressure: QtyAny,
}

impl Environment {
    /// Standard Earth laboratory environment (explicitly chosen, not a
    /// hidden default: the call site names it).
    #[must_use]
    pub fn earth_lab() -> Self {
        Environment {
            gravity: [
                QtyAny::new(0.0, ACCEL_DIMS),
                QtyAny::new(0.0, ACCEL_DIMS),
                QtyAny::new(-9.806_65, ACCEL_DIMS),
            ],
            ambient_temperature: QtyAny::new(293.15, TEMP_DIMS),
            ambient_pressure: QtyAny::new(101_325.0, PRESSURE_DIMS),
        }
    }

    fn check(&self, out: &mut Vec<Violation>) {
        for (i, g) in self.gravity.iter().enumerate() {
            if g.dims != ACCEL_DIMS {
                out.push(Violation {
                    code: "env-gravity-dims",
                    what: format!("gravity component {i} has dimensions {:?}", g.dims.0),
                    fix: "express gravity in m/s² (SI exponents [1,0,-2,0,0,0])".to_string(),
                });
            }
            if !g.value.is_finite() {
                out.push(Violation {
                    code: "env-gravity-nonfinite",
                    what: format!("gravity component {i} is non-finite ({})", g.value),
                    fix: "replace every gravity component with a finite acceleration".to_string(),
                });
            }
        }
        if self.ambient_temperature.dims != TEMP_DIMS {
            out.push(Violation {
                code: "env-temperature-dims",
                what: format!(
                    "ambient temperature has dimensions {:?}",
                    self.ambient_temperature.dims.0
                ),
                fix: "express ambient temperature in kelvin".to_string(),
            });
        }
        if !(self.ambient_temperature.value.is_finite() && self.ambient_temperature.value >= 0.0) {
            out.push(Violation {
                code: "env-temperature-range",
                what: format!(
                    "ambient absolute temperature {} K is non-finite or below absolute zero",
                    self.ambient_temperature.value
                ),
                fix: "use a finite absolute temperature greater than or equal to 0 K".to_string(),
            });
        }
        if self.ambient_pressure.dims != PRESSURE_DIMS {
            out.push(Violation {
                code: "env-pressure-dims",
                what: format!(
                    "ambient pressure has dimensions {:?}",
                    self.ambient_pressure.dims.0
                ),
                fix: "express ambient pressure in pascals".to_string(),
            });
        }
        if !(self.ambient_pressure.value.is_finite() && self.ambient_pressure.value >= 0.0) {
            out.push(Violation {
                code: "env-pressure-range",
                what: format!(
                    "ambient absolute pressure {} Pa is non-finite or negative",
                    self.ambient_pressure.value
                ),
                fix: "use a finite, nonnegative absolute pressure".to_string(),
            });
        }
    }
}

/// A named set of loads applied together.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadCase {
    /// Case name (referenced by combinations).
    pub name: String,
    /// The case's boundary conditions / loads.
    pub bcs: Vec<BoundaryCondition>,
}

/// A factored combination of load cases (code-style, e.g. `1.2D + 1.6L`).
#[derive(Debug, Clone, PartialEq)]
pub struct Combination {
    /// Combination name.
    pub name: String,
    /// `(case name, factor)` terms.
    pub terms: Vec<(String, f64)>,
}

/// A contact/friction pairing between two named regions.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactLaw {
    /// First region.
    pub region_a: String,
    /// Second region.
    pub region_b: String,
    /// The friction model.
    pub model: ContactModel,
}

/// Supported contact models.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContactModel {
    /// No tangential resistance.
    Frictionless,
    /// Coulomb friction with static/kinetic coefficients.
    Coulomb {
        /// Static coefficient.
        mu_s: f64,
        /// Kinetic coefficient (≤ static).
        mu_k: f64,
    },
    /// Fully tied (no relative motion).
    Tied,
}

struct ValidationIndex<K> {
    entries: Vec<(K, usize)>,
}

fn sift_validation_heap<T: Ord, E>(
    values: &mut [T],
    mut root: usize,
    end: usize,
    phase: &'static str,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<(), E> {
    loop {
        checkpoint(phase)?;
        let Some(left) = root.checked_mul(2).and_then(|value| value.checked_add(1)) else {
            return Ok(());
        };
        if left >= end {
            return Ok(());
        }
        let right = left + 1;
        let child = if right < end && values[left] < values[right] {
            right
        } else {
            left
        };
        if values[root] >= values[child] {
            return Ok(());
        }
        values.swap(root, child);
        root = child;
    }
}

pub(crate) fn sort_validation_index<T: Ord, E>(
    values: &mut [T],
    phase: &'static str,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<(), E> {
    let len = values.len();
    for root in (0..len / 2).rev() {
        sift_validation_heap(values, root, len, phase, checkpoint)?;
    }
    for end in (1..len).rev() {
        checkpoint(phase)?;
        values.swap(0, end);
        sift_validation_heap(values, 0, end, phase, checkpoint)?;
    }
    Ok(())
}

impl<K: Ord> ValidationIndex<K> {
    fn build(
        entries: impl ExactSizeIterator<Item = (K, usize)>,
        resource: &'static str,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<Self, ValidationError> {
        let requested = entries.len();
        let mut sorted = Vec::new();
        reserve_validation(&mut sorted, requested, resource)?;
        for entry in entries {
            checkpoint(resource)?;
            sorted.push(entry);
        }
        sort_validation_index(&mut sorted, resource, checkpoint)?;
        Ok(Self { entries: sorted })
    }

    fn first_row(&self, key: &K) -> Option<usize> {
        let position = self
            .entries
            .partition_point(|(candidate, _)| candidate < key);
        self.entries
            .get(position)
            .and_then(|(candidate, row)| (candidate == key).then_some(*row))
    }

    fn contains(&self, key: &K) -> bool {
        self.first_row(key).is_some()
    }
}

fn contact_pair(contact: &ContactLaw) -> (&str, &str) {
    if contact.region_a <= contact.region_b {
        (contact.region_a.as_str(), contact.region_b.as_str())
    } else {
        (contact.region_b.as_str(), contact.region_a.as_str())
    }
}

fn enforce_limit(
    resource: &'static str,
    requested: usize,
    limit: usize,
) -> Result<(), ValidationError> {
    if requested > limit {
        Err(ValidationError::LimitExceeded {
            resource,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}

fn checked_count_add(
    total: &mut usize,
    amount: usize,
    phase: &'static str,
) -> Result<(), ValidationError> {
    *total = total
        .checked_add(amount)
        .ok_or(ValidationError::WorkPlanOverflow { phase })?;
    Ok(())
}

fn checked_finding_slots(
    total: &mut usize,
    records: usize,
    per_record: usize,
    phase: &'static str,
) -> Result<(), ValidationError> {
    let slots = records
        .checked_mul(per_record)
        .ok_or(ValidationError::WorkPlanOverflow { phase })?;
    checked_count_add(total, slots, phase)
}

fn work_units(value: usize, phase: &'static str) -> Result<u128, ValidationError> {
    u128::try_from(value).map_err(|_| ValidationError::WorkPlanOverflow { phase })
}

fn checked_work_add(
    total: &mut u128,
    amount: u128,
    phase: &'static str,
) -> Result<(), ValidationError> {
    *total = total
        .checked_add(amount)
        .ok_or(ValidationError::WorkPlanOverflow { phase })?;
    Ok(())
}

fn checked_work_product(
    left: usize,
    right: usize,
    phase: &'static str,
) -> Result<u128, ValidationError> {
    work_units(left, phase)?
        .checked_mul(work_units(right, phase)?)
        .ok_or(ValidationError::WorkPlanOverflow { phase })
}

fn ceil_log2(value: usize) -> u128 {
    if value <= 1 {
        0
    } else {
        u128::from(usize::BITS - (value - 1).leading_zeros())
    }
}

fn bc_dynamic_scalars(bc: &BoundaryCondition) -> Result<usize, ValidationError> {
    match &bc.value {
        Some(BcValue::Signal(signal)) => {
            signal
                .validation_dynamic_scalars()
                .ok_or(ValidationError::WorkPlanOverflow {
                    phase: "signal scalar count",
                })
        }
        Some(BcValue::Profile(profile)) => Ok(profile.cheb.coeffs().len()),
        Some(BcValue::Uniform(_)) | None => Ok(0),
    }
}

fn bc_flux_checkpoints(bc: &BoundaryCondition) -> usize {
    if bc.kind != BcKind::MassFlowInlet {
        return 0;
    }
    match &bc.value {
        Some(BcValue::Signal(signal)) => signal.net_flux_validation_time_count(),
        _ => 0,
    }
}

fn reserve_validation<T>(
    values: &mut Vec<T>,
    requested: usize,
    resource: &'static str,
) -> Result<(), ValidationError> {
    values
        .try_reserve_exact(requested)
        .map_err(|_| ValidationError::AllocationRefused {
            resource,
            requested,
        })
}

fn reserve_validation_times(times: &mut Vec<f64>, requested: usize) -> Result<(), ValidationError> {
    reserve_validation(times, requested, "net-flux checkpoints")
}

/// The scenario root value.
#[derive(Debug, Clone, PartialEq)]
pub struct Scenario {
    /// Scenario name (IR identity).
    pub name: String,
    /// Study seed (one of the Five Explicits).
    pub seed: u64,
    /// Reference frames (world implicit).
    pub frames: FrameTree,
    /// Always-active boundary conditions.
    pub base_bcs: Vec<BoundaryCondition>,
    /// Named load cases.
    pub cases: Vec<LoadCase>,
    /// Factored combinations over the cases.
    pub combinations: Vec<Combination>,
    /// Stochastic ensembles.
    pub ensembles: Vec<StochasticEnsemble>,
    /// Contact laws.
    pub contacts: Vec<ContactLaw>,
    /// The explicit environment.
    pub environment: Environment,
}

impl Scenario {
    /// A scenario with no loads yet. The environment is a REQUIRED
    /// argument — there is deliberately no default.
    #[must_use]
    pub fn new(name: &str, seed: u64, environment: Environment) -> Self {
        Scenario {
            name: name.to_string(),
            seed,
            frames: FrameTree::new(),
            base_bcs: Vec::new(),
            cases: Vec::new(),
            combinations: Vec::new(),
            ensembles: Vec::new(),
            contacts: Vec::new(),
            environment,
        }
    }

    /// Preflight the complete deterministic semantic-validation work plan.
    ///
    /// Top-level collection caps are checked before their contents are
    /// traversed. Nested collection, signal-scalar, identity-byte, checkpoint,
    /// and work totals use checked arithmetic and are compared with the
    /// caller's explicit limits before validation indexes are allocated.
    pub fn validation_plan(
        &self,
        budget: ValidationBudget,
    ) -> Result<ValidationPlan, ValidationError> {
        let mut checkpoint = |_: &'static str| Ok(());
        self.validation_plan_with_checkpoint(budget, &mut checkpoint)
    }

    fn validation_plan_with_checkpoint(
        &self,
        budget: ValidationBudget,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<ValidationPlan, ValidationError> {
        let frames = self.frames.frames.len();
        let base_bcs = self.base_bcs.len();
        let cases = self.cases.len();
        let combinations = self.combinations.len();
        let ensembles = self.ensembles.len();
        let contacts = self.contacts.len();
        enforce_limit("frames", frames, budget.max_frames)?;
        enforce_limit("base boundary conditions", base_bcs, budget.max_base_bcs)?;
        enforce_limit("load cases", cases, budget.max_cases)?;
        enforce_limit("combinations", combinations, budget.max_combinations)?;
        enforce_limit("ensembles", ensembles, budget.max_ensembles)?;
        enforce_limit("contacts", contacts, budget.max_contacts)?;
        checkpoint("preflight collection caps")?;

        let mut case_bcs = 0usize;
        for case in &self.cases {
            checked_count_add(
                &mut case_bcs,
                case.bcs.len(),
                "case boundary-condition count",
            )?;
            checkpoint("preflight case boundary counts")?;
        }
        enforce_limit("case boundary conditions", case_bcs, budget.max_case_bcs)?;

        let mut combination_terms = 0usize;
        for combination in &self.combinations {
            checked_count_add(
                &mut combination_terms,
                combination.terms.len(),
                "combination term count",
            )?;
            checkpoint("preflight combination term counts")?;
        }
        enforce_limit(
            "combination terms",
            combination_terms,
            budget.max_combination_terms,
        )?;

        let mut identity_bytes = self.name.len();
        let mut signal_scalars = 0usize;
        for frame in &self.frames.frames {
            checked_count_add(
                &mut identity_bytes,
                frame.name.len(),
                "frame identity bytes",
            )?;
            if let FrameMotion::Tilt { angle, .. } = &frame.motion {
                checked_count_add(
                    &mut signal_scalars,
                    angle.validation_dynamic_scalars().ok_or(
                        ValidationError::WorkPlanOverflow {
                            phase: "frame signal scalar count",
                        },
                    )?,
                    "frame signal scalar count",
                )?;
            }
            checkpoint("preflight frames")?;
        }

        let mut base_flux_checkpoint_contribution = 0usize;
        for bc in &self.base_bcs {
            checked_count_add(
                &mut identity_bytes,
                bc.region.len(),
                "base boundary identity bytes",
            )?;
            checked_count_add(
                &mut signal_scalars,
                bc_dynamic_scalars(bc)?,
                "base boundary signal scalars",
            )?;
            checked_count_add(
                &mut base_flux_checkpoint_contribution,
                bc_flux_checkpoints(bc),
                "base flux checkpoints",
            )?;
            checkpoint("preflight base boundary conditions")?;
        }
        for case in &self.cases {
            checked_count_add(&mut identity_bytes, case.name.len(), "case identity bytes")?;
            for bc in &case.bcs {
                checked_count_add(
                    &mut identity_bytes,
                    bc.region.len(),
                    "case boundary identity bytes",
                )?;
                checked_count_add(
                    &mut signal_scalars,
                    bc_dynamic_scalars(bc)?,
                    "case boundary signal scalars",
                )?;
                checkpoint("preflight case boundary conditions")?;
            }
            checkpoint("preflight load cases")?;
        }
        for combination in &self.combinations {
            checked_count_add(
                &mut identity_bytes,
                combination.name.len(),
                "combination identity bytes",
            )?;
            for (case, _) in &combination.terms {
                checked_count_add(
                    &mut identity_bytes,
                    case.len(),
                    "combination reference bytes",
                )?;
                checkpoint("preflight combination terms")?;
            }
            checkpoint("preflight combinations")?;
        }
        for ensemble in &self.ensembles {
            checked_count_add(
                &mut identity_bytes,
                ensemble.name.len(),
                "ensemble identity bytes",
            )?;
            checkpoint("preflight ensembles")?;
        }
        for contact in &self.contacts {
            checked_count_add(
                &mut identity_bytes,
                contact.region_a.len(),
                "contact identity bytes",
            )?;
            checked_count_add(
                &mut identity_bytes,
                contact.region_b.len(),
                "contact identity bytes",
            )?;
            checkpoint("preflight contacts")?;
        }
        enforce_limit("signal scalars", signal_scalars, budget.max_signal_scalars)?;
        enforce_limit("identity bytes", identity_bytes, budget.max_identity_bytes)?;

        let total_bcs =
            base_bcs
                .checked_add(case_bcs)
                .ok_or(ValidationError::WorkPlanOverflow {
                    phase: "finding boundary-condition count",
                })?;
        let mut finding_capacity = FIXED_FINDING_SLOTS;
        for (records, per_record, phase) in [
            (frames, FINDINGS_PER_FRAME, "frame finding capacity"),
            (
                total_bcs,
                FINDINGS_PER_BC,
                "boundary-condition finding capacity",
            ),
            (cases, FINDINGS_PER_CASE, "case finding capacity"),
            (
                combinations,
                FINDINGS_PER_COMBINATION,
                "combination finding capacity",
            ),
            (
                combination_terms,
                FINDINGS_PER_COMBINATION_TERM,
                "combination-term finding capacity",
            ),
            (
                ensembles,
                FINDINGS_PER_ENSEMBLE,
                "ensemble finding capacity",
            ),
            (contacts, FINDINGS_PER_CONTACT, "contact finding capacity"),
        ] {
            checked_finding_slots(&mut finding_capacity, records, per_record, phase)?;
        }
        enforce_limit("validation findings", finding_capacity, budget.max_findings)?;

        let mut flux_checkpoints = base_flux_checkpoint_contribution.checked_add(1).ok_or(
            ValidationError::WorkPlanOverflow {
                phase: "base flux checkpoint count",
            },
        )?;
        let mut flux_sort_work = work_units(flux_checkpoints, "base flux sort work")?
            .checked_mul(ceil_log2(flux_checkpoints))
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "base flux sort work",
            })?;
        let mut flux_evaluation_work =
            checked_work_product(flux_checkpoints, base_bcs, "base flux evaluation work")?;
        for case in &self.cases {
            let mut case_contribution = 0usize;
            for bc in &case.bcs {
                checked_count_add(
                    &mut case_contribution,
                    bc_flux_checkpoints(bc),
                    "case flux checkpoints",
                )?;
                checkpoint("preflight case flux checkpoints")?;
            }
            let set_checkpoints = 1usize
                .checked_add(base_flux_checkpoint_contribution)
                .and_then(|count| count.checked_add(case_contribution))
                .ok_or(ValidationError::WorkPlanOverflow {
                    phase: "effective-set flux checkpoint count",
                })?;
            checked_count_add(
                &mut flux_checkpoints,
                set_checkpoints,
                "total flux checkpoints",
            )?;
            let sort = work_units(set_checkpoints, "case flux sort work")?
                .checked_mul(ceil_log2(set_checkpoints))
                .ok_or(ValidationError::WorkPlanOverflow {
                    phase: "case flux sort work",
                })?;
            checked_work_add(&mut flux_sort_work, sort, "total flux sort work")?;
            let set_bcs =
                base_bcs
                    .checked_add(case.bcs.len())
                    .ok_or(ValidationError::WorkPlanOverflow {
                        phase: "effective-set boundary count",
                    })?;
            checked_work_add(
                &mut flux_evaluation_work,
                checked_work_product(set_checkpoints, set_bcs, "case flux evaluation work")?,
                "total flux evaluation work",
            )?;
            checkpoint("preflight case flux sets")?;
        }
        enforce_limit(
            "flux checkpoints",
            flux_checkpoints,
            budget.max_flux_checkpoints,
        )?;

        let mut record_visits = 1usize;
        for count in [
            frames,
            base_bcs,
            cases,
            case_bcs,
            combinations,
            combination_terms,
            ensembles,
            contacts,
        ] {
            checked_count_add(&mut record_visits, count, "validation record visits")?;
        }
        let mut index_items = frames
            .checked_mul(3)
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "frame index item count",
            })?;
        for count in [cases, combinations, combination_terms, ensembles, contacts] {
            checked_count_add(&mut index_items, count, "validation index item count")?;
        }
        let largest_index = [
            frames,
            cases,
            combinations,
            combination_terms,
            ensembles,
            contacts,
        ]
        .into_iter()
        .max()
        .unwrap_or(0)
        .checked_add(1)
        .ok_or(ValidationError::WorkPlanOverflow {
            phase: "validation index height",
        })?;
        let index_height = ceil_log2(largest_index);
        let index_work_per_item = index_height
            .checked_mul(2)
            .and_then(|work| work.checked_add(2))
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "validation index work multiplier",
            })?;
        let index_work = work_units(index_items, "validation index work")?
            .checked_mul(index_work_per_item)
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "validation index work",
            })?;
        let identity_work = work_units(identity_bytes, "identity comparison work")?
            .checked_mul(index_height.max(1))
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "identity comparison work",
            })?;
        let mut planned_work = work_units(record_visits, "validation record work")?;
        for (amount, phase) in [
            (index_work, "validation index work"),
            (identity_work, "identity comparison work"),
            (
                work_units(signal_scalars, "signal scalar work")?,
                "signal scalar work",
            ),
            (
                checked_work_product(contacts, 3, "contact classification work")?,
                "contact classification work",
            ),
            (flux_sort_work, "flux sort work"),
            (flux_evaluation_work, "flux evaluation work"),
        ] {
            checked_work_add(&mut planned_work, amount, phase)?;
        }
        if planned_work > budget.max_work {
            return Err(ValidationError::WorkExceeded {
                requested: planned_work,
                limit: budget.max_work,
            });
        }
        Ok(ValidationPlan {
            frames,
            base_bcs,
            cases,
            case_bcs,
            combinations,
            combination_terms,
            ensembles,
            contacts,
            signal_scalars,
            flux_checkpoints,
            identity_bytes,
            finding_capacity,
            planned_work,
        })
    }

    /// Validate under [`DEFAULT_VALIDATION_BUDGET`]; empty means admissible.
    /// Resource refusal is retained as a structured violation so legacy
    /// callers cannot mistake an unadmitted scenario for a green one.
    #[must_use]
    pub fn validate(&self) -> Vec<Violation> {
        match self.validation_plan(ValidationBudget::default()) {
            Ok(plan) => {
                let mut checkpoint = |_: &'static str| Ok(());
                match self.validate_admitted(plan.finding_capacity, &mut checkpoint) {
                    Ok(findings) => findings,
                    Err(error) => vec![error.into_violation()],
                }
            }
            Err(error) => vec![error.into_violation()],
        }
    }

    /// Validate under an explicit semantic budget and cancellation context.
    ///
    /// Cancellation is polled before preflight, after the work plan is known,
    /// across semantic phases and record boundaries, and after private
    /// validation but before findings are published. Thus an observed request
    /// returns a typed refusal, never a partial or falsely green finding list.
    pub fn validate_with_budget(
        &self,
        budget: ValidationBudget,
        cx: &Cx<'_>,
    ) -> Result<Vec<Violation>, ValidationError> {
        cx.checkpoint().map_err(|_| ValidationError::Cancelled {
            phase: "initial",
            completed: 0,
            planned: 0,
        })?;
        let mut preflight_checkpoint = |phase: &'static str| {
            cx.checkpoint().map_err(|_| ValidationError::Cancelled {
                phase,
                completed: 0,
                planned: 0,
            })
        };
        let plan = self.validation_plan_with_checkpoint(budget, &mut preflight_checkpoint)?;
        cx.checkpoint().map_err(|_| ValidationError::Cancelled {
            phase: "post-preflight",
            completed: 0,
            planned: plan.planned_work,
        })?;
        let mut checkpoint = |phase: &'static str| {
            cx.checkpoint().map_err(|_| ValidationError::Cancelled {
                phase,
                completed: 0,
                planned: plan.planned_work,
            })
        };
        let findings = self.validate_admitted(plan.finding_capacity, &mut checkpoint)?;
        cx.checkpoint().map_err(|_| ValidationError::Cancelled {
            phase: "pre-publication",
            completed: plan.planned_work,
            planned: plan.planned_work,
        })?;
        Ok(findings)
    }

    fn validate_admitted(
        &self,
        finding_capacity: usize,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<Vec<Violation>, ValidationError> {
        let mut out = Vec::new();
        reserve_validation(&mut out, finding_capacity, "validation findings")?;
        if self.name.is_empty() {
            out.push(Violation {
                code: "scenario-name-empty",
                what: "scenario identity is empty".to_string(),
                fix: "give the scenario a nonempty exact UTF-8 name".to_string(),
            });
        }
        checkpoint("scenario identity")?;
        self.environment.check(&mut out);
        checkpoint("environment")?;
        self.frames.check_with_checkpoint(&mut out, checkpoint)?;
        checkpoint("frames")?;
        let frame_ids = ValidationIndex::build(
            self.frames
                .frames
                .iter()
                .enumerate()
                .map(|(row, frame)| (frame.id.0, row)),
            "scenario frame id index",
            checkpoint,
        )?;
        let first_case_by_name = ValidationIndex::build(
            self.cases
                .iter()
                .enumerate()
                .map(|(row, case)| (case.name.as_str(), row)),
            "case name index",
            checkpoint,
        )?;
        for bc in &self.base_bcs {
            bc.check_with_checkpoint(&mut out, checkpoint)?;
            Self::check_bc_frame(bc, &frame_ids, &mut out);
            checkpoint("base boundary conditions")?;
        }

        for (i, case) in self.cases.iter().enumerate() {
            if case.name.is_empty() {
                out.push(Violation {
                    code: "case-name-empty",
                    what: format!("load case row {i} has an empty identity"),
                    fix: "give every load case a nonempty exact UTF-8 name".to_string(),
                });
            }
            let first = first_case_by_name
                .first_row(&case.name.as_str())
                .unwrap_or(i);
            if first != i {
                out.push(Violation {
                    code: "case-name-duplicate",
                    what: format!(
                        "load case {:?} first appears at row {first} and repeats at row {i}",
                        case.name
                    ),
                    fix: "give every load case a unique name".to_string(),
                });
            }
            for bc in &case.bcs {
                bc.check_with_checkpoint(&mut out, checkpoint)?;
                Self::check_bc_frame(bc, &frame_ids, &mut out);
                checkpoint("case boundary conditions")?;
            }
            checkpoint("load cases")?;
        }

        let first_combo_by_name = ValidationIndex::build(
            self.combinations
                .iter()
                .enumerate()
                .map(|(row, combo)| (combo.name.as_str(), row)),
            "combination name index",
            checkpoint,
        )?;
        for (combo_index, combo) in self.combinations.iter().enumerate() {
            if combo.name.is_empty() {
                out.push(Violation {
                    code: "combo-name-empty",
                    what: format!("combination row {combo_index} has an empty identity"),
                    fix: "give every combination a nonempty exact UTF-8 name".to_string(),
                });
            }
            let first = first_combo_by_name
                .first_row(&combo.name.as_str())
                .unwrap_or(combo_index);
            if first != combo_index {
                out.push(Violation {
                    code: "combo-name-duplicate",
                    what: format!(
                        "combination {:?} first appears at row {first} and repeats at row {combo_index}",
                        combo.name
                    ),
                    fix: "give every combination a unique name".to_string(),
                });
            }

            let first_term_by_case = ValidationIndex::build(
                combo
                    .terms
                    .iter()
                    .enumerate()
                    .map(|(row, (case, _))| (case.as_str(), row)),
                "combination term index",
                checkpoint,
            )?;
            for (term_index, (case, factor)) in combo.terms.iter().enumerate() {
                if case.is_empty() {
                    out.push(Violation {
                        code: "combo-case-empty",
                        what: format!(
                            "combination {:?} term {term_index} has an empty case reference",
                            combo.name
                        ),
                        fix: "reference a nonempty defined load-case name".to_string(),
                    });
                }
                let first = first_term_by_case
                    .first_row(&case.as_str())
                    .unwrap_or(term_index);
                if first != term_index {
                    out.push(Violation {
                        code: "combo-term-duplicate",
                        what: format!(
                            "combination {:?} references case {case:?} at terms {first} and {term_index}",
                            combo.name
                        ),
                        fix: "combine repeated case factors into one term".to_string(),
                    });
                }
                if !first_case_by_name.contains(&case.as_str()) {
                    out.push(Violation {
                        code: "combo-case-missing",
                        what: format!(
                            "combination {:?} references unknown case {case:?}",
                            combo.name
                        ),
                        fix: "reference only defined load cases".to_string(),
                    });
                }
                if !factor.is_finite() {
                    out.push(Violation {
                        code: "combo-factor",
                        what: format!("combination {:?} has non-finite factor", combo.name),
                        fix: "use finite combination factors".to_string(),
                    });
                }
                checkpoint("combination terms")?;
            }
            checkpoint("combinations")?;
        }

        let first_ensemble_by_name = ValidationIndex::build(
            self.ensembles
                .iter()
                .enumerate()
                .map(|(row, ensemble)| (ensemble.name.as_str(), row)),
            "ensemble name index",
            checkpoint,
        )?;
        for (ensemble_index, e) in self.ensembles.iter().enumerate() {
            e.check(&mut out);
            let first = first_ensemble_by_name
                .first_row(&e.name.as_str())
                .unwrap_or(ensemble_index);
            if first != ensemble_index {
                out.push(Violation {
                    code: "ensemble-name-duplicate",
                    what: format!(
                        "ensemble {:?} first appears at row {first} and repeats at row {ensemble_index}",
                        e.name
                    ),
                    fix: "give every ensemble a unique name".to_string(),
                });
            }
            checkpoint("ensembles")?;
        }

        let first_contact_by_pair = ValidationIndex::build(
            self.contacts
                .iter()
                .enumerate()
                .map(|(row, contact)| (contact_pair(contact), row)),
            "contact pair index",
            checkpoint,
        )?;
        let contact_pair_conflicts =
            self.contact_pair_conflicts(&first_contact_by_pair, checkpoint)?;
        for (contact_index, c) in self.contacts.iter().enumerate() {
            if c.region_a.is_empty() || c.region_b.is_empty() {
                out.push(Violation {
                    code: "contact-region-empty",
                    what: format!(
                        "contact row {contact_index} has empty region identity {:?}/{:?}",
                        c.region_a, c.region_b
                    ),
                    fix: "bind both sides to nonempty exact UTF-8 region names".to_string(),
                });
            }
            if c.region_a == c.region_b {
                out.push(Violation {
                    code: "contact-self-pair",
                    what: format!(
                        "contact row {contact_index} pairs region {:?} with itself",
                        c.region_a
                    ),
                    fix: "name two distinct contact regions".to_string(),
                });
            }
            let pair = contact_pair(c);
            let first = first_contact_by_pair
                .first_row(&pair)
                .unwrap_or(contact_index);
            if first != contact_index {
                let (code, fix) = if contact_pair_conflicts
                    .get(contact_index)
                    .copied()
                    .unwrap_or(true)
                {
                    (
                        "contact-pair-conflict",
                        "choose one contact model for the unordered region pair",
                    )
                } else {
                    (
                        "contact-pair-duplicate",
                        "remove the repeated unordered contact pair",
                    )
                };
                out.push(Violation {
                    code,
                    what: format!(
                        "unordered contact pair {pair:?} first appears at row {first} and repeats at row {contact_index}"
                    ),
                    fix: fix.to_string(),
                });
            }
            if let ContactModel::Coulomb { mu_s, mu_k } = c.model
                && !(mu_s.is_finite() && mu_k.is_finite() && mu_k >= 0.0 && mu_s >= mu_k)
            {
                out.push(Violation {
                    code: "contact-coulomb-range",
                    what: format!(
                        "contact {:?}/{:?}: mu_s={mu_s}, mu_k={mu_k}",
                        c.region_a, c.region_b
                    ),
                    fix: "require 0 <= mu_k <= mu_s < inf".to_string(),
                });
            }
            checkpoint("contacts")?;
        }
        self.check_net_flux(&mut out, checkpoint)?;
        if out.len() > finding_capacity {
            return Err(ValidationError::WorkPlanOverflow {
                phase: "validation finding capacity invariant",
            });
        }
        Ok(out)
    }

    fn contact_pair_conflicts(
        &self,
        index: &ValidationIndex<(&str, &str)>,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<Vec<bool>, ValidationError> {
        let contact_count = self.contacts.len();
        let mut conflicts = Vec::new();
        reserve_validation(&mut conflicts, contact_count, "contact conflict flags")?;
        conflicts.resize(contact_count, false);

        let mut start = 0usize;
        while start < index.entries.len() {
            let key = &index.entries[start].0;
            let mut end = start + 1;
            while end < index.entries.len() && &index.entries[end].0 == key {
                checkpoint("contact conflict grouping")?;
                end += 1;
            }
            let first_row = index.entries[start].1;
            let first_model = self.contacts.get(first_row).map(|contact| &contact.model);
            let mut group_conflicts = first_model.is_none();
            for (_, row) in &index.entries[start + 1..end] {
                checkpoint("contact conflict models")?;
                group_conflicts |= self
                    .contacts
                    .get(*row)
                    .map(|contact| Some(&contact.model) != first_model)
                    .unwrap_or(true);
            }
            if group_conflicts {
                for (_, row) in &index.entries[start..end] {
                    checkpoint("contact conflict publication")?;
                    let Some(conflict) = conflicts.get_mut(*row) else {
                        return Err(ValidationError::WorkPlanOverflow {
                            phase: "contact conflict row invariant",
                        });
                    };
                    *conflict = true;
                }
            }
            start = end;
        }
        Ok(conflicts)
    }

    fn check_bc_frame(
        bc: &BoundaryCondition,
        frame_ids: &ValidationIndex<u32>,
        out: &mut Vec<Violation>,
    ) {
        if bc.frame != 0 && !frame_ids.contains(&bc.frame) {
            out.push(Violation {
                code: "bc-frame-missing",
                what: format!(
                    "bc on {:?} references unknown frame {}",
                    bc.region, bc.frame
                ),
                fix: "reference a defined frame id (or 0 for the world)".to_string(),
            });
        }
    }

    /// Net-flux compatibility (the admission check the bead names): for
    /// every effective BC set (base alone, and base + each case), if any
    /// inlet declares `incompressible`, the declared mass flows must
    /// balance to tolerance at every deterministic signal checkpoint OR a
    /// pressure outlet must exist to absorb the imbalance.
    fn check_net_flux(
        &self,
        out: &mut Vec<Violation>,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<(), ValidationError> {
        self.check_net_flux_set(None, out, checkpoint)?;
        for case in &self.cases {
            self.check_net_flux_set(Some(case), out, checkpoint)?;
        }
        Ok(())
    }

    fn check_net_flux_set(
        &self,
        case: Option<&LoadCase>,
        out: &mut Vec<Violation>,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<(), ValidationError> {
        checkpoint("net-flux set")?;
        let case_bcs: &[BoundaryCondition] = case.map_or(&[], |case| case.bcs.as_slice());
        let declares_incompressible = self
            .base_bcs
            .iter()
            .chain(case_bcs.iter())
            .any(|bc| bc.compatibility == Some(Compat::Incompressible));
        if !declares_incompressible {
            return Ok(());
        }
        let has_pressure_outlet = self.base_bcs.iter().chain(case_bcs.iter()).any(|bc| {
            bc.physics == Physics::IncompressibleFlow && bc.kind == BcKind::PressureOutlet
        });
        let label = case.map_or_else(|| "base".to_string(), |case| format!("base+{}", case.name));
        let mut requested_times = 1usize;
        for bc in self.base_bcs.iter().chain(case_bcs.iter()) {
            requested_times = requested_times.checked_add(bc_flux_checkpoints(bc)).ok_or(
                ValidationError::WorkPlanOverflow {
                    phase: "net-flux checkpoint materialization",
                },
            )?;
        }
        let mut validation_times = Vec::new();
        reserve_validation_times(&mut validation_times, requested_times)?;
        validation_times.push(0.0);
        for bc in self.base_bcs.iter().chain(case_bcs.iter()) {
            if bc.kind == BcKind::MassFlowInlet
                && let Some(BcValue::Signal(signal)) = &bc.value
            {
                signal.append_net_flux_validation_times(&mut validation_times);
            }
        }
        debug_assert_eq!(validation_times.len(), requested_times);
        validation_times.sort_by(|a, b| a.total_cmp(b));
        validation_times.dedup_by(|a, b| *a == *b);
        checkpoint("net-flux checkpoints")?;

        let mut first_imbalance = None;
        let mut compatibility_failed = false;
        for time in validation_times {
            let mut net = 0.0f64;
            let mut gross = 0.0f64;
            let mut aggregation_finite = true;
            let mut evaluation_failed = false;
            for bc in self.base_bcs.iter().chain(case_bcs.iter()) {
                checkpoint("net-flux evaluation")?;
                match bc.mass_flow_at(time) {
                    Ok(Some(flux)) => {
                        net += flux;
                        gross += flux.abs();
                        if !flux.is_finite() || !net.is_finite() || !gross.is_finite() {
                            aggregation_finite = false;
                            break;
                        }
                    }
                    Ok(None) => {}
                    Err(error) => {
                        out.push(Violation {
                            code: "flux-evaluation",
                            what: format!(
                                "set {label:?} at t={time:.6e} s: declared mass flow could not be evaluated: {error}"
                            ),
                            fix: "repair the total-flow value or use a uniform/time-signal kg/s declaration that is finite at every compatibility checkpoint"
                                .to_string(),
                        });
                        evaluation_failed = true;
                        break;
                    }
                }
                checkpoint("net-flux evaluation")?;
            }
            if evaluation_failed {
                compatibility_failed = true;
                break;
            }
            if !aggregation_finite {
                out.push(Violation {
                    code: "flux-aggregation-nonfinite",
                    what: format!(
                        "set {label:?} at t={time:.6e} s: declared mass-flow aggregation overflowed or contained a non-finite value"
                    ),
                    fix: "rescale or partition the declared mass flows so their finite net/gross balance can be certified"
                        .to_string(),
                });
                compatibility_failed = true;
                break;
            }
            if !has_pressure_outlet
                && first_imbalance.is_none()
                && gross > 0.0
                && net.abs() > FLUX_REL_TOL * gross
            {
                first_imbalance = Some((time, net, gross));
            }
        }
        if compatibility_failed || has_pressure_outlet {
            // The outlet may absorb every finite per-instant imbalance, but it
            // cannot make malformed or unevaluable declarations admissible.
            return Ok(());
        }
        if let Some((time, net, gross)) = first_imbalance {
            out.push(Violation {
                code: "flux-imbalance",
                what: format!(
                    "set {label:?} at t={time:.6e} s: declared incompressible but net mass flow is \
                     {net:+.6e} kg/s over gross {gross:.6e} kg/s with no pressure outlet"
                ),
                fix: "balance the declared inlet/outlet mass flows at every instant or add a \
                      pressure outlet to absorb the imbalance"
                    .to_string(),
            });
        }
        checkpoint("net-flux set complete")?;
        Ok(())
    }
}

#[cfg(test)]
mod validation_internal_tests {
    use super::{
        Combination, ContactLaw, ContactModel, Environment, LoadCase, Scenario, ValidationBudget,
        ValidationError, reserve_validation_times, sort_validation_index,
    };
    use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
    use crate::ensemble::{SpectrumModel, StochasticEnsemble};
    use crate::frame::{Frame, FrameId, FrameMotion};
    use crate::signal::{Interp, TimeSignal};
    use fs_ga::{Quat, Vec3};
    use fs_qty::{Dims, QtyAny};

    #[test]
    fn checkpoint_capacity_overflow_is_a_typed_refusal() {
        let mut times = Vec::new();
        assert!(matches!(
            reserve_validation_times(&mut times, usize::MAX),
            Err(ValidationError::AllocationRefused {
                resource: "net-flux checkpoints",
                requested: usize::MAX,
            })
        ));
        assert!(times.is_empty());
    }

    #[test]
    fn injected_phase_cancellation_discards_private_findings() {
        let scenario = Scenario::new("phase-cancel", 1, Environment::earth_lab());
        let plan = scenario
            .validation_plan(ValidationBudget::default())
            .expect("minimal semantic plan");
        let mut visited = Vec::new();
        let result = scenario.validate_admitted(plan.finding_capacity, &mut |phase| {
            visited.push(phase);
            if phase == "environment" {
                Err(ValidationError::Cancelled {
                    phase,
                    completed: 0,
                    planned: plan.planned_work,
                })
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(ValidationError::Cancelled {
                phase: "environment",
                completed: 0,
                planned,
            }) if planned == plan.planned_work
        ));
        assert_eq!(visited, ["scenario identity", "environment"]);
    }

    #[test]
    fn injected_preflight_cancellation_refuses_the_plan() {
        let mut scenario = Scenario::new("preflight-cancel", 1, Environment::earth_lab());
        scenario.cases.push(LoadCase {
            name: "load".to_string(),
            bcs: Vec::new(),
        });
        let mut visited = Vec::new();
        let result =
            scenario.validation_plan_with_checkpoint(ValidationBudget::default(), &mut |phase| {
                visited.push(phase);
                if phase == "preflight case boundary counts" {
                    Err(ValidationError::Cancelled {
                        phase,
                        completed: 0,
                        planned: 0,
                    })
                } else {
                    Ok(())
                }
            });

        assert!(matches!(
            result,
            Err(ValidationError::Cancelled {
                phase: "preflight case boundary counts",
                completed: 0,
                planned: 0,
            })
        ));
        assert_eq!(
            visited,
            [
                "preflight collection caps",
                "preflight case boundary counts"
            ]
        );
    }

    #[test]
    fn index_sort_is_total_ordered_and_cancellable() {
        let mut values = vec![("b", 2usize), ("a", 3), ("a", 1), ("c", 0)];
        let mut polls = 0usize;
        sort_validation_index(&mut values, "test index sort", &mut |phase| {
            assert_eq!(phase, "test index sort");
            polls += 1;
            Ok::<(), core::convert::Infallible>(())
        })
        .expect("infallible checkpoint");
        assert_eq!(values, [("a", 1), ("a", 3), ("b", 2), ("c", 0)]);
        assert!(polls >= values.len());

        let mut cancelled = vec![3u32, 2, 1, 0];
        let mut calls = 0usize;
        let result = sort_validation_index(&mut cancelled, "cancel sort", &mut |_| {
            calls += 1;
            if calls == 2 { Err("cancelled") } else { Ok(()) }
        });
        assert_eq!(result, Err("cancelled"));
    }

    fn phase_matrix_flow(region: &str, value: f64) -> BoundaryCondition {
        BoundaryCondition {
            region: region.to_string(),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
            value: Some(BcValue::Signal(TimeSignal::Table {
                times: vec![0.0, 1.0],
                values: vec![value, value],
                dims: Dims([0, 1, -1, 0, 0, 0]),
                interp: Interp::Hold,
            })),
            compatibility: Some(Compat::Incompressible),
            frame: 0,
        }
    }

    fn phase_matrix_scenario() -> Scenario {
        let mut scenario = Scenario::new("phase-matrix", 41, Environment::earth_lab());
        scenario.frames.add(Frame {
            id: FrameId(1),
            name: "fixture".to_string(),
            parent: FrameId(0),
            motion: FrameMotion::Fixed {
                orientation: Quat::identity(),
                translation: Vec3::new(0.0, 0.0, 0.0),
            },
        });
        scenario.base_bcs.push(phase_matrix_flow("base-inlet", 1.0));
        scenario.cases.push(LoadCase {
            name: "load".to_string(),
            bcs: vec![phase_matrix_flow("case-outlet", -1.0)],
        });
        scenario.combinations.push(Combination {
            name: "combo".to_string(),
            terms: vec![("load".to_string(), 1.0)],
        });
        scenario.ensembles.push(StochasticEnsemble {
            name: "ensemble".to_string(),
            seed: 9,
            members: 1,
            duration: QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
            dt: QtyAny::new(0.5, Dims([0, 0, 1, 0, 0, 0])),
            model: SpectrumModel::KanaiTajimi {
                s0: 1.0,
                omega_g: QtyAny::new(2.0, Dims([0, 0, -1, 0, 0, 0])),
                zeta_g: 0.5,
            },
        });
        scenario.contacts.extend([
            ContactLaw {
                region_a: "a".to_string(),
                region_b: "b".to_string(),
                model: ContactModel::Frictionless,
            },
            ContactLaw {
                region_a: "b".to_string(),
                region_b: "a".to_string(),
                model: ContactModel::Tied,
            },
        ]);
        scenario
    }

    #[test]
    fn every_representative_phase_discards_private_findings_on_cancellation() {
        let scenario = phase_matrix_scenario();
        let plan = scenario
            .validation_plan(ValidationBudget::default())
            .expect("representative semantic plan");
        let mut phases = Vec::new();
        scenario
            .validate_admitted(plan.finding_capacity, &mut |phase| {
                if !phases.contains(&phase) {
                    phases.push(phase);
                }
                Ok(())
            })
            .expect("phase discovery validation");

        for required in [
            "frame cycle traversal",
            "signal table times",
            "base boundary conditions",
            "case boundary conditions",
            "combination terms",
            "ensembles",
            "contact conflict models",
            "contacts",
            "net-flux checkpoints",
            "net-flux evaluation",
            "net-flux set complete",
        ] {
            assert!(phases.contains(&required), "missing phase {required}");
        }

        for target in phases {
            let mut injected = false;
            let result = scenario.validate_admitted(plan.finding_capacity, &mut |phase| {
                if !injected && phase == target {
                    injected = true;
                    Err(ValidationError::Cancelled {
                        phase,
                        completed: 0,
                        planned: plan.planned_work,
                    })
                } else {
                    Ok(())
                }
            });
            assert!(injected, "phase {target} was not revisited");
            assert!(matches!(
                result,
                Err(ValidationError::Cancelled {
                    phase,
                    completed: 0,
                    planned,
                }) if phase == target && planned == plan.planned_work
            ));
        }
    }
}

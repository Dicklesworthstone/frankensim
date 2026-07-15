//! The `Scenario` root value: frames + BCs + load cases + combinations +
//! ensembles + environment + contact laws, with whole-scenario validity
//! checking. Violations are STRUCTURED FIXES (code, what, fix) — the
//! agent-facing refusal format — never free-text panics. Environmental
//! fields are explicit constructor arguments: nothing is defaulted
//! silently.

use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
use crate::ensemble::StochasticEnsemble;
use crate::frame::{FrameMotion, FrameTree};
use crate::signal::TimeSignal;
use fs_exec::Cx;
use fs_qty::{Dims, QtyAny};
use std::cmp::Ordering;
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
/// Maximum source bytes copied from one identity into a validation finding.
///
/// Row and term coordinates retain exact provenance when distinct identities
/// share the same bounded prefix.
const DIAGNOSTIC_IDENTITY_PREVIEW_BYTES: usize = 128;

#[derive(Clone, Copy)]
struct DiagnosticIdentity<'a>(&'a str);

impl fmt::Debug for DiagnosticIdentity<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.len() <= DIAGNOSTIC_IDENTITY_PREVIEW_BYTES {
            return fmt::Debug::fmt(self.0, formatter);
        }
        let mut end = DIAGNOSTIC_IDENTITY_PREVIEW_BYTES;
        while !self.0.is_char_boundary(end) {
            end -= 1;
        }
        fmt::Debug::fmt(&self.0[..end], formatter)?;
        write!(formatter, "…<{} bytes total>", self.0.len())
    }
}

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
    max_work: 268_435_456,
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

fn sift_validation_heap_by<T, E>(
    values: &mut [T],
    mut root: usize,
    end: usize,
    phase: &'static str,
    compare: &mut impl FnMut(&T, &T) -> Ordering,
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
        let child = if right < end && compare(&values[left], &values[right]).is_lt() {
            right
        } else {
            left
        };
        if !compare(&values[root], &values[child]).is_lt() {
            return Ok(());
        }
        values.swap(root, child);
        root = child;
    }
}

fn sort_validation_by<T, E>(
    values: &mut [T],
    phase: &'static str,
    mut compare: impl FnMut(&T, &T) -> Ordering,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<(), E> {
    let len = values.len();
    for root in (0..len / 2).rev() {
        sift_validation_heap_by(values, root, len, phase, &mut compare, checkpoint)?;
    }
    for end in (1..len).rev() {
        checkpoint(phase)?;
        values.swap(0, end);
        sift_validation_heap_by(values, 0, end, phase, &mut compare, checkpoint)?;
    }
    Ok(())
}

pub(crate) fn sort_validation_index<T: Ord, E>(
    values: &mut [T],
    phase: &'static str,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<(), E> {
    sort_validation_by(values, phase, |left, right| left.cmp(right), checkpoint)
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

fn heap_sort_work(items: usize, phase: &'static str) -> Result<u128, ValidationError> {
    let per_item = ceil_log2(items)
        .checked_mul(2)
        .and_then(|work| work.checked_add(2))
        .ok_or(ValidationError::WorkPlanOverflow { phase })?;
    work_units(items, phase)?
        .checked_mul(per_item)
        .ok_or(ValidationError::WorkPlanOverflow { phase })
}

/// Ordered lookup envelope: global partition height plus candidate equality,
/// duplicate-neighbor ambiguity detection, and one fixed guard.
fn ordered_lookup_work(
    lookups: usize,
    index_height: u128,
    phase: &'static str,
) -> Result<u128, ValidationError> {
    let per_lookup = index_height
        .checked_add(3)
        .ok_or(ValidationError::WorkPlanOverflow { phase })?;
    work_units(lookups, phase)?
        .checked_mul(per_lookup)
        .ok_or(ValidationError::WorkPlanOverflow { phase })
}

/// Conservative byte work for one string-key index and its later lookups.
///
/// `heap_sort_work` bounds comparator calls; the factor of two charges both
/// key operands at the role's maximum UTF-8 width. Ordered lookups charge both
/// query and candidate operands as well. A zero-byte identity still costs one
/// unit so empty-key comparisons and equality checks are admitted.
fn string_index_identity_work(
    items: usize,
    lookups: usize,
    max_key_bytes: usize,
    index_height: u128,
    phase: &'static str,
) -> Result<u128, ValidationError> {
    let key_width = work_units(max_key_bytes.max(1), phase)?;
    let sort_work = heap_sort_work(items, phase)?
        .checked_mul(2)
        .and_then(|work| work.checked_mul(key_width))
        .ok_or(ValidationError::WorkPlanOverflow { phase })?;
    let lookup_work = ordered_lookup_work(lookups, index_height, phase)?
        .checked_mul(2)
        .and_then(|work| work.checked_mul(key_width))
        .ok_or(ValidationError::WorkPlanOverflow { phase })?;
    sort_work
        .checked_add(lookup_work)
        .ok_or(ValidationError::WorkPlanOverflow { phase })
}

fn repeated_identity_comparison_work(
    items: usize,
    comparisons_per_item: usize,
    max_key_bytes: usize,
    phase: &'static str,
) -> Result<u128, ValidationError> {
    checked_work_product(items, comparisons_per_item, phase)?
        .checked_mul(work_units(max_key_bytes.max(1), phase)?)
        .ok_or(ValidationError::WorkPlanOverflow { phase })
}

fn checked_work_scale(
    work: u128,
    count: usize,
    phase: &'static str,
) -> Result<u128, ValidationError> {
    work.checked_mul(work_units(count, phase)?)
        .ok_or(ValidationError::WorkPlanOverflow { phase })
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

fn bc_flux_evaluation_work(bc: &BoundaryCondition) -> Result<u128, ValidationError> {
    if bc.kind != BcKind::MassFlowInlet {
        return Ok(4);
    }
    match &bc.value {
        Some(BcValue::Signal(TimeSignal::Table { times, .. })) => ceil_log2(times.len())
            .checked_add(4)
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "table flux evaluation work",
            }),
        Some(BcValue::Signal(TimeSignal::Chebfun(profile))) => work_units(
            profile.cheb.coeffs().len(),
            "Chebyshev flux evaluation work",
        )?
        .checked_add(4)
        .ok_or(ValidationError::WorkPlanOverflow {
            phase: "Chebyshev flux evaluation work",
        }),
        Some(BcValue::Signal(TimeSignal::Constant(_) | TimeSignal::Ramp { .. }))
        | Some(BcValue::Uniform(_) | BcValue::Profile(_))
        | None => Ok(4),
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

fn dedup_validation_times<E>(
    times: &mut Vec<f64>,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
) -> Result<(), E> {
    let mut write = 0usize;
    for read in 0..times.len() {
        checkpoint("net-flux checkpoint deduplication")?;
        let value = if times[read] == 0.0 { 0.0 } else { times[read] };
        if write == 0 || value != times[write - 1] {
            times[write] = value;
            write += 1;
        }
    }
    times.truncate(write);
    Ok(())
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
        let mut max_frame_name_bytes = 0usize;
        let mut max_case_name_bytes = 0usize;
        let mut max_combination_name_bytes = 0usize;
        let mut max_term_reference_bytes = 0usize;
        let mut max_ensemble_name_bytes = 0usize;
        let mut max_contact_pair_bytes = 0usize;
        let mut signal_scalars = 0usize;
        for frame in &self.frames.frames {
            max_frame_name_bytes = max_frame_name_bytes.max(frame.name.len().max(1));
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
        let mut base_flux_provider_work = 0u128;
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
            checked_work_add(
                &mut base_flux_provider_work,
                bc_flux_evaluation_work(bc)?,
                "base flux provider work",
            )?;
            checkpoint("preflight base boundary conditions")?;
        }
        for case in &self.cases {
            max_case_name_bytes = max_case_name_bytes.max(case.name.len().max(1));
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
            max_combination_name_bytes =
                max_combination_name_bytes.max(combination.name.len().max(1));
            checked_count_add(
                &mut identity_bytes,
                combination.name.len(),
                "combination identity bytes",
            )?;
            for (case, _) in &combination.terms {
                max_term_reference_bytes = max_term_reference_bytes.max(case.len().max(1));
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
            max_ensemble_name_bytes = max_ensemble_name_bytes.max(ensemble.name.len().max(1));
            checked_count_add(
                &mut identity_bytes,
                ensemble.name.len(),
                "ensemble identity bytes",
            )?;
            checkpoint("preflight ensembles")?;
        }
        for contact in &self.contacts {
            let pair_bytes = contact
                .region_a
                .len()
                .checked_add(contact.region_b.len())
                .and_then(|bytes| bytes.checked_add(2))
                .ok_or(ValidationError::WorkPlanOverflow {
                    phase: "contact pair key width",
                })?;
            max_contact_pair_bytes = max_contact_pair_bytes.max(pair_bytes);
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
        let mut flux_sort_work = heap_sort_work(flux_checkpoints, "base flux sort work")?;
        let mut flux_evaluation_work = checked_work_scale(
            base_flux_provider_work,
            flux_checkpoints,
            "base flux evaluation work",
        )?;
        let mut flux_set_scan_work = checked_work_product(base_bcs, 3, "base flux set scan work")?;
        for case in &self.cases {
            let mut case_contribution = 0usize;
            let mut case_provider_work = 0u128;
            for bc in &case.bcs {
                checked_count_add(
                    &mut case_contribution,
                    bc_flux_checkpoints(bc),
                    "case flux checkpoints",
                )?;
                checked_work_add(
                    &mut case_provider_work,
                    bc_flux_evaluation_work(bc)?,
                    "case flux provider work",
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
            let sort = heap_sort_work(set_checkpoints, "case flux sort work")?;
            checked_work_add(&mut flux_sort_work, sort, "total flux sort work")?;
            let set_bcs =
                base_bcs
                    .checked_add(case.bcs.len())
                    .ok_or(ValidationError::WorkPlanOverflow {
                        phase: "effective-set boundary count",
                    })?;
            checked_work_add(
                &mut flux_set_scan_work,
                checked_work_product(set_bcs, 3, "case flux set scan work")?,
                "total flux set scan work",
            )?;
            let set_provider_work = base_flux_provider_work
                .checked_add(case_provider_work)
                .ok_or(ValidationError::WorkPlanOverflow {
                    phase: "effective-set flux provider work",
                })?;
            checked_work_add(
                &mut flux_evaluation_work,
                checked_work_scale(
                    set_provider_work,
                    set_checkpoints,
                    "case flux evaluation work",
                )?,
                "total flux evaluation work",
            )?;
            checkpoint("preflight case flux sets")?;
        }
        enforce_limit(
            "flux checkpoints",
            flux_checkpoints,
            budget.max_flux_checkpoints,
        )?;
        let flux_checkpoint_work = checked_work_product(
            flux_checkpoints,
            2,
            "flux checkpoint materialization and deduplication work",
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
        let mut identity_work = work_units(identity_bytes, "identity linear scan work")?;
        for (items, lookups, max_key_bytes, phase) in [
            (
                frames,
                frames,
                max_frame_name_bytes,
                "frame name comparison work",
            ),
            (
                cases,
                cases,
                max_case_name_bytes,
                "case name comparison work",
            ),
            (
                combinations,
                combinations,
                max_combination_name_bytes,
                "combination name comparison work",
            ),
            (
                combination_terms,
                combination_terms,
                max_term_reference_bytes,
                "combination term comparison work",
            ),
            (
                ensembles,
                ensembles,
                max_ensemble_name_bytes,
                "ensemble name comparison work",
            ),
            (
                contacts,
                contacts,
                max_contact_pair_bytes,
                "contact pair comparison work",
            ),
        ] {
            checked_work_add(
                &mut identity_work,
                string_index_identity_work(items, lookups, max_key_bytes, index_height, phase)?,
                phase,
            )?;
        }
        let cross_case_key_bytes = max_case_name_bytes
            .max(1)
            .checked_add(max_term_reference_bytes.max(1))
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "combination case-reference key width",
            })?;
        let cross_case_lookup_work = ordered_lookup_work(
            combination_terms,
            index_height,
            "combination case-reference lookup work",
        )?
        .checked_mul(work_units(
            cross_case_key_bytes,
            "combination case-reference lookup work",
        )?)
        .ok_or(ValidationError::WorkPlanOverflow {
            phase: "combination case-reference lookup work",
        })?;
        checked_work_add(
            &mut identity_work,
            cross_case_lookup_work,
            "combination case-reference lookup work",
        )?;
        checked_work_add(
            &mut identity_work,
            repeated_identity_comparison_work(
                contacts,
                5,
                max_contact_pair_bytes,
                "contact direct comparison work",
            )?,
            "contact direct comparison work",
        )?;
        let frame_lookup_count = frames
            .checked_mul(3)
            .and_then(|lookups| lookups.checked_add(total_bcs))
            .ok_or(ValidationError::WorkPlanOverflow {
                phase: "frame lookup count",
            })?;
        let numeric_lookup_work =
            ordered_lookup_work(frame_lookup_count, index_height, "frame lookup work")?;
        let mut planned_work = work_units(record_visits, "validation record work")?;
        for (amount, phase) in [
            (index_work, "validation index work"),
            (identity_work, "identity comparison work"),
            (numeric_lookup_work, "frame lookup work"),
            (
                work_units(signal_scalars, "signal scalar work")?,
                "signal scalar work",
            ),
            (
                checked_work_product(contacts, 3, "contact classification work")?,
                "contact classification work",
            ),
            (flux_set_scan_work, "flux set scan work"),
            (flux_checkpoint_work, "flux checkpoint work"),
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
            let combo_name = DiagnosticIdentity(combo.name.as_str());
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
                        "combination {combo_name:?} first appears at row {first} and repeats at row {combo_index}"
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
                let case_name = DiagnosticIdentity(case.as_str());
                if case.is_empty() {
                    out.push(Violation {
                        code: "combo-case-empty",
                        what: format!(
                            "combination row {combo_index} {combo_name:?} term {term_index} has an empty case reference"
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
                            "combination row {combo_index} {combo_name:?} references case {case_name:?} at terms {first} and {term_index}"
                        ),
                        fix: "combine repeated case factors into one term".to_string(),
                    });
                }
                if !first_case_by_name.contains(&case.as_str()) {
                    out.push(Violation {
                        code: "combo-case-missing",
                        what: format!(
                            "combination row {combo_index} {combo_name:?} term {term_index} references unknown case {case_name:?}"
                        ),
                        fix: "reference only defined load cases".to_string(),
                    });
                }
                if !factor.is_finite() {
                    out.push(Violation {
                        code: "combo-factor",
                        what: format!(
                            "combination row {combo_index} {combo_name:?} term {term_index} has non-finite factor"
                        ),
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
        for _ in 0..contact_count {
            checkpoint("contact conflict initialization")?;
            conflicts.push(false);
        }

        let mut start = 0usize;
        while start < index.entries.len() {
            // A collection of all-unique pairs skips every duplicate-only
            // loop below, so the outer group boundary must poll on its own.
            checkpoint("contact conflict groups")?;
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
        let mut declares_incompressible = false;
        let mut has_pressure_outlet = false;
        for bc in self.base_bcs.iter().chain(case_bcs.iter()) {
            checkpoint("net-flux set classification")?;
            declares_incompressible |= bc.compatibility == Some(Compat::Incompressible);
            has_pressure_outlet |=
                bc.physics == Physics::IncompressibleFlow && bc.kind == BcKind::PressureOutlet;
        }
        if !declares_incompressible {
            return Ok(());
        }
        let label = case.map_or_else(|| "base".to_string(), |case| format!("base+{}", case.name));
        let mut requested_times = 1usize;
        for bc in self.base_bcs.iter().chain(case_bcs.iter()) {
            checkpoint("net-flux checkpoint count")?;
            requested_times = requested_times.checked_add(bc_flux_checkpoints(bc)).ok_or(
                ValidationError::WorkPlanOverflow {
                    phase: "net-flux checkpoint materialization",
                },
            )?;
        }
        let mut validation_times = Vec::new();
        reserve_validation_times(&mut validation_times, requested_times)?;
        checkpoint("net-flux checkpoint materialization")?;
        validation_times.push(0.0);
        for bc in self.base_bcs.iter().chain(case_bcs.iter()) {
            checkpoint("net-flux provider materialization")?;
            if bc.kind == BcKind::MassFlowInlet
                && let Some(BcValue::Signal(signal)) = &bc.value
            {
                signal.append_net_flux_validation_times(&mut validation_times, checkpoint)?;
            }
        }
        debug_assert_eq!(validation_times.len(), requested_times);
        sort_validation_by(
            &mut validation_times,
            "net-flux checkpoint sort",
            |a, b| a.total_cmp(b),
            checkpoint,
        )?;
        dedup_validation_times(&mut validation_times, checkpoint)?;
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
                let evaluated = bc.mass_flow_at_prevalidated_with_checkpoint(time, &mut || {
                    checkpoint("net-flux Chebyshev evaluation")
                })?;
                match evaluated {
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
        Combination, ContactLaw, ContactModel, DIAGNOSTIC_IDENTITY_PREVIEW_BYTES,
        DiagnosticIdentity, Environment, LoadCase, Scenario, ValidationBudget, ValidationError,
        ValidationIndex, ceil_log2, contact_pair, dedup_validation_times, heap_sort_work,
        ordered_lookup_work, reserve_validation_times, sort_validation_index,
        string_index_identity_work,
    };
    use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
    use crate::ensemble::{SpectrumModel, StochasticEnsemble};
    use crate::frame::{Frame, FrameId, FrameMotion};
    use crate::signal::{ChebProfile, Interp, TimeSignal};
    use fs_cheb::Cheb1;
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

        let shared_identity = "x".repeat(256);
        let mut identical = (0..64usize)
            .rev()
            .map(|row| (shared_identity.as_str(), row))
            .collect::<Vec<_>>();
        sort_validation_index(&mut identical, "identical identity sort", &mut |_| {
            Ok::<(), core::convert::Infallible>(())
        })
        .expect("infallible checkpoint");
        assert!(identical.iter().map(|(_, row)| *row).eq(0..identical.len()));
        let height = ceil_log2(identical.len() + 1);
        let charged = string_index_identity_work(
            identical.len(),
            identical.len(),
            shared_identity.len(),
            height,
            "identical identity work",
        )
        .expect("representable identity work");
        let width = u128::try_from(shared_identity.len()).expect("identity width fits u128");
        let sort_work = heap_sort_work(identical.len(), "identical identity work")
            .expect("representable sort work")
            * 2
            * width;
        let lookup_work = ordered_lookup_work(identical.len(), height, "identical identity work")
            .expect("representable lookup work")
            * 2
            * width;
        assert_eq!(charged, sort_work + lookup_work);

        let empty_key_work = string_index_identity_work(
            identical.len(),
            identical.len() * 2,
            0,
            height,
            "empty identity work",
        )
        .expect("empty identities retain comparison work");
        assert!(empty_key_work > 0);
        assert!(matches!(
            string_index_identity_work(
                usize::MAX,
                usize::MAX,
                usize::MAX,
                u128::MAX,
                "identity comparison overflow",
            ),
            Err(ValidationError::WorkPlanOverflow {
                phase: "identity comparison overflow",
            })
        ));
    }

    #[test]
    fn checkpoint_dedup_is_cancellable_and_canonicalizes_signed_zero() {
        let mut values = vec![-0.0, 0.0, 1.0, 1.0];
        let mut polls = 0usize;
        dedup_validation_times(&mut values, &mut |phase| {
            assert_eq!(phase, "net-flux checkpoint deduplication");
            polls += 1;
            Ok::<(), core::convert::Infallible>(())
        })
        .expect("infallible checkpoint");
        assert_eq!(polls, 4);
        assert_eq!(values, [0.0, 1.0]);
        assert_eq!(values[0].to_bits(), 0.0f64.to_bits());

        let mut cancelled = vec![0.0, 1.0, 2.0];
        let mut calls = 0usize;
        let result = dedup_validation_times(&mut cancelled, &mut |_| {
            calls += 1;
            if calls == 2 { Err("cancelled") } else { Ok(()) }
        });
        assert_eq!(result, Err("cancelled"));
    }

    #[test]
    fn chebyshev_net_flux_evaluation_polls_inside_the_recurrence() {
        let mut scenario = Scenario::new("cheb-cancellation", 7, Environment::earth_lab());
        scenario.base_bcs.push(BoundaryCondition {
            region: "cheb-inlet".to_string(),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
            value: Some(BcValue::Signal(TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(0.0, 1.0, vec![2.0, 0.5, 0.25, 0.125]),
                dims: Dims([0, 1, -1, 0, 0, 0]),
            }))),
            compatibility: Some(Compat::Incompressible),
            frame: 0,
        });
        let plan = scenario
            .validation_plan(ValidationBudget::default())
            .expect("Chebyshev scenario plan");

        let mut recurrence_polls = 0usize;
        let result = scenario.validate_admitted(plan.finding_capacity, &mut |phase| {
            if phase == "net-flux Chebyshev evaluation" {
                recurrence_polls += 1;
                if recurrence_polls == 2 {
                    return Err(ValidationError::Cancelled {
                        phase,
                        completed: 0,
                        planned: plan.planned_work,
                    });
                }
            }
            Ok(())
        });
        assert_eq!(recurrence_polls, 2);
        assert!(matches!(
            result,
            Err(ValidationError::Cancelled {
                phase: "net-flux Chebyshev evaluation",
                completed: 0,
                planned,
            }) if planned == plan.planned_work
        ));
    }

    #[test]
    fn repeated_combination_diagnostics_bound_parent_identity_rendering() {
        assert_eq!(format!("{:?}", DiagnosticIdentity("short")), "\"short\"");
        let unicode = "界".repeat(DIAGNOSTIC_IDENTITY_PREVIEW_BYTES);
        let unicode_preview = format!("{:?}", DiagnosticIdentity(unicode.as_str()));
        let unicode_prefix_end =
            DIAGNOSTIC_IDENTITY_PREVIEW_BYTES - (DIAGNOSTIC_IDENTITY_PREVIEW_BYTES % "界".len());
        assert_eq!(
            unicode_preview,
            format!(
                "{:?}…<{} bytes total>",
                &unicode[..unicode_prefix_end],
                unicode.len()
            )
        );

        let combo_name = "x".repeat(DIAGNOSTIC_IDENTITY_PREVIEW_BYTES * 32);
        let case_prefix = "y".repeat(DIAGNOSTIC_IDENTITY_PREVIEW_BYTES);
        let term_count = 64usize;
        let mut scenario = Scenario::new("bounded-diagnostics", 1, Environment::earth_lab());
        scenario.combinations.push(Combination {
            name: combo_name.clone(),
            terms: (0..term_count)
                .map(|term| (format!("{case_prefix}-missing-{term}"), f64::NAN))
                .collect(),
        });

        let findings = scenario.validate();
        let term_findings = findings
            .iter()
            .filter(|finding| matches!(finding.code, "combo-case-missing" | "combo-factor"))
            .collect::<Vec<_>>();
        assert_eq!(term_findings.len(), term_count * 2);
        for finding in term_findings {
            assert!(finding.what.contains("combination row 0"));
            assert!(finding.what.contains("term "));
            assert!(!finding.what.contains(combo_name.as_str()));
            assert!(!finding.what.contains("-missing-"));
            assert!(
                finding.what.len() < DIAGNOSTIC_IDENTITY_PREVIEW_BYTES * 4 + 256,
                "unbounded diagnostic: {} bytes",
                finding.what.len()
            );
        }

        let shared_prefix = "z".repeat(DIAGNOSTIC_IDENTITY_PREVIEW_BYTES);
        let mut collision = Scenario::new("preview-collision", 2, Environment::earth_lab());
        collision.combinations.extend([
            Combination {
                name: format!("{shared_prefix}-alpha"),
                terms: vec![("missing".to_string(), f64::NAN)],
            },
            Combination {
                name: format!("{shared_prefix}-omega"),
                terms: vec![("missing".to_string(), f64::NAN)],
            },
        ]);
        let factor_findings = collision
            .validate()
            .into_iter()
            .filter(|finding| finding.code == "combo-factor")
            .map(|finding| finding.what)
            .collect::<Vec<_>>();
        assert_eq!(factor_findings.len(), 2);
        assert!(factor_findings[0].contains("combination row 0"));
        assert!(factor_findings[1].contains("combination row 1"));
        assert_ne!(factor_findings[0], factor_findings[1]);
    }

    #[test]
    fn unique_contact_groups_poll_cancellation_at_each_group() {
        let mut scenario = Scenario::new("unique-contacts", 3, Environment::earth_lab());
        scenario.contacts.extend([
            ContactLaw {
                region_a: "a".to_string(),
                region_b: "b".to_string(),
                model: ContactModel::Frictionless,
            },
            ContactLaw {
                region_a: "c".to_string(),
                region_b: "d".to_string(),
                model: ContactModel::Frictionless,
            },
            ContactLaw {
                region_a: "e".to_string(),
                region_b: "f".to_string(),
                model: ContactModel::Frictionless,
            },
        ]);
        let index = ValidationIndex::build(
            scenario
                .contacts
                .iter()
                .enumerate()
                .map(|(row, contact)| (contact_pair(contact), row)),
            "unique contact test index",
            &mut |_| Ok(()),
        )
        .expect("fallible contact index allocation");

        let mut initialized = 0usize;
        let mut groups = 0usize;
        let result = scenario.contact_pair_conflicts(&index, &mut |phase| {
            if phase == "contact conflict initialization" {
                initialized += 1;
                return Ok(());
            }
            assert_eq!(phase, "contact conflict groups");
            groups += 1;
            if groups == 2 {
                Err(ValidationError::Cancelled {
                    phase,
                    completed: 0,
                    planned: 0,
                })
            } else {
                Ok(())
            }
        });
        assert_eq!(initialized, scenario.contacts.len());
        assert_eq!(groups, 2);
        assert!(matches!(
            result,
            Err(ValidationError::Cancelled {
                phase: "contact conflict groups",
                completed: 0,
                planned: 0,
            })
        ));
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
            "frame cycle state initialization",
            "frame cycle result initialization",
            "frame cycle traversal",
            "signal table times",
            "base boundary conditions",
            "case boundary conditions",
            "combination terms",
            "ensembles",
            "contact conflict initialization",
            "contact conflict groups",
            "contact conflict models",
            "contacts",
            "net-flux set classification",
            "net-flux checkpoint count",
            "net-flux checkpoint materialization",
            "net-flux provider materialization",
            "net-flux checkpoint sort",
            "net-flux checkpoint deduplication",
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

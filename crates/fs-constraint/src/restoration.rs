//! Restoration-path resource contracts (bead
//! frankensim-constraint-restoration-budget-receipts-x5sev).
//!
//! The elastic solve and the infeasibility diagnosis are synchronous
//! small-fixture paths whose work scales with design dimension,
//! constraint count, skip-mask size, restart/step schedule, deletion
//! filtering, repair candidates, and Monte-Carlo feasibility samples.
//! This module gives those paths the same admitted-resource contract
//! the chance path earned in bead frankensim-oxyjg:
//!
//! - [`RestorationWorkPlan`]: a versioned, checked-arithmetic worst-case
//!   work declaration computed BEFORE allocation/evaluation. Overflow,
//!   counts above versioned caps, and inconsistent totals are typed
//!   refusals, never wraps or silent truncation.
//! - [`RestorationWorkReceipt`]: what actually ran — plan identity,
//!   final [`BudgetConsumption`] from the shared accountant, charged
//!   work units, completed starts, and the memory-authority boundary.
//! - [`RestorationError`]: `Invalid` for pre-admission input faults;
//!   `Refused` when ADMITTED work was stopped by cancellation,
//!   deadline, poll, or cost authority. A refusal always carries its
//!   receipt; no partial success is ever published.
//!
//! Memory authority is an honest boundary, not invented power: a `Cx`
//! assembled by hand carries no operation memory lease (`fs-alloc`), so
//! buffers proportional to dimensions/constraints/repairs are allocated
//! WITHOUT lease admission and the receipt reports
//! [`RestorationMemoryAuthority::NoLeaseNoClaim`]. Only a lease-carrying
//! context may support a stronger claim, and this crate does not mint
//! private leases to fake one.

use crate::ConError;
use fs_exec::{BudgetConsumption, BudgetRefusal};

/// Schema version of [`RestorationWorkPlan`] canonical bytes. Bump on any
/// field-semantics change; plan identity binds this version.
pub const RESTORATION_WORK_PLAN_SCHEMA_VERSION: u32 = 1;

/// Declared work units of one canonical skip-mask entry construction.
pub const RESTORATION_UNITS_SKIP_MASK_ENTRY: u64 = 1;

/// Declared work units of one scalar constraint evaluation (`scalar_at`).
pub const RESTORATION_UNITS_CONSTRAINT_EVALUATION: u64 = 1;

/// Versioned cap on multi-start count per descent schedule.
pub const RESTORATION_MAX_STARTS: u32 = 8;

/// Versioned cap on projected-subgradient steps per start.
pub const RESTORATION_MAX_STEPS_PER_START: u32 = 300;

/// Versioned cap on Monte-Carlo samples per repair feasibility estimate.
pub const RESTORATION_MAX_FEASIBILITY_SAMPLES: u32 = 400;

/// Caller-declared descent/sampling schedule, validated against versioned
/// caps. The [`Default`] values are exactly the historical hardcoded
/// schedule (8 starts, 300 steps, 400 samples), so default plans replay
/// legacy runs bit-for-bit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RestorationWorkLimits {
    /// Multi-start count for one elastic descent.
    pub starts: u32,
    /// Maximum subgradient steps per start.
    pub steps_per_start: u32,
    /// Monte-Carlo samples per repair feasibility estimate.
    pub feasibility_samples: u32,
}

impl Default for RestorationWorkLimits {
    fn default() -> Self {
        Self {
            starts: RESTORATION_MAX_STARTS,
            steps_per_start: RESTORATION_MAX_STEPS_PER_START,
            feasibility_samples: RESTORATION_MAX_FEASIBILITY_SAMPLES,
        }
    }
}

impl RestorationWorkLimits {
    fn checked(self) -> Result<Self, ConError> {
        if self.starts == 0 || self.starts > RESTORATION_MAX_STARTS {
            return Err(ConError::BadParam {
                what: "restoration limits.starts must be in 1..=RESTORATION_MAX_STARTS",
                value: f64::from(self.starts),
            });
        }
        if self.steps_per_start == 0 || self.steps_per_start > RESTORATION_MAX_STEPS_PER_START {
            return Err(ConError::BadParam {
                what: "restoration limits.steps_per_start must be in \
                       1..=RESTORATION_MAX_STEPS_PER_START",
                value: f64::from(self.steps_per_start),
            });
        }
        if self.feasibility_samples == 0
            || self.feasibility_samples > RESTORATION_MAX_FEASIBILITY_SAMPLES
        {
            return Err(ConError::BadParam {
                what: "restoration limits.feasibility_samples must be in \
                       1..=RESTORATION_MAX_FEASIBILITY_SAMPLES",
                value: f64::from(self.feasibility_samples),
            });
        }
        Ok(self)
    }

    fn mix_into(&self, hash: &mut u64) {
        for field in [
            self.starts.to_le_bytes(),
            self.steps_per_start.to_le_bytes(),
            self.feasibility_samples.to_le_bytes(),
        ] {
            for byte in field {
                *hash ^= u64::from(byte);
                *hash = hash.wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
}

/// The inputs a restoration work plan declares, before checked arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RestorationWorkShape {
    /// Admitted host dimension (exact `Rn` dim of the sole variable).
    pub dimensions: u32,
    /// Total constraint count of the evaluated set.
    pub constraints_total: u32,
    /// Count of validated unique in-range skip indices (the canonical
    /// admitted mask cardinality). Skip ORDER never enters the plan:
    /// equivalent orderings share one identity.
    pub skipped_count: u32,
    /// Descent/sampling schedule.
    pub limits: RestorationWorkLimits,
}

/// The admitted worst-case resource contract of a restoration run.
///
/// All aggregates are computed with checked integer arithmetic over the
/// fixed unit weights; a plan that cannot state its own cost refuses
/// instead of wrapping. Deletion-filter and repair allowances bound the
/// data-dependent diagnosis phases at their structural worst cases
/// (`2 * constraints_total + 2` subset solves; `3 * constraints_total`
/// repair estimates), so a diagnosis admits once, up front, for the
/// whole run.
///
/// The plan carries no authority by itself — the diagnose/elastic entry
/// points admit it against the caller's `Cx` budget before the first
/// charged tile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RestorationWorkPlan {
    /// Must equal [`RESTORATION_WORK_PLAN_SCHEMA_VERSION`].
    pub schema_version: u32,
    /// Declared schedule.
    pub limits: RestorationWorkLimits,
    /// Exact host dimension.
    pub dimensions: u32,
    /// Total constraint count.
    pub constraints_total: u32,
    /// Active (non-skipped) constraint count.
    pub active_constraints: u32,
    /// Checked evals in ONE total-violation pass over the active set:
    /// `active * RESTORATION_UNITS_CONSTRAINT_EVALUATION`.
    pub evals_per_total_pass: u64,
    /// Checked evals in ONE finite-difference step over the active set
    /// (both probes, all dimensions): `2 * dimensions * active`.
    pub fd_evals_per_step: u64,
    /// Checked evals in one full descent start: initial pass plus every
    /// capped step's probes and post-step re-evaluation:
    /// `active + steps * (fd + active)`.
    pub descent_evals_per_start: u64,
    /// Checked base-phase units outside the starts loop: canonical
    /// skip-mask build + best-seed pass + final violation pass.
    pub base_phase_units: u64,
    /// Checked worst-case deletion-filter allowance (diagnosis only;
    /// zero for a bare elastic-solve plan).
    pub filter_allowance_units: u64,
    /// Checked worst-case repair-estimate allowance (diagnosis only;
    /// zero for a bare elastic-solve plan).
    pub repair_allowance_units: u64,
    /// Checked grand total: base + starts + allowances.
    pub total_work_units: u64,
}

fn overflow(what: &'static str, value: f64) -> ConError {
    ConError::BadParam { what, value }
}

const fn u32f(value: u32) -> f64 {
    value as f64
}

impl RestorationWorkPlan {
    /// Build a checked plan over the declared shape. Zero dimensions,
    /// an impossible skip count, any count above its versioned cap, and
    /// any arithmetic overflow are typed refusals.
    ///
    /// # Errors
    /// [`ConError::BadParam`] with a teaching message for every boundary.
    pub fn plan(shape: RestorationWorkShape) -> Result<Self, ConError> {
        let limits = shape.limits.checked()?;
        if shape.dimensions == 0 {
            return Err(overflow(
                "restoration work plan needs a positive dimension",
                0.0,
            ));
        }
        if shape.skipped_count > shape.constraints_total {
            return Err(overflow(
                "restoration skipped_count cannot exceed constraints_total",
                u32f(shape.skipped_count),
            ));
        }
        let dims = u64::from(shape.dimensions);
        let ct = u64::from(shape.constraints_total);
        let active = u64::from(shape.constraints_total - shape.skipped_count);
        let steps = u64::from(limits.steps_per_start);
        let starts = u64::from(limits.starts);
        let unit = RESTORATION_UNITS_CONSTRAINT_EVALUATION;

        // One total-violation pass over the active set.
        let evals_per_total_pass = active.checked_mul(unit).ok_or_else(|| {
            overflow(
                "restoration plan overflows evals_per_total_pass",
                u32f(shape.active_constraints()),
            )
        })?;
        // Both finite-difference probes across every dimension.
        let fd_evals_per_step = dims
            .checked_mul(2)
            .and_then(|probes| probes.checked_mul(active))
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows fd_evals_per_step",
                    u32f(shape.dimensions),
                )
            })?;
        // Per step: probes + post-step re-evaluation pass.
        let evals_per_step = fd_evals_per_step
            .checked_add(evals_per_total_pass)
            .ok_or_else(|| {
                overflow("restoration plan overflows per-step work", u32f(shape.dimensions))
            })?;
        // Per start: initial pass + every capped step.
        let descent_evals_per_start = evals_per_step
            .checked_mul(steps)
            .and_then(|steps_work| steps_work.checked_add(evals_per_total_pass))
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows descent work per start",
                    u32f(limits.steps_per_start),
                )
            })?;
        let starts_work = descent_evals_per_start.checked_mul(starts).ok_or_else(|| {
            overflow(
                "restoration plan overflows starts work",
                u32f(limits.starts),
            )
        })?;
        // Base phase outside the starts loop: canonical skip-mask build
        // (every declared entry) + best-seed pass + final violation pass.
        let mask_units =
            ct.checked_mul(RESTORATION_UNITS_SKIP_MASK_ENTRY).ok_or_else(|| {
                overflow(
                    "restoration plan overflows skip-mask work",
                    u32f(shape.constraints_total),
                )
            })?;
        let base_phase_units = mask_units
            .checked_add(evals_per_total_pass)
            .and_then(|work| work.checked_add(evals_per_total_pass))
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows base-phase work",
                    u32f(shape.constraints_total),
                )
            })?;

        // Diagnosis-only allowances. A subset solve rebuilds its own mask
        // over ALL declared entries and descends with at most the full
        // count active, so its upper bound reuses the same formulas at
        // `active == constraints_total`. Deletion filtering performs at
        // most `2*ct + 2` subset solves (index advances + removals +
        // support verification + final verification); repairs perform at
        // most `3*ct` estimates of at most `samples * ct` evals each.
        let subset_fd_step = dims
            .checked_mul(2)
            .and_then(|probes| probes.checked_mul(ct))
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows subset-step allowance",
                    u32f(shape.dimensions),
                )
            })?;
        let subset_evals_per_step = subset_fd_step.checked_add(ct).ok_or_else(|| {
            overflow(
                "restoration plan overflows subset-step allowance",
                u32f(shape.dimensions),
            )
        })?;
        let subset_descent = subset_evals_per_step
            .checked_mul(steps)
            .and_then(|steps_work| steps_work.checked_add(ct))
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows subset-run allowance",
                    u32f(limits.steps_per_start),
                )
            })?;
        let subset_run_upper = mask_units
            .checked_add(ct)
            .and_then(|head| head.checked_add(starts.checked_mul(subset_descent)?))
            .and_then(|body| body.checked_add(ct))
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows subset-run allowance",
                    u32f(shape.constraints_total),
                )
            })?;
        let filter_runs = ct.checked_mul(2).and_then(|runs| runs.checked_add(2)).ok_or_else(
            || {
                overflow(
                    "restoration plan overflows filter-run count",
                    u32f(shape.constraints_total),
                )
            },
        )?;
        let filter_allowance_units = filter_runs.checked_mul(subset_run_upper).ok_or_else(
            || {
                overflow(
                    "restoration plan overflows filter allowance",
                    u32f(shape.constraints_total),
                )
            },
        )?;
        let estimates = ct.checked_mul(3).ok_or_else(|| {
            overflow(
                "restoration plan overflows repair-estimate count",
                u32f(shape.constraints_total),
            )
        })?;
        let estimate_upper = u64::from(limits.feasibility_samples)
            .checked_mul(ct)
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows repair-estimate work",
                    u32f(limits.feasibility_samples),
                )
            })?;
        let repair_allowance_units = estimates
            .checked_mul(estimate_upper)
            .ok_or_else(|| {
                overflow(
                    "restoration plan overflows repair allowance",
                    u32f(shape.constraints_total),
                )
            })?;

        let total_work_units = base_phase_units
            .checked_add(starts_work)
            .and_then(|work| work.checked_add(filter_allowance_units))
            .and_then(|work| work.checked_add(repair_allowance_units))
            .ok_or_else(|| overflow("restoration plan overflows total work", u32f(u32::MAX)))?;

        Ok(Self {
            schema_version: RESTORATION_WORK_PLAN_SCHEMA_VERSION,
            limits,
            dimensions: shape.dimensions,
            constraints_total: shape.constraints_total,
            active_constraints: shape.constraints_total - shape.skipped_count,
            evals_per_total_pass,
            fd_evals_per_step,
            descent_evals_per_start,
            base_phase_units,
            filter_allowance_units,
            repair_allowance_units,
            total_work_units,
        })
    }

    /// Recompute every aggregate from the stored fields and compare.
    ///
    /// # Errors
    /// [`ConError::BadParam`] naming the first mismatched field.
    pub fn verify_consistency(&self) -> Result<(), ConError> {
        let rebuilt = Self::plan(RestorationWorkShape {
            dimensions: self.dimensions,
            constraints_total: self.constraints_total,
            skipped_count: self.constraints_total - self.active_constraints,
            limits: self.limits,
        })?;
        if *self != rebuilt {
            return Err(overflow(
                "restoration work plan fields disagree with their checked aggregates",
                u32f(self.schema_version),
            ));
        }
        Ok(())
    }

    /// Deterministic FNV-1a identity over the canonical little-endian
    /// field order. Two plans declaring the same workload share one
    /// identity regardless of skip ordering; receipts bind it so retained
    /// evidence names the exact contract.
    #[must_use]
    pub fn identity(&self) -> u64 {
        fn mix(mut hash: u64, bytes: &[u8]) -> u64 {
            for byte in bytes {
                hash ^= u64::from(*byte);
                hash = hash.wrapping_mul(0x0100_0000_01b3);
            }
            hash
        }
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        hash = mix(hash, &self.schema_version.to_le_bytes());
        self.limits.mix_into(&mut hash);
        hash = mix(hash, &self.dimensions.to_le_bytes());
        hash = mix(hash, &self.constraints_total.to_le_bytes());
        hash = mix(hash, &self.active_constraints.to_le_bytes());
        hash = mix(hash, &self.evals_per_total_pass.to_le_bytes());
        hash = mix(hash, &self.fd_evals_per_step.to_le_bytes());
        hash = mix(hash, &self.descent_evals_per_start.to_le_bytes());
        hash = mix(hash, &self.base_phase_units.to_le_bytes());
        hash = mix(hash, &self.filter_allowance_units.to_le_bytes());
        hash = mix(hash, &self.repair_allowance_units.to_le_bytes());
        mix(hash, &self.total_work_units.to_le_bytes())
    }
}

impl RestorationWorkShape {
    const fn active_constraints(&self) -> u32 {
        self.constraints_total.saturating_sub(self.skipped_count)
    }
}

/// Honest memory-authority boundary recorded on every receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RestorationMemoryAuthority {
    /// The context carried an operation memory lease; buffers
    /// proportional to dimensions/constraints/repairs MAY be claimed as
    /// lease-admitted by future slices that allocate through it.
    LeaseAdmitted,
    /// The context carried no lease. Allocation happened WITHOUT memory
    /// authority; the report makes NO memory-boundedness claim. This
    /// crate never mints a private lease to upgrade itself.
    NoLeaseNoClaim,
}

/// Retained receipt of one restoration run: which plan ran, how much of
/// the shared budget it consumed, how much work completed, and under
/// which memory authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RestorationWorkReceipt {
    /// Identity of the executed [`RestorationWorkPlan`].
    pub plan_identity: u64,
    /// Schema version of the executed plan.
    pub schema_version: u32,
    /// Final accountant state; `None` only when admission itself was
    /// refused, so no budget contract ever existed to report.
    pub consumption: Option<BudgetConsumption>,
    /// Work units charged as complete (mask entries + evaluations).
    pub work_units_charged: u64,
    /// Descent starts whose full schedule completed.
    pub starts_completed: u32,
    /// Memory-authority boundary observed for this run.
    pub memory: RestorationMemoryAuthority,
}

impl RestorationWorkReceipt {
    pub(crate) fn refused_admission(plan: &RestorationWorkPlan) -> Self {
        Self {
            plan_identity: plan.identity(),
            schema_version: plan.schema_version,
            consumption: None,
            work_units_charged: 0,
            starts_completed: 0,
            memory: RestorationMemoryAuthority::NoLeaseNoClaim,
        }
    }
}

/// Typed outcome of a planned restoration run. `Invalid` carries the
/// ordinary constraint-calculus refusals raised before admission;
/// `Refused` means work WAS admitted and then stopped by the budget
/// authority — the receipt proves exactly how far it got, and NO report
/// is produced (a partial solve never mints feasibility evidence).
#[derive(Debug, Clone, PartialEq)]
pub enum RestorationError {
    /// Typed refusal before any work was admitted.
    Invalid(ConError),
    /// Admitted work stopped by cancellation, deadline, poll, or cost
    /// authority. The receipt retains the exact stop point.
    Refused {
        /// The latched accountant refusal.
        refusal: BudgetRefusal,
        /// Work completed and charged up to the refusal.
        receipt: RestorationWorkReceipt,
    },
}

impl core::fmt::Display for RestorationError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(error) => write!(formatter, "restoration refused: {error}"),
            Self::Refused { refusal, receipt } => write!(
                formatter,
                "restoration stopped after {} starts ({} work units charged): {refusal}",
                receipt.starts_completed, receipt.work_units_charged
            ),
        }
    }
}

impl core::error::Error for RestorationError {}

//! Budget-first admission for exact component-by-component lattice construction.
//!
//! The live exact constructor remains [`crate::qmc::Lattice::cbc`]. This module
//! is the allocation-free prerequisite for moving that implementation behind a
//! typed execution boundary: it validates the count domain, computes checked
//! conservative exact-integer limb/scalar-work/logical-state bounds, and
//! admits only budgets that cover those bounds on the current target. It
//! performs no CBC candidate evaluation and does not allocate.
//!
//! This tranche does not yet claim that `Lattice::cbc` consumes the receipt,
//! supports cancellation or pause/resume, produces a minimality certificate,
//! or has a compact independent checker. The estimator deliberately
//! overcharges all visits at the largest relevant limb widths. Its memory
//! quantity is the maximum of the candidate and product-update phases under
//! the documented exact-capacity state model, including the old/new product
//! overlap created by `ExactNat::mul_assign_factor`. Allocator metadata,
//! allocator growth/rounding, complete call stacks, and process RSS are
//! explicitly outside this tranche's claim. The execution tranche must use
//! flat/exact-capacity storage or charge observed allocator capacity. Work is
//! likewise a deterministic logical debit schedule, not a bound on allocator
//! internals, instructions, energy, or elapsed time.

/// Version of the CBC admission and resource-estimate semantics.
pub const CBC_ADMISSION_SCHEMA_VERSION: u32 = 2;

const LIMB_BYTES: u128 = 4;
const FACTOR_SCRATCH_BYTES: u128 = 4 * LIMB_BYTES;
// Scalar work units are source-level logical charges, not CPU instructions or
// elapsed-time predictions. These constants are part of schema v2. A visit
// debits six residue primitives, ten kernel-numerator primitives, and four
// loop/index primitives. A factor debits one fixed termination primitive plus
// eight decomposition primitives per limb; a Euclidean step debits
// branch/modulo/transfer. Candidate and dimension charges cover their fixed
// loop, selection, and append control. An executor claiming this receipt must
// debit the same schedule or a proven refinement.
const SCALAR_UNITS_PER_LATTICE_VISIT: u128 = 20;
const SCALAR_UNITS_PER_FACTOR: u128 = 1;
const SCALAR_UNITS_PER_FACTOR_LIMB: u128 = 8;
const SCALAR_UNITS_PER_GCD_STEP: u128 = 3;
const SCALAR_UNITS_PER_CANDIDATE: u128 = 16;
const SCALAR_UNITS_PER_DIMENSION: u128 = 8;

/// Exact CBC problem counts. Both values are dimensionless integer counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcProblem {
    point_count: u32,
    dimension: usize,
}

impl CbcProblem {
    /// Validate the structural count domain without allocating.
    ///
    /// # Errors
    /// [`CbcAdmissionError::InvalidPointCount`] when `point_count < 3`, or
    /// [`CbcAdmissionError::InvalidDimension`] when `dimension == 0`.
    #[must_use]
    pub const fn new(point_count: u32, dimension: usize) -> Result<Self, CbcAdmissionError> {
        if point_count < 3 {
            return Err(CbcAdmissionError::InvalidPointCount { point_count });
        }
        if dimension == 0 {
            return Err(CbcAdmissionError::InvalidDimension { dimension });
        }
        Ok(Self {
            point_count,
            dimension,
        })
    }

    /// Number of lattice points `n`.
    #[must_use]
    pub const fn point_count(self) -> u32 {
        self.point_count
    }

    /// Number of generator components.
    #[must_use]
    pub const fn dimension(self) -> usize {
        self.dimension
    }

    /// Compute conservative resource bounds with checked arithmetic and no
    /// allocation or CBC work.
    ///
    /// `lattice_visits` bounds the first-component product update plus every
    /// later candidate/point and chosen-product visit:
    ///
    /// `n + (dimension - 1) * n * (n - 1) + (dimension - 1) * n`.
    ///
    /// `limb_work_units` charges multiplication, worst-case carry propagation,
    /// zero-fill, normalization, and score comparison. `scalar_work_units`
    /// separately charges factor decomposition, residue/kernel evaluation,
    /// Euclidean GCD filtering, candidate control, and dimension control.
    /// Their sum is the budgeted `work_units`. A limb unit is one declared
    /// limb-loop visit; a scalar unit is one declared source-level primitive.
    /// Neither is a claim about CPU instructions.
    ///
    /// # Errors
    /// [`CbcAdmissionError::EstimateOverflow`] naming the first bound that
    /// leaves the `u128` accounting domain, or
    /// [`CbcAdmissionError::TargetCapacityExceeded`] when a modeled `Vec`
    /// length or byte allocation cannot be represented on this target.
    #[must_use]
    #[allow(clippy::too_many_lines)] // One checked derivation keeps the envelope auditable.
    pub fn estimate(self) -> Result<CbcEstimate, CbcAdmissionError> {
        let points = u128::from(self.point_count);
        let dimension =
            u128::try_from(self.dimension).map_err(|_| overflow("dimension conversion"))?;
        let later_dimensions = dimension
            .checked_sub(1)
            .ok_or_else(|| overflow("later dimensions"))?;
        let candidate_upper_bound = self.point_count - 1;
        let candidates = u128::from(candidate_upper_bound);

        // `1 + B2(r/n)` has numerator
        // `7*n^2 + 6*r^2 - 6*r*n <= 7*n^2` over `6*n^2`.
        let point_square = points
            .checked_mul(points)
            .ok_or_else(|| overflow("point-count square"))?;
        let kernel_numerator_upper = 7_u128
            .checked_mul(point_square)
            .ok_or_else(|| overflow("kernel numerator"))?;
        let kernel_numerator_bits = bit_width(kernel_numerator_upper);
        let kernel_factor_limbs = limbs_for_bits(u128::from(kernel_numerator_bits))?;

        let kernel_bits = u128::from(kernel_numerator_bits);
        let max_product_bits = kernel_bits
            .checked_mul(dimension)
            .ok_or_else(|| overflow("product bit width"))?;
        let max_product_limbs = limbs_for_bits(max_product_bits)?;
        // Every multiply/add source is either the initial ExactNat::one or a
        // normalized product from the preceding prefix. The final product is
        // a result, not a source; using d-1 here makes the carry derivation both
        // conservative and semantically exact about which limbs can trigger it.
        let max_source_product_bits = kernel_bits
            .checked_mul(later_dimensions)
            .ok_or_else(|| overflow("source-product bit width"))?;
        let max_source_product_limbs = limbs_for_bits(max_source_product_bits)?.max(1);
        let score_growth_bits = u128::from(ceil_log2(self.point_count));
        let max_score_bits = max_product_bits
            .checked_add(score_growth_bits)
            .ok_or_else(|| overflow("score bit width"))?;
        let max_score_limbs = limbs_for_bits(max_score_bits)?;
        // add_mul_factor first requests source + factor + one spare limb. A
        // multi-term score may then grow to its normalized score-width bound.
        let product_capacity_limbs = max_source_product_limbs
            .checked_add(kernel_factor_limbs)
            .and_then(|limbs| limbs.checked_add(1))
            .ok_or_else(|| overflow("product limb capacity"))?;
        let score_capacity_limbs = product_capacity_limbs.max(max_score_limbs);

        // Immediately before the final update, retained products were created
        // by the preceding update. Dimension one instead starts from `one`.
        let previous_product_capacity_limbs = if later_dimensions == 0 {
            1
        } else {
            let previous_source_prefix = later_dimensions
                .checked_sub(1)
                .ok_or_else(|| overflow("previous source prefix"))?;
            let previous_source_bits = kernel_bits
                .checked_mul(previous_source_prefix)
                .ok_or_else(|| overflow("previous source-product bit width"))?;
            limbs_for_bits(previous_source_bits)?
                .max(1)
                .checked_add(kernel_factor_limbs)
                .and_then(|limbs| limbs.checked_add(1))
                .ok_or_else(|| overflow("previous product limb capacity"))?
        };

        let candidate_count = later_dimensions
            .checked_mul(candidates)
            .ok_or_else(|| overflow("candidate count"))?;
        let candidate_visits = later_dimensions
            .checked_mul(points)
            .and_then(|visits| visits.checked_mul(candidates))
            .ok_or_else(|| overflow("candidate lattice visits"))?;
        let product_update_visits = dimension
            .checked_mul(points)
            .ok_or_else(|| overflow("product-update visits"))?;
        let lattice_visits = candidate_visits
            .checked_add(product_update_visits)
            .ok_or_else(|| overflow("total lattice visits"))?;
        let comparison_count = candidate_count;

        // Model the two actual high-water phases rather than adding mutually
        // exclusive buffers. Candidate evaluation retains preceding-prefix
        // products plus current and best scores. Product update has no score,
        // but mul_assign_factor simultaneously owns the moved old allocation
        // and the replacement allocation. Vec growth/allocator rounding is a
        // no-claim and must be removed by exact-capacity execution or metered.
        let vector_header_bytes = u128::try_from(core::mem::size_of::<Vec<u32>>())
            .map_err(|_| overflow("vector header conversion"))?;
        let best_owner_bytes = u128::try_from(core::mem::size_of::<Option<(Vec<u32>, u32)>>())
            .map_err(|_| overflow("best owner conversion"))?;
        let product_header_bytes = points
            .checked_mul(vector_header_bytes)
            .ok_or_else(|| overflow("product owner array"))?;
        let generator_payload_bytes = dimension
            .checked_mul(LIMB_BYTES)
            .ok_or_else(|| overflow("generator payload"))?;
        let previous_product_payload_bytes = points
            .checked_mul(previous_product_capacity_limbs)
            .and_then(|limbs| limbs.checked_mul(LIMB_BYTES))
            .ok_or_else(|| overflow("previous product payload"))?;
        let final_product_payload_bytes = points
            .checked_mul(product_capacity_limbs)
            .and_then(|limbs| limbs.checked_mul(LIMB_BYTES))
            .ok_or_else(|| overflow("final product payload"))?;
        let old_product_overlap_bytes = previous_product_capacity_limbs
            .checked_mul(LIMB_BYTES)
            .ok_or_else(|| overflow("old product overlap"))?;
        let score_payload_bytes = score_capacity_limbs
            .checked_mul(LIMB_BYTES)
            .ok_or_else(|| overflow("score payload"))?;
        let common_owner_bytes = product_header_bytes
            .checked_add(
                2_u128
                    .checked_mul(vector_header_bytes)
                    .ok_or_else(|| overflow("outer and generator owners"))?,
            )
            .ok_or_else(|| overflow("common owners"))?;
        let common_state_bytes = common_owner_bytes
            .checked_add(generator_payload_bytes)
            .and_then(|bytes| bytes.checked_add(FACTOR_SCRATCH_BYTES))
            .ok_or_else(|| overflow("common logical state"))?;
        // Assignment can transiently retain both the old `best` owner and a
        // new Option owner containing the moved current score. Charge the
        // larger of steady evaluation (Vec + Option) and replacement
        // (Option + Option), including tuple padding on this target.
        let candidate_score_owner_bytes = vector_header_bytes
            .checked_add(best_owner_bytes)
            .ok_or_else(|| overflow("candidate score owners"))?
            .max(
                2_u128
                    .checked_mul(best_owner_bytes)
                    .ok_or_else(|| overflow("replacement score owners"))?,
            );
        let candidate_phase_bytes = if candidate_count == 0 {
            0
        } else {
            common_state_bytes
                .checked_add(previous_product_payload_bytes)
                .and_then(|bytes| bytes.checked_add(candidate_score_owner_bytes))
                .and_then(|bytes| bytes.checked_add(2_u128.checked_mul(score_payload_bytes)?))
                .ok_or_else(|| overflow("candidate-phase logical state"))?
        };
        let update_phase_bytes = common_state_bytes
            .checked_add(final_product_payload_bytes)
            .and_then(|bytes| bytes.checked_add(old_product_overlap_bytes))
            .and_then(|bytes| bytes.checked_add(vector_header_bytes))
            .ok_or_else(|| overflow("update-phase logical state"))?;
        let logical_state_bytes = candidate_phase_bytes.max(update_phase_bytes);

        // Vec APIs take usize element counts and reject allocations larger than
        // isize::MAX bytes. Check each allocation independently. Aggregate
        // logical state is intentionally not compared with isize::MAX: it is
        // not one allocation.
        let usize_limit =
            u128::try_from(usize::MAX).map_err(|_| overflow("target usize conversion"))?;
        let vec_byte_limit =
            u128::try_from(isize::MAX).map_err(|_| overflow("target isize conversion"))?;
        target_bound("point-count element count", points, usize_limit)?;
        target_bound(
            "product owner-array bytes",
            product_header_bytes,
            vec_byte_limit,
        )?;
        target_bound("generator element count", dimension, usize_limit)?;
        target_bound(
            "generator allocation bytes",
            generator_payload_bytes,
            vec_byte_limit,
        )?;
        target_bound(
            "product element capacity",
            product_capacity_limbs,
            usize_limit,
        )?;
        target_bound(
            "product allocation bytes",
            product_capacity_limbs
                .checked_mul(LIMB_BYTES)
                .ok_or_else(|| overflow("product allocation"))?,
            vec_byte_limit,
        )?;
        target_bound(
            "previous-product element capacity",
            previous_product_capacity_limbs,
            usize_limit,
        )?;
        target_bound(
            "previous-product allocation bytes",
            old_product_overlap_bytes,
            vec_byte_limit,
        )?;
        target_bound("score element capacity", score_capacity_limbs, usize_limit)?;
        target_bound(
            "score allocation bytes",
            score_payload_bytes,
            vec_byte_limit,
        )?;

        let multiply_add_units = lattice_visits
            .checked_mul(max_source_product_limbs)
            .and_then(|units| units.checked_mul(kernel_factor_limbs))
            .ok_or_else(|| overflow("multiply-add limb work"))?;
        // add_mul_factor may propagate carry after every source limb. Charge
        // the full accumulator capacity for each such propagation rather than
        // assuming one carry pass per lattice visit.
        let carry_units = lattice_visits
            .checked_mul(max_source_product_limbs)
            .and_then(|units| units.checked_mul(score_capacity_limbs))
            .ok_or_else(|| overflow("carry limb work"))?;
        let comparison_limb_units = comparison_count
            .checked_mul(max_score_limbs)
            .ok_or_else(|| overflow("comparison limb work"))?;
        let zero_fill_limb_units = candidate_count
            .checked_mul(score_capacity_limbs)
            .and_then(|units| {
                product_update_visits
                    .checked_mul(product_capacity_limbs)
                    .and_then(|product_units| units.checked_add(product_units))
            })
            .ok_or_else(|| overflow("zero-fill limb work"))?;
        let normalization_limb_units = candidate_count
            .checked_mul(score_capacity_limbs)
            .and_then(|units| {
                product_update_visits
                    .checked_mul(product_capacity_limbs)
                    .and_then(|product_units| units.checked_add(product_units))
            })
            .ok_or_else(|| overflow("normalization limb work"))?;
        let limb_work_units = multiply_add_units
            .checked_add(carry_units)
            .and_then(|units| units.checked_add(zero_fill_limb_units))
            .and_then(|units| units.checked_add(normalization_limb_units))
            .and_then(|units| units.checked_add(comparison_limb_units))
            .ok_or_else(|| overflow("total limb work"))?;

        let gcd_step_upper_bound = u128::from(ceil_log2(self.point_count))
            .checked_mul(2)
            .and_then(|steps| steps.checked_add(1))
            .ok_or_else(|| overflow("GCD step bound"))?;
        let visit_scalar_units = lattice_visits
            .checked_mul(SCALAR_UNITS_PER_LATTICE_VISIT)
            .ok_or_else(|| overflow("residue/kernel scalar work"))?;
        let factor_scalar_units = lattice_visits
            .checked_mul(kernel_factor_limbs)
            .and_then(|units| units.checked_mul(SCALAR_UNITS_PER_FACTOR_LIMB))
            .and_then(|units| {
                lattice_visits
                    .checked_mul(SCALAR_UNITS_PER_FACTOR)
                    .and_then(|fixed_units| units.checked_add(fixed_units))
            })
            .ok_or_else(|| overflow("factor-decomposition scalar work"))?;
        let gcd_scalar_units = candidate_count
            .checked_mul(gcd_step_upper_bound)
            .and_then(|units| units.checked_mul(SCALAR_UNITS_PER_GCD_STEP))
            .ok_or_else(|| overflow("GCD scalar work"))?;
        let candidate_scalar_units = candidate_count
            .checked_mul(SCALAR_UNITS_PER_CANDIDATE)
            .ok_or_else(|| overflow("candidate scalar work"))?;
        let dimension_scalar_units = dimension
            .checked_mul(SCALAR_UNITS_PER_DIMENSION)
            .and_then(|units| units.checked_add(points))
            .and_then(|units| units.checked_add(dimension))
            .ok_or_else(|| overflow("dimension/initialization scalar work"))?;
        let scalar_work_units = visit_scalar_units
            .checked_add(factor_scalar_units)
            .and_then(|units| units.checked_add(gcd_scalar_units))
            .and_then(|units| units.checked_add(candidate_scalar_units))
            .and_then(|units| units.checked_add(dimension_scalar_units))
            .ok_or_else(|| overflow("total scalar work"))?;
        let work_units = limb_work_units
            .checked_add(scalar_work_units)
            .ok_or_else(|| overflow("total work"))?;

        Ok(CbcEstimate {
            target_pointer_width_bits: usize::BITS,
            candidate_upper_bound,
            candidate_count,
            kernel_numerator_upper,
            kernel_numerator_bits,
            kernel_factor_limbs,
            max_product_bits,
            max_product_limbs,
            max_source_product_limbs,
            max_score_bits,
            max_score_limbs,
            score_capacity_limbs,
            product_capacity_limbs,
            previous_product_capacity_limbs,
            lattice_visits,
            product_update_visits,
            comparison_count,
            limb_work_units,
            scalar_work_units,
            work_units,
            candidate_phase_bytes,
            update_phase_bytes,
            logical_state_bytes,
        })
    }

    /// Admit this problem only when both explicit budgets cover the checked
    /// estimate. Structural/target/estimate refusal precedes work-budget
    /// refusal, which precedes logical-state-budget refusal.
    ///
    /// # Errors
    /// [`CbcAdmissionError::EstimateOverflow`],
    /// [`CbcAdmissionError::TargetCapacityExceeded`],
    /// [`CbcAdmissionError::WorkBudgetExceeded`], or
    /// [`CbcAdmissionError::MemoryBudgetExceeded`].
    #[must_use]
    pub fn admit(self, budget: CbcBudget) -> Result<CbcAdmission, CbcAdmissionError> {
        let estimate = self.estimate()?;
        if estimate.work_units > budget.max_work_units {
            return Err(CbcAdmissionError::WorkBudgetExceeded {
                required: estimate.work_units,
                available: budget.max_work_units,
            });
        }
        if estimate.logical_state_bytes > budget.max_memory_bytes {
            return Err(CbcAdmissionError::MemoryBudgetExceeded {
                required: estimate.logical_state_bytes,
                available: budget.max_memory_bytes,
            });
        }
        Ok(CbcAdmission {
            schema_version: CBC_ADMISSION_SCHEMA_VERSION,
            problem: self,
            budget,
            estimate,
        })
    }
}

/// Explicit CBC admission budgets. Work units are the sum of the declared limb
/// and scalar primitive charges; memory units are logical exact-capacity state
/// bytes under this module's documented phase model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcBudget {
    max_work_units: u128,
    max_memory_bytes: u128,
}

impl CbcBudget {
    /// An effectively unbounded accounting envelope.
    pub const UNBOUNDED: Self = Self {
        max_work_units: u128::MAX,
        max_memory_bytes: u128::MAX,
    };

    /// Construct an explicit work/memory budget.
    #[must_use]
    pub const fn new(max_work_units: u128, max_memory_bytes: u128) -> Self {
        Self {
            max_work_units,
            max_memory_bytes,
        }
    }

    /// Maximum conservative schema-v2 work units.
    #[must_use]
    pub const fn max_work_units(self) -> u128 {
        self.max_work_units
    }

    /// Maximum logical-state bytes admitted by this accounting model.
    #[must_use]
    pub const fn max_memory_bytes(self) -> u128 {
        self.max_memory_bytes
    }
}

/// Checked conservative resource bounds for one CBC problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcEstimate {
    target_pointer_width_bits: u32,
    candidate_upper_bound: u32,
    candidate_count: u128,
    kernel_numerator_upper: u128,
    kernel_numerator_bits: u32,
    kernel_factor_limbs: u128,
    max_product_bits: u128,
    max_product_limbs: u128,
    max_source_product_limbs: u128,
    max_score_bits: u128,
    max_score_limbs: u128,
    score_capacity_limbs: u128,
    product_capacity_limbs: u128,
    previous_product_capacity_limbs: u128,
    lattice_visits: u128,
    product_update_visits: u128,
    comparison_count: u128,
    limb_work_units: u128,
    scalar_work_units: u128,
    work_units: u128,
    candidate_phase_bytes: u128,
    update_phase_bytes: u128,
    logical_state_bytes: u128,
}

impl CbcEstimate {
    /// Target pointer width bound into this estimate and receipt.
    #[must_use]
    pub const fn target_pointer_width_bits(self) -> u32 {
        self.target_pointer_width_bits
    }

    /// Upper bound on candidates examined per later component (`n - 1`).
    #[must_use]
    pub const fn candidate_upper_bound(self) -> u32 {
        self.candidate_upper_bound
    }

    /// Candidate loop iterations charged across all later components.
    #[must_use]
    pub const fn candidate_count(self) -> u128 {
        self.candidate_count
    }

    /// Maximum exact kernel numerator, bounded by `7*n^2`.
    #[must_use]
    pub const fn kernel_numerator_upper(self) -> u128 {
        self.kernel_numerator_upper
    }

    /// Bits needed for the maximum kernel numerator.
    #[must_use]
    pub const fn kernel_numerator_bits(self) -> u32 {
        self.kernel_numerator_bits
    }

    /// Base-2^32 limbs needed for one kernel factor.
    #[must_use]
    pub const fn kernel_factor_limbs(self) -> u128 {
        self.kernel_factor_limbs
    }

    /// Maximum product bit bound at the final prefix.
    #[must_use]
    pub const fn max_product_bits(self) -> u128 {
        self.max_product_bits
    }

    /// Maximum normalized product limbs.
    #[must_use]
    pub const fn max_product_limbs(self) -> u128 {
        self.max_product_limbs
    }

    /// Maximum normalized limbs of any multiplicative source (`d - 1`
    /// factors, or the initial one-limb value).
    #[must_use]
    pub const fn max_source_product_limbs(self) -> u128 {
        self.max_source_product_limbs
    }

    /// Maximum candidate-score bit bound, including summation over `n` points.
    #[must_use]
    pub const fn max_score_bits(self) -> u128 {
        self.max_score_bits
    }

    /// Maximum normalized candidate-score limbs.
    #[must_use]
    pub const fn max_score_limbs(self) -> u128 {
        self.max_score_limbs
    }

    /// Conservative score capacity including factor/spare and summation growth.
    #[must_use]
    pub const fn score_capacity_limbs(self) -> u128 {
        self.score_capacity_limbs
    }

    /// Requested capacity retained per exact product in base-2^32 limbs.
    #[must_use]
    pub const fn product_capacity_limbs(self) -> u128 {
        self.product_capacity_limbs
    }

    /// Retained per-product capacity immediately before the final update.
    #[must_use]
    pub const fn previous_product_capacity_limbs(self) -> u128 {
        self.previous_product_capacity_limbs
    }

    /// Upper bound on point/candidate visits.
    #[must_use]
    pub const fn lattice_visits(self) -> u128 {
        self.lattice_visits
    }

    /// Initial plus chosen product-update visits (`n * d`).
    #[must_use]
    pub const fn product_update_visits(self) -> u128 {
        self.product_update_visits
    }

    /// Upper bound on later-component candidate comparisons.
    #[must_use]
    pub const fn comparison_count(self) -> u128 {
        self.comparison_count
    }

    /// Conservative limb-operation admission units.
    #[must_use]
    pub const fn limb_work_units(self) -> u128 {
        self.limb_work_units
    }

    /// Conservative charged scalar-primitive units.
    #[must_use]
    pub const fn scalar_work_units(self) -> u128 {
        self.scalar_work_units
    }

    /// Total budgeted work units (`limb_work_units + scalar_work_units`).
    #[must_use]
    pub const fn work_units(self) -> u128 {
        self.work_units
    }

    /// Candidate-evaluation phase logical-state envelope.
    #[must_use]
    pub const fn candidate_phase_bytes(self) -> u128 {
        self.candidate_phase_bytes
    }

    /// Chosen-product update phase logical-state envelope, including the old
    /// allocation overlapping its replacement.
    #[must_use]
    pub const fn update_phase_bytes(self) -> u128 {
        self.update_phase_bytes
    }

    /// Maximum logical state in bytes across the documented phases.
    ///
    /// This includes exact requested payloads and modeled owner bytes but
    /// excludes allocator metadata/growth/rounding, complete call stacks, and
    /// process-level RSS effects.
    #[must_use]
    pub const fn logical_state_bytes(self) -> u128 {
        self.logical_state_bytes
    }
}

/// Sealed evidence that one problem fits an explicit CBC resource envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcAdmission {
    schema_version: u32,
    problem: CbcProblem,
    budget: CbcBudget,
    estimate: CbcEstimate,
}

impl CbcAdmission {
    /// Admission schema version.
    #[must_use]
    pub const fn schema_version(self) -> u32 {
        self.schema_version
    }

    /// Admitted problem counts.
    #[must_use]
    pub const fn problem(self) -> CbcProblem {
        self.problem
    }

    /// Explicit admitted budget.
    #[must_use]
    pub const fn budget(self) -> CbcBudget {
        self.budget
    }

    /// Checked resource estimate covered by the budget.
    #[must_use]
    pub const fn estimate(self) -> CbcEstimate {
        self.estimate
    }
}

/// Typed CBC admission refusal. All payloads are fixed-size and allocation-free.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcAdmissionError {
    /// `n` is below the exact constructor's supported domain.
    InvalidPointCount {
        /// Supplied lattice-point count.
        point_count: u32,
    },
    /// The generator would contain no components.
    InvalidDimension {
        /// Supplied dimension.
        dimension: usize,
    },
    /// A conservative bound left the `u128` accounting domain.
    EstimateOverflow {
        /// Stable name of the first overflowing quantity.
        quantity: &'static str,
    },
    /// One modeled allocation length or byte capacity cannot be represented
    /// by this target's `usize`/`isize`-bounded `Vec` API.
    TargetCapacityExceeded {
        /// Stable name of the rejected allocation quantity.
        quantity: &'static str,
        /// Required element count or bytes, as named by `quantity`.
        required: u128,
        /// Target limit in the same unit as `required`.
        limit: u128,
    },
    /// The explicit work budget is insufficient.
    WorkBudgetExceeded {
        /// Required total schema-v2 work units.
        required: u128,
        /// Available units.
        available: u128,
    },
    /// The explicit logical-state budget is insufficient.
    MemoryBudgetExceeded {
        /// Required bytes.
        required: u128,
        /// Available bytes.
        available: u128,
    },
}

impl core::fmt::Display for CbcAdmissionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            Self::InvalidPointCount { point_count } => {
                write!(formatter, "CBC point count {point_count} is below 3")
            }
            Self::InvalidDimension { dimension } => {
                write!(formatter, "CBC dimension {dimension} is below 1")
            }
            Self::EstimateOverflow { quantity } => {
                write!(formatter, "CBC estimate overflowed at {quantity}")
            }
            Self::TargetCapacityExceeded {
                quantity,
                required,
                limit,
            } => write!(
                formatter,
                "CBC {quantity} needs {required} but this target permits at most {limit}"
            ),
            Self::WorkBudgetExceeded {
                required,
                available,
            } => write!(
                formatter,
                "CBC work needs {required} schema-v2 units but budget provides {available}"
            ),
            Self::MemoryBudgetExceeded {
                required,
                available,
            } => write!(
                formatter,
                "CBC state needs {required} bytes but budget provides {available}"
            ),
        }
    }
}

impl std::error::Error for CbcAdmissionError {}

const fn overflow(quantity: &'static str) -> CbcAdmissionError {
    CbcAdmissionError::EstimateOverflow { quantity }
}

const fn bit_width(value: u128) -> u32 {
    u128::BITS - value.leading_zeros()
}

const fn ceil_log2(value: u32) -> u32 {
    u32::BITS - (value - 1).leading_zeros()
}

fn limbs_for_bits(bits: u128) -> Result<u128, CbcAdmissionError> {
    bits.checked_add(31)
        .map(|rounded| rounded / 32)
        .ok_or_else(|| overflow("limb rounding"))
}

fn target_bound(
    quantity: &'static str,
    required: u128,
    limit: u128,
) -> Result<(), CbcAdmissionError> {
    if required > limit {
        return Err(CbcAdmissionError::TargetCapacityExceeded {
            quantity,
            required,
            limit,
        });
    }
    Ok(())
}

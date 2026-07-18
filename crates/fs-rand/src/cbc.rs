//! Budget-first admission for exact component-by-component lattice construction.
//!
//! This module validates the count domain, computes checked exact-integer
//! limb/scalar-work and requested-capacity memory envelopes, and admits only
//! budgets that cover those envelopes on the current target. The resulting
//! schedule is stored in the sealed receipt and consumed by
//! [`crate::cbc_exec::CbcExecutor`]; execution does not mirror its constants.
//!
//! Work is a deterministic logical debit schedule, not a bound on allocator
//! internals, instructions, energy, or elapsed time. Memory means requested
//! `Vec` payload capacity plus actual Rust owner/layout bytes and documented
//! phase overlap. `Vec` may expose more capacity than requested, so allocator
//! rounding, metadata, allocation failure, complete call stacks, and process
//! RSS remain explicit no-claims.

/// Version of the CBC admission and resource-estimate semantics.
pub const CBC_ADMISSION_SCHEMA_VERSION: u32 = 4;

const LIMB_BYTES: u128 = 4;
const FACTOR_SCRATCH_BYTES: u128 = 4 * LIMB_BYTES;
// Scalar work units are source-level logical charges, not CPU instructions or
// elapsed-time predictions. These constants are retained by schema v4. A visit
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
const CERTIFICATE_UNITS_PER_CANDIDATE: u128 = 1;
const CERTIFICATE_FIXED_UNITS_PER_PREFIX: u128 = 6;

/// Execution capability sealed into a CBC admission receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbcExecutionMode {
    /// Construct the exact canonical generator without retaining certificates.
    Construction,
    /// Construct the same generator and retain one certificate per scanned
    /// component under the certificate-capable work/memory envelope.
    Certified,
}

/// Executor-facing debit schedule derived by the admission authority.
///
/// This stays crate-private: callers budget against [`CbcEstimate`], while
/// the executor consumes these already-checked per-unit charges instead of
/// mirroring schema constants or reconstructing limb-width formulas.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CbcExecutionSchedule {
    candidate_visit_units: u128,
    candidate_control_units: u128,
    product_update_visit_units: u128,
    initialization_visit_units: u128,
    prefix_control_units: u128,
    candidate_visit_limb_units: u128,
    candidate_visit_scalar_units: u128,
    candidate_control_limb_units: u128,
    candidate_control_scalar_units: u128,
    product_update_visit_limb_units: u128,
    product_update_visit_scalar_units: u128,
    certificate_candidate_units: u128,
    certificate_score_copy_limb_units: u128,
    certificate_fixed_prefix_units: u128,
}

impl CbcExecutionSchedule {
    fn checked(
        point_count: u32,
        kernel_factor_limbs: u128,
        max_source_product_limbs: u128,
        score_capacity_limbs: u128,
        product_capacity_limbs: u128,
        max_score_limbs: u128,
    ) -> Result<Self, CbcAdmissionError> {
        let visit_limb_units = max_source_product_limbs
            .checked_mul(kernel_factor_limbs)
            .and_then(|multiply_add| {
                max_source_product_limbs
                    .checked_mul(score_capacity_limbs)
                    .and_then(|carry| multiply_add.checked_add(carry))
            })
            .ok_or_else(|| overflow("execution visit limb work"))?;
        let visit_scalar_units = kernel_factor_limbs
            .checked_mul(SCALAR_UNITS_PER_FACTOR_LIMB)
            .and_then(|factor_limbs| {
                SCALAR_UNITS_PER_LATTICE_VISIT
                    .checked_add(SCALAR_UNITS_PER_FACTOR)
                    .and_then(|fixed| fixed.checked_add(factor_limbs))
            })
            .ok_or_else(|| overflow("execution visit scalar work"))?;
        let candidate_visit_units = visit_limb_units
            .checked_add(visit_scalar_units)
            .ok_or_else(|| overflow("execution candidate visit work"))?;
        let gcd_step_upper_bound = u128::from(ceil_log2(point_count))
            .checked_mul(2)
            .and_then(|steps| steps.checked_add(1))
            .ok_or_else(|| overflow("execution GCD step bound"))?;
        let candidate_control_scalar_units = gcd_step_upper_bound
            .checked_mul(SCALAR_UNITS_PER_GCD_STEP)
            .and_then(|gcd_units| gcd_units.checked_add(SCALAR_UNITS_PER_CANDIDATE))
            .ok_or_else(|| overflow("execution candidate scalar control work"))?;
        let candidate_control_limb_units = score_capacity_limbs
            .checked_mul(2)
            .and_then(|score_units| score_units.checked_add(max_score_limbs))
            .ok_or_else(|| overflow("execution candidate limb control work"))?;
        let candidate_control_units = candidate_control_scalar_units
            .checked_add(candidate_control_limb_units)
            .ok_or_else(|| overflow("execution candidate control work"))?;
        let product_update_limb_units = product_capacity_limbs
            .checked_mul(2)
            .and_then(|product_units| visit_limb_units.checked_add(product_units))
            .ok_or_else(|| overflow("execution product-update limb work"))?;
        let product_update_visit_units = product_update_limb_units
            .checked_add(visit_scalar_units)
            .ok_or_else(|| overflow("execution product-update work"))?;
        let initialization_visit_units = product_update_visit_units
            .checked_add(1)
            .ok_or_else(|| overflow("execution initialization work"))?;
        let prefix_control_units = SCALAR_UNITS_PER_DIMENSION
            .checked_add(1)
            .ok_or_else(|| overflow("execution prefix control work"))?;
        Ok(Self {
            candidate_visit_units,
            candidate_control_units,
            product_update_visit_units,
            initialization_visit_units,
            prefix_control_units,
            candidate_visit_limb_units: visit_limb_units,
            candidate_visit_scalar_units: visit_scalar_units,
            candidate_control_limb_units,
            candidate_control_scalar_units,
            product_update_visit_limb_units: product_update_limb_units,
            product_update_visit_scalar_units: visit_scalar_units,
            certificate_candidate_units: CERTIFICATE_UNITS_PER_CANDIDATE,
            certificate_score_copy_limb_units: score_capacity_limbs
                .checked_mul(2)
                .ok_or_else(|| overflow("certificate score-copy work"))?,
            certificate_fixed_prefix_units: CERTIFICATE_FIXED_UNITS_PER_PREFIX,
        })
    }

    fn checked_upper_bound(
        self,
        candidate_visits: u128,
        candidate_count: u128,
        product_update_visits: u128,
        initialization_points: u128,
        prefixes: u128,
    ) -> Result<u128, CbcAdmissionError> {
        candidate_visits
            .checked_mul(self.candidate_visit_units)
            .and_then(|units| {
                candidate_count
                    .checked_mul(self.candidate_control_units)
                    .and_then(|control| units.checked_add(control))
            })
            .and_then(|units| {
                product_update_visits
                    .checked_mul(self.product_update_visit_units)
                    .and_then(|updates| units.checked_add(updates))
            })
            .and_then(|units| units.checked_add(initialization_points))
            .and_then(|units| {
                prefixes
                    .checked_mul(self.prefix_control_units)
                    .and_then(|control| units.checked_add(control))
            })
            .ok_or_else(|| overflow("execution schedule upper bound"))
    }

    fn checked_certificate_upper_bound(
        self,
        candidate_count: u128,
        scanned_prefixes: u128,
        prefix_payload_words: u128,
    ) -> Result<u128, CbcAdmissionError> {
        candidate_count
            .checked_mul(self.certificate_candidate_units)
            .and_then(|units| {
                scanned_prefixes
                    .checked_mul(
                        self.certificate_score_copy_limb_units
                            .checked_add(self.certificate_fixed_prefix_units)?,
                    )
                    .and_then(|prefix_units| units.checked_add(prefix_units))
            })
            .and_then(|units| units.checked_add(prefix_payload_words))
            .ok_or_else(|| overflow("certificate execution schedule"))
    }

    fn checked_limb_upper_bound(
        self,
        candidate_visits: u128,
        candidate_count: u128,
        product_update_visits: u128,
    ) -> Result<u128, CbcAdmissionError> {
        candidate_visits
            .checked_mul(self.candidate_visit_limb_units)
            .and_then(|units| {
                candidate_count
                    .checked_mul(self.candidate_control_limb_units)
                    .and_then(|control| units.checked_add(control))
            })
            .and_then(|units| {
                product_update_visits
                    .checked_mul(self.product_update_visit_limb_units)
                    .and_then(|updates| units.checked_add(updates))
            })
            .ok_or_else(|| overflow("execution limb schedule upper bound"))
    }

    fn checked_scalar_upper_bound(
        self,
        candidate_visits: u128,
        candidate_count: u128,
        product_update_visits: u128,
        initialization_points: u128,
        prefixes: u128,
    ) -> Result<u128, CbcAdmissionError> {
        candidate_visits
            .checked_mul(self.candidate_visit_scalar_units)
            .and_then(|units| {
                candidate_count
                    .checked_mul(self.candidate_control_scalar_units)
                    .and_then(|control| units.checked_add(control))
            })
            .and_then(|units| {
                product_update_visits
                    .checked_mul(self.product_update_visit_scalar_units)
                    .and_then(|updates| units.checked_add(updates))
            })
            .and_then(|units| units.checked_add(initialization_points))
            .and_then(|units| {
                prefixes
                    .checked_mul(self.prefix_control_units)
                    .and_then(|control| units.checked_add(control))
            })
            .ok_or_else(|| overflow("execution scalar schedule upper bound"))
    }

    pub(crate) const fn candidate_visit_units(self) -> u128 {
        self.candidate_visit_units
    }

    pub(crate) const fn candidate_control_units(self) -> u128 {
        self.candidate_control_units
    }

    pub(crate) const fn product_update_visit_units(self) -> u128 {
        self.product_update_visit_units
    }

    pub(crate) const fn initialization_visit_units(self) -> u128 {
        self.initialization_visit_units
    }

    pub(crate) const fn prefix_control_units(self) -> u128 {
        self.prefix_control_units
    }

    pub(crate) const fn certificate_candidate_units(self) -> u128 {
        self.certificate_candidate_units
    }

    pub(crate) fn certificate_prefix_units(self, prefix_len: usize) -> Option<u128> {
        u128::try_from(prefix_len)
            .ok()?
            .checked_add(self.certificate_score_copy_limb_units)?
            .checked_add(self.certificate_fixed_prefix_units)
    }
}

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
    /// `lattice_visits` counts the first-component product update plus every
    /// coprime candidate/point and chosen-product visit:
    ///
    /// `n + (dimension - 1) * n * phi(n) + (dimension - 1) * n`.
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
    /// length/byte allocation cannot be represented on this target or the
    /// simultaneously live logical state exceeds the target address-space
    /// cardinality.
    #[must_use]
    pub fn estimate(self) -> Result<CbcEstimate, CbcAdmissionError> {
        self.estimate_for(CbcExecutionMode::Construction)
    }

    /// Compute the checked resource envelope for one explicit execution mode.
    /// Certified mode includes certificate bookkeeping, retained records, and
    /// record-emission overlap; construction mode cannot later enable them.
    ///
    /// # Errors
    /// The same structural, arithmetic, and target-capacity refusals as
    /// [`Self::estimate`].
    #[must_use]
    #[allow(clippy::too_many_lines)] // One checked derivation keeps the envelope auditable.
    pub fn estimate_for(self, mode: CbcExecutionMode) -> Result<CbcEstimate, CbcAdmissionError> {
        let points = u128::from(self.point_count);
        let dimension =
            u128::try_from(self.dimension).map_err(|_| overflow("dimension conversion"))?;
        let later_dimensions = dimension
            .checked_sub(1)
            .ok_or_else(|| overflow("later dimensions"))?;
        let candidate_upper_bound = self.point_count - 1;
        let candidates = u128::from(candidate_upper_bound);
        let admissible_candidates_per_prefix = u128::from(euler_totient(self.point_count));

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
        let execution_schedule = CbcExecutionSchedule::checked(
            self.point_count,
            kernel_factor_limbs,
            max_source_product_limbs,
            score_capacity_limbs,
            product_capacity_limbs,
            max_score_limbs,
        )?;

        // The executor pre-reserves every resident product at the final
        // admitted capacity before the first update. The moved allocation
        // retained during replacement therefore has this same requested
        // capacity even when its normalized mathematical value is still one
        // limb (notably at dimension one).
        let previous_product_capacity_limbs = product_capacity_limbs;

        let candidate_count = later_dimensions
            .checked_mul(candidates)
            .ok_or_else(|| overflow("candidate count"))?;
        let admissible_candidate_count = later_dimensions
            .checked_mul(admissible_candidates_per_prefix)
            .ok_or_else(|| overflow("admissible candidate count"))?;
        let candidate_visits = later_dimensions
            .checked_mul(points)
            .and_then(|visits| visits.checked_mul(admissible_candidates_per_prefix))
            .ok_or_else(|| overflow("candidate lattice visits"))?;
        let product_update_visits = dimension
            .checked_mul(points)
            .ok_or_else(|| overflow("product-update visits"))?;
        let lattice_visits = candidate_visits
            .checked_add(product_update_visits)
            .ok_or_else(|| overflow("total lattice visits"))?;
        let comparison_count = candidate_count;

        let prefix_payload_words = match mode {
            CbcExecutionMode::Construction => 0,
            CbcExecutionMode::Certified if later_dimensions == 0 => 0,
            CbcExecutionMode::Certified => checked_triangular(dimension)?
                .checked_sub(1)
                .ok_or_else(|| overflow("certificate prefix payload"))?,
        };

        // Model mutually exclusive high-water phases. Candidate evaluation
        // retains current/incumbent scores (plus runner-up in certified mode).
        // Product update has no score, but multiplication simultaneously owns
        // the moved old allocation and its replacement. Allocator rounding is
        // explicitly outside this requested-capacity envelope.
        let exact_nat_owner_bytes = u128::try_from(core::mem::size_of::<crate::qmc::ExactNat>())
            .map_err(|_| overflow("exact-integer owner conversion"))?;
        let executor_inline_bytes =
            u128::try_from(core::mem::size_of::<crate::cbc_exec::CbcExecutor>())
                .map_err(|_| overflow("executor inline-state conversion"))?;
        let product_owner_array_bytes = points
            .checked_mul(exact_nat_owner_bytes)
            .ok_or_else(|| overflow("product owner array"))?;
        let generator_payload_bytes = dimension
            .checked_mul(LIMB_BYTES)
            .ok_or_else(|| overflow("generator payload"))?;
        // The executor reserves every resident product at the final admitted
        // capacity from construction entry onward. Model that actual layout,
        // rather than the smaller mathematical capacity of a prior prefix.
        let resident_product_payload_bytes = points
            .checked_mul(product_capacity_limbs)
            .and_then(|limbs| limbs.checked_mul(LIMB_BYTES))
            .ok_or_else(|| overflow("resident product payload"))?;
        // `ExactNat::mul_assign_factor_with_capacity` temporarily owns the
        // moved old allocation while populating a new allocation requested at
        // the same admitted capacity.
        let product_overlap_bytes = product_capacity_limbs
            .checked_mul(LIMB_BYTES)
            .ok_or_else(|| overflow("product-update overlap"))?;
        let score_payload_bytes = score_capacity_limbs
            .checked_mul(LIMB_BYTES)
            .ok_or_else(|| overflow("score payload"))?;
        let certificate_count = match mode {
            CbcExecutionMode::Construction => 0,
            CbcExecutionMode::Certified => later_dimensions,
        };
        let certificate_owner_bytes = certificate_count
            .checked_mul(
                u128::try_from(core::mem::size_of::<crate::cbc_cert::CbcPrefixCertificate>())
                    .map_err(|_| overflow("certificate owner conversion"))?,
            )
            .ok_or_else(|| overflow("certificate owner array"))?;
        let certificate_prefix_payload_bytes = match mode {
            CbcExecutionMode::Construction => 0,
            CbcExecutionMode::Certified => prefix_payload_words
                .checked_mul(LIMB_BYTES)
                .ok_or_else(|| overflow("certificate prefix payload bytes"))?,
        };
        let certificate_score_payload_bytes = certificate_count
            .checked_mul(2)
            .and_then(|scores| scores.checked_mul(score_payload_bytes))
            .ok_or_else(|| overflow("certificate score payload bytes"))?;
        let certificate_tie_payload_bytes = certificate_count
            .checked_mul(admissible_candidates_per_prefix)
            .and_then(|words| words.checked_mul(LIMB_BYTES))
            .ok_or_else(|| overflow("certificate tie payload bytes"))?;
        let certificate_retained_bytes = certificate_owner_bytes
            .checked_add(certificate_prefix_payload_bytes)
            .and_then(|bytes| bytes.checked_add(certificate_score_payload_bytes))
            .and_then(|bytes| bytes.checked_add(certificate_tie_payload_bytes))
            .ok_or_else(|| overflow("certificate retained state"))?;
        // `size_of::<CbcExecutor>()` binds every inline cursor, schedule,
        // owner, counter, and phase discriminant. Heap-owned ExactNat records,
        // requested payloads, and the fixed factor-decomposition scratch are
        // charged separately.
        let common_state_bytes = executor_inline_bytes
            .checked_add(product_owner_array_bytes)
            .and_then(|bytes| bytes.checked_add(resident_product_payload_bytes))
            .and_then(|bytes| bytes.checked_add(generator_payload_bytes))
            .and_then(|bytes| bytes.checked_add(FACTOR_SCRATCH_BYTES))
            .and_then(|bytes| bytes.checked_add(certificate_retained_bytes))
            .ok_or_else(|| overflow("common logical state"))?;
        // The maximum-sized Phase is already in the executor's inline bytes;
        // candidate evaluation adds current/incumbent payloads and, for a
        // certified receipt, the simultaneously retained runner-up payload.
        let candidate_phase_bytes = if candidate_count == 0 {
            0
        } else {
            let live_score_count: u128 = match mode {
                CbcExecutionMode::Construction => 2,
                CbcExecutionMode::Certified => 3,
            };
            let live_score_payload_bytes = live_score_count
                .checked_mul(score_payload_bytes)
                .ok_or_else(|| overflow("candidate score payloads"))?;
            common_state_bytes
                .checked_add(live_score_payload_bytes)
                .ok_or_else(|| overflow("candidate-phase logical state"))?
        };
        let update_phase_bytes = common_state_bytes
            .checked_add(product_overlap_bytes)
            .ok_or_else(|| overflow("update-phase logical state"))?;
        let logical_state_bytes = candidate_phase_bytes.max(update_phase_bytes);

        // Vec APIs take usize element counts and reject allocations larger than
        // isize::MAX bytes. Check each allocation independently. Aggregate
        // logical state is not one allocation and may exceed isize::MAX, but
        // simultaneously live modeled bytes cannot exceed the target's entire
        // address-space cardinality. This is an impossibility filter, not a
        // promise that the operating system can map every admitted byte.
        let usize_limit =
            u128::try_from(usize::MAX).map_err(|_| overflow("target usize conversion"))?;
        // `2^usize::BITS` is representable in u128 on today's 16/32/64-bit
        // targets. On a hypothetical 128-bit-pointer target the exact
        // cardinality is one beyond this estimator's u128 range, while every
        // representable logical-state estimate is necessarily smaller; `None`
        // therefore means that this particular impossibility check is vacuous,
        // not that admission failed open.
        let address_space_cardinality_bytes = usize_limit.checked_add(1);
        let vec_byte_limit =
            u128::try_from(isize::MAX).map_err(|_| overflow("target isize conversion"))?;
        target_bound("point-count element count", points, usize_limit)?;
        target_bound(
            "product owner-array bytes",
            product_owner_array_bytes,
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
            product_overlap_bytes,
            vec_byte_limit,
        )?;
        target_bound("score element capacity", score_capacity_limbs, usize_limit)?;
        target_bound(
            "score allocation bytes",
            score_payload_bytes,
            vec_byte_limit,
        )?;
        if matches!(mode, CbcExecutionMode::Certified) {
            target_bound(
                "certificate denominator exponent",
                dimension,
                u128::from(u32::MAX),
            )?;
            target_bound(
                "certificate owner-array bytes",
                certificate_owner_bytes,
                vec_byte_limit,
            )?;
            target_bound(
                "certificate tie element capacity",
                admissible_candidates_per_prefix,
                usize_limit,
            )?;
            target_bound(
                "certificate tie allocation bytes",
                admissible_candidates_per_prefix
                    .checked_mul(LIMB_BYTES)
                    .ok_or_else(|| overflow("certificate tie allocation"))?,
                vec_byte_limit,
            )?;
        }
        if let Some(address_space_cardinality_bytes) = address_space_cardinality_bytes {
            target_bound(
                "logical-state address-space bytes",
                logical_state_bytes,
                address_space_cardinality_bytes,
            )?;
        }

        // Derive both the public work decomposition and the executor's
        // per-boundary charges from one typed schedule. This removes copied
        // scalar constants and limb formulas from the execution module.
        let construction_limb_work_units = execution_schedule.checked_limb_upper_bound(
            candidate_visits,
            candidate_count,
            product_update_visits,
        )?;
        let construction_scalar_work_units = execution_schedule.checked_scalar_upper_bound(
            candidate_visits,
            candidate_count,
            product_update_visits,
            points,
            dimension,
        )?;
        let construction_work_units = construction_limb_work_units
            .checked_add(construction_scalar_work_units)
            .ok_or_else(|| overflow("construction work"))?;
        let scheduled_construction_work_units = execution_schedule.checked_upper_bound(
            candidate_visits,
            candidate_count,
            product_update_visits,
            points,
            dimension,
        )?;
        debug_assert_eq!(construction_work_units, scheduled_construction_work_units);
        let certificate_work_units = match mode {
            CbcExecutionMode::Construction => 0,
            CbcExecutionMode::Certified => execution_schedule.checked_certificate_upper_bound(
                candidate_count,
                later_dimensions,
                prefix_payload_words,
            )?,
        };
        let certificate_limb_work_units = certificate_count
            .checked_mul(execution_schedule.certificate_score_copy_limb_units)
            .ok_or_else(|| overflow("certificate limb work"))?;
        let certificate_scalar_work_units = certificate_work_units
            .checked_sub(certificate_limb_work_units)
            .ok_or_else(|| overflow("certificate scalar work"))?;
        let limb_work_units = construction_limb_work_units
            .checked_add(certificate_limb_work_units)
            .ok_or_else(|| overflow("total limb work"))?;
        let scalar_work_units = construction_scalar_work_units
            .checked_add(certificate_scalar_work_units)
            .ok_or_else(|| overflow("total scalar work"))?;
        let work_units = construction_work_units
            .checked_add(certificate_work_units)
            .ok_or_else(|| overflow("total work"))?;

        Ok(CbcEstimate {
            mode,
            executor_schema_version: crate::cbc_exec::CBC_EXECUTOR_SCHEMA_VERSION,
            certificate_schema_version: match mode {
                CbcExecutionMode::Construction => 0,
                CbcExecutionMode::Certified => crate::cbc_cert::CBC_CERTIFICATE_SCHEMA_VERSION,
            },
            target_pointer_width_bits: usize::BITS,
            candidate_upper_bound,
            candidate_count,
            admissible_candidates_per_prefix,
            admissible_candidate_count,
            candidate_visits,
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
            construction_limb_work_units,
            construction_scalar_work_units,
            construction_work_units,
            certificate_work_units,
            limb_work_units,
            scalar_work_units,
            work_units,
            executor_inline_bytes,
            product_owner_array_bytes,
            resident_product_payload_bytes,
            product_overlap_bytes,
            certificate_owner_bytes,
            certificate_prefix_payload_bytes,
            certificate_score_payload_bytes,
            certificate_tie_payload_bytes,
            certificate_retained_bytes,
            candidate_phase_bytes,
            update_phase_bytes,
            logical_state_bytes,
            execution_schedule,
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
        self.admit_for(CbcExecutionMode::Construction, budget)
    }

    /// Admit one explicit execution capability. A construction-only receipt
    /// cannot later enable certificates; a certified receipt covers their
    /// additional schedule and requested-capacity envelope from the start.
    ///
    /// # Errors
    /// The same estimate and budget refusals as [`Self::admit`].
    #[must_use]
    pub fn admit_for(
        self,
        mode: CbcExecutionMode,
        budget: CbcBudget,
    ) -> Result<CbcAdmission, CbcAdmissionError> {
        let estimate = self.estimate_for(mode)?;
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
            mode,
            problem: self,
            budget,
            estimate,
        })
    }
}

/// Explicit CBC admission budgets. Work units are the sum of the declared limb
/// and scalar primitive charges; memory units are requested-capacity state
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

    /// Maximum conservative schema-v4 work units.
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
    mode: CbcExecutionMode,
    executor_schema_version: u32,
    certificate_schema_version: u32,
    target_pointer_width_bits: u32,
    candidate_upper_bound: u32,
    candidate_count: u128,
    admissible_candidates_per_prefix: u128,
    admissible_candidate_count: u128,
    candidate_visits: u128,
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
    construction_limb_work_units: u128,
    construction_scalar_work_units: u128,
    construction_work_units: u128,
    certificate_work_units: u128,
    limb_work_units: u128,
    scalar_work_units: u128,
    work_units: u128,
    executor_inline_bytes: u128,
    product_owner_array_bytes: u128,
    resident_product_payload_bytes: u128,
    product_overlap_bytes: u128,
    certificate_owner_bytes: u128,
    certificate_prefix_payload_bytes: u128,
    certificate_score_payload_bytes: u128,
    certificate_tie_payload_bytes: u128,
    certificate_retained_bytes: u128,
    candidate_phase_bytes: u128,
    update_phase_bytes: u128,
    logical_state_bytes: u128,
    execution_schedule: CbcExecutionSchedule,
}

impl CbcEstimate {
    /// Execution capability covered by this estimate.
    #[must_use]
    pub const fn mode(self) -> CbcExecutionMode {
        self.mode
    }

    /// Executor schema sealed into this estimate.
    #[must_use]
    pub const fn executor_schema_version(self) -> u32 {
        self.executor_schema_version
    }

    /// Certificate schema sealed into this estimate (`0` in construction mode).
    #[must_use]
    pub const fn certificate_schema_version(self) -> u32 {
        self.certificate_schema_version
    }

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

    /// Unit-group cardinality `phi(n)` scored at each scanned prefix.
    #[must_use]
    pub const fn admissible_candidates_per_prefix(self) -> u128 {
        self.admissible_candidates_per_prefix
    }

    /// Coprime candidates that are scored across all scanned prefixes.
    #[must_use]
    pub const fn admissible_candidate_count(self) -> u128 {
        self.admissible_candidate_count
    }

    /// Exact scheduled candidate/point visits after unit-group filtering.
    #[must_use]
    pub const fn candidate_visits(self) -> u128 {
        self.candidate_visits
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

    /// Requested capacity of the retained/moved product allocation during an
    /// update. The executor pre-reserves every product at the final product
    /// capacity; allocator-reported capacity may be larger and is not claimed.
    #[must_use]
    pub const fn previous_product_capacity_limbs(self) -> u128 {
        self.previous_product_capacity_limbs
    }

    /// Exact scheduled product and coprime candidate/point visits.
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

    /// Limb work for exact construction before optional certificate emission.
    #[must_use]
    pub const fn construction_limb_work_units(self) -> u128 {
        self.construction_limb_work_units
    }

    /// Scalar work for exact construction before optional certificates.
    #[must_use]
    pub const fn construction_scalar_work_units(self) -> u128 {
        self.construction_scalar_work_units
    }

    /// Exact logical schedule consumed by a construction-only executor.
    #[must_use]
    pub const fn construction_work_units(self) -> u128 {
        self.construction_work_units
    }

    /// Additional schedule consumed only by certified execution.
    #[must_use]
    pub const fn certificate_work_units(self) -> u128 {
        self.certificate_work_units
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

    /// Inline bytes of the actual executor value, including owners, cursors,
    /// counters, schedule, and the maximum-sized phase representation.
    #[must_use]
    pub const fn executor_inline_bytes(self) -> u128 {
        self.executor_inline_bytes
    }

    /// Bytes in the heap owner array for all exact product values.
    #[must_use]
    pub const fn product_owner_array_bytes(self) -> u128 {
        self.product_owner_array_bytes
    }

    /// Requested heap payload retained by all products from executor entry.
    #[must_use]
    pub const fn resident_product_payload_bytes(self) -> u128 {
        self.resident_product_payload_bytes
    }

    /// Requested new-product payload overlapping one moved old allocation.
    #[must_use]
    pub const fn product_overlap_bytes(self) -> u128 {
        self.product_overlap_bytes
    }

    /// Requested owner-array bytes for retained certificate records.
    #[must_use]
    pub const fn certificate_owner_bytes(self) -> u128 {
        self.certificate_owner_bytes
    }

    /// Requested payload bytes for all retained certificate prefixes.
    #[must_use]
    pub const fn certificate_prefix_payload_bytes(self) -> u128 {
        self.certificate_prefix_payload_bytes
    }

    /// Requested payload bytes for retained winning and runner-up scores.
    #[must_use]
    pub const fn certificate_score_payload_bytes(self) -> u128 {
        self.certificate_score_payload_bytes
    }

    /// Requested payload bytes for retained tie classes.
    #[must_use]
    pub const fn certificate_tie_payload_bytes(self) -> u128 {
        self.certificate_tie_payload_bytes
    }

    /// Total requested certificate owner/payload capacity in the envelope.
    #[must_use]
    pub const fn certificate_retained_bytes(self) -> u128 {
        self.certificate_retained_bytes
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
    mode: CbcExecutionMode,
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

    /// Execution capability covered by this receipt.
    #[must_use]
    pub const fn mode(self) -> CbcExecutionMode {
        self.mode
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

    /// Per-boundary execution charges derived by this admission authority.
    pub(crate) fn execution_schedule(self) -> CbcExecutionSchedule {
        self.estimate.execution_schedule
    }

    pub(crate) fn has_current_authority(self) -> bool {
        self.schema_version == CBC_ADMISSION_SCHEMA_VERSION
            && self.mode == self.estimate.mode
            && self
                .problem
                .estimate_for(self.mode)
                .is_ok_and(|current| current == self.estimate)
            && self.budget.max_work_units >= self.estimate.work_units
            && self.budget.max_memory_bytes >= self.estimate.logical_state_bytes
    }
}

/// Typed CBC admission refusal. All payloads are fixed-size and allocation-free.
#[non_exhaustive]
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
    /// A modeled allocation cannot be represented by this target's
    /// `usize`/`isize`-bounded `Vec` API, or the modeled simultaneously live
    /// state exceeds the target address-space cardinality.
    TargetCapacityExceeded {
        /// Stable name of the rejected target-capacity quantity.
        quantity: &'static str,
        /// Required element count or bytes, as named by `quantity`.
        required: u128,
        /// Target limit in the same unit as `required`.
        limit: u128,
    },
    /// The explicit work budget is insufficient.
    WorkBudgetExceeded {
        /// Required total schema-v4 work units.
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
                "CBC work needs {required} schema-v4 units but budget provides {available}"
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

fn checked_triangular(value: u128) -> Result<u128, CbcAdmissionError> {
    let successor = value
        .checked_add(1)
        .ok_or_else(|| overflow("triangular successor"))?;
    let (left, right) = if value.is_multiple_of(2) {
        (value / 2, successor)
    } else {
        (value, successor / 2)
    };
    left.checked_mul(right)
        .ok_or_else(|| overflow("triangular product"))
}

fn euler_totient(mut value: u32) -> u32 {
    debug_assert!(value >= 1);
    let mut result = value;
    let mut factor = 2_u32;
    while factor <= value / factor {
        if value.is_multiple_of(factor) {
            while value.is_multiple_of(factor) {
                value /= factor;
            }
            result -= result / factor;
        }
        factor = if factor == 2 { 3 } else { factor + 2 };
    }
    if value > 1 {
        result -= result / value;
    }
    result
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

#[cfg(test)]
mod authority_tests {
    use super::{CBC_ADMISSION_SCHEMA_VERSION, CbcBudget, CbcExecutionMode, CbcProblem};
    use crate::cbc_cert::{
        ADMISSIBLE_RULE_UNITS, CbcCertError, CbcPrefixCertificate, TIE_RULE_LOWEST_CANDIDATE,
        verify_consistency_admitted,
    };
    use crate::cbc_exec::{CbcExecError, CbcExecutor};

    fn receipt() -> super::CbcAdmission {
        CbcProblem::new(5, 3)
            .expect("fixture is structural")
            .admit(CbcBudget::UNBOUNDED)
            .expect("fixture admits")
    }

    fn certified_receipt() -> super::CbcAdmission {
        CbcProblem::new(5, 3)
            .expect("fixture is structural")
            .admit_for(CbcExecutionMode::Certified, CbcBudget::UNBOUNDED)
            .expect("certified fixture admits")
    }

    fn certificate_shape() -> CbcPrefixCertificate {
        CbcPrefixCertificate {
            point_count: 5,
            prefix: vec![1, 1],
            winning_score_limbs: vec![1],
            tie_class: vec![1],
            runner_up: None,
            denominator_exponent: 2,
            tie_rule: TIE_RULE_LOWEST_CANDIDATE,
            admissible_rule: ADMISSIBLE_RULE_UNITS,
        }
    }

    #[test]
    fn g0_stale_or_mutated_authority_refuses_before_allocation() {
        CbcExecutor::new(receipt()).expect("fresh receipt matches current authority");

        let mut stale_schema = receipt();
        stale_schema.schema_version = CBC_ADMISSION_SCHEMA_VERSION - 1;
        assert!(matches!(
            CbcExecutor::new(stale_schema),
            Err(CbcExecError::AdmissionAuthorityMismatch)
        ));

        let mut changed_schedule = receipt();
        changed_schedule
            .estimate
            .execution_schedule
            .candidate_visit_units += 1;
        assert!(matches!(
            CbcExecutor::new(changed_schedule),
            Err(CbcExecError::AdmissionAuthorityMismatch)
        ));

        let mut changed_layout = receipt();
        changed_layout.estimate.executor_inline_bytes += 1;
        assert!(matches!(
            CbcExecutor::new(changed_layout),
            Err(CbcExecError::AdmissionAuthorityMismatch)
        ));

        let mut uncovered = receipt();
        uncovered.budget.max_work_units = uncovered.estimate.work_units - 1;
        assert!(matches!(
            CbcExecutor::new(uncovered),
            Err(CbcExecError::AdmissionAuthorityMismatch)
        ));
    }

    #[test]
    fn g0_certificate_checker_rejects_stale_or_mutated_authority() {
        let certificate = certificate_shape();

        let mut stale_schema = certified_receipt();
        stale_schema.schema_version = CBC_ADMISSION_SCHEMA_VERSION - 1;
        assert_eq!(
            verify_consistency_admitted(stale_schema, &certificate),
            Err(CbcCertError::AdmissionMismatch)
        );

        let mut stale_certificate_schema = certified_receipt();
        stale_certificate_schema.estimate.certificate_schema_version -= 1;
        assert_eq!(
            verify_consistency_admitted(stale_certificate_schema, &certificate),
            Err(CbcCertError::AdmissionMismatch)
        );

        let mut changed_schedule = certified_receipt();
        changed_schedule
            .estimate
            .execution_schedule
            .candidate_visit_units += 1;
        assert_eq!(
            verify_consistency_admitted(changed_schedule, &certificate),
            Err(CbcCertError::AdmissionMismatch)
        );
    }
}

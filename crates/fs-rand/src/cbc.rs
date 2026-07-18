//! Budget-first admission for exact component-by-component lattice construction.
//!
//! The live exact constructor remains [`crate::qmc::Lattice::cbc`]. This module
//! is the allocation-free prerequisite for moving that implementation behind a
//! typed execution boundary: it validates the count domain, computes checked
//! conservative exact-integer limb/work/logical-state bounds, and admits only
//! budgets that cover those bounds on the current target. It performs no CBC
//! candidate evaluation and does not allocate.
//!
//! This tranche does not yet claim that `Lattice::cbc` consumes the receipt,
//! supports cancellation or pause/resume, produces a minimality certificate,
//! or has a compact independent checker. The estimator deliberately
//! overcharges all visits at the largest retained limb widths. Its memory
//! quantity covers requested payload and owner bytes for the documented state
//! model; allocator metadata, allocator rounding, thread stacks, and process
//! RSS are explicitly outside this tranche's claim. The execution tranche
//! must use flat/exact-capacity storage or charge observed allocator capacity.

/// Version of the CBC admission and resource-estimate semantics.
pub const CBC_ADMISSION_SCHEMA_VERSION: u32 = 1;

const LIMB_BYTES: u128 = 4;
const FACTOR_SCRATCH_BYTES: u128 = 4 * LIMB_BYTES;

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
    /// `limb_work_units` then charges each visit at the maximum product/factor
    /// multiply width and accumulator carry width, plus a maximum-width score
    /// comparison for every later candidate. This is an admission envelope,
    /// not a prediction of elapsed time.
    ///
    /// # Errors
    /// [`CbcAdmissionError::EstimateOverflow`] naming the first bound that
    /// leaves the `u128` accounting domain.
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

        let max_product_bits = u128::from(kernel_numerator_bits)
            .checked_mul(dimension)
            .ok_or_else(|| overflow("product bit width"))?;
        let max_product_limbs = limbs_for_bits(max_product_bits)?;
        let score_growth_bits = u128::from(ceil_log2(self.point_count));
        let max_score_bits = max_product_bits
            .checked_add(score_growth_bits)
            .ok_or_else(|| overflow("score bit width"))?;
        let max_score_limbs = limbs_for_bits(max_score_bits)?;
        // ExactNat::add_mul_factor reserves source + factor + one spare limb.
        let accumulator_capacity_limbs = max_product_limbs
            .checked_add(kernel_factor_limbs)
            .and_then(|limbs| limbs.checked_add(1))
            .ok_or_else(|| overflow("accumulator limb capacity"))?;
        // ExactNat::mul_assign_factor moves the normalized source out, then
        // requests source + factor + one spare limb for the replacement. A
        // future execution path must retain at least this logical capacity per
        // product even though normalization shortens its visible limb length.
        let product_capacity_limbs = max_product_limbs
            .checked_add(kernel_factor_limbs)
            .and_then(|limbs| limbs.checked_add(1))
            .ok_or_else(|| overflow("product limb capacity"))?;

        let candidate_visits = later_dimensions
            .checked_mul(points)
            .and_then(|visits| visits.checked_mul(candidates))
            .ok_or_else(|| overflow("candidate lattice visits"))?;
        let chosen_update_visits = later_dimensions
            .checked_mul(points)
            .ok_or_else(|| overflow("chosen-product visits"))?;
        let lattice_visits = points
            .checked_add(candidate_visits)
            .and_then(|visits| visits.checked_add(chosen_update_visits))
            .ok_or_else(|| overflow("total lattice visits"))?;
        let comparison_count = later_dimensions
            .checked_mul(candidates)
            .ok_or_else(|| overflow("candidate comparisons"))?;

        let multiply_add_units = lattice_visits
            .checked_mul(max_product_limbs)
            .and_then(|units| units.checked_mul(kernel_factor_limbs))
            .ok_or_else(|| overflow("multiply-add limb work"))?;
        let carry_units = lattice_visits
            .checked_mul(accumulator_capacity_limbs)
            .ok_or_else(|| overflow("carry limb work"))?;
        let comparison_width = max_score_limbs.max(accumulator_capacity_limbs);
        let comparison_limb_units = comparison_count
            .checked_mul(comparison_width)
            .ok_or_else(|| overflow("comparison limb work"))?;
        let limb_work_units = multiply_add_units
            .checked_add(carry_units)
            .and_then(|units| units.checked_add(comparison_limb_units))
            .ok_or_else(|| overflow("total limb work"))?;

        // Logical schoolbook-state envelope: one Vec<u32> owner header per
        // product, an outer products header, the generator header,
        // current/best score headers, retained requested product capacity, two
        // accumulator payloads (also covering product-update overlap),
        // generator payload, and the fixed four-limb factor scratch. This is
        // not an allocator/RSS bound; see the module-level no-claim.
        let vector_header_bytes = u128::try_from(core::mem::size_of::<Vec<u32>>())
            .map_err(|_| overflow("vector header conversion"))?;
        let product_header_bytes = points
            .checked_mul(vector_header_bytes)
            .ok_or_else(|| overflow("product headers"))?;
        let product_payload_bytes = points
            .checked_mul(product_capacity_limbs)
            .and_then(|limbs| limbs.checked_mul(LIMB_BYTES))
            .ok_or_else(|| overflow("product payload"))?;
        let generator_payload_bytes = dimension
            .checked_mul(LIMB_BYTES)
            .ok_or_else(|| overflow("generator payload"))?;
        let accumulator_payload_bytes = 2_u128
            .checked_mul(accumulator_capacity_limbs)
            .and_then(|limbs| limbs.checked_mul(LIMB_BYTES))
            .ok_or_else(|| overflow("accumulator payload"))?;
        let owner_header_bytes = 4_u128
            .checked_mul(vector_header_bytes)
            .ok_or_else(|| overflow("owner headers"))?;
        let resident_bytes = product_header_bytes
            .checked_add(product_payload_bytes)
            .and_then(|bytes| bytes.checked_add(generator_payload_bytes))
            .and_then(|bytes| bytes.checked_add(accumulator_payload_bytes))
            .and_then(|bytes| bytes.checked_add(owner_header_bytes))
            .and_then(|bytes| bytes.checked_add(FACTOR_SCRATCH_BYTES))
            .ok_or_else(|| overflow("resident memory"))?;

        Ok(CbcEstimate {
            candidate_upper_bound,
            kernel_numerator_upper,
            kernel_numerator_bits,
            kernel_factor_limbs,
            max_product_bits,
            max_product_limbs,
            max_score_bits,
            max_score_limbs,
            accumulator_capacity_limbs,
            product_capacity_limbs,
            lattice_visits,
            comparison_count,
            limb_work_units,
            resident_bytes,
        })
    }

    /// Admit this problem only when both explicit budgets cover the checked
    /// estimate. Target-address refusal precedes work refusal, which precedes
    /// memory refusal.
    ///
    /// # Errors
    /// [`CbcAdmissionError::EstimateOverflow`],
    /// [`CbcAdmissionError::AddressSpaceExceeded`],
    /// [`CbcAdmissionError::WorkBudgetExceeded`], or
    /// [`CbcAdmissionError::MemoryBudgetExceeded`].
    #[must_use]
    pub fn admit(self, budget: CbcBudget) -> Result<CbcAdmission, CbcAdmissionError> {
        let estimate = self.estimate()?;
        let addressable =
            u128::try_from(isize::MAX).map_err(|_| overflow("target address-space conversion"))?;
        if estimate.resident_bytes > addressable {
            return Err(CbcAdmissionError::AddressSpaceExceeded {
                required: estimate.resident_bytes,
                addressable,
            });
        }
        if estimate.limb_work_units > budget.max_work_units {
            return Err(CbcAdmissionError::WorkBudgetExceeded {
                required: estimate.limb_work_units,
                available: budget.max_work_units,
            });
        }
        if estimate.resident_bytes > budget.max_memory_bytes {
            return Err(CbcAdmissionError::MemoryBudgetExceeded {
                required: estimate.resident_bytes,
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

/// Explicit CBC admission budgets. Units are conservative limb operations and
/// logical state bytes under this module's documented accounting model,
/// respectively.
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

    /// Maximum conservative limb-operation units.
    #[must_use]
    pub const fn max_work_units(self) -> u128 {
        self.max_work_units
    }

    /// Maximum logical state bytes admitted by this accounting model.
    #[must_use]
    pub const fn max_memory_bytes(self) -> u128 {
        self.max_memory_bytes
    }
}

/// Checked conservative resource bounds for one CBC problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CbcEstimate {
    candidate_upper_bound: u32,
    kernel_numerator_upper: u128,
    kernel_numerator_bits: u32,
    kernel_factor_limbs: u128,
    max_product_bits: u128,
    max_product_limbs: u128,
    max_score_bits: u128,
    max_score_limbs: u128,
    accumulator_capacity_limbs: u128,
    product_capacity_limbs: u128,
    lattice_visits: u128,
    comparison_count: u128,
    limb_work_units: u128,
    resident_bytes: u128,
}

impl CbcEstimate {
    /// Upper bound on candidates examined per later component (`n - 1`).
    #[must_use]
    pub const fn candidate_upper_bound(self) -> u32 {
        self.candidate_upper_bound
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

    /// Conservative temporary score capacity including factor and spare limbs.
    #[must_use]
    pub const fn accumulator_capacity_limbs(self) -> u128 {
        self.accumulator_capacity_limbs
    }

    /// Requested capacity retained per exact product in base-2^32 limbs.
    #[must_use]
    pub const fn product_capacity_limbs(self) -> u128 {
        self.product_capacity_limbs
    }

    /// Upper bound on point/candidate visits.
    #[must_use]
    pub const fn lattice_visits(self) -> u128 {
        self.lattice_visits
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

    /// Conservative logical state envelope in bytes for the documented
    /// exact-capacity layout.
    ///
    /// This includes requested payload and owner bytes but excludes allocator
    /// metadata/rounding, thread stacks, and process-level RSS effects.
    #[must_use]
    pub const fn resident_bytes(self) -> u128 {
        self.resident_bytes
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
    /// The logical state envelope exceeds one target-addressable allocation
    /// domain (`isize::MAX` bytes).
    AddressSpaceExceeded {
        /// Required logical state bytes.
        required: u128,
        /// Target-addressable byte ceiling.
        addressable: u128,
    },
    /// The explicit work budget is insufficient.
    WorkBudgetExceeded {
        /// Required conservative limb-operation units.
        required: u128,
        /// Available units.
        available: u128,
    },
    /// The explicit resident-memory budget is insufficient.
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
            Self::AddressSpaceExceeded {
                required,
                addressable,
            } => write!(
                formatter,
                "CBC state needs {required} logical bytes but this target addresses at most {addressable}"
            ),
            Self::WorkBudgetExceeded {
                required,
                available,
            } => write!(
                formatter,
                "CBC work needs {required} limb units but budget provides {available}"
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

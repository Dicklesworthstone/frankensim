//! G0 admission and resource-bound coverage for exact CBC construction (6ys.20.1).

#![deny(unsafe_code)]

use fs_rand::cbc::{CBC_ADMISSION_SCHEMA_VERSION, CbcAdmissionError, CbcBudget, CbcProblem};

#[test]
fn g0_count_domain_refuses_before_estimation() {
    for point_count in 0..3 {
        assert_eq!(
            CbcProblem::new(point_count, 1),
            Err(CbcAdmissionError::InvalidPointCount { point_count })
        );
    }
    assert_eq!(
        CbcProblem::new(3, 0),
        Err(CbcAdmissionError::InvalidDimension { dimension: 0 })
    );
}

/// Hand calculation for n=5, dimension=3:
///
/// - kernel numerator <= 7*5^2 = 175 (8 bits, one base-2^32 limb),
/// - product <= 24 bits/one limb; five-term score <= 27 bits/one limb,
/// - source, score, and product capacity are one, three, and three limbs,
/// - candidate count = 2*4 = 8; visits = 2*5*4 + 3*5 = 55,
/// - limb work = 55 multiply + 165 carry + 69 zero-fill + 69 normalize
///   + 8 compare = 366,
/// - scalar work = 1,100 visit + 495 factor + 168 GCD + 128 candidate
///   + 32 dimension/initialization = 1,923; total = 2,289,
/// - requested-payload memory uses the actual executor/ExactNat owner layout.
#[test]
fn g0_small_estimate_matches_the_independent_kat() {
    let problem = CbcProblem::new(5, 3).expect("valid finite CBC problem");
    let estimate = problem.estimate().expect("finite estimate");
    assert_eq!(problem.point_count(), 5);
    assert_eq!(problem.dimension(), 3);
    assert_eq!(estimate.candidate_upper_bound(), 4);
    assert_eq!(estimate.kernel_numerator_upper(), 175);
    assert_eq!(estimate.kernel_numerator_bits(), 8);
    assert_eq!(estimate.kernel_factor_limbs(), 1);
    assert_eq!(estimate.max_product_bits(), 24);
    assert_eq!(estimate.max_product_limbs(), 1);
    assert_eq!(estimate.max_source_product_limbs(), 1);
    assert_eq!(estimate.max_score_bits(), 27);
    assert_eq!(estimate.max_score_limbs(), 1);
    assert_eq!(estimate.score_capacity_limbs(), 3);
    assert_eq!(estimate.product_capacity_limbs(), 3);
    assert_eq!(estimate.previous_product_capacity_limbs(), 3);
    assert_eq!(estimate.candidate_count(), 8);
    assert_eq!(estimate.lattice_visits(), 55);
    assert_eq!(estimate.product_update_visits(), 15);
    assert_eq!(estimate.comparison_count(), 8);
    assert_eq!(estimate.limb_work_units(), 366);
    assert_eq!(estimate.scalar_work_units(), 1_923);
    assert_eq!(estimate.work_units(), 2_289);
    assert_eq!(estimate.target_pointer_width_bits(), usize::BITS);

    let common = estimate.executor_inline_bytes()
        + estimate.product_owner_array_bytes()
        + estimate.resident_product_payload_bytes()
        + 3 * 4 // generator payload
        + 4 * 4; // fixed factor scratch
    assert_eq!(estimate.candidate_phase_bytes(), common + 2 * 3 * 4);
    assert_eq!(estimate.update_phase_bytes(), common + 3 * 4);
    assert_eq!(
        estimate.logical_state_bytes(),
        estimate.candidate_phase_bytes()
    );
}

#[test]
fn g0_estimates_are_monotone_in_points_and_dimension() {
    let base = CbcProblem::new(5, 2)
        .expect("base problem")
        .estimate()
        .expect("base estimate");
    let wider = CbcProblem::new(5, 3)
        .expect("wider problem")
        .estimate()
        .expect("wider estimate");
    let denser = CbcProblem::new(7, 3)
        .expect("denser problem")
        .estimate()
        .expect("denser estimate");

    assert!(wider.lattice_visits() > base.lattice_visits());
    assert!(wider.limb_work_units() > base.limb_work_units());
    assert!(wider.scalar_work_units() > base.scalar_work_units());
    assert!(wider.work_units() > base.work_units());
    assert!(wider.logical_state_bytes() > base.logical_state_bytes());
    assert!(denser.candidate_upper_bound() > wider.candidate_upper_bound());
    assert!(denser.lattice_visits() > wider.lattice_visits());
    assert!(denser.limb_work_units() > wider.limb_work_units());
    assert!(denser.scalar_work_units() > wider.scalar_work_units());
    assert!(denser.work_units() > wider.work_units());
    assert!(denser.logical_state_bytes() > wider.logical_state_bytes());
}

/// Hand calculation for n=5, dimension=6, where the multiplicative source —
/// not merely the final product — first becomes multi-limb:
///
/// - final product <= 48 bits/two limbs and source <= 40 bits/two limbs,
/// - score/product capacity = 4 limbs; every resident product is pre-reserved
///   at that same four-limb requested capacity,
/// - candidates = 5*4 = 20; visits = 5*5*4 + 6*5 = 130,
/// - limb work = 260 multiply + 1,040 carry + 200 zero-fill
///   + 200 normalize + 40 compare = 1,740,
/// - scalar work = 2,600 visit + 1,170 factor + 420 GCD + 320 candidate
///   + 59 dimension/initialization = 4,569; total = 6,309.
#[test]
fn g0_multilimb_kat_charges_carry_for_every_source_limb() {
    let boundary = CbcProblem::new(5, 5)
        .expect("valid single-limb boundary problem")
        .estimate()
        .expect("finite boundary estimate");
    let estimate = CbcProblem::new(5, 6)
        .expect("valid multi-limb problem")
        .estimate()
        .expect("finite multi-limb estimate");
    assert_eq!(boundary.max_product_bits(), 40);
    assert_eq!(boundary.max_product_limbs(), 2);
    assert_eq!(boundary.max_source_product_limbs(), 1);
    assert_eq!(boundary.product_capacity_limbs(), 3);
    assert_eq!(boundary.limb_work_units(), 698);
    assert_eq!(boundary.scalar_work_units(), 3_687);
    assert_eq!(boundary.work_units(), 4_385);
    assert_eq!(estimate.max_product_bits(), 48);
    assert_eq!(estimate.max_product_limbs(), 2);
    assert_eq!(estimate.max_source_product_limbs(), 2);
    assert_eq!(estimate.max_score_bits(), 51);
    assert_eq!(estimate.max_score_limbs(), 2);
    assert_eq!(estimate.score_capacity_limbs(), 4);
    assert_eq!(estimate.product_capacity_limbs(), 4);
    assert_eq!(estimate.previous_product_capacity_limbs(), 4);
    assert_eq!(estimate.candidate_count(), 20);
    assert_eq!(estimate.lattice_visits(), 130);
    assert_eq!(estimate.product_update_visits(), 30);
    assert_eq!(estimate.comparison_count(), 20);
    assert_eq!(estimate.limb_work_units(), 1_740);
    assert_eq!(estimate.scalar_work_units(), 4_569);
    assert_eq!(estimate.work_units(), 6_309);
    assert!(estimate.work_units() > boundary.work_units());
}

/// Independent large-factor KAT. Unlike the small/multidimensional cases
/// above, `n=25_000` makes one kernel factor itself two limbs, exercising the
/// second multiplicand axis in every multiplication and factor-decomposition
/// charge without relying on an enormous dimension.
#[test]
fn g0_multilimb_kernel_factor_kat_binds_both_work_axes() {
    let estimate = CbcProblem::new(25_000, 2)
        .expect("valid multi-limb factor problem")
        .estimate()
        .expect("finite multi-limb factor estimate");
    assert_eq!(estimate.kernel_numerator_bits(), 33);
    assert_eq!(estimate.kernel_factor_limbs(), 2);
    assert_eq!(estimate.max_source_product_limbs(), 2);
    assert_eq!(estimate.admissible_candidates_per_prefix(), 10_000);
    assert_eq!(estimate.candidate_visits(), 250_000_000);
    assert_eq!(estimate.lattice_visits(), 250_050_000);
    assert_eq!(estimate.limb_work_units(), 3_501_524_987);
    assert_eq!(estimate.scalar_work_units(), 9_254_599_909);
    assert_eq!(estimate.work_units(), 12_756_124_896);
}

/// Maximum supported point-count KAT. `7*(2^32-1)^2` is 67 bits, so the
/// exact kernel factor occupies three base-2^32 limbs. Dimension one keeps the
/// problem allocation-free while independently pinning the widest live factor,
/// its multiply/carry charges, initialization debit, and update-phase memory.
#[cfg(target_pointer_width = "64")]
#[test]
fn g0_three_limb_kernel_factor_kat_reaches_u32_point_boundary() {
    let estimate = CbcProblem::new(u32::MAX, 1)
        .expect("maximum point count is structurally supported")
        .estimate()
        .expect("dimension-one boundary estimate remains representable");
    assert_eq!(
        estimate.kernel_numerator_upper(),
        129_127_208_455_837_319_175
    );
    assert_eq!(estimate.kernel_numerator_bits(), 67);
    assert_eq!(estimate.kernel_factor_limbs(), 3);
    assert_eq!(estimate.max_source_product_limbs(), 1);
    assert_eq!(estimate.product_capacity_limbs(), 5);
    assert_eq!(estimate.score_capacity_limbs(), 5);
    assert_eq!(estimate.lattice_visits(), u128::from(u32::MAX));
    assert_eq!(estimate.limb_work_units(), 77_309_411_310);
    assert_eq!(estimate.scalar_work_units(), 197_568_495_579);
    assert_eq!(estimate.work_units(), 274_877_906_889);
    assert_eq!(estimate.candidate_phase_bytes(), 0);
    let expected = estimate.executor_inline_bytes()
        + estimate.product_owner_array_bytes()
        + estimate.resident_product_payload_bytes()
        + 4
        + 16
        + estimate.product_overlap_bytes();
    assert_eq!(estimate.update_phase_bytes(), expected);
    assert_eq!(estimate.logical_state_bytes(), expected);
}

/// Bounded property grid over small structural regimes. Dimension growth is
/// monotone. Point-count work intentionally follows `n*phi(n)` and therefore
/// is not monotone across composite/prime transitions; its exact count is
/// checked directly while requested live state remains monotone.
#[test]
fn g0_resource_bounds_are_monotone_over_bounded_grid() {
    for point_count in 3..=31 {
        let mut previous = CbcProblem::new(point_count, 1)
            .expect("grid point count is valid")
            .estimate()
            .expect("small grid estimate is finite");
        for dimension in 2..=12 {
            let current = CbcProblem::new(point_count, dimension)
                .expect("grid dimension is valid")
                .estimate()
                .expect("small grid estimate is finite");
            assert_eq!(
                current.candidate_visits(),
                u128::from(point_count) * current.admissible_candidate_count(),
                "exact coprime visit schedule moved at n={point_count}, d={dimension}"
            );
            assert!(
                current.limb_work_units() >= previous.limb_work_units(),
                "limb work decreased at n={point_count}, d={dimension}"
            );
            assert!(
                current.scalar_work_units() >= previous.scalar_work_units(),
                "scalar work decreased at n={point_count}, d={dimension}"
            );
            assert!(
                current.logical_state_bytes() >= previous.logical_state_bytes(),
                "state decreased at n={point_count}, d={dimension}"
            );
            previous = current;
        }
    }

    for dimension in 1..=12 {
        let mut previous = CbcProblem::new(3, dimension)
            .expect("grid dimension is valid")
            .estimate()
            .expect("small grid estimate is finite");
        for point_count in 4..=31 {
            let current = CbcProblem::new(point_count, dimension)
                .expect("grid point count is valid")
                .estimate()
                .expect("small grid estimate is finite");
            assert_eq!(
                current.candidate_visits(),
                u128::from(point_count) * current.admissible_candidate_count(),
                "exact coprime visit schedule moved at n={point_count}, d={dimension}"
            );
            assert!(
                current.logical_state_bytes() >= previous.logical_state_bytes(),
                "state decreased at n={point_count}, d={dimension}"
            );
            previous = current;
        }
    }
}

#[test]
fn g0_each_vec_capacity_is_checked_before_work_or_budget() {
    let problem = CbcProblem::new(3, usize::MAX).expect("structural counts are nonzero");
    let required = u128::try_from(usize::MAX)
        .expect("usize fits u128")
        .checked_mul(4)
        .expect("four times a supported usize fits u128");
    let limit = u128::try_from(isize::MAX).expect("positive isize maximum fits u128");
    let expected = CbcAdmissionError::TargetCapacityExceeded {
        quantity: "generator allocation bytes",
        required,
        limit,
    };
    assert_eq!(problem.estimate(), Err(expected));
    assert_eq!(problem.admit(CbcBudget::new(0, 0)), Err(expected));
}

#[cfg(target_pointer_width = "64")]
#[test]
fn g0_individual_multilimb_product_bytes_are_checked_on_64_bit() {
    // Generator storage is 4.8e18 bytes and therefore representable, while a
    // single 67-bit-factor product needs more than isize::MAX bytes.
    let problem = CbcProblem::new(u32::MAX, 1_200_000_000_000_000_000)
        .expect("large 64-bit counts remain structurally valid");
    let addressable = u128::try_from(isize::MAX).expect("positive isize maximum fits u128");
    assert!(matches!(
        problem.estimate(),
        Err(CbcAdmissionError::TargetCapacityExceeded {
            quantity: "product allocation bytes",
            required,
            limit,
        }) if required > limit && limit == addressable
    ));
}

#[cfg(target_pointer_width = "32")]
#[test]
fn g0_outer_product_owner_array_is_checked_on_32_bit() {
    let point_count = 200_000_000;
    let problem = CbcProblem::new(point_count, 1).expect("valid count domain");
    let header = u128::try_from(core::mem::size_of::<Vec<u32>>()).expect("header fits u128");
    let required = u128::from(point_count) * header;
    let limit = u128::try_from(isize::MAX).expect("positive isize maximum fits u128");
    assert_eq!(
        problem.estimate(),
        Err(CbcAdmissionError::TargetCapacityExceeded {
            quantity: "product owner-array bytes",
            required,
            limit,
        })
    );
}

#[cfg(target_pointer_width = "64")]
#[test]
fn g0_multilimb_carry_work_overflow_is_fail_closed() {
    // Every modeled Vec remains individually addressable at this dimension;
    // only the checked O(n^2 d^3) carry charge leaves u128.
    let problem = CbcProblem::new(3, 100_000_000_000_000)
        .expect("large 64-bit counts remain structurally valid");
    assert!(matches!(
        problem.estimate(),
        Err(CbcAdmissionError::EstimateOverflow { .. })
    ));
    assert!(matches!(
        problem.admit(CbcBudget::new(0, 0)),
        Err(CbcAdmissionError::EstimateOverflow { .. })
    ));
}

// Aggregate logical state is not one Vec allocation. On 32-bit this fixture's
// phase maximum exceeds isize::MAX, while the outer array, generator, each
// product, and each score allocation remains individually representable.
#[cfg(target_pointer_width = "32")]
#[test]
fn g0_aggregate_state_above_isize_is_not_mistaken_for_one_allocation() {
    let problem = CbcProblem::new(3, 310_000_000).expect("large 32-bit counts remain valid");
    let estimate = problem.estimate().expect("each allocation is addressable");
    let one_allocation_limit =
        u128::try_from(isize::MAX).expect("positive isize maximum fits u128");
    assert!(estimate.logical_state_bytes() > one_allocation_limit);
    problem
        .admit(CbcBudget::UNBOUNDED)
        .expect("aggregate logical bytes are not a single allocation");
}

// Every individual Vec in this fixture remains below isize::MAX, but the
// simultaneously live candidate phase is larger than the entire 32-bit
// address-space cardinality. Per-allocation checks alone therefore cannot
// authenticate the logical-state proposition.
#[cfg(target_pointer_width = "32")]
#[test]
fn g0_aggregate_state_above_address_space_is_refused() {
    let problem = CbcProblem::new(5, 500_000_000).expect("large 32-bit counts remain structural");
    for result in [
        problem.estimate().map(|_| ()),
        problem.admit(CbcBudget::UNBOUNDED).map(|_| ()),
        problem.admit(CbcBudget::new(0, 0)).map(|_| ()),
    ] {
        assert!(matches!(
            result,
            Err(CbcAdmissionError::TargetCapacityExceeded {
                quantity: "logical-state address-space bytes",
                required,
                limit,
            }) if required > limit && limit == 1_u128 << 32
        ));
    }
}

// Exercise the aggregate theorem on the reference 64-bit lane as well. Every
// modeled Vec remains below isize::MAX, while the exact candidate-phase state
// exceeds all 2^64 byte addresses and must refuse before budget comparison.
#[cfg(target_pointer_width = "64")]
#[test]
fn g0_aggregate_state_above_64_bit_address_space_is_refused() {
    let problem = CbcProblem::new(5, 1_800_000_000_000_000_000)
        .expect("large 64-bit counts remain structural");
    for result in [
        problem.estimate().map(|_| ()),
        problem.admit(CbcBudget::UNBOUNDED).map(|_| ()),
        problem.admit(CbcBudget::new(0, 0)).map(|_| ()),
    ] {
        assert!(matches!(
            result,
            Err(CbcAdmissionError::TargetCapacityExceeded {
                quantity: "logical-state address-space bytes",
                required,
                limit,
            }) if required > limit && limit == 1_u128 << 64
        ));
    }
}

#[test]
fn g0_dimension_one_memory_charges_moved_old_and_new_product_overlap() {
    let estimate = CbcProblem::new(5, 1)
        .expect("valid one-dimensional problem")
        .estimate()
        .expect("finite estimate");
    assert_eq!(estimate.candidate_count(), 0);
    assert_eq!(estimate.candidate_phase_bytes(), 0);
    assert_eq!(estimate.previous_product_capacity_limbs(), 3);
    assert_eq!(estimate.product_capacity_limbs(), 3);
    let expected = estimate.executor_inline_bytes()
        + estimate.product_owner_array_bytes()
        + estimate.resident_product_payload_bytes()
        + 4 // generator word
        + 16 // fixed factor scratch
        + estimate.product_overlap_bytes();
    assert_eq!(estimate.update_phase_bytes(), expected);
    assert_eq!(
        estimate.logical_state_bytes(),
        estimate.update_phase_bytes()
    );
}

#[test]
fn g0_work_and_memory_budgets_have_exact_boundaries() {
    let problem = CbcProblem::new(5, 3).expect("valid problem");
    let estimate = problem.estimate().expect("finite estimate");
    let exact = CbcBudget::new(estimate.work_units(), estimate.logical_state_bytes());
    let admission = problem.admit(exact).expect("exact budgets admit");
    assert_eq!(CBC_ADMISSION_SCHEMA_VERSION, 4);
    assert_eq!(admission.schema_version(), CBC_ADMISSION_SCHEMA_VERSION);
    assert_eq!(admission.problem(), problem);
    assert_eq!(admission.budget(), exact);
    assert_eq!(admission.estimate(), estimate);
    assert_eq!(exact.max_work_units(), estimate.work_units());
    assert_eq!(exact.max_memory_bytes(), estimate.logical_state_bytes());

    // Schema v1 admitted this budget because it omitted scalar, zero-fill,
    // normalization, and per-source carry charges. Schemas v2-v4 must not
    // silently preserve that undercount.
    assert_eq!(
        problem.admit(CbcBudget::new(244, estimate.logical_state_bytes())),
        Err(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.work_units(),
            available: 244,
        })
    );

    assert_eq!(
        problem.admit(CbcBudget::new(0, 0)),
        Err(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.work_units(),
            available: 0,
        })
    );
    assert_eq!(
        problem.admit(CbcBudget::new(estimate.work_units() - 1, u128::MAX,)),
        Err(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.work_units(),
            available: estimate.work_units() - 1,
        })
    );
    assert_eq!(
        problem.admit(CbcBudget::new(estimate.work_units(), 0)),
        Err(CbcAdmissionError::MemoryBudgetExceeded {
            required: estimate.logical_state_bytes(),
            available: 0,
        })
    );
    assert_eq!(
        problem.admit(CbcBudget::new(
            estimate.work_units(),
            estimate.logical_state_bytes() - 1,
        )),
        Err(CbcAdmissionError::MemoryBudgetExceeded {
            required: estimate.logical_state_bytes(),
            available: estimate.logical_state_bytes() - 1,
        })
    );
}

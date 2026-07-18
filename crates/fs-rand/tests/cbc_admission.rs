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
/// - score/product capacity reserves 1 product + 1 factor + 1 spare = 3 limbs,
/// - visits = 5 + 2*5*4 + 2*5 = 55; comparisons = 2*4 = 8,
/// - work = 55*(1*1) + 55*3 + 8*3 = 244 limb units,
/// - memory = 112 logical payload/scratch bytes + 9 Vec headers.
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
    assert_eq!(estimate.max_score_bits(), 27);
    assert_eq!(estimate.max_score_limbs(), 1);
    assert_eq!(estimate.accumulator_capacity_limbs(), 3);
    assert_eq!(estimate.product_capacity_limbs(), 3);
    assert_eq!(estimate.lattice_visits(), 55);
    assert_eq!(estimate.comparison_count(), 8);
    assert_eq!(estimate.limb_work_units(), 244);

    let vector_header =
        u128::try_from(core::mem::size_of::<Vec<u32>>()).expect("Vec header size fits u128");
    assert_eq!(estimate.resident_bytes(), 112 + 9 * vector_header);
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
    assert!(wider.resident_bytes() > base.resident_bytes());
    assert!(denser.candidate_upper_bound() > wider.candidate_upper_bound());
    assert!(denser.lattice_visits() > wider.lattice_visits());
    assert!(denser.limb_work_units() > wider.limb_work_units());
    assert!(denser.resident_bytes() > wider.resident_bytes());
}

#[test]
fn g0_hostile_counts_fail_closed_on_checked_estimation() {
    let problem = CbcProblem::new(u32::MAX, usize::MAX).expect("structural counts are nonzero");
    assert!(matches!(
        problem.estimate(),
        Err(CbcAdmissionError::EstimateOverflow { .. })
    ));
    assert!(matches!(
        problem.admit(CbcBudget::UNBOUNDED),
        Err(CbcAdmissionError::EstimateOverflow { .. })
    ));
}

#[test]
fn g0_finite_estimate_beyond_target_address_space_is_not_admitted() {
    let problem =
        CbcProblem::new(u32::MAX, 1_000_000_000).expect("large counts remain structurally valid");
    let estimate = problem
        .estimate()
        .expect("the accounting expression remains within u128");
    let addressable = u128::try_from(isize::MAX).expect("positive isize maximum fits u128");
    assert!(estimate.resident_bytes() > addressable);
    assert_eq!(
        problem.admit(CbcBudget::UNBOUNDED),
        Err(CbcAdmissionError::AddressSpaceExceeded {
            required: estimate.resident_bytes(),
            addressable,
        })
    );
    assert_eq!(
        problem.admit(CbcBudget::new(0, 0)),
        Err(CbcAdmissionError::AddressSpaceExceeded {
            required: estimate.resident_bytes(),
            addressable,
        }),
        "target addressability must refuse before either exhausted budget"
    );
}

#[test]
fn g0_work_and_memory_budgets_have_exact_boundaries() {
    let problem = CbcProblem::new(5, 3).expect("valid problem");
    let estimate = problem.estimate().expect("finite estimate");
    let exact = CbcBudget::new(estimate.limb_work_units(), estimate.resident_bytes());
    let admission = problem.admit(exact).expect("exact budgets admit");
    assert_eq!(admission.schema_version(), CBC_ADMISSION_SCHEMA_VERSION);
    assert_eq!(admission.problem(), problem);
    assert_eq!(admission.budget(), exact);
    assert_eq!(admission.estimate(), estimate);
    assert_eq!(exact.max_work_units(), estimate.limb_work_units());
    assert_eq!(exact.max_memory_bytes(), estimate.resident_bytes());

    assert_eq!(
        problem.admit(CbcBudget::new(0, 0)),
        Err(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.limb_work_units(),
            available: 0,
        })
    );
    assert_eq!(
        problem.admit(CbcBudget::new(estimate.limb_work_units() - 1, u128::MAX,)),
        Err(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.limb_work_units(),
            available: estimate.limb_work_units() - 1,
        })
    );
    assert_eq!(
        problem.admit(CbcBudget::new(estimate.limb_work_units(), 0)),
        Err(CbcAdmissionError::MemoryBudgetExceeded {
            required: estimate.resident_bytes(),
            available: 0,
        })
    );
    assert_eq!(
        problem.admit(CbcBudget::new(
            estimate.limb_work_units(),
            estimate.resident_bytes() - 1,
        )),
        Err(CbcAdmissionError::MemoryBudgetExceeded {
            required: estimate.resident_bytes(),
            available: estimate.resident_bytes() - 1,
        })
    );
}

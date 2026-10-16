//! G0/G3 integration evidence for the admitted synchronous CBC lattice facade
//! (bead 6ys.20.2).

use fs_rand::{
    cbc::{CbcAdmissionError, CbcBudget, CbcProblem},
    qmc::{CbcLatticeError, DEFAULT_CBC_BUDGET, Lattice},
};

const TIE_N: u32 = 5;
const TIE_DIM: usize = 6;
const TIE_GENERATOR: [u32; TIE_DIM] = [1, 2, 1, 2, 1, 2];
const TIE_GENERATOR_LE_BYTES: [u8; TIE_DIM * 4] = [
    1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0,
];

#[test]
fn cbl_001_exact_envelope_constructs_the_pinned_generator() {
    let problem = CbcProblem::new(TIE_N, TIE_DIM).expect("fixture is structural");
    let estimate = problem.estimate().expect("fixture estimate is finite");
    assert_eq!(
        estimate.work_units(),
        6_309,
        "the independent schema-v4 work KAT moved"
    );
    let exact_budget = CbcBudget::new(estimate.work_units(), estimate.logical_state_bytes());

    let lattice = Lattice::try_cbc(TIE_N, TIE_DIM, exact_budget)
        .expect("the exact admitted envelope must complete");
    assert_eq!(lattice.n, TIE_N);
    assert_eq!(lattice.z, TIE_GENERATOR);
    let generator_bytes: Vec<u8> = lattice
        .z
        .iter()
        .flat_map(|component| component.to_le_bytes())
        .collect();
    assert_eq!(generator_bytes, TIE_GENERATOR_LE_BYTES);
}

#[test]
fn cbl_002_one_below_either_explicit_bound_refuses_before_execution() {
    let problem = CbcProblem::new(8, 3).expect("fixture is structural");
    let estimate = problem.estimate().expect("fixture estimate is finite");
    assert_eq!(estimate.work_units(), 3_699);

    let work_error = Lattice::try_cbc(
        8,
        3,
        CbcBudget::new(estimate.work_units() - 1, estimate.logical_state_bytes()),
    )
    .expect_err("one work unit below the sealed estimate must refuse");
    assert_eq!(
        work_error,
        CbcLatticeError::Admission(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.work_units(),
            available: estimate.work_units() - 1,
        })
    );

    let memory_error = Lattice::try_cbc(
        8,
        3,
        CbcBudget::new(estimate.work_units(), estimate.logical_state_bytes() - 1),
    )
    .expect_err("one byte below the sealed state envelope must refuse");
    assert_eq!(
        memory_error,
        CbcLatticeError::Admission(CbcAdmissionError::MemoryBudgetExceeded {
            required: estimate.logical_state_bytes(),
            available: estimate.logical_state_bytes() - 1,
        })
    );
}

#[test]
fn cbl_003_structural_refusals_remain_typed() {
    assert_eq!(
        Lattice::try_cbc(2, 0, CbcBudget::new(0, 0))
            .expect_err("point-count refusal must precede dimension and budget"),
        CbcLatticeError::Admission(CbcAdmissionError::InvalidPointCount { point_count: 2 })
    );
    assert_eq!(
        Lattice::try_cbc(3, 0, CbcBudget::new(0, 0))
            .expect_err("dimension refusal must precede budget"),
        CbcLatticeError::Admission(CbcAdmissionError::InvalidDimension { dimension: 0 })
    );
}

#[test]
fn cbl_004_default_and_exact_envelopes_are_generator_invariant() {
    let cases: &[(u32, usize, &[u32])] = &[(8, 3, &[1, 3, 1]), (127, 5, &[1, 29, 24, 56, 35])];
    for &(point_count, dimension, expected) in cases {
        let problem = CbcProblem::new(point_count, dimension).expect("fixture is structural");
        let estimate = problem.estimate().expect("fixture estimate is finite");
        let exact = Lattice::try_cbc(
            point_count,
            dimension,
            CbcBudget::new(estimate.work_units(), estimate.logical_state_bytes()),
        )
        .expect("exact envelope admits");
        let bounded_default = Lattice::cbc(point_count, dimension);

        assert_eq!(exact.n, bounded_default.n);
        assert_eq!(exact.z, bounded_default.z);
        assert_eq!(bounded_default.z.as_slice(), expected);
    }
}

#[test]
fn cbl_005_legacy_wrapper_has_a_finite_observable_ceiling() {
    assert_eq!(DEFAULT_CBC_BUDGET.max_work_units(), 1_000_000_000);
    assert_eq!(DEFAULT_CBC_BUDGET.max_memory_bytes(), 64 * 1024 * 1024);
    assert_ne!(DEFAULT_CBC_BUDGET, CbcBudget::UNBOUNDED);

    let documented = CbcProblem::new(1_031, 6)
        .expect("documented convergence fixture is structural")
        .estimate()
        .expect("documented convergence fixture has a finite estimate");
    assert!(
        documented.work_units() <= DEFAULT_CBC_BUDGET.max_work_units(),
        "the finite default must cover the documented 1031x6 work fixture"
    );
    assert!(
        documented.logical_state_bytes() <= DEFAULT_CBC_BUDGET.max_memory_bytes(),
        "the finite default must cover the documented 1031x6 state fixture"
    );

    // 5_503 is prime and schema v4 charges this two-dimensional case
    // 1_000_175_665 work units: just beyond the fixed compatibility ceiling,
    // while remaining adjacent enough to the documented fixture to avoid an
    // OOM-class regression test if the wrapper wiring ever moves.
    let outside_point_count = 5_503;
    let outside_dimension = 2;
    let problem = CbcProblem::new(outside_point_count, outside_dimension)
        .expect("outside-envelope fixture is structural");
    let estimate = problem
        .estimate()
        .expect("outside-envelope fixture estimate is finite");
    assert_eq!(estimate.work_units(), 1_000_175_665);
    assert!(
        estimate.work_units() > DEFAULT_CBC_BUDGET.max_work_units(),
        "the refusal fixture must stay beyond the finite compatibility ceiling"
    );
    assert_eq!(
        Lattice::try_cbc(outside_point_count, outside_dimension, DEFAULT_CBC_BUDGET)
            .expect_err("the finite default must refuse before executor allocation"),
        CbcLatticeError::Admission(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.work_units(),
            available: DEFAULT_CBC_BUDGET.max_work_units(),
        })
    );
    assert!(
        std::panic::catch_unwind(|| Lattice::cbc(outside_point_count, outside_dimension)).is_err(),
        "the compatibility wrapper must enforce the same finite envelope"
    );
}

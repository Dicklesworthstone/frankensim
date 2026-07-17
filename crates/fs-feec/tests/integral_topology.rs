//! G0/G4 first-tranche battery for I13.2b exact integral topology algebra.

#![cfg(feature = "moonshot-integral-topology")]

use fs_feec::integral_topology::{
    ExactAlgebraBudget, ExactIntegerMatrix, IntegralTopologyError, IntegralTopologyFailureClass,
    MatrixRole, SmithNormalFormWitness, SmithWitnessStage, TopologyApplicability,
    verify_smith_normal_form, verify_smith_normal_form_with_checkpoint,
};

fn budget(max_extent: usize, max_work: u128) -> ExactAlgebraBudget {
    let matrix_entries = max_extent * max_extent;
    ExactAlgebraBudget::new(
        max_extent,
        max_extent,
        matrix_entries,
        6 * matrix_entries,
        matrix_entries,
        max_work,
    )
}

fn matrix(rows: usize, cols: usize, entries: &[i128]) -> ExactIntegerMatrix {
    ExactIntegerMatrix::try_new(rows, cols, entries.to_vec(), ExactAlgebraBudget::default())
        .expect("fixture matrix")
}

fn identity(size: usize) -> ExactIntegerMatrix {
    let mut entries = vec![0_i128; size * size];
    for index in 0..size {
        entries[index * size + index] = 1;
    }
    matrix(size, size, &entries)
}

fn rank_one_witness() -> (ExactIntegerMatrix, SmithNormalFormWitness) {
    // Row 1 -= 2*row 0, then col 1 -= 2*col 0.
    let source = matrix(2, 2, &[2, 4, 4, 8]);
    let diagonal = matrix(2, 2, &[2, 0, 0, 0]);
    let left = matrix(2, 2, &[1, 0, -2, 1]);
    let left_inverse = matrix(2, 2, &[1, 0, 2, 1]);
    let right = matrix(2, 2, &[1, -2, 0, 1]);
    let right_inverse = matrix(2, 2, &[1, 2, 0, 1]);
    (
        source,
        SmithNormalFormWitness::new(diagonal, left, left_inverse, right, right_inverse),
    )
}

#[test]
fn it_001_complete_integer_witness_verifies_exactly() {
    let (source, witness) = rank_one_witness();
    let verified = verify_smith_normal_form(source.clone(), witness, budget(2, 48))
        .expect("complete exact witness");

    assert_eq!(verified.source(), &source);
    assert_eq!(verified.diagonal().entries(), &[2, 0, 0, 0]);
    assert_eq!(verified.invariant_factors(), &[2]);
    assert_eq!(verified.rank(), 1);
    assert_eq!(verified.scalar_operations(), 48);
    assert_eq!(
        verified.applicability(),
        TopologyApplicability::AbstractAlgebraOnly
    );
    assert_eq!(verified.left_transform().entries(), &[1, 0, -2, 1]);
    assert_eq!(verified.left_inverse().entries(), &[1, 0, 2, 1]);
    assert_eq!(verified.right_transform().entries(), &[1, -2, 0, 1]);
    assert_eq!(verified.right_inverse().entries(), &[1, 2, 0, 1]);
}

#[test]
fn it_002_canonical_diagonal_refuses_every_false_claim() {
    let identity3 = identity(3);
    for (diagonal, expected) in [
        (matrix(3, 3, &[2, 1, 0, 0, 4, 0, 0, 0, 0]), "off-diagonal"),
        (matrix(3, 3, &[-2, 0, 0, 0, 4, 0, 0, 0, 0]), "negative"),
        (matrix(3, 3, &[2, 0, 0, 0, 0, 0, 0, 0, 4]), "zero ordering"),
        (matrix(3, 3, &[4, 0, 0, 0, 6, 0, 0, 0, 0]), "divisibility"),
    ] {
        let witness = SmithNormalFormWitness::new(
            diagonal.clone(),
            identity3.clone(),
            identity3.clone(),
            identity3.clone(),
            identity3.clone(),
        );
        let error = verify_smith_normal_form(diagonal, witness, budget(3, 162))
            .expect_err("noncanonical diagonal must refuse");
        match expected {
            "off-diagonal" => assert!(matches!(
                error,
                IntegralTopologyError::OffDiagonalEntry { .. }
            )),
            "negative" => assert!(matches!(
                error,
                IntegralTopologyError::NegativeInvariantFactor { .. }
            )),
            "zero ordering" => assert!(matches!(
                error,
                IntegralTopologyError::NonzeroAfterZero { .. }
            )),
            "divisibility" => assert!(matches!(
                error,
                IntegralTopologyError::InvariantFactorDivisibility { .. }
            )),
            _ => unreachable!(),
        }
    }
}

#[test]
fn it_003_integer_inverse_is_required_even_when_zero_source_hides_transform() {
    let source = matrix(1, 1, &[0]);
    let witness = SmithNormalFormWitness::new(
        matrix(1, 1, &[0]),
        matrix(1, 1, &[2]),
        matrix(1, 1, &[1]),
        matrix(1, 1, &[1]),
        matrix(1, 1, &[1]),
    );
    let error = verify_smith_normal_form(source, witness, budget(1, 6))
        .expect_err("non-unimodular transform must be refuted");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Refuted);
    assert!(matches!(
        error,
        IntegralTopologyError::WitnessProductMismatch {
            stage: SmithWitnessStage::LeftTimesInverse,
            expected: 1,
            actual: 2,
            ..
        }
    ));
}

#[test]
fn it_004_transform_mutation_cannot_rebind_the_source() {
    let (source, witness) = rank_one_witness();
    let corrupted = SmithNormalFormWitness::new(
        matrix(2, 2, &[2, 0, 0, 2]),
        witness.left().clone(),
        witness.left_inverse().clone(),
        witness.right().clone(),
        witness.right_inverse().clone(),
    );
    assert!(matches!(
        verify_smith_normal_form(source, corrupted, budget(2, 48)),
        Err(IntegralTopologyError::WitnessProductMismatch {
            stage: SmithWitnessStage::DiagonalTransform,
            row: 1,
            col: 1,
            expected: 2,
            actual: 0,
        })
    ));
}

#[test]
fn it_005_shape_and_retained_storage_are_preflighted() {
    let (source, witness) = rank_one_witness();
    let wrong_left = SmithNormalFormWitness::new(
        witness.diagonal().clone(),
        matrix(1, 1, &[1]),
        witness.left_inverse().clone(),
        witness.right().clone(),
        witness.right_inverse().clone(),
    );
    assert!(matches!(
        verify_smith_normal_form(source.clone(), wrong_left, budget(2, 48)),
        Err(IntegralTopologyError::WitnessShape {
            role: MatrixRole::LeftTransform,
            expected_rows: 2,
            expected_cols: 2,
            actual_rows: 1,
            actual_cols: 1,
        })
    ));

    let tight_retained = ExactAlgebraBudget::new(2, 2, 4, 23, 4, 48);
    assert!(matches!(
        verify_smith_normal_form(source, witness, tight_retained),
        Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: 24,
            max: 23,
        })
    ));
}

#[test]
fn it_006_exact_work_cap_and_limit_plus_one_refuse_before_execution() {
    let (source, witness) = rank_one_witness();
    let verified = verify_smith_normal_form(source.clone(), witness.clone(), budget(2, 48))
        .expect("exact scalar-work cap admits");
    assert_eq!(verified.scalar_operations(), 48);
    assert!(matches!(
        verify_smith_normal_form(source, witness, budget(2, 47)),
        Err(IntegralTopologyError::ScalarWorkBudgetExceeded {
            requested: 48,
            max: 47,
        })
    ));
}

#[test]
fn it_007_every_cancellation_poll_is_transactional() {
    let (source, witness) = rank_one_witness();
    let mut poll_count = 0_usize;
    let verified = verify_smith_normal_form_with_checkpoint(
        source.clone(),
        witness.clone(),
        budget(2, 48),
        &mut |_| {
            poll_count += 1;
            true
        },
    )
    .expect("uninterrupted witness");
    assert_eq!(verified.scalar_operations(), 48);

    for stop_at in 0..poll_count {
        let mut observed = 0_usize;
        let result = verify_smith_normal_form_with_checkpoint(
            source.clone(),
            witness.clone(),
            budget(2, 48),
            &mut |_| {
                let keep_running = observed != stop_at;
                observed += 1;
                keep_running
            },
        );
        let error = result.expect_err("every injected cancellation must refuse");
        assert!(
            matches!(error, IntegralTopologyError::Cancelled { .. }),
            "poll {stop_at} must publish only cancellation: {error:?}"
        );
        if stop_at + 1 == poll_count {
            assert!(matches!(
                error,
                IntegralTopologyError::Cancelled {
                    phase: "smith witness finalize",
                    completed_scalar_operations: 48,
                    planned_scalar_operations: 48,
                }
            ));
        }
    }
}

#[test]
fn it_008_checked_arithmetic_refuses_overflow() {
    let shear = i128::MAX;
    let source = matrix(2, 2, &[0, 0, 2, 0]);
    let witness = SmithNormalFormWitness::new(
        matrix(2, 2, &[0, 0, 0, 0]),
        matrix(2, 2, &[1, shear, 0, 1]),
        matrix(2, 2, &[1, -shear, 0, 1]),
        identity(2),
        identity(2),
    );
    let error = verify_smith_normal_form(source, witness, budget(2, 48))
        .expect_err("coefficient explosion must remain unknown");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
    assert!(matches!(
        error,
        IntegralTopologyError::ArithmeticOverflow {
            stage: SmithWitnessStage::LeftTimesSource,
            row: 0,
            col: 0,
            term: 1,
        }
    ));
}

#[test]
fn it_009_empty_exact_complex_is_not_forced_nonempty() {
    let empty = matrix(0, 0, &[]);
    let witness = SmithNormalFormWitness::new(
        empty.clone(),
        empty.clone(),
        empty.clone(),
        empty.clone(),
        empty.clone(),
    );
    let verified =
        verify_smith_normal_form(empty, witness, ExactAlgebraBudget::new(0, 0, 0, 0, 0, 0))
            .expect("zero-dimensional exact algebra is valid");
    assert_eq!(verified.rank(), 0);
    assert!(verified.invariant_factors().is_empty());
    assert_eq!(verified.scalar_operations(), 0);
}

#[test]
fn it_010_matrix_extent_entry_count_and_entry_budget_refuse_exactly() {
    let small = ExactAlgebraBudget::new(2, 2, 4, 24, 4, 48);
    assert!(matches!(
        ExactIntegerMatrix::try_new(3, 1, vec![0; 3], small),
        Err(IntegralTopologyError::MatrixExtentExceeded { .. })
    ));
    assert!(matches!(
        ExactIntegerMatrix::try_new(2, 2, vec![0; 3], small),
        Err(IntegralTopologyError::MatrixEntryCount {
            expected: 4,
            actual: 3,
            ..
        })
    ));
    let entry_tight = ExactAlgebraBudget::new(2, 2, 3, 24, 4, 48);
    assert!(matches!(
        ExactIntegerMatrix::try_new(2, 2, vec![0; 4], entry_tight),
        Err(IntegralTopologyError::MatrixEntryBudgetExceeded {
            requested: 4,
            max: 3,
        })
    ));
}

#[test]
fn it_011_rectangular_and_empty_shapes_remain_exactly_distinct() {
    let rectangular = matrix(2, 3, &[2, 0, 0, 0, 6, 0]);
    let witness = SmithNormalFormWitness::new(
        rectangular.clone(),
        identity(2),
        identity(2),
        identity(3),
        identity(3),
    );
    let verified = verify_smith_normal_form(rectangular, witness, budget(3, 100))
        .expect("rectangular canonical Smith matrix");
    assert_eq!(verified.source().rows(), 2);
    assert_eq!(verified.source().cols(), 3);
    assert_eq!(verified.invariant_factors(), &[2, 6]);
    assert_eq!(verified.scalar_operations(), 100);

    let zero_by_three = matrix(0, 3, &[]);
    let witness = SmithNormalFormWitness::new(
        zero_by_three.clone(),
        matrix(0, 0, &[]),
        matrix(0, 0, &[]),
        identity(3),
        identity(3),
    );
    let verified = verify_smith_normal_form(zero_by_three, witness, budget(3, 54))
        .expect("zero-by-three matrix keeps its shape");
    assert_eq!((verified.source().rows(), verified.source().cols()), (0, 3));
    assert_eq!(verified.scalar_operations(), 54);

    let three_by_zero = matrix(3, 0, &[]);
    let witness = SmithNormalFormWitness::new(
        three_by_zero.clone(),
        identity(3),
        identity(3),
        matrix(0, 0, &[]),
        matrix(0, 0, &[]),
    );
    let verified = verify_smith_normal_form(three_by_zero, witness, budget(3, 54))
        .expect("three-by-zero matrix keeps its shape");
    assert_eq!((verified.source().rows(), verified.source().cols()), (3, 0));
    assert_eq!(verified.scalar_operations(), 54);
}

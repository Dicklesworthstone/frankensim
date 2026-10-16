//! G0 admission and legacy-projection coverage for dense Jacobi eigensolves.
//!
//! These tests exercise only allocation-free shape/work admission and the
//! low-level kernel's panic projection. They do not allocate matrices near the
//! practical cap and make no memory-availability, runtime, accuracy,
//! cancellation, or scientific-certificate claim.

#![deny(unsafe_code)]

use fs_la::eigen::{
    JACOBI_EIGH_ADMISSION_SCHEMA_VERSION, JacobiEighAdmissionError, MAX_EIGEN_WORK_ELEMENTS,
    admit_jacobi_eigh, jacobi_eigh,
};

fn panic_message(payload: Box<dyn core::any::Any + Send>) -> String {
    payload
        .downcast_ref::<&str>()
        .map(|message| (*message).to_owned())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "non-string panic payload".to_owned())
}

/// G0: 4095 is the exact largest dimension satisfying
/// `4 * n * n + 3 * n <= 64 Mi` elements, and every sealed receipt field is
/// an exact projection of that admission decision.
#[test]
fn exact_cap_boundary_returns_a_complete_receipt() {
    let admission = admit_jacobi_eigh(4_095).expect("4095 is the exact admitted boundary");

    assert_eq!(
        admission.schema_version(),
        JACOBI_EIGH_ADMISSION_SCHEMA_VERSION
    );
    assert_eq!(admission.dimension(), 4_095);
    assert_eq!(admission.matrix_entries(), 16_769_025);
    assert_eq!(admission.aggregate_work_elements(), 67_088_385);
    assert_eq!(admission.work_element_cap(), MAX_EIGEN_WORK_ELEMENTS);
    assert_eq!(MAX_EIGEN_WORK_ELEMENTS, 67_108_864);
}

/// G0: the first dimension beyond the practical cap is a typed refusal. The
/// legacy live kernel must project that exact authority before checking input
/// length or allocating, so an empty slice safely detects guard drift.
#[test]
fn first_cap_refusal_is_shared_by_the_live_kernel_guard() {
    let expected = JacobiEighAdmissionError::WorkCapExceeded {
        dimension: 4_096,
        required_elements: 67_121_152,
        cap_elements: MAX_EIGEN_WORK_ELEMENTS,
    };
    assert_eq!(admit_jacobi_eigh(4_096), Err(expected));

    let panic = std::panic::catch_unwind(|| jacobi_eigh(&[], 4_096))
        .expect_err("the legacy kernel must project typed admission as a panic");
    let message = panic_message(panic);
    assert!(
        message.contains("Jacobi aggregate workspace must fit the practical work cap"),
        "unexpected legacy panic: {message}"
    );
    assert!(
        message.contains(&format!("{expected:?}")),
        "legacy guard did not project the typed refusal: {message}"
    );
    assert!(
        !message.contains("a must be n*n"),
        "length checking ran before admission: {message}"
    );
}

/// G0: square-shape overflow and later aggregate-work overflow are separate
/// typed failures on every supported even-width `usize` target.
#[test]
fn shape_and_aggregate_arithmetic_overflows_are_distinct() {
    let first_square_overflow = 1usize << (usize::BITS / 2);
    assert_eq!(
        admit_jacobi_eigh(first_square_overflow),
        Err(JacobiEighAdmissionError::MatrixShapeOverflow {
            dimension: first_square_overflow,
        })
    );

    let aggregate_overflow = first_square_overflow - 1;
    assert!(
        aggregate_overflow.checked_mul(aggregate_overflow).is_some(),
        "fixture must retain a representable square"
    );
    assert_eq!(
        admit_jacobi_eigh(aggregate_overflow),
        Err(JacobiEighAdmissionError::AggregateWorkOverflow {
            dimension: aggregate_overflow,
        })
    );
}

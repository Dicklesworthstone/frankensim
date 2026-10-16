//! G0/G3 canonical-representative coverage for fs-opt SO(3) retractions.
//!
//! This target covers only successful public retraction outputs. It does not
//! claim that arbitrary input or persisted quaternion points are canonical,
//! nor does it cover projection, transport, curve derivatives, Stiefel, or
//! fs-ascent consumer migration.

#![deny(unsafe_code)]

use fs_opt::{Manifold, OptError};

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn assert_close(actual: f64, expected: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= 2e-15,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g0_g3_antipodal_noncommuting_retractions_are_bit_identical() {
    let manifold = Manifold::So3;
    let inverse_sqrt_two = std::f64::consts::FRAC_1_SQRT_2;
    let base = [inverse_sqrt_two, inverse_sqrt_two, 0.0, 0.0];
    let antipode = base.map(|component| -component);
    let body_y_quarter_turn = [0.0, std::f64::consts::FRAC_PI_2, 0.0];

    let canonical = manifold
        .retract(&base, &body_y_quarter_turn)
        .expect("right-composed SO(3) retraction");
    let antipodal = manifold
        .retract(&antipode, &body_y_quarter_turn)
        .expect("antipodal base denotes the same rotation");

    assert_eq!(bits(&canonical), bits(&antipodal));
    for (actual, expected) in canonical.iter().zip([0.5; 4]) {
        assert_close(*actual, expected);
    }
}

#[test]
fn g0_lexicographic_tie_planes_choose_positive_first_nonzero_axis() {
    let manifold = Manifold::So3;
    let zero_step = [0.0; 3];
    let canonical_representatives = [
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];

    for expected in canonical_representatives {
        let antipode = expected.map(|component| -component);
        for base in [expected, antipode] {
            let canonical = manifold
                .retract(&base, &zero_step)
                .expect("unit 180-degree quaternion");
            assert_eq!(bits(&canonical), bits(&expected));
        }
    }
}

#[test]
fn g0_zero_step_coalesces_general_antipodes() {
    let manifold = Manifold::So3;
    let base = [0.5, -0.5, 0.5, -0.5];
    let antipode = base.map(|component| -component);

    let canonical = manifold.retract(&base, &[0.0; 3]).expect("unit base");
    let antipodal = manifold
        .retract(&antipode, &[0.0; 3])
        .expect("unit antipode");

    assert_eq!(bits(&canonical), bits(&antipodal));
    assert_eq!(bits(&canonical), bits(&base));
}

#[test]
fn g0_successful_output_normalizes_every_signed_zero_lane() {
    let canonical = Manifold::So3
        .retract(&[1.0, -0.0, 0.0, -0.0], &[-0.0, 0.0, -0.0])
        .expect("identity rotation with signed-zero spelling");
    assert_eq!(bits(&canonical), bits(&[1.0, 0.0, 0.0, 0.0]));
}

#[test]
fn g0_existing_length_and_nonfinite_refusals_keep_their_provenance() {
    let manifold = Manifold::So3;
    let point = [1.0, 0.0, 0.0, 0.0];
    assert!(matches!(
        manifold.retract(&point[..3], &[0.0; 3]),
        Err(OptError::RetractionLen {
            input: "retraction point",
            expected: 4,
            got: 3,
        })
    ));
    assert!(matches!(
        manifold.retract(&point, &[0.0; 4]),
        Err(OptError::RetractionLen {
            input: "retraction step",
            expected: 3,
            got: 4,
        })
    ));

    let quiet_nan = f64::from_bits(0x7ff8_0000_0000_0022);
    assert!(matches!(
        manifold.retract(&point, &[0.0, quiet_nan, 0.0]),
        Err(OptError::RetractionNonFinite {
            input: "retraction step",
            component: 1,
            bits,
        }) if bits == quiet_nan.to_bits()
    ));
}

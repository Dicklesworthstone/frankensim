//! Force range and original trace-carrier regressions with independent oracles.

use super::*;
use fs_thermochem::SpeciesId;

#[test]
fn tiny_distance_recovers_force_after_an_underflowing_trace_product() {
    // Algebraically one trace weight cancels distance. The other must survive
    // even though multiplying the two weights first is below binary64 range.
    let trace = f64::from_bits((323_u64) << 52); // exactly 2^-700
    assert_eq!(trace * trace, 0.0);
    let expected = 0.75 * trace;
    let actual = solve::pair_force(0.75, trace, trace, trace).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(
        solve::pair_force(-0.75, trace, trace, trace).unwrap(),
        -expected
    );
}

#[test]
fn subnormal_mobility_and_distance_cancel_before_the_last_rounding() {
    let trace = f64::from_bits(2);
    let contrast = 2.0_f64.ln();
    let old_order = ((contrast * trace) * 1.0) / trace;
    assert_eq!(old_order, 0.5);
    // Since weight/distance=1 exactly, the original contrast is the oracle.
    let actual = solve::pair_force(contrast, trace, 1.0, trace).unwrap();
    assert_eq!(actual, contrast);
    assert_eq!(
        solve::pair_force(-contrast, trace, 1.0, trace).unwrap(),
        -contrast
    );
}

#[test]
fn truly_unrepresentable_force_refuses_without_suppressing_a_response() {
    let tiny = f64::from_bits(1);
    assert!(matches!(
        solve::pair_force(0.5, tiny, tiny, 1.0),
        Err(MulticomponentMiscibilityError::Arithmetic(_))
    ));
    assert!(matches!(
        solve::pair_force(f64::MAX, 1.0, 1.0, tiny),
        Err(MulticomponentMiscibilityError::Arithmetic(_))
    ));
    assert_eq!(solve::pair_force(0.0, tiny, tiny, tiny).unwrap(), 0.0);
}

#[test]
fn rounded_unit_carrier_retains_the_original_trace_chemical_potential() {
    let species = ["A", "B", "C"]
        .map(|id| SpeciesId::new(id).unwrap())
        .to_vec();
    let model = MulticomponentRegularSolutionV1::new(
        species,
        vec![0.0; 9],
        vec![0.0; 9],
        Temperature::new(250.0),
        Temperature::new(700.0),
    )
    .unwrap();
    let temperature = Temperature::new(300.0);
    for (left, right) in [
        ([1.0, 1e-20, 2e-20], [1.0, 2e-20, 4e-20]),
        ([0.7, 0.2, 0.1], [0.4, 0.4, 0.2]),
        ([0.7, 0.2, 0.1], [0.70000000000001, 0.19999999999999, 0.1]),
    ] {
        let old = model.evaluate(temperature, &left, || false).unwrap();
        let new = model.evaluate(temperature, &right, || false).unwrap();
        let expected = new.reduced_chemical_potentials()[0].unwrap()
            - old.reduced_chemical_potentials()[0].unwrap();
        let actual =
            solve::ideal_chemical_jump(old.mole_fractions(), new.mole_fractions(), 0).unwrap();
        if left[0] == 1.0 {
            assert_ne!(actual, 0.0);
            assert!((actual - expected).abs() <= 2e-34);
            assert!((actual + 3e-20).abs() <= 2e-34);
        } else {
            assert!((actual - expected).abs() <= 5e-16);
        }
        let reverse =
            solve::ideal_chemical_jump(new.mole_fractions(), old.mole_fractions(), 0).unwrap();
        assert_eq!(reverse, -actual);
    }
}

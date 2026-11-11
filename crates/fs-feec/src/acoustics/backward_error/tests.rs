use super::row_scale;

// The ratio a solver forms for one row.
fn ratio(defect: f64, terms: f64, norm: f64, limit: f64) -> f64 {
    defect / row_scale(terms, norm, limit)
}

// The three rows that the gates refused in the tests of the first run of the fs-feec lib target, measured: (defect, terms, largest
// row scale, limit). The largest row scale stands for the row's normwise scale N_i (N_i is at least the row's own terms).
const MEASURED: [(f64, f64, f64, f64); 3] = [
    (
        4.8127191387251314e-34,
        1.6355755795238293e-33,
        0.447213595499958,
        5e-13,
    ),
    (
        3.0803711672942303e-23,
        7.357340716603305e-22,
        1.1113755284244202e-5,
        5e-12,
    ),
    (
        5.758411920357775e-54,
        1.9796252459503064e-48,
        2.0281292592277076,
        5e-13,
    ),
];

#[test]
fn the_rows_that_were_refused_for_dust_pass_and_were_far_from_the_limit() {
    for (defect, terms, norm, limit) in MEASURED {
        // The componentwise ratio (the old gate) is above the limit by orders of magnitude.
        assert!(defect / terms > 1e3 * limit, "{defect:e} {terms:e}");
        let mixed = ratio(defect, terms, norm, limit);
        assert!(mixed < 0.05 * limit, "{mixed:e} against {limit:e}");
    }
}

// Planted negative 1: a genuinely wrong value in a zero row (no terms) is refused as soon as it exceeds eps * N_i; 1e-6 N_i is
// refused by orders of magnitude, three times eps N_i is refused, and half of eps N_i (roundoff) is not.
#[test]
fn a_wrong_value_in_a_zero_row_is_refused() {
    let limit = 5e-13;
    assert!(ratio(1e-6, 0.0, 1.0, limit) > 1e3 * limit);
    assert!(ratio(3.0 * f64::EPSILON, 0.0, 1.0, limit) > limit);
    assert!(ratio(0.5 * f64::EPSILON, 0.0, 1.0, limit) < limit);
}

// Planted negative 2: a data row (terms of the order of the normwise scale) is judged componentwise as before: 1e-9 relative is
// refused, 1e-14 relative is accepted, and the floor moves the threshold by less than 0.1 % there.
#[test]
fn a_one_part_in_a_billion_error_in_a_data_row_is_refused() {
    let limit = 5e-13;
    assert!(ratio(1e-9, 1.0, 1.0, limit) > limit);
    assert!(ratio(1e-14, 1.0, 1.0, limit) < limit);
    assert!((row_scale(1.0, 1.0, limit) - 1.0).abs() < 1e-3);
}

// A row with no terms and no scale at all (an empty row with exact data 0) has the denominator 0: a nonzero defect there is
// infinite, as in the gates before.
#[test]
fn an_empty_row_has_no_floor() {
    assert_eq!(row_scale(0.0, 0.0, 5e-13), 0.0);
}

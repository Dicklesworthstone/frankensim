use super::*;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn changing_mass_kinetic_increment_matches_independent_totals_in_moderate_frames() {
    for m1 in [0.5, 1.0, 3.0] {
        let m0 = 2.0;
        let p0 = [20.0, -6.0, 4.0];
        let p1 = [15.0, 2.0, -1.0];
        for frame in [[0.0; 3], [10.0, -4.0, 3.0]] {
            let direct = |m: f64, p: [f64; 3]| {
                (0..3)
                    .map(|d| 0.5 * m * (p[d] / m - frame[d]).powi(2))
                    .sum::<f64>()
            };
            near(
                kinetic_change(m0, m1, p0, p1, frame).unwrap(),
                direct(m1, p1) - direct(m0, p0),
                1e-12,
            );
        }
    }
}

#[test]
fn common_drift_does_not_hide_a_small_relative_kinetic_change() {
    let drift = 2.0_f64.powi(40);
    let m0 = 2.0;
    let m1 = 2.5;
    let p0 = [m0 * drift, 0.0, 0.0];
    let p1 = [m1 * drift + 1.0, 2.0, 0.0];
    near(
        kinetic_change(m0, m1, p0, p1, [drift, 0.0, 0.0]).unwrap(),
        1.0,
        1e-15,
    );
    assert!(kinetic_change(m0, 0.0, p0, p1, [0.0; 3]).is_err());
}

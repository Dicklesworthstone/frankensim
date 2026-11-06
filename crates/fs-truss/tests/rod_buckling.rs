//! G1 regression: convergence back to a straight compressed rod is not stability.

use fs_truss::rod_buckling_check;

#[test]
fn g1_zero_and_tiny_initial_bows_cannot_hide_supercritical_compression() {
    let length: f64 = 1.0;
    let area: f64 = 0.01;
    let youngs: f64 = 2e7;
    let inertia = area * area / 12.0;
    let euler = std::f64::consts::PI.powi(2) * youngs * inertia / (4.0 * length * length);
    for bow in [0.0, 1e-12, 2e-3] {
        let (stable, deflection) = rod_buckling_check(length, area, youngs, 0.5 * euler, 1.4, bow);
        assert!(
            stable && deflection < 0.05,
            "subcritical bow={bow:e}, deflection={deflection:e}"
        );
        let (stable, deflection) = rod_buckling_check(length, area, youngs, 0.5 * euler, 2.6, bow);
        assert!(
            !stable,
            "supercritical bow={bow:e}, deflection={deflection:e}"
        );
    }
}

#[test]
fn g0_invalid_compression_inputs_refuse() {
    for args in [
        [0.0, 0.01, 2e7, 1.0, 1.0, 0.0],
        [1.0, -0.01, 2e7, 1.0, 1.0, 0.0],
        [1.0, 0.01, f64::NAN, 1.0, 1.0, 0.0],
        [1.0, 0.01, 2e7, 1.0, 0.0, 0.0],
        [1.0, 0.01, 2e7, 1.0, 1.0, -0.01],
        [1.0, 0.01, 2e7, f64::INFINITY, 1.0, 0.0],
    ] {
        assert_eq!(
            rod_buckling_check(args[0], args[1], args[2], args[3], args[4], args[5]),
            (false, f64::INFINITY)
        );
    }
}

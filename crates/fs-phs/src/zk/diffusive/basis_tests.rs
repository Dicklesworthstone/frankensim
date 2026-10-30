use super::*;

#[test]
fn reusable_basis_preserves_original_coefficients_and_errors_bitwise() {
    for poles in [2, 8, 16] {
        let basis = ZkDiffusiveBasis::new(poles).unwrap();
        assert_eq!(basis.poles_per_field(), poles);
        for [radius, rho, mu, gamma, pr] in [
            [0.001, 1.2, 1.8e-5, 1.4, 0.71],
            [0.002, 0.8, 2.1e-5, 1.25, 0.93],
            [0.0002, 0.15, 1.9e-5, 5.0 / 3.0, 2.0 / 3.0],
            [0.001, 1.2, 1.8e-5, 1.0, 0.71],
        ] {
            let band = [0.1, 10.0];
            let direct = ZkDiffusive::new(radius, rho, mu, gamma, pr, poles, band, 0.02).unwrap();
            let cached = basis.bind(radius, rho, mu, gamma, pr, band, 0.02).unwrap();
            assert_eq!(cached, direct);
            for omega in [0.1, 1.0, 10.0] {
                assert_eq!(
                    cached.response_at(omega).unwrap(),
                    direct.response_at(omega).unwrap()
                );
            }
        }
    }
}

#[test]
fn each_basis_binding_keeps_its_own_physical_band_admission() {
    let basis = ZkDiffusiveBasis::new(2).unwrap();
    let slow = basis
        .bind(1.0, 1.0, 1.0, 1.4, 0.71, [0.1, 1.0], 1e-4)
        .unwrap();
    // The same positive measure does NOT admit an unresolved wider band.
    assert!(
        basis
            .bind(1.0, 1.0, 1.0, 1.4, 0.71, [1.0, 900.0], 1e-4)
            .is_err()
    );
    // Halving rho halves each diffusion time; at twice the frequency the
    // normalized operators scale by two. No eigenproblem depends on rho.
    let fast = basis
        .bind(1.0, 0.5, 1.0, 1.4, 0.71, [0.2, 2.0], 1e-4)
        .unwrap();
    let (zs, ys) = slow.response_at(0.5).unwrap();
    let (zf, yf) = fast.response_at(1.0).unwrap();
    assert_eq!(zf, zs.scale(2.0));
    assert_eq!(yf, ys.scale(2.0));
    for poles in [0, MAX_ZK_DIFFUSIVE_POLES + 1] {
        assert!(ZkDiffusiveBasis::new(poles).is_err());
    }
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            basis
                .bind(invalid, 1.0, 1.0, 1.4, 0.71, [0.1, 1.0], 0.1)
                .is_err()
        );
    }
}

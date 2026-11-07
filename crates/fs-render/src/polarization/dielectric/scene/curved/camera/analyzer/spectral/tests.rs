use super::*;

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);

#[test]
fn g1_spectral_analyzer_interpolates_absolute_rows_and_keeps_authored_knots() {
    let blue = Stokes {
        i: 0.4,
        q: 0.2,
        u: 0.0,
        v: 0.1,
    };
    let red = Stokes {
        i: 0.6,
        q: -0.1,
        u: 0.0,
        v: -0.2,
    };
    let spectrum = SpectralCameraAnalyzer::try_new(X, vec![(400.0, blue), (600.0, red)]).unwrap();
    assert_eq!(spectrum.wavelength_bounds_nm(), (400.0, 600.0));
    assert_eq!(spectrum.at(400.0).unwrap().row(), blue);
    assert_eq!(spectrum.at(600.0).unwrap().row(), red);
    let mid = spectrum.at(500.0).unwrap();
    assert_eq!(mid.reference_axis(), X);
    for (a, b) in [mid.row().i, mid.row().q, mid.row().u, mid.row().v]
        .into_iter()
        .zip([0.5, 0.05, 0.0, -0.05])
    {
        assert!((a - b).abs() < 2e-16);
    }
    // A circular source must see the wavelength-specific absolute response.
    // Normalizing either row before interpolation would change these powers.
    for (w, expected) in [(400.0, 0.5), (500.0, 0.45), (600.0, 0.4)] {
        let a = spectrum.at(w).unwrap().row();
        assert!((a.i + a.v - expected).abs() < 2e-16);
    }
}

#[test]
fn g0_spectral_analyzer_retains_physical_boundary_rows_and_constant_bits() {
    for angle in [0.0, 0.1, 0.7, 1.2, 1.5707963267948966] {
        let a = CameraAnalyzer::try_linear(X, angle).unwrap().row();
        let b = CameraAnalyzer::try_linear(X, angle + 1e-7).unwrap().row();
        let spectrum = SpectralCameraAnalyzer::try_new(X, vec![(380.0, a), (780.0, b)]).unwrap();
        for i in 0..=32 {
            let row = spectrum.at(380.0 + 12.5 * f64::from(i)).unwrap().row();
            let degree = det::hypot(det::hypot(row.q, row.u), row.v);
            assert!(degree <= row.i.min(1.0 - row.i));
        }
        let constant = SpectralCameraAnalyzer::try_new(X, vec![(380.0, a), (780.0, a)]).unwrap();
        assert_eq!(constant.at(527.125).unwrap().row(), a);
    }
    let tiny = Stokes::unpolarized(f64::from_bits(1));
    let spectrum = SpectralCameraAnalyzer::try_new(X, vec![(380.0, tiny), (780.0, tiny)]).unwrap();
    assert_eq!(spectrum.at(580.0).unwrap().row(), tiny);
}

#[test]
fn spectral_analyzer_refuses_missing_support_active_rows_and_invalid_axes() {
    let row = Stokes::unpolarized(0.5);
    let spectrum = SpectralCameraAnalyzer::try_new(X, vec![(400.0, row), (700.0, row)]).unwrap();
    for w in [0.0, -1.0, 399.0, 701.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            spectrum.at(w),
            Err(CameraRadianceError::InvalidPolarization)
        );
    }
    assert!(spectrum.admit_band(400.0, 700.0).is_ok());
    assert!(spectrum.admit_band(380.0, 780.0).is_err());
    assert!(spectrum.admit_band(600.0, 500.0).is_err());
    for knots in [
        vec![],
        vec![(500.0, row)],
        vec![(600.0, row), (500.0, row)],
        vec![(500.0, row), (500.0, row)],
        vec![(400.0, row), (700.0, Stokes::unpolarized(0.0))],
        vec![
            (400.0, row),
            (
                700.0,
                Stokes {
                    i: 0.9,
                    q: 0.2,
                    u: 0.0,
                    v: 0.0,
                },
            ),
        ],
    ] {
        assert!(SpectralCameraAnalyzer::try_new(X, knots).is_err());
    }
    for axis in [Vec3::new(0.0, 0.0, 0.0), Vec3::new(f64::NAN, 1.0, 0.0)] {
        assert!(SpectralCameraAnalyzer::try_new(axis, vec![(400.0, row), (700.0, row)]).is_err());
    }
}

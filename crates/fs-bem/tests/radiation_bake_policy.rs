//! G1/G3 radiation-load regressions on the original dense and iterative owners.
//! Small prescribed panel sources test equation selection and route agreement;
//! they do not qualify a resolved physical mouth or an instrument model.
use fs_bem::helmholtz::iterative::{IterativeError, IterativeOptions};
use fs_bem::helmholtz::{Formulation, HelmholtzError, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_bake::{bake_radiation_load, bake_radiation_load_iterative, lathe_profile};
use fs_bem::radiation_policy::GeometryPolicy;
use fs_math::c64::C64;
use std::f64::consts::PI;

#[test]
fn g1_mouth_area_does_not_set_the_body_fictitious_frequency() {
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let medium = Medium::air();
    let frequencies = [0.5 * medium.sound_speed, PI * medium.sound_speed];
    let whole = vec![true; surface.areas().len()];
    let mut patch = vec![false; surface.areas().len()];
    patch[0] = true;
    let whole = bake_radiation_load(&surface, &whole, &frequencies, medium, 4, "sphere").unwrap();
    let patch = bake_radiation_load(&surface, &patch, &frequencies, medium, 4, "patch").unwrap();

    // A large driven area is not itself a fictitious resonance. Conversely,
    // making the driven area small cannot remove the closed body's spectrum.
    assert!(whole.rows[0].mouth_ka > 0.5);
    assert!(patch.rows[1].mouth_ka < 0.5);
    for bake in [&whole, &patch] {
        assert_eq!(bake.rows[0].formulation, Formulation::PlainCbie);
        assert_eq!(bake.rows[1].formulation, Formulation::BurtonMiller);
        assert!(bake.rows.iter().all(|row| row.radiated_power > 0.0));
    }

    // ka=pi is the smooth sphere's first interior Dirichlet frequency. The
    // exterior breathing impedance is regular there; the unchanged BM
    // discretization must remain inside its existing 320-panel 8% band.
    let analytic =
        C64::new(PI * PI, -PI).scale(medium.density * medium.sound_speed / (1.0 + PI * PI));
    let error = (whole.rows[1].z_specific - analytic).abs() / analytic.abs();
    assert!(
        error < 0.08,
        "sphere impedance at interior resonance: {error}"
    );
}

#[test]
fn g3_iterative_mouth_load_preserves_low_frequency_triangle_accuracy() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let medium = Medium::air();
    let driven: Vec<_> = surface
        .normals()
        .iter()
        .map(|normal| normal[2] > 0.4)
        .collect();
    let frequencies = [0.05 * medium.sound_speed, 0.5 * medium.sound_speed];
    let dense = bake_radiation_load(&surface, &driven, &frequencies, medium, 2, "cap").unwrap();
    // One FMM leaf keeps this small fixture's triangle interactions direct;
    // the production GMRES/preconditioner still solves the original equation.
    let mut options = IterativeOptions {
        tol: 1e-11,
        leaf_capacity: 128,
        ..IterativeOptions::default()
    };
    let iterative =
        bake_radiation_load_iterative(&surface, &driven, &frequencies, medium, 2, "cap", &options)
            .unwrap();
    for (a, b) in dense.rows.iter().zip(&iterative.rows) {
        assert_eq!(a.formulation, Formulation::PlainCbie);
        assert_eq!(b.formulation, Formulation::PlainCbie);
        let error = (a.z_specific - b.z_specific).abs() / a.z_specific.abs();
        let resistance_error = (a.z_specific.re - b.z_specific.re).abs() / a.z_specific.re.abs();
        assert!(error < 1e-9, "complex load discrepancy: {error}");
        assert!(
            resistance_error < 1e-7,
            "radiation resistance discrepancy: {resistance_error}"
        );
        let norm: f64 = a.directivity.coefficients.iter().map(|z| z.norm_sq()).sum();
        let discrepancy: f64 = a
            .directivity
            .coefficients
            .iter()
            .zip(&b.directivity.coefficients)
            .map(|(x, y)| (*x - *y).norm_sq())
            .sum();
        assert!(discrepancy < 1e-16 * norm);
    }

    // The selected equation's own storage limit is still authoritative;
    // failure does not silently switch to a different acoustic formulation.
    options.max_near_entries = 1;
    assert!(matches!(
        bake_radiation_load_iterative(&surface, &driven, &frequencies, medium, 2, "cap", &options),
        Err(IterativeError::NearFieldBudget {
            maximum_entries: 1,
            ..
        })
    ));
}

#[test]
fn g0_bake_requires_actual_closed_geometry_for_spectral_selection() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let centroids = SpherePanels::new(
        surface.centroids().to_vec(),
        surface.normals().to_vec(),
        surface.areas().to_vec(),
    )
    .unwrap();
    let driven = vec![true; surface.areas().len()];
    assert!(matches!(
        bake_radiation_load(&centroids, &driven, &[100.0], Medium::air(), 2, "centroids"),
        Err(HelmholtzError::BadParameter { .. })
    ));
}

#[test]
fn g0_lathed_mouth_cap_shares_the_exact_wall_rim() {
    // The outer cap has three annuli. Computing (0.1 * 3) / 3 would
    // create a distinct radius and an open seam at the original wall.
    let bell = lathe_profile(&[(0.0, 0.01), (0.04, 0.1)], 8).unwrap();
    let surface = SpherePanels::from_triangles(bell.triangles).unwrap();
    let policy = GeometryPolicy::new(&surface).unwrap();
    assert_eq!(policy.first_dirichlet_bounds().len(), 1);
}

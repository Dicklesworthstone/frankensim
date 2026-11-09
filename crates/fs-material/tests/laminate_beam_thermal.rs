//! Analytical narrow-strip thermal bending; no measured-material validation.

use fs_material::laminate::thermal::beam::PreparedLaminateBeam;
use fs_material::laminate::thermal::{
    LaminaThermalStrain, LaminateThermalError, PreparedLaminateThermal,
};
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use std::sync::Arc;

fn isotropic(e: f64, nu: f64, h: f64, density: f64) -> LaminaLayer {
    let q = e / (1.0 - nu * nu);
    let g = e / (2.0 * (1.0 + nu));
    LaminaLayer {
        thickness_m: h,
        orientation_rad: 0.0,
        density_kg_m3: density,
        material: None,
        law: LayerLaw::PlaneStress {
            stiffness_pa: [q, nu * q, 0.0, nu * q, q, 0.0, 0.0, 0.0, g],
            transverse_shear_pa: [g; 2],
        },
    }
}

fn prepared(
    layers: Vec<LaminaLayer>,
    profiles: &[LaminaThermalStrain],
    reference: ReferenceSurface,
    width: f64,
    axis: f64,
) -> PreparedLaminateBeam {
    let source = Arc::new(
        PreparedLaminateThermal::try_new(
            Arc::new(LaminateSection::try_new(layers, reference).unwrap()),
            profiles,
            || false,
        )
        .unwrap(),
    );
    let beam = PreparedLaminateBeam::try_new(Arc::clone(&source), width, axis, || false).unwrap();
    assert!(core::ptr::eq(beam.source(), source.as_ref()));
    beam
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g1_homogeneous_thermal_gradient_uses_free_beam_stiffness_and_signed_curvature() {
    let (e, nu, h, width, mean, gradient) = (70.0e9, 0.29, 0.01, 0.02, 0.002, 0.04);
    let profile = LaminaThermalStrain::from_faces(
        [mean - gradient * h / 2.0, mean - gradient * h / 2.0, 0.0],
        [mean + gradient * h / 2.0, mean + gradient * h / 2.0, 0.0],
    )
    .unwrap();
    let beam = prepared(
        vec![isotropic(e, nu, h, 2700.0)],
        &[profile],
        ReferenceSurface::Midsurface,
        width,
        0.0,
    );
    let ea = e * width * h;
    let ei = e * width * h * h * h / 12.0;
    near(beam.tangent()[0][0], ea, 2e-8);
    near(beam.tangent()[1][1], ei, 2e-12);
    near(beam.tangent()[0][1], 0.0, 1e-10);
    near(beam.thermal_resultants()[0], ea * mean, 1e-10);
    near(beam.thermal_resultants()[1], -ei * gradient, 1e-12);
    let free = beam.free_generalized_strain();
    near(free[0], mean, 2e-18);
    near(free[1], -gradient, 2e-16);
    let response = beam.evaluate(free[0], free[1], || false).unwrap();
    near(response.axial_force_n, 0.0, 1e-10);
    near(response.bending_moment_n_m, 0.0, 1e-12);
    near(response.strain_energy_j_m, 0.0, 1e-25);
    let stress = beam.layer_stress(0, free[0], free[1]).unwrap();
    for value in stress.bottom_pa.into_iter().chain(stress.top_pa) {
        near(value, 0.0, 1e-7);
    }
    let restrained = beam.layer_stress(0, 0.0, 0.0).unwrap();
    near(restrained.bottom_pa[0], -e * profile.faces()[0][0], 1e-7);
    assert_eq!(&restrained.bottom_pa[1..], &[0.0, 0.0]);
    let lateral = beam.layer_strain_faces(0, 0.0, 0.0).unwrap();
    near(lateral[0][1], (1.0 + nu) * profile.faces()[0][1], 2e-18);
}

fn bimetal(nu: [f64; 2], reference: ReferenceSurface, axis: f64) -> PreparedLaminateBeam {
    prepared(
        vec![
            isotropic(70e9, nu[0], 0.004, 2700.0),
            isotropic(110e9, nu[1], 0.006, 8000.0),
        ],
        &[
            LaminaThermalStrain::uniform([0.001, 0.001, 0.0]).unwrap(),
            LaminaThermalStrain::uniform([0.004, 0.004, 0.0]).unwrap(),
        ],
        reference,
        0.03,
        axis,
    )
}

#[test]
fn g1_dissimilar_bimetal_has_coupled_curvature_and_nonzero_self_stress_energy() {
    let beam = bimetal([0.22, 0.38], ReferenceSurface::Midsurface, 0.0);
    let spans = [
        (70e9, -0.005_f64, -0.001_f64, 0.001),
        (110e9, -0.001_f64, 0.005_f64, 0.004),
    ];
    let (mut a, mut b, mut d, mut n, mut m) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (e, lo, hi, alpha) in spans {
        let q0 = hi - lo;
        let q1 = (hi * hi - lo * lo) / 2.0;
        let q2 = (hi.powi(3) - lo.powi(3)) / 3.0;
        a += 0.03 * e * q0;
        b -= 0.03 * e * q1;
        d += 0.03 * e * q2;
        n += 0.03 * e * alpha * q0;
        m -= 0.03 * e * alpha * q1;
    }
    let determinant = a * d - b * b;
    let eps = (d * n - b * m) / determinant;
    let kappa = (a * m - b * n) / determinant;
    let response = beam.evaluate(eps, kappa, || false).unwrap();
    near(beam.tangent()[0][1], b, 1e-10);
    near(beam.free_generalized_strain()[0], eps, 2e-18);
    near(beam.free_generalized_strain()[1], kappa, 2e-15);
    near(response.axial_force_n, 0.0, 1e-10);
    near(response.bending_moment_n_m, 0.0, 1e-12);
    let energy: f64 = spans
        .into_iter()
        .map(|(e, lo, hi, alpha)| {
            let residual = eps - alpha;
            0.5 * 0.03
                * e
                * (residual * residual * (hi - lo) - residual * kappa * (hi * hi - lo * lo)
                    + kappa * kappa * (hi.powi(3) - lo.powi(3)) / 3.0)
        })
        .sum();
    assert!(energy > 1.0);
    near(response.strain_energy_j_m, energy, 2e-12);
    for (layer, ((e, lo, hi, alpha), nu)) in spans.into_iter().zip([0.22, 0.38]).enumerate() {
        let stress = beam.layer_stress(layer, eps, kappa).unwrap();
        let strains = beam.layer_strain_faces(layer, eps, kappa).unwrap();
        for (face, z) in [lo, hi].into_iter().enumerate() {
            let elastic = eps - z * kappa - alpha;
            near(
                [stress.bottom_pa, stress.top_pa][face][0],
                e * elastic,
                2e-7,
            );
            near(strains[face][1], alpha - nu * elastic, 2e-18);
            assert_eq!(&[stress.bottom_pa, stress.top_pa][face][1..], &[0.0, 0.0]);
        }
    }
    // Free lateral stresses make this 1-D response independent of the two
    // Poisson ratios. A global affine plate relaxation would retain additional
    // compatible transverse self-stress, rather than this exact strip energy.
    let zero_nu = bimetal([0.0, 0.0], ReferenceSurface::Midsurface, 0.0);
    near(
        zero_nu
            .evaluate(eps, kappa, || false)
            .unwrap()
            .strain_energy_j_m,
        response.strain_energy_j_m,
        2e-12,
    );
    for i in 0..2 {
        for j in 0..2 {
            near(beam.tangent()[i][j], zero_nu.tangent()[i][j], 2e-8);
        }
    }
}

#[test]
fn g3_source_reference_and_beam_axis_shifts_preserve_physical_energy_stress_and_mass() {
    let axis = 0.0017;
    let shift = -0.0021;
    let a = bimetal([0.22, 0.38], ReferenceSurface::Midsurface, axis);
    let b = bimetal(
        [0.22, 0.38],
        ReferenceSurface::Offset(0.005 + shift),
        axis - shift,
    );
    let (eps, kappa) = (0.002, 0.13);
    let ra = a.evaluate(eps, kappa, || false).unwrap();
    let rb = b.evaluate(eps, kappa, || false).unwrap();
    near(ra.axial_force_n, rb.axial_force_n, 1e-10);
    near(ra.bending_moment_n_m, rb.bending_moment_n_m, 1e-12);
    near(ra.strain_energy_j_m, rb.strain_energy_j_m, 1e-12);
    for i in 0..3 {
        near(
            a.reference_mass_moments()[i],
            b.reference_mass_moments()[i],
            1e-15,
        );
    }
    for layer in 0..2 {
        let sa = a.layer_stress(layer, eps, kappa).unwrap();
        let sb = b.layer_stress(layer, eps, kappa).unwrap();
        for (x, y) in sa
            .bottom_pa
            .into_iter()
            .chain(sa.top_pa)
            .zip(sb.bottom_pa.into_iter().chain(sb.top_pa))
        {
            near(x, y, 2e-7);
        }
    }
    let delta = -0.001;
    let moved = bimetal([0.22, 0.38], ReferenceSurface::Midsurface, axis + delta);
    let rm = moved
        .evaluate(eps - delta * kappa, kappa, || false)
        .unwrap();
    near(rm.strain_energy_j_m, ra.strain_energy_j_m, 1e-12);
    near(
        rm.bending_moment_n_m,
        ra.bending_moment_n_m + delta * ra.axial_force_n,
        1e-12,
    );
    let centroid = a.source().section().i1_kg_m() / a.source().section().i0_kg_m2();
    let centered = bimetal([0.22, 0.38], ReferenceSurface::Midsurface, centroid);
    near(
        centered.reference_mass_moments()[0],
        0.03 * (2700.0 * 0.004 + 8000.0 * 0.006),
        1e-15,
    );
    near(centered.reference_mass_moments()[1], 0.0, 1e-18);
    assert!(centered.reference_mass_moments()[2] > 0.0);
}

#[test]
fn g0_planar_admission_is_per_layer_and_does_not_hide_balanced_shear_coupling() {
    let mut layers = vec![isotropic(90.0, 0.2, 0.5, 2.0); 2];
    for (index, layer) in layers.iter_mut().enumerate() {
        let sign = if index == 0 { 1.0 } else { -1.0 };
        layer.law = LayerLaw::PlaneStress {
            stiffness_pa: [
                90.0,
                20.0,
                sign * 6.0,
                20.0,
                50.0,
                sign * 3.0,
                sign * 6.0,
                sign * 3.0,
                18.0,
            ],
            transverse_shear_pa: [18.0; 2],
        };
    }
    let section = Arc::new(LaminateSection::try_new(layers, ReferenceSurface::Midsurface).unwrap());
    near(section.a_pa_m()[2], 0.0, 0.0);
    let profiles = [LaminaThermalStrain::uniform([0.001, 0.002, 0.0]).unwrap(); 2];
    let source = Arc::new(PreparedLaminateThermal::try_new(section, &profiles, || false).unwrap());
    assert!(PreparedLaminateBeam::try_new(source, 0.1, 0.0, || false).is_err());
    let section = Arc::new(
        LaminateSection::try_new(
            vec![isotropic(90.0, 0.2, 0.5, 2.0); 2],
            ReferenceSurface::Midsurface,
        )
        .unwrap(),
    );
    let shear = [
        LaminaThermalStrain::uniform([0.001, 0.002, 0.0001]).unwrap(),
        LaminaThermalStrain::uniform([0.001, 0.002, -0.0001]).unwrap(),
    ];
    let source = Arc::new(PreparedLaminateThermal::try_new(section, &shear, || false).unwrap());
    assert!(PreparedLaminateBeam::try_new(source, 0.1, 0.0, || false).is_err());
}

#[test]
fn g4_invalid_and_cancelled_trials_leave_original_section_available() {
    let beam = bimetal([0.22, 0.38], ReferenceSurface::Midsurface, 0.0);
    let source = Arc::new(beam.source().clone());
    for width in [0.0, -1.0, f64::INFINITY] {
        assert!(PreparedLaminateBeam::try_new(Arc::clone(&source), width, 0.0, || false).is_err());
    }
    assert!(matches!(
        PreparedLaminateBeam::try_new(Arc::clone(&source), 0.03, 0.0, || true),
        Err(LaminateThermalError::Cancelled)
    ));
    assert!(matches!(
        beam.evaluate(0.0, 0.0, || true),
        Err(LaminateThermalError::Cancelled)
    ));
    assert!(beam.evaluate(f64::NAN, 0.0, || false).is_err());
    assert!(beam.layer_stress(2, 0.0, 0.0).is_err());
    assert_eq!(beam.source(), source.as_ref());
    assert!(beam.evaluate(0.0, 0.0, || false).unwrap().strain_energy_j_m > 0.0);
}

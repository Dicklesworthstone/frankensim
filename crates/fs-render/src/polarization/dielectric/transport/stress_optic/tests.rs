use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass};
use crate::material_state::silica_refraction::corning7980_glass;
use crate::polarization::Mueller;
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_geom::{Point3, Vec3};

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2.0e-14,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn near_power(actual: Stokes, expected: Stokes) {
    near(actual.i, expected.i);
    near(actual.q, expected.q);
    near(actual.u, expected.u);
    near(actual.v, expected.v);
}
fn beam(wavelength: f64, direction: Vec3, power: Stokes) -> PolarizedBeam {
    // The scalar source temperature is explicit, separately from the stress
    // coefficient's unreported laboratory temperature; this is its reference use.
    let glass =
        corning7980_glass(295.0, BeerLambertAbsorption::try_constant(0.7).unwrap()).unwrap();
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        wavelength,
        glass,
        power,
    )
    .unwrap()
}

#[test]
fn tension_positive_source_matches_independent_signed_jones_and_original_beer() {
    let original = beam(
        633.0,
        Vec3::new(0.0, 0.0, 1.0),
        Stokes::try_new(2.0, 0.0, 2.0, 0.0).unwrap(),
    );
    let distance = 0.04;
    let intensity = 2.0 * (-0.7_f64 * distance).exp();
    for sign in [-1.0, 1.0] {
        let sigma = sign * 400_000.0;
        let out = original
            .advance_corning7980_stress_optic(&[sigma, 0.0, 0.0, 0.0, 0.0, 0.0], distance)
            .unwrap();
        // Independent Jones amplitudes E_h=exp(i*delta)/sqrt(2), E_v=1/sqrt(2).
        // Native V=2 Im(E_h conj(E_v)), so horizontal tensile stress gives +V.
        let delta = core::f64::consts::TAU * 3.5e-12 * sigma * distance / 633.0e-9;
        near_power(
            out.power(),
            Stokes {
                i: intensity,
                q: 0.0,
                u: intensity * delta.cos(),
                v: intensity * delta.sin(),
            },
        );
        assert_eq!(out.power().v.is_sign_positive(), sign > 0.0);
        near(out.power().degree(), 1.0);
        assert_eq!(out.ray().origin, Point3::new(0.0, 0.0, distance));
        assert_eq!(out.frame(), original.frame());
        assert_eq!(out.medium(), original.medium());
    }

    let horizontal = beam(
        633.0,
        Vec3::new(0.0, 0.0, 1.0),
        Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap(),
    );
    let tau = 300_000.0;
    let out = horizontal
        .advance_corning7980_stress_optic(&[0.0, 0.0, 0.0, tau, 0.0, 0.0], distance)
        .unwrap();
    let crossed = Mueller::linear_polarizer(core::f64::consts::FRAC_PI_2).apply(out.power());
    let delta = core::f64::consts::TAU * 3.5e-12 * (2.0 * tau) * distance / 633.0e-9;
    near(
        crossed.i,
        (-0.7_f64 * distance).exp() * (0.5 * delta).sin().powi(2),
    );
    assert!(crossed.i > 0.05);
}

#[test]
fn reference_rotation_and_reverse_ray_preserve_native_handedness() {
    let original = beam(
        633.0,
        Vec3::new(0.0, 0.0, 1.0),
        Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(),
    );
    let stress = [400_000.0, 0.0, 0.0, 70_000.0, 0.0, 0.0];
    let direct = original
        .advance_corning7980_stress_optic(&stress, 0.03)
        .unwrap();
    let (s, c) = 0.61_f64.sin_cos();
    let rotated = original
        .with_reference_axis(Vec3::new(c, s, 0.0))
        .unwrap()
        .advance_corning7980_stress_optic(&stress, 0.03)
        .unwrap()
        .with_reference_axis(Vec3::new(1.0, 0.0, 0.0))
        .unwrap();
    near_power(rotated.power(), direct.power());

    // For reverse k with horizontal +x, right-handed vertical is -y.
    // The same WORLD +45-degree input therefore has negative local U.
    let backward = beam(
        633.0,
        Vec3::new(0.0, 0.0, -1.0),
        Stokes::try_new(1.0, 0.0, -1.0, 0.0).unwrap(),
    );
    let backward = backward
        .advance_corning7980_stress_optic(&stress, 0.03)
        .unwrap();
    near(backward.power().i, direct.power().i);
    near(backward.power().q, direct.power().q);
    near(backward.power().u, -direct.power().u);
    near(backward.power().v, -direct.power().v);
    assert_eq!(backward.ray().origin, Point3::new(0.0, 0.0, -0.03));
}

#[test]
fn physical_segment_composition_retains_stress_axis_order_and_cancellation() {
    let original = beam(
        633.0,
        Vec3::new(0.0, 0.0, 1.0),
        Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(),
    );
    let a = [400_000.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let b = [0.0, 0.0, 0.0, 300_000.0, 0.0, 0.0];
    let whole = original.advance_corning7980_stress_optic(&a, 0.05).unwrap();
    let split = original
        .advance_corning7980_stress_optic(&a, 0.02)
        .unwrap()
        .advance_corning7980_stress_optic(&a, 0.03)
        .unwrap();
    near_power(split.power(), whole.power());
    near(split.ray().origin.z, whole.ray().origin.z);
    let cancelled = original
        .advance_corning7980_stress_optic(&a, 0.03)
        .unwrap()
        .advance_corning7980_stress_optic(&a.map(|v| -v), 0.03)
        .unwrap();
    near_power(cancelled.power(), original.advance(0.06).unwrap().power());
    let ab = original
        .advance_corning7980_stress_optic(&a, 0.04)
        .unwrap()
        .advance_corning7980_stress_optic(&b, 0.04)
        .unwrap();
    let ba = original
        .advance_corning7980_stress_optic(&b, 0.04)
        .unwrap()
        .advance_corning7980_stress_optic(&a, 0.04)
        .unwrap();
    assert!((ab.power().q - ba.power().q).abs() > 0.05);
    near(ab.power().i, ba.power().i);
}

#[test]
fn hydrostatic_null_and_source_refusals_preserve_original_material_contract() {
    let original = beam(
        633.0,
        Vec3::new(1.0, 2.0, 3.0),
        Stokes::try_new(1.0, 0.2, -0.3, 0.4).unwrap(),
    );
    for stress in [[0.0; 6], [80.0e6, 80.0e6, 80.0e6, 0.0, 0.0, 0.0]] {
        assert_eq!(
            original
                .advance_corning7980_stress_optic(&stress, 0.2)
                .unwrap(),
            original.advance(0.2).unwrap()
        );
    }
    let other_line = beam(550.0, Vec3::new(0.0, 0.0, 1.0), Stokes::unpolarized(1.0));
    assert_eq!(
        other_line.advance_corning7980_stress_optic(&[0.0; 6], 0.0),
        Err(StressOpticBeamError::Source(
            SilicaStressOpticError::UnsupportedSourceWavelength
        ))
    );
    let other_material = PolarizedBeam::try_new(
        original.ray().origin,
        original.frame(),
        633.0,
        DielectricGlass::new(
            CauchyIor::try_constant(1.45).unwrap(),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        ),
        original.power(),
    )
    .unwrap();
    assert_eq!(
        other_material.advance_corning7980_stress_optic(&[0.0; 6], 0.1),
        Err(StressOpticBeamError::WrongMaterial)
    );
    let mut invalid = [0.0; 6];
    invalid[4] = f64::NAN;
    assert_eq!(
        original.advance_corning7980_stress_optic(&invalid, 0.0),
        Err(StressOpticBeamError::Source(
            SilicaStressOpticError::InvalidStress
        ))
    );
    assert_eq!(
        original.advance_corning7980_stress_optic(&[0.0; 6], -0.1),
        Err(StressOpticBeamError::Source(
            SilicaStressOpticError::InvalidDistance
        ))
    );
}

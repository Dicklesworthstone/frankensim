//! G1/G3 original-material controls for explicit coherent-factor statistics.
//! Authored optical fixtures; these are not measured material qualifications.
use super::*;
use crate::plane_wave::uniaxial::oriented::{
    OrientedUniaxialInterface,
    absorbing::{
        AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, AbsorbingCrystalCoherentIllumination,
        AbsorbingCrystalIllumination, AbsorbingCrystalPartiallyCoherentIllumination,
        CrystalIlluminationError,
    },
};
use fs_math::det;

fn box_limits() -> AbsorbingCrystalBoxLimits {
    AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: 128,
        maximum_workspace_bytes: 1 << 20,
    }
}

#[test]
fn partial_coherence_retains_phase_and_finite_footprint_fringe_contrast() {
    let outside = incident(1.0, MU);
    let (o, e) = axes([1.5; 2], [1.9; 2], [0.08, 0.2], Some(MU), [400.0, 800.0], 1);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside), law, [0.0, 0.0, 1.0],
    ).unwrap();
    let omega = frequency(550.0);
    let p = 0.45;
    let single = interface.prepare_absorbing_field_at_tangential_index(
        omega, p, [C64::ZERO, C64::ONE],
    ).unwrap();
    let kx = single.boundary_response().tangential_wavenumber_rad_m;
    let wavelength = 550e-9;
    let extent = [0.31 * wavelength, 0.4 * wavelength, 0.8 * wavelength];
    let z = 0.1 * wavelength;
    let single_heat = single.slab_absorption(
        extent[0] * extent[1], z, z + extent[2],
    ).unwrap().absorbed_power_w();
    let delta = kx * extent[0];
    let sinc = det::sin(delta) / delta;
    let phase = 0.7;

    // F0 = A + gamma*B, F1 = sqrt(1-|gamma|²)*B gives the positive
    // covariance [[1, conj(gamma)], [gamma, 1]]. Opposite signed kx
    // TE waves have identical transmission/normal decay at this optic axis.
    for coherence in [0.0_f64, 0.4, 1.0] {
        let gamma = C64::new(coherence * det::cos(phase), coherence * det::sin(phase));
        let mut fringe = vec![interface.prepare_absorbing_field_at_tangential_index(
            omega, p, [C64::ZERO, C64::ONE],
        ).unwrap()];
        if coherence > 0.0 {
            fringe.push(interface.prepare_absorbing_field_at_tangential_index(
                omega, -p, [C64::ZERO, gamma],
            ).unwrap());
        }
        let mut remainder = Vec::new();
        if coherence < 1.0 {
            remainder.push(interface.prepare_absorbing_field_at_tangential_index(
                omega, -p,
                [C64::ZERO, C64::new(det::sqrt(1.0 - coherence * coherence), 0.0)],
            ).unwrap());
        }
        let factors = [
            AbsorbingCrystalCoherentIllumination::try_new(&interface, &fringe, 2, || false).unwrap(),
            AbsorbingCrystalCoherentIllumination::try_new(&interface, &remainder, 1, || false).unwrap(),
        ];
        let illumination = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
            &interface, &factors, 2, || false,
        ).unwrap();
        let mut powers = [0.0; 2];
        for (index, sign) in [1.0, -1.0].into_iter().enumerate() {
            let shift = if index == 0 { 0.0 } else { core::f64::consts::PI };
            let center = (shift - phase) / (2.0 * kx);
            let region = AbsorbingCrystalBox::try_new(
                [center - extent[0] / 2.0, -0.2 * wavelength, z], extent,
            ).unwrap();
            let result = illumination.box_absorption(region, box_limits(), || false).unwrap();
            let expected = 2.0 * single_heat * (1.0 + sign * coherence * sinc);
            near(result.absorbed_power_w(), expected, 3e-11 * expected);
            let scale = result.outward_face_power_w().iter().map(|v| v.abs()).sum::<f64>();
            near(result.residual_w(), 0.0, 5e-11 * scale);
            powers[index] = result.absorbed_power_w();
        }
        near((powers[0] - powers[1]) / (powers[0] + powers[1]), coherence * sinc, 3e-11);
    }
}

#[test]
fn complex_unitary_factor_mixing_preserves_original_work_and_all_six_faces() {
    let outside = incident(1.0, MU);
    let (o, e) = axes([1.45, 1.65], [1.8, 2.1], [0.08, 0.2], Some(MU), [400.0, 800.0], 2);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside), law, [0.4, -0.2, 0.9],
    ).unwrap();
    let electric = [
        [C64::new(0.7, 0.1), C64::new(0.2, -0.3)],
        [C64::new(0.1, 0.6), C64::new(0.8, 0.1)],
    ];
    let prepare = |beam: usize, weight: C64| {
        interface.prepare_absorbing_field_at_tangential_index_in_plane(
            frequency(610.0), [0.25, 0.7][beam], [[1.0, 0.0], [0.6, 0.8]][beam],
            electric[beam].map(|value| value * weight),
        ).unwrap()
    };
    let original = [prepare(0, C64::ONE), prepare(1, C64::ONE)];
    let factors = [
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &original[..1], 1, || false).unwrap(),
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &original[1..], 1, || false).unwrap(),
    ];
    // U = [[0.6, 0.8i], [0.8i, 0.6]] is unitary. The two mixed
    // factors each contain spatial interference; those cross terms cancel
    // only AFTER their quadratic powers are formed and summed.
    let mixed = [
        [prepare(0, C64::new(0.6, 0.0)), prepare(1, C64::new(0.0, 0.8))],
        [prepare(0, C64::new(0.0, 0.8)), prepare(1, C64::new(0.6, 0.0))],
    ];
    let mixed_factors = [
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &mixed[0], 2, || false).unwrap(),
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &mixed[1], 2, || false).unwrap(),
    ];
    let original_bank = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &factors, 2, || false,
    ).unwrap();
    let mixed_bank = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &mixed_factors, 2, || false,
    ).unwrap();
    let region = AbsorbingCrystalBox::try_new(
        [0.12e-6, -0.23e-6, 0.15e-6], [0.37e-6, 0.41e-6, 0.93e-6],
    ).unwrap();
    let reference = original_bank.box_absorption(region, box_limits(), || false).unwrap();
    let rotated = mixed_bank.box_absorption(region, box_limits(), || false).unwrap();
    let independent = AbsorbingCrystalIllumination::try_new(&interface, &original, 2, || false).unwrap();
    let extent = region.extent_m();
    let z = region.origin_m()[2];
    let slab = independent.slab_absorption(extent[0] * extent[1], z, z + extent[2], || false).unwrap();
    near(reference.absorbed_power_w(), slab.absorbed_power_w(), 3e-11 * slab.absorbed_power_w());
    near(rotated.absorbed_power_w(), reference.absorbed_power_w(), 4e-11 * reference.absorbed_power_w());
    let faces = reference.outward_face_power_w();
    let scale = faces.iter().map(|v| v.abs()).sum::<f64>();
    for (actual, expected) in rotated.outward_face_power_w().into_iter().zip(faces) {
        near(actual, expected, 4e-11 * scale);
    }
    assert!(faces[0].abs() + faces[2].abs() > 1e-6 * scale);
    near(rotated.residual_w(), 0.0, 5e-11 * scale);
    assert_eq!(original_bank.mode_pair_count(), 8);
    assert_eq!(mixed_bank.mode_pair_count(), 32);
}

#[test]
fn independent_spectral_factors_keep_original_frequency_and_material_loss() {
    let outside = incident(1.0, MU);
    let (o, e) = axes([1.25, 1.85], [1.6, 2.4], [0.1, 0.3], Some(MU), [400.0, 800.0], 4);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside), law, [0.2, 0.3, 0.9],
    ).unwrap();
    let fields = [
        interface.prepare_absorbing_field_at_tangential_index(
            frequency(450.0), 0.2, [C64::ZERO, C64::ONE],
        ).unwrap(),
        interface.prepare_absorbing_field_at_tangential_index(
            frequency(700.0), -0.65, [C64::new(0.3, 0.1), C64::new(0.6, -0.2)],
        ).unwrap(),
    ];
    assert!(AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields, 2, || false).is_err());
    let factors = [
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields[..1], 1, || false).unwrap(),
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields[1..], 1, || false).unwrap(),
    ];
    let illumination = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &factors, 2, || false,
    ).unwrap();
    assert_eq!(illumination.factors()[0].omega_rad_s(), Some(frequency(450.0)));
    assert_eq!(illumination.factors()[1].omega_rad_s(), Some(frequency(700.0)));
    let region = AbsorbingCrystalBox::try_new(
        [0.1e-6, -0.2e-6, 0.15e-6], [0.3e-6, 0.4e-6, 0.7e-6],
    ).unwrap();
    let result = illumination.box_absorption(region, box_limits(), || false).unwrap();
    let independent = AbsorbingCrystalIllumination::try_new(&interface, &fields, 2, || false).unwrap();
    let slab = independent.slab_absorption(0.3e-6 * 0.4e-6, 0.15e-6, 0.85e-6, || false).unwrap();
    near(result.absorbed_power_w(), slab.absorbed_power_w(), 3e-11 * slab.absorbed_power_w());
    let mut expected_faces = [0.0; 6];
    for factor in &factors {
        let value = factor.box_absorption(region, box_limits(), || false).unwrap();
        for (sum, face) in expected_faces.iter_mut().zip(value.outward_face_power_w()) {
            *sum += face;
        }
    }
    let scale = expected_faces.iter().map(|v| v.abs()).sum::<f64>();
    for (actual, expected) in result.outward_face_power_w().into_iter().zip(expected_faces) {
        near(actual, expected, 3e-12 * scale);
    }
    assert!(core::ptr::eq(result.source(), &interface));
    assert_eq!(result.source().crystal_response().material().identity(), law.material().identity());
    near(result.residual_w(), 0.0, 5e-11 * scale);
}

#[test]
fn independent_opposite_factors_remain_bright_with_peak_storage_and_atomic_cancellation() {
    let outside = incident(1.0, MU);
    let (o, e) = axes([1.5; 2], [1.9; 2], [0.08, 0.2], Some(MU), [400.0, 800.0], 3);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let interface = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside), law, [0.0, 0.0, 1.0],
    ).unwrap();
    let other = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside), law, [0.0, 0.0, 1.0],
    ).unwrap();
    let fields = [1.0, -1.0].map(|amplitude| {
        interface.prepare_absorbing_field_at_tangential_index(
            frequency(550.0), 0.3, [C64::ZERO, C64::new(amplitude, 0.0)],
        ).unwrap()
    });
    let factors = [
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields[..1], 1, || false).unwrap(),
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields[1..], 1, || false).unwrap(),
    ];
    assert!(AbsorbingCrystalPartiallyCoherentIllumination::try_new(&other, &factors, 2, || false).is_err());
    assert!(AbsorbingCrystalPartiallyCoherentIllumination::try_new(&interface, &factors, 1, || false).is_err());
    assert!(matches!(
        AbsorbingCrystalPartiallyCoherentIllumination::try_new(&interface, &factors, 2, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    let illumination = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &factors, 2, || false,
    ).unwrap();
    let one = AbsorbingCrystalPartiallyCoherentIllumination::try_new(
        &interface, &factors[..1], 1, || false,
    ).unwrap();
    let coherent = AbsorbingCrystalCoherentIllumination::try_new(&interface, &fields, 2, || false).unwrap();
    let region = AbsorbingCrystalBox::try_new(
        [0.1e-6, -0.2e-6, 0.15e-6], [0.3e-6, 0.4e-6, 0.7e-6],
    ).unwrap();
    let result = illumination.box_absorption(region, box_limits(), || false).unwrap();
    let single = one.box_absorption(region, box_limits(), || false).unwrap();
    assert!(result.absorbed_power_w() > 0.0);
    assert_eq!(result.absorbed_power_w(), 2.0 * single.absorbed_power_w());
    assert_eq!(coherent.box_absorption(region, box_limits(), || false).unwrap().absorbed_power_w(), 0.0);
    assert_eq!(result.workspace_bytes(), single.workspace_bytes());
    assert_eq!(illumination.factor_count(), 2);
    assert_eq!(illumination.mode_pair_count(), 8);
    assert_eq!(result.mode_pair_count(), 8);
    assert!(core::ptr::eq(result.illumination(), &illumination));
    assert_eq!(result.region(), region);
    assert_eq!(result.volume_m3(), region.volume_m3());
    let exact = AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: result.mode_pair_count(),
        maximum_workspace_bytes: result.workspace_bytes(),
    };
    let repeated = illumination.box_absorption(region, exact, || false).unwrap();
    assert_eq!(result.absorbed_power_w(), repeated.absorbed_power_w());
    assert_eq!(result.outward_face_power_w(), repeated.outward_face_power_w());
    for insufficient in [
        AbsorbingCrystalBoxLimits { maximum_mode_pairs: exact.maximum_mode_pairs - 1, ..exact },
        AbsorbingCrystalBoxLimits { maximum_workspace_bytes: exact.maximum_workspace_bytes - 1, ..exact },
    ] {
        let mut polls = 0;
        assert!(illumination.box_absorption(region, insufficient, || { polls += 1; false }).is_err());
        assert_eq!(polls, 1); // refuse before any factor's modal allocation/integration
    }
    let mut total_polls = 0;
    illumination.box_absorption(region, exact, || { total_polls += 1; false }).unwrap();
    for cutoff in [1, 2, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            illumination.box_absorption(region, exact, || { polls += 1; polls >= cutoff }),
            Err(CrystalIlluminationError::Cancelled)
        ));
    }
    assert_eq!(illumination.box_absorption(region, exact, || false).unwrap().absorbed_power_w(), result.absorbed_power_w());
    let empty = AbsorbingCrystalPartiallyCoherentIllumination::try_new(&interface, &[], 0, || false).unwrap();
    let zeros = [
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &[], 0, || false).unwrap(),
        AbsorbingCrystalCoherentIllumination::try_new(&interface, &[], 0, || false).unwrap(),
    ];
    assert!(AbsorbingCrystalPartiallyCoherentIllumination::try_new(&interface, &zeros, 1, || false).is_err());
    let dark_factors = AbsorbingCrystalPartiallyCoherentIllumination::try_new(&interface, &zeros, 2, || false).unwrap();
    for bank in [&empty, &dark_factors] {
        let dark = bank.box_absorption(region, AbsorbingCrystalBoxLimits {
            maximum_mode_pairs: 0, ..box_limits()
        }, || false).unwrap();
        assert_eq!(dark.absorbed_power_w(), 0.0);
        assert_eq!(dark.outward_face_power_w(), [0.0; 6]);
        assert_eq!(dark.mode_pair_count(), 0);
    }
}

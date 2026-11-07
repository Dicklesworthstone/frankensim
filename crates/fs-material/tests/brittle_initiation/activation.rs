//! Endpoint observations must retain actual stress, objectivity and strict
//! affine-path admission without manufacturing an unobserved loading path.

use super::*;

#[test]
fn changed_strength_activates_the_supplied_strain_without_relaxing_path_admission() {
    let material = elastic();
    let held = strain([0.0, 0.75 * FT, 0.0, 0.0, 0.0, 0.0]);
    let original = TensileInitiator::new(&material, FT, resolution()).unwrap();
    assert!(original.activation_at(&held, None).unwrap().is_none());
    // Only the criterion changes. One bisection would not locate an invented
    // zero-to-held crossing; endpoint observation needs no bisection at all.
    let weakened = TensileInitiator::new(
        &material,
        0.5 * FT,
        TensileInitiationSettings {
            maximum_bisections: 1,
            ..resolution()
        },
    )
    .unwrap();
    let observed = weakened.activation_at(&held, None).unwrap().unwrap();
    assert_eq!(observed.strain(), &held);
    assert_eq!(observed.stress_pa(), material.stress(&held, &()));
    assert!(core::ptr::eq(observed.elastic(), &material));
    assert_eq!(observed.tensile_strength_pa(), 0.5 * FT);
    close(observed.principal_stresses_pa()[2], 0.75 * FT, 1e-8);
    same_plane(observed.normal(), [0.0, 1.0, 0.0], 1e-14);
    assert_eq!(observed.principal_multiplicity(), 1);
    assert!(matches!(
        weakened.between(&held, &held, None),
        Err(TensileInitiationError::StartNotElastic { .. })
    ));

    // Exactly representable coefficients isolate equality at the threshold.
    let binary = IsotropicElastic::new(1024.0, 0.0, 0.05).unwrap();
    let detector = TensileInitiator::new(&binary, 1.0, resolution()).unwrap();
    for (stress, active) in [(0.999, false), (1.0, true), (1.001, true)] {
        let observed = detector
            .activation_at(&[stress / 1024.0, 0.0, 0.0, 0.0, 0.0, 0.0], None)
            .unwrap();
        assert_eq!(observed.is_some(), active);
    }
    for stress in [
        [-FT, -FT, -FT, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.2 * FT, 0.0, 0.0],
    ] {
        assert!(
            original
                .activation_at(&strain(stress), None)
                .unwrap()
                .is_none()
        );
    }
}

#[test]
fn endpoint_planes_rotate_and_repeated_principal_tension_requires_a_physical_hint() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    let loaded = strain([0.6 * FT, 0.8 * FT, -0.4 * FT, 0.8 * FT, 0.0, 0.0]);
    let base = detector.activation_at(&loaded, None).unwrap().unwrap();
    let q = rotation([0.0, 1.0, 0.0], 0.7);
    for factor in [1e-120, 1.0, 1e100] {
        let scaled = IsotropicElastic::new(E * factor, NU, 0.05).unwrap();
        let scaled_detector = TensileInitiator::new(&scaled, FT * factor, resolution()).unwrap();
        let observed = scaled_detector
            .activation_at(&rotate(&loaded, &q), None)
            .unwrap()
            .unwrap();
        same_plane(observed.normal(), vector_rotate(base.normal(), &q), 3e-12);
        close(
            observed.principal_stresses_pa()[2] / factor,
            base.principal_stresses_pa()[2],
            3e-5,
        );
        assert!(
            observed.eigen_residual_pa() <= FT * factor * resolution().relative_stress_tolerance
        );
    }

    let repeated = strain([2.0 * FT, 2.0 * FT, -FT, 0.0, 0.0, 0.0]);
    assert!(matches!(
        detector.activation_at(&repeated, None),
        Err(TensileInitiationError::AmbiguousPlane { multiplicity: 2 })
    ));
    assert!(matches!(
        detector.activation_at(&repeated, Some([0.0, 0.0, 1.0])),
        Err(TensileInitiationError::DirectionOutsidePrincipalSpace)
    ));
    let hint = [1.0, 2.0, 3.0];
    let observed = detector
        .activation_at(&repeated, Some(hint))
        .unwrap()
        .unwrap();
    same_plane(
        observed.normal(),
        [1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt(), 0.0],
        1e-12,
    );
    let rotated = detector
        .activation_at(&rotate(&repeated, &q), Some(vector_rotate(hint, &q)))
        .unwrap()
        .unwrap();
    same_plane(
        rotated.normal(),
        vector_rotate(observed.normal(), &q),
        1e-11,
    );
    let hydrostatic = strain([2.0 * FT, 2.0 * FT, 2.0 * FT, 0.0, 0.0, 0.0]);
    let hydrostatic = detector
        .activation_at(&hydrostatic, Some(hint))
        .unwrap()
        .unwrap();
    assert_eq!(hydrostatic.principal_multiplicity(), 3);
    same_plane(
        hydrostatic.normal(),
        hint.map(|v| v / 14.0_f64.sqrt()),
        1e-12,
    );
}

#[test]
fn endpoint_domain_and_spectral_cancellation_refuse_without_a_partial_observation() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    let loaded = strain([2.0 * FT, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(detector.activation_at(&[f64::NAN; 6], None).is_err());
    assert!(matches!(
        detector.activation_at(&[0.02; 6], None),
        Err(TensileInitiationError::Material(
            BrittleError::StrainLimit { .. }
        ))
    ));
    assert!(detector.activation_at(&loaded, Some([0.0; 3])).is_err());
    for cancel_at in [1, 3] {
        let mut polls = 0;
        assert!(matches!(
            detector.activation_at_controlled(&loaded, None, &mut || {
                polls += 1;
                polls == cancel_at
            }),
            Err(TensileInitiationError::Cancelled)
        ));
    }
    assert_eq!(loaded, strain([2.0 * FT, 0.0, 0.0, 0.0, 0.0, 0.0]));
    assert!(detector.activation_at(&loaded, None).unwrap().is_some());
}

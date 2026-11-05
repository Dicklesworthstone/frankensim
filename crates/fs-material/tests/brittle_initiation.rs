//! G0/G1/G3/G4: physical tensile planes, affine onset and conservative handoff.
use fs_material::brittle::initiation::{
    TensileInitiationError, TensileInitiationSettings, TensileInitiator,
};
use fs_material::brittle::{
    BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings, BrittleError,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::tensor::{contract, rotate, rotation};
use fs_material::{SmallStrainLaw, Voigt};

const E: f64 = 3e9;
const NU: f64 = 0.25;
const FT: f64 = 3e6;

fn resolution() -> TensileInitiationSettings {
    TensileInitiationSettings {
        maximum_strain: 0.05,
        relative_stress_tolerance: 1e-10,
        fraction_tolerance: 1e-12,
        maximum_bisections: 64,
    }
}
fn elastic() -> IsotropicElastic {
    IsotropicElastic::new(E, NU, 0.05).unwrap()
}
fn strain(stress: Voigt) -> Voigt {
    [
        (stress[0] - NU * (stress[1] + stress[2])) / E,
        (stress[1] - NU * (stress[0] + stress[2])) / E,
        (stress[2] - NU * (stress[0] + stress[1])) / E,
        (1.0 + NU) * stress[3] / E,
        (1.0 + NU) * stress[4] / E,
        (1.0 + NU) * stress[5] / E,
    ]
}
fn band(elastic: IsotropicElastic, normal: [f64; 3], shear_strength: f64) -> BrittleCrackBand {
    BrittleCrackBand::new(
        elastic,
        BrittleCrackBandParameters {
            tensile_strength_pa: FT,
            shear_strength_pa: shear_strength,
            mode_i_fracture_energy_j_m2: 120.0,
            mode_ii_fracture_energy_j_m2: 360.0,
            bk_exponent: 1.6,
        },
        BrittleCrackBandSettings {
            band_width_m: 0.004,
            normal,
            maximum_strain: 0.05,
            energy_tolerance_j_m3: 1e-6,
        },
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn same_plane(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for i in 0..3 {
        for j in 0..3 {
            close(a[i] * a[j], b[i] * b[j], tolerance);
        }
    }
}
fn vector_rotate(a: [f64; 3], q: &[[f64; 3]; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|j| q[i][j] * a[j]).sum())
}

#[test]
fn g1_affine_uniaxial_onset_seeds_original_elastic_storage_and_native_fracture() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    let start = strain([0.2 * FT, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let end = strain([4.0 * FT, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let event = detector.between(&start, &end, None).unwrap().unwrap();
    assert!(core::ptr::eq(event.elastic(), &material));
    let [before, after] = event.fraction_bracket();
    let expected = 0.8 / 3.8;
    assert!(before <= expected && expected <= after);
    assert!(after - before <= resolution().fraction_tolerance);
    close(
        event.principal_stresses_pa()[2],
        FT,
        2.0 * FT * resolution().relative_stress_tolerance,
    );
    same_plane(event.normal(), [1.0, 0.0, 0.0], 1e-13);
    assert_eq!(event.principal_multiplicity(), 1);
    let receiving = band(material.clone(), event.normal(), 2.0 * FT);
    let handoff = event.seed(&receiving).unwrap();
    assert_eq!(handoff.state.strain(), &start);
    assert_eq!(handoff.state.damage(), 0.0);
    assert_eq!(handoff.state.dissipated_j_m3(), 0.0);
    assert!(
        handoff
            .state
            .cohesive_point()
            .mixed_mode_history()
            .is_none()
    );
    let original_storage = 0.5 * contract(&start, &material.stress(&start, &()));
    close(handoff.elastic_storage_j_m3, original_storage, 1e-12);
    close(handoff.state.stored_energy_j_m3(), original_storage, 1e-10);
    let after = receiving.trial(&end, &handoff.state).unwrap();
    assert!(after.state.damage() > 0.0 && after.state.damage() < 1.0);
    close(
        after.energy.work,
        after.state.stored_energy_j_m3() - original_storage + after.state.dissipated_j_m3(),
        1e-7,
    );
    close(
        after
            .state
            .cohesive_point()
            .mixed_mode_history()
            .unwrap()
            .energy_shear_fraction,
        0.0,
        1e-24,
    );
    assert_eq!(handoff.state.damage(), 0.0);
}

#[test]
fn g0_compression_poisson_expansion_and_subcritical_shear_do_not_nucleate_tension() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    for stress in [
        [-20.0 * FT, -20.0 * FT, -20.0 * FT, 0.0, 0.0, 0.0],
        [-10.0 * FT, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.8 * FT, 0.0, 0.0],
        [0.9 * FT, -5.0 * FT, -8.0 * FT, 0.0, 0.0, 0.0],
    ] {
        for angle in [0.0, 0.3, 1.1, 2.7] {
            let q = rotation([0.0, 0.0, 1.0], angle);
            let end = rotate(&strain(stress), &q);
            assert!(detector.between(&[0.0; 6], &end, None).unwrap().is_none());
            assert!(detector.between(&end, &[0.0; 6], None).unwrap().is_none());
        }
    }
}

#[test]
fn g1_nonproportional_affine_path_matches_independent_quadratic_crossing() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    let start = strain([0.2 * FT, -0.3 * FT, -0.4 * FT, 0.15 * FT, 0.0, 0.0]);
    let end = strain([0.6 * FT, 0.8 * FT, -0.4 * FT, 0.8 * FT, 0.0, 0.0]);
    // det(I - stress_xy/FT) = 1.0175 - 1.595*t + .0175*t^2.
    let root = 2.0 * 1.0175 / (1.595 + (1.595_f64.powi(2) - 4.0 * 0.0175 * 1.0175).sqrt());
    let event = detector.between(&start, &end, None).unwrap().unwrap();
    let [low, high] = event.fraction_bracket();
    assert!(low <= root + 2e-15 && root - 2e-15 <= high);
    let a = 0.2 + 0.4 * root;
    let b = 0.15 + 0.65 * root;
    let length = (b * b + (1.0 - a).powi(2)).sqrt();
    same_plane(event.normal(), [b / length, (1.0 - a) / length, 0.0], 2e-12);
    let receiving = band(material.clone(), event.normal(), 2.0 * FT);
    let handoff = event.seed(&receiving).unwrap();
    let after = receiving.trial(&end, &handoff.state).unwrap();
    assert!(after.state.damage() > 0.0);
    let history = after.state.cohesive_point().mixed_mode_history().unwrap();
    assert!(history.energy_shear_fraction < 1e-18, "{:?}", history);
}

#[test]
fn g3_plane_and_onset_rotate_and_rescale_with_the_original_stress() {
    let start = strain([0.2 * FT, -0.3 * FT, -0.4 * FT, 0.15 * FT, 0.0, 0.0]);
    let end = strain([0.6 * FT, 0.8 * FT, -0.4 * FT, 0.8 * FT, 0.0, 0.0]);
    let material = elastic();
    let base = TensileInitiator::new(&material, FT, resolution())
        .unwrap()
        .between(&start, &end, None)
        .unwrap()
        .unwrap();
    for factor in [1e-120, 1e-6, 1.0, 1e100] {
        let scaled = IsotropicElastic::new(E * factor, NU, 0.05).unwrap();
        let detector = TensileInitiator::new(&scaled, FT * factor, resolution()).unwrap();
        for q in [
            rotation([1.0, 0.0, 0.0], 0.6),
            rotation([0.0, 1.0, 0.0], 1.3),
            rotation([0.0, 0.0, 1.0], -0.7),
        ] {
            let event = detector
                .between(&rotate(&start, &q), &rotate(&end, &q), None)
                .unwrap()
                .unwrap();
            close(
                event.fraction_bracket()[1],
                base.fraction_bracket()[1],
                2e-12,
            );
            same_plane(event.normal(), vector_rotate(base.normal(), &q), 3e-12);
            close(
                event.principal_stresses_pa()[2] / factor,
                base.principal_stresses_pa()[2],
                3e-5,
            );
            assert!(
                event.eigen_residual_pa() <= FT * factor * resolution().relative_stress_tolerance
            );
        }
    }
}

#[test]
fn g3_repeated_principal_tension_requires_a_physical_covariant_direction() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    let end = strain([2.0 * FT, 2.0 * FT, -FT, 0.0, 0.0, 0.0]);
    assert!(matches!(
        detector.between(&[0.0; 6], &end, None),
        Err(TensileInitiationError::AmbiguousPlane { multiplicity: 2 })
    ));
    assert!(matches!(
        detector.between(&[0.0; 6], &end, Some([0.0, 0.0, 1.0])),
        Err(TensileInitiationError::DirectionOutsidePrincipalSpace)
    ));
    let hint = [1.0, 2.0, 3.0];
    let base = detector
        .between(&[0.0; 6], &end, Some(hint))
        .unwrap()
        .unwrap();
    same_plane(
        base.normal(),
        [1.0 / 5.0_f64.sqrt(), 2.0 / 5.0_f64.sqrt(), 0.0],
        1e-12,
    );
    let q = rotation([0.0, 1.0, 0.0], 0.7);
    let rotated = detector
        .between(&[0.0; 6], &rotate(&end, &q), Some(vector_rotate(hint, &q)))
        .unwrap()
        .unwrap();
    same_plane(rotated.normal(), vector_rotate(base.normal(), &q), 1e-11);
    let hydro = strain([2.0 * FT, 2.0 * FT, 2.0 * FT, 0.0, 0.0, 0.0]);
    let event = detector
        .between(&[0.0; 6], &hydro, Some(hint))
        .unwrap()
        .unwrap();
    assert_eq!(event.principal_multiplicity(), 3);
    same_plane(event.normal(), hint.map(|x| x / 14.0_f64.sqrt()), 1e-12);
}

#[test]
fn g0_handoff_refuses_changed_material_plane_and_prior_mixed_mode_onset() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    let end = strain([2.0 * FT, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let event = detector.between(&[0.0; 6], &end, None).unwrap().unwrap();
    assert!(matches!(
        event.seed(&band(material.clone(), [0.0, 1.0, 0.0], 2.0 * FT)),
        Err(TensileInitiationError::PlaneMismatch)
    ));
    assert!(matches!(
        event.seed(&band(
            IsotropicElastic::new(1.1 * E, NU, 0.05).unwrap(),
            event.normal(),
            2.0 * FT
        )),
        Err(TensileInitiationError::MaterialMismatch)
    ));
    // The elastic start lies BELOW Rankine, but its shear on the eventual x
    // plane exceeds this receiving law's separately declared shear strength.
    let start = strain([0.2 * FT, -2.0 * FT, -2.0 * FT, 0.9 * FT, 0.0, 0.0]);
    let end = strain([1.8 * FT, -2.0 * FT, -2.0 * FT, -0.9 * FT, 0.0, 0.0]);
    let event = detector.between(&start, &end, None).unwrap().unwrap();
    same_plane(event.normal(), [1.0, 0.0, 0.0], 2e-12);
    assert!(matches!(
        event.seed(&band(material.clone(), event.normal(), 0.1 * FT)),
        Err(TensileInitiationError::PriorMixedModeInitiation)
    ));
}

#[test]
fn g0_g4_domain_resolution_and_cancellation_refusals_leave_inputs_unchanged() {
    let material = elastic();
    let detector = TensileInitiator::new(&material, FT, resolution()).unwrap();
    let start = [0.0; 6];
    let end = strain([2.0 * FT, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(matches!(
        detector.between(&end, &start, None),
        Err(TensileInitiationError::StartNotElastic { .. })
    ));
    assert!(matches!(
        detector.between(&start, &[0.02; 6], None),
        Err(TensileInitiationError::Material(
            BrittleError::StrainLimit { .. }
        ))
    ));
    assert!(detector.between(&start, &[f64::NAN; 6], None).is_err());
    assert!(detector.between(&start, &end, Some([0.0; 3])).is_err());
    for settings in [
        TensileInitiationSettings {
            maximum_strain: 0.051,
            ..resolution()
        },
        TensileInitiationSettings {
            relative_stress_tolerance: 0.0,
            ..resolution()
        },
        TensileInitiationSettings {
            fraction_tolerance: 0.0,
            ..resolution()
        },
        TensileInitiationSettings {
            maximum_bisections: 0,
            ..resolution()
        },
        TensileInitiationSettings {
            maximum_bisections: 129,
            ..resolution()
        },
    ] {
        assert!(TensileInitiator::new(&material, FT, settings).is_err());
    }
    let budget = TensileInitiator::new(
        &material,
        FT,
        TensileInitiationSettings {
            maximum_bisections: 1,
            ..resolution()
        },
    )
    .unwrap();
    assert!(matches!(
        budget.between(&start, &end, None),
        Err(TensileInitiationError::OnsetBudget)
    ));
    let mut polls = 0;
    assert!(matches!(
        detector.between_controlled(&start, &end, None, &mut || {
            polls += 1;
            polls == 4
        }),
        Err(TensileInitiationError::Cancelled)
    ));
    assert_eq!(start, [0.0; 6]);
    assert_eq!(end, strain([2.0 * FT, 0.0, 0.0, 0.0, 0.0, 0.0]));
    assert!(detector.between(&start, &end, None).unwrap().is_some());
}

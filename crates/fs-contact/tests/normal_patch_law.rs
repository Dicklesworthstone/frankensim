use fs_contact::normal_patch::*;
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn request(law: NormalPatchLaw) -> NormalPatchRequest {
    NormalPatchRequest {
        identity: NormalPatchIdentity {
            model_id: "test/hertz-normal-patch-v1".into(),
            source_id: "test/synthetic-card-v1".into(),
            state_id: "test/loading".into(),
        },
        interface: InterfaceSystemRef::new(
            "test/a->b",
            "test/history",
            "test/source",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap(),
        law,
        indentation_m: 1.0e-5,
        indentation_rate_m_per_s: 2.0e-3,
        step_s: 1.0e-4,
        line_load_n_per_m: 100.0,
        applicability: ApplicabilityInput {
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            yield_strength_pa: 1.0e9,
            characteristic_rate_m_per_s: 1.0,
            temperature_k: 293.15,
            adhesion_energy_j_per_m2: 0.0,
        },
        limits: ApplicabilityLimits {
            max_patch_to_radius: 0.1,
            max_strain: 0.01,
            max_patch_to_depth: 0.1,
            max_patch_to_layer: 0.1,
            max_pressure_to_yield: 0.1,
            max_rate_ratio: 0.1,
            min_temperature_k: 250.0,
            max_temperature_k: 350.0,
        },
        uncertainty: InputUncertainty {
            radius_relative: 0.01,
            modulus_relative: 0.02,
            load_relative: 0.03,
        },
    }
}

#[test]
fn np_001_sphere_reconstructs_force_energy_and_derivative() {
    let req = request(NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    let r = req.evaluate().unwrap();
    let a = (0.02_f64 * 1.0e-5).sqrt();
    let force = 4.0 / 3.0 * 2.0e9 * 0.02_f64.sqrt() * (1.0e-5_f64).powf(1.5);
    assert!((r.normal_force_n_or_n_per_m - force).abs() <= force * 1e-12);
    assert!(
        (r.pressure.resultant_n_or_n_per_m
            - 2.0 / 3.0 * core::f64::consts::PI * a * a * r.pressure.peak_pressure_pa)
            .abs()
            <= force * 1e-12
    );
    assert!((r.pressure.second_moment_m2 - 0.4 * a * a).abs() <= 1e-18);
    assert!((r.reversible_energy_j_or_j_per_m - 0.4 * force * 1.0e-5).abs() <= 1e-18);
    let h = 1.0e-8;
    let mut plus = req.clone();
    plus.indentation_m += h;
    let mut minus = req;
    minus.indentation_m -= h;
    let numerical = (plus.evaluate().unwrap().normal_force_n_or_n_per_m
        - minus.evaluate().unwrap().normal_force_n_or_n_per_m)
        / (2.0 * h);
    assert!((numerical - r.tangent_n_per_m_or_pa).abs() / r.tangent_n_per_m_or_pa < 2.0e-6);
}

#[test]
fn np_002_cylinder_reconstructs_line_resultant_and_scale() {
    let req = request(NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    let r = req.evaluate().unwrap();
    let a = r.patch_radius_or_half_width_m;
    assert!((core::f64::consts::PI * a * r.pressure.peak_pressure_pa / 2.0 - 100.0).abs() < 1e-9);
    assert!((r.pressure.second_moment_m2 - a * a / 4.0).abs() < 1e-18);
    let h = 1.0e-3;
    let mut plus = req.clone();
    plus.line_load_n_per_m += h;
    let mut minus = req.clone();
    minus.line_load_n_per_m -= h;
    let numerical =
        (plus.evaluate().unwrap().approach_m - minus.evaluate().unwrap().approach_m) / (2.0 * h);
    assert!((1.0 / numerical - r.tangent_n_per_m_or_pa).abs() / r.tangent_n_per_m_or_pa < 1e-9);
    let mut scaled = req;
    scaled.law = NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.04,
        reduced_modulus_pa: 2.0e9,
    };
    scaled.line_load_n_per_m = 200.0;
    assert!((scaled.evaluate().unwrap().patch_radius_or_half_width_m / a - 2.0).abs() < 1e-12);
}

#[test]
fn np_003_hunt_crossley_loading_and_unloading_are_passive() {
    let law = NormalPatchLaw::HuntCrossleySphere {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
        dissipation_s_per_m: 10.0,
    };
    let loading = request(law).evaluate().unwrap();
    let mut unloading_request = request(law);
    unloading_request.indentation_rate_m_per_s = -1.0e-3;
    let unloading = unloading_request.evaluate().unwrap();
    assert!(loading.irreversible_work_j > 0.0 && loading.dissipated_power_w > 0.0);
    assert!(unloading.irreversible_work_j > 0.0 && unloading.dissipated_power_w > 0.0);
    assert!(unloading.normal_force_n_or_n_per_m < loading.normal_force_n_or_n_per_m);
}

#[test]
fn np_004_zero_and_hostile_applicability_refuse_without_receipt() {
    let mut zero = request(NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    zero.indentation_m = 0.0;
    let response = zero.evaluate().unwrap();
    assert_eq!(response.normal_force_n_or_n_per_m, 0.0);
    let mut adhesive = zero.clone();
    adhesive.applicability.adhesion_energy_j_per_m2 = 1.0;
    assert!(matches!(
        adhesive.evaluate(),
        Err(NormalPatchError::AdhesionUnsupported { .. })
    ));
    let mut layered = zero;
    layered.indentation_m = 1e-4;
    layered.applicability.layer_thickness_m = 1e-6;
    assert!(matches!(
        layered.evaluate(),
        Err(NormalPatchError::OutsideApplicability {
            ratio: "patch_to_layer",
            ..
        })
    ));
}

#[test]
fn np_005_identity_authority_and_barrier_are_deterministic_and_distinct() {
    let req = request(NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    let a = req.evaluate().unwrap();
    let b = req.evaluate().unwrap();
    assert_eq!(a.request_id, b.request_id);
    assert_eq!(a.receipt_id, b.receipt_id);
    assert_eq!(a.authority, InputAuthority::SyntheticFixture);
    let barrier = ConstraintBarrierReceipt {
        request_id: a.request_id,
        scheme_id: "ipc/barrier".into(),
    };
    assert_eq!(barrier.request_id, a.request_id);
}

#[test]
fn np_006_canonical_identity_covers_limits_uncertainty_and_malformed_extremes() {
    let req = request(NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    let baseline = req.evaluate().unwrap();
    let mut changed_limit = req.clone();
    changed_limit.limits.max_strain *= 2.0;
    assert_ne!(
        changed_limit.evaluate().unwrap().request_id,
        baseline.request_id
    );
    let mut changed_uncertainty = req.clone();
    changed_uncertainty.uncertainty.load_relative *= 2.0;
    assert_ne!(
        changed_uncertainty.evaluate().unwrap().request_id,
        baseline.request_id
    );
    let mut malformed = req.clone();
    malformed.limits.max_rate_ratio = f64::NAN;
    assert!(matches!(
        malformed.evaluate(),
        Err(NormalPatchError::InvalidInput {
            field: "max_rate_ratio"
        })
    ));
    let mut extreme = req;
    extreme.indentation_m = f64::MAX;
    assert!(matches!(
        extreme.evaluate(),
        Err(NormalPatchError::Overflow { .. })
    ));
}

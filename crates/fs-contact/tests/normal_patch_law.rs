use fs_contact::normal_patch::*;
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn request(law: NormalPatchLaw) -> NormalPatchRequest {
    let geometry = match law {
        NormalPatchLaw::HertzCylinderPlane { .. } => NormalPatchGeometry::CylinderPlane,
        NormalPatchLaw::HertzSpherePlane { .. } | NormalPatchLaw::HuntCrossleySphere { .. } => {
            NormalPatchGeometry::SpherePlane
        }
        NormalPatchLaw::HertzEllipticParaboloid { .. } => NormalPatchGeometry::EllipticParaboloid,
    };
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
        geometry,
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

fn point(receipt: NormalPatchReceipt) -> Option<PointNormalPatchReceipt> {
    match receipt {
        NormalPatchReceipt::Point(receipt) => Some(receipt),
        NormalPatchReceipt::Line(_) => None,
    }
}

fn line(receipt: NormalPatchReceipt) -> Option<LineNormalPatchReceipt> {
    match receipt {
        NormalPatchReceipt::Line(receipt) => Some(receipt),
        NormalPatchReceipt::Point(_) => None,
    }
}

#[test]
fn np_001_sphere_reconstructs_force_energy_and_derivative() {
    let req = request(NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    let r = point(req.evaluate().unwrap()).unwrap();
    let a = (0.02_f64 * 1.0e-5).sqrt();
    let force = 4.0 / 3.0 * 2.0e9 * 0.02_f64.sqrt() * (1.0e-5_f64).powf(1.5);
    assert!((r.normal_force_n - force).abs() <= force * 1e-12);
    assert!(
        (r.pressure.resultant_n
            - 2.0 / 3.0 * core::f64::consts::PI * a * a * r.pressure.peak_pressure_pa)
            .abs()
            <= force * 1e-12
    );
    assert!((r.pressure.second_moment_m2 - 0.4 * a * a).abs() <= 1e-18);
    assert!((r.reversible_energy_j - 0.4 * force * 1.0e-5).abs() <= 1e-18);
    let h = 1.0e-8;
    let mut plus = req.clone();
    plus.indentation_m += h;
    let mut minus = req;
    minus.indentation_m -= h;
    let numerical = (point(plus.evaluate().unwrap()).unwrap().normal_force_n
        - point(minus.evaluate().unwrap()).unwrap().normal_force_n)
        / (2.0 * h);
    assert!((numerical - r.tangent_n_per_m).abs() / r.tangent_n_per_m < 2.0e-6);
}

#[test]
fn np_002_cylinder_reconstructs_line_resultant_and_scale() {
    let req = request(NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    let r = line(req.evaluate().unwrap()).unwrap();
    let a = r.patch_half_width_m;
    assert!((core::f64::consts::PI * a * r.pressure.peak_pressure_pa / 2.0 - 100.0).abs() < 1e-9);
    assert!((r.pressure.second_moment_m2 - a * a / 4.0).abs() < 1e-18);
    let log = (4.0 * 0.02 / a).ln();
    let expected_approach = 100.0 * (log + 0.5) / (core::f64::consts::PI * 2.0e9);
    let expected_energy = 100.0_f64.powi(2) * (log + 0.75) / (2.0 * core::f64::consts::PI * 2.0e9);
    assert!((r.approach_m - expected_approach).abs() <= expected_approach * 1e-12);
    assert!((r.reversible_energy_j_per_m - expected_energy).abs() <= expected_energy * 1e-12);
    assert_eq!(r.units, LINE_SI_UNITS);
    let h = 1.0e-3;
    let mut plus = req.clone();
    plus.line_load_n_per_m += h;
    let mut minus = req.clone();
    minus.line_load_n_per_m -= h;
    let numerical = (line(plus.evaluate().unwrap()).unwrap().approach_m
        - line(minus.evaluate().unwrap()).unwrap().approach_m)
        / (2.0 * h);
    assert!((1.0 / numerical - r.tangent_pa).abs() / r.tangent_pa < 1e-9);
    let mut scaled = req;
    scaled.law = NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.04,
        reduced_modulus_pa: 2.0e9,
    };
    scaled.line_load_n_per_m = 200.0;
    assert!((line(scaled.evaluate().unwrap()).unwrap().patch_half_width_m / a - 2.0).abs() < 1e-12);
}

#[test]
fn np_003_hunt_crossley_loading_and_unloading_are_passive() {
    let law = NormalPatchLaw::HuntCrossleySphere {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
        dissipation_s_per_m: 10.0,
    };
    let loading = point(request(law).evaluate().unwrap()).unwrap();
    let mut unloading_request = request(law);
    unloading_request.indentation_rate_m_per_s = -1.0e-3;
    let unloading = point(unloading_request.evaluate().unwrap()).unwrap();
    assert!(loading.irreversible_work_j > 0.0 && loading.dissipated_power_w > 0.0);
    assert!(unloading.irreversible_work_j > 0.0 && unloading.dissipated_power_w > 0.0);
    assert!(unloading.normal_force_n < loading.normal_force_n);
}

#[test]
fn np_004_zero_and_hostile_applicability_refuse_without_receipt() {
    let mut zero = request(NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    zero.indentation_m = 0.0;
    let response = point(zero.evaluate().unwrap()).unwrap();
    assert_eq!(response.normal_force_n, 0.0);
    let mut zero_line = request(NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    zero_line.line_load_n_per_m = 0.0;
    assert_eq!(
        line(zero_line.evaluate().unwrap())
            .unwrap()
            .normal_line_load_n_per_m,
        0.0
    );
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
    let a = point(req.evaluate().unwrap()).unwrap();
    let b = point(req.evaluate().unwrap()).unwrap();
    assert_eq!(a.request_id, b.request_id);
    assert_eq!(a.receipt_id, b.receipt_id);
    assert_eq!(a.authority, InputAuthority::SyntheticFixture);
    assert_eq!(a.interface_system_id, "test/a->b");
    assert_eq!(a.history_id, "test/history");
    assert_eq!(a.input_source_id, "test/source");
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
        changed_limit.evaluate().unwrap().request_id(),
        baseline.request_id()
    );
    let mut changed_uncertainty = req.clone();
    changed_uncertainty.uncertainty.load_relative *= 2.0;
    assert_ne!(
        changed_uncertainty.evaluate().unwrap().request_id(),
        baseline.request_id()
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

#[test]
fn np_007_applicability_boundaries_accept_then_refuse() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let baseline = point(request(law).evaluate().unwrap()).unwrap();

    let mut temperature = request(law);
    temperature.applicability.temperature_k = temperature.limits.min_temperature_k;
    assert!(temperature.evaluate().is_ok());
    temperature.applicability.temperature_k -= 1.0e-9;
    assert!(matches!(
        temperature.evaluate(),
        Err(NormalPatchError::OutsideApplicability {
            ratio: "temperature_k",
            ..
        })
    ));

    let mut rate = request(law);
    rate.indentation_rate_m_per_s =
        rate.applicability.characteristic_rate_m_per_s * rate.limits.max_rate_ratio;
    assert!(rate.evaluate().is_ok());
    rate.indentation_rate_m_per_s *= 1.0 + 1.0e-9;
    assert!(matches!(
        rate.evaluate(),
        Err(NormalPatchError::OutsideApplicability {
            ratio: "rate_ratio",
            ..
        })
    ));

    let mut yield_boundary = request(law);
    yield_boundary.limits.max_pressure_to_yield = baseline.ratios.pressure_to_yield;
    assert!(yield_boundary.evaluate().is_ok());
    yield_boundary.limits.max_pressure_to_yield *= 1.0 - 1.0e-9;
    assert!(matches!(
        yield_boundary.evaluate(),
        Err(NormalPatchError::OutsideApplicability {
            ratio: "pressure_to_yield",
            ..
        })
    ));
}

#[test]
fn np_008_identity_mutations_and_unsupported_geometry_refuse() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let base = request(law);
    let base_id = base.evaluate().unwrap().request_id();
    let mut model = base.clone();
    model.identity.model_id.push_str("/mutated");
    assert_ne!(model.evaluate().unwrap().request_id(), base_id);
    let mut source = base.clone();
    source.identity.source_id.push_str("/mutated");
    assert_ne!(source.evaluate().unwrap().request_id(), base_id);
    let mut state = base.clone();
    state.identity.state_id.push_str("/mutated");
    assert_ne!(state.evaluate().unwrap().request_id(), base_id);
    let mut toroidal = base.clone();
    toroidal.geometry = NormalPatchGeometry::ToroidalOrHighlyElliptical;
    assert!(matches!(
        toroidal.evaluate(),
        Err(NormalPatchError::UnsupportedGeometry {
            geometry: NormalPatchGeometry::ToroidalOrHighlyElliptical,
        })
    ));
    let mut mismatch = base;
    mismatch.geometry = NormalPatchGeometry::CylinderPlane;
    assert!(matches!(
        mismatch.evaluate(),
        Err(NormalPatchError::GeometryLawMismatch {
            geometry: NormalPatchGeometry::CylinderPlane,
        })
    ));
}

#[test]
fn np_009_elliptic_hertz_reconstructs_axes_resultant_energy_and_tangent() {
    let law = NormalPatchLaw::HertzEllipticParaboloid {
        maximum_principal_curvature_m_inverse: 100.0,
        minimum_principal_curvature_m_inverse: 25.0,
        reduced_modulus_pa: 2.0e9,
    };
    let req = request(law);
    let receipt = point(req.evaluate().unwrap()).unwrap();
    let axes = receipt
        .elliptic_patch_axes
        .expect("elliptic law records axes");
    assert!(axes.semi_major_axis_m >= axes.semi_minor_axis_m);
    assert!(axes.aspect_ratio > 1.0);
    assert!((axes.semi_major_axis_m / axes.semi_minor_axis_m - axes.aspect_ratio).abs() < 1e-12);
    assert!(
        (receipt.patch_radius_m - (axes.semi_major_axis_m * axes.semi_minor_axis_m).sqrt()).abs()
            < 1e-18
    );
    assert!(
        (receipt.pressure.resultant_n
            - 2.0 / 3.0
                * core::f64::consts::PI
                * axes.semi_major_axis_m
                * axes.semi_minor_axis_m
                * receipt.pressure.peak_pressure_pa)
            .abs()
            <= receipt.normal_force_n * 1e-12
    );
    assert!(
        (receipt.pressure.second_moment_m2
            - (axes.semi_major_axis_m.powi(2) + axes.semi_minor_axis_m.powi(2)) / 5.0)
            .abs()
            < 1e-18
    );
    assert!(
        (receipt.reversible_energy_j - 0.4 * receipt.normal_force_n * req.indentation_m).abs()
            < 1e-18
    );
    assert_eq!(receipt.irreversible_work_j, 0.0);
    assert_eq!(receipt.dissipated_power_w, 0.0);
    let h = 1.0e-8;
    let mut plus = req.clone();
    plus.indentation_m += h;
    let mut minus = req;
    minus.indentation_m -= h;
    let numerical = (point(plus.evaluate().unwrap()).unwrap().normal_force_n
        - point(minus.evaluate().unwrap()).unwrap().normal_force_n)
        / (2.0 * h);
    assert!((numerical - receipt.tangent_n_per_m).abs() / receipt.tangent_n_per_m < 2.0e-6);
}

#[test]
fn np_010_elliptic_hertz_reduces_to_sphere_and_refuses_out_of_domain_shapes() {
    let sphere = request(NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    });
    let elliptic = request(NormalPatchLaw::HertzEllipticParaboloid {
        maximum_principal_curvature_m_inverse: 50.0,
        minimum_principal_curvature_m_inverse: 50.0,
        reduced_modulus_pa: 2.0e9,
    });
    let sphere_receipt = point(sphere.evaluate().unwrap()).unwrap();
    let elliptic_receipt = point(elliptic.evaluate().unwrap()).unwrap();
    let axes = elliptic_receipt.elliptic_patch_axes.unwrap();
    for (elliptic_value, sphere_value) in [
        (
            elliptic_receipt.normal_force_n,
            sphere_receipt.normal_force_n,
        ),
        (
            elliptic_receipt.tangent_n_per_m,
            sphere_receipt.tangent_n_per_m,
        ),
        (
            elliptic_receipt.pressure.peak_pressure_pa,
            sphere_receipt.pressure.peak_pressure_pa,
        ),
        (
            elliptic_receipt.reversible_energy_j,
            sphere_receipt.reversible_energy_j,
        ),
        (axes.semi_major_axis_m, sphere_receipt.patch_radius_m),
        (axes.semi_minor_axis_m, sphere_receipt.patch_radius_m),
    ] {
        assert!((elliptic_value - sphere_value).abs() <= sphere_value.abs() * 1e-12);
    }
    assert_eq!(axes.aspect_ratio, 1.0);

    let mut reversed = elliptic.clone();
    reversed.law = NormalPatchLaw::HertzEllipticParaboloid {
        maximum_principal_curvature_m_inverse: 25.0,
        minimum_principal_curvature_m_inverse: 50.0,
        reduced_modulus_pa: 2.0e9,
    };
    assert!(matches!(
        reversed.evaluate(),
        Err(NormalPatchError::InvalidInput {
            field: "principal_curvature_order"
        })
    ));

    let mut nonpositive = elliptic.clone();
    nonpositive.law = NormalPatchLaw::HertzEllipticParaboloid {
        maximum_principal_curvature_m_inverse: 50.0,
        minimum_principal_curvature_m_inverse: 0.0,
        reduced_modulus_pa: 2.0e9,
    };
    assert!(matches!(
        nonpositive.evaluate(),
        Err(NormalPatchError::InvalidInput {
            field: "principal_curvature_or_modulus"
        })
    ));

    let mut extreme = elliptic;
    extreme.law = NormalPatchLaw::HertzEllipticParaboloid {
        maximum_principal_curvature_m_inverse: 1.0e16,
        minimum_principal_curvature_m_inverse: 1.0,
        reduced_modulus_pa: 2.0e9,
    };
    assert!(matches!(
        extreme.evaluate(),
        Err(NormalPatchError::EllipticAspectRatioUnsupported { .. })
    ));
}

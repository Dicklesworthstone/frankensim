use fs_contact::normal_patch::*;
use fs_tribo::{InputAuthority, InterfaceMedium, InterfaceSystemRef};

fn law_request(law: NormalPatchLaw) -> NormalPatchRequest {
    let geometry = match law {
        NormalPatchLaw::HertzCylinderPlane { .. } => NormalPatchGeometry::CylinderPlane,
        NormalPatchLaw::HertzSpherePlane { .. } | NormalPatchLaw::HuntCrossleySphere { .. } => {
            NormalPatchGeometry::SpherePlane
        }
        NormalPatchLaw::HertzEllipticParaboloid { .. }
        | NormalPatchLaw::HuntCrossleyEllipticParaboloid { .. } => {
            NormalPatchGeometry::EllipticParaboloid
        }
        NormalPatchLaw::FiniteGapPoint { .. } => NormalPatchGeometry::FiniteGap,
    };
    NormalPatchRequest {
        identity: NormalPatchIdentity {
            model_id: "test/normal-law-v1".into(),
            source_id: "test/synthetic-card-v1".into(),
            state_id: "caller-state-is-mapped".into(),
        },
        interface: InterfaceSystemRef::new(
            "test/actor->reactor",
            "test/history",
            "test/source",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap(),
        law,
        geometry,
        indentation_m: 9.0e-3,
        indentation_rate_m_per_s: 9.0,
        step_s: 9.0,
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

fn request(law: NormalPatchLaw, sample_id: &str, approach_m: f64) -> NormalPatchEmbedRequest {
    NormalPatchEmbedRequest {
        identity: NormalPatchEmbedIdentity {
            solver_id: "test/solver-v1".into(),
            contact_id: "test/contact-7".into(),
            feature_id: "test/feature-a".into(),
            sample_id: sample_id.into(),
        },
        lane: IntegrationLane::SmoothFixed,
        converged: true,
        kinematics: NormalPatchKinematics {
            declared_gap_m: -approach_m,
            approach_m,
            approach_rate_m_per_s: 2.0e-3,
            time_s: 1.0e-3,
            step_s: 1.0e-4,
            iteration: 1,
            normal: [0.0, 0.0, 1.0],
            moment_arm_m: [0.01, 0.0, 0.0],
        },
        law_request: law_request(law),
    }
}

#[test]
fn npe_001_manufactured_point_approach_maps_action_reaction_and_tangent() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let transition = request(law, "sample-1", 1.0e-5)
        .evaluate(&NormalPatchEmbedState::new(0.0, 0.01).unwrap())
        .unwrap();
    assert!(matches!(&transition.port, NormalPatchPort::Point(_)));
    let NormalPatchPort::Point(port) = transition.port else {
        return;
    };
    assert!(matches!(&transition.receipt, NormalPatchReceipt::Point(_)));
    let NormalPatchReceipt::Point(receipt) = transition.receipt else {
        return;
    };
    assert!((port.action_force_n[2] - receipt.normal_force_n).abs() < 1e-12);
    assert_eq!(port.reaction_force_n[2], -port.action_force_n[2]);
    assert_eq!(port.residual_force_n, [0.0; 3]);
    assert!((port.action_moment_n_m[1] + 0.01 * receipt.normal_force_n).abs() < 1e-12);
    assert_eq!(port.tangent_n_per_m, receipt.tangent_n_per_m);
    assert_eq!(transition.applicability, receipt.ratios);
    assert_eq!(transition.uncertainty, receipt.uncertainty);
}

#[test]
fn npe_002_line_port_preserves_per_length_units() {
    let law = NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let transition = request(law, "line-1", 1.0e-5)
        .evaluate(&NormalPatchEmbedState::new(0.0, 0.01).unwrap())
        .unwrap();
    assert!(matches!(&transition.port, NormalPatchPort::Line(_)));
    let NormalPatchPort::Line(port) = transition.port else {
        return;
    };
    assert!(matches!(&transition.receipt, NormalPatchReceipt::Line(_)));
    let NormalPatchReceipt::Line(receipt) = transition.receipt else {
        return;
    };
    assert_eq!(receipt.units, LINE_SI_UNITS);
    assert_eq!(
        port.action_force_n_per_m[2],
        receipt.normal_line_load_n_per_m
    );
    assert_eq!(
        port.reaction_force_n_per_m[2],
        -port.action_force_n_per_m[2]
    );
    assert_eq!(port.residual_force_n_per_m, [0.0; 3]);
    assert_eq!(port.tangent_pa, receipt.tangent_pa);
    assert!((receipt.approach_m / 1.0e-5 - 1.0).abs() < 1e-13);
}

/// G1: prescribed indentation generates the same force, tangent and energy
/// as the original load-controlled law, independently of the old load field.
#[test]
fn cylinder_mechanical_gap_drives_force_tangent_and_elastic_work() {
    let law = NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.01).unwrap();
    let evaluate = |approach| {
        let transition = request(law, "line-gap", approach).evaluate(&state).unwrap();
        let NormalPatchReceipt::Line(receipt) = transition.receipt else { panic!("line") };
        receipt
    };
    for load in [1e-12, 1e-6, 1.0, 100.0, 10_000.0] {
        let mut direct = law_request(law);
        direct.line_load_n_per_m = load;
        direct.indentation_rate_m_per_s = 0.0;
        let NormalPatchReceipt::Line(reference) = direct.evaluate().unwrap() else { panic!("line") };
        let center = evaluate(reference.approach_m);
        assert!((center.normal_line_load_n_per_m / load - 1.0).abs() < 1e-12);
        assert!((center.reversible_energy_j_per_m / reference.reversible_energy_j_per_m - 1.0).abs() < 1e-12);
        let h = center.approach_m * 1e-4;
        let plus = evaluate(center.approach_m + h);
        let minus = evaluate(center.approach_m - h);
        let tangent = (plus.normal_line_load_n_per_m - minus.normal_line_load_n_per_m) / (2.0 * h);
        let force = (plus.reversible_energy_j_per_m - minus.reversible_energy_j_per_m) / (2.0 * h);
        assert!((tangent / center.tangent_pa - 1.0).abs() < 1e-8);
        assert!((force / load - 1.0).abs() < 1e-8);
    }
    let mut sample = request(law, "line-gap", 1e-5);
    let first = sample.evaluate(&state).unwrap();
    sample.law_request.line_load_n_per_m = 1e12;
    let second = sample.evaluate(&state).unwrap();
    assert_eq!(first, second);

    let start = evaluate(1e-7);
    let end = evaluate(1e-5);
    let storage = end.reversible_energy_j_per_m - start.reversible_energy_j_per_m;
    let work = |steps| {
        let mut previous = start.clone();
        let mut sum = 0.0;
        for i in 1..=steps {
            let current = evaluate(1e-7 + (1e-5 - 1e-7) * i as f64 / steps as f64);
            sum += 0.5 * (previous.normal_line_load_n_per_m + current.normal_line_load_n_per_m)
                * (current.approach_m - previous.approach_m);
            previous = current;
        }
        sum
    };
    let coarse = (work(32) - storage).abs();
    let medium = (work(64) - storage).abs();
    let fine = (work(128) - storage).abs();
    assert!(medium < 0.3 * coarse && fine < 0.3 * medium);
    assert!(fine < 4e-6 * storage);
}

/// G0/G4: zero penetration unloads the cylinder; refused inverse samples
/// consume no work key and preserve the original material applicability limits.
#[test]
fn cylinder_gap_control_unloads_and_refuses_outside_its_source_atomically() {
    let law = NormalPatchLaw::HertzCylinderPlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.01).unwrap();
    let checkpoint = state.checkpoint();
    let mut sample = request(law, "line-retry", 0.0);
    sample.kinematics.declared_gap_m = 0.1;
    let zero = sample.evaluate(&state).unwrap();
    let NormalPatchReceipt::Line(receipt) = zero.receipt else { panic!("line") };
    assert_eq!(receipt.normal_line_load_n_per_m, 0.0);
    assert_eq!(receipt.tangent_pa, 0.0);
    assert_eq!(receipt.reversible_energy_j_per_m, 0.0);
    assert!(matches!(request(law, "line-retry", 0.04).evaluate(&state),
        Err(NormalPatchEmbedError::Law(NormalPatchError::OutsideApplicability { .. }))));
    let mut outside = request(law, "line-retry", 1e-5);
    outside.law_request.limits.max_pressure_to_yield = 1e-6;
    assert!(matches!(outside.evaluate(&state),
        Err(NormalPatchEmbedError::Law(NormalPatchError::OutsideApplicability {
            ratio: "pressure_to_yield", .. }))));
    assert_eq!(state.checkpoint(), checkpoint);
    let retry = request(law, "line-retry", 1e-5).evaluate(&state).unwrap();
    assert_eq!(retry, request(law, "line-retry", 1e-5)
        .evaluate(&NormalPatchEmbedState::rollback(&checkpoint)).unwrap());
}

#[test]
fn npe_003_power_work_is_exactly_once_and_rollback_retries() {
    let law = NormalPatchLaw::HuntCrossleySphere {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
        dissipation_s_per_m: 10.0,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.01).unwrap();
    let checkpoint = state.checkpoint();
    let sample = request(law, "work-1", 1.0e-5);
    let transition = sample.evaluate(&state).unwrap();
    assert!(matches!(&transition.port, NormalPatchPort::Point(_)));
    let NormalPatchPort::Point(port) = transition.port else {
        return;
    };
    assert!(matches!(&transition.receipt, NormalPatchReceipt::Point(_)));
    let NormalPatchReceipt::Point(receipt) = transition.receipt.clone() else {
        return;
    };
    assert_eq!(port.dissipated_power_w, receipt.dissipated_power_w);
    assert_eq!(port.irreversible_work_j, receipt.irreversible_work_j);
    assert!(
        (port.irreversible_work_j - port.dissipated_power_w * sample.kinematics.step_s).abs()
            < 1e-18
    );
    assert!(matches!(
        sample.evaluate(&transition.next_state),
        Err(NormalPatchEmbedError::DuplicateWorkKey { .. })
    ));
    let retried = sample
        .evaluate(&NormalPatchEmbedState::rollback(&checkpoint))
        .unwrap();
    assert_eq!(retried.embedding_id, transition.embedding_id);
    assert_eq!(retried.receipt_id, transition.receipt_id);
}

#[test]
fn npe_004_branch_tangent_matches_fixed_branch_difference() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.01).unwrap();
    let h = 1.0e-8;
    let center = request(law, "center", 1.0e-5).evaluate(&state).unwrap();
    let plus = request(law, "plus", 1.0e-5 + h).evaluate(&state).unwrap();
    let minus = request(law, "minus", 1.0e-5 - h).evaluate(&state).unwrap();
    assert!(matches!(&center.port, NormalPatchPort::Point(_)));
    assert!(matches!(&plus.port, NormalPatchPort::Point(_)));
    assert!(matches!(&minus.port, NormalPatchPort::Point(_)));
    let (
        NormalPatchPort::Point(center),
        NormalPatchPort::Point(plus),
        NormalPatchPort::Point(minus),
    ) = (center.port, plus.port, minus.port)
    else {
        return;
    };
    let difference = (plus.action_force_n[2] - minus.action_force_n[2]) / (2.0 * h);
    assert!((difference - center.tangent_n_per_m).abs() / center.tangent_n_per_m < 2.0e-6);
}

#[test]
fn npe_005_refuses_event_stale_future_nonconverged_and_geometry() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.001).unwrap();
    let mut eventful = request(law, "event", 1.0e-5);
    eventful.lane = IntegrationLane::Eventful;
    assert!(matches!(
        eventful.evaluate(&state),
        Err(NormalPatchEmbedError::EventfulLane)
    ));
    let mut nonconverged = request(law, "nonconverged", 1.0e-5);
    nonconverged.converged = false;
    assert!(matches!(
        nonconverged.evaluate(&state),
        Err(NormalPatchEmbedError::Nonconverged)
    ));
    let mut future = request(law, "future", 1.0e-5);
    future.kinematics.time_s = 0.01;
    assert!(matches!(
        future.evaluate(&state),
        Err(NormalPatchEmbedError::FutureState { .. })
    ));
    let accepted = request(law, "accepted", 1.0e-5).evaluate(&state).unwrap();
    assert!(matches!(
        request(law, "stale", 1.0e-5).evaluate(&accepted.next_state),
        Err(NormalPatchEmbedError::StaleState { .. })
    ));
    let mut toroidal = request(law, "toroidal", 1.0e-5);
    toroidal.law_request.geometry = NormalPatchGeometry::ToroidalOrHighlyElliptical;
    assert!(matches!(
        toroidal.evaluate(&state),
        Err(NormalPatchEmbedError::Law(
            NormalPatchError::UnsupportedGeometry { .. }
        ))
    ));
}

#[test]
fn npe_006_identity_mutation_changes_embedded_identity() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.01).unwrap();
    let base = request(law, "identity-base", 1.0e-5)
        .evaluate(&state)
        .unwrap();
    let mut changed = request(law, "identity-base", 1.0e-5);
    changed.identity.feature_id.push_str("/mutated");
    let changed = changed.evaluate(&state).unwrap();
    assert_ne!(changed.embedding_id, base.embedding_id);
    assert_ne!(changed.law_request_id, base.law_request_id);
}

#[test]
fn npe_007_distinct_samples_do_not_share_a_mapped_law_state_identity() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.01).unwrap();
    let first = request(law, "sample-a", 1.0e-5).evaluate(&state).unwrap();
    let second = request(law, "sample-b", 1.0e-5).evaluate(&state).unwrap();

    assert_ne!(first.law_request_id, second.law_request_id);
    assert_ne!(first.receipt_id, second.receipt_id);
    assert_ne!(first.embedding_id, second.embedding_id);
}

#[test]
fn npe_008_delimited_identity_fields_remain_unambiguous() {
    let law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 0.02,
        reduced_modulus_pa: 2.0e9,
    };
    let state = NormalPatchEmbedState::new(0.0, 0.01).unwrap();
    let mut first = request(law, "sample-b", 1.0e-5);
    first.identity.feature_id = "feature/a".into();
    let mut second = request(law, "a/sample-b", 1.0e-5);
    second.identity.feature_id = "feature".into();

    let first = first.evaluate(&state).unwrap();
    let second = second.evaluate(&state).unwrap();
    assert_ne!(first.law_request_id, second.law_request_id);
    assert_ne!(first.receipt_id, second.receipt_id);
    assert_ne!(first.embedding_id, second.embedding_id);
}

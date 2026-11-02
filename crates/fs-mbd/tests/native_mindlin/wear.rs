//! Actual native normal/shear motion supplies source-bound two-surface wear.
use super::*;
use fs_contact::interface_binding::mindlin::fretting::{
    SourceMindlinWear, loading::SourceHertzMindlinWear,
};
use fs_matdb::MaterialCard;
use fs_material::state_point::{resolve_interface_state_point, resolve_material_state_point};
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::mindlin::wear::*;
use fs_qty::Pressure;
use fs_tribo::WearState;

const K: [f64; 2] = [2e-6, 7e-6];
const H: [f64; 2] = [1e9, 4e8];

fn keys() -> [PropertyKey; 2] {
    ["wear-a", "wear-b"].map(|name| PropertyKey::new(name, Dims::NONE))
}
fn hardness_key() -> PropertyKey {
    PropertyKey::new("indentation-hardness", Pressure::DIMS)
}
fn claim(key: PropertyKey, value: f64) -> PropertyClaim {
    PropertyClaim {
        value: PropertyValue::Scalar {
            value,
            dims: key.quantity().dims(),
        },
        key,
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic native fretting fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn law(mu: f64, k: [f64; 2], h: [f64; 2]) -> SourceHertzMindlinWear {
    let bulk = [0, 1].map(|i| {
        let card = common::solid(["wear-A", "wear-B"][i], [1000.0, 2000.0][i], [1e7, 2e7][i]);
        let mut claims = card.claims().clone();
        claims.insert_claim(claim(hardness_key(), h[i])).unwrap();
        MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap()
    });
    let base = common::bind(&bulk[0], &bulk[1], None, 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    claims
        .insert_claim(claim(PropertyKey::new("local-friction", Dims::NONE), mu))
        .unwrap();
    for (key, value) in keys().into_iter().zip(k) {
        claims.insert_claim(claim(key, value)).unwrap();
    }
    let interface = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let normal = PointContactMaterialCards {
        surface_a: &bulk[0],
        surface_b: &bulk[1],
        interface: &interface,
        bulk_selection: [selection.clone(), selection.clone()],
        interface_selection: selection.clone(),
        model_selection: NormalContactModelSelection::SingleSupported,
        authority: InputAuthority::SyntheticFixture,
    }
    .resolve(
        &common::query(),
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.1,
        },
        common::envelope(),
    )
    .unwrap();
    let contact = SourceMindlinContact::resolve(
        normal,
        2e-5,
        ScalarPropertyRequirement::try_new(
            "local-friction",
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
        selection.clone(),
        64,
    )
    .unwrap();
    let requirements = keys().each_ref().map(|key| {
        ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
    });
    let coefficients = resolve_interface_state_point(
        &interface,
        &common::query(),
        &requirements,
        selection.clone(),
    )
    .unwrap();
    let hardness = bulk.each_ref().map(|card| {
        resolve_material_state_point(
            card,
            &common::query(),
            &[ScalarPropertyRequirement::try_with_key(
                &hardness_key(),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            selection.clone(),
        )
        .unwrap()
    });
    SourceHertzMindlinWear::from_fixed_source(
        SourceMindlinWear::new(
            contact,
            coefficients,
            keys(),
            hardness,
            [hardness_key(), hardness_key()],
        )
        .unwrap(),
    )
    .unwrap()
}
fn attachment(law: &SourceHertzMindlinWear, arms: [Vec3; 2]) -> RigidMindlinWearContact<'_> {
    RigidMindlinWearContact {
        law,
        arms_body: arms,
        frame: TangentFrame::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn near(a: f64, b: f64) {
    close(a, b, 3e-10 * a.abs().max(b.abs()).max(1e-300));
}

#[test]
fn native_recoil_and_powered_loading_accumulate_wear_without_changing_mechanics() {
    let law = law(0.6, K, H);
    let source = law.contact();
    let initial = source.initial_state();
    let gross = source.fixed_source().law().gross_slip_displacement_m();
    let loaded = source
        .advance(
            &initial,
            &TangentialStepInput {
                normal_n: initial.normal_force_n(),
                slip_increment_m: [0.7 * gross, 0.0],
                dt_s: 1e-4,
                contact_temperature_k: Some(source.temperature_kelvin()),
                nominal_pressure_pa: None,
                compliance_m_per_n: 0.0,
                solve: controls().solve,
            },
        )
        .unwrap();
    let old_volumes = [1e-15, 3e-15];
    let mut state = law
        .state_from_contact(
            loaded.state,
            old_volumes.map(|v| WearState::new(v).unwrap()),
        )
        .unwrap();
    let (mut bodies, properties, arms) = pair(initial.indentation_m(), 0.0, 0.0);
    let dt = 1e-3;
    let mut cumulative = 0.0;
    for (phase, external) in [
        [Wrench::ZERO; 2],
        [
            Wrench {
                force_world: Vec3::new(3.0, 0.5, 0.0),
                torque_body: Vec3::ZERO,
            },
            Wrench::ZERO,
        ],
        [
            Wrench {
                force_world: Vec3::new(-5.0, -0.1, 0.0),
                torque_body: Vec3::ZERO,
            },
            Wrench::ZERO,
        ],
    ]
    .into_iter()
    .enumerate()
    {
        let original = advance_rigid_mindlin_driven(
            bodies,
            properties,
            contact(source, arms),
            state.contact_state(),
            external,
            dt,
            controls(),
            &mut || false,
        )
        .unwrap();
        let candidate = advance_rigid_mindlin_driven_wear(
            bodies,
            properties,
            attachment(&law, arms),
            &state,
            external,
            dt,
            controls(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(candidate.bodies_after, original.bodies_after);
        assert_eq!(candidate.report, original.report);
        assert_eq!(
            candidate.contact.tangential.state.contact_state(),
            &original.contact.tangential.state
        );
        assert_eq!(
            candidate.contact.normal_response,
            original.contact.normal_response
        );
        assert_eq!(
            candidate.contact.normal_work_j,
            original.contact.normal_work_j
        );
        assert_eq!(
            candidate.contact.tangential.dissipated_j,
            original.contact.tangential.dissipated_j
        );
        let measure = original.contact.tangential.dissipated_j / 0.6;
        assert!(measure > 0.0);
        near(candidate.contact.pressure_weighted_sliding_j, measure);
        cumulative += measure;
        for i in 0..2 {
            near(
                candidate.contact.wear[i].increment_m3,
                K[i] * measure / H[i],
            );
            near(
                candidate.contact.tangential.state.volumes_m3()[i],
                old_volumes[i] + K[i] * cumulative / H[i],
            );
        }
        if phase == 0 {
            assert!(candidate.contact.normal_work_j < 0.0);
            assert!(candidate.contact.tangential.slip_increment_m[0] < 0.0);
            assert!(candidate.contact.tangential.stored_change_j < 0.0);
            assert!(candidate.report.midpoint_impulse_work_j > 0.0);
            let unforced = advance_rigid_mindlin_wear(
                bodies,
                properties,
                attachment(&law, arms),
                &state,
                dt,
                controls(),
                &mut || false,
            )
            .unwrap();
            assert_eq!(unforced, candidate);
        }
        bodies = candidate.bodies_after;
        state = candidate.contact.tangential.state;
    }
}

#[test]
fn native_frictionless_wear_uses_exact_normal_trajectory_and_actual_shear() {
    let law = law(0.0, K, H);
    let state = law.initial_state();
    let d0 = state.contact_state().indentation_m();
    let (bodies, properties, arms) = pair(d0, 0.01, 0.05);
    let result = advance_rigid_mindlin_wear(
        bodies,
        properties,
        attachment(&law, arms),
        &state,
        1e-3,
        controls(),
        &mut || false,
    )
    .unwrap();
    let d1 = result
        .contact
        .tangential
        .state
        .contact_state()
        .indentation_m();
    assert!(d1 > 3.0 * d0);
    let stiffness = law
        .contact()
        .fixed_source()
        .normal_source()
        .stiffness_n_per_m_3_2();
    let mean = 0.4 * stiffness * (d1.powf(2.5) - d0.powf(2.5)) / (d1 - d0);
    let travel = result.report.relative_displacement_m[1].abs();
    near(result.contact.pressure_weighted_sliding_j, mean * travel);
    assert_eq!(result.contact.tangential.dissipated_j, 0.0);
    assert_eq!(result.contact.tangential.traction_n, [0.0; 2]);
    assert!(result.contact.normal_work_j > 0.0);
    let endpoint_mean = 0.5
        * (state.contact_state().normal_force_n() + result.contact.normal_response.normal_force_n);
    assert!((endpoint_mean - mean).abs() > 0.03 * mean);
    for i in 0..2 {
        near(
            result.contact.wear[i].increment_m3,
            K[i] * mean * travel / H[i],
        );
    }
    assert_eq!(state.volumes_m3(), [0.0; 2]);
}

#[test]
fn late_wear_failure_foreign_selected_roles_and_cancel_preserve_bodies_and_history() {
    let law = law(0.6, K, H);
    let state = law.initial_state();
    let saved = state.clone();
    let (bodies, properties, arms) = pair(state.contact_state().indentation_m(), 0.02, 0.0);
    let mut calls = 0;
    let expected = advance_rigid_mindlin_wear(
        bodies,
        properties,
        attachment(&law, arms),
        &state,
        1e-3,
        controls(),
        &mut || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut reached = 0;
    assert!(matches!(
        advance_rigid_mindlin_wear(
            bodies,
            properties,
            attachment(&law, arms),
            &state,
            1e-3,
            controls(),
            &mut || {
                reached += 1;
                reached == calls
            }
        ),
        Err(RigidMindlinWearError::Cancelled)
    ));
    assert_eq!(reached, calls);
    // Equal contact physics and the same original interface still do not permit
    // replacing the selected A/B wear roles in an accepted volume checkpoint.
    let fixed = law.fixed_source();
    let reversed_keys = [keys()[1].clone(), keys()[0].clone()];
    let other = SourceHertzMindlinWear::from_fixed_source(
        SourceMindlinWear::new(
            fixed.contact().clone(),
            fixed.coefficient_source().clone(),
            reversed_keys,
            fixed.hardness_sources().clone(),
            [hardness_key(), hardness_key()],
        )
        .unwrap(),
    )
    .unwrap();
    assert!(matches!(
        advance_rigid_mindlin_wear(
            bodies,
            properties,
            attachment(&other, arms),
            &state,
            1e-3,
            controls(),
            &mut || false
        ),
        Err(RigidMindlinWearError::Wear(_))
    ));
    let huge = self::law(0.6, [K[0], 1e308], [H[0], 1e-300]);
    let huge_state = huge.initial_state();
    let saved_huge = huge_state.clone();
    // Native mechanics succeeds first. B's subsequent wear overflow must not
    // expose that body update or the already finite volume increment for A.
    assert!(
        advance_rigid_mindlin(
            bodies,
            properties,
            contact(huge.contact(), arms),
            huge_state.contact_state(),
            1e-3,
            controls(),
            &mut || false
        )
        .is_ok()
    );
    assert!(matches!(
        advance_rigid_mindlin_wear(
            bodies,
            properties,
            attachment(&huge, arms),
            &huge_state,
            1e-3,
            controls(),
            &mut || false
        ),
        Err(RigidMindlinWearError::Wear(_))
    ));
    assert_eq!(huge_state, saved_huge);
    assert_eq!(state, saved);
    assert_eq!(
        advance_rigid_mindlin_wear(
            bodies,
            properties,
            attachment(&law, arms),
            &state,
            1e-3,
            controls(),
            &mut || false
        )
        .unwrap(),
        expected
    );
}

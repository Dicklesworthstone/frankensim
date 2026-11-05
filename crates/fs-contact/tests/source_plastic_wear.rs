//! Local source-bound wear changes the actual retained finite-gap geometry.

#[path = "source_point/common.rs"]
mod common;

#[path = "source_plastic_wear/load.rs"]
mod load;

#[path = "source_plastic_wear/integration.rs"]
mod integration;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_domain;
use fs_contact::interface_binding::{
    ADHESION_ENERGY_DIMS, ADHESION_ENERGY_PROPERTY, bind_dry_interface_system_card,
    bind_isotropic_elastic_interface, bind_normal_interface_state,
};
use fs_contact::normal_patch::FiniteGapGrid;
use fs_contact::normal_patch::plastic_source::wear::{
    SourcePlasticArchardWear, SourcePlasticWearError,
};
use fs_contact::normal_patch::plastic_source::{
    PlasticContactControls, PlasticContactEnvelope, SourcePlasticContact,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    HardnessLoadStep, HardnessTestContext, InterfaceSystemCard, InterpolationPolicy, MaterialCard,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    UncertaintyModel,
};
use fs_material::hardness::{INDENTATION_HARDNESS_QUANTITY, IndentationHardnessStatePoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
    resolve_material_state_point,
};
use fs_qty::{Dims, QtyAny};
use fs_tribo::InputAuthority;

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 32,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic finite-gap wear fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn claim(key: PropertyKey, value: f64) -> PropertyClaim {
    let observations = key
        .hardness_test()
        .map(HardnessTestContext::observation)
        .into_iter()
        .collect();
    PropertyClaim {
        value: PropertyValue::Scalar {
            value,
            dims: key.quantity().dims(),
        },
        key,
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations,
        provenance: provenance(),
    }
}
fn material(name: &str, hardness: f64) -> (MaterialCard, PropertyKey) {
    let original = common::solid(name, 1000.0, 2.0e7);
    let mut claims = original.claims().clone();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("{name}/original"),
            method: "synthetic force/projected-area indentation".into(),
            artifact: hash_domain("wear-test", name.as_bytes()),
            caveats: "numerical fixture, not experimentally calibrated".into(),
            provenance: provenance(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("indentation-hardness", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "rigid sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic indentation pressure v1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims.insert_claim(claim(key.clone(), hardness)).unwrap();
    (
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap(),
        key,
    )
}
fn fixture(
    k: [f64; 2],
    h: [f64; 2],
    limit: [f64; 2],
    two_cells: bool,
    cx: &Cx<'_>,
) -> SourcePlasticArchardWear {
    let cards = [material("A", h[0]), material("B", h[1])];
    let keys = ["k_a", "k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
    let base = common::bind(&cards[0].0, &cards[1].0, None, 0.0);
    let original = base.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    for (key, value) in keys.clone().into_iter().zip(k) {
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
    let point = common::query();
    let select = MaterialPropertySelection::SingleClaimOnly;
    let solids = cards.each_ref().map(|(card, _)| {
        resolve_isotropic_solid_state_point(card, &point, select.clone()).unwrap()
    });
    let dry = bind_dry_interface_system_card(&interface, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &solids[0], &solids[1]).unwrap();
    let adhesion = resolve_interface_state_point(
        &interface,
        &point,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        select.clone(),
    )
    .unwrap();
    let normal = bind_normal_interface_state(&elastic, &adhesion).unwrap();
    let hardness = cards.each_ref().map(|(card, key)| {
        let source = resolve_material_state_point(
            card,
            &point,
            &[
                ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
                    .unwrap(),
            ],
            select.clone(),
        )
        .unwrap();
        IndentationHardnessStatePoint::try_from_resolved(source, key.clone()).unwrap()
    });
    let coefficients = resolve_interface_state_point(
        &interface,
        &point,
        &keys.each_ref().map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select,
    )
    .unwrap();
    let mut gap = vec![1.0e-3; 25];
    gap[12] = 0.0;
    if two_cells {
        gap[13] = 0.0;
    }
    let contact = SourcePlasticContact::new(
        normal,
        hardness.clone(),
        FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 1.0e-4,
            cell_depth_m: 1.0e-4,
        },
        gap,
        PlasticContactControls {
            maximum_active_set_iterations: 100,
            complementarity_tolerance_m: 1.0e-12,
            boundary_clearance_cells: 1,
        },
        PlasticContactEnvelope {
            maximum_approach_m: 1.0e-4,
            maximum_plastic_depth_m: 1.0e-4,
            half_space_depth_m: [1.0e-2; 2],
            maximum_patch_to_depth: [0.1; 2],
        },
        cx,
    )
    .unwrap();
    SourcePlasticArchardWear::new(contact, coefficients, keys, hardness, limit).unwrap()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2.0e-9 * a.abs().max(b.abs()).max(1.0e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn local_pressure_and_sliding_control_each_body_recession() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], true, cx);
        let loaded = model
            .advance_normal(&model.initial_state(), 2.0e-6, cx)
            .unwrap();
        let before = loaded.candidate.clone();
        let pressure = before.contact_state().contact_state().pressure_pa();
        assert!(pressure[12] > 0.0 && pressure[13] > 0.0);
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-4;
        slip[0] = 10.0; // Sliding on an unloaded cell removes no material.
        let result = model.advance_wear(&before, &slip, cx).unwrap();
        let measure = pressure[12] * 1.0e-4 * 1.0e-8;
        close(result.pressure_weighted_sliding_j, measure);
        for side in 0..2 {
            let expected = [0.01, 0.03][side] * pressure[12] * 1.0e-4 / [1.0e6, 2.0e6][side];
            close(result.recession_increment_m[side][12], expected);
            close(result.wear[side].increment_m3, expected * 1.0e-8);
            close(result.candidate.volumes_m3()[side], expected * 1.0e-8);
            for index in 0..25 {
                if index != 12 {
                    assert_eq!(result.candidate.recession_m()[side][index], 0.0);
                }
            }
        }
        let native = result.candidate.contact_state().contact_state();
        close(
            native.recession_m()[12],
            result.candidate.recession_m()[0][12] + result.candidate.recession_m()[1][12],
        );
        assert!(native.pressure_pa()[12] < pressure[12]);
        assert_eq!(before, loaded.candidate);
        assert_eq!(model.contact().model().undeformed_gap_datum_m(), 0.0);
    });
}

#[test]
fn removal_unloads_retained_plastic_dent_and_survives_reload() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let before = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let native_before = before.contact_state().contact_state();
        assert!(native_before.plastic_depth_m()[12] > 0.0);
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-4;
        let removed = model.advance_wear(&before, &slip, cx).unwrap();
        let native = removed.candidate.contact_state().contact_state();
        assert_eq!(native.plastic_depth_m(), native_before.plastic_depth_m());
        close(
            native.cumulative_plastic_work_j(),
            native_before.cumulative_plastic_work_j(),
        );
        assert!(removed.receipt.parameter_work_j < 0.0);
        close(
            removed.receipt.parameter_work_j,
            removed.receipt.elastic_energy_change_j,
        );
        assert_eq!(removed.receipt.response.plastic_work_j, 0.0);
        close(
            removed.receipt.elastic_energy_j - removed.receipt.previous_elastic_energy_j,
            removed.receipt.elastic_energy_change_j,
        );
        let elastic_closure = native_before.approach_m() - native_before.plastic_depth_m()[12];
        close(
            native.pressure_pa()[12],
            native_before.pressure_pa()[12] * (1.0 - native.recession_m()[12] / elastic_closure),
        );
        let unloaded = model
            .advance_normal(&removed.candidate, -1.0e-5, cx)
            .unwrap();
        assert_eq!(unloaded.receipt.normal.normal_force_n, 0.0);
        assert_eq!(
            unloaded.candidate.volumes_m3(),
            removed.candidate.volumes_m3()
        );
        assert_eq!(
            unloaded.candidate.recession_m(),
            removed.candidate.recession_m()
        );
        let unloaded_sliding = model.advance_wear(&unloaded.candidate, &slip, cx).unwrap();
        assert_eq!(
            unloaded_sliding.candidate.volumes_m3(),
            removed.candidate.volumes_m3()
        );
        assert_eq!(unloaded_sliding.receipt.parameter_work_j, 0.0);
        let reloaded = model
            .advance_normal(&unloaded.candidate, 5.0e-5, cx)
            .unwrap();
        close(
            reloaded
                .candidate
                .contact_state()
                .contact_state()
                .pressure_pa()[12],
            native.pressure_pa()[12],
        );
        assert!(
            model
                .state_from_contact(removed.candidate.contact_state().clone())
                .is_err()
        );
        assert!(
            model
                .state_from_contact(before.contact_state().clone())
                .is_ok()
        );
    });
}

#[test]
fn pressure_recession_feedback_converges_with_sliding_increment_refinement() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let initial = model
            .advance_normal(&model.initial_state(), 1.0e-6, cx)
            .unwrap()
            .candidate;
        let p0 = initial.contact_state().contact_state().pressure_pa()[12];
        let compliance = 1.0e-6 / p0;
        let distance = 2.0e-4;
        let beta = 0.01 / 1.0e6 + 0.03 / 2.0e6;
        let exact = p0 * (-beta * distance / compliance).exp();
        let integrate = |steps: usize| {
            let mut state = initial.clone();
            let mut slip = vec![0.0; 25];
            slip[12] = distance / steps as f64;
            for _ in 0..steps {
                state = model.advance_wear(&state, &slip, cx).unwrap().candidate;
            }
            state.contact_state().contact_state().pressure_pa()[12]
        };
        let coarse = integrate(8);
        let fine = integrate(64);
        assert!(coarse < fine && fine < p0);
        assert!((fine - exact).abs() < 0.2 * (coarse - exact).abs());
    });
}

#[test]
fn source_role_substitution_invalid_slip_late_bounds_and_cancellation_are_atomic() {
    with_cx(|cx, gate| {
        let model = fixture([0.01; 2], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let before = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let swapped = SourcePlasticArchardWear::new(
            model.contact().clone(),
            model.coefficient_source().clone(),
            [
                model.coefficient_keys()[1].clone(),
                model.coefficient_keys()[0].clone(),
            ],
            model.hardness_sources().clone(),
            [1.0e-4; 2],
        )
        .unwrap();
        assert_ne!(swapped.identity(), model.identity());
        assert!(swapped.advance_wear(&before, &[0.0; 25], cx).is_err());
        assert!(
            SourcePlasticArchardWear::new(
                model.contact().clone(),
                model.coefficient_source().clone(),
                model.coefficient_keys().clone(),
                [
                    model.hardness_sources()[1].clone(),
                    model.hardness_sources()[0].clone()
                ],
                [1.0e-4; 2]
            )
            .is_err()
        );
        let warmer = model
            .coefficient_source()
            .at_query_point(
                model.contact().normal_source().elastic().interface().card(),
                &QueryPoint::new().with("T", 310.0).unwrap(),
            )
            .unwrap();
        assert!(
            SourcePlasticArchardWear::new(
                model.contact().clone(),
                warmer,
                model.coefficient_keys().clone(),
                model.hardness_sources().clone(),
                [1.0e-4; 2]
            )
            .is_err()
        );
        assert!(model.advance_wear(&before, &[0.0; 24], cx).is_err());
        for invalid in [f64::NAN, f64::INFINITY, -1.0] {
            let mut slip = vec![0.0; 25];
            slip[12] = invalid;
            assert!(matches!(
                model.advance_wear(&before, &slip, cx),
                Err(SourcePlasticWearError::Invalid(_))
            ));
        }
        let bounded = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4, 1.0e-12], false, cx);
        let bounded_before = bounded
            .advance_normal(&bounded.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let bounded_saved = bounded_before.clone();
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-4;
        assert!(bounded.advance_wear(&bounded_before, &slip, cx).is_err());
        assert_eq!(bounded_before, bounded_saved);
        let zero = model.advance_wear(&before, &[0.0; 25], cx).unwrap();
        assert_eq!(zero.candidate.volumes_m3(), [0.0; 2]);
        assert_eq!(zero.candidate.recession_m(), [&[0.0; 25][..]; 2]);
        gate.request();
        assert!(model.advance_wear(&before, &slip, cx).is_err());
        assert_eq!(before, saved);
    });
}

#[test]
fn thermal_continuation_preserves_each_body_inventory_and_pinned_sources() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let loaded = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-4;
        let worn = model.advance_wear(&loaded, &slip, cx).unwrap().candidate;
        let cards = [material("A", 1.0e6).0, material("B", 2.0e6).0];
        let noop = model
            .refresh_at_temperatures(&worn, cards.each_ref(), [293.15; 2], 293.15, cx)
            .unwrap();
        assert_eq!(noop.source().identity(), model.identity());
        assert_eq!(
            noop.source().coefficient_source(),
            model.coefficient_source()
        );
        assert_eq!(noop.source().hardness_sources(), model.hardness_sources());
        let mut sub_resolution_slip = vec![0.0; 25];
        sub_resolution_slip[12] = f64::MIN_POSITIVE;
        assert!(model.advance_wear(&worn, &sub_resolution_slip, cx).is_err());
        let refreshed = model
            .refresh_at_temperatures(&worn, cards.each_ref(), [300.0, 310.0], 305.0, cx)
            .unwrap();
        let next = refreshed.candidate();
        assert_eq!(next.volumes_m3(), worn.volumes_m3());
        assert_eq!(next.recession_m(), worn.recession_m());
        assert_eq!(
            next.contact_state().contact_state().plastic_depth_m(),
            worn.contact_state().contact_state().plastic_depth_m()
        );
        assert_eq!(
            next.contact_state().contact_state().recession_m(),
            worn.contact_state().contact_state().recession_m()
        );
        for side in 0..2 {
            assert_eq!(
                refreshed.source().hardness_sources()[side].key(),
                model.hardness_sources()[side].key()
            );
            assert_eq!(
                refreshed.source().hardness_sources()[side]
                    .resolved()
                    .card_identity(),
                model.hardness_sources()[side].resolved().card_identity()
            );
        }
        assert_eq!(
            refreshed.source().coefficient_keys(),
            model.coefficient_keys()
        );
        assert_ne!(refreshed.source().identity(), model.identity());
        assert!(model.advance_wear(next, &slip, cx).is_err());
        let continued = refreshed.source().advance_wear(next, &slip, cx).unwrap();
        assert!(continued.candidate.volumes_m3()[0] > worn.volumes_m3()[0]);
        assert!(continued.candidate.volumes_m3()[1] > worn.volumes_m3()[1]);
        let saved = worn.clone();
        assert!(
            model
                .refresh_at_temperatures(&worn, cards.each_ref(), [330.0, 310.0], 305.0, cx)
                .is_err()
        );
        assert_eq!(worn, saved);
    });
}

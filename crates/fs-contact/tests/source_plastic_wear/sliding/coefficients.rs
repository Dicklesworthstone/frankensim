//! Original k(p,v) and k(T,v) use actual irreversible motion, with no calibration claim.

#[path = "coefficients/implicit.rs"]
mod implicit;

use super::*;
use fs_contact::normal_patch::plastic_source::wear::coefficients::LocalWearAxes;
use fs_qty::Temperature;

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);

fn coefficient(side: usize, first: f64, speed: f64, thermal: bool) -> f64 {
    [0.01, 0.03][side]
        * (if thermal {
            1.0 + (first - 293.15) / 20.0
        } else {
            0.5 + first / 1.0e6
        } + speed / 0.002)
}

fn local_source(
    h: [f64; 2],
    thermal: bool,
    two_cells: bool,
    cx: &Cx<'_>,
) -> SourcePlasticSlidingWear {
    let original = source(h, [1.0e-4; 2], two_cells, None, cx);
    let old_wear = original.wear_source();
    let old_normal = old_wear.contact().normal_source();
    let old_card = old_normal.elastic().interface().card();
    let mut claims = old_card.claims().clone();
    let pressure_quantity =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let speed_quantity = QuantitySpec::dimensional(SPEED);
    let keys = ["local_k_a", "local_k_b"].map(|name| PropertyKey::new(name, Dims::NONE));
    let (axis, dims, lo, hi) = if thermal {
        ("T", Temperature::DIMS, 293.15, 313.15)
    } else {
        ("wear_p", Pressure::DIMS, 1.0e5, 2.0e6)
    };
    for (side, key) in keys.iter().enumerate() {
        let mut property = claim(key.clone(), 0.01);
        property.value = PropertyValue::Surface2D {
            axes: [axis.into(), "wear_v".into()],
            axis_dims: [dims, SPEED],
            samples: [[lo, 1.0e-4], [hi, 1.0e-4], [lo, 0.01], [hi, 0.01]]
                .into_iter()
                .map(|point| (point, coefficient(side, point[0], point[1], thermal)))
                .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: Dims::NONE,
        };
        property.validity = property
            .validity
            .with_quantity("wear_p", pressure_quantity, 1.0e5, 2.0e6)
            .with_quantity("wear_v", speed_quantity, 1.0e-4, 0.01);
        property.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(property).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        old_card.surface_a().clone(),
        old_card.surface_b().clone(),
        old_card.context().clone(),
        claims,
        old_card.models().to_vec(),
    )
    .unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(
        &dry,
        old_normal.elastic().surface_a_state(),
        old_normal.elastic().surface_b_state(),
    )
    .unwrap();
    let select = MaterialPropertySelection::SingleClaimOnly;
    let query = common::query();
    let adhesion = resolve_interface_state_point(
        &card,
        &query,
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
    let config = old_wear.contact().model().config();
    let contact = SourcePlasticContact::new_with_linear_solver(
        normal,
        old_wear.contact().hardness_sources().clone(),
        config.grid,
        config.undeformed_gap_m.clone(),
        PlasticContactControls {
            maximum_active_set_iterations: config.maximum_active_set_iterations,
            complementarity_tolerance_m: config.complementarity_tolerance_m,
            boundary_clearance_cells: config.boundary_clearance_cells,
        },
        PlasticContactEnvelope {
            maximum_approach_m: 1.0e-4,
            maximum_plastic_depth_m: 1.0e-4,
            half_space_depth_m: [1.0e-2; 2],
            maximum_patch_to_depth: [0.1; 2],
        },
        old_wear.contact().model().linear_solver(),
        cx,
    )
    .unwrap();
    let coefficients = resolve_interface_state_point(
        &card,
        &query
            .clone()
            .with_quantity("wear_p", pressure_quantity, 1.0e6)
            .unwrap()
            .with_quantity("wear_v", speed_quantity, 0.001)
            .unwrap(),
        &keys.each_ref().map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select.clone(),
    )
    .unwrap();
    let axes = LocalWearAxes {
        pressure: Some("wear_p".into()),
        sliding_speed: Some("wear_v".into()),
    };
    let wear = SourcePlasticArchardWear::new_with_local_coefficients(
        contact,
        coefficients,
        keys,
        old_wear.hardness_sources().clone(),
        [1.0e-4; 2],
        axes,
    )
    .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let friction = SourceThermalCoulomb::resolve(
        dry,
        query,
        requirement("mu_static"),
        requirement("mu_kinetic"),
        select,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
    )
    .unwrap();
    SourcePlasticSlidingWear::new(wear, friction).unwrap()
}

#[test]
fn local_pressure_and_actual_speed_select_original_wear_claims() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], false, true, cx);
        let before = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let mut input = motion();
        input[12] = PlasticSlidingInput {
            released_increment_m: [1.0e-5, 0.0],
            compliance_m_per_n: 0.01,
        };
        input[13] = PlasticSlidingInput {
            released_increment_m: [2.0e-4, 0.0],
            compliance_m_per_n: 0.01,
        };
        let next = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap();
        assert_eq!(next.cells[12].response.regime, FrictionRegime::Sticking);
        let distance = 2.0e-4 - 0.01 * 0.4 * 1.0e6 * 1.0e-8;
        close(next.cells[13].response.slip_increment_m[0], distance);
        let samples = next.coefficient_samples.as_ref().unwrap();
        assert_eq!(samples.sources.len(), 1);
        assert_eq!(samples.cell_source_indices[12], None);
        assert_eq!(samples.cell_source_indices[13], Some(0));
        let source = &samples.sources[0];
        close(source.state_coordinate("wear_p").unwrap(), 1.0e6);
        close(source.state_coordinate("wear_v").unwrap(), distance / 0.1);
        assert_ne!(
            source.state_coordinate("wear_v").unwrap(),
            input[13].released_increment_m[0] / 0.1
        );
        for side in 0..2 {
            let k = coefficient(side, 1.0e6, distance / 0.1, false);
            let depth = k * 1.0e6 * distance / [1.0e6, 2.0e6][side];
            close(next.recession_increment_m[side][13], depth);
            close(next.wear[side].increment_m3, depth * 1.0e-8);
            close(
                next.candidate.wear_state().volumes_m3()[side],
                depth * 1.0e-8,
            );
            assert_eq!(next.recession_increment_m[side][12], 0.0);
            let key = &law.wear_source().coefficient_keys()[side];
            assert_eq!(
                source
                    .property_by_key(key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                law.wear_source()
                    .coefficient_source()
                    .property_by_key(key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            law.wear_source()
                .contact()
                .normal_source()
                .elastic()
                .interface()
                .card()
                .claims()
                .verify_receipt(&source.property_by_key(key).unwrap().answer().receipt)
                .unwrap();
        }
    });
}

#[test]
fn equal_queries_share_receipts_and_duration_changes_k_at_equal_distance() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], false, true, cx);
        let before = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let mut input = motion();
        input[12].released_increment_m = [1.0e-4, 0.0];
        input[13].released_increment_m = [-1.0e-4, 0.0];
        let slow = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap();
        let fast = law
            .advance(&before, &input, 0.05, friction_controls(), cx)
            .unwrap();
        for next in [&slow, &fast] {
            let samples = next.coefficient_samples.as_ref().unwrap();
            assert_eq!(samples.sources.len(), 1);
            assert_eq!(
                samples.cell_source_indices[12],
                samples.cell_source_indices[13]
            );
            close(
                next.recession_increment_m[0][12],
                next.recession_increment_m[0][13],
            );
        }
        close(
            slow.pressure_weighted_sliding_j,
            fast.pressure_weighted_sliding_j,
        );
        close(
            fast.wear[0].increment_m3 / slow.wear[0].increment_m3,
            coefficient(0, 1.0e6, 0.002, false) / coefficient(0, 1.0e6, 0.001, false),
        );
        assert!(fast.wear[0].increment_m3 > slow.wear[0].increment_m3);
        input[13].released_increment_m = [2.0e-4, 0.0];
        let mixed = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap();
        assert_eq!(mixed.coefficient_samples.as_ref().unwrap().sources.len(), 2);
    });
}

#[test]
fn each_original_body_hardness_is_the_local_wear_denominator() {
    with_cx(|cx, _| {
        let mut volumes = Vec::new();
        for hardness_b in [2.0e6, 4.0e6] {
            let law = local_source([1.0e6, hardness_b], false, false, cx);
            let before = law
                .advance_normal(&law.initial_state(), 5.0e-5, cx)
                .unwrap()
                .candidate;
            let mut input = motion();
            input[12].released_increment_m = [1.0e-4, 0.0];
            let next = law
                .advance(&before, &input, 0.1, friction_controls(), cx)
                .unwrap();
            close(next.cells[12].pressure_pa, 1.0e6);
            volumes.push([next.wear[0].increment_m3, next.wear[1].increment_m3]);
        }
        close(volumes[0][0], volumes[1][0]);
        close(volumes[0][1], 2.0 * volumes[1][1]);
    });
}

#[test]
fn zero_slip_needs_no_speed_support_and_active_out_of_domain_is_atomic() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], false, false, cx);
        let before = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let next = law
            .advance(&before, &motion(), 0.1, friction_controls(), cx)
            .unwrap();
        assert!(next.coefficient_samples.unwrap().sources.is_empty());
        assert_eq!(next.wear[0].increment_m3, 0.0);
        assert_eq!(next.wear[1].increment_m3, 0.0);
        let mut input = motion();
        input[12].released_increment_m = [1.0e-4, 0.0];
        assert!(matches!(
            law.advance(&before, &input, 0.001, friction_controls(), cx),
            Err(SourcePlasticSlidingError::Wear(
                SourcePlasticWearError::Material(_)
            ))
        ));
        assert_eq!(before, saved);
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-4;
        assert!(
            law.wear_source()
                .advance_wear(before.wear_state(), &slip, cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn interface_temperature_refresh_retains_local_axes_claims_and_histories() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], true, false, cx);
        let cards = [material("A", 1.0e6).0, material("B", 2.0e6).0];
        let before = law
            .advance_normal(&law.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let mut input = motion();
        input[12].released_increment_m = [1.0e-4, 0.0];
        let worn = law
            .advance(&before, &input, 0.1, friction_controls(), cx)
            .unwrap()
            .candidate;
        let refresh = law
            .refresh_at_temperatures(&worn, [&cards[0], &cards[1]], [293.15; 2], 313.15, cx)
            .unwrap();
        assert_eq!(
            refresh.candidate().wear_state().volumes_m3(),
            worn.wear_state().volumes_m3()
        );
        assert_eq!(
            refresh.candidate().wear_state().recession_m(),
            worn.wear_state().recession_m()
        );
        assert_eq!(
            refresh.candidate().tangential_states(),
            worn.tangential_states()
        );
        assert_eq!(
            refresh.source().wear_source().local_coefficient_axes(),
            law.wear_source().local_coefficient_axes()
        );
        let hot = refresh
            .source()
            .advance(refresh.candidate(), &input, 0.1, friction_controls(), cx)
            .unwrap();
        let cold = law
            .advance(&worn, &input, 0.1, friction_controls(), cx)
            .unwrap();
        close(hot.cells[12].pressure_pa, cold.cells[12].pressure_pa);
        let sampled = &hot.coefficient_samples.as_ref().unwrap().sources[0];
        close(sampled.state_coordinate("T").unwrap(), 313.15);
        close(sampled.state_coordinate("wear_v").unwrap(), 0.001);
        close(
            hot.wear[0].increment_m3 / cold.wear[0].increment_m3,
            coefficient(0, 313.15, 0.001, true) / coefficient(0, 293.15, 0.001, true),
        );
        for key in law.wear_source().coefficient_keys() {
            assert_eq!(
                sampled
                    .property_by_key(key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                law.wear_source()
                    .coefficient_source()
                    .property_by_key(key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    });
}

#[test]
fn local_coefficients_follow_native_maintained_load_and_held_spring_paths() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], false, false, cx);
        let load = law
            .advance_normal_at_load(&law.initial_state(), 0.003, normal_controls(), cx)
            .unwrap();
        let mut input = motion();
        input[12].released_increment_m = [1.0e-5, 0.0];
        let next = law
            .advance_at_load(
                &load.candidate,
                &input,
                0.01,
                friction_controls(),
                0.003,
                normal_controls(),
                cx,
            )
            .unwrap();
        close(
            next.coefficient_samples.as_ref().unwrap().sources[0]
                .state_coordinate("wear_p")
                .unwrap(),
            3.0e5,
        );
        close(next.normal.response.normal.normal_force_n, 0.003);
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 100.0,
            position_m: 5.0e-5,
        };
        let loaded = law
            .advance_normal_with_drive(&law.initial_state(), drive, normal_controls(), cx)
            .unwrap();
        let worn = law
            .advance_with_drive(
                &loaded.candidate,
                &input,
                0.01,
                friction_controls(),
                drive,
                normal_controls(),
                cx,
            )
            .unwrap();
        close(
            worn.coefficient_samples.as_ref().unwrap().sources[0]
                .state_coordinate("wear_p")
                .unwrap(),
            worn.cells[12].pressure_pa,
        );
        assert!(worn.wear[0].increment_m3 > 0.0);
    });
}

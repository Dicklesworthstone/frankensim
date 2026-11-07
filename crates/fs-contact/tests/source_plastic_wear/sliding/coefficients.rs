//! Original k(p,v) and k(T,v) use actual irreversible motion, with no calibration claim.

#[path = "coefficients/implicit.rs"]
mod implicit;

use super::*;
use fs_contact::normal_patch::plastic_source::wear::coefficients::LocalWearAxes;
use fs_contact::normal_patch::plastic_source::wear::integration::{
    WearIntegrationControls, WearIntegrationError,
};
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

fn adaptive_controls() -> WearIntegrationControls {
    WearIntegrationControls {
        absolute_depth_tolerance_m: 1.0e-10,
        absolute_pressure_tolerance_pa: 10.0,
        relative_tolerance: 1.0e-5,
        maximum_accepted_substeps: 4096,
        maximum_step_evaluations: 12000,
    }
}

fn cell_compliance(model: &SourcePlasticArchardWear) -> f64 {
    let width = model.contact().model().config().grid.cell_width_m;
    let modulus = model
        .contact()
        .normal_source()
        .elastic()
        .reduced_modulus_pa();
    4.0 * width * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * modulus)
}

#[test]
fn timed_adaptive_constant_load_preserves_speed_and_exact_archard_volume() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], false, false, cx);
        let model = law.wear_source();
        let force = 0.005;
        let before = model
            .advance_normal_at_load(&model.initial_state(), force, normal_controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let mut slip = [0.0; 25];
        slip[12] = 1.0e-4;
        let half_slip = slip.map(|distance| 0.5 * distance);
        let mut volumes = Vec::new();
        // Include a duration with no finite binary expansion, so elapsed-time
        // endpoint subtraction is exercised independently of dyadic slip splits.
        for duration in [0.05, 0.3] {
            let result = model
                .integrate_wear_at_load_during(
                    &before,
                    &slip,
                    duration,
                    force,
                    normal_controls(),
                    adaptive_controls(),
                    cx,
                )
                .unwrap();
            let first = model
                .integrate_wear_at_load_during(
                    &before,
                    &half_slip,
                    0.5 * duration,
                    force,
                    normal_controls(),
                    adaptive_controls(),
                    cx,
                )
                .unwrap();
            let second = model
                .integrate_wear_at_load_during(
                    &first.candidate,
                    &half_slip,
                    0.5 * duration,
                    force,
                    normal_controls(),
                    adaptive_controls(),
                    cx,
                )
                .unwrap();
            // G1/G3: constant pressure and speed give an exact independent
            // Archard integral; splitting both distance and time cannot alter k.
            let pressure = force / 1.0e-8;
            close(
                result
                    .candidate
                    .contact_state()
                    .contact_state()
                    .pressure_pa()[12],
                pressure,
            );
            close(result.receipt.pressure_weighted_sliding_j, force * slip[12]);
            for side in 0..2 {
                let k = coefficient(side, pressure, slip[12] / duration, false);
                let expected = k * force * slip[12] / [1.0e6, 2.0e6][side];
                close(result.receipt.removed_volume_m3[side], expected);
                close(second.candidate.volumes_m3()[side], expected);
                close(
                    first.receipt.removed_volume_m3[side] + second.receipt.removed_volume_m3[side],
                    expected,
                );
            }
            assert_eq!(result.receipt.accepted_substeps, 2);
            volumes.push(result.receipt.removed_volume_m3[0]);
        }
        assert!(volumes[0] > volumes[1]);
        assert_eq!(before, saved);
    });
}

#[test]
fn timed_adaptive_pressure_feedback_converges_at_held_approach_and_spring_drive() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], false, false, cx);
        let model = law.wear_source();
        let c = cell_compliance(model);
        let area = 1.0e-8;
        let p0 = 0.5e6;
        let beta = 0.01 / 1.0e6 + 0.03 / 2.0e6;
        let speed = 0.002;
        let a = 0.5 + speed / 0.002;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 2.0 * c * p0,
        };
        for spring in [false, true] {
            let before = if spring {
                model
                    .advance_normal_with_drive(&model.initial_state(), drive, normal_controls(), cx)
                    .unwrap()
                    .candidate
            } else {
                model
                    .advance_normal(&model.initial_state(), c * p0, cx)
                    .unwrap()
                    .candidate
            };
            let saved = before.clone();
            let effective_compliance = if spring { 2.0 * c } else { c };
            let mut slip = [0.0; 25];
            slip[12] = effective_compliance / (beta * a);
            let duration = slip[12] / speed;
            let solve = |controls| {
                if spring {
                    model.integrate_wear_with_drive_during(
                        &before,
                        &slip,
                        duration,
                        drive,
                        normal_controls(),
                        controls,
                        cx,
                    )
                } else {
                    model.integrate_wear_during(&before, &slip, duration, controls, cx)
                }
                .unwrap()
            };
            let loose = solve(WearIntegrationControls {
                absolute_depth_tolerance_m: 1.0e-8,
                absolute_pressure_tolerance_pa: 1000.0,
                relative_tolerance: 1.0e-3,
                ..adaptive_controls()
            });
            let fine = solve(adaptive_controls());
            // G1: k_i=k0_i*(a+p/P) gives dp/ds=-beta*p*(a+p/P)/C_eff.
            // This closed-form nonlinear decay includes actual speed in a and
            // the spring compliance A/K in C_eff, independently of the driver.
            let decay = (-beta * a * slip[12] / effective_compliance).exp();
            let exact_pressure = a * p0 * decay / (a + p0 / 1.0e6 * (1.0 - decay));
            let exact_recession = effective_compliance * (p0 - exact_pressure);
            let next = fine.candidate.contact_state().contact_state();
            let loose_pressure = loose
                .candidate
                .contact_state()
                .contact_state()
                .pressure_pa()[12];
            assert!((next.pressure_pa()[12] - exact_pressure).abs() < 0.01 * exact_pressure);
            assert!(
                (next.pressure_pa()[12] - exact_pressure).abs()
                    < (loose_pressure - exact_pressure).abs()
            );
            assert!(fine.receipt.accepted_substeps > loose.receipt.accepted_substeps);
            assert!(fine.receipt.rejected_intervals > 0);
            for side in 0..2 {
                let fraction = [0.01 / 1.0e6, 0.03 / 2.0e6][side] / beta;
                let exact_volume = area * fraction * exact_recession;
                assert!(
                    (fine.receipt.removed_volume_m3[side] - exact_volume).abs()
                        < 0.01 * exact_volume
                );
                close(
                    fine.receipt.removed_volume_m3[side],
                    fine.candidate.volumes_m3()[side],
                );
            }
            close(
                fine.candidate.volumes_m3()[1] / fine.candidate.volumes_m3()[0],
                1.5,
            );
            assert_eq!(
                next.plastic_depth_m(),
                before.contact_state().contact_state().plastic_depth_m()
            );
            assert_eq!(fine.receipt.plastic_work_j, 0.0);
            if spring {
                let initial_extension =
                    drive.position_m - before.contact_state().contact_state().approach_m();
                let final_extension = drive.position_m - next.approach_m();
                close(
                    fine.receipt.normal_work_j,
                    0.5 * drive.stiffness_n_per_m
                        * (initial_extension * initial_extension
                            - final_extension * final_extension),
                );
            } else {
                assert_eq!(fine.receipt.normal_work_j, 0.0);
            }
            assert_eq!(fine.source_identity, model.identity());
            assert_eq!(before, saved);
        }
    });
}

#[test]
fn timed_adaptive_duration_limits_zero_paths_and_failures_are_transactional() {
    with_cx(|cx, _| {
        let law = local_source([1.0e6, 2.0e6], false, false, cx);
        let model = law.wear_source();
        let force = 0.005;
        let before = model
            .advance_normal_at_load(&model.initial_state(), force, normal_controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let stiffness = 1.0e-8 / cell_compliance(model);
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: stiffness,
            position_m: before.contact_state().contact_state().approach_m() + force / stiffness,
        };
        let smallest_duration = f64::from_bits(1);
        for duration in [0.1, smallest_duration, 0.0, -1.0, f64::NAN, f64::INFINITY] {
            let responses = [
                model.integrate_wear_during(&before, &[0.0; 25], duration, adaptive_controls(), cx),
                model.integrate_wear_at_load_during(
                    &before,
                    &[0.0; 25],
                    duration,
                    force,
                    normal_controls(),
                    adaptive_controls(),
                    cx,
                ),
                model.integrate_wear_with_drive_during(
                    &before,
                    &[0.0; 25],
                    duration,
                    drive,
                    normal_controls(),
                    adaptive_controls(),
                    cx,
                ),
            ];
            for response in responses {
                if duration.is_finite() && duration > 0.0 {
                    let result = response.unwrap();
                    assert_eq!(result.candidate, before);
                    assert_eq!(result.receipt.removed_volume_m3, [0.0; 2]);
                    assert_eq!(result.receipt.accepted_substeps, 1);
                } else {
                    assert!(matches!(response, Err(WearIntegrationError::Invalid(_))));
                }
            }
        }
        let mut slip = [0.0; 25];
        slip[12] = 1.0e-4;
        assert!(matches!(
            model.integrate_wear_during(&before, &slip, smallest_duration, adaptive_controls(), cx),
            Err(WearIntegrationError::ResolutionLimit { .. })
        ));
        for duration in [0.001, 10.0] {
            assert!(matches!(
                model.integrate_wear_during(&before, &slip, duration, adaptive_controls(), cx),
                Err(WearIntegrationError::Step(
                    SourcePlasticWearError::Material(_)
                ))
            ));
        }
        // The untimed API still requires a rate-independent law. A supplied
        // duration enables the speed-dependent law; it never fills a missing one.
        assert!(matches!(
            model.integrate_wear(&before, &slip, adaptive_controls(), cx),
            Err(WearIntegrationError::Step(SourcePlasticWearError::Invalid(
                _
            )))
        ));
        slip[12] = cell_compliance(model) / ((0.01 / 1.0e6 + 0.03 / 2.0e6) * 1.5);
        assert!(matches!(
            model.integrate_wear_during(
                &before,
                &slip,
                slip[12] / 0.002,
                WearIntegrationControls {
                    maximum_accepted_substeps: 2,
                    ..adaptive_controls()
                },
                cx,
            ),
            Err(WearIntegrationError::BudgetExceeded {
                accepted_substeps: 2,
                ..
            })
        ));
        assert_eq!(before, saved);
    });
}

#[test]
fn timed_adaptive_refinement_preserves_supported_speed_endpoints() {
    with_cx(|cx, _| {
        for speed in [1.0e-4_f64, 0.01] {
            for duration in [0.3_f64, 1.7] {
                for just_inside in [false, true] {
                    // Keep pressure strictly inside its source support for the
                    // full path, including the larger upper-bound speed.
                    let hardness = if speed == 0.01 {
                        [1.0e9, 2.0e9]
                    } else {
                        [1.0e6, 2.0e6]
                    };
                    let law = local_source(hardness, false, false, cx);
                    let model = law.wear_source();
                    let c = cell_compliance(model);
                    let force = 0.005;
                    let p0 = force / 1.0e-8;
                    let drive = FiniteGapPlasticDrive {
                        stiffness_n_per_m: 1.0e-8 / c,
                        position_m: 2.0 * c * p0,
                    };
                    let mut slip = [0.0; 25];
                    slip[12] = if !just_inside {
                        speed * duration
                    } else if speed == 0.01 {
                        (speed * duration).next_down()
                    } else {
                        (speed * duration).next_up()
                    };
                    assert!((1.0e-4..=0.01).contains(&(slip[12] / duration)));
                    for boundary in 0..3 {
                        let before = match boundary {
                            0 => {
                                model
                                    .advance_normal(&model.initial_state(), c * p0, cx)
                                    .unwrap()
                                    .candidate
                            }
                            1 => {
                                model
                                    .advance_normal_at_load(
                                        &model.initial_state(),
                                        force,
                                        normal_controls(),
                                        cx,
                                    )
                                    .unwrap()
                                    .candidate
                            }
                            _ => {
                                model
                                    .advance_normal_with_drive(
                                        &model.initial_state(),
                                        drive,
                                        normal_controls(),
                                        cx,
                                    )
                                    .unwrap()
                                    .candidate
                            }
                        };
                        let saved = before.clone();
                        let coarse = match boundary {
                            0 => model.advance_wear_during(&before, &slip, duration, cx),
                            1 => model
                                .advance_wear_at_load_during(
                                    &before,
                                    &slip,
                                    duration,
                                    force,
                                    normal_controls(),
                                    cx,
                                )
                                .map(|next| next.wear),
                            _ => model
                                .advance_wear_with_drive_during(
                                    &before,
                                    &slip,
                                    duration,
                                    drive,
                                    normal_controls(),
                                    cx,
                                )
                                .map(|next| next.wear),
                        };
                        // Each complete mean-speed query is admitted by the
                        // unmodified single-step path at the original bounds.
                        assert!(coarse.is_ok(), "{coarse:?}");
                        let refined = match boundary {
                            0 => model.integrate_wear_during(&before, &slip, duration, adaptive_controls(), cx),
                            1 => model.integrate_wear_at_load_during(&before, &slip, duration, force, normal_controls(), adaptive_controls(), cx),
                            _ => model.integrate_wear_with_drive_during(&before, &slip, duration, drive, normal_controls(), adaptive_controls(), cx),
                        }.unwrap_or_else(|error| panic!(
                            "supported speed lost: boundary={boundary}, speed={speed:e}, duration={duration}, inside={just_inside}: {error}"
                        ));
                        if boundary != 1 {
                            assert!(refined.receipt.accepted_substeps > 2);
                        }
                        assert!(
                            refined
                                .candidate
                                .contact_state()
                                .contact_state()
                                .pressure_pa()[12]
                                > 1.0e5
                        );
                        assert_eq!(refined.source_identity, model.identity());
                        assert_eq!(before, saved);
                    }
                }
            }
        }
    });
}

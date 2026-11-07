//! G1/G3/G4: evolving original normal loads with retained shear and finite heat.
use super::*;

fn normal(load: f64) -> NormalPatchView {
    let original = reference(1.0, PartialSlipState::zero());
    let n = original.normal_patch();
    NormalPatchView::new(
        n.patch_id(),
        n.normal_model_id(),
        n.source_id(),
        n.authority(),
        load,
        n.semi_axis_longitudinal_m(),
        n.semi_axis_lateral_m(),
        n.pressure_second_moment_m2(),
    )
    .unwrap()
}

fn uniform(load: f64) -> MovingPressureProfile {
    MovingPressureProfile::new(
        &normal(load),
        vec![MovingPressureCell {
            span_m: [0.0, 1.0],
            normal_pressure_pa: load,
        }],
    )
    .unwrap()
}

fn loaded_step(
    contact: &MovingPartialSlip,
    thermal: &DepositionSimulation<'static>,
    target: f64,
    pressure: &MovingPressureProfile,
) -> MovingPartialSlipCandidate<'static, ()> {
    contact
        .step_candidate_with_pressure_profile(
            thermal,
            target,
            contact.position(),
            pressure,
            &work(),
            controls(),
            &vec![0.0; thermal.cells().len()],
            300.0,
            |t, regions, _, _| Ok::<_, Infallible>(source(t, regions, 1.0)),
            || false,
        )
        .unwrap()
}

#[test]
fn unloading_releases_original_shear_and_reloading_preserves_the_remaining_strains() {
    let (thermal, contact) = fixture(&[-3.0, 5.0], 1.0, [0.004, 0.01]);
    let old = contact.parcels().to_vec();
    let low = uniform(0.2);
    let unloaded = loaded_step(&contact, &thermal, 0.1, &low);
    let r = unloaded.report();
    // Independent return: k=[60,40], new capacities=[0.06,0.14].
    // e_new=[0.001,0.0035], D=sum(C*(e_old-e_new)),
    // D_num=sum(k*(e_old-e_new)^2/2). Normal displacement is not prescribed.
    near(r.force_on_a_n[0], -0.2, 2e-14);
    near(r.storage_before_j, 0.00248, 2e-14);
    near(r.storage_after_j, 0.000275, 2e-14);
    near(r.heat.dissipated_j, 0.00109, 2e-14);
    near(r.numerical_dissipation_j, 0.001115, 2e-14);
    near(r.heat.thermal.stored_energy_change_j, 0.00109, 3e-10);
    near(r.relative_slip_m, 0.0, 0.0);
    near(r.mechanical_work_j, 0.0, 0.0);
    near(r.elastic_outflow_j, 0.0, 0.0);
    near(r.total_energy_residual_j, 0.0, 5e-10);
    assert_eq!(r.normal_load_before_n, 1.0);
    assert_eq!(r.normal_load_after_n, 0.2);
    assert_eq!(unloaded.contact().normal_patch(), low.normal_patch());
    assert_eq!(contact.parcels(), old);
    assert!(contact.pressure_profile().is_none());
    let (unloaded, heated, _) = unloaded.into_parts();
    let strain = unloaded.parcels()[0]
        .state()
        .region_elastic_displacement_m();
    let branches = unloaded.parcels()[0].state().sliding_regions();
    let reloaded = loaded_step(&unloaded, &heated, 0.2, &uniform(2.0));
    near(reloaded.report().force_on_a_n[0], -0.2, 2e-14);
    near(reloaded.report().heat.dissipated_j, 0.0, 0.0);
    near(reloaded.report().numerical_dissipation_j, 0.0, 0.0);
    near(reloaded.report().storage_after_j, 0.000275, 2e-14);
    assert_eq!(reloaded.report().normal_load_before_n, 0.2);
    assert_eq!(reloaded.report().normal_load_after_n, 2.0);
    assert_eq!(
        reloaded.contact().parcels()[0]
            .state()
            .region_elastic_displacement_m(),
        strain
    );
    assert_eq!(
        reloaded.contact().parcels()[0].state().sliding_regions(),
        branches
    );
    // The ordinary continuation now uses the accepted 2 N receipt, without
    // requiring the caller to resupply it or rebuilding the original brush.
    let (reloaded, heated, _) = reloaded.into_parts();
    let continued = step(&reloaded, &heated, 0.3, position(0.0, -0.001, 0.0), 1.0);
    near(continued.report().force_on_a_n[0], -0.3, 2e-14);
    near(continued.report().heat.dissipated_j, 0.0, 0.0);
    near(continued.report().storage_after_j, 0.000525, 2e-14);
    near(continued.report().mechanical_work_j, 0.0003, 2e-14);
    assert_eq!(continued.report().normal_load_before_n, 2.0);
    assert_eq!(continued.report().normal_load_after_n, 2.0);
    assert_eq!(
        continued
            .contact()
            .reference()
            .normal_patch()
            .normal_load_n(),
        1.0
    );
}

#[test]
fn pressure_redistribution_at_constant_load_releases_only_the_unloaded_material() {
    let (thermal, contact) = fixture(&[-3.0, 0.0, 0.5, 1.0, 5.0], 1.0, [0.004, 0.01]);
    let pressure = MovingPressureProfile::new(
        &normal(1.0),
        vec![
            MovingPressureCell {
                span_m: [0.0, 0.5],
                normal_pressure_pa: 0.4,
            },
            MovingPressureCell {
                span_m: [0.5, 1.0],
                normal_pressure_pa: 1.6,
            },
        ],
    )
    .unwrap();
    let next = loaded_step(&contact, &thermal, 0.1, &pressure);
    let r = next.report();
    // The low-pressure half returns to [0.002,0.007] m. The other half
    // retains [0.004,0.01] m. No additional tangential displacement occurs.
    near(r.force_on_a_n[0], -0.52, 2e-14);
    near(r.storage_after_j, 0.00179, 2e-14);
    near(r.heat.dissipated_j, 0.00054, 2e-14);
    near(r.numerical_dissipation_j, 0.00015, 2e-14);
    near(r.total_energy_residual_j, 0.0, 5e-10);
    near(
        r.regions.iter().map(|s| s.region.normal_load_n).sum(),
        1.0,
        2e-14,
    );
    for region in &r.regions {
        assert_eq!(region.region.accepted_state(), contact.parcels()[0].state());
        if region.region.material_span_m[0] >= 0.5 {
            assert_eq!(
                region.response.next_state.region_elastic_displacement_m(),
                [[0.004, 0.0], [0.01, 0.0]]
            );
            near(region.response.dissipation.heat_j, 0.0, 0.0);
        }
    }
    near(
        r.heat
            .patches
            .iter()
            .filter(|p| p.cells[0] == CellId(1))
            .map(|p| p.dissipated_j)
            .sum(),
        0.00054,
        2e-14,
    );
    near(
        r.heat
            .patches
            .iter()
            .filter(|p| p.cells[0] == CellId(2))
            .map(|p| p.dissipated_j)
            .sum(),
        0.0,
        0.0,
    );
}

#[test]
fn gradual_normal_unloading_refines_physical_heat_and_separate_return_loss() {
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let (mut thermal, mut contact) = fixture(&[-3.0, 5.0], 1.0, [0.004, 0.01]);
        let mut heat = 0.0;
        let mut numerical = 0.0;
        for i in 1..=steps {
            let t = i as f64 / steps as f64;
            let next = loaded_step(&contact, &thermal, t, &uniform(1.0 - 0.8 * t));
            let (c, h, r) = next.into_parts();
            heat += r.heat.dissipated_j;
            numerical += r.numerical_dissipation_j;
            near(r.mechanical_work_j, 0.0, 0.0);
            near(r.total_energy_residual_j, 0.0, 5e-10);
            contact = c;
            thermal = h;
        }
        // Under quasi-static pressure release at held total displacement,
        // integral(C de_plastic) equals the lost original spring energy.
        let exact_loss = 0.00248 - 0.000275;
        near(contact.stored_energy_j().unwrap(), 0.000275, 2e-14);
        near(heat + numerical, exact_loss, 3e-14);
        assert!(heat > 0.0 && numerical > 0.0);
        errors.push(exact_loss - heat);
    }
    assert!(
        errors[1] < 0.55 * errors[0] && errors[2] < 0.55 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn normal_continuation_refuses_foreign_lineage_geometry_or_load_before_source() {
    let (thermal, contact) = fixture(&[-3.0, 5.0], 1.0, [0.004, 0.01]);
    let n = normal(0.2);
    let calls = Cell::new(0);
    for variant in 0..7 {
        let changed = NormalPatchView::new(
            if variant == 0 {
                "foreign patch"
            } else {
                n.patch_id()
            },
            if variant == 1 {
                "foreign model"
            } else {
                n.normal_model_id()
            },
            if variant == 2 {
                "foreign source"
            } else {
                n.source_id()
            },
            if variant == 3 {
                NormalPatchAuthority::Estimated
            } else {
                n.authority()
            },
            n.normal_load_n(),
            n.semi_axis_longitudinal_m() * if variant == 4 { 2.0 } else { 1.0 },
            n.semi_axis_lateral_m(),
            n.pressure_second_moment_m2() * if variant == 5 { 2.0 } else { 1.0 },
        )
        .unwrap();
        let pressure = MovingPressureProfile::new(
            &changed,
            vec![MovingPressureCell {
                span_m: [0.0, 1.0],
                normal_pressure_pa: if variant == 6 { 0.3 } else { 0.2 },
            }],
        )
        .unwrap();
        let result = contact.step_candidate_with_pressure_profile(
            &thermal,
            0.1,
            contact.position(),
            &pressure,
            &work(),
            controls(),
            &[0.0; 2],
            300.0,
            |t, r, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, Infallible>(source(t, r, 1.0))
            },
            || false,
        );
        assert!(result.is_err(), "variant {variant}");
    }
    assert_eq!(calls.get(), 0);
    assert!(contact.pressure_profile().is_none());
    near(contact.stored_energy_j().unwrap(), 0.00248, 2e-14);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn changed_pressure_late_source_failure_cancellation_and_budget_leave_a_retryable_state() {
    let (thermal, contact) = fixture(&[-3.0, 5.0], 1.0, [0.004, 0.01]);
    let pressure = uniform(0.2);
    let before = contact.parcels().to_vec();
    let calls = Cell::new(0);
    let failed = contact.step_candidate_with_pressure_profile(
        &thermal,
        0.1,
        contact.position(),
        &pressure,
        &work(),
        controls(),
        &[0.0; 2],
        300.0,
        |t, r, _, _| {
            calls.set(calls.get() + 1);
            if calls.get() >= 3 {
                Err("late normal-path refusal")
            } else {
                Ok(source(t, r, 1.0))
            }
        },
        || false,
    );
    assert!(matches!(
        failed,
        Err(MovingPartialSlipError::Source("late normal-path refusal"))
    ));
    let cancelled = contact.step_candidate_with_pressure_profile(
        &thermal,
        0.1,
        contact.position(),
        &pressure,
        &work(),
        controls(),
        &[0.0; 2],
        300.0,
        |t, r, _, _| Ok::<_, Infallible>(source(t, r, 1.0)),
        || true,
    );
    assert!(matches!(
        cancelled,
        Err(MovingPartialSlipError::Thermal(DepositionError::Cancelled))
    ));
    let split = MovingPressureProfile::new(
        &normal(0.2),
        vec![
            MovingPressureCell {
                span_m: [0.0, 0.5],
                normal_pressure_pa: 0.2,
            },
            MovingPressureCell {
                span_m: [0.5, 1.0],
                normal_pressure_pa: 0.2,
            },
        ],
    )
    .unwrap();
    let mut limit = controls();
    limit.maximum_regions = 1;
    let limited = contact.step_candidate_with_pressure_profile(
        &thermal,
        0.1,
        contact.position(),
        &split,
        &work(),
        limit,
        &[0.0; 2],
        300.0,
        |_, _, _, _| -> Result<MovingPartialSlipSource<()>, Infallible> {
            panic!("pressure budget must precede source")
        },
        || false,
    );
    assert!(matches!(limited, Err(MovingPartialSlipError::RegionBudget)));
    assert_eq!(contact.parcels(), before);
    assert_eq!(contact.normal_patch().normal_load_n(), 1.0);
    assert!(contact.pressure_profile().is_none());
    assert_eq!(thermal.time_s(), 0.0);
    assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    let retry = loaded_step(&contact, &thermal, 0.1, &pressure);
    near(retry.report().heat.dissipated_j, 0.00109, 2e-14);
    assert_eq!(retry.contact().pressure_profile(), Some(&pressure));
}

//! Supplied spatial pressure, material memory and original finite heat.
use super::*;

fn cell(a: f64, b: f64, pressure: f64) -> MovingPressureCell {
    MovingPressureCell {
        span_m: [a, b],
        normal_pressure_pa: pressure,
    }
}

fn edge_core() -> Vec<MovingPressureCell> {
    // The two edges carry 0.1 N each; the centre carries 0.8 N.
    vec![
        cell(0.0, 0.25, 0.4),
        cell(0.25, 0.75, 1.6),
        cell(0.75, 1.0, 0.4),
    ]
}

fn loaded(
    edges: &[f64],
    mu: f64,
    state: PartialSlipState,
    cells: Vec<MovingPressureCell>,
) -> (DepositionSimulation<'static>, MovingPartialSlip) {
    let (thermal, route) = fixture_route(edges);
    let original = reference(mu, state);
    let pressure = MovingPressureProfile::new(original.normal_patch(), cells).unwrap();
    let contact = MovingPartialSlip::new_with_pressure_profile(
        &thermal,
        route,
        original,
        position(0.0, 0.0, 0.0),
        pressure,
        || Ok(()),
    )
    .unwrap();
    (thermal, contact)
}

#[test]
fn uniform_pressure_cells_reproduce_the_original_force_history_and_heat() {
    let edges = [-3.0, 0.0, 0.4, 0.8, 1.0, 4.0];
    let (thermal, uniform) = fixture(&edges, 1.0, [0.0; 2]);
    let end = position(0.02, 0.0, 0.02);
    let baseline = step(&uniform, &thermal, 0.1, end, 1.0);
    for cells in [
        vec![cell(0.0, 1.0, 1.0)],
        vec![
            cell(0.0, 0.25, 1.0),
            cell(0.25, 0.75, 1.0),
            cell(0.75, 1.0, 1.0),
        ],
    ] {
        let (thermal, contact) = loaded(&edges, 1.0, initial([0.0; 2]), cells.clone());
        assert_eq!(contact.pressure_profile().unwrap().cells(), cells);
        assert_eq!(
            contact.pressure_profile().unwrap().normal_patch(),
            contact.reference().normal_patch()
        );
        let candidate = step(&contact, &thermal, 0.1, end, 1.0);
        let a = baseline.report();
        let b = candidate.report();
        assert_eq!(baseline.contact().parcels(), candidate.contact().parcels());
        for j in 0..3 {
            near(a.force_on_a_n[j], b.force_on_a_n[j], 2e-13);
            near(
                a.moment_on_a_about_origin_nm[j],
                b.moment_on_a_about_origin_nm[j],
                2e-13,
            );
        }
        for (x, y) in [
            (a.mechanical_work_j, b.mechanical_work_j),
            (a.storage_after_j, b.storage_after_j),
            (a.numerical_dissipation_j, b.numerical_dissipation_j),
            (a.elastic_outflow_j, b.elastic_outflow_j),
            (a.heat.dissipated_j, b.heat.dissipated_j),
        ] {
            near(x, y, 2e-13);
        }
        near(
            b.regions.iter().map(|r| r.region.normal_load_n).sum(),
            1.0,
            2e-14,
        );
        for r in &b.regions {
            assert_eq!(r.region.normal_pressure_pa, 1.0);
            assert_eq!(
                r.response.checkpoint.normal_patch(),
                contact.reference().normal_patch()
            );
        }
        for (id, before) in baseline.thermal().cells() {
            near(
                before.specific_h_j_kg,
                candidate.thermal().cells()[id].specific_h_j_kg,
                2e-10,
            );
        }
    }
}

#[test]
fn pressure_edges_microslip_and_heat_while_the_core_sticks() {
    let edges = [-2.0, 0.0, 0.25, 0.75, 1.0, 3.0];
    let (thermal, contact) = loaded(&edges, 1.0, initial([0.0; 2]), edge_core());
    let next = step(&contact, &thermal, 0.1, position(0.0, -0.004, 0.0), 1.0);
    let r = next.report();
    assert_eq!(r.regions.len(), 3);
    for (i, region) in r.regions.iter().enumerate() {
        assert_eq!(region.region.normal_pressure_pa, [0.4, 1.6, 0.4][i]);
        near(region.region.normal_load_n, [0.1, 0.8, 0.1][i], 1e-15);
        assert_eq!(
            region.response.state,
            if i == 1 {
                PartialSlipStateKind::Sticking
            } else {
                PartialSlipStateKind::PartialSlip
            }
        );
        assert_eq!(
            region.response.checkpoint.normal_patch().authority(),
            contact.reference().normal_patch().authority()
        );
    }
    // Independent two-spring calculation: k=[60,40]. Edge capacities are
    // [0.12,0.28] N, so the weak spring returns to 0.002 m; the core stays
    // elastic at 0.004 m. Edge area is 1/2 and core area is 1/2.
    near(r.force_on_a_n[0], -0.34, 2e-13);
    near(r.mechanical_work_j, 0.00136, 2e-14);
    near(r.storage_after_j, 0.00062, 2e-14);
    near(r.numerical_dissipation_j, 0.00062, 2e-14);
    near(r.heat.dissipated_j, 0.00012, 2e-14);
    near(r.heat.thermal.stored_energy_change_j, 0.00012, 3e-10);
    near(r.total_energy_residual_j, 0.0, 5e-10);
    for (id, expected) in [(CellId(1), 0.00006), (CellId(2), 0.0), (CellId(3), 0.00006)] {
        near(
            r.heat
                .patches
                .iter()
                .filter(|p| p.cells[0] == id)
                .map(|p| p.dissipated_j)
                .sum(),
            expected,
            2e-14,
        );
    }
    let (_, uniform) = fixture(&edges, 1.0, [0.0; 2]);
    let old = step(&uniform, &thermal, 0.1, position(0.0, -0.004, 0.0), 1.0);
    near(old.report().force_on_a_n[0], -0.4, 2e-13);
    near(old.report().heat.dissipated_j, 0.0, 2e-14);
}

#[test]
fn material_memory_crosses_pressure_boundaries_without_reset_or_thermalizing_outflow() {
    let (thermal, contact) = loaded(
        &[-2.0, 0.0, 0.25, 0.75, 1.0, 3.0],
        1.0,
        initial([0.0; 2]),
        edge_core(),
    );
    let first = step(&contact, &thermal, 0.1, position(0.0, -0.004, 0.0), 1.0);
    let (contact, thermal, _) = first.into_parts();
    let old_right = contact
        .parcels()
        .iter()
        .find(|p| p.span_m == [0.75, 1.0])
        .unwrap()
        .state()
        .clone();
    // No additional relative travel. The footprint advances over A material;
    // a previously elastic piece enters the low-pressure leading region, and
    // a previously microslipping piece enters the high-pressure core.
    let next = step(&contact, &thermal, 0.2, position(0.0, -0.004, 0.25), 1.0);
    let r = next.report();
    near(r.relative_slip_m, 0.0, 0.0);
    near(r.mechanical_work_j, 0.0, 2e-14);
    near(r.elastic_outflow_j, 0.00011, 2e-14);
    near(r.storage_before_j, 0.00062, 2e-14);
    near(r.storage_after_j, 0.00042, 2e-14);
    // The lowered local normal bound releases weak-spring energy. Its real
    // plastic part alone heats the cells; return-map loss and outflow remain
    // separate, independently closed mechanical accounts.
    near(r.heat.dissipated_j, 0.00006, 2e-14);
    near(r.numerical_dissipation_j, 0.00003, 2e-14);
    near(r.heat.thermal.stored_energy_change_j, 0.00006, 3e-10);
    near(r.total_energy_residual_j, 0.0, 5e-10);
    let carried = r
        .regions
        .iter()
        .find(|r| r.region.material_span_m == [0.75, 1.0])
        .unwrap();
    assert_eq!(carried.region.accepted_state(), &old_right);
    assert_eq!(
        carried.response.next_state.region_elastic_displacement_m(),
        old_right.region_elastic_displacement_m()
    );
    assert_eq!(
        carried.response.next_state.sliding_regions(),
        old_right.sliding_regions()
    );
    assert_eq!(carried.region.normal_pressure_pa, 1.6);
    let entering = r
        .regions
        .iter()
        .find(|r| r.region.material_span_m == [1.0, 1.25])
        .unwrap();
    assert_eq!(entering.region.accepted_state(), &PartialSlipState::zero());
    assert_eq!(entering.region.normal_pressure_pa, 0.4);
    assert_eq!(
        next.contact().pressure_profile(),
        contact.pressure_profile()
    );
}

#[test]
fn local_pressure_and_endpoint_temperature_feed_the_original_friction_source() {
    let (thermal, contact) = loaded(&[-3.0, 5.0], 1.0, initial([0.0; 2]), edge_core());
    let dt = 0.02;
    let capacity = 800.0;
    let beta = 10000.0;
    let qa = 2.0 / (1.0 + dt * 2.0 / capacity);
    let qb = 6.0 / (1.0 + dt * 6.0 / capacity);
    let share = qa / (qa + qb);
    let response = |heat: f64| {
        let temperature = 300.0 + share * heat / capacity;
        edge_core().iter().fold((0.0, 0.0), |(force, loss), c| {
            let fraction = c.span_m[1] - c.span_m[0];
            let q = c.normal_pressure_pa
                / (1.0 + 0.25 * c.normal_pressure_pa + beta * (temperature - 300.0));
            (
                force + fraction * q,
                loss + fraction * (0.1 * q - 0.01375 * q * q),
            )
        })
    };
    // Independent scalar energy balance: Ta=300+share*D/C and each fully
    // sliding pair of original springs gives D_i=q_i*du-0.01375*q_i^2.
    let mut lo = 0.0;
    let mut hi = 0.1;
    for _ in 0..100 {
        let middle = f64::midpoint(lo, hi);
        if response(middle).1 > middle {
            lo = middle;
        } else {
            hi = middle;
        }
    }
    let heat = f64::midpoint(lo, hi);
    let force = response(heat).0;
    let next = contact
        .step_candidate(
            &thermal,
            dt,
            position(0.0, -0.1, 0.0),
            &work(),
            controls(),
            &[0.0; 2],
            300.0,
            |t, regions, _, _| {
                let mut s = source(t, regions, 1.0);
                s.laws = regions
                    .iter()
                    .map(|r| {
                        let temperature = t.temperature_k(r.surface_a_cell).unwrap();
                        law(1.0
                            / (1.0 + 0.25 * r.normal_pressure_pa + beta * (temperature - 300.0)))
                    })
                    .collect();
                Ok::<_, Infallible>(s)
            },
            || false,
        )
        .unwrap();
    near(next.report().force_on_a_n[0], -force, 3e-8);
    near(next.report().heat.dissipated_j, heat, 4e-9);
    near(
        next.thermal().temperature_k(CellId(0)).unwrap(),
        300.0 + share * heat / capacity,
        3e-11,
    );
    near(
        next.report().heat.thermal.stored_energy_change_j,
        heat,
        4e-9,
    );
    near(next.report().total_energy_residual_j, 0.0, 5e-10);
    assert!(
        force < 0.95 * response(0.0).0,
        "temperature feedback must change actual friction"
    );
}

#[test]
fn pressure_coverage_load_and_source_authority_refuse_without_repair() {
    let (thermal, route) = fixture_route(&[-3.0, 5.0]);
    let original = reference(1.0, initial([0.0; 2]));
    for cells in [
        vec![],
        vec![cell(0.0, 1.0, 0.0)],
        vec![cell(0.0, 1.0, f64::NAN)],
        vec![cell(0.0, 0.4, 1.0), cell(0.5, 1.0, 1.0)],
        vec![cell(0.0, 0.6, 1.0), cell(0.5, 1.0, 1.0)],
    ] {
        assert!(matches!(
            MovingPressureProfile::new(original.normal_patch(), cells),
            Err(MovingPartialSlipError::Invalid(_))
        ));
    }
    for cells in [
        vec![cell(0.0, 1.0, 1.1)],
        vec![cell(0.1, 0.9, 1.25)],
        vec![cell(-0.5, 1.5, 0.5)],
    ] {
        let profile = MovingPressureProfile::new(original.normal_patch(), cells).unwrap();
        let rejected = MovingPartialSlip::new_with_pressure_profile(
            &thermal,
            route.clone(),
            original.clone(),
            position(0.0, 0.0, 0.0),
            profile,
            || Ok(()),
        );
        assert!(matches!(rejected, Err(MovingPartialSlipError::Invalid(_))));
    }
    let n = original.normal_patch();
    let foreign = NormalPatchView::new(
        n.patch_id(),
        n.normal_model_id(),
        "different normal source",
        n.authority(),
        n.normal_load_n(),
        n.semi_axis_longitudinal_m(),
        n.semi_axis_lateral_m(),
        n.pressure_second_moment_m2(),
    )
    .unwrap();
    let profile = MovingPressureProfile::new(&foreign, edge_core()).unwrap();
    assert!(matches!(
        MovingPartialSlip::new_with_pressure_profile(
            &thermal,
            route,
            original,
            position(0.0, 0.0, 0.0),
            profile,
            || Ok(()),
        ),
        Err(MovingPartialSlipError::Invalid(_))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
}

#[test]
fn pressure_resolution_and_late_source_refusal_preserve_retryable_history() {
    let (thermal, contact) = loaded(&[-3.0, 5.0], 1.0, initial([0.0; 2]), edge_core());
    let snapshot = contact.parcels().to_vec();
    let profile = contact.pressure_profile().cloned();
    let calls = Cell::new(0);
    let mut limited = controls();
    limited.maximum_regions = 2;
    let rejected = contact.step_candidate(
        &thermal,
        0.1,
        position(0.0, -0.004, 0.0),
        &work(),
        limited,
        &[0.0; 2],
        300.0,
        |t, r, _, _| {
            calls.set(calls.get() + 1);
            Ok::<_, Infallible>(source(t, r, 1.0))
        },
        || false,
    );
    assert!(matches!(
        rejected,
        Err(MovingPartialSlipError::RegionBudget)
    ));
    assert_eq!(calls.get(), 0);
    let rejected = contact.step_candidate(
        &thermal,
        0.1,
        position(0.0, -0.004, 0.0),
        &work(),
        controls(),
        &[0.0; 2],
        300.0,
        |t, r, _, _| {
            calls.set(calls.get() + 1);
            if calls.get() >= 3 {
                Err("late pressure-source refusal")
            } else {
                Ok(source(t, r, 1.0))
            }
        },
        || false,
    );
    assert!(matches!(
        rejected,
        Err(MovingPartialSlipError::Source(
            "late pressure-source refusal"
        ))
    ));
    assert_eq!(contact.parcels(), snapshot);
    assert_eq!(contact.pressure_profile(), profile.as_ref());
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    let retry = step(&contact, &thermal, 0.1, position(0.0, -0.004, 0.0), 1.0);
    near(retry.report().force_on_a_n[0], -0.34, 2e-13);
    near(retry.report().heat.dissipated_j, 0.00012, 2e-14);
}

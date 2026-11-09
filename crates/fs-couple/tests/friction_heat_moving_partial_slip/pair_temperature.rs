//! Independent local traction/heat oracles for one original A cell crossing
//! distinct B cells. A common A temperature cannot stand in for either pair's
//! actual material or massless interface temperature.

use super::*;

fn fixture(b_edges: &[f64]) -> (DepositionSimulation<'static>, MovingPartialSlip) {
    let grids = [&[-3.0, 5.0][..], b_edges];
    let count = b_edges.len();
    let mut events = Vec::new();
    let mut surfaces: [Vec<_>; 2] = std::array::from_fn(|_| Vec::new());
    for (side, grid) in grids.into_iter().enumerate() {
        for span in grid.windows(2) {
            let id = CellId(events.len() as u64);
            let length = span[1] - span[0];
            let centre = f64::midpoint(span[0], span[1]);
            let y = if side == 0 { [-1.0, 0.0] } else { [0.0, 1.0] };
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: side as u64,
                    material_id: MaterialId(0),
                    mass_kg: 100.0 * length,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [centre, y[0], 0.0],
                        end_m: [centre, y[1], 0.0],
                        volume_m3: length,
                        surface_area_m2: 4.0 * length + 2.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
            surfaces[side].push(StripSurfaceCell {
                cell: id,
                span_m: [span[0], span[1]],
            });
        }
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: count,
            max_links: 1,
            max_events: count,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let route = MovingFrictionRoute::new(
        &thermal,
        StripFrame {
            origin_m: [0.0; 3],
            tangent: [1.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            width_m: 1.0,
        },
        surfaces,
        [0.0, 1.0],
        MovingFrictionBudget {
            maximum_surface_cells: count,
            maximum_overlap_tests: count * count,
            maximum_patches: count * count,
        },
        || Ok(()),
    )
    .unwrap();
    let contact = MovingPartialSlip::new(
        &thermal,
        route,
        reference(1.0, PartialSlipState::zero()),
        position(0.0, 0.0, 0.0),
        || Ok(()),
    )
    .unwrap();
    (thermal, contact)
}

fn paired_source(
    thermal: &DepositionSimulation<'_>,
    regions: &[MovingSlipRegion],
    patches: &[SweptContactPatch],
) -> MovingPartialSlipSource<()> {
    let mut result = source(thermal, regions, 1.0);
    result.laws = regions
        .iter()
        .map(|r| {
            assert_eq!(
                patches[r.endpoint_interface_index].cells,
                [r.surface_a_cell, r.surface_b_cell],
            );
            assert_eq!(r.surface_a_cell, CellId(0));
            law(if r.surface_b_cell == CellId(1) {
                0.25
            } else {
                1.0
            })
        })
        .collect();
    result
}

fn paired_step(
    contact: &MovingPartialSlip,
    thermal: &DepositionSimulation<'static>,
    target_s: f64,
    offset: f64,
) -> MovingPartialSlipCandidate<'static, ()> {
    contact
        .step_candidate(
            thermal,
            target_s,
            position(offset, 0.0, offset),
            &work(),
            controls(),
            &vec![0.0; thermal.cells().len()],
            300.0,
            |t, regions, patches, _| Ok::<_, Infallible>(paired_source(t, regions, patches)),
            || false,
        )
        .unwrap()
}

#[test]
fn opposing_material_cells_resolve_distinct_traction_and_original_heat_residence() {
    let (thermal, contact) = fixture(&[-3.0, 0.6, 5.0]);
    let candidate = paired_step(&contact, &thermal, 0.02, 0.1);
    let report = candidate.report();
    assert_eq!(report.regions.len(), 2);
    for (index, r) in report.regions.iter().enumerate() {
        near(r.region.area_fraction, 0.5, 1e-14);
        assert_eq!(r.region.surface_b_cell, CellId(index as u64 + 1));
        assert_eq!(
            report.heat.patches[r.region.endpoint_interface_index].cells,
            [r.region.surface_a_cell, r.region.surface_b_cell],
        );
        let mu = if index == 0 { 0.25 } else { 1.0 };
        let strain = r.response.next_state.region_elastic_displacement_m();
        near(strain[0][0], 0.005 * mu, 1e-14);
        near(strain[1][0], 0.0175 * mu, 1e-14);
    }
    near(report.force_on_a_n[0], -0.625, 1e-13);
    near(report.mechanical_work_j, 0.0625, 1e-13);
    near(report.storage_after_j, 0.00365234375, 1e-13);
    near(report.numerical_dissipation_j, 0.00365234375, 1e-13);
    near(report.heat.dissipated_j, 0.0551953125, 1e-13);
    // A's right half spends 10% of this interval's area-time on the left B
    // cell before crossing its edge. Endpoint source selection must not move
    // that already-generated work to the right endpoint receiver.
    near(report.heat.patches[0].mean_area_m2, 0.55, 1e-14);
    near(report.heat.patches[1].mean_area_m2, 0.45, 1e-14);
    near(report.heat.patches[0].dissipated_j, 0.0163828125, 1e-13);
    near(report.heat.patches[1].dissipated_j, 0.0388125, 1e-13);
    near(
        report.heat.thermal.stored_energy_change_j,
        report.heat.dissipated_j,
        3e-10,
    );
    near(report.total_energy_residual_j, 0.0, 5e-10);
    assert_eq!(contact.parcels().len(), 1);
    assert_eq!(contact.parcels()[0].state(), &PartialSlipState::zero());
}

#[test]
fn interface_and_opposing_body_feedback_match_independent_condensed_heat_equations() {
    check_interface_feedback(None);
}

#[test]
fn imperfect_contact_resistance_changes_partial_slip_through_the_same_thermal_feedback() {
    check_interface_feedback(Some([[0.25, 0.5], [1.0, 0.125]]));
}

fn check_interface_feedback(resistance: Option<[[f64; 2]; 2]>) {
    let (thermal, contact) = fixture(&[-3.0, 0.6, 5.0]);
    let dt = 0.02;
    // Independent series resistance: original half-depth d=0.5 m and the
    // swept areas 0.55/0.45 m². Conductivities k_A=1 and k_B in {1,3} W/(m K)
    // remain separate from each interface's generation-plane resistance.
    let r = resistance.unwrap_or([[0.0; 2]; 2]);
    let ga = [0.55 / (0.5 + r[0][0]), 0.45 / (0.5 + r[1][0])];
    let gb = [0.55 / (0.5 + r[0][1]), 0.45 / (0.5 / 3.0 + r[1][1])];
    let cb = [360.0, 440.0];
    let bulk_weight = std::array::from_fn::<_, 2, _>(|i| dt * gb[i] / (cb[i] + dt * gb[i]));
    let a_weight = ga.map(|g| dt * g / (800.0 + dt * (ga[0] + ga[1])));
    let l: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            dt * (if i == j {
                ga[i] + gb[i] * (1.0 - bulk_weight[i])
            } else {
                0.0
            } - ga[i] * a_weight[j])
        })
    });
    let determinant = l[0][0] * l[1][1] - l[0][1] * l[1][0];
    let plastic_heat = |mu: f64| 0.1 * mu - 0.01375 * mu * mu;
    let temperatures = |mu: [f64; 2]| {
        let d = [
            0.5 * plastic_heat(mu[0]) + 0.05 * plastic_heat(mu[1]),
            0.45 * plastic_heat(mu[1]),
        ];
        [
            (l[1][1] * d[0] - l[0][1] * d[1]) / determinant,
            (l[0][0] * d[1] - l[1][0] * d[0]) / determinant,
        ]
    };
    // Eliminate the original three finite caloric temperatures analytically.
    // Two nested monotone scalar brackets then solve only the independently
    // authored constitutive equations, not the production Newton residual.
    let residual = |mu: [f64; 2], i: usize| {
        mu[i] * (1.0 + (1.0 + 1000.0 * bulk_weight[i]) * temperatures(mu)[i]) - 1.0
    };
    let left_mu = |right: f64| {
        let (mut lo, mut hi) = (0.0, 1.0);
        for _ in 0..90 {
            let mid = f64::midpoint(lo, hi);
            if residual([mid, right], 0) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        f64::midpoint(lo, hi)
    };
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..90 {
        let mid = f64::midpoint(lo, hi);
        if residual([left_mu(mid), mid], 1) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let right = f64::midpoint(lo, hi);
    let mu = [left_mu(right), right];
    let interface = temperatures(mu);
    let candidate = contact
        .step_candidate(
            &thermal,
            dt,
            position(0.1, 0.0, 0.1),
            &work(),
            controls(),
            &[0.0; 3],
            300.0,
            |t, regions, patches, temperatures| {
                let mut result = paired_source(t, regions, patches);
                result.conductivity_w_m_k.insert(CellId(1), 1.0);
                result.area_resistance_m2_k_w = resistance.map(|r| r.to_vec());
                result.laws = regions
                    .iter()
                    .map(|r| {
                        let interface = temperatures[r.endpoint_interface_index];
                        let body_b = t.temperature_k(r.surface_b_cell).unwrap();
                        law(1.0 / (1.0 + interface - 300.0 + 1000.0 * (body_b - 300.0)))
                    })
                    .collect();
                Ok::<_, Infallible>(result)
            },
            || false,
        )
        .unwrap();
    let report = candidate.report();
    near(report.force_on_a_n[0], -0.5 * (mu[0] + mu[1]), 2e-8);
    near(
        report.heat.dissipated_j,
        0.5 * (plastic_heat(mu[0]) + plastic_heat(mu[1])),
        3e-9,
    );
    for i in 0..2 {
        assert_eq!(report.heat.patches[i].area_resistance_m2_k_w, r[i]);
        near(
            report.heat.patches[i].interface_temperature_k,
            300.0 + interface[i],
            3e-8,
        );
        near(
            candidate
                .thermal()
                .temperature_k(CellId(i as u64 + 1))
                .unwrap(),
            300.0 + bulk_weight[i] * interface[i],
            3e-11,
        );
    }
    near(
        candidate.thermal().temperature_k(CellId(0)).unwrap(),
        300.0 + a_weight[0] * interface[0] + a_weight[1] * interface[1],
        3e-11,
    );
    assert!(
        (mu[0] - mu[1]).abs() > 0.05,
        "two pair temperatures must change traction differently"
    );
    assert!(
        mu[0] < 0.8 && mu[1] < 0.8,
        "finite heat feedback must be active"
    );
    near(report.total_energy_residual_j, 0.0, 5e-10);
}

#[test]
fn opposing_boundary_crossings_keep_material_strains_and_kinetic_history() {
    let (thermal, contact) = fixture(&[-3.0, 0.6, 5.0]);
    let (contact, thermal, first) = paired_step(&contact, &thermal, 0.02, 0.1).into_parts();
    let old = contact.parcels().to_vec();
    let next = paired_step(&contact, &thermal, 0.04, 0.2);
    let report = next.report();
    assert_eq!(report.regions.len(), 3);
    let expected = [[0.0, 0.4], [0.4, 0.5], [0.5, 1.0]];
    for (r, span) in report.regions.iter().zip(expected) {
        near(r.region.material_span_m[0], span[0], 1e-14);
        near(r.region.material_span_m[1], span[1], 1e-14);
    }
    assert_eq!(report.regions[0].region.accepted_state(), old[0].state());
    assert_eq!(report.regions[1].region.accepted_state(), old[0].state());
    assert_eq!(report.regions[2].region.accepted_state(), old[1].state());
    assert_eq!(report.regions[1].region.surface_b_cell, CellId(2));
    assert_ne!(
        report.regions[1].region.accepted_state(),
        report.regions[2].region.accepted_state()
    );
    near(report.storage_before_j, first.storage_after_j, 1e-14);
    near(report.force_on_a_n[0], -0.7, 1e-13);
    near(report.mechanical_work_j, 0.07, 1e-13);
    near(report.storage_after_j, 0.004296875, 1e-13);
    near(report.heat.dissipated_j, 0.06896875, 1e-13);
    near(report.numerical_dissipation_j, 0.00038671875, 1e-13);
    near(report.elastic_outflow_j, 0.0, 1e-14);
    near(report.total_energy_residual_j, 0.0, 5e-10);
    assert_eq!(contact.parcels(), old);
    assert!(
        next.contact()
            .parcels()
            .iter()
            .all(|p| p.state().accepted_steps() == 2)
    );
}

#[test]
fn constant_source_mechanics_are_invariant_to_opposing_thermal_grid_refinement() {
    let (coarse_heat, coarse) = fixture(&[-3.0, 5.0]);
    let (fine_heat, fine) = fixture(&[-3.0, 0.3, 0.6, 0.9, 5.0]);
    let end = position(0.1, 0.0, 0.1);
    let a = step(&coarse, &coarse_heat, 0.02, end, 1.0);
    let b = step(&fine, &fine_heat, 0.02, end, 1.0);
    assert_eq!(a.contact().parcels(), b.contact().parcels());
    assert_eq!(a.report().regions.len(), 1);
    assert_eq!(b.report().regions.len(), 4);
    for (x, y) in [
        (a.report().force_on_a_n[0], b.report().force_on_a_n[0]),
        (a.report().mechanical_work_j, b.report().mechanical_work_j),
        (a.report().storage_after_j, b.report().storage_after_j),
        (
            a.report().numerical_dissipation_j,
            b.report().numerical_dissipation_j,
        ),
        (a.report().heat.dissipated_j, b.report().heat.dissipated_j),
    ] {
        near(x, y, 1e-13);
    }
    near(
        a.report().heat.thermal.stored_energy_change_j,
        b.report().heat.thermal.stored_energy_change_j,
        5e-10,
    );
}

#[test]
fn opposing_cell_resolution_spends_the_original_region_budget_before_source_evaluation() {
    let (thermal, contact) = fixture(&[-3.0, 0.6, 5.0]);
    let parcels = contact.parcels().to_vec();
    let calls = Cell::new(0);
    let mut controls = controls();
    controls.maximum_regions = 1;
    let refused = contact.step_candidate(
        &thermal,
        0.02,
        position(0.1, 0.0, 0.1),
        &work(),
        controls,
        &[0.0; 3],
        300.0,
        |t, regions, patches, _| {
            calls.set(calls.get() + 1);
            Ok::<_, Infallible>(paired_source(t, regions, patches))
        },
        || false,
    );
    assert!(matches!(refused, Err(MovingPartialSlipError::RegionBudget)));
    assert_eq!(calls.get(), 0);
    assert_eq!(contact.parcels(), parcels);
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    let retried = paired_step(&contact, &thermal, 0.02, 0.1);
    near(retried.report().force_on_a_n[0], -0.625, 1e-13);
}

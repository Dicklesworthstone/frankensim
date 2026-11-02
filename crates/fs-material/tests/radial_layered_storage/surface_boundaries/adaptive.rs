//! Fixed physical faces carried through the original adaptive/event owner.
use super::*;
use fs_material::thermal::transport::stack::transient::implicit::LayeredThermalAdaptiveError;
use fs_thermal::enthalpy::{TimedEvent, conduction::implicit::ImplicitConductionError};

fn adaptive_controls(span_s: f64, temperature_tolerance_k: f64) -> LayeredThermalAdaptiveControls {
    let mut controls = LayeredThermalAdaptiveControls::default();
    controls.initial_dt_s = span_s;
    controls.step.maximum_dt_s = span_s;
    controls.temperature_absolute_tolerance_k = temperature_tolerance_k;
    controls.enthalpy_absolute_tolerance_j_kg = CP * temperature_tolerance_k;
    controls.relative_tolerance = 0.0;
    controls
}

fn scheduled<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    temperature_k: f64,
    boundary: Boundary,
    changes: &[(f64, Boundary)],
) -> DepositionSimulation<'a> {
    let mut events = bound
        .birth_events(0.0, &[temperature_k], conduction(), || Ok(()))
        .unwrap();
    let Event::Birth(cell) = &mut events[0].event else {
        panic!("original birth")
    };
    cell.boundary = boundary;
    events.extend(changes.iter().map(|&(time_s, boundary)| TimedEvent {
        time_s,
        event: Event::Boundary(bound.cell_ids()[0], boundary),
    }));
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: changes.len() + 1,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn g1_adaptive_physical_faces_refine_toward_the_independent_cooling_exponential() {
    let card = constant("adaptive finite material face", 0.4);
    let body = prepare(
        &[&card],
        &[0.04],
        &[1],
        &[],
        LayerStackGeometry::Planar { area_m2: 0.02 },
        false,
    )
    .unwrap();
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let old = initial(&bound, &[500.0], &[0.0], 80.0, 300.0);
    let before = old.cells().clone();
    let patches = [
        patch(LayeredFace::Inner, LayeredBoundaryChannel::Ambient, 0.02),
        patch(LayeredFace::Outer, LayeredBoundaryChannel::Ambient, 0.02),
    ];
    let capacity = RHO * CP * body.total_volume_m3();
    let conductance = 2.0 / (0.02 / (0.4 * 0.02) + 1.0 / (80.0 * 0.02));
    let exact = 300.0 + 200.0 * (-conductance * 5.0 / capacity).exp();
    let mut errors = Vec::new();
    let mut comparisons = Vec::new();
    for tolerance in [0.05, 0.001] {
        let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
        let (next, report) = bound
            .advance_adaptive_to_with_surface_boundaries(
                &old,
                5.0,
                &patches,
                conduction(),
                adaptive_controls(5.0, tolerance),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        let t = next.temperature_k(bound.cell_ids()[0]).unwrap();
        errors.push((t - exact).abs());
        comparisons.push(report.accepted_comparisons);
        close(
            report.solve.energy.external_heat_j,
            capacity * (t - 500.0),
            2e-8,
        );
        close(
            next.energy_j() - old.energy_j(),
            report.solve.energy.external_heat_j,
            2e-8,
        );
        let entropy = report.solve.entropy;
        close(
            entropy.thermal_storage_change_j_k,
            capacity * (t / 500.0).ln(),
            2e-10,
        );
        close(
            entropy.reservoir_entropy_input_j_k,
            report.solve.energy.external_heat_j / 300.0,
            2e-10,
        );
        assert!(!entropy.has_unclassified_sources);
        assert_eq!(entropy.unclassified_source_at_cell_j_k, 0.0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        if tolerance == 0.001 {
            assert!(report.accuracy_rejections > 0);
        }
    }
    assert!(
        errors[1] < 0.5 * errors[0] && errors[1] < 0.01,
        "{errors:?}"
    );
    assert!(comparisons[1] > comparisons[0]);
    // The existing API and an explicit empty map keep precisely the same path.
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let (legacy, legacy_report) = bound
        .advance_adaptive_to(
            &old,
            5.0,
            conduction(),
            adaptive_controls(5.0, 0.05),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let (empty, empty_report) = bound
        .advance_adaptive_to_with_surface_boundaries(
            &old,
            5.0,
            &[],
            conduction(),
            adaptive_controls(5.0, 0.05),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(legacy.cells(), empty.cells());
    assert_eq!(legacy_report, empty_report);
    assert!(exact - legacy.temperature_k(bound.cell_ids()[0]).unwrap() > 10.0);
    assert_eq!(old.cells(), &before);
}

#[test]
fn g1_adaptive_events_refresh_baths_films_and_charge_only_accepted_half_steps() {
    let card = constant("scheduled physical faces", 0.4);
    let body = prepare(
        &[&card],
        &[0.04],
        &[1],
        &[],
        LayerStackGeometry::Planar { area_m2: 0.02 },
        false,
    )
    .unwrap();
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial_boundary = Boundary {
        ambient_k: 250.0,
        htc_w_m2_k: 80.0,
        source_w: 5.0,
        reservoir: Some(Reservoir {
            area_m2: 0.02,
            conductance_w_k: 2.0,
            temperature_k: 700.0,
        }),
    };
    let changed = Boundary {
        ambient_k: 350.0,
        htc_w_m2_k: 20.0,
        source_w: -2.0,
        reservoir: Some(Reservoir {
            area_m2: 0.02,
            conductance_w_k: 0.3,
            temperature_k: 600.0,
        }),
    };
    // The target event changes the NEXT interval; it supplies no retroactive heat.
    let at_target = Boundary {
        ambient_k: 1000.0,
        htc_w_m2_k: 0.0,
        source_w: 0.0,
        reservoir: Some(Reservoir {
            area_m2: 0.02,
            conductance_w_k: 0.5,
            temperature_k: 300.0,
        }),
    };
    let old = scheduled(
        &bound,
        450.0,
        initial_boundary,
        &[(2.0, changed), (5.0, at_target)],
    );
    let before = old.cells().clone();
    let patches = [
        patch(LayeredFace::Inner, LayeredBoundaryChannel::Reservoir, 0.02),
        patch(LayeredFace::Outer, LayeredBoundaryChannel::Ambient, 0.02),
    ];
    let capacity = RHO * CP * body.total_volume_m3();
    let material_resistance = 0.02 / (0.4 * 0.02);
    let mut expected = 450.0;
    let mut energy = 0.0;
    let mut bath_entropy = 0.0;
    let mut source_entropy = 0.0;
    let mut generation = 0.0;
    // With loose time-error limits, the native event-clipped comparisons accept
    // exactly these four backward-Euler half steps. The coarse work is excluded.
    for (dt, boundary) in [
        (1.0, initial_boundary),
        (1.0, initial_boundary),
        (1.5, changed),
        (1.5, changed),
    ] {
        let reservoir = boundary.reservoir.unwrap();
        let ga = 1.0 / (material_resistance + 1.0 / (boundary.htc_w_m2_k * 0.02));
        let gr = 1.0 / (material_resistance + 1.0 / reservoir.conductance_w_k);
        expected = (capacity * expected
            + dt * (boundary.source_w + ga * boundary.ambient_k + gr * reservoir.temperature_k))
            / (capacity + dt * (ga + gr));
        let qa = ga * (boundary.ambient_k - expected);
        let qr = gr * (reservoir.temperature_k - expected);
        energy += dt * (boundary.source_w + qa + qr);
        bath_entropy += dt * (qa / boundary.ambient_k + qr / reservoir.temperature_k);
        source_entropy += dt * boundary.source_w / expected;
        generation += dt
            * (qa * (1.0 / expected - 1.0 / boundary.ambient_k)
                + qr * (1.0 / expected - 1.0 / reservoir.temperature_k));
    }
    let controls = adaptive_controls(5.0, 1000.0);
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let mut calls = 0;
    let (next, report) = bound
        .advance_adaptive_to_with_surface_boundaries(
            &old,
            5.0,
            &patches,
            conduction(),
            controls,
            &mut workspace,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    close(
        next.temperature_k(bound.cell_ids()[0]).unwrap(),
        expected,
        2e-9,
    );
    close(report.solve.energy.external_heat_j, energy, 2e-8);
    close(
        report.solve.energy.stored_energy_change_j,
        capacity * (expected - 450.0),
        2e-8,
    );
    assert_eq!(report.accepted_comparisons, 2);
    assert_eq!(report.solve.energy.substeps, 4);
    assert_eq!(report.solve.energy.events, 2);
    assert_eq!(report.accuracy_rejections, 0);
    assert_eq!(next.cells()[&bound.cell_ids()[0]].boundary, at_target);
    let entropy = report.solve.entropy;
    close(entropy.reservoir_entropy_input_j_k, bath_entropy, 2e-10);
    close(
        entropy.unclassified_source_at_cell_j_k,
        source_entropy,
        2e-12,
    );
    close(entropy.conduction_generation_j_k, generation, 2e-10);
    close(
        entropy.thermal_storage_change_j_k,
        capacity * (expected / 450.0).ln(),
        2e-10,
    );
    assert!(entropy.has_unclassified_sources);
    let (continued, _) = bound
        .advance_adaptive_to_with_surface_boundaries(
            &next,
            5.5,
            &patches,
            conduction(),
            adaptive_controls(0.5, 1000.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let gr = 1.0 / (material_resistance + 1.0 / 0.5);
    let expected_continued = 300.0 + (expected - 300.0) / (1.0 + 0.25 * gr / capacity).powi(2);
    close(
        continued.temperature_k(bound.cell_ids()[0]).unwrap(),
        expected_continued,
        2e-9,
    );

    // Cancellation at the final checkpoint discards all intermediate accepted
    // candidates and consumed native events, rather than publishing a prefix.
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let mut late = 0;
    assert!(matches!(
        bound.advance_adaptive_to_with_surface_boundaries(
            &old,
            5.0,
            &patches,
            conduction(),
            controls,
            &mut workspace,
            || {
                late += 1;
                if late == calls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(LayeredThermalAdaptiveError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(old.cells(), &before);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn g0_target_events_cannot_publish_incompatible_coverage_or_unsupported_material_faces() {
    let card = constant("adaptive terminal event admission", 0.4);
    let body = prepare(
        &[&card],
        &[0.04],
        &[1],
        &[],
        LayerStackGeometry::Planar { area_m2: 0.02 },
        false,
    )
    .unwrap();
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let boundary = Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 20.0,
        source_w: 0.0,
        reservoir: Some(Reservoir {
            area_m2: 0.02,
            conductance_w_k: 0.5,
            temperature_k: 350.0,
        }),
    };
    let patches = [
        patch(LayeredFace::Inner, LayeredBoundaryChannel::Reservoir, 0.02),
        patch(LayeredFace::Outer, LayeredBoundaryChannel::Ambient, 0.02),
    ];
    let incompatible = [
        Boundary {
            reservoir: None,
            ..boundary
        },
        Boundary {
            reservoir: Some(Reservoir {
                area_m2: 0.01,
                ..boundary.reservoir.unwrap()
            }),
            ..boundary
        },
        Boundary {
            ambient_k: 1000.0,
            htc_w_m2_k: 1e6,
            ..boundary
        },
    ];
    for time in [2.0, 5.0] {
        for (i, bad) in incompatible.iter().copied().enumerate() {
            let old = scheduled(&bound, 450.0, boundary, &[(time, bad)]);
            let before = old.cells().clone();
            let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
            let result = bound.advance_adaptive_to_with_surface_boundaries(
                &old,
                5.0,
                &patches,
                conduction(),
                adaptive_controls(5.0, 1000.0),
                &mut workspace,
                || Ok(()),
            );
            if i == 2 {
                assert!(matches!(
                    result,
                    Err(LayeredThermalAdaptiveError::Implicit(
                        ImplicitConductionError::Thermal(DepositionError::MaterialDomain)
                    ))
                ));
            } else {
                assert!(result.is_err());
            }
            assert_eq!(old.cells(), &before);
            assert_eq!(old.time_s(), 0.0);
        }
    }
    let old = scheduled(&bound, 450.0, boundary, &[]);
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    assert!(
        bound
            .advance_adaptive_to_with_surface_boundaries(
                &old,
                0.0,
                &[LayeredSurfacePatch {
                    area_m2: 0.01,
                    ..patches[0]
                }],
                conduction(),
                adaptive_controls(1.0, 1.0),
                &mut workspace,
                || Ok(()),
            )
            .is_err()
    );
}

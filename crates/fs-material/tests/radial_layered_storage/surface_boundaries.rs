//! Physical face/film integration on the original material-backed time owner.
use super::*;
use fs_material::thermal::transport::stack::transient::{
    LayeredBoundaryChannel, LayeredFace, LayeredSurfacePatch,
};
use fs_thermal::enthalpy::{Boundary, Reservoir};

fn patch(face: LayeredFace, channel: LayeredBoundaryChannel, area_m2: f64) -> LayeredSurfacePatch {
    LayeredSurfacePatch {
        face,
        channel,
        area_m2,
    }
}

// The existing parent fixture supplies sources, geometry and caloric laws.
// This one-cell variant additionally authors an original reservoir boundary.
fn reservoir_initial<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    temperature: f64,
    boundary: Boundary,
) -> DepositionSimulation<'a> {
    let mut events = bound
        .birth_events(0.0, &[temperature], conduction(), || Ok(()))
        .unwrap();
    let Event::Birth(cell) = &mut events[0].event else {
        panic!("original cell birth")
    };
    cell.boundary = boundary;
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1,
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
fn g1_face_series_resistance_uses_actual_planar_cylindrical_and_spherical_areas() {
    let card = constant("physical terminal faces", 8.0);
    for geometry in [
        LayerStackGeometry::Planar { area_m2: 0.02 },
        shape(true, 0.01),
        shape(false, 0.01),
    ] {
        let body = prepare(&[&card], &[0.03], &[1], &[], geometry, false).unwrap();
        let laws = body.thermal_laws(true).unwrap();
        let bound = body.bind(&laws, ids()).unwrap();
        let cell = body.cell_geometries().next().unwrap();
        let center = body.cell_temperature_positions_m().next().unwrap();
        let radii = [cell.start_m[0], cell.end_m[0]];
        let areas = radii.map(|r| area(geometry, r));
        let reservoir_area = 0.25 * (areas[0] + areas[1]);
        let source = 7.0;
        let boundary = Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 80.0,
            source_w: source,
            reservoir: Some(Reservoir {
                area_m2: reservoir_area,
                conductance_w_k: 40.0 * reservoir_area,
                temperature_k: 700.0,
            }),
        };
        let old = reservoir_initial(&bound, 500.0, boundary);
        let old_cells = old.cells().clone();
        let patches = [
            patch(
                LayeredFace::Inner,
                LayeredBoundaryChannel::Ambient,
                0.75 * areas[0],
            ),
            patch(
                LayeredFace::Inner,
                LayeredBoundaryChannel::Reservoir,
                0.25 * areas[0],
            ),
            patch(
                LayeredFace::Outer,
                LayeredBoundaryChannel::Ambient,
                0.75 * areas[1],
            ),
            patch(
                LayeredFace::Outer,
                LayeredBoundaryChannel::Reservoir,
                0.25 * areas[1],
            ),
        ];
        let dt = 0.02;
        let step = LayeredThermalStep::new(&bound, &old, dt, conduction(), || Ok(()))
            .unwrap()
            .with_surface_boundaries(&patches, || Ok(()))
            .unwrap();
        let trial = step.evaluate(&[0.0], || Ok(())).unwrap();
        assert_eq!(trial.payload().surface_responses().len(), 4);
        let mut channel_g = [0.0; 2];
        let mut total_outward = 0.0;
        for surface in trial.payload().surface_responses() {
            let face = if surface.patch.face == LayeredFace::Inner {
                0
            } else {
                1
            };
            let (channel, fraction, htc, bath) = match surface.patch.channel {
                LayeredBoundaryChannel::Ambient => (0, 0.75, 80.0, 300.0),
                LayeredBoundaryChannel::Reservoir => (1, 0.25, 40.0, 700.0),
            };
            let film_g = htc * fraction * areas[face];
            // Independent textbook resistance integral; a partial patch is
            // an explicitly parallel normal streamtube, not lateral spreading.
            let material_r = metric(geometry, center.min(radii[face]), center.max(radii[face]))
                / (8.0 * fraction);
            let expected_g = 1.0 / (material_r + 1.0 / film_g);
            let expected_q = expected_g * (500.0 - bath);
            let response = surface.response;
            close(response.conductance_w_k, expected_g, 2e-10);
            close(response.power_to_reservoir_w, expected_q, 2e-8);
            close(response.endpoint_derivative_w_k[0], expected_g, 2e-10);
            close(response.endpoint_derivative_w_k[1], -expected_g, 2e-10);
            close(
                response.material_face_temperature_k,
                bath + expected_q / film_g,
                2e-8,
            );
            assert_eq!(surface.patch.area_m2, fraction * areas[face]);
            assert_eq!(surface.cell_id, bound.cell_ids()[0]);
            channel_g[channel] += expected_g;
            total_outward += expected_q;
        }
        close(trial.residual_j()[0], dt * (total_outward - source), 2e-8);
        close(
            trial.audit().external_heat_j,
            dt * (source - total_outward),
            2e-8,
        );
        let capacity = RHO * CP * body.total_volume_m3();
        close(
            trial.payload().diagonal()[0],
            1.0 + dt * (channel_g[0] + channel_g[1]) / capacity,
            2e-11,
        );

        // Both original bath accounts and the authored source occur exactly once.
        let expected_t = (capacity * 500.0
            + dt * (source + channel_g[0] * 300.0 + channel_g[1] * 700.0))
            / (capacity + dt * (channel_g[0] + channel_g[1]));
        let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
        let (next, _, report) = step
            .solve_with_entropy(
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        close(
            next.temperature_k(bound.cell_ids()[0]).unwrap(),
            expected_t,
            2e-9,
        );
        close(
            report.energy.stored_energy_change_j,
            capacity * (expected_t - 500.0),
            2e-8,
        );
        close(
            report.energy.external_heat_j,
            dt * (source
                + channel_g[0] * (300.0 - expected_t)
                + channel_g[1] * (700.0 - expected_t)),
            2e-8,
        );
        assert!(report.energy.residual_j.abs() <= 1e-8);
        let entropy = report.entropy.unwrap();
        close(
            entropy.thermal_storage_change_j_k,
            capacity * (expected_t / 500.0).ln(),
            2e-10,
        );
        close(
            entropy.reservoir_entropy_input_j_k,
            dt * (channel_g[0] * (300.0 - expected_t) / 300.0
                + channel_g[1] * (700.0 - expected_t) / 700.0),
            2e-10,
        );
        close(
            entropy.conduction_generation_j_k,
            dt * (channel_g[0] * (expected_t - 300.0).powi(2) / (expected_t * 300.0)
                + channel_g[1] * (expected_t - 700.0).powi(2) / (expected_t * 700.0)),
            2e-10,
        );
        close(
            entropy.unclassified_source_at_cell_j_k,
            dt * source / expected_t,
            2e-12,
        );
        assert!(entropy.has_unclassified_sources);
        close(
            entropy.temporal_defect_j_k,
            entropy.thermal_storage_change_j_k
                - entropy.reservoir_entropy_input_j_k
                - entropy.conduction_generation_j_k
                - entropy.unclassified_source_at_cell_j_k,
            2e-12,
        );

        let original = LayeredThermalStep::new(&bound, &old, dt, conduction(), || Ok(()))
            .unwrap()
            .evaluate(&[0.0], || Ok(()))
            .unwrap();
        let empty = LayeredThermalStep::new(&bound, &old, dt, conduction(), || Ok(()))
            .unwrap()
            .with_surface_boundaries(&[], || Ok(()))
            .unwrap()
            .evaluate(&[0.0], || Ok(()))
            .unwrap();
        assert_eq!(original.residual_j(), empty.residual_j());
        assert_eq!(original.payload().diagonal(), empty.payload().diagonal());
        assert!(empty.payload().surface_responses().is_empty());
        assert_eq!(old.cells(), &old_cells);
    }
}

#[test]
fn g1_nonlinear_face_tangent_and_native_energy_entropy_use_the_same_heat() {
    let card = material(
        "affine face conductivity",
        [(200.0, 2.0), (450.0, 4.5), (800.0, 8.0)],
    );
    let geometry = shape(true, 0.01);
    let body = prepare(&[&card], &[0.03], &[3], &[], geometry, false).unwrap();
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let sources = [2.0, 0.0, -1.0];
    let old = initial(&bound, &[600.0, 480.0, 350.0], &sources, 50.0, 300.0);
    let old_cells = old.cells().clone();
    let ids = bound.cell_ids();
    let cells: Vec<_> = body.cell_geometries().collect();
    let centers: Vec<_> = body.cell_temperature_positions_m().collect();
    let patches = [
        patch(
            LayeredFace::Inner,
            LayeredBoundaryChannel::Ambient,
            area(geometry, cells[0].start_m[0]),
        ),
        patch(
            LayeredFace::Outer,
            LayeredBoundaryChannel::Ambient,
            area(geometry, cells[2].end_m[0]),
        ),
    ];
    let dt = 0.02;
    let step = LayeredThermalStep::new(&bound, &old, dt, conduction(), || Ok(()))
        .unwrap()
        .with_surface_boundaries(&patches, || Ok(()))
        .unwrap();
    let energy = [0.01, -0.02, 0.03];
    let direction = [0.3, -0.2, 0.4];
    let trial = step.evaluate(&energy, || Ok(())).unwrap();
    for surface in trial.payload().surface_responses() {
        let i = if surface.patch.face == LayeredFace::Inner {
            0
        } else {
            2
        };
        assert_eq!(surface.cell_id, ids[i]);
        let radius = if i == 0 {
            cells[i].start_m[0]
        } else {
            cells[i].end_m[0]
        };
        let multiplier = 1.0 / metric(geometry, radius.min(centers[i]), radius.max(centers[i]));
        let film = 1.0 / (50.0 * surface.patch.area_m2);
        let cell_t = surface.cell_temperature_k;
        // k(T)=0.01*T gives a quadratic equation for u=T_face-T_bath.
        // This closed form is independent of the material/film root iteration.
        let a = 0.005 * multiplier * film;
        let b = 1.0 + 0.01 * multiplier * film * 300.0;
        let d = a * (cell_t * cell_t - 300.0_f64.powi(2));
        let offset = 2.0 * d / (b + (b * b + 4.0 * a * d).sqrt());
        let face_t = 300.0 + offset;
        let q = offset / film;
        let tangent = multiplier * 0.01 * cell_t / (1.0 + film * multiplier * 0.01 * face_t);
        close(surface.response.material_face_temperature_k, face_t, 2e-8);
        close(surface.response.power_to_reservoir_w, q, 2e-8);
        close(surface.response.endpoint_derivative_w_k[0], tangent, 2e-9);
    }
    let mut action = [0.0; 3];
    trial
        .payload()
        .apply_into(&direction, &mut action, || Ok(()))
        .unwrap();
    let epsilon = 1e-4;
    let plus = step
        .evaluate(
            &std::array::from_fn::<_, 3, _>(|i| energy[i] + epsilon * direction[i]),
            || Ok(()),
        )
        .unwrap();
    let minus = step
        .evaluate(
            &std::array::from_fn::<_, 3, _>(|i| energy[i] - epsilon * direction[i]),
            || Ok(()),
        )
        .unwrap();
    for i in 0..3 {
        close(
            action[i],
            (plus.residual_j()[i] - minus.residual_j()[i]) / (2.0 * epsilon),
            2e-6 * action[i].abs().max(1.0),
        );
    }

    let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
    let (next, linearization, report) = step
        .solve_with_entropy(
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    let t = linearization.cell_temperatures_k();
    let surface_heat: f64 = linearization
        .surface_responses()
        .iter()
        .map(|surface| surface.response.power_to_reservoir_w)
        .sum();
    // Any ulp-sized uncovered middle-cell remainder keeps its original channel.
    let middle_heat = old.exposed_area_m2(ids[1]).unwrap() * 50.0 * (t[1] - 300.0);
    let external = dt * (sources.iter().sum::<f64>() - surface_heat - middle_heat);
    close(report.energy.external_heat_j, external, 2e-8);
    close(next.energy_j() - old.energy_j(), external, 2e-8);
    assert!(report.energy.residual_j.abs() <= 1e-8);
    let entropy = report.entropy.unwrap();
    close(
        entropy.reservoir_entropy_input_j_k,
        -dt * (surface_heat + middle_heat) / 300.0,
        2e-10,
    );
    let source_entropy: f64 = sources.iter().zip(t).map(|(q, t)| q / t).sum();
    close(
        entropy.unclassified_source_at_cell_j_k,
        dt * source_entropy,
        2e-12,
    );
    let cell_capacity = RHO * CP * body.total_volume_m3() / 3.0;
    let storage_entropy: f64 = t
        .iter()
        .zip([600.0, 480.0, 350.0])
        .map(|(new, old)| cell_capacity * (new / old).ln())
        .sum();
    close(entropy.thermal_storage_change_j_k, storage_entropy, 2e-10);
    let mut generation = middle_heat * (1.0 / 300.0 - 1.0 / t[1]);
    for surface in linearization.surface_responses() {
        generation += surface.response.power_to_reservoir_w
            * (1.0 / surface.bath_temperature_k - 1.0 / surface.cell_temperature_k);
    }
    for (i, link) in linearization.link_responses().iter().enumerate() {
        generation += link.power_a_to_b_w * (1.0 / t[i + 1] - 1.0 / t[i]);
    }
    close(entropy.conduction_generation_j_k, dt * generation, 2e-10);
    assert!(entropy.conduction_generation_j_k >= 0.0);
    assert_eq!(old.cells(), &old_cells);
}

#[test]
fn g1_outside_baths_require_supported_material_faces_and_zero_g_is_insulation() {
    let card = constant("supported material face", 8.0);
    let geometry = LayerStackGeometry::Planar { area_m2: 0.02 };
    let body = prepare(&[&card], &[0.04], &[1], &[], geometry, false).unwrap();
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let patches = [
        patch(LayeredFace::Inner, LayeredBoundaryChannel::Ambient, 0.02),
        patch(LayeredFace::Outer, LayeredBoundaryChannel::Ambient, 0.02),
    ];
    for bath in [100.0, 1000.0] {
        let old = initial(&bound, &[300.0], &[0.0], 4.0, bath);
        let trial = LayeredThermalStep::new(&bound, &old, 0.1, conduction(), || Ok(()))
            .unwrap()
            .with_surface_boundaries(&patches, || Ok(()))
            .unwrap()
            .evaluate(&[0.0], || Ok(()))
            .unwrap();
        let expected_q = (300.0 - bath) / (0.02 / (8.0 * 0.02) + 1.0 / (4.0 * 0.02));
        for surface in trial.payload().surface_responses() {
            close(surface.response.power_to_reservoir_w, expected_q, 2e-9);
            close(
                surface.response.material_face_temperature_k,
                bath + expected_q / 0.08,
                2e-8,
            );
            assert!((200.0..=800.0).contains(&surface.response.material_face_temperature_k));
        }
        let unsupported = initial(&bound, &[300.0], &[0.0], 1e6, bath);
        let unsupported_step =
            LayeredThermalStep::new(&bound, &unsupported, 0.1, conduction(), || Ok(()))
                .unwrap()
                .with_surface_boundaries(&patches, || Ok(()))
                .unwrap();
        assert!(matches!(
            unsupported_step.evaluate(&[0.0], || Ok(())),
            Err(DepositionError::MaterialDomain)
        ));

        let insulated = initial(&bound, &[300.0], &[0.0], 0.0, bath);
        let trial = LayeredThermalStep::new(&bound, &insulated, 0.1, conduction(), || Ok(()))
            .unwrap()
            .with_surface_boundaries(&patches, || Ok(()))
            .unwrap()
            .evaluate(&[0.0], || Ok(()))
            .unwrap();
        assert_eq!(trial.residual_j(), &[0.0]);
        assert_eq!(trial.payload().diagonal(), &[1.0]);
        for surface in trial.payload().surface_responses() {
            assert_eq!(surface.response.power_to_reservoir_w, 0.0);
            assert_eq!(surface.response.endpoint_derivative_w_k, [0.0, 0.0]);
            assert_eq!(surface.response.conductance_w_k, 0.0);
            assert_eq!(surface.response.material_face_temperature_k, 300.0);
        }
    }
}

#[test]
fn g0_surface_mapping_refusals_and_late_cancellation_leave_original_state_untouched() {
    let card = constant("surface admission", 8.0);
    let body = prepare(
        &[&card],
        &[0.04],
        &[1],
        &[],
        LayerStackGeometry::Planar { area_m2: 0.02 },
        false,
    )
    .unwrap();
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let old = initial(&bound, &[400.0], &[0.0], 50.0, 300.0);
    let before = old.cells().clone();
    let inner = patch(LayeredFace::Inner, LayeredBoundaryChannel::Ambient, 0.02);
    let outer = patch(LayeredFace::Outer, LayeredBoundaryChannel::Ambient, 0.02);
    let make = || LayeredThermalStep::new(&bound, &old, 0.1, conduction(), || Ok(())).unwrap();
    let invalid = [
        vec![inner],        // Incomplete original channel coverage.
        vec![inner, inner], // Duplicate face/channel.
        vec![
            LayeredSurfacePatch {
                area_m2: 0.021,
                ..inner
            },
            LayeredSurfacePatch {
                area_m2: 0.019,
                ..outer
            },
        ], // A face is overfilled.
        vec![patch(
            LayeredFace::Inner,
            LayeredBoundaryChannel::Reservoir,
            0.02,
        )],
        vec![inner; 5], // Explicit bounded map.
    ];
    for patches in invalid {
        assert!(make().with_surface_boundaries(&patches, || Ok(())).is_err());
    }
    for bad_area in [0.0, -0.01, f64::INFINITY, f64::NAN] {
        assert!(
            make()
                .with_surface_boundaries(
                    &[
                        LayeredSurfacePatch {
                            area_m2: bad_area,
                            ..inner
                        },
                        outer
                    ],
                    || Ok(()),
                )
                .is_err()
        );
    }
    assert!(matches!(
        make().with_surface_boundaries(&[inner, outer], || Err(DepositionError::Cancelled),),
        Err(DepositionError::Cancelled)
    ));
    let step = make()
        .with_surface_boundaries(&[inner, outer], || Ok(()))
        .unwrap();
    let mut calls = 0;
    let reference = step
        .evaluate(&[0.0], || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut late = 0;
    assert!(matches!(
        step.evaluate(&[0.0], || {
            late += 1;
            if late == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    let retry = step.evaluate(&[0.0], || Ok(())).unwrap();
    assert_eq!(reference.residual_j(), retry.residual_j());
    assert_eq!(
        reference.payload().surface_responses(),
        retry.payload().surface_responses()
    );
    assert_eq!(old.cells(), &before);
    assert_eq!(old.time_s(), 0.0);
}

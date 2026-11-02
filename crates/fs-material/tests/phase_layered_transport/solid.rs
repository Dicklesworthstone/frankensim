//! Phase-changing solid centres use regular radial symmetry and native heat.
use super::*;

fn solid_shapes() -> [LayerStackGeometry; 2] {
    [
        LayerStackGeometry::Cylindrical {
            inner_radius_m: 0.0,
            length_m: 0.7,
        },
        LayerStackGeometry::Spherical {
            inner_radius_m: 0.0,
        },
    ]
}

fn solid<'a>(
    sources: &[&'a PreparedPhaseConductivity<'a>],
    widths: &[f64],
    counts: &[usize],
    contacts: &[f64],
    shape: LayerStackGeometry,
) -> PreparedPhaseLayeredThermalBody<'a> {
    let mut radius = 0.0;
    let layers = sources
        .iter()
        .enumerate()
        .map(|(i, &conductivity)| {
            let mass_kg = volume(shape, radius, radius + widths[i]);
            radius += widths[i];
            PhaseLayerSource {
                id: 100 + i as u64,
                thickness_m: widths[i],
                mass_kg,
                cell_count: counts[i],
                conductivity,
            }
        })
        .collect::<Vec<_>>();
    match shape {
        LayerStackGeometry::Cylindrical { length_m, .. } => {
            PreparedPhaseLayeredThermalBody::try_new_solid_cylindrical(
                &layers,
                contacts,
                length_m,
                || Ok(()),
            )
            .unwrap()
        }
        LayerStackGeometry::Spherical { .. } => {
            PreparedPhaseLayeredThermalBody::try_new_solid_spherical(&layers, contacts, || Ok(()))
                .unwrap()
        }
        LayerStackGeometry::Planar { .. } => unreachable!(),
    }
}

#[test]
fn solid_centres_retain_physical_inventory_and_refuse_an_inner_bath() {
    let a = material(1.0, 4.0);
    let b = material(2.0, 6.0);
    let ca = caloric(&a, 400.0, 100.0);
    let cb = caloric(&b, 400.0, 150.0);
    let ka = prepare(&a, &ca);
    let kb = prepare(&b, &cb);
    for shape in solid_shapes() {
        let body = solid(&[&ka, &kb], &[0.4, 0.6], &[2, 3], &[0.1], shape);
        assert!(body.has_solid_center());
        assert_eq!(body.geometry(), shape);
        close(body.total_volume_m3(), volume(shape, 0.0, 1.0), 1e-13);
        close(body.total_mass_kg(), body.total_volume_m3(), 1e-13);
        let cells = body.cell_geometries().collect::<Vec<_>>();
        let positions = body.cell_temperature_positions_m().collect::<Vec<_>>();
        assert_eq!(cells[0].start_m[0], 0.0);
        close(positions[0], centroid(shape, 0.0, cells[0].end_m[0]), 1e-15);
        close(
            cells[0].surface_area_m2,
            area(shape, cells[0].end_m[0]),
            1e-15,
        );
        for (cell, &position) in cells.iter().zip(&positions) {
            assert!(cell.start_m[0] < position && position < cell.end_m[0]);
            close(
                cell.volume_m3,
                volume(shape, cell.start_m[0], cell.end_m[0]),
                1e-13,
            );
        }
        // Zero radial origin is accepted only by the explicit symmetry API.
        assert!(
            PreparedPhaseLayeredThermalBody::try_new(body.layers(), &[0.1], shape, || Ok(()))
                .is_err()
        );
        assert!(body.bind(ids(), [Some(0.0), None]).is_err());
        let bound = body.bind(ids(), [None, Some(0.1)]).unwrap();
        let before = initial(&bound, &[20.0; 5], [None, Some(500.0)]);
        for (&id, cell) in before.cells() {
            let exposed = before.exposed_area_m2(id).unwrap();
            assert!(exposed >= 0.0);
            assert!(exposed <= 8.0 * f64::EPSILON * cell.geometry.surface_area_m2);
            assert_eq!(cell.boundary.htc_w_m2_k, 0.0);
        }
        assert!(
            bound
                .birth_events(0.0, &[20.0; 5], [Some(300.0), Some(500.0)], || Ok(()))
                .is_err()
        );
        match shape {
            LayerStackGeometry::Cylindrical { length_m, .. } => {
                for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                    assert!(
                        PreparedPhaseLayeredThermalBody::try_new_solid_cylindrical(
                            body.layers(),
                            &[0.1],
                            invalid,
                            || Ok(())
                        )
                        .is_err()
                    );
                }
                assert!(matches!(
                    PreparedPhaseLayeredThermalBody::try_new_solid_cylindrical(
                        body.layers(),
                        &[0.1],
                        length_m,
                        || Err(DepositionError::Cancelled)
                    ),
                    Err(DepositionError::Cancelled)
                ));
            }
            LayerStackGeometry::Spherical { .. } => {
                assert!(matches!(
                    PreparedPhaseLayeredThermalBody::try_new_solid_spherical(
                        body.layers(),
                        &[0.1],
                        || Err(DepositionError::Cancelled)
                    ),
                    Err(DepositionError::Cancelled)
                ));
            }
            LayerStackGeometry::Planar { .. } => unreachable!(),
        }
    }
}

#[test]
fn solid_core_latent_heating_uses_original_conductivity_and_outer_face_resistance() {
    let card = material(1.0, 4.0);
    let curve = caloric(&card, 400.0, 100.0);
    let law = prepare(&card, &curve);
    for shape in solid_shapes() {
        let body = solid(&[&law], &[1.0], &[1], &[], shape);
        let bound = body.bind(ids(), [None, Some(0.2)]).unwrap();
        let before = initial(&bound, &[20.0], [None, Some(500.0)]);
        let original = before.cells().clone();
        let id = bound.cell_ids()[0];
        let mass = volume(shape, 0.0, 1.0);
        let centre = centroid(shape, 0.0, 1.0);
        let bulk = metric(shape, centre, 1.0);
        let film = 0.2 / area(shape, 1.0);
        let dt = 0.01;
        let power = |heat: f64| {
            let k = 1.0 + 3.0 * (20.0 + heat / mass) / 100.0;
            100.0 / (bulk / k + film)
        };
        let expected = root(0.0, 80.0 * mass, |heat| heat - dt * power(heat));
        let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = chain
            .advance_to(&before, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        close(
            next.cells()[&id].specific_h_j_kg,
            20.0 + expected / mass,
            2e-8,
        );
        close(next.temperature_k(id).unwrap(), 400.0, 1e-12);
        close(report.energy.stored_energy_change_j, expected, 1e-8);
        close(report.energy.external_heat_j, expected, 1e-8);
        close(report.energy.residual_j, 0.0, 1e-8);
        close(
            report.entropy.thermal_storage_change_j_k,
            expected / 400.0,
            1e-10,
        );
        close(
            report.entropy.reservoir_entropy_input_j_k,
            expected / 500.0,
            1e-10,
        );
        close(
            report.entropy.conduction_generation_j_k,
            expected * (1.0 / 400.0 - 1.0 / 500.0),
            1e-10,
        );
        close(report.entropy.temporal_defect_j_k, 0.0, 1e-10);
        assert_eq!(before.cells(), &original);
    }
}

// Independent dimensional steady balance: cumulative generated power through
// each radial face, integrated against dr/(k A). No production geometry helper.
fn steady_reference(shape: LayerStackGeometry, n: usize) -> (Vec<f64>, Vec<f64>) {
    let dimension = if matches!(shape, LayerStackGeometry::Cylindrical { .. }) {
        2.0
    } else {
        3.0
    };
    let radii = (0..=n)
        .map(|i| (i as f64 / n as f64).powf(1.0 / dimension))
        .collect::<Vec<_>>();
    let centres = radii
        .windows(2)
        .map(|r| centroid(shape, r[0], r[1]))
        .collect::<Vec<_>>();
    let mut temperature = vec![0.0; n];
    temperature[n - 1] =
        350.0 + 400.0 * volume(shape, 0.0, 1.0) * metric(shape, centres[n - 1], 1.0);
    for i in (0..n - 1).rev() {
        temperature[i] = temperature[i + 1]
            + 400.0 * volume(shape, 0.0, radii[i + 1]) * metric(shape, centres[i], centres[i + 1]);
    }
    (centres, temperature)
}

#[test]
fn native_radial_heating_computes_a_melt_front_with_spatial_refinement_and_one_heat_account() {
    let card = material(1.0, 1.0);
    let curve = caloric(&card, 400.0, 100.0);
    let law = prepare(&card, &curve);
    for shape in solid_shapes() {
        let dimension: f64 = if matches!(shape, LayerStackGeometry::Cylindrical { .. }) {
            2.0
        } else {
            3.0
        };
        let exact_front = (1.0 - dimension / 4.0).sqrt();
        let mut previous_errors: Option<[f64; 2]> = None;
        for n in [8, 16, 32] {
            let body = solid(&[&law], &[1.0], &[n], &[], shape);
            let bound = body.bind(ids(), [None, Some(0.0)]).unwrap();
            let mut events = bound
                .birth_events(0.0, &vec![-50.0; n], [None, Some(350.0)], || Ok(()))
                .unwrap();
            for timed in &mut events {
                if let Event::Birth(cell) = &mut timed.event {
                    cell.boundary.source_w = 400.0 * cell.geometry.volume_m3;
                }
            }
            let before = DepositionSimulation::new(
                events,
                &bound.materials(),
                DepositionBudget {
                    max_cells: n,
                    max_links: n - 1,
                    max_events: 2 * n - 1,
                    max_substeps: 100_000,
                },
            )
            .unwrap()
            .advance_to(
                0.0,
                StepConfig {
                    max_dt_s: 1.0,
                    energy_tolerance_j: 1e-9,
                },
                || Ok(()),
            )
            .unwrap()
            .0;
            let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
            let mut workspace = ImplicitConductionWorkspace::new(n).unwrap();
            let (next, report) = chain
                .advance_to(&before, 10.0, controls(0.1), &mut workspace, || Ok(()))
                .unwrap();
            let (centres, reference) = steady_reference(shape, n);
            let mut temperature = Vec::with_capacity(n);
            let mut stored = 0.0;
            let mut expected_stored = 0.0;
            let mut field_error: f64 = 0.0;
            for (i, &id) in bound.cell_ids().iter().enumerate() {
                let cell = &next.cells()[&id];
                let t = next.temperature_k(id).unwrap();
                close(t, reference[i], 2e-5);
                let exact_t = 350.0 + 400.0 * (1.0 - centres[i] * centres[i]) / (2.0 * dimension);
                field_error = field_error.max((t - exact_t).abs());
                let expected_h = if reference[i] > 400.0 {
                    reference[i] - 400.0 + 100.0
                } else {
                    reference[i] - 400.0
                };
                stored += cell.mass_kg * (cell.specific_h_j_kg + 50.0);
                expected_stored += cell.mass_kg * (expected_h + 50.0);
                let phase = curve
                    .state_at_specific_enthalpy(cell.specific_h_j_kg)
                    .unwrap();
                assert_eq!(
                    phase.liquid_mass_fraction(),
                    if t > 400.0 { 1.0 } else { 0.0 }
                );
                temperature.push(t);
            }
            close(stored, expected_stored, 1e-4);
            close(report.energy.stored_energy_change_j, stored, 1e-8);
            close(report.energy.external_heat_j, stored, 1e-8);
            close(report.energy.residual_j, 0.0, 1e-8);
            assert!(temperature.windows(2).all(|pair| pair[0] > pair[1]));
            let i = temperature
                .windows(2)
                .position(|pair| pair[0] > 400.0 && pair[1] < 400.0)
                .unwrap();
            let front = centres[i]
                + (400.0 - temperature[i]) * (centres[i + 1] - centres[i])
                    / (temperature[i + 1] - temperature[i]);
            let front_error = (front - exact_front).abs();
            if let Some(previous) = previous_errors {
                assert!(
                    field_error < 0.7 * previous[0],
                    "{field_error} vs {previous:?}"
                );
                assert!(
                    front_error < 0.3 * previous[1],
                    "{front_error} vs {previous:?}"
                );
            }
            previous_errors = Some([field_error, front_error]);
            assert_eq!(before.time_s(), 0.0);
            assert!(
                before
                    .cells()
                    .values()
                    .all(|cell| cell.specific_h_j_kg == -50.0)
            );
        }
    }
}

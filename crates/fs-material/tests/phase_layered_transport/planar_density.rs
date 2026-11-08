//! Original density changes the native implicit flux, not reference inventories.
use super::*;
use fs_material::phase::{HeatCapacityKnot, SolidLiquidPhase};
use fs_thermal::enthalpy::{MaterialId, ThermalLaw};

fn liquid_curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    // Manufactured affine Cp/rho(T), not empirical material accuracy.
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 10.0,
                bulk_density_kg_m3: 1000.0,
            },
            HeatCapacityKnot {
                temperature_k: 400.0,
                specific_heat_capacity_j_kg_k: 10.0,
                bulk_density_kg_m3: 500.0,
            },
        ],
    )
    .unwrap()
}
fn column<'a>(
    laws: &[&'a PreparedPhaseConductivity<'a>],
    widths: &[f64],
    h: &[f64],
    area: f64,
    contacts: &[f64],
) -> PreparedPhaseLayeredThermalBody<'a> {
    let layers: Vec<_> = laws
        .iter()
        .enumerate()
        .map(|(i, &law)| PhaseLayerSource {
            id: 100 + i as u64,
            thickness_m: widths[i],
            mass_kg: area
                * widths[i]
                * law
                    .caloric_curve()
                    .values_at_specific_enthalpy(h[i])
                    .unwrap()
                    .bulk_density_kg_m3,
            cell_count: 1,
            conductivity: law,
        })
        .collect();
    PreparedPhaseLayeredThermalBody::try_new(
        &layers,
        contacts,
        LayerStackGeometry::Planar { area_m2: area },
        || Ok(()),
    )
    .unwrap()
}
fn source_initial<'a>(
    bound: &BoundPhaseLayeredThermalBody<'_, 'a>,
    h: f64,
    specific_power: f64,
) -> DepositionSimulation<'a> {
    let mut events = bound.birth_events(0.0, &[h], [None; 2], || Ok(())).unwrap();
    if let Event::Birth(cell) = &mut events[0].event {
        cell.boundary.source_w = cell.mass_kg * specific_power;
    }
    DepositionSimulation::new(
        events,
        &bound.materials(),
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1000,
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
    .0
}

#[test]
fn two_liquid_sources_share_one_implicit_heat_transfer_at_current_thickness() {
    let (card_a, card_b) = (material(1.0, 1.0), material(2.0, 2.0));
    let (curve_a, curve_b) = (liquid_curve(&card_a), liquid_curve(&card_b));
    let (ka, kb) = (prepare(&card_a, &curve_a), prepare(&card_b, &curve_b));
    let (area, contact, dt) = (0.2, 0.002, 2.0);
    let body = column(&[&ka, &kb], &[0.01, 0.015], &[500.0, 0.0], area, &[contact]);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    let before = initial(&bound, &[500.0, 0.0], [None; 2]);
    let moving = bound
        .planar_liquid_density_column(&before, || Ok(()))
        .unwrap();
    let (ma, mb) = (
        before.cells()[&fs_thermal::enthalpy::CellId(11)].mass_kg,
        before.cells()[&fs_thermal::enthalpy::CellId(12)].mass_kg,
    );
    let equilibrium_heat = 50.0 / (1.0 / (10.0 * ma) + 1.0 / (10.0 * mb));
    let expected = root(0.0, equilibrium_heat, |q| {
        let (ha, hb) = (500.0 - q / ma, q / mb);
        let (rho_a, rho_b) = (1000.0 - 0.5 * ha, 1000.0 - 0.5 * hb);
        let resistance = ma / (2.0 * area * rho_a) + contact + mb / (4.0 * area * rho_b);
        q - dt * area * ((ha - hb) / 10.0) / resistance
    });
    let frozen_g = area / (0.005 + contact + 0.0075 / 2.0);
    let frozen_heat =
        dt * frozen_g * 50.0 / (1.0 + dt * frozen_g * (1.0 / (10.0 * ma) + 1.0 / (10.0 * mb)));
    assert!((expected - frozen_heat).abs() > 0.1);
    let chain = moving.implicit_chain(&before, || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (after, report) = chain
        .advance_to(&before, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    let shape = moving.current_geometry(&after, || Ok(())).unwrap();
    let states: Vec<_> = after.cells().values().collect();
    close(states[0].specific_h_j_kg, 500.0 - expected / ma, 2e-8);
    close(states[1].specific_h_j_kg, expected / mb, 2e-8);
    close(report.energy.external_heat_j, 0.0, 1e-12);
    close(report.energy.stored_energy_change_j, 0.0, 1e-8);
    close(report.energy.residual_j, 0.0, 1e-8);
    assert_eq!(report.energy.substeps, 1);
    assert!(report.nonlinear_iterations > 0);
    close(shape.footprint_area_m2, area, 0.0);
    close(shape.cells[0].interval_m[0], 0.0, 0.0);
    close(
        shape.cells[0].interval_m[1],
        shape.cells[1].interval_m[0],
        0.0,
    );
    let mut thickness = 0.0;
    for (i, cell) in states.into_iter().enumerate() {
        assert_eq!(cell.mass_kg, before.cells()[&cell.id].mass_kg);
        assert_eq!(cell.geometry, before.cells()[&cell.id].geometry);
        let original = if i == 0 { &curve_a } else { &curve_b };
        assert!(after.uses_material_law(cell.material_id, original));
        let volume = cell.mass_kg / (1000.0 - 0.5 * cell.specific_h_j_kg);
        close(shape.cells[i].volume_m3, volume, 1e-14);
        thickness += volume / area;
    }
    close(shape.total_thickness_m, thickness, 1e-14);
    close(shape.total_mass_kg, ma + mb, 1e-14);
}

#[test]
fn both_following_end_baths_use_the_same_current_bulk_half_thickness() {
    let card = material(1.0, 1.0);
    let curve = liquid_curve(&card);
    let law = prepare(&card, &curve);
    let (area, dt) = (0.1, 2.0);
    let body = column(&[&law], &[0.02], &[0.0], area, &[]);
    let bound = body.bind(ids(), [Some(0.002), Some(0.005)]).unwrap();
    let before = initial(&bound, &[0.0], [Some(360.0), Some(390.0)]);
    let moving = bound
        .planar_liquid_density_column(&before, || Ok(()))
        .unwrap();
    let mass = before.mass_kg();
    let expected_h = root(0.0, 900.0, |h| {
        let temperature = 300.0 + h / 10.0;
        let half = mass / (2.0 * area * (1000.0 - 0.5 * h));
        mass * h
            - dt * area
                * ((360.0 - temperature) / (half + 0.002) + (390.0 - temperature) / (half + 0.005))
    });
    let chain = moving.implicit_chain(&before, || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (after, report) = chain
        .advance_to(&before, dt, controls(dt), &mut workspace, || Ok(()))
        .unwrap();
    let cell = after.cells().values().next().unwrap();
    close(cell.specific_h_j_kg, expected_h, 2e-8);
    close(report.energy.external_heat_j, mass * expected_h, 1e-7);
    close(
        report.energy.stored_energy_change_j,
        report.energy.external_heat_j,
        1e-8,
    );
    close(report.energy.residual_j, 0.0, 1e-8);
    assert_eq!(cell.boundary, before.cells()[&cell.id].boundary);
    let geometry = moving.current_geometry(&after, || Ok(())).unwrap();
    close(
        geometry.total_thickness_m,
        mass / (area * (1000.0 - 0.5 * expected_h)),
        1e-12,
    );
    assert!(geometry.total_thickness_m > 0.02);
}

#[test]
fn isolated_inventory_uses_native_adaptive_sdirk_without_a_dummy_heat_channel() {
    let card = material(1.0, 1.0);
    let curve = liquid_curve(&card);
    let law = prepare(&card, &curve);
    let body = column(&[&law], &[0.02], &[0.0], 0.1, &[]);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    let before = source_initial(&bound, 0.0, 500.0);
    let moving = bound
        .planar_liquid_density_column(&before, || Ok(()))
        .unwrap();
    let chain = moving.implicit_chain(&before, || Ok(())).unwrap();
    let policy = AdaptiveImplicitConductionControls {
        step: controls(1.0),
        initial_dt_s: 1.0,
        maximum_trials: 16,
        enthalpy_absolute_tolerance_j_kg: 1e-7,
        temperature_absolute_tolerance_k: 1e-8,
        phase_fraction_absolute_tolerance: Some(1e-10),
        relative_tolerance: 1e-7,
    };
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (after, report) = chain
        .advance_sdirk2_adaptive_to(&before, 1.0, policy, &mut workspace, || Ok(()))
        .unwrap();
    let cell = after.cells().values().next().unwrap();
    close(cell.specific_h_j_kg, 500.0, 1e-8);
    close(after.temperature_k(cell.id).unwrap(), 350.0, 1e-9);
    close(
        report.solve.energy.external_heat_j,
        500.0 * cell.mass_kg,
        1e-8,
    );
    close(report.solve.energy.residual_j, 0.0, 1e-8);
    assert!(after.links().is_empty());
    assert!(cell.boundary.reservoir.is_none());
    assert_eq!(cell.geometry, before.cells()[&cell.id].geometry);
    let current = moving.current_geometry(&after, || Ok(())).unwrap();
    close(current.cells[0].density_kg_m3, 750.0, 1e-8);
    close(current.total_thickness_m, 0.02 * (1000.0 / 750.0), 1e-12);
    assert!(report.accepted_comparisons > 0);
}

#[test]
fn geometry_inventory_and_original_caloric_identity_are_required() {
    let card = material(1.0, 1.0);
    let curve = liquid_curve(&card);
    let law = prepare(&card, &curve);
    for shape in [
        LayerStackGeometry::Spherical {
            inner_radius_m: 0.1,
        },
        LayerStackGeometry::Cylindrical {
            inner_radius_m: 0.1,
            length_m: 1.0,
        },
    ] {
        let body = super::body(&[&law], &[0.01], &[1], &[], shape);
        let bound = body.bind(ids(), [None; 2]).unwrap();
        let thermal = initial(&bound, &[0.0], [None; 2]);
        assert!(
            bound
                .planar_liquid_density_column(&thermal, || Ok(()))
                .is_err()
        );
    }
    let body = column(&[&law], &[0.02], &[0.0], 0.1, &[]);
    let mut wrong_layers = body.layers().to_vec();
    wrong_layers[0].mass_kg *= 1.01;
    let wrong = PreparedPhaseLayeredThermalBody::try_new(
        &wrong_layers,
        &[],
        LayerStackGeometry::Planar { area_m2: 0.1 },
        || Ok(()),
    )
    .unwrap();
    let wrong_bound = wrong.bind(ids(), [None; 2]).unwrap();
    let wrong_thermal = initial(&wrong_bound, &[0.0], [None; 2]);
    assert!(
        wrong_bound
            .planar_liquid_density_column(&wrong_thermal, || Ok(()))
            .is_err()
    );
    let bound = body.bind(ids(), [None; 2]).unwrap();
    let copied_curve = curve.clone();
    let foreign = DepositionSimulation::new(
        bound
            .birth_events(0.0, &[0.0], [None; 2], || Ok(()))
            .unwrap(),
        &[(MaterialId(21), &copied_curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 10,
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
    assert!(
        bound
            .planar_liquid_density_column(&foreign, || Ok(()))
            .is_err()
    );
}

#[test]
fn native_trials_cannot_cross_into_solid_and_failed_advances_publish_nothing() {
    let card = material(1.0, 1.0);
    let curve = caloric(&card, 300.0, 100.0);
    let law = prepare(&card, &curve);
    let body = column(&[&law], &[0.02], &[110.0], 0.1, &[]);
    let bound = body.bind(ids(), [None; 2]).unwrap();
    for h in [-10.0, 50.0] {
        let thermal = initial(&bound, &[h], [None; 2]);
        assert!(
            bound
                .planar_liquid_density_column(&thermal, || Ok(()))
                .is_err()
        );
    }
    let before = source_initial(&bound, 110.0, -20.0);
    let moving = bound
        .planar_liquid_density_column(&before, || Ok(()))
        .unwrap();
    let chain = moving.implicit_chain(&before, || Ok(())).unwrap();
    let original_cells = before.cells().clone();
    let original_geometry = moving.current_geometry(&before, || Ok(())).unwrap();
    let mut limits = controls(1.0);
    limits.maximum_halvings = 0;
    limits.maximum_newton_iterations = 8;
    limits.maximum_backtracks = 8;
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    assert!(
        chain
            .advance_to(&before, 1.0, limits, &mut workspace, || Ok(()))
            .is_err()
    );
    assert_eq!(before.cells(), &original_cells);
    assert_eq!(
        moving.current_geometry(&before, || Ok(())).unwrap(),
        original_geometry
    );
    let mut polls = 0;
    assert!(
        chain
            .advance_to(&before, 0.1, controls(0.001), &mut workspace, || {
                polls += 1;
                if polls == 100 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(polls, 100);
    assert_eq!(before.time_s(), 0.0);
    assert_eq!(before.cells(), &original_cells);
    assert_eq!(
        moving.current_geometry(&before, || Ok(())).unwrap(),
        original_geometry
    );
}

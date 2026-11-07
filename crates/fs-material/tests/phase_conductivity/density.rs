//! Latent density changes drive the same native implicit resistance and heat.
use super::*;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{
    ConductionObserver, ConductiveReservoir, PlanarDensityLayer,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, Reservoir, StepConfig, TimedEvent,
};

#[test]
fn original_phase_volume_tangent_clips_consistently_and_drives_implicit_latent_heat() {
    let card = card(vec![claim(
        phase_values(&[(0.0, 1.0), (1.0, 1.0)]),
        [0.0, 1.0],
    )]);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 300.0, 0.0, 1.0),
            (100.0, 400.0, 0.0, 1.0),
            (300.0, 400.0, 1.0, 0.5),
            (500.0, 600.0, 1.0, 0.5),
        ]
        .into_iter()
        .map(|(h, t, f, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: rho,
        })
        .collect(),
    )
    .unwrap();
    let path = prepare(&card, &curve, [0.0, 500.0], || Ok(())).unwrap();
    let volume_slope = |path: &PreparedPhaseConductivity<'_>, h| {
        ImplicitConductionLaw::specific_volume_derivative_at_specific_enthalpy(path, h)
            .unwrap()
            .unwrap()
    };
    close(volume_slope(&path, 100.0), 0.005);
    close(volume_slope(&path, 200.0), 0.005);
    close(volume_slope(&path, 300.0), 0.0);
    let solid = prepare(&card, &curve, [0.0, 100.0], || Ok(())).unwrap();
    let latent = prepare(&card, &curve, [100.0, 300.0], || Ok(())).unwrap();
    close(volume_slope(&solid, 100.0), 0.0);
    close(volume_slope(&latent, 300.0), 0.005);
    assert!(
        ImplicitConductionLaw::specific_volume_derivative_at_specific_enthalpy(&solid, 101.0)
            .is_err()
    );

    let initial = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(0),
                material_id: MaterialId(0),
                body_id: 0,
                segment_id: 0,
                mass_kg: 1.0,
                specific_h_j_kg: 100.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 2.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: 1.0,
                        conductance_w_k: 0.0,
                        temperature_k: 500.0,
                    }),
                },
            }),
        }],
        &[(MaterialId(0), &curve)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 100,
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
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(0), &path)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_material_density_layers(
        &initial,
        &[PlanarDensityLayer {
            cell_id: CellId(0),
            footprint_area_m2: 1.0,
        }],
        || Ok(()),
    )
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&initial, &observer, &[(MaterialId(0), &path)], || Ok(()))
            .unwrap();
    let (melted, report) = chain
        .advance_to(
            &initial,
            1.0,
            ImplicitConductionControls {
                maximum_dt_s: 1.0,
                maximum_halvings: 0,
                cell_tolerance_j: 1e-10,
                total_tolerance_j: 1e-9,
                ..ImplicitConductionControls::default()
            },
            &mut ImplicitConductionWorkspace::new(1).unwrap(),
            || Ok(()),
        )
        .unwrap();

    // x=h-100; v=1+x/200; G=2/v, T=400. BE requires x=200/v.
    // The constant-density result x=200 would overpredict complete melting.
    let expected_h = 100.0 * 5.0_f64.sqrt();
    let actual_h = melted.cells()[&CellId(0)].specific_h_j_kg;
    assert!((actual_h - expected_h).abs() < 1e-9);
    assert_eq!(melted.temperature_k(CellId(0)).unwrap(), 400.0);
    close(
        curve
            .values_at_specific_enthalpy(actual_h)
            .unwrap()
            .liquid_mass_fraction,
        (expected_h - 100.0) / 200.0,
    );
    assert!((report.energy.external_heat_j - (expected_h - 100.0)).abs() < 1e-9);
    assert!(report.energy.residual_j.abs() < 1e-9);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(melted.cells()[&CellId(0)].mass_kg, 1.0);
    assert_eq!(
        melted.cells()[&CellId(0)].geometry,
        initial.cells()[&CellId(0)].geometry
    );
}

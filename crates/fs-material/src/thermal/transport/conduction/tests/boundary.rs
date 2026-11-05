//! G1/G4: native Cp/Cv storage, integrated material faces and actual bath heat.

use super::super::boundary::IntegratedCaloricBoundary;
use super::*;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
};
use fs_thermal::enthalpy::{
    Reservoir,
    conduction::{ConductiveAmbient, ConductiveReservoir},
};

const CELL: CellId = CellId(0);
const MATERIAL: MaterialId = MaterialId(11);
const MASS: f64 = 0.5;

fn single_cell<'a>(
    law: &'a dyn ThermalLaw,
    initial_energy: f64,
    boundary: Boundary,
) -> DepositionSimulation<'a> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 1,
                body_id: 1,
                material_id: MATERIAL,
                mass_kg: MASS,
                specific_h_j_kg: initial_energy,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [MASS / (RHO * AREA), 0.0, 0.0],
                    volume_m3: MASS / RHO,
                    surface_area_m2: 3.0 * AREA,
                },
                boundary,
            }),
        }],
        &[(MATERIAL, law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 8,
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

fn controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-6,
        total_tolerance_j: 1e-6,
        ..ImplicitConductionControls::default()
    }
}

#[test]
fn g1_original_cp_cv_boundary_matches_independent_film_root_and_complete_tangents() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_k =
        IsochoricConductionLaw::try_new(&cv_source.card, &cv, &cv_source.conductivity, SELECTION)
            .unwrap();
    let cp_k =
        IsobaricConductionLaw::try_new(&cp_source.card, &cp, &cp_source.conductivity, SELECTION)
            .unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&cv_k).unwrap(),
        CaloricConductionLayer::from_isobaric(&cp_k).unwrap(),
    ];
    let geometry = ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: 0.015,
        surface_resistance_m2_k_w: 0.004,
    };
    for (side, layer) in layers.iter().enumerate() {
        let boundary = IntegratedCaloricBoundary::new(layer, geometry, layer_controls()).unwrap();
        assert!(core::ptr::eq(boundary.caloric_law(), layer.caloric_law()));
        let k = [0.04, 0.02][side];
        for [t, bath] in [[620.0, 340.0], [340.0, 620.0], [500.0, 500.0]] {
            // For k(T)=k*T, p*Ts²+Ts=Tbath+p*T² exactly, with p=R''*k/(2L).
            // This independent quadratic never calls the production layer root.
            let p = geometry.surface_resistance_m2_k_w * k / (2.0 * geometry.half_distance_m);
            let rhs = bath + p * t * t;
            let face = 2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt());
            let expected = geometry.area_m2 / geometry.surface_resistance_m2_k_w * (face - bath);
            let response = boundary
                .response_at_specific_energy(energy(side, t), bath, &mut || Ok(()))
                .unwrap();
            near(response.material_face_temperature_k, face, 2e-10);
            near(response.power_to_reservoir_w, expected, 3e-9);
            let denominator =
                geometry.half_distance_m + geometry.surface_resistance_m2_k_w * k * face;
            near(
                response.endpoint_derivative_w_k[0],
                AREA * k * t / denominator,
                2e-10,
            );
            near(
                response.endpoint_derivative_w_k[1],
                -AREA * k * face / denominator,
                2e-10,
            );
            let state = boundary
                .state_at_specific_enthalpy(energy(side, t), bath, geometry, &mut || Ok(()))
                .unwrap();
            let capacity = if side == 0 {
                200.0 + 0.5 * t
            } else {
                400.0 + t
            };
            near(
                state.power_enthalpy_derivative_w_kg_j,
                response.endpoint_derivative_w_k[0] / capacity,
                2e-12,
            );
            let upper = boundary
                .response_at_specific_energy(energy(side, t) + 0.25, bath, &mut || Ok(()))
                .unwrap();
            let lower = boundary
                .response_at_specific_energy(energy(side, t) - 0.25, bath, &mut || Ok(()))
                .unwrap();
            near(
                state.power_enthalpy_derivative_w_kg_j,
                (upper.power_to_reservoir_w - lower.power_to_reservoir_w) / 0.5,
                2e-8,
            );
        }
    }
}

#[test]
fn g1_native_reservoir_and_ambient_heat_use_original_storage_and_each_bath_once() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_k =
        IsochoricConductionLaw::try_new(&cv_source.card, &cv, &cv_source.conductivity, SELECTION)
            .unwrap();
    let cp_k =
        IsobaricConductionLaw::try_new(&cp_source.card, &cp, &cp_source.conductivity, SELECTION)
            .unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&cv_k).unwrap(),
        CaloricConductionLayer::from_isobaric(&cp_k).unwrap(),
    ];
    let conductivity: [&dyn ConductivityLaw; 2] = [&cv_k, &cp_k];
    let caloric_conduction: [&dyn ImplicitConductionLaw; 2] = [&cv_k, &cp_k];
    let (initial, dt) = (320.0, 8.0);
    let mut temperatures = [0.0; 2];
    for side in 0..2 {
        // Deliberately unrelated authored coefficients must be REPLACED by
        // the selected physical faces, never added as a second heat source.
        let baths = Boundary {
            ambient_k: 680.0,
            htc_w_m2_k: 999.0,
            source_w: 0.0,
            reservoir: Some(Reservoir {
                area_m2: AREA,
                conductance_w_k: 999.0,
                temperature_k: 360.0,
            }),
        };
        let thermal = single_cell(layers[side].caloric_law(), energy(side, initial), baths);
        let reservoir_geometry = ImplicitBoundaryGeometry {
            area_m2: AREA,
            half_distance_m: 0.02,
            surface_resistance_m2_k_w: 0.0,
        };
        let ambient_geometry = ImplicitBoundaryGeometry {
            area_m2: thermal.exposed_area_m2(CELL).unwrap(),
            half_distance_m: 0.01,
            surface_resistance_m2_k_w: 0.0,
        };
        let reservoir =
            IntegratedCaloricBoundary::new(&layers[side], reservoir_geometry, layer_controls())
                .unwrap();
        let ambient =
            IntegratedCaloricBoundary::new(&layers[side], ambient_geometry, layer_controls())
                .unwrap();
        let observer = ConductionObserver::new(
            &thermal,
            &[(MATERIAL, conductivity[side])],
            &[],
            &[ConductiveReservoir {
                cell_id: CELL,
                half_distance_m: 0.02,
                surface_resistance_m2_k_w: 0.0,
            }],
            || Ok(()),
        )
        .unwrap()
        .with_ambient_faces(
            &thermal,
            &[ConductiveAmbient {
                cell_id: CELL,
                half_distance_m: 0.01,
                surface_resistance_m2_k_w: 0.0,
            }],
            || Ok(()),
        )
        .unwrap();
        let chain = ImplicitConductionChain::new(
            &thermal,
            &observer,
            &[(MATERIAL, caloric_conduction[side])],
            || Ok(()),
        )
        .unwrap()
        .with_boundary_laws(&thermal, &[(CELL, &reservoir)], &[(CELL, &ambient)], || {
            Ok(())
        })
        .unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let mut calls = 0;
        let (next, report) = chain
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || {
                calls += 1;
                Ok(())
            })
            .unwrap();
        let k = [0.04, 0.02][side];
        let gr = reservoir_geometry.area_m2 / reservoir_geometry.half_distance_m * k / 2.0;
        let ga = ambient_geometry.area_m2 / ambient_geometry.half_distance_m * k / 2.0;
        let (a, b) = if side == 0 {
            (200.0, 0.25)
        } else {
            (400.0, 0.5)
        };
        let quadratic = MASS * b + dt * (gr + ga);
        let linear = MASS * a;
        let rhs = MASS * (a * initial + b * initial * initial)
            + dt * (gr * 360.0_f64.powi(2) + ga * 680.0_f64.powi(2));
        let expected = 2.0 * rhs / (linear + (linear * linear + 4.0 * quadratic * rhs).sqrt());
        let t = next.temperature_k(CELL).unwrap();
        temperatures[side] = t;
        near(t, expected, 2e-8);
        assert!(t > 475.0 && t < 650.0);
        let qr = dt * gr * (360.0_f64.powi(2) - t * t);
        let qa = dt * ga * (680.0_f64.powi(2) - t * t);
        near(report.energy.external_heat_j, qr + qa, 1e-6);
        near(
            report.energy.stored_energy_change_j,
            MASS * (energy(side, t) - energy(side, initial)),
            2e-6,
        );
        near(report.energy.residual_j, 0.0, 1e-6);
        near(
            report.entropy.reservoir_entropy_input_j_k,
            qr / 360.0 + qa / 680.0,
            2e-8,
        );
        near(
            report.entropy.thermal_storage_change_j_k,
            MASS * (a * (t / initial).ln() + 2.0 * b * (t - initial)),
            2e-8,
        );
        assert_eq!(report.energy.substeps, 1);
        assert_eq!(next.cells()[&CELL].boundary, baths);
        assert_eq!(
            next.material_storage_basis(MATERIAL).unwrap(),
            [
                ThermalStorageBasis::InternalEnergy,
                ThermalStorageBasis::Enthalpy
            ][side]
        );

        // Final cancellation after all physical work publishes no candidate;
        // retrying the identical interval and scratch reproduces the endpoint.
        let mut cancelled_calls = 0;
        assert!(matches!(
            chain.advance_to(&thermal, dt, controls(dt), &mut workspace, || {
                cancelled_calls += 1;
                if cancelled_calls == calls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(
            thermal.cells()[&CELL].specific_h_j_kg,
            energy(side, initial)
        );
        let (retry, retry_report) = chain
            .advance_to(&thermal, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        assert_eq!(retry.cells(), next.cells());
        assert_eq!(retry_report, report);
    }
    assert!(temperatures[0] > temperatures[1] + 70.0);
}

#[test]
fn g4_external_bath_is_not_a_material_sample_and_failed_source_steps_preserve_history() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conductivity =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION)
            .unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conductivity).unwrap();
    let geometry = ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: 0.02,
        surface_resistance_m2_k_w: 0.1,
    };
    let boundary = IntegratedCaloricBoundary::new(&layer, geometry, layer_controls()).unwrap();
    // The bath is outside [300,700] K, while the finite-film material face
    // remains supported. Its analytic root is independent of source queries.
    let response = boundary
        .response_at_specific_energy(energy(0, 500.0), 1000.0, &mut || Ok(()))
        .unwrap();
    let p = geometry.surface_resistance_m2_k_w * 0.04 / (2.0 * geometry.half_distance_m);
    let rhs = 1000.0 + p * 500.0_f64.powi(2);
    let face = 2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt());
    near(response.material_face_temperature_k, face, 2e-10);
    assert!(response.power_to_reservoir_w < 0.0 && face > 500.0 && face < 700.0);
    let perfect = IntegratedCaloricBoundary::new(
        &layer,
        ImplicitBoundaryGeometry {
            surface_resistance_m2_k_w: 0.0,
            ..geometry
        },
        layer_controls(),
    )
    .unwrap();
    assert!(
        perfect
            .response_at_specific_energy(energy(0, 500.0), 1000.0, &mut || Ok(()))
            .is_err()
    );
    assert!(
        boundary
            .response_at_specific_energy(layer.energy_bounds_j_kg()[1] + 1.0, 500.0, &mut || Ok(()))
            .is_err()
    );
    assert!(
        boundary
            .state_at_specific_enthalpy(
                energy(0, 500.0),
                1000.0,
                ImplicitBoundaryGeometry {
                    area_m2: 2.0 * AREA,
                    ..geometry
                },
                &mut || Ok(())
            )
            .is_err()
    );

    let thermal = single_cell(
        &cv,
        energy(0, 690.0),
        Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            source_w: 0.0,
            reservoir: Some(Reservoir {
                area_m2: AREA,
                conductance_w_k: 0.0,
                temperature_k: 1000.0,
            }),
        },
    );
    let observer = ConductionObserver::new(
        &thermal,
        &[(MATERIAL, &conductivity)],
        &[],
        &[ConductiveReservoir {
            cell_id: CELL,
            half_distance_m: geometry.half_distance_m,
            surface_resistance_m2_k_w: geometry.surface_resistance_m2_k_w,
        }],
        || Ok(()),
    )
    .unwrap();
    let chain =
        ImplicitConductionChain::new(&thermal, &observer, &[(MATERIAL, &conductivity)], || Ok(()))
            .unwrap()
            .with_boundary_laws(&thermal, &[(CELL, &boundary)], &[], || Ok(()))
            .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    assert!(
        chain
            .advance_to(
                &thermal,
                1000.0,
                controls(1000.0),
                &mut workspace,
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CELL].specific_h_j_kg, energy(0, 690.0));
    let (retry, _) = chain
        .advance_to(&thermal, 0.1, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();
    assert!(retry.temperature_k(CELL).unwrap() > 690.0);
    assert!(retry.temperature_k(CELL).unwrap() < 700.0);
}

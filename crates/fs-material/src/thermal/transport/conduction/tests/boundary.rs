//! G1/G4: native Cp/Cv storage, integrated material faces and actual bath heat.

use super::super::boundary::IntegratedCaloricBoundary;
use super::*;
use crate::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
};
use crate::thermal::conduction::phase::PreparedPhaseConductivity;
use crate::thermal::source::ThermalMaterialPoint;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
};
use fs_thermal::enthalpy::{
    Reservoir,
    conduction::{ConductiveAmbient, ConductiveReservoir, PlanarDensityLayer},
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

fn phase_point(source: &Fixture) -> ThermalMaterialPoint<'_> {
    ThermalMaterialPoint {
        material: &source.card,
        context: &source.lower,
        temperature_axis: "T",
        liquid_fraction_axis: None,
    }
}

fn exact_phase_curve(source: &Fixture, phase: SolidLiquidPhase) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        source.card.content_hash(),
        phase,
        0.0,
        &[300.0, 500.0, 700.0].map(|t| HeatCapacityKnot {
            temperature_k: t,
            specific_heat_capacity_j_kg_k: 400.0 + t,
            bulk_density_kg_m3: 1600.0 - 2.0 * t,
        }),
    )
    .unwrap()
}

#[test]
fn g1_equilibrium_layers_keep_exact_cp_source_bounds_and_inward_phase_tangents() {
    let source = Fixture::new(false);
    for phase in [SolidLiquidPhase::Solid, SolidLiquidPhase::Liquid] {
        let curve = exact_phase_curve(&source, phase);
        let resolved = ResolvedConductionLayer::new(
            &curve, phase_point(&source), &source.conductivity, SELECTION,
            phase, [350.0, 500.0], || Ok(()),
        ).unwrap();
        let layer = CaloricConductionLayer::from_equilibrium_phase(&resolved);
        assert!(core::ptr::addr_eq(layer.caloric_law(), &curve as &dyn ThermalLaw));
        assert!(core::ptr::eq(layer.source_line(), resolved.source_line().as_ref()));
        assert!(layer.supports_material_card(source.card.content_hash()));
        assert_eq!(layer.caloric_law().storage_basis(), ThermalStorageBasis::Enthalpy);
        assert_eq!(layer.energy_bounds_j_kg(), [energy(1, 350.0), energy(1, 500.0)]);
        assert_eq!(layer.conductivity_curve(), resolved.conductivity_curve());
        for (initial_t, bath, target_t) in [
            (400.0_f64, 450.0, 420.0),
            (450.0, 350.0, 420.0),
        ] {
            // Independent integral of Cp(T)/(T_bath-T), Cp(T)=400+T.
            let x = (400.0 + bath) * ((bath - initial_t) / (bath - target_t)).ln()
                - (target_t - initial_t);
            let exchange = layer.linear_heat_exchange(energy(1, initial_t), bath, x).unwrap();
            let heat = energy(1, target_t) - energy(1, initial_t);
            let mean = bath - heat / x;
            assert_eq!(exchange.storage_basis, ThermalStorageBasis::Enthalpy);
            assert_eq!(exchange.initial_specific_energy_j_kg, energy(1, initial_t));
            near(exchange.final_specific_energy_j_kg, energy(1, target_t), 2e-7);
            near(exchange.temperature_k, target_t, 2e-10);
            near(exchange.specific_heat_j_kg, heat, 2e-7);
            near(exchange.heat_residual_j_kg, 0.0, 2e-9);
            near(exchange.mean_temperature_k, mean, 2e-10);
            near(exchange.inlet_energy_gain, (bath - target_t) / (bath - initial_t), 3e-14);
            near(exchange.remaining_temperature_difference_k, bath - target_t, 3e-11);
            near(exchange.mean_log_exposure_derivative_k, target_t - mean, 3e-10);
        }
        // The original caloric chart extends farther than this k-source layer.
        // Fixed-bath convenience must not silently widen the prepared support.
        assert!(layer.linear_heat_exchange(energy(1, 400.0), 700.0, 10_000.0).is_err());
        assert!(layer.linear_heat_exchange(energy(1, 300.0), 450.0, 0.0).is_err());
        let end = layer.state_at_specific_energy(energy(1, 500.0)).unwrap();
        assert_eq!(end.temperature_k, 500.0);
        near(end.temperature_slope_kg_k_j, 1.0 / 900.0, 1e-16);
        near(end.conductivity_slope_kg_w_m_k_j, 0.02 / 900.0, 1e-18);
        assert!(layer.state_at_specific_energy(energy(1, 350.0) - 1.0).is_err());
        assert!(layer.state_at_specific_energy(energy(1, 500.0) + 1.0).is_err());
    }

    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        source.card.content_hash(),
        [(0.0, 300.0, 0.0, 1000.0), (160_000.0, 500.0, 0.0, 600.0),
         (210_000.0, 500.0, 1.0, 600.0), (410_000.0, 700.0, 1.0, 200.0)]
            .map(|(h, t, f, rho)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h, temperature_k: t,
                liquid_mass_fraction: f, bulk_density_kg_m3: rho,
            }).to_vec(),
    ).unwrap();
    // The original solidus has a zero RIGHT derivative into its latent plateau.
    // The selected all-solid layer must retain its positive inward derivative.
    assert_eq!(curve.temperature_derivative_at_specific_enthalpy(160_000.0).unwrap(), 0.0);
    assert_eq!(curve.temperature_at_specific_enthalpy(185_000.0).unwrap(), 500.0);
    for (phase, temperatures, bounds, end_slope) in [
        (SolidLiquidPhase::Solid, [350.0, 500.0], [40_000.0, 160_000.0], 1.0 / 800.0),
        (SolidLiquidPhase::Liquid, [500.0, 650.0], [210_000.0, 360_000.0], 1.0 / 1000.0),
    ] {
        let resolved = ResolvedConductionLayer::new(
            &curve, phase_point(&source), &source.conductivity, SELECTION,
            phase, temperatures, || Ok(()),
        ).unwrap();
        let layer = CaloricConductionLayer::from_equilibrium_phase(&resolved);
        assert_eq!(layer.energy_bounds_j_kg(), bounds);
        near(layer.state_at_specific_energy(bounds[1]).unwrap().temperature_slope_kg_k_j,
            end_slope, 1e-16);
        assert!(layer.state_at_specific_energy(185_000.0).is_err());
        assert!(layer.linear_heat_exchange(185_000.0, 550.0, 1.0).is_err());
    }
    assert!(ResolvedConductionLayer::new(
        &curve, phase_point(&source), &source.conductivity, SELECTION,
        SolidLiquidPhase::SolidLiquid, [350.0, 650.0], || Ok(()),
    ).is_err());
    assert!(ResolvedConductionLayer::new(
        &curve, phase_point(&source), &source.conductivity, SELECTION,
        SolidLiquidPhase::Solid, [450.0, 550.0], || Ok(()),
    ).is_err());
}

#[test]
fn g1_original_phase_calorics_expand_and_contract_through_integrated_film_heat() {
    let source = Fixture::new(false);
    for (phase, initial, bath, temperatures) in [
        (SolidLiquidPhase::Solid, 400.0, 460.0, [350.0, 500.0]),
        (SolidLiquidPhase::Liquid, 620.0, 560.0, [500.0, 680.0]),
    ] {
        let curve = exact_phase_curve(&source, phase);
        let resolved = ResolvedConductionLayer::new(
            &curve, phase_point(&source), &source.conductivity, SELECTION,
            phase, temperatures, || Ok(()),
        ).unwrap();
        let layer = CaloricConductionLayer::from_equilibrium_phase(&resolved);
        let point_law = PreparedPhaseConductivity::new(
            &curve, phase_point(&source), &source.conductivity, SELECTION,
            layer.energy_bounds_j_kg(), || Ok(()),
        ).unwrap();
        let density = |t: f64| 1600.0 - 2.0 * t;
        let distance = |t: f64| MASS / (2.0 * AREA * density(t));
        let geometry = ImplicitBoundaryGeometry {
            area_m2: AREA,
            half_distance_m: distance(initial),
            surface_resistance_m2_k_w: 0.02,
        };
        let film = IntegratedCaloricBoundary::new(&layer, geometry, layer_controls()).unwrap();
        // Independent exact Kirchhoff/film root for k(T)=0.02*T, with the
        // original affine rho(T) and the same conserved mass/fixed footprint.
        let power = |t: f64, moving: bool| {
            let d = if moving { distance(t) } else { distance(initial) };
            let p = geometry.surface_resistance_m2_k_w * 0.02 / (2.0 * d);
            let rhs = bath + p * t * t;
            let face = 2.0 * rhs / (1.0 + (1.0 + 4.0 * p * rhs).sqrt());
            AREA / geometry.surface_resistance_m2_k_w * (face - bath)
        };
        let t = initial + 12.0 * (bath - initial).signum();
        let h = energy(1, t);
        let state = film.state_at_specific_enthalpy_with_geometry_tangent(
            h, bath,
            ImplicitBoundaryGeometry { half_distance_m: distance(t), ..geometry },
            2.0 / (density(t) * (400.0 + t)), &mut || Ok(()),
        ).unwrap().unwrap();
        near(state.conductance_w_k * (t - bath), power(t, true), 2e-9);
        near(
            state.power_enthalpy_derivative_w_kg_j,
            (power(temperature(1, h + 0.25), true)
                - power(temperature(1, h - 0.25), true)) / 0.5,
            2e-8,
        );

        let volume = MASS / density(initial);
        let original = DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CELL, material_id: MATERIAL, body_id: 1, segment_id: 1,
                    mass_kg: MASS, specific_h_j_kg: energy(1, initial),
                    geometry: CellGeometry {
                        start_m: [0.0; 3], end_m: [volume / AREA, 0.0, 0.0],
                        volume_m3: volume, surface_area_m2: 2.0 * AREA,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0, htc_w_m2_k: 0.0, source_w: 0.0,
                        reservoir: Some(Reservoir {
                            area_m2: AREA, conductance_w_k: 0.0, temperature_k: bath,
                        }),
                    },
                }),
            }],
            &[(MATERIAL, &curve)],
            DepositionBudget { max_cells: 1, max_links: 0, max_events: 1, max_substeps: 8 },
        ).unwrap().advance_to(
            0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 }, || Ok(()),
        ).unwrap().0;
        let observer = ConductionObserver::new(
            &original, &[(MATERIAL, &point_law)], &[],
            &[ConductiveReservoir {
                cell_id: CELL, half_distance_m: geometry.half_distance_m,
                surface_resistance_m2_k_w: geometry.surface_resistance_m2_k_w,
            }], || Ok(()),
        ).unwrap().with_material_density_layers(
            &original, &[PlanarDensityLayer { cell_id: CELL, footprint_area_m2: AREA }],
            || Ok(()),
        ).unwrap();
        let chain = ImplicitConductionChain::new(
            &original, &observer, &[(MATERIAL, &point_law)], || Ok(()),
        ).unwrap().with_boundary_laws(
            &original, &[(CELL, &film)], &[], || Ok(()),
        ).unwrap();
        let dt = 200.0;
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = chain.advance_to(
            &original, dt, controls(dt), &mut workspace, || Ok(()),
        ).unwrap();

        let residual = |t: f64, moving: bool| {
            MASS * (energy(1, t) - energy(1, initial)) + dt * power(t, moving)
        };
        let (mut lo, mut hi) = (initial.min(bath), initial.max(bath));
        assert!(residual(lo, true) < 0.0 && residual(hi, true) > 0.0);
        for _ in 0..70 {
            let t = 0.5 * (lo + hi);
            if residual(t, true) > 0.0 { hi = t; } else { lo = t; }
        }
        let final_t = next.temperature_k(CELL).unwrap();
        near(final_t, 0.5 * (lo + hi), 2e-8);
        near(report.energy.external_heat_j, -dt * power(final_t, true), 1e-6);
        near(report.energy.stored_energy_change_j,
            MASS * (energy(1, final_t) - energy(1, initial)), 1e-6);
        near(report.energy.residual_j, 0.0, 1e-6);
        assert!(residual(final_t, false).abs() > 1.0);
        assert_eq!(next.mass_kg(), MASS);
        assert_eq!(next.cells()[&CELL].geometry, original.cells()[&CELL].geometry);
        assert_eq!(original.cells()[&CELL].specific_h_j_kg, energy(1, initial));
    }
}

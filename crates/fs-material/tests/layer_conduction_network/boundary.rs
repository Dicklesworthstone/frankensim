//! Material-source boundary fluxes in the original implicit enthalpy history.
use super::*;
use fs_material::thermal::conduction::network::ImplicitLayerBoundary;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
};
use fs_thermal::enthalpy::conduction::{ConductiveAmbient, ConductiveReservoir};

const AREA: f64 = 0.004;
const DISTANCE: f64 = 0.02;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} differs from {expected} by more than {tolerance}"
    );
}

fn geometry(film: f64) -> ImplicitBoundaryGeometry {
    ImplicitBoundaryGeometry {
        area_m2: AREA,
        half_distance_m: DISTANCE,
        surface_resistance_m2_k_w: film,
    }
}

fn baths(ambient_k: f64, reservoir_k: Option<f64>, ambient: bool) -> Boundary {
    Boundary {
        ambient_k,
        // Selected material faces replace these unrelated authored values.
        htc_w_m2_k: if ambient { 999.0 } else { 0.0 },
        source_w: 0.0,
        reservoir: reservoir_k.map(|temperature_k| Reservoir {
            area_m2: AREA,
            conductance_w_k: 999.0,
            temperature_k,
        }),
    }
}

fn one_cell<'a>(
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    temperature: f64,
    boundary: Boundary,
    mut later_events: Vec<TimedEvent>,
) -> DepositionSimulation<'a> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(1),
            segment_id: 1,
            body_id: 1,
            material_id: MaterialId(1),
            mass_kg: 0.01,
            specific_h_j_kg: 100.0 * (temperature - 200.0),
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [0.1, 0.0, 0.0],
                volume_m3: 1e-5,
                surface_area_m2: AREA + boundary.reservoir.map_or(0.0, |r| r.area_m2),
            },
            boundary,
        }),
    }];
    events.append(&mut later_events);
    let budget = DepositionBudget {
        max_cells: 1,
        max_links: 0,
        max_events: events.len(),
        max_substeps: 100,
    };
    DepositionSimulation::new(events, &[(MaterialId(1), caloric)], budget)
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

fn observer<'a>(
    thermal: &DepositionSimulation<'_>,
    law: &'a dyn ConductivityLaw,
    shape: ImplicitBoundaryGeometry,
    ambient: bool,
) -> ConductionObserver<'a> {
    let reservoirs: Vec<_> = thermal.cells()[&CellId(1)]
        .boundary
        .reservoir
        .iter()
        .map(|_| ConductiveReservoir {
            cell_id: CellId(1),
            half_distance_m: shape.half_distance_m,
            surface_resistance_m2_k_w: shape.surface_resistance_m2_k_w,
        })
        .collect();
    let observer = ConductionObserver::new(
        thermal,
        &[(MaterialId(1), law)],
        &[],
        &reservoirs,
        || Ok(()),
    )
    .unwrap();
    if ambient {
        observer
            .with_ambient_faces(
                thermal,
                &[ConductiveAmbient {
                    cell_id: CellId(1),
                    half_distance_m: shape.half_distance_m,
                    surface_resistance_m2_k_w: shape.surface_resistance_m2_k_w,
                }],
                || Ok(()),
            )
            .unwrap()
    } else {
        observer
    }
}

// Independent Fourier integral for k(T)=2+0.01T. The cell has m*cp=1 J/K,
// so backward Euler reduces to a scalar quadratic with this positive root.
fn perfect_contact_endpoint(initial: f64, bath: f64, dt: f64) -> f64 {
    let scale = dt * AREA / DISTANCE;
    let a = 0.005 * scale;
    let b = 1.0 + 2.0 * scale;
    let c = -(initial + scale * (2.0 * bath + 0.005 * bath * bath));
    -2.0 * c / (b + (b * b - 4.0 * a * c).sqrt())
}

#[test]
fn g1_source_boundaries_match_analytic_enthalpy_in_chain_and_graph() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let law = transport(&card, &caloric, &context);
    let shape = geometry(0.0);
    let face = ImplicitLayerBoundary::new(&layer, shape, config()).unwrap();
    let bindings = [(CellId(1), &face as &dyn ImplicitBoundaryLaw)];
    for ambient in [false, true] {
        for (initial, bath) in [(600.0, 300.0), (300.0, 700.0), (450.0, 450.0)] {
            let old = one_cell(
                &caloric,
                initial,
                baths(bath, (!ambient).then_some(bath), ambient),
                vec![],
            );
            let observed = observer(&old, &law, shape, ambient);
            let (reservoirs, ambients) = if ambient {
                (&[][..], &bindings[..])
            } else {
                (&bindings[..], &[][..])
            };
            let chain =
                ImplicitConductionChain::new(&old, &observed, &[(MaterialId(1), &law)], || Ok(()))
                    .unwrap()
                    .with_boundary_laws(&old, reservoirs, ambients, || Ok(()))
                    .unwrap();
            let graph = ImplicitConductionGraph::new(
                &old,
                &observed,
                &[(MaterialId(1), &law)],
                ImplicitConductionLinearControls::default(),
                || Ok(()),
            )
            .unwrap()
            .with_boundary_laws(&old, reservoirs, ambients, || Ok(()))
            .unwrap();
            let dt = 0.3;
            let expected = perfect_contact_endpoint(initial, bath, dt);
            let mut work = ImplicitConductionWorkspace::new(1).unwrap();
            for result in [
                chain.advance_to(&old, dt, controls(dt), &mut work, || Ok(())),
                graph.advance_to(&old, dt, controls(dt), &mut work, || Ok(())),
            ] {
                let (next, report) = result.unwrap();
                near(next.temperature_k(CellId(1)).unwrap(), expected, 1e-8);
                near(report.energy.external_heat_j, expected - initial, 1e-8);
                near(
                    report.energy.stored_energy_change_j,
                    expected - initial,
                    1e-8,
                );
                near(report.energy.residual_j, 0.0, 1e-8);
                near(
                    report.entropy.reservoir_entropy_input_j_k,
                    (expected - initial) / bath,
                    1e-10,
                );
                near(
                    report.entropy.thermal_storage_change_j_k,
                    (expected / initial).ln(),
                    1e-10,
                );
                assert!(report.entropy.conduction_generation_j_k >= -1e-14);
                assert!(!report.entropy.has_unclassified_sources);
                assert_eq!(report.energy.substeps, 1);
            }
            assert_eq!(old.time_s(), 0.0);
            near(old.temperature_k(CellId(1)).unwrap(), initial, 0.0);
        }
    }
}

#[test]
fn g1_two_distinct_baths_keep_separate_source_integrals_and_entropy() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let law = transport(&card, &caloric, &context);
    let shape = geometry(0.0);
    let old = one_cell(&caloric, 480.0, baths(700.0, Some(300.0), true), vec![]);
    let observed = observer(&old, &law, shape, true);
    let face = ImplicitLayerBoundary::new(&layer, shape, config()).unwrap();
    let bindings = [(CellId(1), &face as &dyn ImplicitBoundaryLaw)];
    let chain = ImplicitConductionChain::new(&old, &observed, &[(MaterialId(1), &law)], || Ok(()))
        .unwrap()
        .with_boundary_laws(&old, &bindings, &bindings, || Ok(()))
        .unwrap();
    let mut work = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = chain
        .advance_to(&old, 0.25, controls(0.25), &mut work, || Ok(()))
        .unwrap();
    // At 500 K, independent integrals give -240 W to the 300 K bath and
    // +320 W from the 700 K bath: 20 J net input over the quarter second.
    near(next.temperature_k(CellId(1)).unwrap(), 500.0, 1e-8);
    near(report.energy.external_heat_j, 20.0, 1e-8);
    near(report.energy.residual_j, 0.0, 1e-8);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        -60.0 / 300.0 + 80.0 / 700.0,
        1e-10,
    );
    near(
        report.entropy.conduction_generation_j_k,
        -60.0 * (1.0 / 500.0 - 1.0 / 300.0) + 80.0 * (1.0 / 500.0 - 1.0 / 700.0),
        1e-10,
    );
}

#[test]
fn g1_finite_films_admit_external_baths_through_supported_material_faces() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let law = transport(&card, &caloric, &context);
    for (initial, endpoint, surface, bath) in
        [(650.0, 600.0, 400.0, 100.0), (350.0, 400.0, 600.0, 1100.0)]
    {
        let power = AREA / DISTANCE * (endpoint - surface) * (2.0 + 0.005 * (endpoint + surface));
        let film = (surface - bath) / power;
        let shape = geometry(film * AREA);
        let face = ImplicitLayerBoundary::new(&layer, shape, config()).unwrap();
        let h = 100.0 * (endpoint - 200.0);
        let sample = face
            .state_at_specific_enthalpy(h, bath, shape, &mut || Ok(()))
            .unwrap();
        near(sample.conductance_w_k * (endpoint - bath), power, 1e-8);
        let derivative = AREA / DISTANCE * (2.0 + 0.01 * endpoint)
            / (1.0 + AREA / DISTANCE * (2.0 + 0.01 * surface) * film)
            / 100.0;
        near(sample.power_enthalpy_derivative_w_kg_j, derivative, 1e-11);
        let delta = 0.01;
        let plus = face
            .state_at_specific_enthalpy(h + delta, bath, shape, &mut || Ok(()))
            .unwrap();
        let minus = face
            .state_at_specific_enthalpy(h - delta, bath, shape, &mut || Ok(()))
            .unwrap();
        near(
            (plus.conductance_w_k * (endpoint + delta / 100.0 - bath)
                - minus.conductance_w_k * (endpoint - delta / 100.0 - bath))
                / (2.0 * delta),
            derivative,
            1e-8,
        );
        let old = one_cell(&caloric, initial, baths(300.0, Some(bath), false), vec![]);
        let observed = observer(&old, &law, shape, false);
        let chain =
            ImplicitConductionChain::new(&old, &observed, &[(MaterialId(1), &law)], || Ok(()))
                .unwrap()
                .with_boundary_laws(&old, &[(CellId(1), &face)], &[], || Ok(()))
                .unwrap();
        let dt = (initial - endpoint) / power;
        let mut work = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = chain
            .advance_to(&old, dt, controls(dt), &mut work, || Ok(()))
            .unwrap();
        near(next.temperature_k(CellId(1)).unwrap(), endpoint, 1e-7);
        near(report.energy.external_heat_j, endpoint - initial, 1e-7);
        near(report.energy.residual_j, 0.0, 1e-8);
        let perfect = ImplicitLayerBoundary::new(&layer, geometry(0.0), config()).unwrap();
        assert!(
            perfect
                .state_at_specific_enthalpy(h, bath, geometry(0.0), &mut || Ok(()))
                .is_err()
        );
    }
}

#[test]
fn g3_bath_events_refresh_the_original_boundary_on_exact_native_intervals() {
    let card = source();
    let caloric = curve(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context);
    let law = transport(&card, &caloric, &context);
    let changed = baths(700.0, None, true);
    let old = one_cell(
        &caloric,
        600.0,
        baths(300.0, None, true),
        vec![TimedEvent {
            time_s: 0.2,
            event: Event::Boundary(CellId(1), changed),
        }],
    );
    let shape = geometry(0.0);
    let observed = observer(&old, &law, shape, true);
    let face = ImplicitLayerBoundary::new(&layer, shape, config()).unwrap();
    let chain = ImplicitConductionChain::new(&old, &observed, &[(MaterialId(1), &law)], || Ok(()))
        .unwrap()
        .with_boundary_laws(&old, &[], &[(CellId(1), &face)], || Ok(()))
        .unwrap();
    let mut work = ImplicitConductionWorkspace::new(1).unwrap();
    let mut polls = 0;
    let (next, report) = chain
        .advance_to(&old, 0.5, controls(0.5), &mut work, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let middle = perfect_contact_endpoint(600.0, 300.0, 0.2);
    let expected = perfect_contact_endpoint(middle, 700.0, 0.3);
    near(next.temperature_k(CellId(1)).unwrap(), expected, 1e-8);
    assert_eq!(report.energy.substeps, 2);
    assert_eq!(next.time_s(), 0.5);
    assert_eq!(next.cells()[&CellId(1)].boundary, changed);
    near(
        report.entropy.reservoir_entropy_input_j_k,
        (middle - 600.0) / 300.0 + (expected - middle) / 700.0,
        1e-10,
    );
    let mut calls = 0;
    assert!(matches!(
        chain.advance_to(&old, 0.5, controls(0.5), &mut work, || {
            calls += 1;
            if calls == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(
            fs_thermal::enthalpy::conduction::implicit::ImplicitConductionError::Thermal(
                DepositionError::Cancelled
            )
        )
    ));
    assert_eq!(old.time_s(), 0.0);
    near(old.temperature_k(CellId(1)).unwrap(), 600.0, 0.0);
    let retried = chain
        .advance_to(&old, 0.5, controls(0.5), &mut work, || Ok(()))
        .unwrap();
    assert_eq!(retried.0.cells(), next.cells());
    assert_eq!(retried.1, report);
}

#[test]
fn g0_boundary_retains_source_branch_identity_geometry_and_upper_endpoint_tangent() {
    let card = source();
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let caloric = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 200.0, 0.0),
            (40_000.0, 600.0, 0.0),
            (60_000.0, 600.0, 1.0),
            (100_000.0, 1000.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1000.0,
        })
        .collect(),
    )
    .unwrap();
    let solid = ResolvedConductionLayer::new(
        &caloric,
        ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        [200.0, 600.0],
        || Ok(()),
    )
    .unwrap();
    let shape = geometry(0.0);
    let face = ImplicitLayerBoundary::new(&solid, shape, config()).unwrap();
    let sample = face
        .state_at_specific_enthalpy(40_000.0, 300.0, shape, &mut || Ok(()))
        .unwrap();
    assert_eq!(
        caloric
            .temperature_derivative_at_specific_enthalpy(40_000.0)
            .unwrap(),
        0.0
    );
    near(sample.power_enthalpy_derivative_w_kg_j, 0.016, 1e-13);
    assert!(
        face.state_at_specific_enthalpy(40_000.0_f64.next_up(), 300.0, shape, &mut || Ok(()))
            .is_err()
    );
    for changed in [
        ImplicitBoundaryGeometry {
            area_m2: 2.0 * AREA,
            ..shape
        },
        ImplicitBoundaryGeometry {
            half_distance_m: 2.0 * DISTANCE,
            ..shape
        },
        ImplicitBoundaryGeometry {
            surface_resistance_m2_k_w: 0.1,
            ..shape
        },
    ] {
        assert!(
            face.state_at_specific_enthalpy(30_000.0, 300.0, changed, &mut || Ok(()))
                .is_err()
        );
    }
    assert!(matches!(
        face.state_at_specific_enthalpy(30_000.0, 300.0, shape, &mut || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));

    let original = curve(&card);
    let equal_but_separate = curve(&card);
    let other_layer = layer(&card, &equal_but_separate, &context);
    let foreign = ImplicitLayerBoundary::new(&other_layer, shape, config()).unwrap();
    let old = one_cell(&original, 500.0, baths(300.0, Some(300.0), false), vec![]);
    let law = transport(&card, &original, &context);
    let observed = observer(&old, &law, shape, false);
    assert!(
        ImplicitConductionChain::new(&old, &observed, &[(MaterialId(1), &law)], || Ok(()))
            .unwrap()
            .with_boundary_laws(&old, &[(CellId(1), &foreign)], &[], || Ok(()))
            .is_err()
    );
    assert_eq!(old.time_s(), 0.0);
}

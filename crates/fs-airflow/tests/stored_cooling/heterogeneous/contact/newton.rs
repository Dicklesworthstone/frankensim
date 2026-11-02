//! Real stored-solid/coolant solves checked against independent surface moments.
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::contact::newton::{
    ContactNewtonStoredCoolingControls, StoredContactRadiation,
};
use fs_conduction::radiation::{
    GrayDiffuseEnclosure, OpenDiffuseEnclosure, RadiationSurface, ViewFactorEvidence,
    ViewFactorMatrix, ViewFactorTolerance,
};
use fs_conduction::{PointwiseEnclosureTolerance, STEFAN_BOLTZMANN_W_M2_K4};

const NAMES: [&str; 2] = ["radiator-left", "radiator-right"];
const EPSILON: f64 = 0.9;
const DT: f64 = 0.2;

fn newton_controls() -> ContactNewtonStoredCoolingControls {
    let mut thermal = controls();
    thermal.solid.newton.relative_tolerance = 0.0;
    ContactNewtonStoredCoolingControls {
        thermal,
        max_residual_evaluations: 256,
        max_face_pairs: 2,
        residual_tolerance_j: 1e-8,
        max_total_enthalpy_solves: 20_000,
    }
}

fn smooth_contact() -> InterfaceSystemCard {
    let base = card(false, 450.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if let PropertyValue::Curve { knots, .. } = &mut claim.value {
            *knots = vec![(250.0, 0.002), (450.0, 0.02)];
        }
        claims.insert_claim(claim).unwrap();
    }
    InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

fn unequal_history(
    storage: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    charts: &[ResolvedHeatCapacityCurve; 2],
) -> Vec<f64> {
    storage
        .vertex_material_ids()
        .iter()
        .map(|&id| {
            charts[id]
                .state_at_temperature([340.0, 300.0][id])
                .unwrap()
                .specific_enthalpy_j_kg()
        })
        .collect()
}

fn radiating_boundary(
    mesh: &ConductionMesh,
    rows: &[CorrelationRobinBoundary],
) -> std::result::Result<ThermalBoundary, ConductionError> {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            NAMES[0],
            |face| face.outward_normal[2] < -0.99 && face.vertices.iter().all(|&v| v < 8),
            ThermalBc::adiabatic(),
        )?
        .region(
            NAMES[1],
            |face| face.outward_normal[2] < -0.99 && face.vertices.iter().all(|&v| v >= 8),
            ThermalBc::adiabatic(),
        )?
        .region(
            "body",
            |face| face.outward_normal[2] > 0.99,
            rows[0].boundary_condition().clone(),
        )?
        .adiabatic_remainder()
        .finish()
}

fn surfaces(mesh: &ConductionMesh) -> Vec<RadiationSurface> {
    let emissivity = radiation_patch().emissivity().clone();
    NAMES
        .into_iter()
        .enumerate()
        .map(|(body, name)| {
            RadiationSurface::new(
                mesh,
                name,
                |face| {
                    face.outward_normal[2] < -0.99
                        && face.vertices.iter().all(|&v| v as usize / 8 == body)
                },
                emissivity.clone(),
            )
            .unwrap()
        })
        .collect()
}

// Exact P1 fourth moment from simplex monomials, independent of the production
// face quadrature and radiosity factors: mean(T_h^4) = sum_(a+b+c=4) T^abc / 15.
fn emission(mesh: &ConductionMesh, temperature: &[f64], body: usize) -> (f64, f64) {
    let mut integral = 0.0;
    let mut area = 0.0;
    for face in mesh.boundary().iter().filter(|face| {
        face.outward_normal[2] < -0.99 && face.vertices.iter().all(|&v| v as usize / 8 == body)
    }) {
        let t = face.vertices.map(|v| temperature[v as usize]);
        let mut moment = 0.0;
        for a in 0..=4 {
            for b in 0..=4 - a {
                moment += t[0].powi(a) * t[1].powi(b) * t[2].powi(4 - a - b);
            }
        }
        integral += face.area * moment / 15.0;
        area += face.area;
    }
    (area, STEFAN_BOLTZMANN_W_M2_K4 * integral / area)
}

fn stored_change(
    state: &HeterogeneousStoredCoolingState<'_, '_, '_, '_>,
    old: &[f64],
    tank_h: f64,
) -> f64 {
    let solid: f64 = state
        .storage()
        .reference_nodal_masses_kg()
        .iter()
        .zip(state.solid_specific_enthalpy_j_kg())
        .zip(old)
        .map(|((&m, &h), &h0)| m * (h - h0))
        .sum();
    solid
        + state.reservoir().mass_kg()
            * (state.reservoir().state().specific_enthalpy_j_kg() - tank_h)
}

#[test]
fn simultaneous_contact_matches_the_existing_constant_resistance_cooling_step() {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = adjacent_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let source = card(false, 450.0);
        let contacts = contacts(cx, &mesh, &source);
        let old = unequal_history(&storage, &charts);
        let mut legacy = HeterogeneousStoredCoolingState::try_new(
            cx,
            uniform.reservoir().clone(),
            &storage,
            old.clone(),
        )
        .unwrap();
        let mut coupled = HeterogeneousStoredCoolingState::try_new(
            cx,
            uniform.reservoir().clone(),
            &storage,
            old,
        )
        .unwrap();
        let heat = left_heat(&mesh, 20.0);
        let problem = StoredSolidProblem {
            mesh: &mesh,
            material,
            source: &heat,
            element_materials: None,
            interfaces: None,
        };
        network
            .advance_stored_solid_with_pump_and_contacts(
                cx,
                &mut legacy,
                DT,
                0.0,
                &[320.0],
                pump,
                problem,
                &contacts,
                None,
                contact_controls(),
                |_, rows| boundary(&mesh, rows),
            )
            .unwrap();
        let result = network
            .advance_stored_solid_with_pump_and_contact_newton(
                cx,
                &mut coupled,
                DT,
                0.0,
                &[320.0],
                pump,
                problem,
                &contacts,
                StoredContactRadiation::None,
                newton_controls(),
                |_, rows| boundary(&mesh, rows),
            )
            .unwrap();
        assert!(
            legacy
                .solid_specific_enthalpy_j_kg()
                .iter()
                .zip(coupled.solid_specific_enthalpy_j_kg())
                .all(|(a, b)| (a - b).abs() < 2e-6)
        );
        assert!(
            (legacy.reservoir().state().temperature_k()
                - coupled.reservoir().state().temperature_k())
            .abs()
                < 1e-9
        );
        assert_eq!(
            result.total_enthalpy_solves,
            result.step.liquid.total_solid_solves
        );
        assert_eq!(
            result.total_contact_iterations,
            result.total_enthalpy_solves
        );
    });
}

#[test]
fn pointwise_ambient_radiation_contact_and_coolant_retain_separate_heat_accounts() {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = adjacent_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let source = smooth_contact();
        let contacts = contacts(cx, &mesh, &source);
        let old = unequal_history(&storage, &charts);
        let tank_h = uniform.reservoir().state().specific_enthalpy_j_kg();
        let mut state = HeterogeneousStoredCoolingState::try_new(
            cx,
            uniform.reservoir().clone(),
            &storage,
            old.clone(),
        )
        .unwrap();
        let heat = left_heat(&mesh, 20.0);
        let patches: Vec<_> = NAMES
            .into_iter()
            .map(|name| {
                AmbientRadiationPatch::new(name, radiation_patch().emissivity().clone(), 250.0)
                    .unwrap()
            })
            .collect();
        let result = network
            .advance_stored_solid_with_pump_and_contact_newton(
                cx,
                &mut state,
                DT,
                0.0,
                &[320.0],
                pump,
                StoredSolidProblem {
                    mesh: &mesh,
                    material,
                    source: &heat,
                    element_materials: None,
                    interfaces: None,
                },
                &contacts,
                StoredContactRadiation::Ambient(&patches),
                newton_controls(),
                |_, rows| radiating_boundary(&mesh, rows),
            )
            .unwrap();
        let response = &result.step.liquid.circuit.state.solid.state;
        let solid = &response.thermal.conduction;
        let expected: f64 = (0..2)
            .map(|body| {
                let (area, emitted) = emission(&mesh, &solid.temperature, body);
                area * EPSILON * (emitted - STEFAN_BOLTZMANN_W_M2_K4 * 250.0_f64.powi(4))
            })
            .sum();
        assert!((response.thermal.physical.radiation_out_w - expected).abs() < 1e-10);
        assert!(expected > 0.0);
        assert!(solid.pointwise_radiation.is_some());
        assert!(response.thermal.radiation.is_none());
        assert!(response.thermal.combined_boundary.is_none());
        assert_eq!(response.external_convection_out_w, 0.0);
        assert!((stored_change(&state, &old, tank_h) - DT * (20.0 - expected)).abs() < 2e-6);
        assert!((DT * response.coolant_out_w - result.step.liquid.returned_heat_j).abs() < 2e-6);
        assert_eq!(
            state.solid_specific_enthalpy_j_kg(),
            solid.specific_enthalpy_j_kg.as_slice()
        );
        for contact in response.thermal.physical_contacts.states() {
            let expected_r = 0.002 + (contact.interface_temperature_k - 250.0) * 0.018 / 200.0;
            assert!((contact.resistance.value_m2_k_per_w() - expected_r).abs() < 1e-14);
        }
        assert_eq!(
            result.total_enthalpy_solves,
            result.step.liquid.total_solid_solves
        );
    });
}

#[test]
fn enclosure_reflection_redistributes_heat_and_only_black_bath_exchange_is_external() {
    for open in [false, true] {
        fixture(1.0, |cx, _, network, pump, uniform, _, material| {
            let (mesh, ids) = adjacent_bodies();
            let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
            let storage = prepare(cx, &mesh, &ids, &charts);
            let source = smooth_contact();
            let contacts = contacts(cx, &mesh, &source);
            let old = unequal_history(&storage, &charts);
            let tank_h = uniform.reservoir().state().specific_enthalpy_j_kg();
            let mut state = HeterogeneousStoredCoolingState::try_new(
                cx,
                uniform.reservoir().clone(),
                &storage,
                old.clone(),
            )
            .unwrap();
            let heat = left_heat(&mesh, 20.0);
            // Prescribed equal-area algebra fixtures, not view factors inferred
            // from the displayed mesh separation or surface orientation.
            let closed = GrayDiffuseEnclosure::new(
                surfaces(&mesh),
                ViewFactorMatrix::infinite_parallel_plates(0.003).unwrap(),
            )
            .unwrap();
            let partial = OpenDiffuseEnclosure::new(
                cx,
                surfaces(&mesh),
                vec![vec![0.0, 0.25], vec![0.25, 0.0]],
                ViewFactorEvidence::Analytic {
                    geometry: "declared partial-view thermal fixture".into(),
                },
                ViewFactorTolerance::default(),
            )
            .unwrap();
            let radiation = if open {
                StoredContactRadiation::OpenEnclosure {
                    enclosure: &partial,
                    ambient_temperature_k: 420.0,
                    tolerance: PointwiseEnclosureTolerance::default(),
                }
            } else {
                StoredContactRadiation::ClosedEnclosure {
                    enclosure: &closed,
                    tolerance: PointwiseEnclosureTolerance::default(),
                }
            };
            let result = network
                .advance_stored_solid_with_pump_and_contact_newton(
                    cx,
                    &mut state,
                    DT,
                    0.0,
                    &[320.0],
                    pump,
                    StoredSolidProblem {
                        mesh: &mesh,
                        material,
                        source: &heat,
                        element_materials: None,
                        interfaces: None,
                    },
                    &contacts,
                    radiation,
                    newton_controls(),
                    |_, rows| radiating_boundary(&mesh, rows),
                )
                .unwrap();
            let response = &result.step.liquid.circuit.state.solid.state;
            let solid = &response.thermal.conduction;
            let report = solid.pointwise_enclosure_radiation.as_ref().unwrap();
            let (area, e0) = emission(&mesh, &solid.temperature, 0);
            let (_, e1) = emission(&mesh, &solid.temperature, 1);
            let view = if open { 0.25 } else { 1.0 };
            let eb = STEFAN_BOLTZMANN_W_M2_K4 * 420.0_f64.powi(4);
            let reflection = (1.0 - EPSILON) * view;
            let rhs = [e0, e1].map(|e| EPSILON * e + (1.0 - EPSILON) * (1.0 - view) * eb);
            // Independent two-by-two radiosity elimination.
            let j = [
                (rhs[0] + reflection * rhs[1]) / (1.0 - reflection * reflection),
                (rhs[1] + reflection * rhs[0]) / (1.0 - reflection * reflection),
            ];
            let q = [
                area * (j[0] - view * j[1] - (1.0 - view) * eb),
                area * (j[1] - view * j[0] - (1.0 - view) * eb),
            ];
            for (patch, expected) in report.patches.iter().zip(q) {
                assert!((patch.outward_heat_w - expected).abs() < 1e-10);
            }
            let external = q[0] + q[1];
            if open {
                assert!(external < 0.0); // Hot black bath supplies energy.
                assert!((report.ambient.as_ref().unwrap().outward_heat_w + external).abs() < 1e-10);
            } else {
                assert!(q[0] > 0.0 && q[1] < 0.0);
                assert!(external.abs() < 1e-10);
                assert!(report.ambient.is_none());
            }
            assert!((response.thermal.physical.radiation_out_w - external).abs() < 1e-10);
            assert!((stored_change(&state, &old, tank_h) - DT * (20.0 - external)).abs() < 2e-6);
            assert!((result.step.solid_external_heat_j - DT * (20.0 - external)).abs() < 1e-8);
            assert!(result.step.interface_energy_residual_j.abs() < 2e-6);
            assert_eq!(
                result.total_enthalpy_solves,
                result.step.liquid.total_solid_solves
            );
        });
    }
}

#[test]
fn newton_radiation_refusals_publish_neither_solid_nor_tank_history() {
    for case in 0..4 {
        fixture(1.0, |cx, gate, network, pump, uniform, _, material| {
            let (mesh, ids) = adjacent_bodies();
            let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
            let storage = prepare(cx, &mesh, &ids, &charts);
            let source = if case == 0 {
                card(false, 300.0)
            } else {
                smooth_contact()
            };
            let contacts = contacts(cx, &mesh, &source);
            let old = unequal_history(&storage, &charts);
            let tank = uniform.reservoir().clone();
            let mut state =
                HeterogeneousStoredCoolingState::try_new(cx, tank.clone(), &storage, old.clone())
                    .unwrap();
            let heat = left_heat(&mesh, 20.0);
            let mut config = newton_controls();
            if case == 1 {
                config.max_total_enthalpy_solves = 1;
            }
            let patches: Vec<_> = NAMES
                .into_iter()
                .map(|name| {
                    AmbientRadiationPatch::new(name, radiation_patch().emissivity().clone(), 250.0)
                        .unwrap()
                })
                .collect();
            let result = network.advance_stored_solid_with_pump_and_contact_newton(
                cx,
                &mut state,
                DT,
                0.0,
                &[320.0],
                pump,
                StoredSolidProblem {
                    mesh: &mesh,
                    material,
                    source: &heat,
                    element_materials: None,
                    interfaces: None,
                },
                &contacts,
                StoredContactRadiation::Ambient(&patches),
                config,
                |_, rows| {
                    if case == 2 {
                        gate.request();
                    }
                    if case == 3 {
                        boundary(&mesh, rows)
                    }
                    // Missing declared radiator traces.
                    else {
                        radiating_boundary(&mesh, rows)
                    }
                },
            );
            assert!(result.is_err(), "case {case} unexpectedly published");
            if case == 1 {
                assert!(matches!(
                    result,
                    Err(StoredCoolingError::ThermalWorkBudget { .. })
                ));
            }
            assert_eq!(state.reservoir(), &tank);
            assert_eq!(state.solid_specific_enthalpy_j_kg(), old.as_slice());
        });
    }
}

//! Actual contact/FEM/pump/tank composition, using the existing cooling fixture.
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::contact::ContactStoredCoolingControls;
use fs_conduction::{ThermalInterfaces, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, TEMPERATURE_DIMS};
use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContact, TemperatureContacts};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_blake3::ContentHash;

fn adjacent_bodies() -> (ConductionMesh, Vec<usize>) {
    let (a, mut positions) = fs_conduction::fixtures::box_grid([1; 3], [0.03, 0.1, 0.01]);
    let (b, other) = fs_conduction::fixtures::box_grid([1; 3], [0.03, 0.1, 0.01]);
    let offset = positions.len() as u32;
    positions.extend(other.into_iter().map(|[x, y, z]| [x + 0.03, y, z]));
    let mut ids = vec![0; a.tets.len()]; ids.extend(vec![1; b.tets.len()]);
    let mut tets = a.tets;
    tets.extend(b.tets.into_iter().map(|t| t.map(|v| v + offset)));
    (ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap(), ids)
}

fn card(varying: bool, upper: f64) -> InterfaceSystemCard {
    let knots = if upper == 300.0 { vec![(250.0, 0.002), (300.0, 0.002)] }
        else if varying { vec![(250.0, 0.002), (300.0, 0.002), (301.0, 0.02), (450.0, 0.02)] }
        else { vec![(250.0, 0.002), (450.0, 0.002)] };
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: TEMPERATURE_DIMS,
            knots, dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 250.0, upper).with("normal_pressure", 1e6, 1e6),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic contact curve, not measured data".into(),
            license: "test-only".into(), artifact: None },
    }).unwrap();
    let side = |name: &str| SurfaceSpec { material: MaterialStateId { chemistry: name.into(),
        phase: "solid".into(), process: "synthetic fixture".into(), revision: 0 }, texture_frame: "fixed normal".into() };
    InterfaceSystemCard::assemble(side("body A"), side("body B"), SystemContext {
        medium: "dry".into(), third_body: None, environment: "declared surroundings".into(), history: "fixed geometry".into(),
    }, claims, vec![]).unwrap()
}

fn contacts<'m, 'c>(cx: &Cx<'_>, mesh: &'m ConductionMesh, source: &'c InterfaceSystemCard)
    -> TemperatureContacts<'m, 'c>
{
    let base = ThermalBoundaryBuilder::new(mesh).adiabatic_remainder().finish().unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap().with("normal_pressure", 1e6).unwrap();
    let pairs = ThermalInterfaces::coincident_face_pairs(mesh).unwrap();
    assert_eq!(pairs.len(), 2);
    let law = TemperatureContact::new("joint", pairs, source, &point, "T", SelectionPolicy::SingleClaimOnly).unwrap();
    TemperatureContacts::new(cx, mesh, &base, vec![], vec![law], 2).unwrap()
}

fn contact_controls() -> ContactStoredCoolingControls {
    ContactStoredCoolingControls { thermal: controls(), contact: ContactIterationConfig {
        max_iterations: 64, relaxation: 0.5, residual_tolerance_j: 1e-8,
    }, max_total_enthalpy_solves: 200_000 }
}

fn left_heat(mesh: &ConductionMesh, watts: f64) -> ScalarField {
    let mut values = vec![0.0; mesh.vertex_count()];
    values[..8].fill(watts / 0.00003);
    ScalarField::Nodal(values)
}

#[test]
fn constant_source_contact_matches_the_existing_frozen_contact_transaction() {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = adjacent_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let source = card(false, 450.0); let contacts = contacts(cx, &mesh, &source);
        let base = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let frozen = contacts.evaluate(cx, &base, &vec![300.0; mesh.vertex_count()]).unwrap();
        let mut fixed = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(),
            &storage, initial(&storage, &charts)).unwrap();
        let mut sourced = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(),
            &storage, initial(&storage, &charts)).unwrap();
        let heat = left_heat(&mesh, 20.0);
        let problem = StoredSolidProblem { mesh: &mesh, material, source: &heat, element_materials: None, interfaces: None };
        let a = network.advance_stored_solid_with_pump(cx, &mut fixed, 1.0, 0.0, &[300.0], pump,
            StoredSolidProblem { interfaces: Some(frozen.frozen_interfaces()), ..problem },
            controls(), |_, rows| boundary(&mesh, rows)).unwrap();
        let b = network.advance_stored_solid_with_pump_and_contacts(cx, &mut sourced, 1.0, 0.0,
            &[300.0], pump, problem, &contacts, None, contact_controls(), |_, rows| boundary(&mesh, rows)).unwrap();
        assert!((fixed.reservoir().state().temperature_k() - sourced.reservoir().state().temperature_k()).abs() < 1e-8);
        assert!(a.liquid.circuit.state.solid.state.temperature.iter()
            .zip(&b.step.liquid.circuit.state.solid.state.thermal.conduction.temperature)
            .all(|(x, y)| (x-y).abs() < 1e-8));
        assert!(b.total_enthalpy_solves >= b.step.liquid.total_solid_solves);
        assert_eq!(b.total_enthalpy_solves, b.total_contact_iterations);
        assert!(b.step.combined_energy_residual_j.abs() < 2e-6);
    });
}

fn pulse(varying: bool) -> (f64, Vec<f64>) {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = adjacent_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let source = card(varying, 450.0); let contacts = contacts(cx, &mesh, &source);
        let old = initial(&storage, &charts);
        let tank_h = uniform.reservoir().state().specific_enthalpy_j_kg();
        let mut state = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(), &storage, old.clone()).unwrap();
        let mut walls = vec![300.0]; let mut temperatures = Vec::new();
        let mut pulse_tank = 0.0;
        for (end, watts) in [(1.0, 20.0), (2.0, 0.0)] {
            let heat = left_heat(&mesh, watts);
            let result = network.advance_stored_solid_with_pump_and_contacts(cx, &mut state, end, 0.0,
                &walls, pump, StoredSolidProblem { mesh: &mesh, material, source: &heat,
                    element_materials: None, interfaces: None }, &contacts, None,
                contact_controls(), |_, rows| boundary(&mesh, rows)).unwrap();
            let response = &result.step.liquid.circuit.state.solid.state;
            assert!(response.thermal.physical.residual_norm_j <= contact_controls().contact.residual_tolerance_j);
            assert!(response.thermal.physical.energy_residual_j.abs() <= controls().solid.energy_tolerance_j);
            assert_eq!(response.thermal.physical_contacts.states().len(), 2);
            for local in response.thermal.physical_contacts.states() {
                assert!(local.interface_temperature_k > 300.0);
                if varying { assert!(local.resistance.value_m2_k_per_w() > 0.002); }
                else { assert!((local.resistance.value_m2_k_per_w()-0.002).abs() < 1e-14); }
            }
            assert_eq!(state.solid_specific_enthalpy_j_kg(), response.thermal.conduction.specific_enthalpy_j_kg.as_slice());
            let solid_j: f64 = storage.reference_nodal_masses_kg().iter().zip(&old)
                .zip(state.solid_specific_enthalpy_j_kg()).map(|((&m, &h0), &h)| m*(h-h0)).sum();
            let tank_j = state.reservoir().mass_kg() * (state.reservoir().state().specific_enthalpy_j_kg()-tank_h);
            assert!((solid_j + tank_j - 20.0).abs() < 4e-6);
            assert!(result.step.interface_energy_residual_j.abs() < 2e-6);
            if end == 1.0 { pulse_tank = state.reservoir().state().temperature_k(); }
            else { assert!(state.reservoir().state().temperature_k() > pulse_tank); }
            temperatures = response.thermal.conduction.temperature.clone();
            walls = result.step.liquid.circuit.state.solid.regions.iter().map(|r| r.mean_wall_temperature_k).collect();
        }
        (state.reservoir().state().temperature_k(), temperatures)
    })
}

#[test]
fn contact_material_alone_changes_the_coupled_history_without_creating_heat() {
    let constant = pulse(false); let varying = pulse(true);
    assert!((constant.0-varying.0).abs() > 1e-10);
    assert!(constant.1.iter().zip(&varying.1).any(|(a,b)| (a-b).abs() > 1e-7));
}

#[test]
fn source_contacts_and_radiation_only_surfaces_share_one_complete_cooling_step() {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = adjacent_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let source = card(true, 450.0); let contacts = contacts(cx, &mesh, &source);
        let mut state = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(),
            &storage, initial(&storage, &charts)).unwrap();
        let heat = left_heat(&mesh, 20.0); let patches = [radiation_patch()];
        let rad = AmbientRadiationConfig { max_iterations: 48, relaxation: 1.0,
            temperature_tolerance_k: 1e-9, balance_tolerance_w: 1e-8, balance_relative_tolerance: 1e-10 };
        let result = network.advance_stored_solid_with_pump_and_contacts(cx, &mut state, 1.0, 0.0,
            &[300.0], pump, StoredSolidProblem { mesh: &mesh, material, source: &heat,
                element_materials: None, interfaces: None }, &contacts, Some((&patches, rad)),
            contact_controls(), |_, rows| ThermalBoundaryBuilder::new(&mesh)
                .region("exposed", |f| f.outward_normal[2] < -0.99, ThermalBc::adiabatic())?
                .region("body", |f| f.outward_normal[2] > 0.99, rows[0].boundary_condition().clone())?
                .adiabatic_remainder().finish()).unwrap();
        let response = &result.step.liquid.circuit.state.solid.state;
        assert_eq!(response.external_convection_out_w, 0.0);
        assert_eq!(response.coolant_robin_fluxes.len(), 1);
        assert_eq!(response.coolant_robin_fluxes[0].region, "body");
        assert!(response.thermal.physical.radiation_out_w > 0.0);
        assert!(response.thermal.radiation.is_some());
        assert!(response.thermal.combined_boundary.is_some());
        assert!((response.coolant_out_w - result.step.liquid.returned_heat_j).abs() < 2e-6);
        assert!((result.step.solid_external_heat_j - (20.0-response.thermal.physical.radiation_out_w)).abs() < 1e-8);
        assert!(result.step.combined_energy_residual_j.abs() < 2e-6);
        assert!(result.total_enthalpy_solves >= result.total_contact_iterations);
        assert!(result.total_contact_iterations >= result.step.liquid.total_solid_solves);
    });
}

#[test]
fn contact_heat_can_change_latent_fraction_without_reconstructing_enthalpy_from_temperature() {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = adjacent_bodies();
        let sensible = chart_with_origin(1.0, 0.0);
        let latent = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([174; 32]),
            [(0.0,250.0,0.0), (25000.0,300.0,0.0), (125000.0,300.0,1.0), (200000.0,450.0,1.0)]
                .into_iter().map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
                    temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: 1000.0 }).collect()).unwrap();
        let storage = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &[
            ReferenceEnthalpyMaterial { curve: sensible.curve(), reference_density_kg_m3: 1000.0 },
            ReferenceEnthalpyMaterial { curve: &latent, reference_density_kg_m3: 1000.0 },
        ], &ids, EnthalpyBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count() }).unwrap();
        let h0 = sensible.state_at_temperature(300.0).unwrap().specific_enthalpy_j_kg();
        let old: Vec<_> = storage.vertex_material_ids().iter().map(|&id| if id == 0 { h0 } else { 75000.0 }).collect();
        let source = card(true, 450.0); let contacts = contacts(cx, &mesh, &source);
        let mut state = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(), &storage, old).unwrap();
        let heat = left_heat(&mesh, 20.0);
        let result = network.advance_stored_solid_with_pump_and_contacts(cx, &mut state, 1.0, 0.0,
            &[300.0], pump, StoredSolidProblem { mesh: &mesh, material, source: &heat,
                element_materials: None, interfaces: None }, &contacts, None,
            contact_controls(), |_, rows| boundary(&mesh, rows)).unwrap();
        let solid = &result.step.liquid.circuit.state.solid.state.thermal.conduction;
        assert!(solid.temperature[8..].iter().all(|t| (*t-300.0).abs() < 1e-10));
        assert!(solid.liquid_mass_fraction[8..].iter().any(|f| *f > 0.5));
        assert!(solid.specific_enthalpy_j_kg[8..].iter().any(|h| *h > 75000.0));
        assert_eq!(state.solid_specific_enthalpy_j_kg(), solid.specific_enthalpy_j_kg.as_slice());
        assert!((solid.stored_energy_change_j + result.step.liquid.stored_enthalpy_change_j - 20.0).abs() < 2e-6);
    });
}

#[test]
fn contact_source_work_ownership_and_cancel_failures_publish_neither_history() {
    for case in 0..5 {
        fixture(1.0, |cx, gate, network, pump, uniform, _, material| {
            let (mesh, ids) = adjacent_bodies();
            let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
            let storage = prepare(cx, &mesh, &ids, &charts);
            let source = card(true, if case == 0 { 300.0 } else { 450.0 });
            let contacts = contacts(cx, &mesh, &source);
            let old = initial(&storage, &charts); let tank = uniform.reservoir().clone();
            let mut state = HeterogeneousStoredCoolingState::try_new(cx, tank.clone(), &storage, old.clone()).unwrap();
            let heat = left_heat(&mesh, 20.0);
            let base = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
            let frozen = contacts.evaluate(cx, &base, &vec![300.0; mesh.vertex_count()]).unwrap();
            let mut limits = contact_controls(); if case == 1 { limits.max_total_enthalpy_solves = 1; }
            let result = network.advance_stored_solid_with_pump_and_contacts(cx, &mut state, 1.0, 0.0,
                &[300.0], pump, StoredSolidProblem { mesh: &mesh, material, source: &heat, element_materials: None,
                    interfaces: if case == 2 { Some(frozen.frozen_interfaces()) } else { None } },
                &contacts, None, limits, |_, rows| {
                    if case == 3 { gate.request(); }
                    if case == 4 {
                        ThermalBoundaryBuilder::new(&mesh)
                            .region("body", |f| f.outward_normal[2] > 0.99, rows[0].boundary_condition().clone())?
                            .region("unowned", |f| f.outward_normal[2] < -0.99, ThermalBc::robin(5.0,280.0)?)?
                            .adiabatic_remainder().finish()
                    } else { boundary(&mesh, rows) }
                });
            assert!(result.is_err(), "case {case} unexpectedly published a step");
            if case == 2 || case == 4 { assert!(matches!(result, Err(StoredCoolingError::Invalid(_)))); }
            assert_eq!(state.reservoir(), &tank);
            assert_eq!(state.solid_specific_enthalpy_j_kg(), old.as_slice());
        });
    }
}

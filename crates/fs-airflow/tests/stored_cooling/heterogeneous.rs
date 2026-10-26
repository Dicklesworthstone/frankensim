//! Exercise the actual per-material FEM through the original coolant transaction.
#[path = "heterogeneous/contact.rs"]
mod contact;
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::stored_solid::{
    HeterogeneousStoredCoolingState, StoredMaterialCoolingState,
    radiation::RadiatingStoredCoolingControls,
};
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, SurfaceEmissivity,
    ThermalBc, SURFACE_EMISSIVITY_PROPERTY};
use fs_conduction::transient::enthalpy::{EnthalpyError,
    heterogeneous::{HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial}};
use fs_rep_mesh::TetComplex;

// Two separated components have the original fixture's TOTAL volume/top area.
// Distinct materials never share an enthalpy DOF or an invented mixture chart.
fn two_bodies() -> (ConductionMesh, Vec<usize>) {
    let (left, mut positions) = fs_conduction::fixtures::box_grid([1; 3], [0.03, 0.1, 0.01]);
    let (right, other) = fs_conduction::fixtures::box_grid([1; 3], [0.03, 0.1, 0.01]);
    let offset = positions.len() as u32;
    positions.extend(other.into_iter().map(|[x, y, z]| [x + 0.04, y, z]));
    let mut ids = vec![0; left.tets.len()];
    ids.extend(vec![1; right.tets.len()]);
    let mut tets = left.tets;
    tets.extend(right.tets.into_iter().map(|t| t.map(|v| v + offset)));
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    (mesh, ids)
}

fn chart_with_origin(scale: f64, origin: f64) -> ResolvedHeatCapacityCurve {
    let card = solid_card(scale);
    ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &key("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
        density: &key("density", Dims([-3, 1, 0, 0, 0, 0])), temperature_axis: "T",
        lower: &point(250.0), upper: &point(450.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, origin, 0.01).unwrap()
}

fn prepare<'m, 'c>(cx: &Cx<'_>, mesh: &'m ConductionMesh, ids: &[usize],
    charts: &'c [ResolvedHeatCapacityCurve; 2]) -> HeterogeneousEnthalpyBackwardEuler<'m, 'c>
{
    let materials = charts.iter().map(|chart| ReferenceEnthalpyMaterial {
        curve: chart.curve(), reference_density_kg_m3: 1000.0,
    }).collect::<Vec<_>>();
    HeterogeneousEnthalpyBackwardEuler::new(cx, mesh, &materials, ids,
        EnthalpyBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count() }).unwrap()
}

fn initial(storage: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    charts: &[ResolvedHeatCapacityCurve; 2]) -> Vec<f64>
{
    storage.vertex_material_ids().iter().map(|&id|
        charts[id].state_at_temperature(300.0).unwrap().specific_enthalpy_j_kg()).collect()
}

#[test]
fn a_single_assigned_material_reproduces_the_original_uniform_transaction() {
    fixture(1.0, |cx, _, network, pump, uniform, mesh, material| {
        let storage = HeterogeneousEnthalpyBackwardEuler::new(cx, mesh,
            &[ReferenceEnthalpyMaterial { curve: uniform.storage().phase_curve(),
                reference_density_kg_m3: 1000.0 }], &vec![0; mesh.element_count()],
            EnthalpyBudget { max_vertices: mesh.vertex_count(), max_elements: mesh.element_count() }).unwrap();
        let mut assigned = HeterogeneousStoredCoolingState::try_new(cx,
            uniform.reservoir().clone(), &storage, uniform.solid_specific_enthalpy_j_kg().to_vec()).unwrap();
        let heat = ScalarField::uniform("20 W", 20.0 / 0.00006).unwrap();
        let problem = StoredSolidProblem { mesh, material, source: &heat,
            element_materials: None, interfaces: None };
        let a = network.advance_stored_solid_with_pump(cx, uniform, 1.0, 0.0, &[300.0],
            pump, problem, controls(), |_, rows| boundary(mesh, rows)).unwrap();
        let b = network.advance_stored_solid_with_pump(cx, &mut assigned, 1.0, 0.0, &[300.0],
            pump, problem, controls(), |_, rows| boundary(mesh, rows)).unwrap();
        assert_eq!(uniform.reservoir(), assigned.reservoir());
        assert_eq!(uniform.solid_specific_enthalpy_j_kg(), assigned.solid_specific_enthalpy_j_kg());
        assert_eq!(a.liquid.circuit.state.solid.state.temperature, b.liquid.circuit.state.solid.state.temperature);
        assert_eq!(a.combined_energy_residual_j, b.combined_energy_residual_j);
    });
}

#[test]
fn initial_history_is_checked_against_each_assigned_chart_not_the_first_material() {
    fixture(1.0, |cx, _, _, _, uniform, _, _| {
        let (mesh, ids) = two_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let mut h = initial(&storage, &charts);
        let state = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(), &storage, h.clone()).unwrap();
        assert!(std::ptr::eq(state.storage(), &storage));
        assert_eq!(state.storage().vertex_material_ids(), storage.vertex_material_ids());
        let vertex = storage.vertex_material_ids().iter().position(|&id| id == 1).unwrap();
        // This value is perfectly valid for material zero, but outside the
        // second material's deliberately disjoint specific-enthalpy interval.
        h[vertex] = h[0];
        let refused = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(), &storage, h);
        assert!(matches!(refused, Err(StoredCoolingError::Solid(EnthalpyError::Phase { vertex: v, .. })) if v == vertex));
    });
}

fn pulse(scale: f64, origin: f64) -> (f64, Vec<f64>) {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = two_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(scale, origin)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let old = initial(&storage, &charts);
        let tank_h = uniform.reservoir().state().specific_enthalpy_j_kg();
        let mut state = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(), &storage, old.clone()).unwrap();
        let mut walls = vec![300.0];
        let mut pulse_temperature = 0.0;
        let mut final_temperature = Vec::new();
        for (end, watts) in [(1.0, 20.0), (2.0, 0.0)] {
            let heat = ScalarField::uniform("pulse", watts / 0.00006).unwrap();
            let result = network.advance_stored_solid_with_pump(cx, &mut state, end, 0.0, &walls, pump,
                StoredSolidProblem { mesh: &mesh, material, source: &heat, element_materials: None, interfaces: None },
                controls(), |_, rows| boundary(&mesh, rows)).unwrap();
            let solid = &result.liquid.circuit.state.solid.state;
            assert!(result.interface_energy_residual_j.abs() <= controls().joint_energy_tolerance_j);
            assert!(result.combined_energy_residual_j.abs() <= controls().joint_energy_tolerance_j);
            assert_eq!(state.reservoir().time_s(), end);
            assert_eq!(state.solid_specific_enthalpy_j_kg(), solid.specific_enthalpy_j_kg.as_slice());
            for (v, &h) in state.solid_specific_enthalpy_j_kg().iter().enumerate() {
                let local = storage.phase_curve_at(v).unwrap().state_at_specific_enthalpy(h).unwrap();
                assert_eq!(local.temperature_k(), solid.temperature[v]);
            }
            let solid_j: f64 = storage.reference_nodal_masses_kg().iter().zip(&old)
                .zip(state.solid_specific_enthalpy_j_kg()).map(|((&m, &h0), &h)| m * (h - h0)).sum();
            let tank_j = state.reservoir().mass_kg() * (state.reservoir().state().specific_enthalpy_j_kg() - tank_h);
            assert!((solid_j + tank_j - 20.0).abs() < 4e-6);
            if end == 1.0 { pulse_temperature = state.reservoir().state().temperature_k(); }
            else {
                assert!(state.reservoir().state().temperature_k() > pulse_temperature);
                assert!(solid.stored_energy_change_j < 0.0);
            }
            walls = result.liquid.circuit.state.solid.regions.iter().map(|r| r.mean_wall_temperature_k).collect();
            final_temperature = solid.temperature.clone();
        }
        (state.reservoir().state().temperature_k(), final_temperature)
    })
}

#[test]
fn material_substitution_changes_cooling_but_enthalpy_origins_do_not() {
    let a = pulse(2.0, 0.0);
    let shifted = pulse(2.0, 1e6);
    assert!((a.0 - shifted.0).abs() < 2e-7);
    assert!(a.1.iter().zip(&shifted.1).all(|(x, y)| (x - y).abs() < 2e-7));
    let larger_capacity = pulse(4.0, 0.0);
    assert!(larger_capacity.0 < a.0);
    // The second component's only changed input is its caloric source data.
    assert!(larger_capacity.1[8..].iter().sum::<f64>() < a.1[8..].iter().sum::<f64>());
    assert!(a.1[..8].iter().sum::<f64>() > a.1[8..].iter().sum::<f64>());
}

fn radiation_patch() -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::with_quantity(SURFACE_EMISSIVITY_PROPERTY, QuantitySpec::dimensional(Dims::NONE)),
        value: PropertyValue::Scalar { value: 0.9, dims: Dims::NONE },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 450.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "synthetic finish".into(), license: "test-only".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "test finish".into(), phase: "solid".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    AmbientRadiationPatch::new("exposed", SurfaceEmissivity::from_card("exposed", &card,
        300.0, SelectionPolicy::SingleClaimOnly).unwrap(), 250.0).unwrap()
}

#[test]
fn heterogeneous_radiation_retains_original_charts_and_separate_heat_destinations() {
    fixture(1.0, |cx, _, network, pump, uniform, _, material| {
        let (mesh, ids) = two_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let mut state = HeterogeneousStoredCoolingState::try_new(cx, uniform.reservoir().clone(), &storage, initial(&storage, &charts)).unwrap();
        let heat = ScalarField::uniform("20 W", 20.0 / 0.00006).unwrap();
        let result = network.advance_stored_solid_with_pump_and_radiation(cx, &mut state,
            1.0, 0.0, &[300.0], pump, StoredSolidProblem { mesh: &mesh, material, source: &heat,
                element_materials: None, interfaces: None }, &[radiation_patch()],
            RadiatingStoredCoolingControls { thermal: controls(), radiation: AmbientRadiationConfig {
                max_iterations: 48, relaxation: 1.0, temperature_tolerance_k: 1e-9,
                balance_tolerance_w: 1e-8, balance_relative_tolerance: 1e-10 }, max_total_enthalpy_solves: 65536 },
            |_, rows| ThermalBoundaryBuilder::new(&mesh)
                .region("exposed", |f| f.outward_normal[2] < -0.99, ThermalBc::adiabatic())?
                .region("body", |f| f.outward_normal[2] > 0.99, rows[0].boundary_condition().clone())?
                .adiabatic_remainder().finish()).unwrap();
        let response = &result.step.liquid.circuit.state.solid.state;
        assert_eq!(response.external_convection_out_w, 0.0);
        assert_eq!(response.coolant_robin_fluxes.len(), 1);
        assert_eq!(response.thermal.convective_robin_fluxes.len(), 1);
        assert!(response.thermal.radiation.nonlinear_radiation_out_w > 0.0);
        assert!(result.step.combined_energy_residual_j.abs() < 2e-6);
        assert!(result.step.interface_energy_residual_j.abs() < 2e-6);
        assert!(response.thermal.physical_residual_norm_j <= response.thermal.physical_residual_tolerance_j);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), response.thermal.conduction.specific_enthalpy_j_kg.as_slice());
    });
}

#[test]
fn heterogeneous_failures_leave_all_materials_and_the_reservoir_unpublished() {
    fixture(1.0, |cx, gate, network, pump, uniform, _, material| {
        let (mesh, ids) = two_bodies();
        let charts = [chart_with_origin(1.0, 0.0), chart_with_origin(2.0, 1e6)];
        let storage = prepare(cx, &mesh, &ids, &charts);
        let old = initial(&storage, &charts);
        let tank = uniform.reservoir().clone();
        let mut state = StoredMaterialCoolingState::try_new(cx, tank.clone(), &storage, old.clone()).unwrap();
        let heat = ScalarField::uniform("20 W", 20.0 / 0.00006).unwrap();
        let problem = StoredSolidProblem { mesh: &mesh, material, source: &heat, element_materials: None, interfaces: None };
        let mut budget = controls(); budget.pump.max_total_solid_solves = 1;
        assert!(network.advance_stored_solid_with_pump(cx, &mut state, 1.0, 0.0, &[300.0], pump,
            problem, budget, |_, rows| boundary(&mesh, rows)).is_err());
        assert_eq!(state.reservoir(), &tank); assert_eq!(state.solid_specific_enthalpy_j_kg(), old.as_slice());
        assert!(network.advance_stored_solid_with_pump(cx, &mut state, 1.0, 0.0, &[300.0], pump,
            problem, controls(), |_, rows| { gate.request(); boundary(&mesh, rows) }).is_err());
        assert_eq!(state.reservoir(), &tank); assert_eq!(state.solid_specific_enthalpy_j_kg(), old.as_slice());
    });
}

//! Manufactured finite core/winding exchange; no measured magnetic fit claim.
use super::*;
use crate::circuit_enthalpy::{
    contacts::{heat::{MaterialHeatNetwork, PreparedHeatNetwork}, material::ElectricalContactNetwork},
    material::ConductionLink,
};
use fs_conduction::{
    GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity, ViewFactorMatrix,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
};
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_exec::{CancelGate, Cx};
use fs_material::{state_point::THERMAL_CONDUCTIVITY_DIMS, thermal::interface::ThermalInterfacePoint};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};

const CORE: CellId = CellId(10);
const COIL: CellId = CellId(11);

fn with_heat<R>(emissivity: f64, resistance: f64, core_upper: f64,
    run: impl FnOnce(&Cx<'_>, &PreparedHeatNetwork<'_, '_>,
        &DepositionSimulation<'_>, &ThermalHystereticCircuit<'_, '_>) -> R) -> R {
    let (original, _) = material(1.0, 0.005);
    let mut claims = original.claims().clone();
    for (name, dims, value) in [("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, emissivity)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "manufactured core/winding exchange".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(200.0), upper: &point(2000.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    let mut events = [CORE, COIL].into_iter().enumerate().map(|(i, id)| TimedEvent {
        time_s: 0.0, event: Event::Birth(DepositedCell {
            id, segment_id: i as u64, body_id: i as u64, material_id: MaterialId(0), mass_kg: 1.0,
            specific_h_j_kg: if id == CORE { 90.0 } else { 100.0 },
            geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 1.0, 0.0, 0.0],
                volume_m3: 1.0, surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }) }).collect::<Vec<_>>();
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 0, a: CORE, b: COIL, area_m2: 1.0, conductance_w_k: 999.0 }) });
    let thermal = DepositionSimulation::new(events, &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 100 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.01, energy_tolerance_j: 1e-9 }, || Ok(())).unwrap().0;
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [CORE, COIL].map(|cell| MaterialCell { cell, material_id: MaterialId(0), curve: curve.curve(),
        point: ThermalMaterialPoint { material: &card, context: &context, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)) });
    let resistors = [CellResistor { destination: ResistorCell { branch: 1, cell: COIL },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 1.0, area_m2: 1.0 }];
    let network = MaterialNetwork { cells: &cells, resistors: &resistors,
        links: &[ConductionLink { link: 0, distance_a_m: 0.5, distance_b_m: 0.5 }] };
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let sources = network.resolve_with_bases(&thermal, &bases, || false).unwrap().resistors;
    let core = core(true);
    let start = initial(&core, 0.0, vec![5.0], vec![]);
    let circuit = HystereticCircuit::new(rl_graph(), 2, core, sources, start).unwrap();
    let driver = ThermalHystereticCircuit::bind(circuit, &thermal, &network, &bases, CORE,
        [280.0, core_upper], || false).unwrap();
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(RP, RD),
        value: PropertyValue::Scalar { value: resistance, dims: RD },
        validity: ValidityDomain::unconstrained().with("aT", 200.0, 2000.0).with("bT", 200.0, 2000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "manufactured core/coil contact".into(),
            license: "CC0-1.0".into(), artifact: None } }).unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec { material: card.id().clone(), texture_frame: "core".into() },
        SurfaceSpec { material: card.id().clone(), texture_frame: "coil".into() },
        SystemContext { medium: "fixed".into(), third_body: None, environment: "fixture".into(), history: "fixed".into() },
        claims, vec![]).unwrap();
    let ip = QueryPoint::new().with("aT", 290.0).unwrap().with("bT", 300.0).unwrap();
    let contacts = [MaterialContactLink { link: 0, sides: [CORE, COIL].map(|cell|
        ContactCell { cell, material_id: MaterialId(0), curve: curve.curve() }),
        point: ThermalInterfacePoint { interface: &interface, materials: [&card; 2], context: &ip,
            temperature_axes: ["aT", "bT"], liquid_fraction_axes: [None; 2], mean_temperature_axis: None },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
    let electrical = ElectricalContactNetwork { bulk: &network, contacts: &[], electrical_bases: &bases };
    let (complex, positions) = fs_conduction::fixtures::box_grid([1; 3], [1.0; 3]);
    let mesh = fs_conduction::mesh::ConductionMesh::new(complex, positions).unwrap();
    let surfaces = [0, 1].map(|i| RadiationSurface::new(&mesh, format!("face-{i}"),
        |face| fs_conduction::fixtures::on_box_face(face.centroid[0], i as f64),
        SurfaceEmissivity::from_thermal_state("fixture", cells[i].point,
            curve.curve().state_at_specific_enthalpy(100.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly).unwrap().0).unwrap()).to_vec();
    // Cube faces supply area; the view-factor fixture is the parallel-plate limit.
    let enclosure = GrayDiffuseEnclosure::new(surfaces, ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap()).unwrap();
    let surfaces = [CORE, COIL].into_iter().enumerate().map(|(surface, cell)| RadiatingCell {
        surface, cell, material_id: MaterialId(0), curve: curve.curve(), material: cells[surface].point,
        selection: MaterialPropertySelection::SingleClaimOnly }).collect::<Vec<_>>();
    let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
        power_tolerance_w: 1e-9, radiosity_tolerance_w_m2: 1e-9 };
    let assembly = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &contacts, radiation: Some(&radiation) };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey { seed: 43, kernel_id: 1, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        let prepared = assembly.prepare(&cx, &thermal).unwrap();
        run(&cx, &prepared, &thermal, &driver)
    })
}

#[test]
fn g1_g3_contact_and_radiation_change_actual_current_and_hysteresis_heat() {
    let run = |emissivity, resistance| with_heat(emissivity, resistance, 330.0,
        |cx, prepared, thermal, driver| {
            let next = driver.step_candidate_with_heat_network(cx, prepared, thermal, 0.05,
                &[5.0], heat_controls(), controls(), &[0.0; 2]).unwrap();
            let [tc, tw] = [CORE, COIL].map(|id| next.thermal().temperature_k(id).unwrap());
            let heat = 0.05 * ((tw - tc) / (1.0 + resistance)
                + fs_conduction::STEFAN_BOLTZMANN_W_M2_K4 * (tw.powi(4) - tc.powi(4)) / (2.0 / emissivity - 1.0));
            near(tc - thermal.temperature_k(CORE).unwrap(), next.report().electrical.core.core_heat_j + heat, 3e-8);
            near(tw - thermal.temperature_k(COIL).unwrap(), next.report().electrical.circuit.heat_j - heat, 3e-8);
            near(next.report().deposited_joule_j,
                next.report().electrical.circuit.heat_j + next.report().electrical.core.core_heat_j, 1e-10);
            near(next.report().total_energy_residual_j, 0.0, 2e-8);
            assert_eq!(next.report().source.thermal_contacts.len(), 1);
            assert!(next.report().source.radiation.is_some());
            assert!(next.report().electrical.core.core_heat_j > 0.0);
            // A second accepted interval must retain magnetic memory, not restart
            // the original loop when the heat assembly is reused.
            let old_state = next.circuit().circuit().clone();
            let continued = next.circuit().step_candidate_with_heat_network(cx, prepared,
                next.thermal(), 0.1, &[5.0], heat_controls(), controls(), &[0.0; 2]).unwrap();
            assert_eq!(next.circuit().circuit(), &old_state);
            assert!(continued.circuit().circuit().flux_linkage_wb() > old_state.flux_linkage_wb());
            (continued.thermal().temperature_k(COIL).unwrap(), continued.circuit().circuit().current_a(),
                continued.report().electrical.core.core_heat_j)
        });
    let weak = run(0.1, 2.0);
    let contact = run(0.1, 0.001);
    let both = run(0.9, 0.001);
    assert!(both.0 < contact.0 && contact.0 < weak.0);
    assert!(both.1 > contact.1 && contact.1 > weak.1);
    assert!((both.2 - weak.2).abs() > 1e-10);
}

#[test]
fn g4_late_core_domain_refusal_preserves_memory_heat_and_retry() {
    with_heat(0.8, 0.001, 290.01, |cx, prepared, thermal, driver| {
        let before = driver.circuit().clone();
        let cells = thermal.cells().clone();
        assert!(driver.step_candidate_with_heat_network(cx, prepared, thermal, 0.05,
            &[5.0], heat_controls(), controls(), &[0.0; 2]).is_err());
        assert_eq!(driver.circuit(), &before);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
        let retry = || driver.step_candidate_with_heat_network(cx, prepared, thermal, 1e-6,
            &[5.0], heat_controls(), controls(), &[0.0; 2]).unwrap();
        let a = retry(); let b = retry();
        assert_eq!(a.thermal().cells(), b.thermal().cells());
        assert_eq!(a.circuit().circuit(), b.circuit().circuit());
    });
}

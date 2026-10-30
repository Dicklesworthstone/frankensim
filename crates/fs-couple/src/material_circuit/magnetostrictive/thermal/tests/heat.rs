//! Manufactured actuator/finite-sink heat exchange, not measured core data.
use super::*;
use crate::circuit_enthalpy::{contacts::{heat::MaterialHeatNetwork,
    material::ElectricalContactNetwork}, material::ConductionLink};
use fs_conduction::{GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity,
    ViewFactorMatrix, EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_exec::{CancelGate, Cx};
use fs_material::{state_point::THERMAL_CONDUCTIVITY_DIMS, thermal::interface::ThermalInterfacePoint};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};

const SINK: CellId = CellId(42);
const AREA: f64 = 0.01;

fn with_heat<R>(emissivity: f64, resistance: f64, upper: f64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &PreparedHeatNetwork<'_, '_>,
        &DepositionSimulation<'_>, &ThermalMagnetostrictiveCircuit<'_, '_>, f64) -> R) -> R {
    let (original, _) = winding(0.01, ElectricalPropertyBasis::Resistivity, upper);
    let mut claims = original.claims().clone();
    for (name, dims, value) in [("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, emissivity)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 600.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "manufactured thermal actuator".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    let material = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&material, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(290.0), upper: &point(600.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    let cell = initial(&curve, 0.0).cells()[&COIL].clone();
    let thermal = DepositionSimulation::new(vec![
        TimedEvent { time_s: 0.0, event: Event::Birth(cell.clone()) },
        TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell { id: SINK, ..cell }) },
        TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink { id: 0, a: COIL, b: SINK,
            area_m2: AREA, conductance_w_k: 999.0 }) },
    ], &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 2, max_links: 1, max_events: 3, max_substeps: 100,
    }).unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.01, energy_tolerance_j: 1e-10 },
        || Ok(())).unwrap().0;
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let (_, resistors) = bindings(&material, &curve, &context);
    let cells = [COIL, SINK].map(|cell| MaterialCell { cell, material_id: MaterialId(0),
        curve: curve.curve(), point: ThermalMaterialPoint { material: &material, context: &context,
            temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)) });
    let network = MaterialNetwork { cells: &cells, resistors: &resistors,
        links: &[ConductionLink { link: 0, distance_a_m: 0.05, distance_b_m: 0.05 }] };
    let rod = rod(0.2, 10.0);
    let voltage = RESISTANCE * rod.bias_current_a() + 0.2;
    let driver = ThermalMagnetostrictiveCircuit::new(rod, 0.0, 0.0, 0.0, voltage, &thermal,
        &network, ElectricalPropertyBasis::Resistivity, || false).unwrap();
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(RP, RD),
        value: PropertyValue::Scalar { value: resistance, dims: RD },
        validity: ValidityDomain::unconstrained().with("aT", 290.0, 600.0).with("bT", 290.0, 600.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: Provenance { source: "manufactured thermal contact".into(),
            license: "CC0-1.0".into(), artifact: None } }).unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec { material: material.id().clone(), texture_frame: "coil".into() },
        SurfaceSpec { material: material.id().clone(), texture_frame: "sink".into() },
        SystemContext { medium: "fixed".into(), third_body: None, environment: "fixture".into(), history: "fixed".into() },
        claims, vec![]).unwrap();
    let ip = QueryPoint::new().with("aT", 300.0).unwrap().with("bT", 300.0).unwrap();
    let contacts = [MaterialContactLink { link: 0, sides: [COIL, SINK].map(|cell|
        ContactCell { cell, material_id: MaterialId(0), curve: curve.curve() }),
        point: ThermalInterfacePoint { interface: &interface, materials: [&material; 2], context: &ip,
            temperature_axes: ["aT", "bT"], liquid_fraction_axes: [None; 2], mean_temperature_axis: None },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
    let electrical = ElectricalContactNetwork { bulk: &network, contacts: &[],
        electrical_bases: &[(0, ElectricalPropertyBasis::Resistivity)] };
    let (complex, positions) = fs_conduction::fixtures::box_grid([1; 3], [0.1; 3]);
    let mesh = fs_conduction::mesh::ConductionMesh::new(complex, positions).unwrap();
    let surfaces = [0, 1].map(|i| RadiationSurface::new(&mesh, format!("face-{i}"),
        |face| fs_conduction::fixtures::on_box_face(face.centroid[0], 0.1 * i as f64),
        SurfaceEmissivity::from_thermal_state("fixture", cells[i].point,
            curve.curve().state_at_specific_enthalpy(100.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly).unwrap().0).unwrap()).to_vec();
    // Algebraic parallel-plate limit; cube faces carry area, not finite-cube view factors.
    let enclosure = GrayDiffuseEnclosure::new(surfaces, ViewFactorMatrix::infinite_parallel_plates(AREA).unwrap()).unwrap();
    let surfaces = [COIL, SINK].into_iter().enumerate().map(|(surface, cell)| RadiatingCell {
        surface, cell, material_id: MaterialId(0), curve: curve.curve(), material: cells[surface].point,
        selection: MaterialPropertySelection::SingleClaimOnly }).collect::<Vec<_>>();
    let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
        power_tolerance_w: 1e-9, radiosity_tolerance_w_m2: 1e-9 };
    let assembly = MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &contacts, radiation: Some(&radiation) };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey { seed: 19, kernel_id: 1, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        let prepared = assembly.prepare(&cx, &thermal).unwrap();
        run(&cx, &gate, &prepared, &thermal, &driver, voltage)
    })
}

#[test]
fn actual_contact_and_radiation_cooling_changes_biased_current_and_motion() {
    let run = |emissivity, resistance| with_heat(emissivity, resistance, 600.0,
        |cx, _, prepared, thermal, driver, voltage| {
            let dt = 0.02;
            let next = driver.step_candidate_with_heat_network(cx, prepared, thermal, dt,
                voltage, 0.0, controls(), &[0.0; 2]).unwrap();
            let [ta, tb] = [COIL, SINK].map(|id| next.thermal().temperature_k(id).unwrap());
            let q = dt * (AREA / (0.1 + resistance) * (ta - tb)
                + fs_conduction::STEFAN_BOLTZMANN_W_M2_K4 * AREA * (ta.powi(4) - tb.powi(4)) / (2.0 / emissivity - 1.0));
            near(CAPACITY * (tb - 300.0), q, 2e-10);
            near(CAPACITY * (ta - 300.0), next.report().electrical.coil_heat_j - q, 2e-10);
            near(next.report().deposited_joule_j, next.report().electrical.coil_heat_j, 1e-12);
            near(next.report().total_energy_residual_j, 0.0, 2e-10);
            assert_eq!(next.report().source.thermal_contacts.len(), 1);
            for surface in next.report().source.radiation.as_ref().unwrap().surfaces() {
                assert_eq!(surface.phase, next.report().source.material.cells.iter()
                    .find(|c| c.cell == surface.cell).unwrap().phase);
            }
            (ta, next.report().electrical.midpoint_current_a, next.driver().circuit().extension_m())
        });
    let weak = run(0.1, 0.2);
    let contact = run(0.1, 0.001);
    let both = run(0.9, 0.001);
    assert!(both.0 < contact.0 && contact.0 < weak.0);
    assert!(both.1 > contact.1 && contact.1 > weak.1);
    assert!((both.2 - weak.2).abs() > 1e-12);
}

#[test]
fn complete_heat_source_failure_and_retry_leave_all_original_histories_intact() {
    with_heat(0.8, 0.001, 300.001, |cx, _, prepared, thermal, driver, voltage| {
        assert!(driver.step_candidate_with_heat_network(cx, prepared, thermal, 0.02,
            voltage, 0.0, controls(), &[0.0; 2]).is_err());
        near(thermal.temperature_k(COIL).unwrap(), 300.0, 1e-12);
        near(thermal.temperature_k(SINK).unwrap(), 300.0, 1e-12);
        near(driver.time_s(), 0.0, 0.0);
        let retry = || driver.step_candidate_with_heat_network(cx, prepared, thermal, 1e-7,
            voltage, 0.0, controls(), &[0.0; 2]).unwrap();
        let a = retry(); let b = retry();
        assert_eq!(a.thermal().temperature_k(COIL).unwrap().to_bits(), b.thermal().temperature_k(COIL).unwrap().to_bits());
        assert_eq!(a.driver().current_a().unwrap().to_bits(), b.driver().current_a().unwrap().to_bits());
    });
}

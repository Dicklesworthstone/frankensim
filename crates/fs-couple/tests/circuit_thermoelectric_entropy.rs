//! Source-backed thermoelectric entropy on the original four-store circuit.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use common::{boundary, config, point};
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, ResistorCell};
use fs_couple::circuit_enthalpy::entropy::CircuitEntropyAudit;
use fs_couple::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::circuit_enthalpy::thermoelectric::{self, ThermoelectricLeg};
use fs_couple::material_circuit::MaterialCircuit;
use fs_evidence::ValidityDomain;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::thermoelectric::{ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS};
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;

const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (1, ElectricalPropertyBasis::Resistivity), (3, ElectricalPropertyBasis::Resistivity),
];
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:.16e} != {b:.16e}; tolerance {tolerance:e}");
}
fn balance(a: &CircuitEntropyAudit) {
    close(a.thermal_storage_change_j_k, a.reservoir_entropy_input_j_k
        + a.conduction_generation_j_k + a.electrical_loss_generation_j_k
        + a.unclassified_source_at_cell_j_k + a.reversible_transfer_entropy_j_k
        + a.temporal_defect_j_k + a.algebraic_residual_j_k + a.accounting_residual_j_k, 2e-14);
}
fn source(s: f64, slope: f64, resistance: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = common::material_with_conductivity(resistance, 10.0, 0.0, 1.0);
    let mut claims = base.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, s-150.0*slope), (2000.0, s+1650.0*slope)], dims: SEEBECK_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 200.0, 2000.0)
            .with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "manufactured entropy test, not measured thermoelectric data".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(200.0), upper: &point(2000.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, caloric)
}
fn sources(slope: f64) -> ([MaterialCard; 4], [ResolvedHeatCapacityCurve; 4]) {
    let [(a,ca), (b,cb), (c,cc), (d,cd)] = [source(0.0,0.0,1.0), source(0.0,0.0,1.0),
        source(0.02,slope,2.0), source(-0.01,0.0,3.0)];
    ([a,b,c,d], [ca,cb,cc,cd])
}
fn initial(curves: &[ResolvedHeatCapacityCurve; 4], temperatures: [f64; 4], baths: [Boundary; 4])
    -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for i in 0..4 {
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
            mass_kg: 1.0, specific_h_j_kg: 10.0*(temperatures[i]-200.0),
            geometry: CellGeometry { start_m: [i as f64,0.0,0.0], end_m: [i as f64+1.0,0.0,0.0],
                volume_m3: 1.0, surface_area_m2: 6.0 }, boundary: baths[i],
        }) });
    }
    for (id, [a,b]) in [[0,2], [2,1], [1,3], [3,0]].into_iter().enumerate() {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: id as u64, a: CellId(a), b: CellId(b), area_m2: 1.0, conductance_w_k: 999.0,
        }) });
    }
    let laws = (0..4).map(|i| (MaterialId(i as u64), curves[i].curve() as &dyn ThermalLaw)).collect::<Vec<_>>();
    DepositionSimulation::new(events, &laws, DepositionBudget {
        max_cells: 4, max_links: 4, max_events: 8, max_substeps: 100,
    }).unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0
}
fn cells<'a>(cards: &'a [MaterialCard; 4], curves: &'a [ResolvedHeatCapacityCurve; 4],
    context: &'a QueryPoint, key: &'a PropertyKey) -> [MaterialCell<'a>; 4] {
    [0,1,2,3].map(|i| MaterialCell {
        cell: CellId(i as u64), material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
    })
}
fn resistors() -> [CellResistor; 2] {
    [(1,2), (3,3)].map(|(branch, cell)| CellResistor {
        destination: ResistorCell { branch, cell: CellId(cell) }, selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0, area_m2: 1.0,
    })
}
fn links() -> [ConductionLink; 4] {
    [0,1,2,3].map(|link| ConductionLink { link, distance_a_m: 0.5, distance_b_m: 0.5 })
}
fn legs() -> [ThermoelectricLeg; 2] {
    [ThermoelectricLeg { source_branch: 0, resistor_branch: 1, terminals: [CellId(0),CellId(1)],
        thermal_links: [0,1], selection: MaterialPropertySelection::SingleClaimOnly },
     ThermoelectricLeg { source_branch: 2, resistor_branch: 3, terminals: [CellId(1),CellId(0)],
        thermal_links: [2,3], selection: MaterialPropertySelection::SingleClaimOnly }]
}
fn circuit(network: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>, current: Option<f64>) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph { node_count: 5, branches: vec![
        (1,2,Branch::VoltageSource { port: 0 }), (2,3,Branch::Resistor { ohms: 999.0 }),
        (3,4,Branch::VoltageSource { port: 1 }), (4,0,Branch::Resistor { ohms: 999.0 }),
        (1,0, if current.is_some() { Branch::CurrentSource { port: 2 } } else { Branch::VoltageSource { port: 2 } }),
    ], transformers: vec![] }, network.resolve(thermal, || false).unwrap().resistors,
        vec![], vec![], &[0.0,0.0,current.unwrap_or(0.0)]).unwrap()
}

#[test]
fn generator_entropy_reuses_original_history_and_does_not_count_kelvin_transfer_as_a_source() {
    let (cards, curves) = sources(0.0);
    let context = point(300.0); let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards,&curves,&context,&key); let resistors = resistors(); let links = links();
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let t0 = [400.0,300.0,350.0,350.0];
    let old = initial(&curves,t0,[boundary();4]);
    let electrical = circuit(&material,&old,None);
    let plain = thermoelectric::step_candidate(&material,&old,0.1,&electrical,&[0.0;3],&legs(),
        &BASES,config(),&[0.0;4],||false).unwrap();
    let (next,a) = thermoelectric::entropy::step_candidate(&material,&old,0.1,&electrical,&[0.0;3],&legs(),
        &BASES,config(),&[0.0;4],||false).unwrap();
    assert_eq!(next.circuit(),plain.circuit());
    assert_eq!(next.thermal().cells(),plain.thermal().cells());
    assert_eq!(next.report().coupled.circuit_trials,plain.report().coupled.circuit_trials);
    assert_eq!(next.report().coupled.source,plain.report().coupled.source);
    assert!(!a.has_unclassified_sources);
    assert_eq!(a.unclassified_source_at_cell_j_k,0.0);
    close(a.reversible_transfer_entropy_j_k,0.0,2e-17);
    let t = [0,1,2,3].map(|i|next.thermal().temperature_k(CellId(i)).unwrap());
    let storage: f64 = (0..4).map(|i|10.0*(t[i]/t0[i]).ln()).sum();
    close(a.thermal_storage_change_j_k,storage,1e-13);
    let current = next.report().coupled.source.legs[0].transfer.current_a;
    close(a.electrical_loss_generation_j_k,0.1*current*current*(2.0/t[2]+3.0/t[3]),1e-13);
    assert!(a.conduction_generation_j_k>0.0 && a.electrical_loss_generation_j_k>0.0);
    close(next.thermal().energy_j(),old.energy_j(),2e-6);
    balance(&a);
}

#[test]
fn variable_seebeck_bulk_entropy_closes_for_both_current_directions() {
    let (cards,curves)=sources(2e-4);
    let context=point(300.0); let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,&curves,&context,&key); let resistors=resistors(); let links=links();
    let material=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let old=initial(&curves,[400.0,300.0,350.0,350.0],[boundary();4]);
    for current in [-0.5,0.5] {
        let electrical=circuit(&material,&old,Some(current));
        let (next,a)=thermoelectric::entropy::step_candidate(&material,&old,0.1,&electrical,
            &[0.0,0.0,current],&legs(),&BASES,config(),&[0.0;4],||false).unwrap();
        let leg=&next.report().coupled.source.legs[0];
        let [ta,tb]=leg.material.temperature_k(); let tc=leg.material.bulk_temperature_k();
        let [qa,qb]=leg.transfer.terminal_heat_w; let qc=leg.transfer.bulk_heat_w;
        assert!(qc*current>0.0);
        assert!((qa/ta+qb/tb).abs()>1e-3, "omitting bulk entropy would fail this balance");
        close(qa/ta+qc/tc+qb/tb,0.0,2e-17);
        assert!(!a.has_unclassified_sources);
        assert_eq!(a.unclassified_source_at_cell_j_k,0.0);
        close(a.reversible_transfer_entropy_j_k,0.0,2e-17);
        assert!(a.electrical_loss_generation_j_k>0.0);
        balance(&a);
    }
}

#[test]
fn thermoelectric_transfer_does_not_hide_authored_sources_or_real_bath_entropy() {
    let (cards,curves)=sources(2e-4);
    let context=point(300.0); let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,&curves,&context,&key); let resistors=resistors(); let links=links();
    let material=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let baths=[Boundary {ambient_k:450.0,htc_w_m2_k:0.2,source_w:1e-10,..boundary()},
        Boundary {ambient_k:280.0,htc_w_m2_k:0.1,..boundary()},boundary(),boundary()];
    let old=initial(&curves,[400.0,300.0,350.0,350.0],baths);
    let electrical=circuit(&material,&old,Some(0.5));
    let (next,a)=thermoelectric::entropy::step_candidate(&material,&old,0.1,&electrical,
        &[0.0,0.0,0.5],&legs(),&BASES,config(),&[0.0;4],||false).unwrap();
    let ta=next.thermal().temperature_k(CellId(0)).unwrap();
    let tb=next.thermal().temperature_k(CellId(1)).unwrap();
    assert!(a.has_unclassified_sources);
    assert!(a.unclassified_source_at_cell_j_k>0.0);
    close(a.unclassified_source_at_cell_j_k,1e-11/ta,1e-28);
    // Each junction has two original 1 m2 links: exposed surface is 4 m2.
    close(a.reservoir_entropy_input_j_k,0.1*(0.8*(450.0-ta)/450.0+0.4*(280.0-tb)/280.0),1e-14);
    close(a.reversible_transfer_entropy_j_k,0.0,2e-17);
    balance(&a);
}

#[test]
fn requested_entropy_preserves_last_accepted_state_on_late_cancellation_and_budget_failure() {
    let (cards,curves)=sources(2e-4);
    let context=point(300.0); let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,&curves,&context,&key); let resistors=resistors(); let links=links();
    let material=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let old=initial(&curves,[400.0,300.0,350.0,350.0],[boundary();4]);
    let electrical=circuit(&material,&old,None);
    let saved=old.cells().clone(); let saved_circuit=electrical.clone();
    let mut polls=0;
    let (success,_)=thermoelectric::entropy::step_candidate(&material,&old,0.1,&electrical,&[0.0;3],
        &legs(),&BASES,config(),&[0.0;4],||{polls+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(thermoelectric::entropy::step_candidate(&material,&old,0.1,&electrical,&[0.0;3],
        &legs(),&BASES,config(),&[0.0;4],||{actual+=1;actual==polls}),Err(CircuitEnthalpyError::Cancelled)));
    let mut limited=config(); limited.maximum_circuit_trials=success.report().coupled.circuit_trials-1;
    assert!(matches!(thermoelectric::entropy::step_candidate(&material,&old,0.1,&electrical,&[0.0;3],
        &legs(),&BASES,limited,&[0.0;4],||false),Err(CircuitEnthalpyError::CircuitBudget)));
    assert_eq!(old.cells(),&saved); assert_eq!(old.time_s(),0.0); assert_eq!(electrical,saved_circuit);
}

#[test]
fn thermoelectric_numerical_entropy_converges_under_clock_refinement() {
    let (cards,curves)=sources(2e-4);
    let context=point(300.0); let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,&curves,&context,&key); let resistors=resistors(); let links=links();
    let material=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let t0=[400.0,300.0,350.0,350.0];
    let mut defects=Vec::new();
    for count in [16,32,64] {
        let mut thermal=initial(&curves,t0,[boundary();4]);
        let mut electrical=circuit(&material,&thermal,None);
        let energy0=thermal.energy_j(); let (mut storage,mut defect)=(0.0,0.0);
        for i in 1..=count {
            let (next,a)=thermoelectric::entropy::step_candidate(&material,&thermal,i as f64/count as f64,
                &electrical,&[0.0;3],&legs(),&BASES,config(),&[0.0;4],||false).unwrap();
            balance(&a); assert!(!a.has_unclassified_sources);
            storage+=a.thermal_storage_change_j_k; defect+=a.temporal_defect_j_k;
            close(next.thermal().energy_j(),energy0,2e-6);
            let (t,c,_)=next.into_parts(); thermal=t; electrical=c;
        }
        let exact:f64=(0..4).map(|i|10.0*(thermal.temperature_k(CellId(i as u64)).unwrap()/t0[i]).ln()).sum();
        close(storage,exact,2e-12); assert!(defect>0.0); defects.push(defect);
    }
    for pair in defects.windows(2) { assert!((0.45..0.55).contains(&(pair[1]/pair[0])),"{defects:?}"); }
}

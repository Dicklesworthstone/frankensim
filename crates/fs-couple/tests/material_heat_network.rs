//! Same material assembly through real electrical, contact and radiation owners.
#[path = "material_heat_network/layers.rs"]
mod layers;
#[path = "../../fs-conduction/tests/radiation_material/common.rs"]
mod common;
use common::*;
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_conduction::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_couple::circuit_enthalpy::{self, CircuitEnthalpyInput, CircuitEnthalpySolve,
    ResistorCell};
use fs_couple::circuit_enthalpy::contacts;
use fs_couple::circuit_enthalpy::contacts::material::*;
use fs_couple::circuit_enthalpy::contacts::heat::*;
use fs_couple::circuit_enthalpy::material::*;
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_exec::{CancelGate, Cx};
use fs_material::conductor::{ElectricalPropertyBasis, ELECTRICAL_RESISTIVITY_DIMS};
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_DIMS as ED,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as EP};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_matdb::*;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::Temperature;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

fn source(e: f64, r: f64) -> MaterialCard {
    let card = card(e, e);
    let mut claims = card.claims().clone();
    for (name, dims, a, b) in [("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, r, 2.0*r),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0, 2.0)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0,a),(500.0,b)], dims },
            validity: fs_evidence::ValidityDomain::unconstrained().with("T",300.0,500.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic combined transport".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap()
}
fn interface(cards: &[MaterialCard;2], rt: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [(RP,RD,rt),(EP,ED,0.05)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name,dims),
            value: PropertyValue::Curve { abscissa:"aT".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(300.0,value),(500.0,2.0*value)],dims },
            validity:fs_evidence::ValidityDomain::unconstrained().with("aT",300.0,500.0)
                .with("bT",300.0,500.0).with("pressure",100000.0,100000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic electrical/thermal interface".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {material:cards[0].id().clone(),texture_frame:"a".into()},
        SurfaceSpec {material:cards[1].id().clone(),texture_frame:"b".into()},
        SystemContext {medium:"fixed contact".into(),third_body:None,
            environment:"authored".into(),history:"fixed".into()},claims,vec![]).unwrap()
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {maximum_cells:2,maximum_circuit_trials:300,maximum_newton_iterations:32,
        difference_step_j:1e-4,cell_energy_tolerance_j:1e-7,thermal_energy_tolerance_j:1e-7,total_energy_tolerance_j:1e-7,
        circuit:MaterialCircuitSolve {energy_tolerance_j:1e-7,power_tolerance_w:1e-7},
        newton:NewtonKrylovConfig {absolute_tolerance:1e-9,relative_tolerance:1e-14,
            forcing_minimum:1e-4,linear_restart:4,..NewtonKrylovConfig::default()}}
}

// Borrow the SAME laws/cards into each existing owner. The helper returns no
// borrowed state, and cannot advance or clone an independent physical timeline.
fn with_case<R>(e: f64, rt: f64, run: impl FnOnce(&Cx<'_>, &CancelGate,
    &MaterialHeatNetwork<'_>, DepositionSimulation<'_>, MaterialCircuit) -> R) -> R {
    let cards = [source(e,1.0),source(0.6,2.0)];
    let curves = [curve(&cards[0],100.0),curve(&cards[1],100.0)];
    let context = context();
    let ids = [CellId(30),CellId(10)];
    let key = PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0,1].map(|i|MaterialCell {cell:ids[i],material_id:MaterialId(i as u64),curve:&curves[i],
        point:binding(&cards[i],&context),conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let resistors = [0,1].map(|i|CellResistor {destination:ResistorCell {branch:i+1,cell:ids[i]},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:1.0});
    let links = [ConductionLink {link:7,distance_a_m:0.5,distance_b_m:0.5}];
    let bulk = MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let card = interface(&cards,rt);
    let p = QueryPoint::new().with("aT",400.0).unwrap().with("bT",400.0).unwrap()
        .with("aM",0.0).unwrap().with("bM",0.0).unwrap().with("pressure",100000.0).unwrap();
    let interface_point = ThermalInterfacePoint {interface:&card,materials:[&cards[0],&cards[1]],context:&p,
        temperature_axes:["aT","bT"],liquid_fraction_axes:[Some("aM"),Some("bM")],mean_temperature_axis:None};
    let electrical_contacts = [MaterialElectricalContact {link:7,branch:1,point:interface_point,
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:0.3}];
    let bases = [(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let electrical = ElectricalContactNetwork {bulk:&bulk,contacts:&electrical_contacts,electrical_bases:&bases};
    let thermal_contacts = [MaterialContactLink {link:7,point:interface_point,
        sides:[0,1].map(|i|ContactCell {cell:ids[i],material_id:MaterialId(i as u64),curve:&curves[i]}),
        selection:ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)}];
    let surfaces = [0,1].map(|i|RadiatingCell {surface:i,cell:ids[i],material_id:MaterialId(i as u64),curve:&curves[i],
        material:binding(&cards[i],&context),selection:MaterialPropertySelection::SingleClaimOnly});
    let enclosure = enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
    let radiation = MaterialRadiation {enclosure:&enclosure,cells:&surfaces,
        power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
    let assembly = MaterialHeatNetwork {electrical:&electrical,thermal_contacts:&thermal_contacts,radiation:Some(&radiation)};
    let mut events = [0,1].map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:ids[i],segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:1.0,
        specific_h_j_kg:if i==0 {150.0} else {-50.0},geometry:CellGeometry {start_m:[i as f64,0.0,0.0],
            end_m:[i as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:2.5},
        boundary:Boundary {ambient_k:400.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})}).to_vec();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:7,a:ids[0],b:ids[1],area_m2:0.1,conductance_w_k:999.0})});
    let thermal = DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),
        (MaterialId(1),&curves[1] as &dyn ThermalLaw)],DepositionBudget {
            max_cells:2,max_links:1,max_events:3,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-7},||Ok(())).unwrap().0;
    let sample = electrical.resolve(&thermal,||false).unwrap();
    let circuit = MaterialCircuit::new_with_contacts(CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:0.01}),(1,2,Branch::Resistor {ohms:0.0}),
        (2,0,Branch::Resistor {ohms:0.0})],transformers:vec![]},sample.bulk.resistors,sample.contacts,
        vec![],vec![1.0],&[]).unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate,arena,fs_exec::StreamKey {seed:73,kernel_id:51,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        run(&cx,&gate,&assembly,thermal,circuit)
    })
}

#[test]
fn composed_sources_match_the_existing_manual_physics_and_do_not_accumulate() {
    with_case(0.2,0.2,|cx,_,assembly,thermal,circuit| {
        let selected = assembly.prepare(cx,&thermal).unwrap();
        let first = selected.resolve(cx,&thermal).unwrap();
        let second = selected.resolve(cx,&thermal).unwrap();
        assert_eq!(first.bulk.conductances_w_k,second.bulk.conductances_w_k);
        assert_eq!(first.bulk.sources_w,second.bulk.sources_w);
        let tc = &first.bulk.source.thermal_contacts[0];
        near(tc.conductance_w_k,1.0/(1.0/tc.bulk_conductance_w_k+tc.interface.value_m2_k_per_w()/0.1));
        let radiation = first.bulk.source.radiation.as_ref().unwrap();
        assert!(radiation.sources_w()[&CellId(30)]<0.0 && radiation.sources_w()[&CellId(10)]>0.0);
        let expected = assembly.electrical.prepare(&thermal,||false).unwrap();
        let destinations = assembly.electrical.bulk.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
        let routes = assembly.electrical.heat_destinations(&thermal).unwrap();
        let original = contacts::step_candidate(&thermal,0.0005,
            CircuitEnthalpyInput {circuit:&circuit,ports:&[],destinations:&destinations},&routes,controls(),&[0.0;2],|state| {
                let mut source = expected.resolve(state,||false).map_err(|e|e.to_string())?;
                apply_material_contacts(state,assembly.thermal_contacts,&mut source.bulk.conductances_w_k,||false).map_err(|e|e.to_string())?;
                let radiation = assembly.radiation.unwrap().resolve(cx,state).map_err(|e|e.to_string())?;
                radiation.add_to(&mut source.bulk.sources_w).map_err(|e|e.to_string())?;
                Ok::<_,String>(source)
            },||false).unwrap();
        let actual = selected.step_candidate(cx,&thermal,0.0005,&circuit,&[],controls(),&[0.0;2]).unwrap();
        assert_eq!(original.thermal().cells(),actual.thermal().cells());
        assert_eq!(original.circuit(),actual.circuit());
        near(actual.thermal().energy_j()+actual.report().electrical.stored_after_j,150.0);
        assert_eq!(thermal.time_s(),0.0);
        assert_eq!(circuit.charges_c(),&[1.0]);
    });
}

#[test]
fn data_only_emissivity_and_contact_changes_reach_both_bodies_and_electrical_history() {
    let run = |e,rt| with_case(e,rt,|cx,_,assembly,mut thermal,mut circuit| {
        let selected = assembly.prepare(cx,&thermal).unwrap();
        for n in 1..=40 {
            let proposal = selected.step_candidate(cx,&thermal,f64::from(n)*0.0005,&circuit,&[],controls(),&[0.0;2]).unwrap();
            let report = proposal.report();
            assert_eq!(report.electrical.contacts.len(),1);
            assert_eq!(report.source.material.thermal_contacts.len(),1);
            for surface in report.source.material.radiation.as_ref().unwrap().surfaces() {
                let cell = report.source.material.material.cells.iter().find(|c|c.cell==surface.cell).unwrap();
                assert_eq!(cell.phase,surface.phase);
            }
            near(proposal.thermal().energy_j()+report.electrical.stored_after_j,150.0);
            let (t,c,_) = proposal.into_parts(); thermal=t; circuit=c;
        }
        [thermal.temperature_k(CellId(30)).unwrap(),thermal.temperature_k(CellId(10)).unwrap(),circuit.charges_c()[0]]
    });
    let original = run(0.2,0.2);
    for changed in [run(0.8,0.2),run(0.2,2.0)] {
        for i in 0..3 { assert!((original[i]-changed[i]).abs()>1e-8); }
    }
}

#[test]
fn absent_heat_paths_reduce_to_contact_only_and_cancellation_keeps_history() {
    with_case(0.2,0.2,|cx,gate,assembly,thermal,circuit| {
        let empty = MaterialHeatNetwork {electrical:assembly.electrical,thermal_contacts:&[],radiation:None};
        let selected = empty.prepare(cx,&thermal).unwrap();
        let actual = selected.step_candidate(cx,&thermal,0.0005,&circuit,&[],controls(),&[0.0;2]).unwrap();
        let original = assembly.electrical.step_candidate(&thermal,0.0005,&circuit,&[],controls(),&[0.0;2],||false).unwrap();
        assert_eq!(original.thermal().cells(),actual.thermal().cells());
        assert_eq!(original.circuit(),actual.circuit());
        gate.request();
        assert!(matches!(assembly.prepare(cx,&thermal),Err(HeatNetworkError::Cancelled)));
        assert!(matches!(selected.resolve(cx,&thermal),Err(HeatNetworkError::Cancelled)));
        assert!(matches!(selected.step_candidate(cx,&thermal,0.0005,&circuit,&[],controls(),&[0.0;2]),
            Err(circuit_enthalpy::CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.time_s(),0.0); assert_eq!(circuit.charges_c(),&[1.0]);
    });
}

#[test]
fn separately_authored_radiation_cannot_change_the_same_bodys_pressure() {
    with_case(0.2,0.2,|cx,_,assembly,thermal,_| {
        let original = assembly.radiation.unwrap();
        let other = context().with("pressure",105000.0).unwrap();
        let cells = original.cells.iter().map(|c|RadiatingCell {surface:c.surface,cell:c.cell,
            material_id:c.material_id,curve:c.curve,selection:c.selection.clone(),
            material:fs_material::thermal::source::ThermalMaterialPoint {context:&other,..c.material}}).collect::<Vec<_>>();
        let radiation = MaterialRadiation {enclosure:original.enclosure,cells:&cells,
            power_tolerance_w:original.power_tolerance_w,radiosity_tolerance_w_m2:original.radiosity_tolerance_w_m2};
        // Each source is valid alone. They disagree on the state of the SAME body.
        assert!(radiation.resolve(cx,&thermal).is_ok());
        let wrong = MaterialHeatNetwork {electrical:assembly.electrical,thermal_contacts:assembly.thermal_contacts,radiation:Some(&radiation)};
        assert!(matches!(wrong.prepare(cx,&thermal),Err(HeatNetworkError::Invalid(_))));
        assert!(assembly.prepare(cx,&thermal).is_ok());
    });
}

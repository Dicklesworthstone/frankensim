#![allow(dead_code)]
//! Synthetic data; the original material, enclosure, beam and circuit owners.
#[path = "../strip_enthalpy/common.rs"]
mod strip;
pub use strip::*;
use std::collections::BTreeMap;
use fs_conduction::{GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity,
    ViewFactorMatrix, EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_couple::circuit_enthalpy::ResistorCell;
use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
use fs_couple::circuit_enthalpy::contacts::material::{ElectricalContactNetwork, MaterialElectricalContact};
use fs_couple::circuit_enthalpy::material::{MaterialNetwork, MaterialCell, CellResistor, ConductionLink};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_couple::strip_enthalpy::bridge::{BridgeEnthalpyConfig, BridgeThermalBinding};
use fs_couple::strip_enthalpy::bridge::electrical::ElectricalConfig;
use fs_couple::strip_enthalpy::LayerThermalBinding;
use fs_exec::{Cx, CancelGate};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as EP,
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as ED};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS, MaterialPropertySelection};
use fs_material::thermal::{interface::ThermalInterfacePoint, source::ThermalMaterialPoint};
use fs_material::thermovisco::{Limits, source::{ResolvedThermoMaxwell, ThermoMaxwellSourceQuery}};
use fs_matdb::*;
use fs_qty::{Dims, Density, Pressure, Time};
use fs_phs::circuit::{CircuitGraph, Branch};
use fs_solid::bridge::{self as solid, Bridge, Support};
use fs_thermal::enthalpy::*;

pub const CAP: f64=0.01;
pub const CHARGE: f64=0.1;
pub const DT: f64=0.02;
pub fn controls()->ElectricalConfig {
    ElectricalConfig {bridge:BridgeEnthalpyConfig {feedback:config(),maximum_deflection_m:0.009},
        circuit:MaterialCircuitSolve {energy_tolerance_j:1e-10,power_tolerance_w:1e-9},
        maximum_circuit_trials:1000}
}
pub fn mechanics()->solid::Config {
    solid::Config {section:config().mechanics,maximum_deflection_m:0.009}
}
pub fn map(m:&[Material;2],single:bool)->Vec<BridgeThermalBinding<'_>> {
    (0..if single {1}else{2}).map(|i|BridgeThermalBinding {axial_cell:20,
        section:LayerThermalBinding {layer:LAYERS[i],cell:CELLS[i],
            material:&m[i].mechanical,caloric:&m[i].caloric}}).collect()
}
fn claim(name:&str,dims:Dims,value:f64)->PropertyClaim {
    PropertyClaim {key:PropertyKey::new(name,dims),value:PropertyValue::Scalar {value,dims},
        validity:fs_evidence::ValidityDomain::unconstrained(),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance {source:"synthetic heat/beam coupling".into(),
            license:"CC0-1.0".into(),artifact:None}}
}
fn heat_material(cp:f64,e:f64)->Material {
    let old=material_with_transport(cp,3e8,0.05,0.0,vec![(280.0,0.002),(320.0,0.002)],
        None,vec![(280.0,0.08),(320.0,0.12)]);
    let mut claims=old.card.claims().clone();
    claims.insert_claim(claim(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS,e)).unwrap();
    let card=MaterialCard::assemble(old.card.id().clone(),claims,vec![]).unwrap();
    // Re-resolve BOTH constitutive and caloric data after adding the optical
    // property: no retained curve may cite the previous card's content identity.
    let mechanical=ResolvedThermoMaxwell::try_new(&card,ThermoMaxwellSourceQuery {
        equilibrium_modulus:&PropertyKey::new("equilibrium",Pressure::DIMS),
        branches:&[(PropertyKey::new("branch",Pressure::DIMS),PropertyKey::new("time",Time::DIMS))],
        modulus_scale:&PropertyKey::new("modulus_scale",Dims([0;6])),
        relaxation_scale:&PropertyKey::new("shift",Dims([0;6])),
        free_strain:&PropertyKey::new("free",Dims([0;6])),temperature_axis:"T",
        lower:&point(280.0),upper:&point(320.0),selection:SelectionPolicy::SingleClaimOnly,
    },Limits {max_strain:0.04,max_strain_rate_s:10.0},||false).unwrap();
    let caloric=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("cp",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",
        lower:&point(280.0),upper:&point(320.0),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,-20.0*cp,1e-8).unwrap();
    Material {card,mechanical,caloric}
}

pub fn arrival(id:u64,axial:u64,start:f64,length:f64,material:usize)->TimedEvent {
    let mut c=cell(CellId(id),material);
    c.segment_id=axial; c.body_id=1;
    c.mass_kg=DENSITY*AREA*length;
    c.geometry.start_m=[start,0.0,geometry(0).centroid_z_m];
    c.geometry.end_m=[start+length,0.0,geometry(0).centroid_z_m];
    c.geometry.volume_m3=AREA*length;
    TimedEvent {time_s:DT,event:Event::Birth(c)}
}

// `single` makes the second thermal/radiative receiver passive mechanically.
// Every returned borrow stays inside the callback's original card/law lifetime.
pub fn with_case<R>(e:f64,rt:f64,extra_heat:bool,single:bool,later:Vec<TimedEvent>,
    run:impl FnOnce(&[Material;2],&Cx<'_>,&CancelGate,&MaterialHeatNetwork<'_>,
        &DepositionSimulation<'_>,&Bridge<'_>,&MaterialCircuit)->R)->R {
    let m=[heat_material(5.0,e),heat_material(10.0,0.6)];
    let mut events=Vec::new();
    for i in 0..2 {
        let mut c=cell(CELLS[i],i); c.segment_id=20;
        c.body_id=if single && i==1 {2}else{1};
        c.specific_h_j_kg=if i==0 {50.0}else{-100.0};
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(c)});
    }
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:5,a:CELLS[0],b:CELLS[1],area_m2:0.001,conductance_w_k:99.0})});
    events.extend(later);
    let thermal=DepositionSimulation::new(events,&[
        (MaterialId(0),m[0].caloric.curve() as &dyn ThermalLaw),
        (MaterialId(1),m[1].caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:8,max_links:8,max_events:16,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.01,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    let mut beam=Bridge::new(&[(20,LENGTH)],Support::CantileverLeft,
        solid::Budget {max_cells:4,max_layers_per_cell:2,max_branch_states:16}).unwrap();
    for i in 0..if single {1}else{2} {
        beam=beam.deposit_bonded_material_layer(LAYERS[i],geometry(i),&m[i].mechanical,&m[i].caloric,
            &BTreeMap::from([(20,thermal.temperature_k(CELLS[i]).unwrap())]),&mut ||false).unwrap();
    }
    let context=point(300.0);let key=PropertyKey::new("conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=[0,1].map(|i|MaterialCell {cell:CELLS[i],material_id:MaterialId(i as u64),curve:m[i].caloric.curve(),
        point:ThermalMaterialPoint {material:&m[i].card,context:&context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let resistors=[CellResistor {destination:ResistorCell {branch:1,cell:CELLS[0]},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:LENGTH,area_m2:AREA}];
    let links=[ConductionLink {link:5,distance_a_m:0.0005,distance_b_m:0.0005}];
    let bulk=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let mut claims=ClaimSet::new();
    claims.insert_claim(claim(RP,RD,rt)).unwrap();
    claims.insert_claim(claim(EP,ED,1.0)).unwrap();
    let card=InterfaceSystemCard::assemble(
        SurfaceSpec {material:m[0].card.id().clone(),texture_frame:"a".into()},
        SurfaceSpec {material:m[1].card.id().clone(),texture_frame:"b".into()},
        SystemContext {medium:"fixed synthetic interface".into(),third_body:None,
            environment:"test".into(),history:"fixed".into()},claims,vec![]).unwrap();
    let ip=QueryPoint::new().with("Ta",310.0).unwrap().with("Tb",290.0).unwrap();
    let interface=ThermalInterfacePoint {interface:&card,materials:[&m[0].card,&m[1].card],context:&ip,
        temperature_axes:["Ta","Tb"],liquid_fraction_axes:[None,None],mean_temperature_axis:None};
    let contacts=[MaterialElectricalContact {link:5,branch:1,point:interface,
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:0.25}];
    let bases=[(1,ElectricalPropertyBasis::Resistivity)];
    let electrical=ElectricalContactNetwork {bulk:&bulk,contacts:&contacts,electrical_bases:&bases};
    let thermal_contacts=[MaterialContactLink {link:5,point:interface,
        sides:[0,1].map(|i|ContactCell {cell:CELLS[i],material_id:MaterialId(i as u64),curve:m[i].caloric.curve()}),
        selection:ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)}];
    let surfaces=[0,1].map(|i|RadiatingCell {surface:i,cell:CELLS[i],material_id:MaterialId(i as u64),
        curve:m[i].caloric.curve(),material:cells[i].point,selection:MaterialPropertySelection::SingleClaimOnly});
    let (complex,positions)=box_grid([1,1,1],[0.1,0.01,0.001]);
    let mesh=ConductionMesh::new(complex,positions).unwrap();
    let radiation_surfaces=[0,1].map(|i| {
        let phase=m[i].caloric.curve().state_at_specific_enthalpy(thermal.cells()[&CELLS[i]].specific_h_j_kg).unwrap();
        let e=SurfaceEmissivity::from_thermal_state("synthetic",cells[i].point,phase,
            MaterialPropertySelection::SingleClaimOnly).unwrap().0;
        RadiationSurface::new(&mesh,format!("surface-{i}"),|f|on_box_face(f.centroid[2],i as f64*0.001),e).unwrap()
    });
    // Infinite-plate view factors are an algebraic limiting fixture, not ray geometry.
    let enclosure=GrayDiffuseEnclosure::new(radiation_surfaces.into_iter().collect(),
        ViewFactorMatrix::infinite_parallel_plates(0.001).unwrap()).unwrap();
    let radiation=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,
        power_tolerance_w:1e-10,radiosity_tolerance_w_m2:1e-8};
    let network=MaterialHeatNetwork {electrical:&electrical,
        thermal_contacts:if extra_heat {&thermal_contacts}else{&[]},
        radiation:if extra_heat {Some(&radiation)}else{None}};
    let source=electrical.resolve(&thermal,||false).unwrap();
    let circuit=MaterialCircuit::new_with_contacts(CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::Capacitor {farads:CAP}),(1,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},
        source.bulk.resistors,source.contacts,vec![],vec![CHARGE],&[]).unwrap();
    let gate=CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:73,kernel_id:51,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        run(&m,&cx,&gate,&network,&thermal,&beam,&circuit)
    })
}

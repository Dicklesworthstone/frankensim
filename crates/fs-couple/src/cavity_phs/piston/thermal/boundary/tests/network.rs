//! Synthetic cards through the actual circuit, thermal network and gas owners.
mod warping;
use super::*;
use crate::cavity_phs::piston::thermal::boundary::network::*;
use crate::cavity_phs::piston::thermal::boundary::network::material::*;
use crate::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use crate::circuit_enthalpy::{CircuitEnthalpyError,CircuitEnthalpySolve,ResistorCell};
use crate::circuit_enthalpy::material::*;
use crate::material_circuit::{MaterialCircuit,MaterialCircuitSolve};
use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
use fs_material::phase::{ResolvedHeatCapacityCurve,HeatCapacitySourceQuery,SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection,THERMAL_CONDUCTIVITY_DIMS,SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_qty::{Temperature,Density};
use fs_phs::circuit::{CircuitGraph,Branch};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

fn point(t:f64)->QueryPoint {QueryPoint::new().with("T",t).unwrap().with("lot",7.0).unwrap()}
fn material(k:f64,r:f64,slope:f64)->(MaterialCard,ResolvedHeatCapacityCurve) {
    let mut claims=ClaimSet::new();
    for (name,dims,a,b) in [
        ("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS,r*0.0025*(1.0-200.0*slope),r*0.0025*(1.0+400.0*slope)),
        ("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS,k,2.0*k),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,500.0,500.0),
        ("density",Density::DIMS,4000.0,4000.0)] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(200.0,a),(800.0,b)],dims},
            validity:ValidityDomain::unconstrained().with("T",200.0,800.0).with("lot",7.0,7.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"synthetic network/gas material".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"authored test".into(),revision:0},claims,vec![]).unwrap();
    let curve=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",lower:&point(200.0),upper:&point(800.0),
        selection:SelectionPolicy::SingleClaimOnly},SolidLiquidPhase::Solid,0.0,1e-5).unwrap();
    (card,curve)
}
fn cell(i:usize)->CellId {CellId(if i==0 {30}else{10})}
fn bindings<'a>(cards:&'a [MaterialCard;2],curves:&'a [ResolvedHeatCapacityCurve;2],context:&'a QueryPoint,key:&'a PropertyKey)
    ->[MaterialCell<'a>;2] {
    [0,1].map(|i|MaterialCell {cell:cell(i),material_id:MaterialId(i as u64),curve:curves[i].curve(),
        point:ThermalMaterialPoint {material:&cards[i],context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((key,MaterialPropertySelection::SingleClaimOnly))})
}
fn resistors()->[CellResistor;2] {[0,1].map(|i|CellResistor {
    destination:ResistorCell {branch:i+1,cell:cell(i)},selection:MaterialPropertySelection::SingleClaimOnly,
    length_m:0.01,area_m2:0.000025})}
fn initial<'a>(curves:&'a [ResolvedHeatCapacityCurve;2],ambient:f64)->DepositionSimulation<'a> {
    let mut events:Vec<_>=[0,1].map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:cell(i),segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:0.001,
        specific_h_j_kg:100000.0,geometry:CellGeometry {start_m:[i as f64,0.0,0.0],end_m:[i as f64+0.01,0.0,0.0],
            volume_m3:0.001/4000.0,surface_area_m2:0.04},boundary:Boundary {ambient_k:400.0,
                htc_w_m2_k:ambient,reservoir:None,source_w:0.0}})}).to_vec();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {id:7,a:cell(0),b:cell(1),
        area_m2:0.002,conductance_w_k:999.0})});
    DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],DepositionBudget {
        max_cells:2,max_links:1,max_events:3,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-7},||Ok(())).unwrap().0
}
fn circuit(binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'_>)->MaterialCircuit {
    MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Capacitor {farads:0.01}),(1,2,Branch::Resistor {ohms:0.0}),
        (2,0,Branch::Resistor {ohms:0.0})],transformers:vec![]},
        binding.resolve(thermal,||false).unwrap().resistors,vec![],vec![1.0],&[]).unwrap()
}
fn solve()->NetworkGasSolve {NetworkGasSolve {network:CircuitEnthalpySolve {
    maximum_cells:2,maximum_circuit_trials:300,maximum_newton_iterations:32,difference_step_j:1e-4,
    cell_energy_tolerance_j:1e-7,thermal_energy_tolerance_j:1e-7,total_energy_tolerance_j:1e-7,
    circuit:MaterialCircuitSolve {energy_tolerance_j:1e-7,power_tolerance_w:1e-7},
    newton:NewtonKrylovConfig {absolute_tolerance:1e-8,relative_tolerance:0.0,
        forcing_minimum:1e-4,linear_restart:4,..NewtonKrylovConfig::default()}},
    maximum_boundary_solves:300,total_energy_tolerance_j:1e-6}}
fn port(dt:f64)->NetworkGasDrive {NetworkGasDrive {cell:cell(0),input:PistonInput {duration_s:dt,..drive(0.0)},
    piston:controls(),exchange:exchange()}}
fn geometry()->GasConductionGeometry {GasConductionGeometry {area_m2:0.01,gas_distance_m:0.001,cell_distance_m:0.02}}
fn near(a:f64,b:f64) {assert!((a-b).abs()<2e-6,"{a} != {b}");}

#[test]
fn conductivity_data_changes_finite_supply_two_bodies_and_actual_gas_motion() {
    let run=|k| {
        let (a,ca)=material(k,2.0,0.002);let (b,cb)=material(0.4,4.0,0.001);
        let cards=[a,b];let curves=[ca,cb];let context=point(400.0);
        let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells=bindings(&cards,&curves,&context,&key);let resistors=resistors();
        let links=[ConductionLink {link:7,distance_a_m:0.001,distance_b_m:0.001}];
        let binding=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
        let mut thermal=initial(&curves,0.0);let mut circuit=circuit(&binding,&thermal);let mut gas=piston();
        let old=gas.clone();let h0=thermal.energy_j();let mut external=0.0;
        for i in 1..=40 {
            let target=f64::from(i)*0.0005;let drive=port(target-thermal.time_s());
            let candidate=gas.material_network_candidate(&binding,&thermal,target,&circuit,&[],drive,solve(),geometry(),&[0.0;2],||false).unwrap();
            assert!(candidate.step().total_energy_residual_j.abs()<=solve().total_energy_tolerance_j);
            let (t,c,g,r)=candidate.into_parts();
            near(r.coupled.external_network_heat_j,0.0);
            near(r.network.non_electrical_heat_j,-r.boundary.exchanged_heat_j);
            let body=r.network.source.cells.iter().find(|s|s.cell==cell(0)).unwrap();
            assert_eq!(body.conductivity.as_ref(),Some(&r.interface));
            near(body.phase.temperature_k(),r.boundary.boundary_temperature_k);
            external+=r.boundary.piston.exterior_pressure_work_j;
            near(t.energy_j()-h0+r.network.electrical.stored_after_j-50.0
                +g.gas().mass_kg()*(g.gas().state().specific_internal_energy_j_kg()-old.gas().state().specific_internal_energy_j_kg())
                +0.5*g.mass_kg()*g.velocity_m_s().powi(2)+external,0.0);
            thermal=t;circuit=c;gas=g;
        }
        assert!(gas.position_m()>old.position_m());assert!(circuit.charges_c()[0]<1.0);
        (thermal.temperature_k(cell(0)).unwrap(),thermal.temperature_k(cell(1)).unwrap(),circuit.charges_c()[0],gas)
    };
    let a=run(0.2);let b=run(2.0);
    assert!((a.0-b.0).abs()>1e-4 && (a.1-b.1).abs()>1e-4);
    assert!((a.2-b.2).abs()>1e-9);
    assert!((a.3.position_m()-b.3.position_m()).abs()>1e-9);
    assert!((a.3.gas().state().state().sound_speed-b.3.gas().state().state().sound_speed).abs()>1e-6);
}

#[test]
fn isolated_gas_patch_retains_exact_original_circuit_and_motion_owners() {
    let (a,ca)=material(0.2,2.0,0.002);let (b,cb)=material(0.4,4.0,0.001);
    let cards=[a,b];let curves=[ca,cb];let context=point(400.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=bindings(&cards,&curves,&context,&key);let resistors=resistors();
    let links=[ConductionLink {link:7,distance_a_m:0.001,distance_b_m:0.001}];
    let binding=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let thermal=initial(&curves,0.0);let electrical=circuit(&binding,&thermal);let gas=piston();
    let independent=binding.step_candidate(&thermal,0.0005,&electrical,&[],solve().network,&[0.0;2],||false).unwrap();
    let mut drive=port(0.0005);drive.input.heat_j=0.01;
    let candidate=gas.material_network_candidate(&binding,&thermal,0.0005,&electrical,&[],drive,solve(),
        GasConductionGeometry {area_m2:0.0,..geometry()},&[0.0;2],||false).unwrap();
    assert_eq!(candidate.network().thermal().cells(),independent.thermal().cells());
    assert_eq!(candidate.network().circuit(),independent.circuit());
    let mut expected=gas.clone();expected.step(drive.input,drive.piston,||false).unwrap();
    assert_eq!(candidate.boundary().piston(),&expected);
    assert_eq!(candidate.boundary().step().exchanged_heat_j,0.0);
}

#[test]
fn late_refusal_cancellation_and_work_exhaustion_publish_no_owner() {
    let (a,ca)=material(0.2,2.0,0.002);let (b,cb)=material(0.4,4.0,0.001);
    let cards=[a,b];let curves=[ca,cb];let context=point(400.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=bindings(&cards,&curves,&context,&key);let resistors=resistors();
    let links=[ConductionLink {link:7,distance_a_m:0.001,distance_b_m:0.001}];
    let binding=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let thermal=initial(&curves,0.0);let electrical=circuit(&binding,&thermal);let gas=piston();
    let saved=gas.clone();let cells0=thermal.cells().clone();let circuit0=electrical.clone();
    let mut polls=0;
    let accepted=gas.material_network_candidate(&binding,&thermal,0.0005,&electrical,&[],port(0.0005),solve(),geometry(),
        &[0.0;2],||{polls+=1;false}).unwrap();
    let total=accepted.step().boundary_solves;assert!(total>1);
    let mut actual=0;
    assert!(matches!(gas.material_network_candidate(&binding,&thermal,0.0005,&electrical,&[],port(0.0005),solve(),
        geometry(),&[0.0;2],||{actual+=1;actual==polls}),Err(CircuitEnthalpyError::Cancelled)));
    assert!(gas.material_network_candidate(&binding,&thermal,0.0005,&electrical,&[],port(0.0005),
        NetworkGasSolve {maximum_boundary_solves:total-1,..solve()},geometry(),&[0.0;2],||false).is_err());
    assert!(gas.material_network_candidate(&binding,&thermal,0.0005,&electrical,&[],port(0.0006),solve(),
        geometry(),&[0.0;2],||false).is_err());
    assert!(gas.material_network_candidate(&binding,&thermal,0.0005,&electrical,&[],port(0.0005),solve(),
        GasConductionGeometry {area_m2:0.05,..geometry()},&[0.0;2],||false).is_err());
    let convecting=initial(&curves,1.0);
    assert!(gas.material_network_candidate(&binding,&convecting,0.0005,&electrical,&[],port(0.0005),solve(),
        geometry(),&[0.0;2],||false).is_err());
    assert_eq!(gas,saved);assert_eq!(thermal.cells(),&cells0);assert_eq!(electrical,circuit0);
}

#[test]
fn independent_network_and_gas_sources_are_not_counted_as_internal_exchange() {
    let (a,ca)=material(0.2,2.0,0.002);let (b,cb)=material(0.4,4.0,0.001);
    let cards=[a,b];let curves=[ca,cb];let context=point(400.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=bindings(&cards,&curves,&context,&key);let resistors=resistors();
    let links=[ConductionLink {link:7,distance_a_m:0.001,distance_b_m:0.001}];
    let binding=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let thermal=initial(&curves,0.0);let electrical=circuit(&binding,&thermal);let gas=piston();
    let destinations:Vec<_>=resistors.iter().map(|r|r.destination).collect();
    let mut drive=port(0.0005);drive.input.heat_j=0.01;
    let candidate=gas.network_candidate(&thermal,0.0005,crate::circuit_enthalpy::CircuitEnthalpyInput {
        circuit:&electrical,ports:&[],destinations:&destinations},drive,solve(),&[0.0;2],
        |state| {
            let mut coefficients=binding.resolve(state,||false)?;
            coefficients.sources_w.insert(cell(1),30.0);
            Ok::<_,MaterialNetworkError>(coefficients)
        },|_,_,_|Ok::<_,Infallible>(ConductanceSample {conductance_w_k:0.2,source:()}),||false).unwrap();
    near(candidate.step().external_network_heat_j,0.015);
    let (next,circuit,piston,report)=candidate.into_parts();
    near(report.network.non_electrical_heat_j+report.boundary.exchanged_heat_j,0.015);
    let energy=next.energy_j()-thermal.energy_j()+(circuit.charges_c()[0].powi(2)-1.0)/0.02
        +piston.gas().mass_kg()*(piston.gas().state().specific_internal_energy_j_kg()-gas.gas().state().specific_internal_energy_j_kg())
        +0.5*piston.mass_kg()*piston.velocity_m_s().powi(2)+100000.0*(piston.gas().volume_m3()-gas.gas().volume_m3());
    near(energy,0.025);
}

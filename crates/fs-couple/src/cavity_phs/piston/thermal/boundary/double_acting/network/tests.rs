//! Synthetic source cards exercise the actual finite circuit, cells and paired gas.
use super::*;
use crate::cavity_phs::piston::{GasPiston, tests::{controls, gas}};
use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, source::ThermalMaterialPoint};
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry, resolve_contact};
use crate::circuit_enthalpy::{CircuitEnthalpySolve, ResistorCell};
use crate::circuit_enthalpy::material::*;
use crate::material_circuit::MaterialCircuitSolve;
use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
use fs_material::phase::{ResolvedHeatCapacityCurve, HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_qty::{Temperature, Density};
use fs_phs::circuit::{CircuitGraph, Branch};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T",t).unwrap().with("lot",7.0).unwrap() }
fn card(k: f64, r: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name,dims,a,b) in [
        ("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS,r*0.0025*0.6,r*0.0025*1.8),
        ("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS,k,k),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,500.0,500.0),
        ("density",Density::DIMS,4000.0,4000.0)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name,dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0,a),(800.0,b)], dims },
            validity: ValidityDomain::unconstrained().with("T",200.0,800.0).with("lot",7.0,7.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic paired-network specimen".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "authored".into(), revision: 0 },claims,vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density",Density::DIMS), temperature_axis: "T",
        lower: &point(200.0), upper: &point(800.0), selection: SelectionPolicy::SingleClaimOnly },
        SolidLiquidPhase::Solid,0.0,1e-5).unwrap();
    (card,curve)
}
fn cell(i: usize) -> CellId { CellId(if i==0 {30} else {10}) }
fn cells<'a>(cards: &'a [MaterialCard;2], curves: &'a [ResolvedHeatCapacityCurve;2], p: &'a QueryPoint, key: &'a PropertyKey)
    -> [MaterialCell<'a>;2] {
    [0,1].map(|i| MaterialCell { cell: cell(i), material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: p, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((key,MaterialPropertySelection::SingleClaimOnly)) })
}
fn resistors() -> [CellResistor;2] {
    [0,1].map(|i| CellResistor { destination: ResistorCell { branch: i+1, cell: cell(i) },
        selection: MaterialPropertySelection::SingleClaimOnly, length_m: 0.01, area_m2: 0.000025 })
}
fn initial<'a>(curves: &'a [ResolvedHeatCapacityCurve;2], temperatures: [f64;2]) -> DepositionSimulation<'a> {
    let events = [0,1].map(|i| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: cell(i), segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64), mass_kg: 0.001,
        specific_h_j_kg: 500.0*(temperatures[i]-200.0), geometry: CellGeometry { start_m: [i as f64,0.0,0.0],
            end_m: [i as f64+0.01,0.0,0.0], volume_m3: 0.001/4000.0, surface_area_m2: 0.04 },
        boundary: Boundary { ambient_k: 400.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 } }) }).to_vec();
    DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],DepositionBudget {
        max_cells: 2,max_links: 1,max_events: 3,max_substeps: 100 }).unwrap()
        .advance_to(0.0,StepConfig { max_dt_s: 0.001,energy_tolerance_j: 1e-7 },||Ok(())).unwrap().0
}
fn circuit(binding: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>, charge: f64) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph { node_count: 3,branches: vec![
        (1,0,Branch::Capacitor { farads: 0.01 }),(1,2,Branch::Resistor { ohms: 0.0 }),
        (2,0,Branch::Resistor { ohms: 0.0 })],transformers: vec![] },
        binding.resolve(thermal,||false).unwrap().resistors,vec![],vec![charge],&[]).unwrap()
}
fn solve() -> NetworkGasSolve { NetworkGasSolve { network: CircuitEnthalpySolve {
    maximum_cells: 2,maximum_circuit_trials: 300,maximum_newton_iterations: 32,difference_step_j: 1e-4,
    cell_energy_tolerance_j: 1e-7,thermal_energy_tolerance_j: 1e-7,total_energy_tolerance_j: 1e-7,
    circuit: MaterialCircuitSolve { energy_tolerance_j: 1e-7,power_tolerance_w: 1e-7 },
    newton: NewtonKrylovConfig { absolute_tolerance: 1e-8,relative_tolerance: 0.0,forcing_minimum: 1e-4,
        linear_restart: 4,..NewtonKrylovConfig::default() } },maximum_boundary_solves: 300,total_energy_tolerance_j: 1e-6 } }
fn drive(dt: f64) -> DoubleNetworkDrive {
    DoubleNetworkDrive { cells: [cell(0),cell(1)],
        input: DoublePistonInput { duration_s: dt,heat_j: [0.0;2],external_force_n: 0.0 },
        piston: DoublePistonSolve { piston: controls(),opposing_gas: controls().gas,volume_tolerance_m3: 1e-15 },
        exchange: DoubleHeatExchangeSolve { faces: [HeatExchangeSolve { heat_bracket_j: [-1.0,1.0],
            source_tolerance_j: 1e-10,energy_tolerance_j: 1e-7,maximum_iterations: 128 };2],maximum_joint_trials: 20000 } }
}
fn geometry() -> [GasConductionGeometry;2] { [GasConductionGeometry {
    area_m2: 0.01,gas_distance_m: 0.001,cell_distance_m: 0.02 };2] }
fn piston() -> DoubleActingGasPiston {
    DoubleActingGasPiston::new(GasPiston::new(gas(0.75,0.001),1.0,0.01,0.0).unwrap(),gas(0.25,0.001),0.01).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs()<2e-6,"{a} != {b}"); }

#[test]
fn finite_electrical_store_drives_both_gases_and_material_swap_changes_current() {
    let run = |k| {
        let (a,ca)=card(k,2.0); let (b,cb)=card(0.4,4.0);
        let cards=[a,b]; let curves=[ca,cb]; let point=point(400.0);
        let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
        let cells=cells(&cards,&curves,&point,&key); let resistors=resistors();
        let binding=MaterialNetwork { cells: &cells,resistors: &resistors,links: &[] };
        let mut thermal=initial(&curves,[400.0;2]); let mut circuit=circuit(&binding,&thermal,1.0);
        let mut piston=piston(); let old=piston.clone(); let initial_heat=thermal.energy_j();
        for i in 1..=16 {
            let target=f64::from(i)*0.0005; let dt=target-thermal.time_s();
            let candidate=piston.material_network_candidate(&binding,&thermal,target,&circuit,&[],
                drive(dt),solve(),geometry(),&[0.0;2],||false).unwrap();
            let (t,c,p,r)=candidate.into_parts();
            for i in 0..2 {
                let material=r.network.source.cells.iter().find(|s|s.cell==r.cells[i]).unwrap();
                assert_eq!(material.conductivity.as_ref(),Some(&r.interfaces[i]));
                near(r.boundary.boundary_temperature_k[i],t.temperature_k(r.cells[i]).unwrap());
                assert!(r.boundary.source_residual_j[i].abs()<=drive(dt).exchange.faces[i].source_tolerance_j);
            }
            let gas_energy=p.primary().gas().mass_kg()*(p.primary().gas().state().specific_internal_energy_j_kg()
                -old.primary().gas().state().specific_internal_energy_j_kg())
                +p.opposing().mass_kg()*(p.opposing().state().specific_internal_energy_j_kg()
                -old.opposing().state().specific_internal_energy_j_kg());
            near(t.energy_j()-initial_heat+r.network.electrical.stored_after_j-50.0+gas_energy
                +0.5*p.primary().mass_kg()*p.primary().velocity_m_s().powi(2),0.0);
            near(r.coupled.external_network_heat_j,0.0);
            assert!(r.coupled.total_energy_residual_j.abs()<=solve().total_energy_tolerance_j);
            assert!(r.coupled.boundary_solves>0 && r.coupled.exchange_trials>=r.coupled.boundary_solves);
            thermal=t; circuit=c; piston=p;
        }
        (thermal.temperature_k(cell(0)).unwrap(),circuit.charges_c()[0],piston)
    };
    let a=run(0.2); let b=run(2.0);
    assert!((a.0-b.0).abs()>1e-6);
    assert!((a.1-b.1).abs()>1e-10);
    assert!((a.2.primary().velocity_m_s()-b.2.primary().velocity_m_s()).abs()>1e-9);
    assert_ne!(a.2.opposing().state().state().sound_speed,b.2.opposing().state().state().sound_speed);
}

#[test]
fn both_faces_debit_one_existing_cell_once_and_refuse_overlapping_coverage() {
    let (a,ca)=card(0.2,2.0); let (b,cb)=card(0.4,4.0);
    let cards=[a,b]; let curves=[ca,cb]; let point=point(400.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,&curves,&point,&key); let resistors=resistors();
    let binding=MaterialNetwork { cells: &cells,resistors: &resistors,links: &[] };
    let thermal=initial(&curves,[500.0,400.0]); let circuit=circuit(&binding,&thermal,0.0); let old=piston();
    let drive=DoubleNetworkDrive { cells: [cell(0);2],..drive(0.0005) };
    let c=old.material_network_candidate(&binding,&thermal,0.0005,&circuit,&[],drive,solve(),geometry(),&[0.0;2],||false).unwrap();
    let (next,_,_,r)=c.into_parts();
    assert!(r.boundary.exchanged_heat_j.iter().all(|q|*q>0.0));
    let net=r.boundary.exchanged_heat_j.iter().sum::<f64>();
    near(next.energy_j()-thermal.energy_j(),-net);
    near(r.debited_heat_j,net);
    assert_eq!(next.cells()[&cell(1)],thermal.cells()[&cell(1)]);
    let mut patches=geometry(); for p in &mut patches { p.area_m2=0.03; }
    assert!(old.material_network_candidate(&binding,&thermal,0.0005,&circuit,&[],drive,solve(),patches,&[0.0;2],||false).is_err());
}

#[test]
fn independent_external_supply_and_final_cancellation_keep_all_histories_atomic() {
    let (a,ca)=card(0.2,2.0); let (b,cb)=card(0.4,4.0);
    let cards=[a,b]; let curves=[ca,cb]; let point=point(400.0);
    let key=PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=cells(&cards,&curves,&point,&key); let resistors=resistors();
    let binding=MaterialNetwork { cells: &cells,resistors: &resistors,links: &[] };
    let thermal=initial(&curves,[500.0,300.0]); let old_cells=thermal.cells().clone();
    let circuit=circuit(&binding,&thermal,0.1); let old=piston(); let saved=old.clone();
    let destinations:Vec<_>=resistors.iter().map(|r|r.destination).collect();
    let resolve=|s:&DepositionSimulation<'_>| {
        let mut c=binding.resolve(s,||false)?; c.sources_w.insert(cell(0),10.0);
        Ok::<_,MaterialNetworkError>(c)
    };
    let contacts=|s:&DepositionSimulation<'_>,m:&MaterialNetworkSample,g:&[ResolvedThermochemicalGasState;2]|
        -> Result<[ConductanceSample<fs_material::state_point::ResolvedMaterialStatePoint>;2],MaterialNetworkError> { Ok([
        resolve_contact(&binding,cell(0),s,m,&g[0],geometry()[0])?,
        resolve_contact(&binding,cell(1),s,m,&g[1],geometry()[1])?,
    ]) };
    let electrical=CircuitEnthalpyInput { circuit: &circuit,ports: &[],destinations: &destinations };
    let mut polls=0;
    let candidate=old.network_candidate(&thermal,0.0005,electrical,drive(0.0005),solve(),&[0.0;2],resolve,contacts,
        || { polls+=1; false }).unwrap();
    near(candidate.step().external_network_heat_j,0.005);
    let mut count=0;
    assert!(matches!(old.network_candidate(&thermal,0.0005,electrical,drive(0.0005),solve(),&[0.0;2],resolve,contacts,
        || { count+=1; count==polls }),Err(CircuitEnthalpyError::Cancelled)));
    assert!(old.network_candidate(&thermal,0.0005,electrical,drive(0.0005),NetworkGasSolve {
        maximum_boundary_solves:1,..solve() },&[0.0;2],resolve,contacts,||false).is_err());
    assert_eq!(old,saved); assert_eq!(thermal.cells(),&old_cells);
    assert_eq!(thermal.time_s(),0.0); assert_eq!(circuit.charges_c(),&[0.1]);
}

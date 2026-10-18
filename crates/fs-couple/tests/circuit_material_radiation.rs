//! Radiation-only data changes propagate through real circuit/thermal evolution.
//! Reuse the source/trace fixtures used by the radiation owner's own tests.
#[path="../../fs-conduction/tests/radiation_material/common.rs"] mod common;
use common::*;
use fs_conduction::transient::deposition::radiation::*;
use fs_couple::circuit_enthalpy::{self, CircuitEnthalpyInput, CircuitEnthalpySolve,
    CircuitThermalCoefficients, ResistorCell};
use fs_couple::circuit_enthalpy::material::{MaterialCell, MaterialNetwork, CellResistor};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS;
use fs_thermal::enthalpy::*;
use fs_matdb::*;
use fs_qty::Temperature;
use fs_phs::circuit::{CircuitGraph,Branch};
use fs_solver::NewtonKrylovConfig;

fn electrical_card(e:f64,r:f64)->MaterialCard {
    let card=common::card(e,e);let mut claims=card.claims().clone();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new("electrical_resistivity",ELECTRICAL_RESISTIVITY_DIMS),
        value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(300.0,r),(500.0,2.0*r)],dims:ELECTRICAL_RESISTIVITY_DIMS},
        validity:fs_evidence::ValidityDomain::unconstrained().with("T",300.0,500.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance {source:"synthetic temperature-dependent resistor".into(),
            license:"CC0-1.0".into(),artifact:None},
    }).unwrap();
    MaterialCard::assemble(card.id().clone(),claims,vec![]).unwrap()
}
#[test]
fn emissivity_only_changes_both_body_temperatures_and_finite_electrical_history() {
    let run=|emissivity_a:f64| {
        let cards=[electrical_card(emissivity_a,1.0),electrical_card(0.6,2.0)];
        let curves=[curve(&cards[0],100.0),curve(&cards[1],100.0)];let context=context();
        let ids=[CellId(30),CellId(10)];let destinations=[0,1].map(|i|ResistorCell {branch:i+1,cell:ids[i]});
        let cells=[0,1].map(|i|MaterialCell {cell:ids[i],material_id:MaterialId(i as u64),curve:&curves[i],
            point:binding(&cards[i],&context),conductivity:None});
        let resistors=[0,1].map(|i|CellResistor {destination:destinations[i],selection:MaterialPropertySelection::SingleClaimOnly,
            length_m:1.0,area_m2:1.0});
        let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
        let surfaces=[0,1].map(|i|RadiatingCell {surface:i,cell:ids[i],material_id:MaterialId(i as u64),curve:&curves[i],
            material:binding(&cards[i],&context),selection:MaterialPropertySelection::SingleClaimOnly});
        let enclosure=enclosure([emissivity(&cards[0],&curves[0],150.0),emissivity(&cards[1],&curves[1],-50.0)]);
        let radiation=MaterialRadiation {enclosure:&enclosure,cells:&surfaces,power_tolerance_w:1e-8,radiosity_tolerance_w_m2:1e-8};
        let events=[0,1].map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id:ids[i],segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:1.0,
            specific_h_j_kg:if i==0 {150.0} else {-50.0},geometry:CellGeometry {
                start_m:[i as f64,0.0,0.0],end_m:[i as f64+1.0,0.0,0.0],volume_m3:1.0,surface_area_m2:1.5},
            boundary:Boundary {ambient_k:350.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0},
        })}).to_vec();
        let mut thermal=DepositionSimulation::new(events,&[(MaterialId(0),&curves[0] as &dyn ThermalLaw),
            (MaterialId(1),&curves[1] as &dyn ThermalLaw)],DepositionBudget {
            max_cells:2,max_links:0,max_events:2,max_substeps:100}).unwrap()
            .advance_to(0.0,StepConfig {max_dt_s:0.001,energy_tolerance_j:1e-7},||Ok(())).unwrap().0;
        let mut circuit=MaterialCircuit::new(CircuitGraph {node_count:3,branches:vec![
            (1,0,Branch::Capacitor {farads:0.01}),(1,2,Branch::Resistor {ohms:0.0}),
            (2,0,Branch::Resistor {ohms:0.0})],transformers:vec![]},
            network.resolve(&thermal,||false).unwrap().resistors,vec![],vec![1.0],&[]).unwrap();
        let controls=CircuitEnthalpySolve {maximum_cells:2,maximum_circuit_trials:300,maximum_newton_iterations:32,
            difference_step_j:1e-4,cell_energy_tolerance_j:1e-7,thermal_energy_tolerance_j:1e-7,total_energy_tolerance_j:1e-7,
            circuit:MaterialCircuitSolve {energy_tolerance_j:1e-7,power_tolerance_w:1e-7},
            newton:NewtonKrylovConfig {absolute_tolerance:1e-9,relative_tolerance:0.0,forcing_minimum:1e-4,
                linear_restart:4,..NewtonKrylovConfig::default()}};
        let gate=fs_exec::CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx=fs_exec::Cx::new(&gate,arena,fs_exec::StreamKey {seed:73,kernel_id:51,tile:0,iteration:0},
                fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
            for n in 1..=40 {
                let candidate=circuit_enthalpy::step_candidate(&thermal,n as f64*0.0005,
                    CircuitEnthalpyInput {circuit:&circuit,ports:&[],destinations:&destinations},controls,&[0.0;2],|state| {
                        let mut coefficients=network.resolve(state,||false).map_err(|e|e.to_string())?;
                        let emission=radiation.resolve(&cx,state).map_err(|e|e.to_string())?;
                        emission.add_to(&mut coefficients.sources_w).map_err(|e|e.to_string())?;
                        Ok::<_,String>(CircuitThermalCoefficients {resistors:coefficients.resistors,sources_w:coefficients.sources_w,
                            conductances_w_k:coefficients.conductances_w_k,source:(coefficients.source,emission)})
                    },||gate.is_requested()).unwrap();
                let report=candidate.report();
                near(candidate.thermal().energy_j()+report.electrical.stored_after_j,150.0);
                for sample in report.source.1.surfaces() {
                    let original=report.source.0.cells.iter().find(|c|c.cell==sample.cell).unwrap();
                    assert_eq!(sample.phase,original.phase);
                }
                let (t,c,_)=candidate.into_parts();thermal=t;circuit=c;
            }
        });
        (thermal.temperature_k(ids[0]).unwrap(),thermal.temperature_k(ids[1]).unwrap(),circuit.charges_c()[0])
    };
    let low=run(0.2);let high=run(0.8);
    assert!((low.0-high.0).abs()>0.1 && (low.1-high.1).abs()>0.1);
    assert!((low.2-high.2).abs()>1e-6,"radiation must feed back into resistance and finite electrical storage");
}

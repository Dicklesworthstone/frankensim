//! G0/G1/G3: material atmosphere -> actual finite-cell heat, without another stepper.
use std::collections::BTreeMap;
use fs_airflow::graph::thermal::medium::convection::{
    ForcedConvectionSettings, evaluate_forced_from_material, deposition::*,
};
use fs_convection::CorrelationId;
use fs_conduction::transient::deposition::*;
use fs_conduction::transient::deposition::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId};
use fs_qty::{Length, Temperature};
use fs_qty::semantic::{Composition, CompositionBasis};

const IDS: [CellId;2] = [CellId(30), CellId(10)];
const AREA: f64 = 0.01;
const MASS: f64 = 0.001;
const CP: f64 = 1000.0;
fn near(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a} != {b}"); }
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),
        phase:"solid".into(),process:"convection test, not measured data".into(),revision:0},
        ClaimSet::new(),vec![]).unwrap();
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),[280.0,320.0].map(|t|
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:CP*(t-300.0),temperature_k:t,
            liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0}).to_vec()).unwrap()
}
fn initial(curve:&EquilibriumEnthalpyPhaseCurve, temperatures:[f64;2], htc:f64,
    coverage:f64) -> DepositionSimulation<'_> {
    let mut events=IDS.into_iter().enumerate().map(|(i,id)|TimedEvent {time_s:0.0,
        event:Event::Birth(DepositedCell {id,segment_id:i as u64,body_id:i as u64,
            material_id:MaterialId(0),mass_kg:MASS,specific_h_j_kg:CP*(temperatures[i]-300.0),
            geometry:CellGeometry {start_m:[0.0,i as f64,0.0],end_m:[0.1,i as f64,0.0],
                volume_m3:1e-6,surface_area_m2:AREA},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:htc,reservoir:None,source_w:0.0}})}).collect::<Vec<_>>();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:5,a:IDS[0],b:IDS[1],area_m2:coverage,conductance_w_k:0.0})});
    DepositionSimulation::new(events,&[(MaterialId(0),curve as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:1.0,energy_tolerance_j:1e-10},||Ok(())).unwrap().0
}
fn binding(fluid:&dyn FluidTransport,cell:CellId,reference:f64)->ConvectiveCell<'_> {
    ConvectiveCell {cell,fluid,reference_temperature:Temperature::new(reference),
        characteristic_length:Length::new(0.1),correlation:CorrelationId::FlatPlateLaminarAverage,
        speed_m_s:1.0,settings:ForcedConvectionSettings::default()}
}

#[test]
fn changing_only_pressure_or_composition_changes_actual_exposed_cell_power() {
    let chart=curve();let heat=initial(&chart,[310.0,310.0],0.0,0.002);
    let spec=GasSpec::dry_air_ussa1976();
    let low=GasState::try_new(&spec,300.0,100000.0).unwrap();
    let high=GasState::try_new(&spec,300.0,200000.0).unwrap();
    let run=|fluid:&dyn FluidTransport|evaluate_cell_convection(&heat,&[binding(fluid,IDS[0],300.0)],1,||false).unwrap();
    let a=run(&low);let b=run(&high);
    near(b.sources_w()[&IDS[0]]/a.sources_w()[&IDS[0]],2.0_f64.sqrt(),1e-10);
    let expected=evaluate_forced_from_material(CorrelationId::FlatPlateLaminarAverage,&low,
        1.0,Length::new(0.1),ForcedConvectionSettings::default()).unwrap();
    near(a.sources_w()[&IDS[0]],expected.coefficient().value.value()*0.008*(-10.0),1e-12);
    near(a.samples()[0].exposed_area_m2,0.008,1e-15);
    let mix=GasState::try_new_mixture(&[spec,GasSpec::water_vapor_nist()],
        &Composition::new(CompositionBasis::MoleFraction,vec![0.98,0.02]).unwrap(),
        300.0,100000.0,Some(1)).unwrap();
    assert!((run(&mix).sources_w()[&IDS[0]]-a.sources_w()[&IDS[0]]).abs()>1e-8);
    assert_eq!(heat.time_s(),0.0);
}

#[test]
fn implicit_heating_and_cooling_use_actual_endpoint_temperatures_and_original_storage() {
    let chart=curve();let heat=initial(&chart,[310.0,290.0],0.0,0.002);
    let fluid=FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap();
    let cells=IDS.map(|id|binding(&fluid,id,300.0));
    let source=evaluate_cell_convection(&heat,&cells,2,||false).unwrap();
    let dt=0.2;
    let step=BackwardEulerStep::new(&heat,dt,||Ok(())).unwrap();
    let increments=step.cell_ids().iter().map(|&id| {
        let sample=source.samples().iter().find(|s|s.cell==id).unwrap();
        let g=sample.coefficient.value.value()*sample.exposed_area_m2;
        let old=heat.temperature_k(id).unwrap();
        MASS*CP*((300.0+(old-300.0)/(1.0+dt*g/(MASS*CP)))-old)
    }).collect::<Vec<_>>();
    let trial=step.evaluate(&increments,|endpoint| {
        let sources=evaluate_cell_convection(endpoint,&cells,2,||false)?;
        Ok::<_,ConvectionSourceError>(StepCoefficients {sources_w:sources.sources_w().clone(),
            conductances_w_k:BTreeMap::new(),payload:sources})
    },||Ok(())).unwrap();
    let (next,actual,audit)=trial.accept(1e-10,1e-10,||Ok(())).unwrap();
    for sample in actual.samples() {
        near(sample.wall_temperature_k,next.temperature_k(sample.cell).unwrap(),1e-12);
        assert!(sample.heat_into_cell_w.abs()<source.sources_w()[&sample.cell].abs());
    }
    assert!(actual.sources_w()[&IDS[0]]<0.0 && actual.sources_w()[&IDS[1]]>0.0);
    near(audit.external_heat_j,dt*actual.total_heat_into_cells_w(),1e-12);
    near(next.energy_j()-heat.energy_j(),audit.external_heat_j,1e-10);
    assert_eq!(heat.temperature_k(IDS[0]).unwrap(),310.0);
    assert_eq!(next.cells()[&IDS[0]].boundary.htc_w_m2_k,0.0);
    assert_eq!(next.cells()[&IDS[0]].boundary.source_w,0.0);
}

#[test]
fn liquid_and_counted_fluid_reuse_one_original_transport_snapshot_per_cell() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self)->Result<FluidTransportState,fs_material::MaterialError> {
            self.0.set(self.0.get()+1);
            FluidTransportState::try_new(1000.0,4000.0,0.5,0.001)
        }
    }
    let chart=curve();let heat=initial(&chart,[300.0;2],0.0,0.0);
    let fluid=Counted(std::cell::Cell::new(0));
    let result=evaluate_cell_convection(&heat,&[binding(&fluid,IDS[0],300.0)],1,||false).unwrap();
    assert_eq!(fluid.0.get(),1);assert_eq!(result.sources_w()[&IDS[0]],0.0);
    near(result.samples()[0].nusselt.groups()["Re"],100000.0,1e-8);
    near(result.samples()[0].nusselt.groups()["Pr"],8.0,1e-12);
    let mut sources=BTreeMap::from([(IDS[0],2.0)]);
    result.add_to(&mut sources,||false).unwrap();
    assert_eq!(fluid.0.get(),1);assert_eq!(sources[&IDS[0]],2.0);
}

#[test]
fn existing_contact_coverage_and_stable_cell_order_remain_authoritative() {
    let chart=curve();let fluid=FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap();
    let uncovered=initial(&chart,[310.0;2],0.0,0.0);
    let covered=initial(&chart,[310.0;2],0.0,AREA/2.0);
    let bindings=IDS.map(|id|binding(&fluid,id,300.0));
    let a=evaluate_cell_convection(&uncovered,&bindings,2,||false).unwrap();
    let b=evaluate_cell_convection(&covered,&bindings,2,||false).unwrap();
    let c=evaluate_cell_convection(&covered,&[bindings[1],bindings[0]],2,||false).unwrap();
    assert_eq!(b.sources_w(),c.sources_w());assert_eq!(b.samples()[0].cell,CellId(10));
    for id in IDS {near(b.sources_w()[&id],0.5*a.sources_w()[&id],1e-12);}
}

#[test]
fn duplicate_robin_missing_cells_invalid_reservoir_and_work_limits_refuse() {
    let chart=curve();let heat=initial(&chart,[300.0;2],0.0,0.0);
    let robin=initial(&chart,[300.0;2],5.0,0.0);
    let fluid=FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap();
    let b=binding(&fluid,IDS[0],290.0);
    assert!(evaluate_cell_convection(&robin,&[b],1,||false).is_err());
    assert!(evaluate_cell_convection(&heat,&[b,b],2,||false).is_err());
    assert!(evaluate_cell_convection(&heat,&[b],0,||false).is_err());
    assert!(evaluate_cell_convection(&heat,&[ConvectiveCell {cell:CellId(999),..b}],1,||false).is_err());
    assert!(evaluate_cell_convection(&heat,&[ConvectiveCell {reference_temperature:Temperature::new(0.0),..b}],1,||false).is_err());
    assert!(evaluate_cell_convection(&heat,&[ConvectiveCell {speed_m_s:-1.0,..b}],1,||false).is_err());
}

#[test]
fn final_cancellation_and_source_overflow_do_not_publish_partial_heat() {
    let chart=curve();let heat=initial(&chart,[310.0;2],0.0,0.0);
    let before=heat.cells().clone();
    let fluid=FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap();
    let cells=IDS.map(|id|binding(&fluid,id,300.0));
    let mut count=0;
    let sources=evaluate_cell_convection(&heat,&cells,2,||{count+=1;false}).unwrap();
    let mut at=0;
    assert!(matches!(evaluate_cell_convection(&heat,&cells,2,||{at+=1;at==count}),
        Err(ConvectionSourceError::Cancelled)));
    let mut destination=BTreeMap::from([(IDS[0],f64::NAN),(IDS[1],2.0)]);
    assert!(sources.add_to(&mut destination,||false).is_err());
    assert!(destination[&IDS[0]].is_nan());assert_eq!(destination[&IDS[1]],2.0);
    let mut destination=BTreeMap::from([(IDS[0],1.0),(IDS[1],2.0)]);
    let saved=destination.clone();let mut polls=0;
    sources.add_to(&mut destination,||{polls+=1;false}).unwrap();
    let mut destination=saved.clone();let mut at=0;
    assert!(matches!(sources.add_to(&mut destination,||{at+=1;at==polls}),Err(ConvectionSourceError::Cancelled)));
    assert_eq!(destination,saved);assert_eq!(heat.cells(),&before);
}

#[test]
fn dittus_boelter_direction_follows_the_solved_wall_not_a_stale_flag() {
    let chart=curve();let heat=initial(&chart,[310.0,290.0],0.0,0.0);
    let fluid=FluidTransportState::try_new(1000.0,4000.0,0.5,0.001).unwrap();
    let cells=IDS.map(|id|ConvectiveCell {
        correlation:CorrelationId::DittusBoelter,
        settings:ForcedConvectionSettings {length_over_hydraulic_diameter:Some(20.0),
            ..ForcedConvectionSettings::default()},..binding(&fluid,id,300.0)});
    let source=evaluate_cell_convection(&heat,&cells,2,||false).unwrap();
    let hot=source.samples().iter().find(|s|s.cell==IDS[0]).unwrap();
    let cold=source.samples().iter().find(|s|s.cell==IDS[1]).unwrap();
    near(hot.coefficient.value.value()/cold.coefficient.value.value(),8.0_f64.powf(0.1),1e-10);
    assert!(hot.heat_into_cell_w<0.0 && cold.heat_into_cell_w>0.0);
}

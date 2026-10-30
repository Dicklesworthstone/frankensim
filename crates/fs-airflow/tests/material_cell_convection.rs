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

#[path = "material_cell_convection/thermal_source.rs"]
mod thermal_source;

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
        let g=sample.evaluation.as_ref().unwrap().coefficient.value.value()*sample.exposed_area_m2;
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
    near(result.samples()[0].evaluation.as_ref().unwrap().nusselt.groups()["Re"],100000.0,1e-8);
    near(result.samples()[0].evaluation.as_ref().unwrap().nusselt.groups()["Pr"],8.0,1e-12);
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
    near(hot.evaluation.as_ref().unwrap().coefficient.value.value()
        /cold.evaluation.as_ref().unwrap().coefficient.value.value(),8.0_f64.powf(0.1),1e-10);
    assert!(hot.heat_into_cell_w<0.0 && cold.heat_into_cell_w>0.0);
}

fn natural(fluid:&dyn FluidTransport,cell:CellId)->MaterialConvectiveCell<'_> {
    MaterialConvectiveCell {cell,fluid,reference_temperature:Temperature::new(300.0),
        characteristic_length:Length::new(0.1),
        law:CellConvectionLaw::VerticalPlate {gravity_m_s2:9.81,expansion_per_k:0.002}}
}

// Independent Churchill-Chu algebra on explicitly declared source properties.
fn natural_h(fluid:FluidTransportState,delta:f64,beta:f64)->f64 {
    let rho=fluid.density_kg_m3();let cp=fluid.specific_heat_j_kg_k();
    let k=fluid.thermal_conductivity_w_m_k();let mu=fluid.dynamic_viscosity_pa_s();
    let ra=(9.81*beta*delta*0.1_f64.powi(3)/((mu/rho)*(k/(rho*cp)))).abs();
    let pr=mu*cp/k;
    (0.825+0.387*ra.powf(1.0/6.0)/(1.0+(0.492/pr).powf(9.0/16.0)).powf(8.0/27.0)).powi(2)*k/0.1
}

fn natural_endpoint(fluid:FluidTransportState,old:f64,heater:f64,area:f64,dt:f64)->f64 {
    let mut lower=old.min(300.0);
    let mut upper=old.max(300.0)+dt*heater/(MASS*CP);
    for _ in 0..80 {
        let t=0.5*lower+0.5*upper;
        let power=heater+area*natural_h(fluid,t-300.0,0.002)*(300.0-t);
        let residual=MASS*CP*(t-old)-dt*power;
        if residual>0.0 {upper=t;} else {lower=t;}
    }
    0.5*lower+0.5*upper
}

#[test]
fn natural_hot_cold_and_material_changes_follow_independent_rayleigh_heat() {
    let chart=curve();let heat=initial(&chart,[310.0,290.0],0.0,0.002);
    let first=FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap();
    let second=FluidTransportState::try_new(2.0,1200.0,0.03,3e-5).unwrap();
    let run=|fluid:&FluidTransportState| {
        let cells=IDS.map(|id|natural(fluid,id));
        evaluate_material_cell_convection(&heat,&cells,2,||false).unwrap()
    };
    let a=run(&first);let b=run(&second);
    for (result,fluid) in [(&a,first),(&b,second)] {
        for sample in result.samples() {
            let delta=sample.wall_temperature_k-300.0;
            let h=natural_h(fluid,delta,0.002);
            near(sample.evaluation.as_ref().unwrap().coefficient.value.value(),h,1e-10);
            near(sample.heat_into_cell_w,-0.008*h*delta,1e-10);
            assert_eq!(sample.transport,fluid);
        }
        near(result.sources_w()[&IDS[0]],-result.sources_w()[&IDS[1]],1e-12);
        assert!(result.sources_w()[&IDS[0]]<0.0 && result.sources_w()[&IDS[1]]>0.0);
    }
    assert!((a.sources_w()[&IDS[0]]-b.sources_w()[&IDS[0]]).abs()>1e-3);
    let reversed=MaterialConvectiveCell {law:CellConvectionLaw::VerticalPlate {
        gravity_m_s2:9.81,expansion_per_k:-0.002},..natural(&first,IDS[0])};
    let reverse=evaluate_material_cell_convection(&heat,&[reversed],1,||false).unwrap();
    near(reverse.sources_w()[&IDS[0]],a.sources_w()[&IDS[0]],1e-12);
}

#[test]
fn natural_boundaries_close_original_backward_euler_storage_at_the_actual_endpoint() {
    let chart=curve();let heat=initial(&chart,[310.0,290.0],0.0,0.002);
    let fluid=FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap();
    let cells=IDS.map(|id|natural(&fluid,id));
    let before=evaluate_material_cell_convection(&heat,&cells,2,||false).unwrap();
    let dt=0.2;let step=BackwardEulerStep::new(&heat,dt,||Ok(())).unwrap();
    let increments=step.cell_ids().iter().map(|&id| {
        let old=heat.temperature_k(id).unwrap();
        MASS*CP*(natural_endpoint(fluid,old,0.0,0.008,dt)-old)
    }).collect::<Vec<_>>();
    let trial=step.evaluate(&increments,|endpoint| {
        let sources=evaluate_material_cell_convection(endpoint,&cells,2,||false)?;
        Ok::<_,ConvectionSourceError>(StepCoefficients {sources_w:sources.sources_w().clone(),
            conductances_w_k:BTreeMap::new(),payload:sources})
    },||Ok(())).unwrap();
    let (next,actual,audit)=trial.accept(1e-9,1e-9,||Ok(())).unwrap();
    for sample in actual.samples() {
        let old=heat.temperature_k(sample.cell).unwrap();
        near(next.temperature_k(sample.cell).unwrap(),natural_endpoint(fluid,old,0.0,0.008,dt),1e-10);
        near(sample.wall_temperature_k,next.temperature_k(sample.cell).unwrap(),1e-12);
        let original=before.samples().iter().find(|s|s.cell==sample.cell).unwrap();
        assert!(sample.evaluation.as_ref().unwrap().coefficient.value.value()
            <original.evaluation.as_ref().unwrap().coefficient.value.value());
        near(MASS*CP*(sample.wall_temperature_k-old),dt*sample.heat_into_cell_w,1e-9);
    }
    near(next.energy_j()-heat.energy_j(),audit.external_heat_j,1e-9);
    near(audit.external_heat_j,dt*actual.total_heat_into_cells_w(),1e-12);
    assert_eq!(heat.time_s(),0.0);
}

#[test]
fn natural_equilibrium_retains_source_admission_and_can_start_a_heated_step() {
    struct Counted {calls:std::cell::Cell<usize>,state:FluidTransportState}
    impl FluidTransport for Counted {
        fn transport_state(&self)->Result<FluidTransportState,fs_material::MaterialError> {
            self.calls.set(self.calls.get()+1);Ok(self.state)
        }
    }
    let chart=curve();let heat=initial(&chart,[300.0;2],0.0,0.0);
    let fluid=Counted {calls:std::cell::Cell::new(0),
        state:FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap()};
    let cells=IDS.map(|id|natural(&fluid,id));
    let zero=evaluate_material_cell_convection(&heat,&cells,2,||false).unwrap();
    assert_eq!(fluid.calls.get(),2);
    assert_eq!(zero.total_heat_into_cells_w(),0.0);
    assert!(zero.samples().iter().all(|s|s.evaluation.is_none() && s.heat_into_cell_w==0.0));
    for cell in [MaterialConvectiveCell {characteristic_length:Length::new(0.0),..cells[0]},
        MaterialConvectiveCell {law:CellConvectionLaw::VerticalPlate {gravity_m_s2:-1.0,
            expansion_per_k:0.002},..cells[0]},
        MaterialConvectiveCell {law:CellConvectionLaw::VerticalPlate {gravity_m_s2:9.81,
            expansion_per_k:f64::NAN},..cells[0]}]
    {assert!(evaluate_material_cell_convection(&heat,&[cell],1,||false).is_err());}
    struct Missing;
    impl FluidTransport for Missing {
        fn transport_state(&self)->Result<FluidTransportState,fs_material::MaterialError> {
            Err(fs_material::MaterialError::State {what:"missing source".into()})
        }
    }
    assert!(evaluate_material_cell_convection(&heat,&[natural(&Missing,IDS[0])],1,||false).is_err());
    let dt=0.2;let heater=1.0;
    let step=BackwardEulerStep::new(&heat,dt,||Ok(())).unwrap();
    let expected=natural_endpoint(fluid.state,300.0,heater,AREA,dt);
    let increments=vec![MASS*CP*(expected-300.0);2];
    let coefficients=|endpoint:&DepositionSimulation<'_>| {
        let sources=evaluate_material_cell_convection(endpoint,&cells,2,||false)?;
        let mut powers=BTreeMap::from(IDS.map(|id|(id,heater)));
        sources.add_to(&mut powers,||false)?;
        Ok::<_,ConvectionSourceError>(StepCoefficients {sources_w:powers,
            conductances_w_k:BTreeMap::new(),payload:sources})
    };
    // The exact initial guess is a valid residual query, not a fabricated
    // positive temperature offset chosen just to make the correlation evaluate.
    step.evaluate(&[0.0;2],coefficients,||Ok(())).unwrap();
    let trial=step.evaluate(&increments,coefficients,||Ok(())).unwrap();
    let (next,sources,audit)=trial.accept(1e-9,1e-9,||Ok(())).unwrap();
    assert!(next.temperature_k(IDS[0]).unwrap()>300.0);
    assert!(sources.samples().iter().all(|s|s.evaluation.is_some() && s.heat_into_cell_w<0.0));
    near(next.temperature_k(IDS[0]).unwrap(),expected,1e-10);
    near(next.energy_j()-heat.energy_j(),dt*(2.0*heater+sources.total_heat_into_cells_w()),1e-9);
    near(audit.external_heat_j,next.energy_j()-heat.energy_j(),1e-9);
}

#[test]
fn mixed_cell_laws_preserve_forced_values_boundary_ownership_and_rollback() {
    let chart=curve();let heat=initial(&chart,[310.0,290.0],0.0,0.002);
    let fluid=FluidTransportState::try_new(1.0,1000.0,0.02,2e-5).unwrap();
    let forced=binding(&fluid,IDS[0],300.0);
    let cells=[MaterialConvectiveCell::from(forced),natural(&fluid,IDS[1])];
    let mixed=evaluate_material_cell_convection(&heat,&cells,2,||false).unwrap();
    let old=evaluate_cell_convection(&heat,&[forced],1,||false).unwrap();
    assert_eq!(mixed.sources_w()[&IDS[0]],old.sources_w()[&IDS[0]]);
    let reordered=evaluate_material_cell_convection(&heat,&[cells[1],cells[0]],2,||false).unwrap();
    assert_eq!(mixed.sources_w(),reordered.sources_w());
    assert_eq!(mixed.samples()[0].cell,IDS[1]);
    let before=heat.cells().clone();
    assert!(evaluate_material_cell_convection(&heat,&[cells[0],natural(&fluid,IDS[0])],2,||false).is_err());
    let owned=initial(&chart,[310.0,290.0],5.0,0.0);
    assert!(evaluate_material_cell_convection(&owned,&cells,2,||false).is_err());
    for law in [CellConvectionLaw::VerticalPlate {gravity_m_s2:0.0,expansion_per_k:0.002},
        CellConvectionLaw::VerticalPlate {gravity_m_s2:9.81,expansion_per_k:0.0},
        CellConvectionLaw::VerticalPlate {gravity_m_s2:9.81,expansion_per_k:1e30}]
    {
        let invalid=MaterialConvectiveCell {law,..cells[1]};
        assert!(evaluate_material_cell_convection(&heat,&[cells[0],invalid],2,||false).is_err());
    }
    let mut count=0;
    evaluate_material_cell_convection(&heat,&cells,2,||{count+=1;false}).unwrap();
    let mut at=0;
    assert!(matches!(evaluate_material_cell_convection(&heat,&cells,2,||{at+=1;at==count}),
        Err(ConvectionSourceError::Cancelled)));
    let mut destination=BTreeMap::from([(IDS[0],1.0),(IDS[1],2.0)]);
    let saved=destination.clone();
    assert!(mixed.add_to(&mut destination,||true).is_err());
    assert_eq!(destination,saved);assert_eq!(heat.cells(),&before);
}

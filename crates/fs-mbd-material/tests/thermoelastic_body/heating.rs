use super::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_matdb::SelectionPolicy;
use fs_mbd_material::thermal::heating::{HeatingError, ThermalRegion};
use fs_mbd_material::thermal::heating::constrained::{ConstrainedHeatedBody, ConstrainedHeatError, ConstrainedRegion};
use fs_thermal::enthalpy::*;
use fs_thermal::observers::ObserverCouplingError;
use std::collections::BTreeMap;

fn source(cp: f64, alpha_upper: f64, which: usize) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, dims, value, upper) in [
        (DENSITY_PROPERTY, Density::DIMS, if which == 0 {6.} else {12.}, 400.),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, cp, 400.),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, if which == 0 {1000.} else {2000.}, 400.),
        (POISSON_RATIO_PROPERTY, Dims([0;6]), 0.25, 400.),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, 1e-4, alpha_upper),
    ] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(key,dims),
            value: PropertyValue::Curve { abscissa:"T".into(), abscissa_dims:Temperature::DIMS,
                knots:vec![(300.,value),(upper,value)], dims },
            validity:ValidityDomain::unconstrained().with("T",300.,upper).with("pressure",100000.,100000.),
            uncertainty:UncertaintyModel::Unstated, interpolation:InterpolationPolicy::LinearInside,
            observations:vec![], provenance:Provenance {source:format!("synthetic heated tet {which}/{key}; not measurements"),
                license:"CC0-1.0".into(), artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:KEYS[which].into(),phase:"solid".into(),
        process:"synthetic heat fixture".into(),revision:0},claims,vec![]).unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    let point=|t|QueryPoint::new().with("T",t).unwrap().with("pressure",100000.).unwrap();
    ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),temperature_axis:"T",
        lower:&point(300.),upper:&point(400.),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.,1e-8).unwrap()
}
fn boundary(w: f64) -> Boundary { Boundary {ambient_k:300.,htc_w_m2_k:0.,reservoir:None,source_w:w} }
fn step() -> StepConfig { StepConfig {max_dt_s:0.25,energy_tolerance_j:1e-8} }
fn make<'a>(reference:&'a ResolvedBody,cards:&'a [MaterialCard],curves:&'a [ResolvedHeatCapacityCurve],
    watts:f64,contact:bool,mut events:Vec<TimedEvent>,supports:Vec<usize>,cx:&Cx<'_>,
) -> ConstrainedHeatedBody<'a> {
    let mesh=TetrahedralMass::bind(reference,&NODES,&TETS,&KEYS,5,||false).unwrap();
    let regions=(0..2).map(|i|ConstrainedRegion {storage:ThermalRegion {
        key:KEYS[i],card:&cards[i],caloric:&curves[i],cell_id:CellId(i as u64),
        material_id:MaterialId(i as u64),segment_id:i as u64,body_id:1,
        geometry:CellGeometry {start_m:[0.,0.25,if i==0 {0.25}else{-0.25}],
            end_m:[0.5,0.25,if i==0 {0.25}else{-0.25}],volume_m3:1./6.,
            surface_area_m2:1.5+3.0_f64.sqrt()/2.},
        initial_specific_enthalpy_j_kg:0.,boundary:boundary(if i==0 {watts}else{0.}),
        expansion_selection:MaterialPropertySelection::SingleClaimOnly,fixed_point_m:Vec3::default(),
    },stiffness_selection:MaterialPropertySelection::SingleClaimOnly}).collect();
    if contact {events.insert(0,TimedEvent {time_s:0.,event:Event::Contact(ThermalLink {
        id:10,a:CellId(0),b:CellId(1),area_m2:0.5,conductance_w_k:100.,
    })});}
    ConstrainedHeatedBody::new(mesh,regions,supports,config(),events,
        DepositionBudget {max_cells:2,max_links:1,max_events:32,max_substeps:128},step(),cx).unwrap()
}

#[test]
fn actual_heat_substeps_produce_the_same_constrained_field_as_a_fresh_fem_solve() {
    let cards=[source(10.,400.,0),source(10.,400.,1)];let reference=reference(&cards);
    let curves=cards.each_ref().map(caloric);
    with_cx(|cx| {
        let body=make(&reference,&cards,&curves,20.,true,vec![],SUPPORT.to_vec(),cx);
        let (plain,expected)=body.thermal().advance_to(1.,step(),||Ok(())).unwrap();
        let (next,report)=body.advance_to(1.,step(),128,cx).unwrap();
        assert_eq!(next.thermal().cells(),plain.cells());assert_eq!(report.thermal,expected);
        assert_eq!(report.shape_solves,expected.substeps);assert!(expected.substeps>1);
        let t=[plain.temperature_k(CellId(0)).unwrap(),plain.temperature_k(CellId(1)).unwrap()];
        assert!(t[0]>t[1] && t[1]>300.);
        let direct=body.mesh().solve_thermal(&inputs(&cards,&t),&SUPPORT,config(),cx).unwrap();
        assert_eq!(&direct,next.configuration());
        assert!(next.configuration().stress().elastic_energy_j()>0.);
        for (old,new) in reference.regions().iter().zip(next.configuration().body().regions()) {
            assert_eq!(old.inertia().mass().to_bits(),new.inertia().mass().to_bits());
        }
        near(next.thermal().energy_j(),20.);near(body.thermal().energy_j(),0.);
    });
}

#[test]
fn heat_capacity_only_changes_compatible_strain_stress_and_actual_rotor_response() {
    let mut responses=Vec::new();
    for cp in [10.,20.] {
        let cards=[source(cp,400.,0),source(cp,400.,1)];let reference=reference(&cards);
        let curves=cards.each_ref().map(caloric);
        with_cx(|cx| {
            let body=make(&reference,&cards,&curves,10.,false,vec![],SUPPORT.to_vec(),cx);
            let (next,report)=body.advance_to(1.,step(),128,cx).unwrap();
            near(next.thermal().temperature_k(CellId(0)).unwrap(),300.+10./cp);
            near(next.thermal().temperature_k(CellId(1)).unwrap(),300.);
            let shared=1e-4*(10./cp)/3.;
            near(next.configuration().displacement().displacement_m()[1][0],shared);
            near(next.configuration().stress().elements()[1].stress_mandel_pa()[0],2000./0.75*shared);
            near(report.thermal.external_heat_j,10.);
            responses.push(acceleration(next.configuration().body().inertia()));
        });
    }
    assert!(responses[0]<responses[1]);
}

#[test]
fn fully_restrained_heating_updates_stress_without_free_scaling() {
    let cards=[source(10.,400.,0),source(10.,400.,1)];let reference=reference(&cards);
    let curves=cards.each_ref().map(caloric);
    with_cx(|cx| {
        let body=make(&reference,&cards,&curves,10.,false,vec![],(0..15).collect(),cx);
        let (next,report)=body.advance_to(1.,step(),128,cx).unwrap();
        assert_eq!(next.configuration().body().inertia(),reference.inertia());
        assert_eq!(next.configuration().body().nodes_m(),NODES);
        near(next.configuration().stress().elements()[0].stress_mandel_pa()[0],-0.2);
        near(report.thermal.external_heat_j,10.);
        assert!(next.configuration().stress().elastic_energy_j()>0.);
    });
}

#[test]
fn a_bad_intermediate_source_state_cannot_hide_behind_later_cooling() {
    let cards=[source(10.,300.5,0),source(10.,400.,1)];let reference=reference(&cards);
    let curves=cards.each_ref().map(caloric);
    with_cx(|cx| {
        let body=make(&reference,&cards,&curves,10.,false,
            vec![TimedEvent {time_s:1.,event:Event::Boundary(CellId(0),boundary(-10.))}],SUPPORT.to_vec(),cx);
        let before=body.configuration().clone();
        let (plain,_)=body.thermal().advance_to(2.,step(),||Ok(())).unwrap();
        near(plain.temperature_k(CellId(0)).unwrap(),300.);
        assert!(matches!(body.advance_to(2.,step(),128,cx),Err(ConstrainedHeatError::Shape(_))));
        assert_eq!(body.configuration(),&before);near(body.thermal().time_s(),0.);
    });
}

#[derive(Clone,Debug,Default,PartialEq)]
struct Heater { joules:f64,calls:usize,fail:bool }
impl DepositionObserver for Heater {
    type Error=DepositionError;
    fn sources(&mut self,_:&DepositionSimulation<'_>,check:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error> {check()?;Ok(BTreeMap::from([(CellId(0),10.)]))}
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self,_:f64,dt:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        check()?;self.calls+=1;self.joules+=10.*dt;
        if self.fail && self.calls==3 {return Err(DepositionError::InvalidInput("late heat-owner refusal"));} Ok(())
    }
}

#[test]
fn generated_heat_and_its_history_commit_or_roll_back_with_actual_fem() {
    let cards=[source(10.,400.,0),source(10.,400.,1)];let reference=reference(&cards);
    let curves=cards.each_ref().map(caloric);
    with_cx(|cx| {
        let body=make(&reference,&cards,&curves,0.,true,vec![],SUPPORT.to_vec(),cx);
        let source=Heater::default();
        let (next,accepted,report)=body.advance_with(&source,1.,step(),128,cx).unwrap();
        near(accepted.joules,10.);near(report.thermal.external_heat_j,10.);near(next.thermal().energy_j(),10.);
        assert_eq!(accepted.calls,report.shape_solves);assert_eq!(accepted.calls,report.thermal.substeps);
        assert!(next.configuration().stress().elastic_energy_j()>0.);
        assert!(next.thermal().cells().values().all(|c|c.boundary.source_w==0.));
        let fail=Heater {fail:true,..Heater::default()};
        assert!(matches!(body.advance_with(&fail,1.,step(),128,cx),Err(ObserverCouplingError::First(_))));
        assert_eq!(source,Heater::default());assert_eq!(fail.calls,0);near(body.thermal().energy_j(),0.);
    });
}

#[test]
fn aligned_restart_no_op_and_work_refusal_preserve_complete_configuration() {
    let cards=[source(10.,400.,0),source(10.,400.,1)];let reference=reference(&cards);
    let curves=cards.each_ref().map(caloric);
    with_cx(|cx| {
        let body=make(&reference,&cards,&curves,10.,false,
            vec![TimedEvent {time_s:1.,event:Event::Boundary(CellId(0),boundary(0.))}],SUPPORT.to_vec(),cx);
        let (same,zero)=body.advance_to(0.,step(),0,cx).unwrap();
        assert_eq!(same.configuration(),body.configuration());assert_eq!(zero.shape_solves,0);
        let (whole,a)=body.advance_to(2.,step(),128,cx).unwrap();
        let (mid,b)=body.advance_to(1.,step(),128,cx).unwrap();
        let (end,c)=mid.advance_to(2.,step(),128,cx).unwrap();
        assert_eq!(whole.thermal().cells(),end.thermal().cells());assert_eq!(whole.configuration(),end.configuration());
        near(a.thermal.external_heat_j,b.thermal.external_heat_j+c.thermal.external_heat_j);
        assert_eq!(a.shape_solves,b.shape_solves+c.shape_solves);
        assert!(matches!(body.advance_to(2.,step(),1,cx),
            Err(ConstrainedHeatError::Heating(HeatingError::Thermal(DepositionError::Budget(_))))));
        near(body.thermal().time_s(),0.);
    });
}

//! Manufactured source data; independent rigid momentum and two-temperature
//! storage/entropy references, with optional retained Maxwell memory.
use super::*;
use fs_couple::friction_heat::third_body::rigid::thermal::*;
use fs_evidence::ValidityDomain;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, ResolvedShearLayer,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    SHEAR_LAYER_MODULUS_PROPERTY, SHEAR_LAYER_YIELD_PROPERTY, PLASTIC_SHEAR_VISCOSITY_PROPERTY};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    InterpolationPolicy, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{Density, Dims, DynViscosity, Pressure, Temperature, Time};

const CP: &str = "specific_heat_capacity";
const RHO: &str = "density";
const CP_DIMS: Dims = Dims([2,0,-2,-1,0,0]);
const ARM: &str = "manufactured_memory_shear_modulus";
const TAU: &str = "manufactured_memory_time";

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    source: ResolvedShearLayer,
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T",temperature).unwrap().with("omega",1000.0).unwrap()
}
fn fixture(slope: f64, arm_slope: Option<f64>, yield_pa: f64, kink: bool) -> Fixture {
    let mut claims=ClaimSet::new();
    let mut requirements=Vec::new();
    let mut roles=vec![
        (SHEAR_LAYER_MODULUS_PROPERTY,Pressure::DIMS,10.0,Some(slope),true),
        (SHEAR_LAYER_YIELD_PROPERTY,Pressure::DIMS,yield_pa,None,true),
        (PLASTIC_SHEAR_VISCOSITY_PROPERTY,DynViscosity::DIMS,0.3,None,true),
        (CP,CP_DIMS,1.0,None,false),(RHO,Density::DIMS,CAP[2],None,false),
    ];
    if let Some(slope)=arm_slope {
        roles.push((ARM,Pressure::DIMS,3.0,Some(slope),true));
        roles.push((TAU,Time::DIMS,0.4,None,true));
    }
    for (name,dims,value,slope,mechanical) in roles {
        let data=if let Some(slope)=slope {
            PropertyValue::Curve { abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,value-20.0*slope),
                    (300.0,value+if kink && name==SHEAR_LAYER_MODULUS_PROPERTY {0.1}else{0.0}),
                    (320.0,value+20.0*slope)],dims }
        } else { PropertyValue::Scalar {value,dims} };
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),value:data,
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
            interpolation:if slope.is_some(){InterpolationPolicy::LinearInside}else{InterpolationPolicy::ConstantWithinValidity},
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"manufactured rigid thermoelastic layer".into(),license:"test-only".into(),artifact:None}
        }).unwrap();
        if mechanical {
            requirements.push(ScalarPropertyRequirement::try_new(name,dims,
                if name==SHEAR_LAYER_YIELD_PROPERTY {ScalarAdmissibility::NonNegative}
                else {ScalarAdmissibility::StrictlyPositive}).unwrap());
        }
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"manufactured".into(),
        phase:"solid".into(),process:"rigid layer test".into(),revision:0},claims,vec![]).unwrap();
    let source=resolve_material_state_point(&card,&point(300.0),&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let mut source=ResolvedShearLayer::new(source,1.0,1.0,0.1).unwrap();
    if arm_slope.is_some(){source=source.with_relaxation(vec![(ARM.into(),TAU.into())]).unwrap();}
    let caloric=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new(CP,CP_DIMS),density:&PropertyKey::new(RHO,Density::DIMS),
        temperature_axis:"T",lower:&point(280.0),upper:&point(320.0),
        selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,-20.0,1e-8).unwrap();
    Fixture {card,caloric,source}
}
fn bind(f: &Fixture) -> ThermoelasticLayerMaterial<'_> {
    ThermoelasticLayerMaterial::try_new(f.source.clone(),&f.card,&f.caloric,"T",1000.0).unwrap()
}
fn inventory(f: &Fixture) -> DepositionSimulation<'_> {
    let original=thermal();
    let mut events=Vec::new();
    for cell in original.cells().values() {
        let mut cell=cell.clone();
        if cell.id==LAYER {cell.material_id=MaterialId(1);}
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(cell)});
    }
    for link in original.links().values() {
        events.push(TimedEvent {time_s:0.0,event:Event::Contact(*link)});
    }
    DepositionSimulation::new(events,&[(MaterialId(0),&CALORIC as &dyn ThermalLaw),
        (MaterialId(1),f.caloric.curve() as &dyn ThermalLaw)],DepositionBudget {
            max_cells:3,max_links:2,max_events:5,max_substeps:32,
        }).unwrap().advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-10},||Ok(())).unwrap().0
}
fn transport(source: &ResolvedShearLayer, power_w: f64) -> RigidLayerThermalTransport<ResolvedShearLayer> {
    RigidLayerThermalTransport {conductivity_w_m_k:K,area_resistance_m2_k_w:[0.0;2],
        transport:StepCoefficients {sources_w:if power_w==0.0 {BTreeMap::new()}
            else {BTreeMap::from([(LAYER,power_w)])},conductances_w_k:BTreeMap::new(),payload:source.clone()}}
}

// Closed-form collinear force at a GIVEN temperature. At zero initial elastic
// and memory strain, the Maxwell midpoint gives k_eff=k+ka/(1+dt/(2tau)).
// The two physical face wrenches give W=17/6, including their spin mobility.
fn reference(theta: f64, slope: f64, arm_slope: Option<f64>, yield_n: f64, dt: f64)
    -> (f64,f64,f64,f64,f64) {
    let k=10.0+slope*theta;
    let arm=arm_slope.map_or(0.0,|s|3.0+s*theta);
    let factor=1.0+dt/(2.0*0.4);
    let effective=k+arm/factor;
    let a=1.0+effective*dt*dt*(17.0/6.0)/4.0;
    let trial=effective*dt*0.8/(2.0*a);
    let b=effective*dt/(2.0*0.3);
    let force=if trial<=yield_n {trial}else{(effective*dt*0.8/2.0+b*yield_n)/(a+b)};
    let elastic=2.0*force/effective;
    let memory=if arm_slope.is_some(){elastic/factor}else{0.0};
    let plastic=dt/0.3*(force-yield_n).max(0.0);
    let dissipated=force*plastic+dt*arm*(0.5*memory).powi(2)/0.4;
    let reversible=0.5*(300.0+theta)*(slope*elastic*elastic+arm_slope.unwrap_or(0.0)*memory*memory);
    (force,dissipated,reversible,elastic,memory)
}
fn endpoint(slope: f64, arm_slope: Option<f64>, yield_n: f64, power_w: f64, dt: f64) -> f64 {
    let g=[0.8,2.0/3.0];
    let capacity=CAP[2]+dt*g[0]*CAP[0]/(CAP[0]+dt*g[0])+dt*g[1]*CAP[1]/(CAP[1]+dt*g[1]);
    let (mut lower,mut upper)=(-2.0,5.0);
    for _ in 0..100 {
        let theta=0.5*(lower+upper);
        let (_,d,q,_,_)=reference(theta,slope,arm_slope,yield_n,dt);
        if capacity*theta-d-q-dt*power_w>0.0 {upper=theta;}else{lower=theta;}
    }
    0.5*(lower+upper)
}

#[test]
fn g1_affine_material_and_memory_softening_close_actual_rigid_momentum_heat_and_entropy() {
    let dt=0.1;
    let mut forces=Vec::new();
    for (slope,arm_slope,yield_n,power) in [(-0.02,None,0.2,4.0),(0.02,None,0.2,4.0),
        (-0.02,Some(-0.01),0.2,4.0),(-0.02,None,100.0,0.0)] {
        let f=fixture(slope,arm_slope,yield_n,false);
        let material=bind(&f);
        let thermal=inventory(&f);
        let route=ThirdBodyHeatRoute::from_links(&thermal,LAYER,[7,8],K,[0.0;2]).unwrap();
        let (_,law)=material.at_temperature(300.0).unwrap();
        let history=RigidLayerHistory::new(&law,geometry(),0.0).unwrap();
        let old=bodies(false,Vec3::new(0.0,0.8,0.0));
        let candidate=step_thermoelastic_heated_rigid_layer_candidate(&old,properties(false),&history,
            &thermal,&route,&material,dt,10.0,law_solve(),heat_controls(),[0.0;2],&[0.0;3],
            |_,source|Ok::<_,&str>(transport(source,power)),||false).unwrap();
        let theta=endpoint(slope,arm_slope,yield_n,power,dt);
        let (force,dissipated,reversible,elastic,memory)=reference(theta,slope,arm_slope,yield_n,dt);
        let report=candidate.report();
        close(report.heat.layer_temperature_k,300.0+theta,3e-8);
        close(report.mechanical.applied_force_on_a_n.y,force,3e-9);
        close(candidate.history().deformation().elastic_slip_m()[0],elastic,3e-9);
        if arm_slope.is_some(){close(candidate.history().deformation().relaxation_slip_m()[0][0],memory,3e-9);}
        close(report.heat.dissipated_j,dissipated,4e-10);
        close(report.heat.reversible_heat_j,reversible,4e-10);
        close(report.heat.external_heat_j,dt*power,1e-12);
        close(report.combined_energy_residual_j,0.0,5e-10);
        let account=report.thermoelastic.unwrap();
        let internal=0.5*(10.0-300.0*slope)*elastic*elastic
            +arm_slope.map_or(0.0,|s|0.5*(3.0-300.0*s)*memory*memory);
        close(account.internal_energy_change_j,internal,4e-10);
        close(account.reversible_heat_j,reversible,4e-10);
        close(account.coefficient_free_energy_change_j,0.0,0.0);
        let entropy=report.heat.entropy.unwrap();
        close(entropy.friction_generation_j_k,dissipated/(300.0+theta),2e-12);
        close(entropy.mechanical_storage_change_j_k,-reversible/(300.0+theta),2e-12);
        close(entropy.reversible_caloric_exchange_j_k,reversible/(300.0+theta),2e-12);
        close(entropy.thermal_storage_change_j_k+entropy.mechanical_storage_change_j_k,
            entropy.reservoir_entropy_input_j_k+entropy.network_conduction_generation_j_k
                +entropy.interface_conduction_generation_j_k+entropy.friction_generation_j_k
                +entropy.unclassified_source_at_cell_j_k+entropy.temporal_defect_j_k
                +entropy.interface_residual_entropy_j_k,2e-12);
        assert_eq!(entropy.has_unclassified_sources,power!=0.0);
        for side in 0..2 {
            let sign=if side==0{1.0}else{-1.0};
            vector_close(candidate.bodies()[side].linear_momentum_world().sub(old[side].linear_momentum_world()),
                Vec3::new(0.0,sign*dt*force,0.0),3e-10);
            vector_close(candidate.bodies()[side].angular_momentum_body(),Vec3::new(0.0,0.0,dt*force),3e-10);
        }
        vector_close(report.mechanical.linear_momentum_change_ns,Vec3::ZERO,3e-10);
        vector_close(report.mechanical.angular_momentum_change_nms,Vec3::ZERO,3e-10);
        close(report.source.stiffness_n_m(),10.0+slope*theta,2e-9);
        assert_eq!(report.source.resolved().card_identity(),f.card.content_hash());
        if power==0.0 {assert_eq!(dissipated,0.0);assert!(reversible<0.0);assert!(theta<0.0);}
        forces.push(force);
        assert_eq!(history.time_s(),0.0);assert_eq!(thermal.time_s(),0.0);
    }
    assert!((forces[0]-forces[1]).abs()>1e-5,"source modulus must change solved face forces");
}

#[test]
fn g1_warm_continuation_retains_original_memory_and_actual_two_temperature_storage() {
    let f=fixture(-0.02,Some(-0.01),0.2,false);
    let material=bind(&f);
    let mut thermal=inventory(&f);
    let route=ThirdBodyHeatRoute::from_links(&thermal,LAYER,[7,8],K,[0.0;2]).unwrap();
    let (_,law)=material.at_temperature(300.0).unwrap();
    let mut history=RigidLayerHistory::new(&law,geometry(),0.0).unwrap();
    let mut bodies=bodies(false,Vec3::new(0.0,0.8,0.0));
    for step in 1..=3 {
        let slopes=material.thermal_slopes();
        let old=history.law().thermodynamics_with_relaxation(history.deformation(),
            slopes.equilibrium_n_m_k,&slopes.relaxation_n_m_k).unwrap();
        let candidate=step_thermoelastic_heated_rigid_layer_candidate(&bodies,properties(false),&history,
            &thermal,&route,&material,step as f64*0.05,10.0,law_solve(),heat_controls(),[0.0;2],&[0.0;3],
            |_,source|Ok::<_,&str>(transport(source,4.0)),||false).unwrap();
        let new=candidate.history().law().thermodynamics_with_relaxation(candidate.history().deformation(),
            slopes.equilibrium_n_m_k,&slopes.relaxation_n_m_k).unwrap();
        let account=candidate.report().thermoelastic.unwrap();
        close(account.internal_energy_change_j,new.internal_energy_j-old.internal_energy_j,2e-11);
        close(account.free_energy_change_j,new.free_energy_j-old.free_energy_j,2e-11);
        close(account.mechanical_entropy_change_j_k,new.entropy_j_k-old.entropy_j_k,2e-12);
        close(account.reversible_heat_j+candidate.report().heat.layer_temperature_k
            *account.mechanical_entropy_change_j_k,0.0,2e-10);
        if step>1 {assert_ne!(account.coefficient_free_energy_change_j,0.0);}
        assert_eq!(candidate.history().deformation().relaxation_slip_m().len(),1);
        close(candidate.report().combined_energy_residual_j,0.0,5e-10);
        let (next_bodies,next_history,next_thermal,_)=candidate.into_parts();
        bodies=next_bodies;history=next_history;thermal=next_thermal;
    }
}

#[test]
fn g4_nonaffine_foreign_and_late_source_failures_publish_no_partial_history() {
    let kink=fixture(-0.02,None,0.2,true);
    assert!(ThermoelasticLayerMaterial::try_new(kink.source.clone(),&kink.card,&kink.caloric,"T",1000.0).is_err());
    let f=fixture(-0.02,Some(-0.01),0.2,false);
    let foreign=fixture(0.02,Some(-0.01),0.2,false);
    assert!(ThermoelasticLayerMaterial::try_new(f.source.clone(),&f.card,&foreign.caloric,"T",1000.0).is_err());
    let material=bind(&f);
    assert!(material.at_temperature(321.0).is_err());
    let thermal=inventory(&f);
    let route=ThirdBodyHeatRoute::from_links(&thermal,LAYER,[7,8],K,[0.0;2]).unwrap();
    let (_,law)=material.at_temperature(300.0).unwrap();
    let history=RigidLayerHistory::new(&law,geometry(),0.0).unwrap();
    let old=bodies(false,Vec3::new(0.0,0.8,0.0));
    let saved=history.deformation().clone();
    let calls=Cell::new(0);
    let refused=step_thermoelastic_heated_rigid_layer_candidate(&old,properties(false),&history,
        &thermal,&route,&material,0.1,10.0,law_solve(),heat_controls(),[0.0;2],&[0.0;3],
        |_,source|{calls.set(calls.get()+1);if calls.get()==3 {Err("late transport source")}
            else {Ok(transport(source,4.0))}},||false);
    assert!(matches!(refused,Err(HeatedRigidLayerError::Source(ThermoelasticLayerError::Transport("late transport source")))));
    assert_eq!(history.deformation(),&saved);assert_eq!(history.time_s(),0.0);assert_eq!(thermal.time_s(),0.0);
    let mismatched=super::thermal();
    let wrong_route=ThirdBodyHeatRoute::from_links(&mismatched,LAYER,[7,8],K,[0.0;2]).unwrap();
    assert!(step_thermoelastic_heated_rigid_layer_candidate(&old,properties(false),&history,
        &mismatched,&wrong_route,&material,0.1,10.0,law_solve(),heat_controls(),[0.0;2],&[0.0;3],
        |_,source|Ok::<_,&str>(transport(source,4.0)),||false).is_err());
}

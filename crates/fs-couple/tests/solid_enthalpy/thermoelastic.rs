//! G1/G3/G4: finite heat and the SAME native load-bearing elastic trajectory.
//! Synthetic source values; these tests do not qualify a measured material.
use super::{with_cx,near,source,NODES,TETS,FREE};
use fs_blake3::ContentHash;
use fs_couple::solid_enthalpy::{SolidEnthalpy,SolidEquilibriumConfig};
use fs_couple::solid_enthalpy::thermoelastic::*;
use fs_evidence::ValidityDomain;
use fs_material::phase::{EnthalpyPhaseKnot,EquilibriumEnthalpyPhaseCurve,SolidLiquidPhase};
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS,SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_qty::{Density,Dims,Pressure,QuantitySpec,Temperature,semantic::{QuantityKind,SemanticType,ValueForm}};
use fs_solid::linear3::{TetAssemblyBudget,TetStaticSolveConfig};
use fs_solid::loads3::VolumeLoad3;
use fs_thermal::enthalpy::*;
use std::collections::BTreeMap;

const IDS:[CellId;2]=[CellId(30),CellId(10)];
const VOLUME:f64=1.0/6.0;
const PAIR_NODES:[[f64;3];5]=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,-1.,0.]];
const PAIR_TETS:[[usize;4];2]=[[0,1,2,3],[0,1,3,4]];
fn quantity()->QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature,ValueForm::Static))
}
fn context(t:f64)->QueryPoint {QueryPoint::new().with_quantity("T",quantity(),t).unwrap()}
fn fixture(young:f64,slope:f64,alpha:f64,cp:f64,full_phase:bool)
    ->(MaterialCard,EquilibriumEnthalpyPhaseCurve) {
    let mut claims=ClaimSet::new();
    for (name,dims,value,derivative) in [
        ("young_modulus",Pressure::DIMS,young,slope),
        ("poisson_ratio",Dims([0;6]),0.0,0.0),
        ("linear_thermal_expansion_coefficient",INVERSE_TEMPERATURE_DIMS,alpha,0.0),
        ("density",Density::DIMS,6.0,0.0),
        ("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS,cp,0.0),
    ] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,value-20.0*derivative),(500.0,value+200.0*derivative)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",quantity(),280.0,500.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"authored spatial thermoelastic fixture".into(),
                license:"test-only".into(),artifact:None}}).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic elastic solid".into(),
        phase:"solid".into(),process:"numerical source".into(),revision:0},claims,vec![]).unwrap();
    let knot=|t,h,f,rho|EnthalpyPhaseKnot {temperature_k:t,specific_enthalpy_j_kg:h,
        liquid_mass_fraction:f,bulk_density_kg_m3:rho};
    let curve=if full_phase {
        EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),vec![
            knot(280.0,0.0,0.0,6.0),knot(400.0,120.0*cp,0.0,6.0),
            knot(400.0,120.0*cp+10.0,1.0,3.0),knot(500.0,220.0*cp+10.0,1.0,3.0),
        ]).unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(),SolidLiquidPhase::Solid,
            vec![knot(280.0,0.0,0.0,6.0),knot(500.0,220.0*cp,0.0,6.0)]).unwrap()
    };
    (card,curve)
}
fn mechanical()->SolidEquilibriumConfig {
    SolidEquilibriumConfig {static_solve:TetStaticSolveConfig::default(),
        mesh_frame:ContentHash([1;32]),maximum_nodal_displacement_m:0.1}
}
fn solve()->ThermoelasticSolve {
    ThermoelasticSolve {maximum_iterations:64,heat_relaxation:1.0,heat_tolerance_j:1e-11,
        energy_tolerance_j:1e-9,stress_tolerance_pa:1e-8}
}
fn thermal_controls()->StepConfig {StepConfig {max_dt_s:0.025,energy_tolerance_j:1e-10}}
fn thermal<'a>(curves:&[&'a EquilibriumEnthalpyPhaseCurve],cp:&[f64],power:&[f64],conductance:f64)
    ->DepositionSimulation<'a> {
    let mut events:Vec<_>=curves.iter().enumerate().map(|(i,_)|TimedEvent {time_s:0.0,
        event:Event::Birth(DepositedCell {id:IDS[i],segment_id:i as u64,body_id:3,
            material_id:MaterialId(i as u64),mass_kg:1.0,specific_h_j_kg:20.0*cp[i],
            geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:VOLUME,surface_area_m2:3.0},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:power[i]}})}).collect();
    if curves.len()==2 {events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:9,a:IDS[0],b:IDS[1],area_m2:0.5,conductance_w_k:conductance})});}
    let materials:Vec<_>=curves.iter().enumerate().map(|(i,curve)|
        (MaterialId(i as u64),*curve as &dyn ThermalLaw)).collect();
    DepositionSimulation::new(events,&materials,DepositionBudget {
        max_cells:2,max_links:1,max_events:4,max_substeps:1000,
    }).unwrap().advance_to(0.0,thermal_controls(),||Ok(())).unwrap().0
}
fn root(mut low:f64,mut high:f64,residual:impl Fn(f64)->f64)->f64 {
    assert!(residual(low)<0.0 && residual(high)>0.0);
    for _ in 0..100 {let mid=0.5*(low+high);
        if residual(mid)<0.0 {low=mid;} else {high=mid;}
    }
    0.5*(low+high)
}

// Independent nu=0 invariants with one common axial DOF and two clamped
// transverse directions. No native material/source/mesh helper enters this oracle.
fn element(t:f64,u:f64,young:f64,slope:f64,alpha:f64)->(f64,f64) {
    let e=young+slope*(t-300.0); let free=alpha*(t-300.0);
    let stress=e*(u-free);
    let psi=0.5*VOLUME*e*((u-free).powi(2)+2.0*free*free);
    let entropy=VOLUME*alpha*e*(u-3.0*free)-(slope/e)*psi;
    (stress,psi+t*entropy)
}
fn displacement(t:[f64;2],force:f64,young:[f64;2],slope:[f64;2],alpha:[f64;2])->f64 {
    let k=[0,1].map(|i|VOLUME*(young[i]+slope[i]*(t[i]-300.0)));
    (force+(0..2).map(|i|k[i]*alpha[i]*(t[i]-300.0)).sum::<f64>())/(k[0]+k[1])
}
#[allow(clippy::too_many_arguments)]
fn pair_endpoint(old_t:[f64;2],old_u:f64,force:f64,young:[f64;2],slope:[f64;2],
    alpha:[f64;2],capacity:[f64;2],power:[f64;2],conductance:f64,dt:f64)->([f64;2],f64,[f64;2]) {
    let evaluate=|t:[f64;2]| {
        let u=displacement(t,force,young,slope,alpha);
        let heat=[0,1].map(|i| {
            let (s0,u0)=element(old_t[i],old_u,young[i],slope[i],alpha[i]);
            let (s1,u1)=element(t[i],u,young[i],slope[i],alpha[i]);
            0.5*VOLUME*(s0+s1)*(u-old_u)-(u1-u0)
        });
        // The declared four forward-Euler thermal substeps use the same
        // conservative interface transfer; heat is a constant interval power.
        let mut predicted=old_t;
        for _ in 0..4 {let before=predicted;
            predicted=[0,1].map(|i|before[i]+dt/4.0/capacity[i]
                *(power[i]+heat[i]/dt+conductance*(before[1-i]-before[i])));
        }
        ([0,1].map(|i|capacity[i]*(t[i]-predicted[i])),u,heat)
    };
    // Independent two-variable Newton elimination, distinct from the runtime's
    // heat fixed point. This is a small analytic-fixture reference only.
    let mut t=old_t;
    for _ in 0..30 {
        let (r,u,q)=evaluate(t);
        if r[0].abs().max(r[1].abs())<1e-12 {return(t,u,q);}
        let mut j=[[0.0;2];2];
        for col in 0..2 {let mut probe=t;probe[col]+=1e-4;let next=evaluate(probe).0;
            for row in 0..2 {j[row][col]=(next[row]-r[row])/1e-4;}
        }
        let determinant=j[0][0]*j[1][1]-j[0][1]*j[1][0];
        t[0]-=(j[1][1]*r[0]-j[0][1]*r[1])/determinant;
        t[1]-=(j[0][0]*r[1]-j[1][0]*r[0])/determinant;
        assert!(t.iter().all(|v|(280.0..360.0).contains(v)));
    }
    panic!("independent two-cell thermoelastic reference did not converge");
}

#[test]
fn g1_free_and_constrained_heating_use_the_same_finite_thermoelastic_potential() {
    with_cx(|cx| {
        for (alpha,clamped) in [(0.0003,false),(0.0006,false),(0.0003,true)] {
            let (card,curve)=fixture(1000.0,2.0,alpha,0.2,false);
            let thermal=thermal(&[&curve],&[0.2],&[2.0],0.0);let context=context(300.0);
            let source=SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
                vec![source(30,&card,&curve,&context)],TetAssemblyBudget::standard(),
                0.05,mechanical().mesh_frame,cx).unwrap();
            let loads=source.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
            let fixed=if clamped {(0..12).collect::<Vec<_>>()} else {FREE.to_vec()};
            let body=ThermoelasticSolid::try_new(&source,&thermal,&loads,&fixed,&[[280.0,360.0]],
                mechanical(),solve(),cx).unwrap();
            let next=body.advance(&thermal,&(),0.1,&loads,thermal_controls(),solve(),cx).unwrap();
            let expected_t=if clamped {
                root(280.0,360.0,|t|0.2*(t-300.0)+element(t,0.0,1000.0,2.0,alpha).1-0.2)
            } else {301.0};
            near(next.thermal.temperature_k(IDS[0]).unwrap(),expected_t,2e-8);
            let u=next.solid.equilibrium().solution.displacement.displacement_m();
            let free=alpha*(expected_t-300.0);
            for (node,point) in NODES.iter().enumerate() {for k in 0..3 {
                near(u[node][k],if clamped {0.0} else {free*point[k]},2e-10);
            }}
            let stress=next.solid.equilibrium().solution.stress.elements()[0].stress_mandel_pa();
            for k in 0..3 {near(stress[k],if clamped {
                -(1000.0+2.0*(expected_t-300.0))*free
            } else {0.0},2e-8);}
            near(next.report.mechanical_work_j,0.0,1e-12);
            near(next.report.external_heat_j,0.2,1e-9);
            near(next.report.combined_residual_j,0.0,1e-9);
            near(next.thermal.energy_j()-thermal.energy_j()+next.report.internal_energy_change_j,0.2,1e-9);
            if clamped {assert!(expected_t>301.07);assert!(next.report.reversible_heat_j>0.01);}
            else {near(next.report.reversible_heat_j,0.0,1e-10);}
            assert_eq!(thermal.time_s(),0.0);assert_eq!(body.time_s(),0.0);
        }
    });
}

#[test]
fn g1_heterogeneous_temperature_loads_and_source_swaps_share_one_compatible_trajectory() {
    with_cx(|cx| {
        let mut outcomes=Vec::new();
        for scale in [1.0,2.0] {
            let young=[1000.0,1600.0];let slope=[2.0,-1.0];
            let alpha=[0.0002*scale,0.0003*scale];let capacity=[0.2,0.3];let power=[1.0,0.1];
            let first=fixture(young[0],slope[0],alpha[0],capacity[0],false);
            let second=fixture(young[1],slope[1],alpha[1],capacity[1],false);
            let cards=[first.0,second.0];let curves=[first.1,second.1];let context=context(300.0);
            let mut thermal=thermal(&[&curves[0],&curves[1]],&capacity,&power,0.02);
            let source=SolidEnthalpy::try_new(&thermal,&PAIR_NODES,&PAIR_TETS,
                vec![source(30,&cards[0],&curves[0],&context),source(10,&cards[1],&curves[1],&context)],
                TetAssemblyBudget::standard(),0.05,mechanical().mesh_frame,cx).unwrap();
            let fixed:Vec<_>=(0..15).filter(|d|*d!=3).collect();
            let initial_load=source.geometry().assemble(&[],VolumeLoad3::None,&[(3,0.8)],cx).unwrap();
            let mut body=ThermoelasticSolid::try_new(&source,&thermal,&initial_load,&fixed,
                &[[280.0,360.0];2],mechanical(),solve(),cx).unwrap();
            let mut old_t=[300.0;2];let mut old_u=displacement(old_t,0.8,young,slope,alpha);let mut old_force=0.8;
            let initial_energy=thermal.energy_j()+body.potentials().iter().map(|p|p.internal_energy_j).sum::<f64>();
            let mut work=0.0;let mut heat=0.0;
            for (interval,force) in [1.0,0.6,0.8].into_iter().enumerate() {
                let target=(interval+1) as f64*0.1;let dt=target-thermal.time_s();
                let expected=pair_endpoint(old_t,old_u,force,young,slope,alpha,capacity,power,0.02,dt);
                let loads=source.geometry().assemble(&[],VolumeLoad3::None,&[(3,force)],cx).unwrap();
                let next=body.advance(&thermal,&(),target,&loads,thermal_controls(),solve(),cx).unwrap();
                for i in 0..2 {
                    let t=next.thermal.temperature_k(IDS[i]).unwrap();
                    near(t,expected.0[i],3e-8);
                    near(next.report.cells[i].reversible_heat_j,expected.2[i],2e-9);
                    assert_eq!(next.report.cells[i].cell,IDS[i]);
                    assert_eq!(next.solid.state().samples()[i].elastic.card_identity(),cards[i].content_hash());
                    near(next.solid.state().samples()[i].reference_mass_kg,1.0,1e-14);
                    let sigma=element(t,expected.1,young[i],slope[i],alpha[i]).0;
                    near(next.solid.equilibrium().solution.stress.elements()[i].stress_mandel_pa()[0],sigma,2e-7);
                }
                near(next.solid.equilibrium().solution.displacement.displacement_m()[1][0],expected.1,2e-10);
                near(next.report.mechanical_work_j,0.5*(old_force+force)*(expected.1-old_u),2e-10);
                near(next.report.external_heat_j,dt*(power[0]+power[1]),2e-9);
                near(next.report.mechanical_residual_j,0.0,1e-9);
                near(next.report.combined_residual_j,0.0,1e-9);
                work+=next.report.mechanical_work_j;heat+=next.report.external_heat_j;
                near(next.thermal.energy_j()+next.solid.potentials().iter().map(|p|p.internal_energy_j).sum::<f64>(),
                    initial_energy+work+heat,4e-9);
                let (next_thermal,next_body,_,_)=next.into_parts();thermal=next_thermal;body=next_body;
                old_t=expected.0;old_u=expected.1;old_force=force;
            }
            assert!((old_t[0]-old_t[1]).abs()>1.0);
            outcomes.push((old_t,old_u));
        }
        assert!((outcomes[0].1-outcomes[1].1).abs()>1e-4);
        assert!((outcomes[0].0[0]-outcomes[1].0[0]).abs()>1e-3);
    });
}

#[derive(Clone,Debug,PartialEq)]
struct Steps {count:usize,cancel_at:usize}
impl DepositionObserver for Steps {
    type Error=DepositionError;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self,_:f64,_:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        _:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        self.count+=1;
        if self.count==self.cancel_at {Err(DepositionError::Cancelled)} else {Ok(())}
    }
}

#[test]
fn g4_complete_phase_chart_is_retained_and_failed_coupling_publishes_no_history() {
    with_cx(|cx| {
        let (card,curve)=fixture(1000.0,2.0,0.0003,0.2,true);
        let thermal=thermal(&[&curve],&[0.2],&[2.0],0.0);let context=context(300.0);
        let source=SolidEnthalpy::try_new(&thermal,&NODES,&TETS,
            vec![source(30,&card,&curve,&context)],TetAssemblyBudget::standard(),
            0.05,mechanical().mesh_frame,cx).unwrap();
        let fixed:Vec<_>=(0..12).collect();let loads=source.geometry().assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        // The SAME chart retains its later latent plateau and liquid branch.
        let body=ThermoelasticSolid::try_new(&source,&thermal,&loads,&fixed,&[[280.0,360.0]],
            mechanical(),solve(),cx).unwrap();
        assert!(ThermoelasticSolid::try_new(&source,&thermal,&loads,&fixed,&[[280.0,420.0]],
            mechanical(),solve(),cx).is_err());
        let before_cells=thermal.cells().clone();let before_u=body.equilibrium().solution.displacement.clone();
        let observer=Steps {count:0,cancel_at:usize::MAX};
        let accepted=body.advance(&thermal,&observer,0.1,&loads,thermal_controls(),solve(),cx).unwrap();
        assert!(accepted.observer.count>=4);assert_eq!(observer.count,0);
        let mut budget=solve();budget.maximum_iterations=1;
        assert!(matches!(body.advance(&thermal,&observer,0.1,&loads,thermal_controls(),budget,cx),
            Err(ThermoelasticFeedbackError::NotConverged {..})));
        let late_cancel=Steps {count:0,cancel_at:2};
        assert!(body.advance(&thermal,&late_cancel,0.1,&loads,thermal_controls(),solve(),cx).is_err());
        assert_eq!(late_cancel.count,0);
        let narrow=ThermoelasticSolid::try_new(&source,&thermal,&loads,&fixed,&[[280.0,300.01]],
            mechanical(),solve(),cx).unwrap();
        assert!(narrow.advance(&thermal,&observer,0.1,&loads,thermal_controls(),solve(),cx).is_err());
        let (hot,_) = thermal.advance_to(0.025,thermal_controls(),||Ok(())).unwrap();
        assert!(body.advance(&hot,&observer,0.1,&loads,thermal_controls(),solve(),cx).is_err());
        assert_eq!(thermal.cells(),&before_cells);assert_eq!(thermal.time_s(),0.0);
        assert_eq!(body.equilibrium().solution.displacement,before_u);assert_eq!(body.time_s(),0.0);
        let retry=body.advance(&thermal,&observer,0.1,&loads,thermal_controls(),solve(),cx).unwrap();
        assert_eq!(retry.thermal.cells(),accepted.thermal.cells());
        assert_eq!(retry.solid.equilibrium().solution.displacement,accepted.solid.equilibrium().solution.displacement);
        assert_eq!(retry.observer,accepted.observer);
    });
}

#[path = "thermoelastic/dynamic.rs"]
mod dynamic;

#[path = "thermoelastic/varying_coefficients.rs"]
mod varying_coefficients;

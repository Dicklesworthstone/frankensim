//! G1/G3/G4: original caloric heating, compatible cut-solid response and rollback.
use std::collections::BTreeMap;
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_couple::cut_thermoelastic::*;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticEquilibriumSettings3, ThermoelasticExchangeLimits3,
};
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityFormulation3, ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_ivl::Interval;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::elastic::thermal::{
    AffineThermoelastic, IsothermalThermoelasticInput, IsothermalThermoelasticState,
    ThermoelasticSymmetry,
};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve,
    SolidLiquidPhase,
};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solid::cutfront::dynamics3::{CutElasticityRef3, CutSecondOrderProblem3};
use fs_solver::LinearOp;
use fs_time::galpha::SecondOrderProblem;
use fs_thermal::enthalpy::*;


const CP: f64 = 20.0;
const DENSITY: f64 = 7.0;
const POINTS: usize = 100_000;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn proceed() -> ControlFlow<()> { ControlFlow::Continue(()) }
struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64;3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo:[f64;3], hi:[f64;3]) -> Interval {
        Interval::new(lo[2]-0.73, hi[2]-0.73)
    }
    fn derivative_enclose(&self, _: [f64;3], _: [f64;3], axis:HeightAxis) -> Interval {
        let x=if axis==HeightAxis::Z {1.0}else{0.0};
        Interval::new(x,x)
    }
}
enum Model { Cartesian(CutElasticity3), Adaptive(AdaptiveElasticity3) }
macro_rules! native {
    ($owner:expr, $method:ident $(, $arg:expr)* $(,)?) => {
        match $owner {
            Model::Cartesian(op) => op.$method($($arg),*),
            Model::Adaptive(op) => op.$method($($arg),*),
        }
    };
}
impl Model {
    fn source(&self)->CutElasticityRef3<'_> {
        match self {
            Self::Cartesian(op)=>CutElasticityRef3::Cartesian(op),
            Self::Adaptive(op)=>CutElasticityRef3::Adaptive(op),
        }
    }
    fn mutable(&mut self)->CutElasticityMut3<'_> {
        match self {
            Self::Cartesian(op)=>CutElasticityMut3::Cartesian(op),
            Self::Adaptive(op)=>CutElasticityMut3::Adaptive(op),
        }
    }
    fn operator(&self)->&dyn LinearOp {
        match self {Self::Cartesian(op)=>op,Self::Adaptive(op)=>op}
    }
}
fn build(adaptive:bool, source:&IsothermalThermoelasticState)->Model {
    let domain=HexCell::try_new([0.0;3],[1.0;3]).unwrap();
    let mut callback=|_|proceed();
    let mut control=QuadratureControl3::new(
        QuadratureOptions3 {depth:0,..Default::default()},&mut callback).unwrap();
    let options=ElasticityOptions3 {formulation:ElasticityFormulation3::BBar,..Default::default()};
    let field=ElasticMaterialField3::ThermoelasticUniform(source);
    let clamp=|p:[f64;3]|p[0]==0.0;
    if adaptive {
        let tree=Octree3::uniform(1,4,4096).unwrap();
        let tree=tree.refined(&[*tree.leaves().iter().next().unwrap()],proceed).unwrap();
        Model::Adaptive(AdaptiveElasticity3::build_with_material_field(
            domain,&tree,&Top,field,&clamp,options,&mut control).unwrap())
    } else {
        Model::Cartesian(CutElasticity3::build_with_material_field(
            domain,[2,1,1],&Top,field,&clamp,options,&mut control).unwrap())
    }
}
fn point(t:f64)->QueryPoint {
    QueryPoint::new().with_quantity("T",QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,ValueForm::Static)),t).unwrap()
}
fn card()->MaterialCard {
    let mut claims=ClaimSet::new();
    for (name,dims,value,slope) in [
        (YOUNG_MODULUS_PROPERTY,Pressure::DIMS,1000.0,-1.0),
        (POISSON_RATIO_PROPERTY,Dims::NONE,0.2,0.0002),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,Dims([0,0,0,-1,0,0]),0.0005,0.000001),
        ("density",Density::DIMS,DENSITY,0.0),
        ("specific_heat_capacity",CP_DIMS,CP,0.0),
        ("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS,2.0,0.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(name,dims),
            value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,value-20.0*slope),(300.0,value),(360.0,value+60.0*slope)],dims},
            validity:ValidityDomain::unconstrained().with_quantity("T",
                *point(300.0).axis_quantities().get("T").unwrap(),280.0,360.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {
                source:"authored cut-solid reciprocal heat fixture; not measurements".into(),
                license:"test-only".into(),artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {
        chemistry:"synthetic thermal cut solid".into(),phase:"solid".into(),
        process:"analytic source fixture".into(),revision:0,
    },claims,vec![]).unwrap()
}
fn law(card:&MaterialCard)->AffineThermoelastic {
    let requirement=|name,dims,admissibility|ScalarPropertyRequirement::try_new(name,dims,admissibility).unwrap();
    let elastic=resolve_material_state_point(card,&point(300.0),&[
        requirement(YOUNG_MODULUS_PROPERTY,Pressure::DIMS,ScalarAdmissibility::StrictlyPositive),
        requirement(POISSON_RATIO_PROPERTY,Dims::NONE,ScalarAdmissibility::Finite),
    ],MaterialPropertySelection::SingleClaimOnly).unwrap();
    let expansion=resolve_material_state_point(card,&point(300.0),&[
        requirement(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,Dims([0,0,0,-1,0,0]),
            ScalarAdmissibility::Finite),
    ],MaterialPropertySelection::SingleClaimOnly).unwrap();
    AffineThermoelastic::new(card,&elastic,&expansion,ThermoelasticSymmetry::Isotropic,
        ThermoelasticSymmetry::Isotropic,[280.0,360.0],0.08).unwrap()
}
fn source(law:&AffineThermoelastic,t:f64)->IsothermalThermoelasticState {
    law.at_temperature(t).unwrap().isothermal_elasticity(IsothermalThermoelasticInput {
        reference_density_kg_m3:DENSITY,mechanical_material_identity:ContentHash([71;32]),
    }).unwrap()
}
fn caloric(card:&MaterialCard)->ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new_exact(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",CP_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",
        lower:&point(280.0),upper:&point(360.0),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.0).unwrap()
}
fn settings(n:usize)->ThermoelasticEquilibriumSettings3 {
    ThermoelasticEquilibriumSettings3 {
        refresh:Default::default(),maximum_native_dofs:8*n,relative_residual_tolerance:1e-11,
        maximum_iterations:4000,poll_iterations:32,maximum_free_residual_n:1e-7,
        maximum_kinematic_strain:0.08,maximum_points:POINTS,maximum_mass_terms:POINTS,
        maximum_prestress_entries:Some(POINTS),moment_origin_m:[0.0;3],
    }
}
fn config()->CutThermoelasticConfig3 {
    CutThermoelasticConfig3 {
        maximum_cells:1024,maximum_iterations:32,heat_relaxation:1.0,
        heat_tolerance_j:2e-8,energy_tolerance_j:1e-7,maximum_numerical_work_j:10.0,
        exchange:ThermoelasticExchangeLimits3 {
            maximum_cells:1024,maximum_points:POINTS,maximum_workspace_bytes:2_000_000,
        },
    }
}
fn step()->StepConfig {StepConfig {max_dt_s:0.1,energy_tolerance_j:2e-8}}
fn binding<'a>(card:&'a MaterialCard,context:&'a QueryPoint)->ThermalMaterialPoint<'a> {
    ThermalMaterialPoint {material:card,context,temperature_axis:"T",liquid_fraction_axis:None}
}
fn cells<'a>(model:&Model,card:&'a MaterialCard,context:&'a QueryPoint,
    curve:&'a EquilibriumEnthalpyPhaseCurve,law:&AffineThermoelastic)->Vec<CutThermoelasticCell3<'a>> {
    (0..native!(model,cells)).map(|i|CutThermoelasticCell3 {
        cell:CellId(i as u64+10),curve,point:binding(card,context),law:law.clone(),
    }).collect()
}
fn forces(model:&Model)->Vec<f64> {
    let mut f=vec![0.0;model.operator().n()];
    let nodes=native!(model,nodes);
    let tip:Vec<_>=nodes.iter().enumerate().filter(|(_,p)|p[0]==1.0).map(|(i,_)|i).collect();
    for i in tip.iter().copied() {f[3*i+2]=-0.02/(tip.len() as f64);}
    f
}
fn timeline<'a>(model:&Model,curve:&'a EquilibriumEnthalpyPhaseCurve,contact:bool)->DepositionSimulation<'a> {
    let h=curve.knots().iter().find(|k|k.temperature_k==300.0).unwrap().specific_enthalpy_j_kg;
    let mut events:Vec<_>=(0..native!(model,cells)).map(|i| {
        let bounds=native!(model,reference_cell_bounds,i).unwrap();
        let volume=native!(model,reference_cell_volume_m3,i).unwrap();
        let mass=DENSITY*volume;
        let rate=if bounds.lo()[0]<0.5 {1.0}else{3.0};
        TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id:CellId(i as u64+10),segment_id:i as u64,body_id:1,material_id:MaterialId(0),
            mass_kg:mass,specific_h_j_kg:h,
            geometry:CellGeometry {start_m:bounds.lo(),end_m:bounds.hi(),volume_m3:volume,surface_area_m2:3.0},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:mass*CP*rate},
        })}
    }).collect();
    if contact {
        assert_eq!(native!(model,cells),2);
        events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
            id:7,a:CellId(10),b:CellId(11),area_m2:0.73,conductance_w_k:0.0,
        })});
    }
    DepositionSimulation::new(events,&[(MaterialId(0),curve as &dyn ThermalLaw)],
        DepositionBudget {max_cells:1024,max_links:1024,max_events:2048,max_substeps:1000}).unwrap()
        .advance_to(0.0,step(),||Ok(())).unwrap().0
}
#[derive(Clone,Debug,Default,PartialEq)]
struct Observe {steps:usize,duration_s:f64,fail_after_s:Option<f64>}
impl DepositionObserver for Observe {
    type Error=DepositionError;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self,_:f64,dt:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        check()?;self.steps+=1;self.duration_s+=dt;
        if self.fail_after_s.is_some_and(|t|self.duration_s>=t) {
            return Err(DepositionError::InvalidInput("authored late observer failure"));
        }
        Ok(())
    }
}
fn near(a:f64,b:f64,tol:f64) {
    assert!(a.is_finite()&&b.is_finite()&&(a-b).abs()<=tol*a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}");
}
fn with_cx(run:impl FnOnce(&Cx<'_>,&CancelGate)) {
    let gate=CancelGate::new();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:73,kernel_id:112,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        run(&cx,&gate)
    })
}
fn action(model:&Model)->Vec<f64> {
    let x:Vec<_>=(0..model.operator().n()).map(|i|((i+1) as f64).sin()).collect();
    let mut y=vec![0.0;x.len()];model.operator().apply(&x,&mut y);y
}

#[test]
fn original_finite_heating_changes_compatible_raw_and_hanging_equilibrium_with_closed_feedback() {
    for adaptive in [false,true] {
        with_cx(|cx,gate| {
            let card=card();let law=law(&card);let curve=caloric(&card);let context=point(300.0);
            let mut model=build(adaptive,&source(&law,300.0));
            let mut thermal=timeline(&model,curve.curve(),false);
            let loads=forces(&model);
            let original_cells=thermal.cells().clone();
            let original_energy=thermal.energy_j();
            let authored: f64=thermal.cells().values().map(|c|c.boundary.source_w).sum();
            let initial_mass=thermal.mass_kg();
            let mut solid=CutThermoelasticSolid3::try_new(model.source(),&thermal,
                cells(&model,&card,&context,curve.curve(),&law),&loads,
                settings(model.operator().n()),config(),cx).unwrap();
            let old_u=solid.equilibrium().solution.coefficients().to_vec();
            let occupancy=vec![1.0;native!(&model,cells)];
            let original_diagonal=solid.equilibrium().mass.diagonal(&occupancy,proceed).unwrap();
            let mut observer=Observe::default();
            let one_way=thermal.advance_to(1.0,step(),||Ok(())).unwrap().0;
            let report=solid.advance(model.mutable(),&mut thermal,&mut observer,1.0,&loads,
                step(),config(),cx).unwrap();
            assert_eq!(solid.time_s(),thermal.time_s());
            assert_eq!(solid.time_s(),1.0);
            assert!(report.iterations>1,"reversible source must actually be iterated");
            assert_eq!(observer.steps,report.thermal_audit.substeps);
            near(observer.duration_s,1.0,1e-14);
            assert!(solid.equilibrium().solution.coefficients().iter().zip(old_u)
                .any(|(a,b)|(a-b).abs()>1e-6));
            assert!(report.reversible_heat_j.abs()>1e-5);
            assert!(thermal.cells().iter().any(|(id,c)|
                (c.specific_h_j_kg-one_way.cells()[id].specific_h_j_kg).abs()>1e-7));
            assert!(report.heat_defect_j<=config().heat_tolerance_j);
            assert!(report.combined_residual_j.abs()<=config().energy_tolerance_j);
            near(report.external_heat_j,authored,1e-9);
            near(thermal.energy_j()-original_energy+report.internal_energy_change_j
                -report.external_work_j+report.numerical_work_j,authored,2e-9);
            near(report.native_work_j-report.physical_bulk_work_j,report.numerical_work_j,1e-12);
            near(report.physical_bulk_work_j-report.internal_energy_change_j,
                report.reversible_heat_j,1e-12);
            assert_eq!(thermal.mass_kg(),initial_mass);
            for (id,c) in thermal.cells() {
                assert_eq!(c.mass_kg,original_cells[id].mass_kg);
                assert_eq!(c.geometry,original_cells[id].geometry);
            }
            native!(&model,validate_material_mass,&solid.equilibrium().mass).unwrap();
            native!(&model,validate_eigenstrain,&solid.equilibrium().eigenstrain).unwrap();
            assert!(solid.equilibrium().prestress.is_some());
            assert_eq!(solid.equilibrium().mass.diagonal(&occupancy,proceed).unwrap(),original_diagonal);
            {
                let accepted=solid.equilibrium();
                let dynamics=CutSecondOrderProblem3::new_with_eigenstrain(
                    model.source(),&accepted.mass,&occupancy,&accepted.eigenstrain,
                    0.08,POINTS,model.operator().n(),gate).unwrap()
                    .with_reference_prestress(accepted.solution.coefficients(),POINTS,POINTS).unwrap();
                let q=dynamics.restrict_force(accepted.solution.coefficients()).unwrap();
                let expected=dynamics.restrict_force(&loads).unwrap();
                let mut resistance=vec![0.0;q.len()];
                dynamics.internal_force(&q,&mut resistance);
                for (a,b) in resistance.iter().zip(expected) {near(*a,b,2e-7);}
                let direction:Vec<_>=(0..model.operator().n()).map(|i|(i as f64).sin()).collect();
                let mut packet=vec![0.0;direction.len()];let mut held=packet.clone();
                accepted.prestress.as_ref().unwrap().apply_into(&direction,&mut packet,proceed).unwrap();
                dynamics.reference_prestress().unwrap().apply_into(&direction,&mut held,proceed).unwrap();
                assert_eq!(packet,held);
            }
            // Continue the SAME thermal/observer history; accepted retries were not emitted twice.
            let prior_steps=observer.steps;
            let report2=solid.advance(model.mutable(),&mut thermal,&mut observer,2.0,&loads,
                step(),config(),cx).unwrap();
            assert_eq!(observer.steps,prior_steps+report2.thermal_audit.substeps);
            assert_eq!(solid.time_s(),2.0);
            near(observer.duration_s,2.0,1e-14);
        });
    }
}

#[test]
fn valid_picard_rejection_and_late_mechanical_failure_restore_every_accepted_owner() {
    with_cx(|cx,gate| {
        let card=card();let law=law(&card);let curve=caloric(&card);let context=point(300.0);
        let mut model=build(true,&source(&law,300.0));
        let mut thermal=timeline(&model,curve.curve(),false);
        let loads=forces(&model);
        let mut solid=CutThermoelasticSolid3::try_new(model.source(),&thermal,
            cells(&model,&card,&context,curve.curve(),&law),&loads,
            settings(model.operator().n()),config(),cx).unwrap();
        let original=thermal.cells().clone();let original_action=action(&model);
        let original_u=solid.equilibrium().solution.coefficients().to_vec();
        let mut observer=Observe::default();
        let once=CutThermoelasticConfig3 {maximum_iterations:1,heat_tolerance_j:1e-25,..config()};
        assert!(matches!(solid.advance(model.mutable(),&mut thermal,&mut observer,1.0,&loads,
            step(),once,cx),Err(CutThermoelasticError3::NotConverged {..})));
        assert_eq!(thermal.time_s(),0.0);assert_eq!(solid.time_s(),0.0);
        assert_eq!(thermal.cells(),&original);assert_eq!(observer,Observe::default());
        assert_eq!(action(&model),original_action);
        assert_eq!(solid.equilibrium().solution.coefficients(),original_u);
        native!(&model,validate_material_mass,&solid.equilibrium().mass).unwrap();
        native!(&model,validate_eigenstrain,&solid.equilibrium().eigenstrain).unwrap();
        let mut late=Observe {fail_after_s:Some(0.95),..Default::default()};
        let untouched_late=late.clone();
        assert!(matches!(solid.advance(model.mutable(),&mut thermal,&mut late,1.0,&loads,
            step(),config(),cx),Err(CutThermoelasticError3::Observer(_))));
        assert_eq!(late,untouched_late);assert_eq!(thermal.cells(),&original);
        assert_eq!(action(&model),original_action);assert_eq!(solid.time_s(),0.0);
        let excessive:Vec<_>=loads.iter().map(|f|f*1e10).collect();
        assert!(solid.advance(model.mutable(),&mut thermal,&mut observer,1.0,&excessive,
            step(),config(),cx).is_err());
        assert_eq!(thermal.cells(),&original);assert_eq!(observer,Observe::default());
        assert_eq!(action(&model),original_action);assert_eq!(solid.time_s(),0.0);
        gate.request();
        assert!(matches!(solid.advance(model.mutable(),&mut thermal,&mut observer,1.0,&loads,
            step(),config(),cx),Err(CutThermoelasticError3::Cancelled)));
        assert_eq!(thermal.cells(),&original);assert_eq!(observer,Observe::default());
        assert_eq!(action(&model),original_action);
    });
}

#[test]
fn numerical_stabilization_work_is_reported_and_gated_without_becoming_heat() {
    with_cx(|cx,_| {
        let card=card();let law=law(&card);let curve=caloric(&card);let context=point(300.0);
        let mut model=build(true,&source(&law,300.0));
        let mut thermal=timeline(&model,curve.curve(),false);
        let loads=forces(&model);let original=thermal.cells().clone();
        let mut solid=CutThermoelasticSolid3::try_new(model.source(),&thermal,
            cells(&model,&card,&context,curve.curve(),&law),&loads,
            settings(model.operator().n()),config(),cx).unwrap();
        let mut observer=Observe::default();
        let strict=CutThermoelasticConfig3 {maximum_numerical_work_j:0.0,..config()};
        assert!(matches!(solid.advance(model.mutable(),&mut thermal,&mut observer,1.0,&loads,
            step(),strict,cx),Err(CutThermoelasticError3::NumericalWork {..})));
        assert_eq!(thermal.cells(),&original);assert_eq!(solid.time_s(),0.0);
        let report=solid.advance(model.mutable(),&mut thermal,&mut observer,1.0,&loads,
            step(),config(),cx).unwrap();
        assert!(report.numerical_work_j.abs()>1e-14);
        near(report.reversible_heat_j,report.cells.iter()
            .map(|c|c.exchange.retained_reversible_heat_j().to_f64()).sum(),1e-12);
        near(report.physical_bulk_work_j-report.internal_energy_change_j,
            report.reversible_heat_j,1e-12);
    });
}

/// The original observer alone panics when displaced; its candidate clones do not.
#[derive(Debug)]
struct DropObserver {state:Observe,panic_on_drop:bool}
impl Clone for DropObserver {
    fn clone(&self)->Self {Self {state:self.state.clone(),panic_on_drop:false}}
}
impl Drop for DropObserver {
    fn drop(&mut self) {if self.panic_on_drop {panic!("authored displaced observer destructor");}}
}
impl DepositionObserver for DropObserver {
    type Error=DepositionError;
    fn event(&mut self,e:&Event,t:&DepositionSimulation<'_>)->Result<(),Self::Error>{self.state.event(e,t)}
    fn step(&mut self,s:f64,dt:f64,old:&BTreeMap<CellId,f64>,t:&DepositionSimulation<'_>,
        c:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error>{
        self.state.step(s,dt,old,t,c)
    }
}
#[test]
fn observer_destructor_panic_occurs_only_after_all_four_accepted_owners_are_replaced() {
    with_cx(|cx,_| {
        let card=card();let law=law(&card);let curve=caloric(&card);let context=point(300.0);
        let mut model=build(false,&source(&law,300.0));
        let mut thermal=timeline(&model,curve.curve(),false);let loads=forces(&model);
        let mut solid=CutThermoelasticSolid3::try_new(model.source(),&thermal,
            cells(&model,&card,&context,curve.curve(),&law),&loads,
            settings(model.operator().n()),config(),cx).unwrap();
        let mut observer=DropObserver {state:Observe::default(),panic_on_drop:true};
        let caught=std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            solid.advance(model.mutable(),&mut thermal,&mut observer,1.0,&loads,
                step(),config(),cx).unwrap();
        }));
        assert!(caught.is_err());
        assert_eq!(solid.time_s(),1.0);assert_eq!(thermal.time_s(),1.0);
        assert!(!observer.panic_on_drop);near(observer.state.duration_s,1.0,1e-14);
        native!(&model,validate_material_mass,&solid.equilibrium().mass).unwrap();
        native!(&model,validate_eigenstrain,&solid.equilibrium().eigenstrain).unwrap();
        for (i,&t) in solid.temperatures_k().iter().enumerate() {
            assert_eq!(native!(&model,cell_thermoelastic_state,i).unwrap().source().temperature_k(),t);
            assert_eq!(thermal.temperature_k(CellId(i as u64+10)).unwrap(),t);
        }
    });
}

#[path = "cut_thermoelastic/network.rs"]
mod network;

#[path = "cut_thermoelastic/dynamic.rs"]
mod dynamic;

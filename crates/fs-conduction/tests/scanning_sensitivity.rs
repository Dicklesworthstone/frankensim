//! Scan controls through original nonlinear enthalpy/radiation adjoints.
mod support;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory,ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun,EnthalpyHistoryMaterials,EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::{ScanningEnthalpyRun,ScanRadiation};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::ScanSensitivityBudget;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget,MovingHeatSource};
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::transient::enthalpy::{EnthalpyBudget,EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh,ConductivityModel,ConductivityTable,ScalarField,ThermalBc,
    ThermalBoundaryBuilder,LinearConfig,AmbientRadiationConfig,SurfaceEmissivity,
    COMPONENT_POWER_DIMS,TEMPERATURE_DIMS,EMISSIVITY_DIMS,SURFACE_EMISSIVITY_PROPERTY};
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve,EnthalpyPhaseKnot};
use fs_matdb::{ClaimSet,PropertyClaim,PropertyKey,PropertyValue,UncertaintyModel,InterpolationPolicy,
    Provenance,MaterialCard,MaterialStateId,SelectionPolicy};
use fs_qty::{Dims,QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal,Interp};
use fs_solver::NewtonKrylovConfig;
use support::{with_cx,with_cancelled_cx};
fn constant(v:f64,d:Dims)->TimeSignal{TimeSignal::Constant(QtyAny::new(v,d))}
fn mesh()->ConductionMesh{ConductionMesh::new(TetComplex::from_tets(4,vec![[0,1,2,3]]),
    vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).unwrap()}
fn chart()->EquilibriumEnthalpyPhaseCurve{
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([74;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot{specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect()).unwrap()
}
fn config()->EnthalpyStepConfig{EnthalpyStepConfig{newton:NewtonKrylovConfig{
    absolute_tolerance:1e-12,relative_tolerance:1e-12,linear_restart:4,max_linear_cycles:16,
    forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
    max_newton_iterations:64,energy_tolerance_j:1e-9}}
fn linear()->LinearConfig{LinearConfig{tolerance:1e-11,max_iterations:1000,restart:4}}
fn budget()->ScanSensitivityBudget{ScanSensitivityBudget{max_source_vertices:32,max_radiation_feedback_entries:1024}}
fn heater(p:[f64;7])->MovingHeatSource{MovingHeatSource::new(
    [p[0],p[1],p[2]].map(|v|constant(v,Dims([1,0,0,0,0,0]))),[p[3],p[4],p[5]],
    constant(p[6],COMPONENT_POWER_DIMS),0.,MovingHeatBudget{max_vertices:4,max_elements:1,
        min_active_vertices:2,power_relative_tolerance:1e-12}).unwrap()}
fn history<'m>(cx:&Cx<'_>,m:&'m ConductionMesh)->ThermalBoundaryHistory<'m>{
    ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(m).remainder("surface",ThermalBc::adiabatic()).unwrap(),
        &Environment::earth_lab(),vec![],0.).unwrap()
}
fn radiation()->ScanRadiation{
    let mut claims=ClaimSet::new();claims.insert_claim(PropertyClaim{
        key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar{value:0.2,dims:EMISSIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.),uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
        provenance:Provenance{source:"synthetic scan gradient fixture".into(),license:"test-only".into(),artifact:None}}).unwrap();
    let card=MaterialCard::assemble(MaterialStateId{chemistry:"gray fixture".into(),phase:"solid".into(),
        process:"synthetic".into(),revision:0},claims,vec![]).unwrap();
    ScanRadiation{patches:vec![RadiationRegionHistory{region:"surface".into(),
        emissivity:SurfaceEmissivity::from_card("surface",&card,300.,SelectionPolicy::SingleClaimOnly).unwrap(),
        ambient_temperature:constant(293.,TEMPERATURE_DIMS)}],config:AmbientRadiationConfig::default()}
}
fn run<'h,'m,'c>(cx:&Cx<'_>,history:&'h ThermalBoundaryHistory<'m>,material:&'h ConductivityModel,
    curve:&'c EquilibriumEnthalpyPhaseCurve,old:&[f64],sources:Vec<MovingHeatSource>,rad:bool)
    ->ScanningEnthalpyRun<'h,'m,'c>{
    let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem{boundary:history,material,
        element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
        EnthalpyHistoryMaterials{table:&[ReferenceEnthalpyMaterial{curve,reference_density_kg_m3:1.}],
            of_element:&[0],phase_conductivity:None},old,
        EnthalpyHistoryBudget{spatial:EnthalpyBudget{max_vertices:4,max_elements:1},max_steps:8}).unwrap();
    ScanningEnthalpyRun::new(cx,base,sources,4,rad.then(radiation)).unwrap()
}
fn flatten(g:&fs_conduction::bc::history::march::moving::sensitivity::MovingHeatGradient)->[f64;7]{
    [g.center_m[0],g.center_m[1],g.center_m[2],g.half_width_m[0],g.half_width_m[1],g.half_width_m[2],g.absorbed_power_w]
}
fn close(a:f64,b:f64){assert!((a-b).abs()<3e-5*(1.+a.abs().max(b.abs())),"{a:e} != {b:e}");}
#[test]
fn g3_scan_gradients_cross_nonlinear_conductivity_and_radiation_without_resolving_each_control(){
    let m=mesh();let curve=chart();let mat=ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280.,1.),(380.,3.)]).unwrap());
    with_cx(|cx|{
        let hist=history(cx,&m);let p=[0.2,0.25,0.3,1.4,1.3,1.2,3.];let old=[10.,20.,30.,40.];
        let seed=[0.3,-0.1,0.7,0.5];
        for rad in [false,true]{
            let mut evolving=run(cx,&hist,&mat,&curve,&old,vec![heater(p)],rad);
            let endpoint=evolving.advance_linearized_to(cx,0.002,config(),budget()).unwrap();
            let g=endpoint.pullback_temperature(cx,&seed,linear()).unwrap();let a=flatten(&g.sources[0]);
            assert!(g.storage.transport.relative_residual<=linear().tolerance);
            let mut direct=run(cx,&hist,&mat,&curve,&old,vec![heater(p)],rad);
            assert_eq!(endpoint.report().solution().specific_enthalpy_j_kg,
                direct.advance_to(cx,0.002,config()).unwrap().solution().specific_enthalpy_j_kg);
            let objective=|p|{
                let mut trial=run(cx,&hist,&mat,&curve,&old,vec![heater(p)],rad);
                let step=trial.advance_to(cx,0.002,config()).unwrap();
                seed.iter().zip(&step.solution().temperature).map(|(s,t)|s*t).sum::<f64>()
            };
            for i in 0..7{let(mut plus,mut minus)=(p,p);let eps=2e-5;plus[i]+=eps;minus[i]-=eps;
                close(a[i],(objective(plus)-objective(minus))/(2.*eps));}
        }
    });
}
#[test]
fn g3_fixed_grid_reverse_sweep_includes_earlier_heat_and_independent_sources(){
    let m=mesh();let curve=chart();let mat=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx|{
        let hist=history(cx,&m);let p=[0.2,0.25,0.3,1.4,1.3,1.2,3.];let old=[10.,20.,30.,40.];
        let second=[0.6,0.1,0.2,1.5,1.4,1.3,1.5];let seed=[0.3,-0.1,0.7,0.5];
        let mut evolving=run(cx,&hist,&mat,&curve,&old,vec![heater(p),heater(second)],false);
        let first=evolving.advance_linearized_to(cx,0.002,config(),budget()).unwrap();
        let last=evolving.advance_linearized_to(cx,0.005,config(),budget()).unwrap();
        let last_g=last.pullback(cx,&seed,linear()).unwrap();
        let first_g=first.pullback(cx,&last_g.storage.transport.previous_specific_enthalpy,linear()).unwrap();
        assert_eq!(last_g.sources.len(),2);
        let objective=|p|{
            let mut trial=run(cx,&hist,&mat,&curve,&old,vec![heater(p),heater(second)],false);
            trial.advance_to(cx,0.002,config()).unwrap();
            let end=trial.advance_to(cx,0.005,config()).unwrap();
            seed.iter().zip(&end.solution().specific_enthalpy_j_kg).map(|(s,h)|s*h).sum::<f64>()
        };
        let a=flatten(&first_g.sources[0]);let b=flatten(&last_g.sources[0]);
        for i in 0..7{let(mut plus,mut minus)=(p,p);let eps=2e-5;plus[i]+=eps;minus[i]-=eps;
            close(a[i]+b[i],(objective(plus)-objective(minus))/(2.*eps));}
    });
}
#[test]
fn g1_latent_plateau_has_zero_temperature_but_nonzero_energy_source_sensitivity(){
    let m=mesh();let curve=chart();let mat=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx|{
        let hist=history(cx,&m);let p=[0.2,0.25,0.3,1.4,1.3,1.2,3.];
        let mut evolving=run(cx,&hist,&mat,&curve,&[150.;4],vec![heater(p)],false);
        let endpoint=evolving.advance_linearized_to(cx,0.002,config(),budget()).unwrap();
        assert!(endpoint.report().solution().temperature.iter().all(|&t|t==310.));
        let t=endpoint.pullback_temperature(cx,&[1.;4],linear()).unwrap();
        assert!(flatten(&t.sources[0]).iter().all(|&v|v==0.));
        let energy=endpoint.pullback(cx,evolving.conduction().storage().reference_nodal_masses_kg(),linear()).unwrap();
        close(energy.sources[0].absorbed_power_w,0.002);
        for i in 0..3{close(energy.sources[0].center_m[i],0.);close(energy.sources[0].half_width_m[i],0.);}
    });
}
#[test]
fn g0_derivative_refusals_are_atomic_and_event_binding_keeps_the_primal_side(){
    let m=mesh();let curve=chart();let mat=ConductivityModel::isotropic_declared(2.).unwrap();
    let hist=with_cx(|cx|history(cx,&m));let p=[0.2,0.25,0.3,1.4,1.3,1.2,3.];
    let mut evolving=with_cx(|cx|run(cx,&hist,&mat,&curve,&[20.;4],vec![heater(p)],false));
    with_cancelled_cx(|cx|assert!(evolving.advance_linearized_to(cx,0.002,config(),budget()).is_err()));
    with_cx(|cx|{
        assert!(evolving.advance_linearized_to(cx,0.002,config(),ScanSensitivityBudget{max_source_vertices:0,..budget()}).is_err());
        assert_eq!(evolving.conduction().accepted_steps(),0);assert_eq!(evolving.conduction().specific_enthalpy_j_kg(),&[20.;4]);
        let mut off=p;off[6]=0.;
        let mut kink=run(cx,&hist,&mat,&curve,&[100.;4],vec![heater(off)],false);
        assert!(kink.advance_linearized_to(cx,0.002,config(),budget()).is_err());
        assert_eq!(kink.conduction().accepted_steps(),0);
        assert_eq!(kink.conduction().specific_enthalpy_j_kg(),&[100.;4]);
        kink.advance_to(cx,0.002,config()).unwrap();
        let held=MovingHeatSource::new([p[0],p[1],p[2]].map(|v|constant(v,Dims([1,0,0,0,0,0]))),
            [p[3],p[4],p[5]],TimeSignal::Table{times:vec![0.,0.01],values:vec![3.,0.],
                dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,heater(p).budget()).unwrap();
        let mut events=run(cx,&hist,&mat,&curve,&[20.;4],vec![held],false);
        let endpoint=events.advance_linearized_event_step(cx,0.02,0.02,config(),budget()).unwrap();
        assert_eq!(endpoint.end_time_s(),0.01);assert_eq!(endpoint.report().deposition[0].absorbed_power_w,3.);
        assert!(endpoint.pullback(cx,&[],linear()).is_err());
        assert_eq!(events.conduction().accepted_steps(),1);
    });
}

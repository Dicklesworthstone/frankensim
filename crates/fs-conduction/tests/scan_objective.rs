//! Scan controls through original nonlinear enthalpy/radiation adjoints.
mod support;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory,ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun,EnthalpyHistoryMaterials,EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::{ScanningEnthalpyRun,ScanRadiation};
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
fn close(a:f64,b:f64){assert!((a-b).abs()<3e-5*(1.+a.abs().max(b.abs())),"{a:e} != {b:e}");}

use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    evaluate_scan_objective, ScanObservable, ScanObservation, ScanObjectiveBudget,
    ScanObjectiveSampling,
};

fn objective_budget() -> ScanObjectiveBudget {
    ScanObjectiveBudget { max_steps: 8, max_observations: 32, max_retained_vertices: 32,
        max_retained_source_vertices: 128, max_schedule_controls: 64,
        max_retained_radiation_entries: 8192 }
}
fn powers(values: [f64; 3], interp: Interp) -> MovingHeatSource {
    let template = heater([0.2,0.25,0.3,1.4,1.3,1.2,1.]);
    MovingHeatSource::new([0.2,0.25,0.3].map(|v|constant(v,Dims([1,0,0,0,0,0]))),
        [1.4,1.3,1.2], TimeSignal::Table { times: vec![0.,0.002,0.006],
            values: values.to_vec(), dims: COMPONENT_POWER_DIMS, interp }, 0., template.budget()).unwrap()
}
fn observation(endpoint: usize, vertex: usize, observable: ScanObservable, target: f64, scale: f64)
    -> ScanObservation { ScanObservation { endpoint, vertex, observable, target, scale } }

#[test]
fn g3_full_objective_gradient_includes_intermediate_data_earlier_heat_and_both_sources() {
    let m=mesh();let curve=chart();
    let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.),(380.,3.)]).unwrap());
    with_cx(|cx| {
        let hist=history(cx,&m);let times=[0.001,0.002,0.004,0.006];
        let observations=[observation(0,0,ScanObservable::Temperature,301.2,2.),
            observation(1,2,ScanObservable::SpecificEnthalpy,31.,4.),
            observation(3,1,ScanObservable::Temperature,304.,0.5),
            observation(3,3,ScanObservable::SpecificEnthalpy,42.,3.)];
        for rad in [false,true] { for interp in [Interp::Linear,Interp::Hold] {
            let evaluate=|a,b,old:&[f64]| evaluate_scan_objective(cx,
                run(cx,&hist,&material,&curve,old,vec![powers(a,interp),powers(b,interp)],rad),
                &times,&observations,ScanObjectiveSampling::EventLeft,config(),linear(),objective_budget()).unwrap();
            let a=[2.,4.,1.];let b=[0.5,1.,2.];let old=[10.,20.,30.,40.];
            let actual=evaluate(a,b,&old);
            assert_eq!(actual.endpoint_count,4);
            assert!(actual.max_transpose_relative_residual<=linear().tolerance);
            let objective=|a,b,old:&[f64]| {
                let mut forward=run(cx,&hist,&material,&curve,old,vec![powers(a,interp),powers(b,interp)],rad);
                let mut loss=0.;
                for (i,&time) in times.iter().enumerate() {
                    let solved=forward.advance_event_step(cx,time,time-forward.conduction().time_s(),config()).unwrap();
                    assert_eq!(solved.end_time_s,time);
                    for observation in observations.iter().filter(|o|o.endpoint==i) {
                        let predicted=match observation.observable {
                            ScanObservable::Temperature=>solved.report.solution().temperature[observation.vertex],
                            ScanObservable::SpecificEnthalpy=>solved.report.solution().specific_enthalpy_j_kg[observation.vertex],
                        };
                        loss+=0.5*((predicted-observation.target)/observation.scale).powi(2);
                    }
                }
                loss
            };
            close(actual.value,objective(a,b,&old));
            for source in 0..2 { for index in 0..3 {
                let(mut ap,mut am,mut bp,mut bm)=(a,a,b,b);let eps=1e-4;
                if source==0 {ap[index]+=eps;am[index]-=eps;} else {bp[index]+=eps;bm[index]-=eps;}
                close(actual.sources[source].absorbed_power_values_w[index],
                    (objective(ap,bp,&old)-objective(am,bm,&old))/(2.*eps));
            }}
            for vertex in 0..4 {
                let(mut plus,mut minus)=(old,old);let eps=1e-4;plus[vertex]+=eps;minus[vertex]-=eps;
                close(actual.initial_specific_enthalpy[vertex],
                    (objective(a,b,&plus)-objective(a,b,&minus))/(2.*eps));
            }
        }}
    });
}

#[test]
fn g1_plateau_targets_recover_authored_pulse_gradients_without_temperature_signal() {
    let m=mesh();let curve=chart();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let hist=history(cx,&m);
        let source=MovingHeatSource::new([constant(0.5,Dims([1,0,0,0,0,0])),
            constant(0.5,Dims([1,0,0,0,0,0])),constant(0.5,Dims([1,0,0,0,0,0]))],
            [2.;3],TimeSignal::Table {times:vec![0.,0.002,0.004],values:vec![2.,5.,0.],
                dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,heater([0.5,0.5,0.5,2.,2.,2.,1.]).budget()).unwrap();
        let solved=evaluate_scan_objective(cx,run(cx,&hist,&material,&curve,&[150.;4],vec![source],false),
            &[0.002,0.004],&[observation(0,0,ScanObservable::SpecificEnthalpy,150.,1.),
                observation(1,0,ScanObservable::SpecificEnthalpy,150.,1.),
                observation(1,0,ScanObservable::Temperature,311.,1.)],
            ScanObjectiveSampling::EventLeft,config(),linear(),objective_budget()).unwrap();
        assert!((solved.value-0.5*(0.024_f64.powi(2)+0.084_f64.powi(2)+1.)).abs()<1e-10);
        let power=&solved.sources[0].absorbed_power_values_w;
        assert!((power[0]-0.001296).abs()<1e-10);
        assert!((power[1]-0.001008).abs()<1e-10);
        assert_eq!(power[2],0.);
        assert!((solved.initial_specific_enthalpy[0]-0.108).abs()<1e-9);
        assert!(solved.final_step.solution().temperature.iter().all(|&t|t==310.));
    });
}

#[test]
fn g0_objective_rejects_bad_grid_observations_storage_limits_and_cancellation() {
    let m=mesh();let curve=chart();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let hist=with_cx(|cx|history(cx,&m));
    let fresh=|cx:&Cx<'_>|run(cx,&hist,&material,&curve,&[20.;4],vec![powers([2.,4.,1.],Interp::Hold)],false);
    let observations=[observation(0,0,ScanObservable::Temperature,300.,1.)];
    let input=with_cx(|cx|fresh(cx));
    with_cancelled_cx(|cx|assert!(evaluate_scan_objective(cx,input,&[0.001],&observations,
        ScanObjectiveSampling::EventLeft,config(),linear(),objective_budget()).is_err()));
    with_cx(|cx| {
        let evaluate=|times:&[f64],obs:&[ScanObservation],budget| evaluate_scan_objective(cx,fresh(cx),times,obs,
            ScanObjectiveSampling::EventLeft,config(),linear(),budget);
        assert!(evaluate(&[0.004],&observations,objective_budget()).is_err());
        for times in [vec![],vec![f64::NAN],vec![0.],vec![0.001,0.001]] {
            assert!(evaluate(&times,&observations,objective_budget()).is_err());
        }
        for obs in [observation(1,0,ScanObservable::Temperature,300.,1.),
            observation(0,4,ScanObservable::Temperature,300.,1.),
            observation(0,0,ScanObservable::Temperature,300.,0.)] {
            assert!(evaluate(&[0.001],&[obs],objective_budget()).is_err());
        }
        for budget in [ScanObjectiveBudget{max_retained_vertices:0,..objective_budget()},
            ScanObjectiveBudget{max_retained_source_vertices:0,..objective_budget()},
            ScanObjectiveBudget{max_schedule_controls:0,..objective_budget()}] {
            assert!(evaluate(&[0.001],&observations,budget).is_err());
        }
        assert!(evaluate(&[0.001],&observations,objective_budget()).is_ok());
    });
}

//! Independent thermal response constraints through the original scan adjoints.
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
mod response_rows {
    use super::*;
    use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
        evaluate_scan_objective, ScanObjectiveBudget, ScanObjectiveSampling, ScanObservable, ScanObservation,
        responses::{evaluate_scan_objective_with_responses, ScanResponse, ScanResponseBudget},
    };
    fn objective_budget() -> ScanObjectiveBudget {
        ScanObjectiveBudget { max_steps: 2, max_observations: 8, max_retained_vertices: 8,
            max_retained_source_vertices: 8, max_schedule_controls: 16,
            max_retained_radiation_entries: 2048 }
    }
    fn response_budget() -> ScanResponseBudget {
        ScanResponseBudget { max_responses: 6, max_gradient_entries: 128, max_transpose_solves: 16 }
    }
    fn observations() -> Vec<ScanObservation> {
        vec![ScanObservation { endpoint: 0, vertex: 1, observable: ScanObservable::Temperature,
            target: 302.0, scale: 2.0 },
            ScanObservation { endpoint: 1, vertex: 2, observable: ScanObservable::SpecificEnthalpy,
                target: 32.0, scale: 10.0 }]
    }
    fn components(g: &fs_conduction::bc::history::march::moving::sensitivity::schedule::MovingHeatScheduleGradient)
        -> [f64; 7]
    { [g.center_values_m[0][0], g.center_values_m[1][0], g.center_values_m[2][0],
        g.half_width_m[0], g.half_width_m[1], g.half_width_m[2], g.absorbed_power_values_w[0]] }

    #[test]
    fn g3_response_rows_reuse_the_objective_primal_and_match_nonlinear_radiative_differences() {
        let m = mesh(); let curve = chart();
        let mat = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.),(380.,3.)]).unwrap());
        with_cx(|cx| {
            let hist = history(cx, &m); let old = [10.,20.,30.,40.];
            let p = [0.2,0.25,0.3,1.4,1.3,1.2,3.]; let endpoints = [0.002,0.005];
            let requests = [ScanResponse { endpoint: 0, vertex: 0, observable: ScanObservable::Temperature },
                ScanResponse { endpoint: 1, vertex: 3, observable: ScanObservable::SpecificEnthalpy }];
            for rad in [false,true] {
                let make = |p| run(cx, &hist, &mat, &curve, &old, vec![heater(p)], rad);
                let (actual, rows) = evaluate_scan_objective_with_responses(cx, make(p), &endpoints,
                    &observations(), ScanObjectiveSampling::RightEndpoint, config(), linear(), objective_budget(),
                    &requests, response_budget()).unwrap();
                let original = evaluate_scan_objective(cx, make(p), &endpoints, &observations(),
                    ScanObjectiveSampling::RightEndpoint, config(), linear(), objective_budget()).unwrap();
                assert_eq!(actual.value, original.value);
                assert_eq!(actual.predictions, original.predictions);
                assert_eq!(actual.sources, original.sources);
                assert_eq!(actual.final_step.solution().specific_enthalpy_j_kg,
                    original.final_step.solution().specific_enthalpy_j_kg);
                let sample = |p| {
                    let mut trial = make(p);
                    let first = trial.advance_to(cx, endpoints[0], config()).unwrap();
                    let last = trial.advance_to(cx, endpoints[1], config()).unwrap();
                    [first.solution().temperature[0], last.solution().specific_enthalpy_j_kg[3]]
                };
                assert_eq!(rows.iter().map(|r| r.value).collect::<Vec<_>>(), sample(p));
                for coordinate in 0..7 {
                    let (mut plus, mut minus) = (p,p); let eps = 2e-5;
                    plus[coordinate] += eps; minus[coordinate] -= eps;
                    let (a,b) = (sample(plus),sample(minus));
                    for (i,row) in rows.iter().enumerate() {
                        let expected = (a[i]-b[i])/(2.0*eps);
                        let got = components(&row.sources[0])[coordinate];
                        assert!((got-expected).abs() < 2e-6*(1.0+expected.abs()), "{got:e} != {expected:e}");
                        assert!(row.max_transpose_relative_residual <= linear().tolerance);
                    }
                }
            }
        });
    }

    #[test]
    fn g1_latent_responses_keep_heat_sensitivity_and_do_not_assign_it_to_future_controls() {
        let m = mesh(); let curve = chart(); let mat = ConductivityModel::isotropic_declared(2.).unwrap();
        with_cx(|cx| {
            let hist = history(cx,&m); let p = [0.2,0.25,0.3,1.4,1.3,1.2,3.];
            let source = MovingHeatSource::new([p[0],p[1],p[2]].map(|v|constant(v,Dims([1,0,0,0,0,0]))),
                [p[3],p[4],p[5]],TimeSignal::Table { times:vec![0.,0.002,0.005],values:vec![3.,0.5,0.],
                    dims:COMPONENT_POWER_DIMS,interp:Interp::Hold },0.,heater(p).budget()).unwrap();
            let mut requests = vec![ScanResponse { endpoint:0,vertex:0,observable:ScanObservable::Temperature },
                ScanResponse { endpoint:0,vertex:0,observable:ScanObservable::SpecificEnthalpy }];
            requests.extend((0..4).map(|vertex|ScanResponse { endpoint:1,vertex,observable:ScanObservable::SpecificEnthalpy }));
            let (_,rows) = evaluate_scan_objective_with_responses(cx,
                run(cx,&hist,&mat,&curve,&[150.;4],vec![source],false),&[0.002,0.005],&observations(),
                ScanObjectiveSampling::EventLeft,config(),linear(),objective_budget(),&requests,response_budget()).unwrap();
            assert_eq!(rows[0].value,310.);
            assert!(rows[0].sources[0].absorbed_power_values_w.iter().all(|&g|g==0.));
            assert!(rows[0].initial_specific_enthalpy.iter().all(|&g|g==0.));
            assert!(rows[1].sources[0].absorbed_power_values_w[0]>0.);
            assert_eq!(rows[1].sources[0].absorbed_power_values_w[1],0.);
            for (control,expected) in [0.002,0.003,0.].into_iter().enumerate() {
                let total: f64 = rows[2..].iter().map(|r|r.sources[0].absorbed_power_values_w[control]/24.).sum();
                assert!((total-expected).abs()<1e-10);
            }
        });
    }

    #[test]
    fn g0_response_shape_work_storage_and_cancellation_refuse_without_returning_partial_rows() {
        let m = mesh(); let curve = chart(); let mat = ConductivityModel::isotropic_declared(2.).unwrap();
        let hist = with_cx(|cx|history(cx,&m));
        let request = ScanResponse { endpoint:1,vertex:0,observable:ScanObservable::SpecificEnthalpy };
        let check = |cx: &Cx<'_>,request,b| {
            evaluate_scan_objective_with_responses(cx,run(cx,&hist,&mat,&curve,&[20.;4],
                vec![heater([0.2,0.25,0.3,1.4,1.3,1.2,3.])],false),&[0.002,0.005],&observations(),
                ScanObjectiveSampling::RightEndpoint,config(),linear(),objective_budget(),&[request],b)
        };
        with_cx(|cx| {
            for b in [ScanResponseBudget { max_responses:0,..response_budget() },
                ScanResponseBudget { max_gradient_entries:0,..response_budget() },
                ScanResponseBudget { max_transpose_solves:3,..response_budget() }]
            { assert!(check(cx,request,b).is_err()); }
            assert!(check(cx,ScanResponse { vertex:4,..request },response_budget()).is_err());
            assert!(check(cx,ScanResponse { endpoint:2,..request },response_budget()).is_err());
        });
        // Construct under a live scope: the cancelled call under test, not its
        // run constructor, is responsible for the refusal.
        let ready = with_cx(|cx|run(cx,&hist,&mat,&curve,&[20.;4],
            vec![heater([0.2,0.25,0.3,1.4,1.3,1.2,3.])],false));
        with_cancelled_cx(|cx| assert!(evaluate_scan_objective_with_responses(cx,ready,
            &[0.002,0.005],&observations(),ScanObjectiveSampling::RightEndpoint,config(),linear(),
            objective_budget(),&[request],response_budget()).is_err()));
    }
}

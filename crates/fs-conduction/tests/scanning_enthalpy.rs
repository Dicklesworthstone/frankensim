//! End-to-end scanning, phase, source superposition and radiation tests.
mod support;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::{ThermalBoundaryHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::moving::motion::ScanMotionLimit;
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun, EnthalpyHistoryMaterials, EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::{ScanningEnthalpyRun, ScanRadiation, ScanTransport};
use fs_conduction::bc::history::march::enthalpy::adaptive::PhaseErrorTolerance;
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, ScalarField,
    ThermalBc, ThermalBoundaryBuilder, COMPONENT_POWER_DIMS, TEMPERATURE_DIMS};
use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal, Interp};
use fs_solver::NewtonKrylovConfig;
use support::{with_cx, with_cancelled_cx};

const LENGTH: Dims = Dims([1,0,0,0,0,0]);
fn constant(value: f64, dims: Dims) -> TimeSignal { TimeSignal::Constant(QtyAny::new(value,dims)) }
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(8,vec![[0,1,2,3],[4,5,6,7]]),
        vec![[0.,0.,0.],[0.5,0.,0.],[0.,0.5,0.],[0.,0.,0.5],
            [1.,0.,0.],[1.5,0.,0.],[1.,0.5,0.],[1.,0.,0.5]]).unwrap()
}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([91;32]),
        [(0.,300.,0.),(100.,310.,0.),(200.,310.,1.),(1000.,390.,1.)]
        .into_iter().map(|(h,t,f)|EnthalpyPhaseKnot { specific_enthalpy_j_kg:h,
            temperature_k:t,liquid_mass_fraction:f,bulk_density_kg_m3:1. }).collect()).unwrap()
}
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11,
        relative_tolerance:1e-12,linear_restart:8,max_linear_cycles:16,
        forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default() },
        max_newton_iterations:64,energy_tolerance_j:1e-8 }
}
fn heater(x: TimeSignal, power: TimeSignal) -> MovingHeatSource {
    MovingHeatSource::new([x,constant(0.1,LENGTH),constant(0.1,LENGTH)], [0.4,0.6,0.6],
        power,0.,MovingHeatBudget {max_vertices:8,max_elements:2,min_active_vertices:1,
            power_relative_tolerance:1e-12}).unwrap()
}
fn history<'m>(cx:&Cx<'_>,mesh:&'m ConductionMesh)->ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(mesh)
        .remainder("surface",ThermalBc::adiabatic()).unwrap(),
        &Environment::earth_lab(),Vec::new(),0.).unwrap()
}
fn base<'h,'m,'c>(cx:&Cx<'_>,history:&'h ThermalBoundaryHistory<'m>,
    material:&'h ConductivityModel,curve:&'c EquilibriumEnthalpyPhaseCurve,max_steps:usize)
    ->EnthalpyHistoryRun<'h,'m,'c>
{
    EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:history,material,
        element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(2.))},
        EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial {curve,reference_density_kg_m3:1.}],
            of_element:&[0,0],phase_conductivity:None },&[150.;8],
        EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:2},max_steps}).unwrap()
}

#[test]
fn g1_motion_limited_scan_retains_original_endpoints_and_step_budgets() {
    let mesh=mesh();let curve=curve();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=history(cx,&mesh);
        let source=heater(TimeSignal::Ramp {t_start:0.,t_end:1.,
            from:QtyAny::new(0.1,LENGTH),to:QtyAny::new(0.5,LENGTH)},constant(0.1,COMPONENT_POWER_DIMS));
        let make=||ScanningEnthalpyRun::new(cx,base(cx,&history,&material,&curve,20),vec![source.clone()],1,None).unwrap();
        let mut limited=make();let mut manual=make();
        let policy=ScanMotionLimit {max_normalized_travel:0.125,max_chebyshev_coefficients:0};
        limited.set_motion_limit(cx,policy).unwrap();
        let before=limited.conduction().specific_enthalpy_j_kg().to_vec();
        assert!(limited.advance_to(cx,1.,config()).is_err());
        assert_eq!(limited.conduction().specific_enthalpy_j_kg(),before.as_slice());
        assert_eq!(limited.conduction().accepted_steps(),0);
        let mut energy=0.;
        while limited.conduction().time_s()<0.5 {
            let step=limited.advance_event_step(cx,0.5,1.,config()).unwrap();
            assert!(step.end_time_s-step.start_time_s<=0.125);
            let reference=manual.advance_to(cx,step.end_time_s,config()).unwrap();
            assert_eq!(step.report.solution().specific_enthalpy_j_kg,reference.solution().specific_enthalpy_j_kg);
            energy+=step.report.solution().stored_energy_change_j;
        }
        assert!((energy-0.5*(0.1+2.*mesh.total_volume())).abs()<1e-8);
        assert!(limited.set_motion_limit(cx,policy).is_err());
        assert_eq!(limited.motion_limit(),Some(policy));
        assert_eq!(limited.conduction().accepted_steps(),manual.conduction().accepted_steps());
    });
}

#[test]
fn g4_adaptive_scan_applies_the_tightest_motion_limit_and_resumes_without_partial_state() {
    let mesh=mesh();let curve=curve();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let history=with_cx(|cx|history(cx,&mesh));
    let path=|speed|TimeSignal::Ramp {t_start:0.,t_end:1.,from:QtyAny::new(0.1,LENGTH),to:QtyAny::new(0.1+speed,LENGTH)};
    let sources=vec![heater(path(0.2),constant(0.1,COMPONENT_POWER_DIMS)),
        heater(path(0.4),constant(0.1,COMPONENT_POWER_DIMS))];
    let mut run=with_cx(|cx|ScanningEnthalpyRun::new(cx,base(cx,&history,&material,&curve,2),sources,2,None).unwrap());
    let policy=ScanMotionLimit {max_normalized_travel:0.125,max_chebyshev_coefficients:0};
    with_cancelled_cx(|cx|assert!(run.set_motion_limit(cx,policy).is_err()));
    assert_eq!(run.motion_limit(),None);
    with_cx(|cx|{run.set_motion_limit(cx,policy).unwrap();});
    let adaptive=AdaptiveHistoryConfig {absolute_temperature_k:1.,relative_temperature:0.,
        min_step_s:1e-6,max_step_s:1.,max_attempts:8};
    let phase=PhaseErrorTolerance {absolute_enthalpy_j_kg:100.,relative_enthalpy:0.,absolute_liquid_fraction:1.};
    let initial=run.conduction().specific_enthalpy_j_kg().to_vec();
    with_cancelled_cx(|cx|assert!(run.advance_adaptive_step(cx,0.5,1.,config(),adaptive,phase).is_err()));
    assert_eq!(run.conduction().specific_enthalpy_j_kg(),initial.as_slice());
    with_cx(|cx| {
        let pair=run.advance_adaptive_step(cx,0.5,1.,config(),adaptive,phase).unwrap();
        assert!((pair.end_time_s-0.125).abs()<1e-13);
        assert_eq!(pair.rejected_attempts,0);
        assert_eq!(run.conduction().accepted_steps(),2);
        assert!(run.advance_event_step(cx,0.5,1.,config()).is_err());
        run.extend_step_budget(3).unwrap();
        let next=run.advance_event_step(cx,0.5,1.,config()).unwrap();
        assert!(next.end_time_s-next.start_time_s<=0.125);
        assert_eq!(run.motion_limit(),Some(policy));
    });
}

#[test]
fn g1_scan_position_and_power_events_drive_the_original_phase_solver() {
    let mesh=mesh();let curve=curve();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=history(cx,&mesh);
        let source=heater(TimeSignal::Table {times:vec![0.,0.2],values:vec![0.1,1.1],dims:LENGTH,interp:Interp::Hold},
            TimeSignal::Table {times:vec![0.,0.1,0.3],values:vec![0.,0.5,0.],dims:COMPONENT_POWER_DIMS,interp:Interp::Hold});
        let manual=base(cx,&history,&material,&curve,4);
        let mut run=ScanningEnthalpyRun::new(cx,base(cx,&history,&material,&curve,4),vec![source.clone()],1,None).unwrap();
        let mut old=vec![150.;8];let mut previous=0.;let mut stored=0.;
        for end in [0.1,0.2,0.3,0.4] {
            let boundary=history.at(cx,end).unwrap();
            let sampled=source.before(cx,&mesh,end).unwrap();
            let field=ScalarField::Nodal((0..8).map(|i|2.+sampled.source.at(i)).collect());
            let expected=manual.storage().advance(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,
                material:&material,source:&field,element_materials:None},None,&old,end-previous,config()).unwrap();
            let actual=run.advance_event_step(cx,0.4,1.,config()).unwrap();
            assert_eq!(actual.end_time_s,end);
            assert_eq!(actual.report.solution().specific_enthalpy_j_kg,expected.specific_enthalpy_j_kg);
            assert_eq!(actual.report.deposition[0].center_m,sampled.center_m);
            assert!(actual.report.solution().temperature.iter().all(|&t|t==310.));
            stored+=actual.report.solution().stored_energy_change_j;
            old=expected.specific_enthalpy_j_kg;previous=end;
        }
        assert!((stored-(0.5*0.2+2.*mesh.total_volume()*0.4)).abs()<1e-8);
        assert!(run.conduction().liquid_mass_fraction().iter().any(|&f|f>0.5));
        assert_eq!(run.conduction().accepted_steps(),4);
    });
}

#[test]
fn g1_independent_sources_superpose_and_adaptation_preserves_both_half_step_reports() {
    let mesh=mesh();let curve=curve();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=history(cx,&mesh);
        let sources=vec![heater(constant(0.1,LENGTH),TimeSignal::Ramp {t_start:0.,t_end:1.,
            from:QtyAny::new(0.5,COMPONENT_POWER_DIMS),to:QtyAny::new(1.5,COMPONENT_POWER_DIMS)}),
            heater(constant(1.1,LENGTH),constant(0.25,COMPONENT_POWER_DIMS))];
        let mut run=ScanningEnthalpyRun::new(cx,base(cx,&history,&material,&curve,2),sources,2,None).unwrap();
        let pair=run.advance_adaptive_step(cx,0.1,0.1,config(),AdaptiveHistoryConfig {
            absolute_temperature_k:0.001,relative_temperature:0.,min_step_s:1e-6,max_step_s:0.1,max_attempts:20},
            PhaseErrorTolerance {absolute_enthalpy_j_kg:0.01,relative_enthalpy:0.,absolute_liquid_fraction:0.0001}).unwrap();
        assert!(pair.normalized_error<=1. && pair.rejected_attempts>0);
        assert_eq!(pair.first_half.deposition.len(),2);
        let expected=(pair.midpoint_time_s-pair.start_time_s)*pair.first_half.solution().source_w
            +(pair.end_time_s-pair.midpoint_time_s)*pair.second_half.solution().source_w;
        let actual=pair.first_half.solution().stored_energy_change_j+pair.second_half.solution().stored_energy_change_j;
        assert!((actual-expected).abs()<1e-8);
        assert_eq!(run.conduction().accepted_steps(),2);
        assert!(run.advance_to(cx,0.2,config()).is_err());
        run.extend_step_budget(3).unwrap();
        run.advance_to(cx,0.2,config()).unwrap();
    });
}

#[test]
fn g4_failed_source_and_cancelled_steps_leave_all_caloric_fields_unchanged() {
    let mesh=mesh();let curve=curve();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let history=with_cx(|cx|history(cx,&mesh));
    let source=heater(TimeSignal::Table {times:vec![0.,1.],values:vec![0.1,100.],dims:LENGTH,interp:Interp::Hold},
        constant(0.5,COMPONENT_POWER_DIMS));
    let mut run=with_cx(|cx|ScanningEnthalpyRun::new(cx,base(cx,&history,&material,&curve,3),vec![source.clone()],1,None).unwrap());
    with_cx(|cx|{run.advance_to(cx,0.1,config()).unwrap();});
    let h=run.conduction().specific_enthalpy_j_kg().to_vec();
    let f=run.conduction().liquid_mass_fraction().to_vec();
    with_cancelled_cx(|cx|assert!(run.advance_to(cx,0.2,config()).is_err()));
    with_cx(|cx| {
        assert!(run.advance_to(cx,1.,config()).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(),h.as_slice());
        assert_eq!(run.conduction().liquid_mass_fraction(),f.as_slice());
        assert_eq!(run.conduction().time_s(),0.1);
        assert_eq!(run.conduction().accepted_steps(),1);
        assert!(ScanningEnthalpyRun::new(cx,base(cx,&history,&material,&curve,3),vec![source],0,None).is_err());
        run.advance_to(cx,0.2,config()).unwrap();
    });
}

#[test]
fn g1_scanning_keeps_optional_radiation_and_its_physical_heat_balance() {
    use fs_conduction::{SurfaceEmissivity,AmbientRadiationConfig,EMISSIVITY_DIMS,SURFACE_EMISSIVITY_PROPERTY};
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet,PropertyClaim,PropertyKey,PropertyValue,UncertaintyModel,InterpolationPolicy,
        Provenance,MaterialCard,MaterialStateId,SelectionPolicy};
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar {value:0.8,dims:EMISSIVITY_DIMS},validity:ValidityDomain::unconstrained().with("T",280.,380.),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:Vec::new(),provenance:Provenance {source:"synthetic scan fixture".into(),license:"internal-test-use".into(),artifact:None}}).unwrap();
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"gray fixture".into(),phase:"solid".into(),
        process:"synthetic".into(),revision:0},claims,Vec::new()).unwrap();
    let emissivity=SurfaceEmissivity::from_card("surface",&card,310.,SelectionPolicy::SingleClaimOnly).unwrap();
    let mesh=mesh();let curve=curve();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=history(cx,&mesh);
        let radiation=ScanRadiation {patches:vec![RadiationRegionHistory {region:"surface".into(),emissivity,
            ambient_temperature:constant(300.,TEMPERATURE_DIMS)}],config:AmbientRadiationConfig::default()};
        let mut run=ScanningEnthalpyRun::new(cx,base(cx,&history,&material,&curve,2),
            vec![heater(constant(0.1,LENGTH),constant(0.5,COMPONENT_POWER_DIMS))],1,Some(radiation)).unwrap();
        let step=run.advance_event_step(cx,0.001,0.001,config()).unwrap();
        let ScanTransport::Radiation(report)=step.report.transport else {panic!("radiation was omitted")};
        assert_eq!(report.convective_out_w,0.);
        assert!(report.radiation.nonlinear_radiation_out_w>0.);
        assert!(report.physical_energy_residual_j.abs()<=config().energy_tolerance_j);
        assert!(report.physical_residual_norm_j<=report.physical_residual_tolerance_j);
        assert!((report.conduction.source_w-(0.5+2.*mesh.total_volume())).abs()<1e-10);
    });
}

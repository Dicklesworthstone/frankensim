//! Original material-contact solves under schedules, radiation and time adaptation.
mod support;

use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun};
use fs_conduction::bc::history::march::enthalpy::contact::ContactEnthalpyHistoryRun;
use fs_conduction::bc::history::march::enthalpy::adaptive::PhaseErrorTolerance;
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContact, TemperatureContacts};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::{HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial};
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionMesh, ConductionProblem,
    ConductivityModel, ConductivityTable, ScalarField, SurfaceEmissivity, ThermalBc,
    ThermalBoundary, ThermalBoundaryBuilder, ThermalInterfaces, TEMPERATURE_DIMS, EMISSIVITY_DIMS,
    SURFACE_EMISSIVITY_PROPERTY, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const INITIAL: [f64; 8] = [80.,80.,80.,80.,450.,450.,450.,450.];
fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(8,vec![[0,1,2,3],[4,5,6,7]]),
        vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.],
            [0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).unwrap()
}
fn phase(hot: bool) -> EquilibriumEnthalpyPhaseCurve {
    let knots=if hot { [(-200.,280.,0.),(500.,350.,0.),(600.,350.,1.),(900.,380.,1.)] }
        else { [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)] };
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([if hot {21} else {20};32]),
        knots.into_iter().map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h,
            temperature_k:t,liquid_mass_fraction:f,bulk_density_kg_m3:1. }).collect()).unwrap()
}
fn state(name: &str) -> MaterialStateId {
    MaterialStateId { chemistry:name.into(),phase:"solid".into(),process:"synthetic test".into(),revision:0 }
}
fn claim(name: &str, dims: Dims, value: PropertyValue, interpolation: InterpolationPolicy) -> PropertyClaim {
    PropertyClaim { key:PropertyKey::new(name,dims),value,
        validity:ValidityDomain::unconstrained().with("T",280.,380.),
        uncertainty:UncertaintyModel::Unstated, interpolation, observations:Vec::new(),
        provenance:Provenance {source:"synthetic contact history fixture".into(),license:"internal-test-use".into(),artifact:None} }
}
fn card() -> InterfaceSystemCard {
    let dims=AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS; let mut claims=ClaimSet::new();
    claims.insert_claim(claim(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,dims,
        PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:TEMPERATURE_DIMS,
            knots:vec![(280.,0.1),(380.,0.7)],dims},InterpolationPolicy::LinearInside)).unwrap();
    InterfaceSystemCard::assemble(SurfaceSpec {material:state("A"),texture_frame:"reference".into()},
        SurfaceSpec {material:state("B"),texture_frame:"reference".into()},
        SystemContext {medium:"dry".into(),third_body:None,environment:"vacuum".into(),history:"fixed".into()},
        claims,Vec::new()).unwrap()
}
fn patch(signal: TimeSignal) -> RadiationRegionHistory {
    let mut claims=ClaimSet::new(); claims.insert_claim(claim(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS,
        PropertyValue::Scalar {value:0.2,dims:EMISSIVITY_DIMS},InterpolationPolicy::ConstantWithinValidity)).unwrap();
    let card=MaterialCard::assemble(state("surface"),claims,Vec::new()).unwrap();
    RadiationRegionHistory {region:"exterior".into(),ambient_temperature:signal,
        emissivity:SurfaceEmissivity::from_card("exterior",&card,320.,SelectionPolicy::SingleClaimOnly).unwrap()}
}
fn boundary<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(mesh)
        .region("exterior",|f| f.centroid[2]<0. && f.centroid[1]==0.,ThermalBc::adiabatic()).unwrap()
        .adiabatic_remainder(),&Environment::earth_lab(),Vec::new(),0.).unwrap()
}
fn contacts<'m,'i>(cx: &Cx<'_>, mesh: &'m ConductionMesh, boundary: &ThermalBoundary,
    card: &'i InterfaceSystemCard) -> TemperatureContacts<'m,'i>
{
    let law=TemperatureContact::new("joint",ThermalInterfaces::coincident_face_pairs(mesh).unwrap(),card,
        &QueryPoint::new().with("T",320.).unwrap(),"T",SelectionPolicy::SingleClaimOnly).unwrap();
    TemperatureContacts::new(cx,mesh,boundary,Vec::new(),vec![law],1).unwrap()
}
fn source() -> HeatSourceHistory {
    HeatSourceHistory::Scaled {field:ScalarField::Nodal(vec![0.,0.,0.,0.,60.,60.,60.,60.]),
        scale:TimeSignal::Ramp {t_start:0.,t_end:0.04,from:QtyAny::new(1.,Dims([0;6])),to:QtyAny::new(2.,Dims([0;6]))}}
}
fn source_at(t: f64) -> ScalarField {
    let q=60.*(1.+(t/0.04).clamp(0.,1.)); ScalarField::Nodal(vec![0.,0.,0.,0.,q,q,q,q])
}
fn table<'c>(cold: &'c EquilibriumEnthalpyPhaseCurve, hot: &'c EquilibriumEnthalpyPhaseCurve)
    -> [ReferenceEnthalpyMaterial<'c>;2]
{
    [ReferenceEnthalpyMaterial {curve:cold,reference_density_kg_m3:1.},
        ReferenceEnthalpyMaterial {curve:hot,reference_density_kg_m3:2.}]
}
fn base<'h,'m,'c>(cx:&Cx<'_>, history:&'h ThermalBoundaryHistory<'m>, material:&'h ConductivityModel,
    cold:&'c EquilibriumEnthalpyPhaseCurve,hot:&'c EquilibriumEnthalpyPhaseCurve,
    interfaces:Option<&'h ThermalInterfaces>,max_steps:usize) -> EnthalpyHistoryRun<'h,'m,'c>
{
    EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:history,material,element_materials:None,
        interfaces,source:source()},EnthalpyHistoryMaterials {table:&table(cold,hot),of_element:&[0,1],phase_conductivity:None},
        &INITIAL,EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:2},max_steps}).unwrap()
}
fn controls() -> EnthalpyStepConfig {
    EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
        linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
        max_newton_iterations:64,energy_tolerance_j:1e-8}
}
fn coupling() -> ContactIterationConfig {
    ContactIterationConfig {max_iterations:64,relaxation:0.5,residual_tolerance_j:1e-8}
}
fn material() -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.6),(380.,3.6)]).unwrap())
}

#[test]
fn g1_contact_history_matches_original_scheduled_endpoints_with_and_without_radiation() {
    let mesh=mesh();let cold=phase(false);let hot=phase(true);let material=material();let card=card();
    with_cx(|cx| {
        let history=boundary(cx,&mesh);let fixed=history.at(cx,0.).unwrap();
        let manual=contacts(cx,&mesh,&fixed,&card);
        let storage=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&table(&cold,&hot),&[0,1],
            EnthalpyBudget {max_vertices:8,max_elements:2}).unwrap();
        let scheduled=patch(TimeSignal::Table {times:vec![0.,0.02],values:vec![300.,350.],dims:TEMPERATURE_DIMS,interp:Interp::Hold});
        for radiation in [false,true] {
            let mut run=ContactEnthalpyHistoryRun::new(cx,base(cx,&history,&material,&cold,&hot,None,8),
                contacts(cx,&mesh,&fixed,&card),coupling()).unwrap();
            if radiation { run=run.with_radiation(cx,vec![scheduled.clone()],AmbientRadiationConfig::default()).unwrap(); }
            let mut old=INITIAL.to_vec();let mut previous=0.;
            for time in [0.01,0.02,0.04] {
                let source=source_at(time);
                let problem=ConductionProblem {mesh:&mesh,boundary:&fixed,material:&material,source:&source,element_materials:None};
                let expected=if radiation {
                    let reservoir=if time<=0.02 {300.} else {350.};
                    let sampled=AmbientRadiationPatch::new("exterior",scheduled.emissivity.clone(),reservoir).unwrap();
                    manual.advance_enthalpy_with_ambient_radiation(cx,&storage,problem,&old,time-previous,
                        controls(),coupling(),&[sampled],AmbientRadiationConfig::default()).unwrap()
                } else { manual.advance_enthalpy(cx,&storage,problem,&old,time-previous,controls(),coupling()).unwrap() };
                let actual=run.advance_event_step(cx,time,1.,controls()).unwrap();
                assert_eq!(actual.end_time_s,time);
                assert_eq!(actual.report.conduction.specific_enthalpy_j_kg,expected.conduction.specific_enthalpy_j_kg);
                assert_eq!(actual.report.physical_contacts.states(),expected.physical_contacts.states());
                assert_eq!(actual.report.physical.energy_residual_j,expected.physical.energy_residual_j);
                assert!(actual.report.physical.residual_norm_j<=coupling().residual_tolerance_j);
                assert!(actual.report.physical.energy_residual_j.abs()<=controls().energy_tolerance_j);
                assert_eq!(actual.report.radiation.is_some(),radiation);
                let local=&actual.report.physical_contacts.states()[0];
                assert!((local.resistance.value_m2_k_per_w()-(0.1+0.006*(local.interface_temperature_k-280.))).abs()<1e-12);
                old=expected.conduction.specific_enthalpy_j_kg;previous=time;
            }
            assert_eq!(run.conduction().specific_enthalpy_j_kg(),old.as_slice());
            assert_eq!(run.conduction().accepted_steps(),3);
        }
    });
}

#[test]
fn g1_adaptive_contact_radiation_keeps_the_actual_pair_and_no_partial_rejected_state() {
    let mesh=mesh();let cold=phase(false);let hot=phase(true);let material=material();let card=card();
    let history=with_cx(|cx| boundary(cx,&mesh));
    let scheduled=patch(TimeSignal::Constant(QtyAny::new(300.,TEMPERATURE_DIMS)));
    let mut run=with_cx(|cx| {
        let fixed=history.at(cx,0.).unwrap();
        ContactEnthalpyHistoryRun::new(cx,base(cx,&history,&material,&cold,&hot,None,2),
            contacts(cx,&mesh,&fixed,&card),coupling()).unwrap()
            .with_radiation(cx,vec![scheduled.clone()],AmbientRadiationConfig::default()).unwrap()
    });
    let adaptive=AdaptiveHistoryConfig {absolute_temperature_k:0.002,relative_temperature:0.,min_step_s:1e-7,max_step_s:0.04,max_attempts:24};
    let phase=PhaseErrorTolerance {absolute_enthalpy_j_kg:0.02,relative_enthalpy:0.,absolute_liquid_fraction:0.001};
    with_cancelled_cx(|cx| assert!(run.advance_adaptive_step(cx,0.04,0.04,controls(),adaptive,phase).is_err()));
    with_cx(|cx| {
        let strict=AdaptiveHistoryConfig {absolute_temperature_k:1e-12,max_attempts:1,..adaptive};
        assert!(run.advance_adaptive_step(cx,0.04,0.04,controls(),strict,phase).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(),&INITIAL);
        assert_eq!(run.conduction().time_s(),0.);assert_eq!(run.conduction().accepted_steps(),0);
        let selected=run.advance_adaptive_step(cx,0.04,0.04,controls(),adaptive,phase).unwrap();
        assert!(selected.normalized_error<=1.);assert_eq!(run.conduction().accepted_steps(),2);
        let fixed=history.at(cx,0.).unwrap();let mut old=INITIAL.to_vec();let mut from=0.;
        for (to,actual) in [(selected.midpoint_time_s,&selected.first_half),(selected.end_time_s,&selected.second_half)] {
            let source=source_at(to);let patch=scheduled.at(to).unwrap();
            let expected=run.contacts().advance_enthalpy_with_ambient_radiation(cx,run.conduction().storage(),
                ConductionProblem {mesh:&mesh,boundary:&fixed,material:&material,element_materials:None,source:&source},
                &old,to-from,controls(),coupling(),&[patch],AmbientRadiationConfig::default()).unwrap();
            assert_eq!(actual.conduction.specific_enthalpy_j_kg,expected.conduction.specific_enthalpy_j_kg);
            assert_eq!(actual.physical_contacts.states(),expected.physical_contacts.states());
            assert!(actual.physical.energy_residual_j.abs()<=controls().energy_tolerance_j);
            old=expected.conduction.specific_enthalpy_j_kg;from=to;
        }
        assert_eq!(run.conduction().specific_enthalpy_j_kg(),old.as_slice());
        assert!(run.advance_to(cx,0.04,controls()).is_err());
        run.extend_step_budget(3).unwrap();
        run.advance_to(cx,0.04,controls()).unwrap();
        assert_eq!(run.conduction().accepted_steps(),3);
        assert!(run.radiation().is_some());
    });
}

#[test]
fn g0_contact_admission_and_invalid_future_reservoir_preserve_original_history() {
    let mesh=mesh();let other_mesh=self::mesh();let cold=phase(false);let hot=phase(true);let material=material();let card=card();
    with_cx(|cx| {
        let history=boundary(cx,&mesh);let fixed=history.at(cx,0.).unwrap();
        let make=|| base(cx,&history,&material,&cold,&hot,None,4);
        let other=boundary(cx,&other_mesh).at(cx,0.).unwrap();
        assert!(ContactEnthalpyHistoryRun::new(cx,make(),contacts(cx,&other_mesh,&other,&card),coupling()).is_err());
        let contacts=contacts(cx,&mesh,&fixed,&card);
        let snapshot=contacts.evaluate(cx,&fixed,make().temperature()).unwrap();
        let duplicate=base(cx,&history,&material,&cold,&hot,Some(snapshot.frozen_interfaces()),4);
        assert!(ContactEnthalpyHistoryRun::new(cx,duplicate,contacts,coupling()).is_err());
        let mut run=ContactEnthalpyHistoryRun::new(cx,make(),self::contacts(cx,&mesh,&fixed,&card),coupling()).unwrap()
            .with_radiation(cx,vec![patch(TimeSignal::Table {times:vec![0.,0.03],values:vec![300.,-1.],
                dims:TEMPERATURE_DIMS,interp:Interp::Hold})],AmbientRadiationConfig::default()).unwrap();
        run.advance_to(cx,0.01,controls()).unwrap();
        let old=run.conduction().specific_enthalpy_j_kg().to_vec();
        let temperature=run.conduction().temperature().to_vec();let fraction=run.conduction().liquid_mass_fraction().to_vec();
        assert!(run.advance_to(cx,0.03,controls()).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(),old.as_slice());
        assert_eq!(run.conduction().temperature(),temperature.as_slice());
        assert_eq!(run.conduction().liquid_mass_fraction(),fraction.as_slice());
        assert_eq!(run.conduction().time_s(),0.01);assert_eq!(run.conduction().accepted_steps(),1);
        run.advance_to(cx,0.02,controls()).unwrap();
    });
}

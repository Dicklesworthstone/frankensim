//! Latent error coordinates, radiative material retention and accepted-pair work.
mod support;

use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem, ThermalHistoryRun};
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun};
use fs_conduction::bc::history::march::enthalpy::adaptive::PhaseErrorTolerance;
use fs_conduction::bc::history::march::enthalpy::radiation::RadiativeEnthalpyHistoryRun;
use fs_conduction::bc::history::march::radiation::{RadiationRegionHistory, RadiativeThermalHistoryRun};
use fs_conduction::transient::{VolumetricHeatCapacity, backward_euler::{BackwardEuler, StepConfig}};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionMesh, ConductionProblem,
    ConductivityModel, LinearConfig, ScalarField, SurfaceEmissivity, ThermalBc, ThermalBoundaryBuilder,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0,1,2,3]]),
        vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).unwrap()
}
fn boundary<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .remainder("surface", ThermalBc::adiabatic()).unwrap(), &Environment::earth_lab(), Vec::new(), 0.).unwrap()
}
fn problem<'h,'m>(history: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel,
    source: HeatSourceHistory) -> ThermalHistoryProblem<'h,'m>
{
    ThermalHistoryProblem { boundary: history, material, source, element_materials: None, interfaces: None }
}
fn curve(offset: f64) -> EquilibriumEnthalpyPhaseCurve {
    let mut identity=[74;32]; identity[..8].copy_from_slice(&offset.to_le_bytes());
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash(identity),
        [(-200.,280.,0.),(100.,310.,0.),(500.,310.,1.),(1200.,380.,1.)].into_iter()
        .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h+offset, temperature_k:t,
            liquid_mass_fraction:f, bulk_density_kg_m3:1. }).collect()).unwrap()
}
fn ramp(from: f64, to: f64, end: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Ramp { t_start:0., t_end:end, from:QtyAny::new(from,dims), to:QtyAny::new(to,dims) }
}
fn controls() -> EnthalpyStepConfig {
    EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11, relative_tolerance:1e-12,
        linear_restart:4, max_linear_cycles:16, forcing_minimum:1e-12, forcing_maximum:1e-3,
        ..Default::default() }, max_newton_iterations:64, energy_tolerance_j:1e-8 }
}
fn adaptive() -> AdaptiveHistoryConfig {
    AdaptiveHistoryConfig { absolute_temperature_k:1e-5, relative_temperature:0., min_step_s:1e-7,
        max_step_s:1., max_attempts:24 }
}
fn phase() -> PhaseErrorTolerance {
    PhaseErrorTolerance { absolute_enthalpy_j_kg:1e9, relative_enthalpy:0., absolute_liquid_fraction:0.001 }
}
fn run<'h,'m,'c>(cx:&Cx<'_>, history:&'h ThermalBoundaryHistory<'m>, material:&'h ConductivityModel,
    curve:&'c EquilibriumEnthalpyPhaseCurve, offset:f64, source:HeatSourceHistory)
    -> EnthalpyHistoryRun<'h,'m,'c>
{
    EnthalpyHistoryRun::new(cx, problem(history,material,source), EnthalpyHistoryMaterials {
        table:&[ReferenceEnthalpyMaterial { curve, reference_density_kg_m3:1. }], of_element:&[0],
        phase_conductivity:None }, &[150.+offset;4], EnthalpyHistoryBudget {
            spatial:EnthalpyBudget { max_vertices:4,max_elements:1 }, max_steps:4 }).unwrap()
}
fn emissivity() -> SurfaceEmissivity {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim { key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar { value:0.8,dims:EMISSIVITY_DIMS },
        validity:ValidityDomain::unconstrained().with("T",280.,380.), uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity, observations:Vec::new(),
        provenance:Provenance { source:"synthetic adaptive radiation fixture".into(), license:"internal-test-use".into(),artifact:None }}).unwrap();
    let card=MaterialCard::assemble(MaterialStateId { chemistry:"gray fixture".into(),phase:"solid".into(),
        process:"synthetic".into(),revision:0 },claims,Vec::new()).unwrap();
    SurfaceEmissivity::from_card("surface",&card,320.,SelectionPolicy::SingleClaimOnly).unwrap()
}

#[test]
fn g1_phase_error_rejects_an_inaccurate_isothermal_plateau_and_ignores_enthalpy_origin() {
    let mesh=mesh(); let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=boundary(cx,&mesh); let mut accepted_ends=Vec::new();
        for offset in [0.,10000.] {
            let chart=curve(offset);
            let mut model=run(cx,&history,&material,&chart,offset,
                HeatSourceHistory::Uniform(ramp(0.,40.,1.,VOLUMETRIC_SOURCE_DIMS)));
            let solved=model.advance_adaptive_step(cx,1.,1.,controls(),adaptive(),phase()).unwrap();
            assert!(solved.rejected_attempts>0);
            assert!(solved.normalized_error<=1.);
            assert!(model.temperature().iter().all(|&t| t==310.));
            // h'=40*t. The accepted equal-half backward-Euler path is h0+30*dt^2.
            let expected=150.+offset+30.*solved.end_time_s*solved.end_time_s;
            assert!(model.specific_enthalpy_j_kg().iter().all(|&h| (h-expected).abs()<1e-8));
            assert!(model.liquid_mass_fraction().iter().all(|&f| (f-(expected-offset-100.)/400.).abs()<1e-10));
            assert_eq!(model.accepted_steps(),2);
            accepted_ends.push(solved.end_time_s);
            for half in [&solved.first_half,&solved.second_half] { assert!(half.energy_residual_j.abs()<=controls().energy_tolerance_j); }
        }
        assert_eq!(accepted_ends[0],accepted_ends[1]);
    });
}

#[test]
fn g0_rejected_phase_trials_and_cancellation_publish_no_partial_half_step() {
    let mesh=mesh();let material=ConductivityModel::isotropic_declared(2.).unwrap();let chart=curve(0.);
    let history=with_cx(|cx| boundary(cx,&mesh));
    let mut model=with_cx(|cx| run(cx,&history,&material,&chart,0.,
        HeatSourceHistory::Uniform(ramp(0.,40.,1.,VOLUMETRIC_SOURCE_DIMS))));
    with_cancelled_cx(|cx| assert!(model.advance_adaptive_step(cx,1.,1.,controls(),adaptive(),phase()).is_err()));
    with_cx(|cx| {
        let limited=AdaptiveHistoryConfig { max_attempts:1,..adaptive() };
        assert!(model.advance_adaptive_step(cx,1.,1.,controls(),limited,phase()).is_err());
        assert_eq!(model.specific_enthalpy_j_kg(),&[150.;4]);
        assert_eq!(model.liquid_mass_fraction(),&[0.125;4]);
        assert_eq!(model.temperature(),&[310.;4]);
        assert_eq!(model.time_s(),0.); assert_eq!(model.accepted_steps(),0);
        // Independently tighten enthalpy rather than the fraction allowance.
        let enthalpy_only=PhaseErrorTolerance { absolute_enthalpy_j_kg:0.05,absolute_liquid_fraction:1.,..phase() };
        let accepted=model.advance_adaptive_step(cx,1.,1.,controls(),adaptive(),enthalpy_only).unwrap();
        assert!(accepted.rejected_attempts>0); assert!(accepted.normalized_error<=1.);
    });
}

#[test]
fn g1_radiative_phase_adaptation_retains_latent_heat_and_reservoir_time() {
    let mesh=mesh(); let material=ConductivityModel::isotropic_declared(2.).unwrap(); let chart=curve(0.);
    with_cx(|cx| {
        let history=boundary(cx,&mesh);
        let patch=RadiationRegionHistory { region:"surface".into(),emissivity:emissivity(),
            ambient_temperature:ramp(310.,340.,0.02,TEMPERATURE_DIMS) };
        let mut model=RadiativeEnthalpyHistoryRun::new(cx,run(cx,&history,&material,&chart,0.,
            HeatSourceHistory::Fixed(ScalarField::Uniform(0.))),vec![patch.clone()],AmbientRadiationConfig::default()).unwrap();
        let step=model.advance_adaptive_step(cx,0.02,0.02,controls(),adaptive(),phase()).unwrap();
        assert!(step.rejected_attempts>0); assert!(step.normalized_error<=1.);
        assert!(model.conduction().temperature().iter().all(|&t| t==310.));
        let flux=|time:f64| 0.8*fs_conduction::STEFAN_BOLTZMANN_W_M2_K4*((310.+1500.*time).powi(4)-310.0_f64.powi(4));
        let supplied=flux(step.midpoint_time_s)*step.midpoint_time_s
            +flux(step.end_time_s)*(step.end_time_s-step.midpoint_time_s);
        for (node,&h) in model.conduction().specific_enthalpy_j_kg().iter().enumerate() {
            let weight=if node==0 { 0.5 } else { (1.+3.0_f64.sqrt()/2.)/3. };
            assert!((h-(150.+weight*supplied/(1./24.))).abs()<1e-7);
        }
        for half in [&step.first_half,&step.second_half] {
            assert!(half.physical_energy_residual_j.abs()<=controls().energy_tolerance_j);
            assert_eq!(half.convective_out_w,0.); assert!(half.convective_robin_fluxes.is_empty());
            assert_eq!(half.radiation.patches[0].patch.emissivity(),&patch.emissivity);
        }
        assert_eq!(model.conduction().accepted_steps(),2);
        assert_eq!(model.patches(),&[patch]);
    });
}

#[test]
fn g1_radiative_temperature_adaptation_matches_the_actual_two_endpoint_solves() {
    let mesh=mesh();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=boundary(cx,&mesh);let zero=ScalarField::Uniform(0.);
        let base=ThermalHistoryRun::new(cx,problem(&history,&material,HeatSourceHistory::Fixed(zero.clone())),
            VolumetricHeatCapacity::declared(15.).unwrap(),&[350.;4],2).unwrap();
        let patch=RadiationRegionHistory { region:"surface".into(),emissivity:emissivity(),
            ambient_temperature:TimeSignal::Table { times:vec![0.,0.004],values:vec![300.,340.],
                dims:TEMPERATURE_DIMS,interp:Interp::Hold }};
        let mut model=RadiativeThermalHistoryRun::new(cx,base,vec![patch.clone()],AmbientRadiationConfig::default()).unwrap();
        let config=StepConfig { linear:LinearConfig { tolerance:1e-12,..LinearConfig::default() },energy_tolerance_j:1e-7 };
        let selected=model.advance_adaptive_step(cx,0.01,0.01,config,adaptive()).unwrap();
        assert!(selected.end_time_s<=0.004);assert!(selected.normalized_error<=1.);
        let engine=BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(15.).unwrap()).unwrap();
        let boundary=history.at(cx,0.).unwrap();
        let original=ConductionProblem { mesh:&mesh,boundary:&boundary,material:&material,source:&zero,element_materials:None };
        let fixed=AmbientRadiationPatch::new("surface",patch.emissivity.clone(),300.).unwrap();
        let first=engine.advance_prescribed_with_ambient_radiation(cx,original,None,&[350.;4],
            selected.midpoint_time_s,config,None,&[fixed.clone()],AmbientRadiationConfig::default()).unwrap();
        let second=engine.advance_prescribed_with_ambient_radiation(cx,original,None,&first.conduction.temperature,
            selected.end_time_s-selected.midpoint_time_s,config,None,&[fixed],AmbientRadiationConfig::default()).unwrap();
        assert_eq!(selected.first_half.conduction.temperature,first.conduction.temperature);
        assert_eq!(selected.second_half.conduction.temperature,second.conduction.temperature);
        assert_eq!(selected.second_half.radiation,second.radiation);
        assert_eq!(model.conduction().temperature(),second.conduction.temperature.as_slice());
        assert_eq!(model.conduction().accepted_steps(),2);
    });
}

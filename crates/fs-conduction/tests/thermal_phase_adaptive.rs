//! Local enthalpy/fraction accuracy and atomic adaptive phase publication.
mod support;

use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun};
use fs_conduction::bc::history::march::enthalpy::adaptive::PhaseErrorTolerance;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh, ConductivityModel, ThermalBc, ThermalBoundaryBuilder, VOLUMETRIC_SOURCE_DIMS};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::TimeSignal;
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

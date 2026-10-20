//! Liquid fraction is a constitutive observable, not scaled temperature.
mod support;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun, EnthalpyHistoryMaterials, EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityBudget, ScanSensitivityError};
use fs_conduction::bc::history::march::moving::{MovingHeatSource, MovingHeatBudget};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh, ConductivityModel, ScalarField, ThermalBc, ThermalBoundaryBuilder,
    LinearConfig, COMPONENT_POWER_DIMS};
use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::{scenario::Environment, signal::TimeSignal};
use fs_solver::NewtonKrylovConfig;
use support::{with_cx, with_cancelled_cx};

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0,1,2,3]]),
        vec![[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]).unwrap()
}
fn chart(kink: bool) -> EquilibriumEnthalpyPhaseCurve {
    let mut knots = vec![(-200.,280.,0.), (100.,310.,0.)];
    if kink { knots.push((150.,310.,0.25)); }
    knots.extend([(200.,310.,1.), (900.,380.,1.)]);
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([91;32]), knots.into_iter()
        .map(|(h,t,f)| EnthalpyPhaseKnot {specific_enthalpy_j_kg:h, temperature_k:t,
            liquid_mass_fraction:f, bulk_density_kg_m3:1.}).collect()).unwrap()
}
fn source(power: f64) -> MovingHeatSource {
    MovingHeatSource::new([0.2,0.25,0.3].map(|v|
        TimeSignal::Constant(QtyAny::new(v, Dims([1,0,0,0,0,0])))), [1.4,1.3,1.2],
        TimeSignal::Constant(QtyAny::new(power, COMPONENT_POWER_DIMS)), 0.,
        MovingHeatBudget {max_vertices:4, max_elements:1, min_active_vertices:2,
            power_relative_tolerance:1e-12}).unwrap()
}
fn history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .remainder("surface", ThermalBc::adiabatic()).unwrap(),
        &Environment::earth_lab(), vec![], 0.).unwrap()
}
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-12,
        relative_tolerance:1e-12, linear_restart:4, max_linear_cycles:16,
        forcing_minimum:1e-12, forcing_maximum:1e-3, ..Default::default()},
        max_newton_iterations:64, energy_tolerance_j:1e-9}
}
fn linear() -> LinearConfig { LinearConfig {tolerance:1e-11, max_iterations:1000, restart:4} }
fn sensitivity() -> ScanSensitivityBudget {
    ScanSensitivityBudget {max_source_vertices:4, max_radiation_feedback_entries:1024}
}
fn run<'h,'m,'c>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, curve: &'c EquilibriumEnthalpyPhaseCurve,
    old: &[f64], power: f64) -> ScanningEnthalpyRun<'h,'m,'c>
{
    let base = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {boundary:history, material,
        element_materials:None, interfaces:None, source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
        EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve, reference_density_kg_m3:1.}],
            of_element:&[0], phase_conductivity:None}, old,
        EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:4, max_elements:1}, max_steps:16}).unwrap();
    ScanningEnthalpyRun::new(cx, base, vec![source(power)], 1, None).unwrap()
}
fn close(a: f64, b: f64) { assert!((a-b).abs() < 2e-7*(1.+a.abs().max(b.abs())), "{a:e} != {b:e}"); }

#[test]
fn g3_fraction_source_and_history_gradients_match_fresh_physical_differences() {
    let mesh=mesh(); let curve=chart(false); let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=history(cx,&mesh); let seed=[0.3,-0.1,0.7,0.5];
        for old in [[20.;4], [130.,150.,160.,170.], [300.;4]] {
            let mut model=run(cx,&history,&material,&curve,&old,3.);
            let binding=model.advance_linearized_to(cx,0.002,config(),sensitivity()).unwrap();
            let actual=binding.pullback_liquid_mass_fraction(cx,&seed,linear()).unwrap();
            let objective=|old:&[f64],power| {
                let mut candidate=run(cx,&history,&material,&curve,old,power);
                let step=candidate.advance_to(cx,0.002,config()).unwrap();
                seed.iter().zip(&step.solution().liquid_mass_fraction).map(|(s,f)|s*f).sum::<f64>()
            };
            let eps=1e-4;
            close(actual.sources[0].absorbed_power_w,
                (objective(&old,3.+eps)-objective(&old,3.-eps))/(2.*eps));
            for vertex in 0..4 {
                let (mut plus,mut minus)=(old,old); plus[vertex]+=eps; minus[vertex]-=eps;
                close(actual.storage.transport.previous_specific_enthalpy[vertex],
                    (objective(&plus,3.)-objective(&minus,3.))/(2.*eps));
            }
            if old[0]==130. {
                assert!(actual.sources[0].absorbed_power_w>0.);
                let temperature=binding.pullback_temperature(cx,&seed,linear()).unwrap();
                assert_eq!(temperature.sources[0].absorbed_power_w,0.);
            } else { assert_eq!(actual.sources[0].absorbed_power_w,0.); }
        }
    });
}

#[test]
fn g0_fraction_only_corner_refuses_without_disabling_transport_or_unobserved_nodes() {
    let mesh=mesh(); let curve=chart(true); let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let history=history(cx,&mesh);
        let mut model=run(cx,&history,&material,&curve,&[150.,160.,160.,160.],0.);
        let binding=model.advance_linearized_to(cx,0.002,config(),sensitivity()).unwrap();
        binding.pullback(cx,&[1.,0.,0.,0.],linear()).unwrap();
        binding.pullback_temperature(cx,&[1.,0.,0.,0.],linear()).unwrap();
        assert!(matches!(binding.pullback_liquid_mass_fraction(cx,&[1.,0.,0.,0.],linear()),
            Err(ScanSensitivityError::Adjoint(EnthalpyAdjointError::ChartKink {vertex:0,..}))));
        let projected=binding.liquid_mass_fraction_pullback(cx,&[0.,1.,0.,0.]).unwrap();
        close(projected[1],0.015);
        assert_eq!(binding.liquid_mass_fraction_pullback(cx,&[0.;4]).unwrap(),vec![0.;4]);
        assert_eq!(model.conduction().accepted_steps(),1);
        assert_eq!(model.conduction().specific_enthalpy_j_kg(),&[150.,160.,160.,160.]);
    });
}

#[test]
fn g0_fraction_pullback_rejects_bad_seeds_and_honors_cancellation() {
    let mesh=mesh(); let curve=chart(false); let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let history=with_cx(|cx|history(cx,&mesh));
    let binding=with_cx(|cx|run(cx,&history,&material,&curve,&[150.;4],3.)
        .advance_linearized_to(cx,0.002,config(),sensitivity()).unwrap());
    with_cx(|cx| {
        assert!(binding.liquid_mass_fraction_pullback(cx,&[]).is_err());
        assert!(binding.liquid_mass_fraction_pullback(cx,&[f64::NAN;4]).is_err());
        assert!(binding.liquid_mass_fraction_pullback(cx,&[f64::INFINITY;4]).is_err());
        assert_eq!(binding.liquid_mass_fraction_pullback(cx,&[0.;4]).unwrap(),vec![0.;4]);
    });
    with_cancelled_cx(|cx| assert!(binding.pullback_liquid_mass_fraction(cx,&[1.;4],linear()).is_err()));
}

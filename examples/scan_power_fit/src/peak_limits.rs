//! One spatial/time hot-spot constraint, rather than one constraint per node.
use super::*;
use fs_ascent::{SqpRunReport,SqpStop};
use fs_ascent::scan_design::{ScanGeometryControl,constraints::{ConstrainedScanStudy,
    ScanPeakBound,ScanPeakResponse,ScanResponseBudget,ScanConstraintEvaluation}};

struct PeakFit {
    report:SqpRunReport,powers:Vec<f64>,constraints:ScanConstraintEvaluation,
    temperature:Vec<f64>,evaluations:usize,gradient_error:f64,
}

fn fit_peak(cx:&Cx<'_>, maximum:f64, smoothing:f64, audit:bool) -> Result<PeakFit> {
    if !maximum.is_finite() || maximum<=150. || maximum>=190.
        || !smoothing.is_finite() || smoothing<=0. || smoothing>1.
    {return Err("require 150 < maximum < 190 J/kg and 0 < smoothing <= 1 J/kg".into());}
    let (complex,positions)=fs_conduction::fixtures::box_grid([1;3],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let material=ConductivityModel::isotropic_declared(2.)?;
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([123;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),vec![],0.)?;
    let source=MovingHeatSource::new([0.2,0.3,0.4].map(|v|TimeSignal::Constant(
        QtyAny::new(v,Dims([1,0,0,0,0,0])))),[1.8;3],TimeSignal::Table {
            times:vec![0.,1.,2.],values:vec![1.,1.,0.],dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,
        MovingHeatBudget {max_vertices:8,max_elements:6,min_active_vertices:1,power_relative_tolerance:1e-11})?;
    let factory=|cx:&Cx<'_>,sources:Vec<MovingHeatSource>|
        ->std::result::Result<ScanningEnthalpyRun<'_, '_, '_>,ScanSensitivityError> {
        let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&boundary,material:&material,
            element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
            EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}],
                of_element:&[0;6],phase_conductivity:None},&[150.;8],
            EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:6},max_steps:2})?;
        Ok(ScanningEnthalpyRun::new(cx,base,sources,1,None)?)
    };
    let objective=ScanPowerObjective {endpoints_s:vec![1.,2.],observations:vec![
        ScanObservation {endpoint:0,vertex:7,observable:ScanObservable::SpecificEnthalpy,target:152.,scale:10.},
        ScanObservation {endpoint:1,vertex:7,observable:ScanObservable::SpecificEnthalpy,target:157.,scale:10.}],
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
            linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
            max_newton_iterations:64,energy_tolerance_j:1e-8},
        transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
        budget:ScanObjectiveBudget {max_steps:2,max_observations:2,max_retained_vertices:16,
            max_retained_source_vertices:16,max_schedule_controls:16,max_retained_radiation_entries:0},
        source_reference_time_s:0.};
    let geometry=if audit {vec![
        ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:0.2,scale_m:0.25},
        ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:1.8}]} else {vec![]};
    let study=ScanPowerStudy::new_joint(cx,vec![source],vec![
        ScanPowerControl {source:0,value:0,reference_w:1.},ScanPowerControl {source:0,value:1,reference_w:1.}],
        geometry,objective,factory,1)?;
    let bounds=vec![ScanPeakBound {response:ScanPeakResponse {first_endpoint:0,last_endpoint:1,
        vertices:(0..8).collect(),observable:ScanObservable::SpecificEnthalpy,smoothing_scale:smoothing},
        maximum,scale:10.}];
    // Two primal endpoints plus two peak transpose calls, not sixteen sweeps.
    let mut constrained=ConstrainedScanStudy::new_with_peak_bounds(cx,study,vec![],vec![],vec![],bounds,
        ScanResponseBudget {max_responses:1,max_gradient_entries:32,max_transpose_solves:4},16,8)?;
    let initial=constrained.point().to_vec();
    assert_eq!(constrained.peak_bounds().len(),1);
    let paused=constrained.advance_optimizer(cx,1e-8,120)?;
    if paused.stop!=SqpStop::EvaluationLimit || constrained.point()!=initial.as_slice() {
        return Err("peak-constrained checkpoint changed on budget exhaustion".into());
    }
    if constrained.constraints_at(cx,&initial).is_ok() {return Err("cheap checks silently omitted peak limits".into());}
    constrained.extend_evaluation_budget(127)?;
    let mut gradient_error=0.0_f64;
    if audit {
        let before=constrained.evaluations();
        let (_,rows)=constrained.evaluate_at(cx,&initial)?;
        assert_eq!(constrained.evaluations(),before+1);
        for i in 0..initial.len() {
            let (mut plus,mut minus)=(initial.clone(),initial.clone());plus[i]+=1e-5;minus[i]-=1e-5;
            let (_,p)=constrained.evaluate_at(cx,&plus)?;let (_,m)=constrained.evaluate_at(cx,&minus)?;
            let fd=(p.residuals[0]-m.residuals[0])/2e-5;
            gradient_error=gradient_error.max((fd-rows.jacobian[i]).abs());
        }
        let used=constrained.evaluations();let gate=CancelGate::new_clock_free();gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled=Cx::new(&gate,arena,StreamKey {seed:123,kernel_id:995,tile:0,iteration:0},
                Budget::INFINITE,ExecMode::Deterministic);
            assert!(constrained.advance_optimizer(&cancelled,1e-8,1).is_err());
        });
        assert_eq!(constrained.evaluations(),used);assert_eq!(constrained.point(),initial.as_slice());
    }
    let report=constrained.advance_optimizer(cx,1e-8,if audit {0} else {120})?;
    let before=constrained.evaluations();let sources=constrained.sources(cx)?;
    let powers=sources[0].absorbed_power_values(3)?;
    assert_eq!(powers[2],0.);assert_eq!(constrained.evaluations(),before);
    constrained.extend_evaluation_budget(128)?;
    let (replay,rows)=constrained.replay_with_constraints(cx)?;
    if (replay.physical.value-report.solution.f).abs()>1e-10
        || (rows.residuals[0]-constrained.sample().ci[0]).abs()>1e-10
    {return Err("peak-constrained replay disagrees with the accepted point".into());}
    let final_max=replay.physical.final_step.solution().specific_enthalpy_j_kg.iter()
        .copied().fold(f64::NEG_INFINITY,f64::max);
    // Heating is nonnegative and all nodes stay on a uniform-temperature latent
    // plateau in this fixture, so its largest sampled h is at the final endpoint.
    assert!((final_max-rows.sampled_peak_values[0]).abs()<1e-8);
    Ok(PeakFit {report,powers,constraints:rows,
        temperature:replay.physical.final_step.solution().temperature.clone(),
        evaluations:constrained.evaluations(),gradient_error})
}

pub(super) fn run_cli(args:&[String]) -> Result<()> {
    let (maximum,smoothing)=match args {[] => (155.,0.05),[a,b] => (a.parse()?,b.parse()?),
        _=>return Err("--peak-limits [maximum_enthalpy_J_per_kg smoothing_J_per_kg]".into())};
    let fit=with_context(|cx|fit_peak(cx,maximum,smoothing,false))?;
    println!("# status={:?} objective={:.12e} kkt={:?} evaluations={}",
        fit.report.stop,fit.report.solution.f,fit.report.solution.kkt,fit.evaluations);
    println!("# 8 vertices x 2 endpoints: sampled_max_J_per_kg={:?} envelope_J_per_kg={:?} smoothing_margin_J_per_kg={:?}",
        fit.constraints.sampled_peak_values,fit.constraints.peak_values,fit.constraints.peak_smoothing_gaps);
    println!("# limit_J_per_kg={maximum} smoothing_J_per_kg={smoothing} final_temperature_K={:?}",fit.temperature);
    println!("time_s,absorbed_power_w\n0,{:.12e}\n1,{:.12e}\n2,0",fit.powers[0],fit.powers[1]);
    if !fit.report.solution.converged || fit.constraints.max_scaled_violation()>1e-8 {
        return Err("sampled peak or KKT tolerance not met; schedule is only a candidate".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_one_peak_row_limits_all_selected_nodes_and_times_while_replaying_the_fit() {
        for smoothing in [0.05,0.2] {
            let fit=with_context(|cx|fit_peak(cx,155.,smoothing,false)).unwrap();
            assert!(fit.report.solution.converged);
            assert_eq!(fit.constraints.residuals.len(),1);
            assert!(fit.constraints.max_scaled_violation()<=1e-8);
            assert!(fit.constraints.peak_values[0]<=155.+1e-7);
            assert!(fit.constraints.sampled_peak_values[0]<=fit.constraints.peak_values[0]);
            assert!(fit.constraints.peak_smoothing_gaps[0]<=smoothing*16_f64.ln()+1e-12);
            assert!(fit.report.solution.f>1e-4);
            assert!(fit.temperature.iter().all(|&t|t==310.));
            assert_eq!(fit.powers[2],0.);assert!(fit.evaluations<=128);
        }
    }
    #[test]
    fn g3_peak_constraint_chains_joint_geometry_and_power_gradients_without_changing_accepted_state() {
        let fit=with_context(|cx|fit_peak(cx,155.,0.05,true)).unwrap();
        assert!(fit.gradient_error<2e-6,"{}",fit.gradient_error);
        assert_eq!(fit.report.solution.iters,0);
        assert_eq!(fit.constraints.jacobian.len(),4);
    }
    #[test]
    fn g0_peak_cli_rejects_invalid_scales_and_limits_before_physics() {
        for args in [vec!["155"],vec!["nan","0.05"],vec!["150","0.05"],vec!["155","0"],vec!["155","-1"]] {
            assert!(run_cli(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err());
        }
    }
}

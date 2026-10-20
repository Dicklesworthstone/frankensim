//! Output-state limits, including heat limits hidden by an isothermal plateau.
use super::*;
use fs_ascent::{SqpRunReport,SqpStop};
use fs_ascent::scan_design::{ScanGeometryControl,constraints::{ConstrainedScanStudy,
    ScanStateBound,ScanResponse,ScanResponseBudget,ScanConstraintEvaluation}};

struct StateFit {
    report:SqpRunReport,powers:Vec<f64>,constraints:ScanConstraintEvaluation,
    predictions:Vec<f64>,temperature:Vec<f64>,evaluations:usize,gradient_error:f64,
}

fn fit(cx:&Cx<'_>, first_limit:f64, final_limit:f64, temperature:bool, audit:bool) -> Result<StateFit> {
    let old = if temperature {20.} else {150.};
    let base_value = if temperature {302.} else {old};
    if !first_limit.is_finite() || !final_limit.is_finite()
        || first_limit<=base_value || final_limit<=first_limit
        || (!temperature && final_limit>=200.) || (temperature && final_limit>=310.)
    { return Err("state limits must increase strictly inside the selected caloric branch".into()); }
    let (complex,positions)=fs_conduction::fixtures::box_grid([1;3],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let material=ConductivityModel::isotropic_declared(2.)?;
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([118;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),vec![],0.)?;
    let source=MovingHeatSource::new(std::array::from_fn(|_|TimeSignal::Constant(
        QtyAny::new(0.5,Dims([1,0,0,0,0,0])))),[2.;3],TimeSignal::Table {
            times:vec![0.,1.,2.],values:vec![1.,1.,0.],dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,
        MovingHeatBudget {max_vertices:8,max_elements:6,min_active_vertices:1,power_relative_tolerance:1e-11})?;
    let factory=|cx:&Cx<'_>,sources:Vec<MovingHeatSource>|
        ->std::result::Result<ScanningEnthalpyRun<'_, '_, '_>,ScanSensitivityError> {
        let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&boundary,material:&material,
            element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
            EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}],
                of_element:&[0;6],phase_conductivity:None},&[old;8],
            EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:6},max_steps:2})?;
        Ok(ScanningEnthalpyRun::new(cx,base,sources,1,None)?)
    };
    let objective=ScanPowerObjective {endpoints_s:vec![1.,2.],observations:vec![
        ScanObservation {endpoint:0,vertex:0,observable:ScanObservable::SpecificEnthalpy,target:old+2.,scale:10.},
        ScanObservation {endpoint:1,vertex:0,observable:ScanObservable::SpecificEnthalpy,target:old+7.,scale:10.}],
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
            linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
            max_newton_iterations:64,energy_tolerance_j:1e-8},
        transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
        budget:ScanObjectiveBudget {max_steps:2,max_observations:2,max_retained_vertices:16,
            max_retained_source_vertices:16,max_schedule_controls:16,max_retained_radiation_entries:0},
        source_reference_time_s:0.};
    let geometry=if audit {vec![
        ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:0.5,scale_m:0.25},
        ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:2.}]} else {vec![]};
    let study=ScanPowerStudy::new_joint(cx,vec![source],vec![
        ScanPowerControl {source:0,value:0,reference_w:1.},ScanPowerControl {source:0,value:1,reference_w:1.}],
        geometry,objective,factory,1)?;
    let observable=if temperature {ScanObservable::Temperature} else {ScanObservable::SpecificEnthalpy};
    let scale=if temperature {1.} else {10.};
    let bounds=vec![ScanStateBound {response:ScanResponse {endpoint:0,vertex:0,observable},
        lower:Some(base_value),upper:Some(first_limit),scale},
        ScanStateBound {response:ScanResponse {endpoint:1,vertex:0,observable},
            lower:None,upper:Some(final_limit),scale}];
    let mut constrained=ConstrainedScanStudy::new_with_state_bounds(cx,study,vec![],vec![],bounds,
        ScanResponseBudget {max_responses:2,max_gradient_entries:64,max_transpose_solves:5},16)?;
    let initial=constrained.point().to_vec();
    let paused=constrained.advance_optimizer(cx,1e-8,120)?;
    if paused.stop!=SqpStop::EvaluationLimit || constrained.point()!=initial.as_slice() {
        return Err("state-constrained checkpoint changed on budget exhaustion".into());
    }
    // A cheap input-only check must not silently claim all state limits passed.
    if constrained.constraints_at(cx,&initial).is_ok() || constrained.evaluations()!=1 {
        return Err("state constraints were omitted or charged as input-only checks".into());
    }
    constrained.extend_evaluation_budget(127)?;
    let mut gradient_error=0.0_f64;
    if audit {
        let before=constrained.evaluations();
        let (_,rows)=constrained.evaluate_at(cx,&initial)?;
        if constrained.evaluations()!=before+1 {return Err("response limits duplicated the physical evaluation".into());}
        for i in 0..initial.len() {
            let(mut plus,mut minus)=(initial.clone(),initial.clone());plus[i]+=1e-5;minus[i]-=1e-5;
            let (_,p)=constrained.evaluate_at(cx,&plus)?;let (_,m)=constrained.evaluate_at(cx,&minus)?;
            for row in 0..rows.residuals.len() {
                let fd=(p.residuals[row]-m.residuals[row])/2e-5;
                gradient_error=gradient_error.max((fd-rows.jacobian[row*initial.len()+i]).abs());
            }
        }
        if constrained.point()!=initial.as_slice() {return Err("trial response evaluation changed accepted design".into());}
        let used=constrained.evaluations();let gate=CancelGate::new_clock_free();gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled=Cx::new(&gate,arena,StreamKey {seed:118,kernel_id:994,tile:0,iteration:0},
                Budget::INFINITE,ExecMode::Deterministic);
            assert!(constrained.advance_optimizer(&cancelled,1e-8,1).is_err());
        });
        assert_eq!(constrained.evaluations(),used);assert_eq!(constrained.point(),initial.as_slice());
    }
    let report=constrained.advance_optimizer(cx,1e-8,if audit {0} else {120})?;
    let used=constrained.evaluations();let sources=constrained.sources(cx)?;
    let powers=sources[0].absorbed_power_values(3)?;
    if powers[2]!=0. || used!=constrained.evaluations() {return Err("export changed shutdown or spent work".into());}
    constrained.extend_evaluation_budget(128)?;
    let (replay,constraints)=constrained.replay_with_constraints(cx)?;
    if (replay.physical.value-report.solution.f).abs()>1e-10 {return Err("state-constrained replay objective differs".into());}
    for (&a,&b) in constraints.residuals.iter().zip(&constrained.sample().ci) {
        if (a-b).abs()>1e-10 {return Err("state-constrained replay inequalities differ".into());}
    }
    Ok(StateFit {report,powers,constraints,predictions:replay.physical.predictions,
        temperature:replay.physical.final_step.solution().temperature.clone(),
        evaluations:constrained.evaluations(),gradient_error})
}

pub(super) fn run_cli(args:&[String]) -> Result<()> {
    let (first,last)=match args {[] => (151.,155.),[a,b] => (a.parse()?,b.parse()?),
        _=>return Err("--state-limits [first_enthalpy_max_J_per_kg final_enthalpy_max_J_per_kg]".into())};
    let fit=with_context(|cx|fit(cx,first,last,false,false))?;
    println!("# status={:?} objective={:.12e} kkt={:?} evaluations={}",
        fit.report.stop,fit.report.solution.f,fit.report.solution.kkt,fit.evaluations);
    println!("# target_enthalpy_J_per_kg=[152,157] predicted={:?} state_values_J_per_kg={:?}",fit.predictions,fit.constraints.state_values);
    println!("# state_caps_J_per_kg=[{first},{last}] scaled_violation={:.12e} final_temperature_K={:?}",
        fit.constraints.max_scaled_violation(),fit.temperature);
    println!("time_s,absorbed_power_w\n0,{:.12e}\n1,{:.12e}\n2,0",fit.powers[0],fit.powers[1]);
    if !fit.report.solution.converged || fit.constraints.max_scaled_violation()>1e-8 {
        return Err("sampled state limits or KKT tolerance not met; output is only a candidate".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_enthalpy_limits_remain_effective_on_an_isothermal_melting_plateau() {
        let fit=with_context(|cx|fit(cx,151.,155.,false,false)).unwrap();
        assert!(fit.report.solution.converged);
        assert!((fit.powers[0]-1.).abs()<2e-5 && (fit.powers[1]-4.).abs()<2e-5);
        assert!((fit.report.solution.f-0.025).abs()<1e-7);
        assert!(fit.temperature.iter().all(|&t|t==310.));
        assert!(fit.constraints.max_scaled_violation()<=1e-8);
        assert!(fit.report.solution.nu[1]>0. && fit.report.solution.nu[2]>0.);
        assert!(fit.evaluations<=128);
    }
    #[test]
    fn g1_temperature_limits_use_the_caloric_chain_not_enthalpy_units() {
        let fit=with_context(|cx|fit(cx,302.1,302.5,true,false)).unwrap();
        assert!(fit.report.solution.converged);
        assert!((fit.powers[0]-1.).abs()<3e-5 && (fit.powers[1]-4.).abs()<3e-5);
        assert!((fit.constraints.state_values[0]-302.1).abs()<1e-7);
        assert!((fit.constraints.state_values[1]-302.5).abs()<1e-7);
        assert!(fit.constraints.max_scaled_violation()<=1e-8);
    }
    #[test]
    fn g3_joint_state_jacobians_and_failed_trials_preserve_the_checkpoint() {
        let fit=with_context(|cx|fit(cx,151.5,156.,false,true)).unwrap();
        assert!(fit.gradient_error<2e-7);
        assert_eq!(fit.powers,vec![1.,1.,0.]);
        assert_eq!(fit.report.solution.iters,0);
    }
    #[test]
    fn g0_malformed_state_limit_commands_refuse() {
        for (first,last) in [(150.,155.),(155.,151.),(151.,200.),(f64::NAN,155.)] {
            assert!(with_context(|cx|fit(cx,first,last,false,false)).is_err());
        }
        assert!(run_cli(&["151".into()]).is_err());
    }
}

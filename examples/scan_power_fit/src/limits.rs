//! Feasible best-fit heating when the requested trajectory exceeds a design cap.
use super::*;
use fs_ascent::{SqpRunReport,SqpStop};
use fs_ascent::scan_design::{ScanGeometryControl,constraints::{ConstrainedScanStudy,
    ScanControlBound,ScanEnergyCap,ScanConstraintEvaluation}};

struct LimitedFit {
    report:SqpRunReport,powers:Vec<f64>,constraints:ScanConstraintEvaluation,
    predictions:Vec<f64>,loss:f64,evaluations:usize,gradient_error:f64,
}

fn fit(cx:&Cx<'_>,maximum_w:f64,maximum_j:f64,audit_geometry:bool) -> Result<LimitedFit> {
    if !maximum_w.is_finite() || maximum_w<=0.01 || !maximum_j.is_finite() || maximum_j<=0.0 {
        return Err("--limits requires maximum power > 0.01 W and positive maximum absorbed joules".into());
    }
    let (complex,positions)=fs_conduction::fixtures::box_grid([1;3],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let material=ConductivityModel::isotropic_declared(2.)?;
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([111;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),vec![],0.)?;
    let source=MovingHeatSource::new(std::array::from_fn(|_|TimeSignal::Constant(
        QtyAny::new(0.5,Dims([1,0,0,0,0,0])))),[2.;3],TimeSignal::Table {
            times:vec![0.,1.,2.],values:vec![1.,1.,0.],dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,
        MovingHeatBudget {max_vertices:8,max_elements:6,min_active_vertices:1,power_relative_tolerance:1e-11})?;
    let mut sources=vec![source.clone()];
    if audit_geometry { sources.push(source.with_absorbed_power_values(&[0.25;3],0.,3)?); }
    let factory=|cx:&Cx<'_>,sources:Vec<MovingHeatSource>|
        ->std::result::Result<ScanningEnthalpyRun<'_, '_, '_>,ScanSensitivityError> {
        let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&boundary,material:&material,
            element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
            EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}],
                of_element:&[0;6],phase_conductivity:None},&[150.;8],
            EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:6},max_steps:2})?;
        Ok(ScanningEnthalpyRun::new(cx,base,sources,2,None)?)
    };
    let objective=ScanPowerObjective {endpoints_s:vec![1.,2.],observations:vec![
        ScanObservation {endpoint:0,vertex:0,observable:ScanObservable::SpecificEnthalpy,target:152.,scale:10.},
        ScanObservation {endpoint:1,vertex:0,observable:ScanObservable::SpecificEnthalpy,target:157.,scale:10.}],
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
            linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
            max_newton_iterations:64,energy_tolerance_j:1e-8},
        transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
        budget:ScanObjectiveBudget {max_steps:2,max_observations:2,max_retained_vertices:16,
            max_retained_source_vertices:32,max_schedule_controls:32,max_retained_radiation_entries:0},
        source_reference_time_s:0.};
    let geometry=if audit_geometry {vec![
        ScanGeometryControl::CenterValue {source:0,axis:0,value:0,reference_m:0.5,scale_m:0.25},
        ScanGeometryControl::HalfWidth {source:0,axis:0,reference_m:2.}]} else {vec![]};
    let study=ScanPowerStudy::new_joint(cx,sources,vec![
        ScanPowerControl {source:0,value:0,reference_w:1.},ScanPowerControl {source:0,value:1,reference_w:1.}],
        geometry,objective,factory,1)?;
    let mut bounds=vec![
        ScanControlBound {coordinate:0,lower:Some(0.01),upper:Some(maximum_w),scale:1.},
        ScanControlBound {coordinate:1,lower:Some(0.01),upper:Some(maximum_w),scale:1.}];
    if audit_geometry { bounds.extend([
        ScanControlBound {coordinate:2,lower:Some(-0.2),upper:Some(0.8),scale:0.5},
        ScanControlBound {coordinate:3,lower:Some(1.5),upper:Some(2.5),scale:1.}]); }
    let mut constrained=ConstrainedScanStudy::new(cx,study,bounds,
        vec![ScanEnergyCap {start_time_s:0.,end_time_s:2.,maximum_j,scale_j:1.}],32)?;
    let initial=constrained.point().to_vec();
    let paused=constrained.advance_optimizer(cx,1e-8,120)?;
    if paused.stop!=SqpStop::EvaluationLimit || constrained.point()!=initial.as_slice() {
        return Err("constrained checkpoint changed after evaluation-budget exhaustion".into());
    }
    let mut gradient_error=0.0_f64;
    if audit_geometry {
        let rows=constrained.constraints_at(cx,&initial)?;
        if (rows.absorbed_energy_j[0]-2.5).abs()>1e-13 {return Err("energy cap dropped an unselected heater".into());}
        for i in 0..initial.len() {
            let(mut plus,mut minus)=(initial.clone(),initial.clone());plus[i]+=1e-5;minus[i]-=1e-5;
            let p=constrained.constraints_at(cx,&plus)?;let m=constrained.constraints_at(cx,&minus)?;
            for row in 0..rows.residuals.len() {
                let fd=(p.residuals[row]-m.residuals[row])/2e-5;
                gradient_error=gradient_error.max((fd-rows.jacobian[row*initial.len()+i]).abs());
            }
        }
        let last=rows.residuals.len()-1;
        if rows.jacobian[last*4+2]!=0. || rows.jacobian[last*4+3]!=0. {
            return Err("changing a domain-absorbed footprint spuriously changes total joules".into());
        }
        assert!(constrained.constraints_at(cx,&[]).is_err());
        assert_eq!(constrained.evaluations(),1);
        let gate=CancelGate::new_clock_free();gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled=Cx::new(&gate,arena,StreamKey {seed:111,kernel_id:991,tile:0,iteration:0},
                Budget::INFINITE,ExecMode::Deterministic);
            assert!(constrained.advance_optimizer(&cancelled,1e-8,1).is_err());
        });
        assert_eq!(constrained.point(),initial.as_slice());
        assert_eq!(constrained.evaluations(),1);
    }
    constrained.extend_evaluation_budget(127)?;
    let report=constrained.advance_optimizer(cx,1e-8,if audit_geometry {0} else {120})?;
    let used=constrained.evaluations();
    let sources=constrained.sources(cx)?;
    let powers=sources[0].absorbed_power_values(3)?;
    assert_eq!(powers[2],0.);assert_eq!(used,constrained.evaluations());
    let constraints=constrained.constraints_at(cx,constrained.point())?;
    let actual_energy=sources.iter().map(|s|s.integrate_absorbed_power(cx,0.,2.,3).map(|e|e.energy_j))
        .collect::<std::result::Result<Vec<_>,_>>()?.iter().sum::<f64>();
    if (constraints.absorbed_energy_j[0]-actual_energy).abs()>1e-11 {
        return Err("exported source energy differs from constrained design".into());
    }
    constrained.extend_evaluation_budget(128)?;
    let replay=constrained.replay(cx)?;
    if (replay.physical.value-report.solution.f).abs()>1e-10 {return Err("constrained replay objective differs".into());}
    Ok(LimitedFit {report,powers,constraints,predictions:replay.physical.predictions,
        loss:replay.physical.value,evaluations:constrained.evaluations(),gradient_error})
}

pub(super) fn run_cli(args:&[String]) -> Result<()> {
    let (maximum_w,maximum_j)=match args {
        []=>(3.5,5.),[w,j]=>(w.parse()?,j.parse()?),
        _=>return Err("--limits [maximum_power_W maximum_absorbed_energy_J]".into()),
    };
    let result=with_context(|cx|fit(cx,maximum_w,maximum_j,false))?;
    println!("# status={:?} objective={:.12e} kkt={:?} evaluations={}",
        result.report.stop,result.loss,result.report.solution.kkt,result.evaluations);
    println!("# target_enthalpy_J_per_kg=[152,157] predicted={:?}",result.predictions);
    println!("# max_power_W={maximum_w} max_absorbed_energy_J={maximum_j} actual_energy_J={:?} scaled_violation={:.12e}",
        result.constraints.absorbed_energy_j,result.constraints.max_scaled_violation());
    println!("time_s,absorbed_power_w\n0,{:.12e}\n1,{:.12e}\n2,0",result.powers[0],result.powers[1]);
    if !result.report.solution.converged || result.constraints.max_scaled_violation()>1e-8 {
        return Err("constrained fit did not meet its KKT/feasibility tolerance; output is only a candidate".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_energy_limited_fit_retains_target_mismatch_and_active_constraint() {
        let fit=with_context(|cx|fit(cx,3.5,5.,false)).unwrap();
        assert!(fit.report.solution.converged);
        assert!((fit.powers[0]-2.).abs()<2e-5 && (fit.powers[1]-3.).abs()<2e-5);
        assert!((fit.loss-0.02).abs()<1e-7);
        assert!(fit.constraints.max_scaled_violation()<=1e-8);
        assert!(*fit.report.solution.nu.last().unwrap()>0.0);
        assert!(fit.evaluations<=128);
    }
    #[test]
    fn g1_peak_power_limit_uses_native_constrained_optimum_not_clipping() {
        let fit=with_context(|cx|fit(cx,3.,20.,false)).unwrap();
        assert!(fit.report.solution.converged);
        assert!((fit.powers[0]-3.).abs()<2e-5 && (fit.powers[1]-3.).abs()<2e-5);
        assert!((fit.loss-0.01).abs()<1e-7);
        assert!(fit.constraints.max_scaled_violation()<=1e-8);
    }
    #[test]
    fn g3_joint_constraint_jacobians_include_fixed_heaters_and_preserve_exact_shutdowns() {
        let fit=with_context(|cx|fit(cx,3.5,5.,true)).unwrap();
        assert!(fit.gradient_error<1e-8);
        assert_eq!(fit.powers[2],0.);
        assert!((fit.constraints.absorbed_energy_j[0]-2.5).abs()<1e-13);
    }
    #[test]
    fn g0_bad_limit_commands_refuse_before_a_physical_solve() {
        for (w,j) in [(0.,5.),(3.,-1.),(f64::NAN,5.),(3.,f64::INFINITY)] {
            assert!(with_context(|cx|fit(cx,w,j,false)).is_err());
        }
        assert!(run_cli(&["3".into()]).is_err());
    }
}

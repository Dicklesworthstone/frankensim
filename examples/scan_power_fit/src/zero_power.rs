//! Inclusive zero-power design on the original synthetic enthalpy fixture.
use super::*;
use fs_ascent::{SqpRunReport, SqpStop};
use fs_ascent::scan_design::{ScanPowerCoordinates, ScanGeometryControl,
    constraints::{ConstrainedScanStudy, ScanControlBound, ScanEnergyCap,
        ScanStateBound, ScanPeakBound, ScanResponse, ScanResponseBudget,
        ScanPeakResponse, ScanConstraintEvaluation}};

const USAGE: &str = "--heater-off [first_target_J_per_kg second_target_J_per_kg [initial_first_W initial_second_W]]";

#[derive(Debug, Clone, Copy)]
struct Input { targets: [f64; 2], initial: [f64; 2] }
impl Input {
    fn parse(args: &[String]) -> Result<Self> {
        let mut out = Self { targets: [154.0, 152.0], initial: [0.0, 0.0] };
        if ![0, 2, 4].contains(&args.len()) { return Err(USAGE.into()); }
        for (i, arg) in args.iter().enumerate() {
            if i < 2 { out.targets[i] = arg.parse()?; }
            else { out.initial[i-2] = arg.parse()?; }
        }
        out.validate()?;
        Ok(out)
    }
    fn validate(self) -> Result<()> {
        if self.targets.iter().any(|&t| !t.is_finite() || t <= 100.0 || t >= 190.0)
            || self.initial.iter().any(|&p| !p.is_finite() || p < 0.0)
            || self.initial.iter().sum::<f64>() >= 40.0
        { return Err("synthetic fixture needs targets inside (100,190) J/kg and nonnegative initial watts summing to less than 40".into()); }
        Ok(())
    }
}

struct OffFit {
    report: SqpRunReport,
    powers: Vec<f64>,
    predictions: Vec<f64>,
    temperature: Vec<f64>,
    constraints: ScanConstraintEvaluation,
    evaluations: usize,
    checkpoint_work: Option<ScanCheckpointWork>,
    gradient_error: f64,
}

fn fit(cx: &Cx<'_>, input: Input, scales: [f64; 2], checkpointed: bool, audit: bool) -> Result<OffFit> {
    input.validate()?;
    let endpoints = 8usize;
    let (complex, positions) = fs_conduction::fixtures::box_grid([1;3], [1.;3]);
    let mesh = ConductionMesh::new(complex, positions)?;
    let material = ConductivityModel::isotropic_declared(2.0)?;
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([139;32]),
        [(-200.,280.,0.), (100.,310.,0.), (200.,310.,1.), (900.,380.,1.)].into_iter()
            .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
                liquid_mass_fraction:f, bulk_density_kg_m3:1. }).collect())?;
    let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface", ThermalBc::adiabatic())?, &Environment::earth_lab(), vec![], 0.0)?;
    let source = MovingHeatSource::new(std::array::from_fn(|_| TimeSignal::Constant(
        QtyAny::new(0.5, Dims([1,0,0,0,0,0])))), [2.;3], TimeSignal::Table {
            times:vec![0.,1.,2.], values:vec![input.initial[0],input.initial[1],0.],
            dims:COMPONENT_POWER_DIMS, interp:Interp::Hold }, 0.0,
        MovingHeatBudget { max_vertices:8, max_elements:6, min_active_vertices:1,
            power_relative_tolerance:1e-11 })?;
    let factory = |cx: &Cx<'_>, sources: Vec<MovingHeatSource>|
        -> std::result::Result<ScanningEnthalpyRun<'_, '_, '_>, ScanSensitivityError> {
        let base = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
            boundary:&boundary, material:&material, element_materials:None, interfaces:None,
            source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)) },
            EnthalpyHistoryMaterials { table:&[ReferenceEnthalpyMaterial { curve:&curve,
                reference_density_kg_m3:1. }], of_element:&[0;6], phase_conductivity:None }, &[150.;8],
            EnthalpyHistoryBudget { spatial:EnthalpyBudget { max_vertices:8, max_elements:6 },
                max_steps:endpoints })?;
        Ok(ScanningEnthalpyRun::new(cx, base, sources, 1, None)?)
    };
    let retained = if checkpointed {2} else {endpoints};
    let objective = ScanPowerObjective { endpoints_s:(1..=endpoints).map(|i| i as f64 / 4.0).collect(),
        observations:vec![ScanObservation { endpoint:3, vertex:0, observable:ScanObservable::SpecificEnthalpy,
            target:input.targets[0], scale:10. }, ScanObservation { endpoint:7, vertex:0,
            observable:ScanObservable::SpecificEnthalpy, target:input.targets[1], scale:10. }],
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig { newton:NewtonKrylovConfig { absolute_tolerance:1e-11, relative_tolerance:1e-12,
            linear_restart:8, max_linear_cycles:16, forcing_minimum:1e-12, forcing_maximum:1e-3,
            ..Default::default() }, max_newton_iterations:64, energy_tolerance_j:1e-8 },
        transpose:LinearConfig { tolerance:1e-10, max_iterations:1000, restart:8 },
        budget:ScanObjectiveBudget { max_steps:endpoints, max_observations:2,
            max_retained_vertices:8*retained, max_retained_source_vertices:8*retained,
            max_schedule_controls:16, max_retained_radiation_entries:0 }, source_reference_time_s:0.0 };
    let geometry = if audit { vec![
        ScanGeometryControl::CenterValue { source:0, axis:0, value:0, reference_m:0.5, scale_m:0.25 },
        ScanGeometryControl::HalfWidth { source:0, axis:0, reference_m:2.0 }] } else { vec![] };
    let mut study = ScanPowerStudy::new_joint_with_power_coordinates(cx, vec![source], vec![
        ScanPowerControl { source:0, value:0, reference_w:scales[0] },
        ScanPowerControl { source:0, value:1, reference_w:scales[1] }], geometry,
        ScanPowerCoordinates::NonnegativeLinear, objective, factory, 1)?;
    if checkpointed {
        study.set_checkpointing(cx, ScanCheckpointBudget { segment_steps:2,
            max_checkpoint_vertices:40, max_forward_solves:16 })?;
    }
    if audit {
        // An unconstrained solver cannot silently omit the new physical domain.
        assert!(study.begin_optimizer(cx, 8).is_err());
        assert_eq!(study.evaluations(), 0);
        assert!(study.powers_at(&[0.,0.,0.]).is_err());
    }
    let mut constrained = if audit {
        ConstrainedScanStudy::new_with_peak_bounds(cx, study,
            vec![ScanControlBound { coordinate:0, lower:None, upper:Some(20.0), scale:2.0 }],
            vec![ScanEnergyCap { start_time_s:0.0, end_time_s:2.0, maximum_j:40.0, scale_j:3.0 }],
            vec![ScanStateBound { response:ScanResponse { endpoint:7, vertex:0,
                observable:ScanObservable::LiquidMassFraction }, lower:None, upper:Some(0.9), scale:0.1 }],
            vec![ScanPeakBound { response:ScanPeakResponse { first_endpoint:0, last_endpoint:7,
                vertices:(0..8).collect(), observable:ScanObservable::SpecificEnthalpy, smoothing_scale:0.05 },
                maximum:190.0, scale:10.0 }],
            ScanResponseBudget { max_responses:2, max_gradient_entries:64, max_transpose_solves:24 },64,16)?
    } else { ConstrainedScanStudy::new(cx, study, vec![], vec![], 4)? };
    let initial = constrained.point().to_vec();
    let paused = constrained.advance_optimizer(cx, 1e-9, 100)?;
    if paused.stop != SqpStop::EvaluationLimit || constrained.point() != initial.as_slice() {
        return Err("nonnegative study changed its accepted point while exhausted".into());
    }
    constrained.extend_evaluation_budget(127)?;
    let mut gradient_error = 0.0_f64;
    if audit {
        for point in [initial.clone(), vec![2.0/scales[0],1.0/scales[1],0.2,0.1]] {
            let (base, rows) = constrained.evaluate_at(cx, &point)?;
            assert_eq!(&rows.residuals[..2], &[-point[0],-point[1]]);
            assert_eq!(rows.nonnegative_power_values_w, base.powers_w);
            for i in 0..point.len() {
                let epsilon = 1e-4;
                let mut plus = point.clone(); plus[i] += epsilon;
                let (p, pr) = constrained.evaluate_at(cx, &plus)?;
                let (fd, derivatives) = if i < 2 && point[i] == 0.0 {
                    let mut twice = point.clone(); twice[i] += 2.0*epsilon;
                    let (q, qr) = constrained.evaluate_at(cx, &twice)?;
                    ((-3.0*base.physical.value+4.0*p.physical.value-q.physical.value)/(2.0*epsilon),
                        rows.residuals.iter().zip(&pr.residuals).zip(&qr.residuals)
                            .map(|((&a,&b),&c)|(-3.0*a+4.0*b-c)/(2.0*epsilon)).collect::<Vec<_>>())
                } else {
                    let mut minus = point.clone(); minus[i] -= epsilon;
                    let (m, mr) = constrained.evaluate_at(cx, &minus)?;
                    ((p.physical.value-m.physical.value)/(2.0*epsilon),
                        pr.residuals.iter().zip(&mr.residuals).map(|(&a,&b)|(a-b)/(2.0*epsilon)).collect())
                };
                gradient_error = gradient_error.max((fd-base.gradient[i]).abs());
                for (row, &fd) in derivatives.iter().enumerate() {
                    gradient_error = gradient_error.max((fd-rows.jacobian[row*point.len()+i]).abs());
                }
            }
        }
        let spent = constrained.evaluations();
        let mut negative = initial.clone(); negative[0] = -1e-12;
        assert!(constrained.evaluate_at(cx, &negative).is_err());
        assert_eq!(constrained.evaluations(), spent);
        let gate = CancelGate::new_clock_free(); gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(&gate, arena, StreamKey {seed:139,kernel_id:999,tile:0,iteration:0},
                Budget::INFINITE, ExecMode::Deterministic);
            assert!(constrained.advance_optimizer(&cancelled,1e-9,1).is_err());
        });
        assert_eq!(constrained.evaluations(), spent);
        assert_eq!(constrained.point(), initial.as_slice());
    }
    let report = constrained.advance_optimizer(cx, 1e-9, if audit {0} else {100})?;
    let spent = constrained.evaluations();
    let sources = constrained.sources(cx)?;
    let powers = sources[0].absorbed_power_values(3)?;
    if powers[2] != 0.0 || constrained.evaluations() != spent {
        return Err("source export changed the shutdown or physical evaluation count".into());
    }
    constrained.extend_evaluation_budget(128)?;
    let (replay, constraints) = constrained.replay_with_constraints(cx)?;
    if (replay.physical.value-report.solution.f).abs() > 1e-10 {
        return Err("off-state physical replay changed the accepted objective".into());
    }
    if constraints.residuals != constrained.sample().ci {
        return Err("off-state physical replay changed accepted constraints".into());
    }
    Ok(OffFit {report,powers,predictions:replay.physical.predictions,
        temperature:replay.physical.final_step.solution().temperature.clone(),constraints,
        evaluations:constrained.evaluations(),checkpoint_work:replay.checkpoint_work,gradient_error})
}

pub(super) fn run_cli(args: &[String]) -> Result<()> {
    let input = Input::parse(args)?;
    let result = with_context(|cx|fit(cx,input,[1.0,1.0],true,false))?;
    println!("# status={:?} objective={:.12e} kkt={:?} evaluations={}", result.report.stop,
        result.report.solution.f,result.report.solution.kkt,result.evaluations);
    println!("# targets_J_per_kg={:?} predictions_J_per_kg={:?} final_temperature_K={:?}",
        input.targets,result.predictions,result.temperature);
    println!("# powers_W={:?} nonnegative_domain_multipliers={:?} checkpoint_work={:?}",
        result.powers,result.report.solution.nu,result.checkpoint_work);
    println!("time_s,absorbed_power_w\n0,{:.17e}\n1,{:.17e}\n2,0",result.powers[0],result.powers[1]);
    if !result.report.solution.converged || result.constraints.max_scaled_violation()>1e-9 {
        return Err("nonnegative fit has not met its KKT/feasibility checks; output is only a candidate".into());
    }
    // A constrained optimum can miss impossible cooling targets. Report the
    // actual mismatch; never relabel it as target recovery or snap small watts.
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_off_state_optima_and_activation_from_zero_use_actual_enthalpy() {
        for (targets,initial,expected,loss) in [([154.,152.],[0.,0.],[3.,0.],0.01),
            ([152.,157.],[0.,0.],[2.,5.],0.0),([149.,148.],[1.,1.],[0.,0.],0.025)] {
            let fitted=with_context(|cx|fit(cx,Input {targets,initial},[0.25,4.],true,false)).unwrap();
            assert!(fitted.report.solution.converged, "{:?}",fitted.report);
            for i in 0..2 { assert!((fitted.powers[i]-expected[i]).abs()<2e-5); }
            assert_eq!(fitted.powers[2],0.0);
            assert!((fitted.report.solution.f-loss).abs()<1e-8);
            assert!(fitted.temperature.iter().all(|&t|t==310.0));
            assert!(fitted.constraints.max_scaled_violation()<=1e-9);
            assert!(fitted.evaluations<=128);
            if targets==[154.,152.] {assert!(fitted.report.solution.nu[1]>0.0);}
        }
    }
    #[test]
    fn g5_short_tape_preserves_inclusive_boundary_fit_and_resume() {
        let input=Input {targets:[154.,152.],initial:[1.,1.]};
        let full=with_context(|cx|fit(cx,input,[1.,1.],false,false)).unwrap();
        let short=with_context(|cx|fit(cx,input,[1.,1.],true,false)).unwrap();
        assert!(full.report.solution.converged && short.report.solution.converged);
        assert_eq!(full.powers,short.powers);assert_eq!(full.predictions,short.predictions);
        assert_eq!(full.constraints,short.constraints);assert_eq!(full.evaluations,short.evaluations);
        let work=short.checkpoint_work.unwrap();
        assert_eq!(work.original_forward_solves,8);assert_eq!(work.recomputed_forward_solves,8);
        assert!(full.checkpoint_work.is_none());
    }
    #[test]
    fn g3_zero_power_and_joint_response_gradients_do_not_use_the_logarithmic_chain() {
        for checkpointed in [false,true] {
            let fitted=with_context(|cx|fit(cx,Input {targets:[154.,152.],initial:[0.,0.]},
                [0.25,4.0],checkpointed,true)).unwrap();
            assert!(fitted.gradient_error<2e-5,"{}",fitted.gradient_error);
            assert_eq!(fitted.powers,vec![0.0,0.0,0.0]);
            // Automatic rows precede declared power, energy, state and peak rows.
            assert_eq!(fitted.constraints.nonnegative_power_values_w,vec![0.0,0.0]);
            assert_eq!(&fitted.constraints.jacobian[..8],&[-1.,0.,0.,0.,0.,-1.,0.,0.]);
            assert!(fitted.constraints.jacobian[3*4]>0.0 && fitted.constraints.jacobian[3*4+1]>0.0);
        }
    }
    #[test]
    fn g0_zero_is_valid_but_negative_and_nonfinite_heater_inputs_are_not() {
        assert!(Input::parse(&[]).is_ok());
        for initial in [[-1.,0.],[f64::NAN,0.],[0.,f64::INFINITY],[20.,20.]] {
            assert!(Input {targets:[154.,152.],initial}.validate().is_err());
        }
        assert!(Input::parse(&["154".into()]).is_err());
        assert!(with_context(|cx|fit(cx,Input {targets:[154.,152.],initial:[0.,0.]},
            [0.,1.],false,false)).is_err());
    }
}

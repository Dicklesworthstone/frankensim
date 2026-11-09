//! Authored reservoir histories and emissivity through whole-scan inverse design.
use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    ScanRadiationScheduleGradient,
    checkpoint::{ScanCheckpointBudget, evaluate_scan_objective_checkpointed_with_peaks},
    responses::{ScanResponse, ScanResponseBudget, peaks::{
        ScanPeakResponse, evaluate_scan_objective_with_peaks}},
    value::{ScanValueBudget, evaluate_scan_values, responses::{
        ScanValueResponseBudget, evaluate_scan_values_with_responses}},
};
use fs_scenario::signal::LinearChirp;

fn radiative_run<'h,'m,'c>(cx:&Cx<'_>,history:&'h ThermalBoundaryHistory<'m>,
    material:&'h ConductivityModel,curve:&'c EquilibriumEnthalpyPhaseCurve,
    old:&[f64],sources:Vec<MovingHeatSource>,radiation:ScanRadiation)
    ->ScanningEnthalpyRun<'h,'m,'c>{
    let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem{boundary:history,material,
        element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
        EnthalpyHistoryMaterials{table:&[ReferenceEnthalpyMaterial{curve,reference_density_kg_m3:1.}],
            of_element:&[0],phase_conductivity:None},old,
        EnthalpyHistoryBudget{spatial:EnthalpyBudget{max_vertices:4,max_elements:1},max_steps:8}).unwrap();
    ScanningEnthalpyRun::new(cx,base,sources,4,Some(radiation)).unwrap()
}

fn sourced_radiation(ambient_temperature:TimeSignal,emissivity:f64)->ScanRadiation{
    let mut claims=ClaimSet::new();claims.insert_claim(PropertyClaim{
        key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar{value:emissivity,dims:EMISSIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.),uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
        provenance:Provenance{source:"synthetic reservoir schedule fixture".into(),license:"test-only".into(),artifact:None}}).unwrap();
    let card=MaterialCard::assemble(MaterialStateId{chemistry:"gray fixture".into(),phase:"solid".into(),
        process:"synthetic".into(),revision:0},claims,vec![]).unwrap();
    ScanRadiation{patches:vec![RadiationRegionHistory{region:"surface".into(),
        emissivity:SurfaceEmissivity::from_card("surface",&card,300.,SelectionPolicy::SingleClaimOnly).unwrap(),
        ambient_temperature}],config:AmbientRadiationConfig::default()}
}

fn reservoir(kind:usize,values:&[f64])->TimeSignal{
    match kind{
        0=>constant(values[0],TEMPERATURE_DIMS),
        1=>TimeSignal::Ramp{t_start:0.,t_end:0.006,
            from:QtyAny::new(values[0],TEMPERATURE_DIMS),to:QtyAny::new(values[1],TEMPERATURE_DIMS)},
        2|3=>TimeSignal::Table{times:vec![0.,0.002,0.006],values:values.to_vec(),
            dims:TEMPERATURE_DIMS,interp:if kind==2{Interp::Linear}else{Interp::Hold}},
        _=>unreachable!("finite fixture profile family"),
    }
}
fn observations()->[ScanObservation;2]{
    [observation(0,1,ScanObservable::Temperature,302.,2.),
        observation(3,3,ScanObservable::SpecificEnthalpy,42.,4.)]
}
fn requests()->([ScanResponse;1],[ScanPeakResponse;1]){
    ([ScanResponse{endpoint:1,vertex:0,observable:ScanObservable::Temperature}],
        [ScanPeakResponse{first_endpoint:0,last_endpoint:3,vertices:vec![0,3],
            observable:ScanObservable::Temperature,smoothing_scale:0.5}])
}
fn response_budget()->ScanResponseBudget{
    ScanResponseBudget{max_responses:2,max_gradient_entries:128,max_transpose_solves:10}
}
fn checkpoints()->ScanCheckpointBudget{
    ScanCheckpointBudget{segment_steps:2,max_checkpoint_vertices:12,max_forward_solves:8}
}
fn components(gradient:&ScanRadiationScheduleGradient)->Vec<f64>{
    let mut values=gradient.reservoir_values_k.clone();values.push(gradient.emissivity);values
}

#[test]
fn g3_radiative_schedules_match_whole_history_response_and_peak_differences(){
    let m=mesh();let curve=chart();let material=ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280.,1.),(380.,3.)]).unwrap());
    with_cx(|cx|{
        let hist=history(cx,&m);let old=[10.,20.,30.,40.];
        let p=[0.2,0.25,0.3,1.4,1.3,1.2,3.];let times=[0.001,0.002,0.004,0.006];
        let (points,peaks)=requests();
        for (kind,sampling) in [(0,ScanObjectiveSampling::RightEndpoint),
            (1,ScanObjectiveSampling::EventLeft),(2,ScanObjectiveSampling::EventLeft),
            (3,ScanObjectiveSampling::EventLeft),(3,ScanObjectiveSampling::RightEndpoint)]{
            let values=match kind{0=>vec![293.],1=>vec![291.,297.],_=>vec![291.,295.,298.]};
            let make=|values:&[f64],emissivity|radiative_run(cx,&hist,&material,&curve,&old,
                vec![heater(p)],sourced_radiation(reservoir(kind,values),emissivity));
            let (objective,rows,envelopes)=evaluate_scan_objective_with_peaks(cx,make(&values,0.2),
                &times,&observations(),sampling,config(),linear(),objective_budget(),
                &points,&peaks,response_budget(),8).unwrap();
            let replay=evaluate_scan_objective_checkpointed_with_peaks(cx,make(&values,0.2),
                &times,&observations(),sampling,config(),linear(),objective_budget(),
                &points,&peaks,response_budget(),8,checkpoints()).unwrap();
            assert_eq!(objective.radiation,replay.objective.radiation);
            assert_eq!(rows[0].radiation,replay.responses[0].radiation);
            assert_eq!(envelopes[0].radiation,replay.peaks[0].radiation);
            assert_eq!(objective.sources,replay.objective.sources);
            assert_eq!(rows[0].sources,replay.responses[0].sources);
            assert_eq!(envelopes[0].sources,replay.peaks[0].sources);
            assert_eq!(replay.work.original_forward_solves,4);
            assert_eq!(replay.work.recomputed_forward_solves,4);
            assert_eq!(objective.radiation[0].region,"surface");
            assert_eq!(objective.radiation[0].reservoir_values_k.len(),values.len());
            let analytic=[components(&objective.radiation[0]),components(&rows[0].radiation[0]),
                components(&envelopes[0].radiation[0])];
            // Value-only outputs use the SAME primal and peak reduction without
            // building a tangent during the finite-difference probes.
            let evaluate=|values:&[f64],emissivity|{
                let (objective,rows,peaks)=evaluate_scan_values_with_responses(cx,make(values,emissivity),
                    &times,&observations(),sampling,config(),
                    ScanValueBudget{max_steps:4,max_observations:2,max_schedule_controls:7},
                    &points,&peaks,ScanValueResponseBudget{max_responses:2,max_peak_samples:8}).unwrap();
                [objective.value,rows[0].value,peaks[0].value]
            };
            assert_eq!(evaluate(&values,0.2),[objective.value,rows[0].value,envelopes[0].value]);
            for coordinate in 0..=values.len(){
                let (mut plus,mut minus)=(values.clone(),values.clone());
                let (mut ep,mut em)=(0.2,0.2);let delta=if coordinate==values.len(){1e-5}else{1e-3};
                if coordinate==values.len(){ep+=delta;em-=delta;}
                else{plus[coordinate]+=delta;minus[coordinate]-=delta;}
                let (a,b)=(evaluate(&plus,ep),evaluate(&minus,em));
                for row in 0..3{close(analytic[row][coordinate],(a[row]-b[row])/(2.*delta));}
            }
            if kind==3&&sampling==ScanObjectiveSampling::EventLeft{
                // The original interval-ending side excludes every future hold.
                assert_eq!(objective.radiation[0].reservoir_values_k[2],0.);
                assert_eq!(rows[0].radiation[0].reservoir_values_k[1..],[0.,0.]);
                assert_eq!(envelopes[0].radiation[0].reservoir_values_k[2],0.);
            }
        }
    });
}

#[test]
fn g0_radiative_controls_charge_objective_and_response_scalar_allowances(){
    let m=mesh();let curve=chart();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx|{
        let hist=history(cx,&m);let old=[10.,20.,30.,40.];
        let p=[0.2,0.25,0.3,1.4,1.3,1.2,3.];let times=[0.001,0.002,0.004,0.006];
        let (points,peaks)=requests();
        let make=||radiative_run(cx,&hist,&material,&curve,&old,vec![heater(p)],
            sourced_radiation(reservoir(3,&[291.,295.,298.]),0.2));
        let call=|budget,rows|evaluate_scan_objective_with_peaks(cx,make(),&times,&observations(),
            ScanObjectiveSampling::EventLeft,config(),linear(),budget,&points,&peaks,rows,8);
        // Seven source coordinates + three reservoir values + one emissivity.
        let exact=ScanObjectiveBudget{max_schedule_controls:11,..objective_budget()};
        let rows=ScanResponseBudget{max_gradient_entries:30,..response_budget()};
        assert!(call(ScanObjectiveBudget{max_schedule_controls:10,..exact},rows).is_err());
        assert!(call(exact,ScanResponseBudget{max_gradient_entries:29,..rows}).is_err());
        assert!(call(exact,rows).is_ok());
        assert!(evaluate_scan_objective_checkpointed_with_peaks(cx,make(),&times,&observations(),
            ScanObjectiveSampling::EventLeft,config(),linear(),exact,&points,&peaks,
            ScanResponseBudget{max_gradient_entries:29,..rows},8,checkpoints()).is_err());
        let ordinary=evaluate_scan_objective(cx,run(cx,&hist,&material,&curve,&old,vec![heater(p)],false),
            &times,&observations(),ScanObjectiveSampling::EventLeft,config(),linear(),
            ScanObjectiveBudget{max_schedule_controls:7,..objective_budget()}).unwrap();
        assert!(ordinary.radiation.is_empty());
    });
}

#[test]
fn g0_value_only_radiation_keeps_chirps_and_valid_constitutive_corners(){
    let m=mesh();let curve=chart();let material=ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx|{
        let hist=history(cx,&m);let p=[0.2,0.25,0.3,1.4,1.3,1.2,3.];
        let chirp=TimeSignal::Chirp(LinearChirp::new(0.,0.006,1.,2.,0.,
            QtyAny::new(293.,TEMPERATURE_DIMS),QtyAny::new(1.,TEMPERATURE_DIMS)).unwrap());
        let make=||radiative_run(cx,&hist,&material,&curve,&[20.;4],vec![heater(p)],
            sourced_radiation(chirp.clone(),0.2));
        let observed=[observation(0,0,ScanObservable::Temperature,300.,1.)];
        let values=ScanValueBudget{max_steps:1,max_observations:1,max_schedule_controls:7};
        assert!(evaluate_scan_values(cx,make(),&[0.001],&observed,
            ScanObjectiveSampling::RightEndpoint,config(),values).is_ok());
        assert!(evaluate_scan_objective(cx,make(),&[0.001],&observed,
            ScanObjectiveSampling::RightEndpoint,config(),linear(),objective_budget()).is_err());
        let mut off=p;off[6]=0.;
        let corner=||radiative_run(cx,&hist,&material,&curve,&[100.;4],vec![heater(off)],
            sourced_radiation(constant(310.,TEMPERATURE_DIMS),0.2));
        let result=evaluate_scan_values(cx,corner(),&[0.001],&observed,
            ScanObjectiveSampling::RightEndpoint,config(),values).unwrap();
        assert_eq!(result.final_step.solution().temperature,[310.;4]);
        assert!(evaluate_scan_objective(cx,corner(),&[0.001],&observed,
            ScanObjectiveSampling::RightEndpoint,config(),linear(),objective_budget()).is_err());
    });
}

//! CSV-driven inverse fitting on an explicitly synthetic reference specimen.
//! The library supports caller-owned meshes/materials; this executable fixes
//! one model so external measurement parsing never silently chooses physics.
use std::io::Read;
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerCoordinates, ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::constraints::ConstrainedScanStudy;
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryRun, EnthalpyHistoryMaterials, EnthalpyHistoryBudget};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityError,
    objective::{ScanObservable, ScanObjectiveSampling, ScanObjectiveBudget}};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::spatial::{ScanReadout, ScanSpatialObservation};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::loss::{ScanObservationWithLoss, ScanResidualLoss};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::temporal::ScanTemporalObservation;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::checkpoint::ScanCheckpointBudget;
use fs_conduction::bc::history::march::moving::{MovingHeatSource, MovingHeatBudget};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh, ConductivityModel, ScalarField, ThermalBc,
    ThermalBoundaryBuilder, LinearConfig, COMPONENT_POWER_DIMS};
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot};
use fs_qty::{Dims, QtyAny};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{TimeSignal, Interp};
use fs_solver::NewtonKrylovConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const MAX_BYTES: usize = 128 * 1024;
const MAX_ROWS: usize = 256;
const MAX_STEPS: usize = 4096;
const MAX_SAMPLES: usize = 65_536;
const HEADER: &str = "time_s,observable,target,scale,readout,selection";
const USAGE: &str = "spatial_measurement_fit measurements.csv [initial_first_W initial_second_W] [--huber positive_threshold] [--substeps 1..128]; time_s may be start..end for a time mean";

struct Options {
    path: String,
    initial: [f64;2],
    loss: ScanResidualLoss,
    substeps: usize,
}
fn options(args: &[String]) -> Result<Options> {
    let mut positional=Vec::new();let mut threshold=None;let mut substeps=None;
    let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--huber" if threshold.is_none() => {
                i+=1;threshold=Some(finite(args.get(i).ok_or(USAGE)?)?);
            }
            "--substeps" if substeps.is_none() => {
                i+=1;substeps=Some(args.get(i).ok_or(USAGE)?.parse::<usize>()?);
            }
            value if !value.starts_with("--") => positional.push(value),
            _ => return Err(USAGE.into()),
        }
        i+=1;
    }
    if !matches!(positional.len(),1|3) {return Err(USAGE.into());}
    let loss=threshold.map_or(ScanResidualLoss::Quadratic,|threshold|ScanResidualLoss::Huber {threshold});
    loss.validate()?;
    let initial=if positional.len()==3 {[finite(positional[1])?,finite(positional[2])?]} else {[1.,1.]};
    if initial.iter().any(|p|*p<0.) {return Err("initial watts must be nonnegative".into());}
    let substeps=substeps.unwrap_or(1);
    if !(1..=128).contains(&substeps) {return Err("substeps must be in 1..128".into());}
    Ok(Options {path:positional[0].to_string(),initial,loss,substeps})
}

struct Measurements<'m> {
    endpoints: Vec<f64>,
    rows: Vec<ScanSpatialObservation<'m>>,
    starts: Vec<Option<usize>>,
    contributions: usize,
}
fn finite(text: &str) -> Result<f64> {
    let value: f64=text.parse()?;
    if !value.is_finite() { return Err("numeric input must be finite".into()); }
    Ok(value)
}
fn indices(text: &str, maximum: usize) -> Result<Vec<usize>> {
    let mut out=Vec::new();
    for part in text.split(';') {
        if out.len() >= maximum { return Err("readout selection exceeds this mesh's size".into()); }
        out.push(part.trim().parse()?);
    }
    Ok(out)
}
fn readout<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh, kind: &str, selection: &str)
    -> Result<ScanReadout<'m>>
{
    Ok(match kind {
        "vertex" => ScanReadout::vertex(cx,mesh,selection.parse()?)?,
        "tet" => {
            let p:Vec<_>=selection.split(':').map(str::trim).collect();
            if p.len()!=5 { return Err("tet selection is element:b0:b1:b2:b3".into()); }
            ScanReadout::in_tetrahedron(cx,mesh,p[0].parse()?,
                [finite(p[1])?,finite(p[2])?,finite(p[3])?,finite(p[4])?])?
        }
        "volume" => ScanReadout::volume_mean(cx,mesh,&indices(selection,mesh.element_count())?,mesh.element_count())?,
        "weighted" => {
            let mut weights=Vec::new();
            for p in selection.split(';') {
                if weights.len()>=mesh.vertex_count() { return Err("sensor row exceeds this mesh's size".into()); }
                let (v,w)=p.trim().split_once(':').ok_or("weighted selection is vertex:weight;vertex:weight")?;
                weights.push((v.trim().parse()?,finite(w.trim())?));
            }
            ScanReadout::weighted_mean(cx,mesh,&weights,mesh.vertex_count())?
        }
        _ => return Err("readout must be vertex, tet, volume or weighted".into()),
    })
}
fn parse<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh, text: &str, substeps: usize) -> Result<Measurements<'m>> {
    if text.len()>MAX_BYTES { return Err("measurement file exceeds 128 KiB".into()); }
    if !(1..=128).contains(&substeps) {return Err("substeps must be in 1..128".into());}
    let mut lines=text.lines();
    if lines.next().map(str::trim)!=Some(HEADER) { return Err(format!("required CSV header: {HEADER}").into()); }
    let mut rows=Vec::new(); let mut times=Vec::new();
    for (line_index,line) in lines.enumerate() {
        cx.checkpoint().map_err(|_| "measurement parsing cancelled")?;
        let line=line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let row=(|| -> Result<_> {
            if rows.len()>=MAX_ROWS { return Err("measurement count exceeds 256".into()); }
            let fields:Vec<_>=line.split(',').map(str::trim).collect();
            if fields.len()!=6 || fields.iter().any(|p|p.is_empty()) {
                return Err("each unquoted numeric measurement needs exactly six fields".into());
            }
            let (start,time)=match fields[0].split_once("..") {
                Some((start,end)) => (Some(finite(start.trim())?),finite(end.trim())?),
                None => (None,finite(fields[0])?),
            };
            if time<=0. || time>2. { return Err("reference specimen measurements require 0 < time_s <= 2".into()); }
            if start.is_some_and(|start|start<=0. || start>=time) {
                return Err("time mean needs 0 < start < end <= 2 seconds; initial-state windows are unsupported".into());
            }
            let observable=match fields[1] {
                "temperature_k" => ScanObservable::Temperature,
                "specific_enthalpy_j_kg" => ScanObservable::SpecificEnthalpy,
                "liquid_mass_fraction" => ScanObservable::LiquidMassFraction,
                _ => return Err("observable must explicitly name temperature_k, specific_enthalpy_j_kg or liquid_mass_fraction".into()),
            };
            let target=finite(fields[2])?; let scale=finite(fields[3])?;
            if scale<=0. { return Err("measurement scale must be positive in observable units".into()); }
            if observable==ScanObservable::LiquidMassFraction && !(0.0..=1.0).contains(&target) {
                return Err("liquid mass fraction is a fraction in [0,1], not a percentage".into());
            }
            Ok(((start,time),ScanSpatialObservation {endpoint:0,readout:readout(cx,mesh,fields[4],fields[5])?,observable,target,scale}))
        })().map_err(|e|format!("measurement line {}: {e}",line_index+2))?;
        times.push(row.0); rows.push(row.1);
    }
    if rows.is_empty() { return Err("at least one measurement is required".into()); }
    // Preserve exact representable times and every source switch. Inserting
    // measurements cannot make the grid step over the one-second power change.
    let mut endpoints=Vec::new();
    for &(start,end) in &times {endpoints.extend(start);endpoints.push(end);}
    endpoints.extend([1.0,2.0]);
    endpoints.sort_by(f64::total_cmp); endpoints.dedup();
    let count=endpoints.len().checked_mul(substeps).filter(|&n|n<=MAX_STEPS)
        .ok_or("measurement time grid exceeds 4096 physical endpoints")?;
    let mut refined=Vec::new();refined.try_reserve_exact(count)?;
    let mut previous=0.0;
    for &end in &endpoints {
        for part in 1..=substeps {
            cx.checkpoint().map_err(|_|"measurement grid construction cancelled")?;
            let time=if part==substeps {end} else {previous+(end-previous)*(part as f64/substeps as f64)};
            let old=refined.last().copied().unwrap_or(0.0);
            if time<=old || (part<substeps && time>=end) {
                return Err("requested substeps are not representable between exact measurement/event times".into());
            }
            refined.push(time);
        }
        previous=end;
    }
    endpoints=refined;
    let mut starts=Vec::with_capacity(rows.len());let mut contributions=0usize;
    for (row,(start,time)) in rows.iter_mut().zip(times) {
        row.endpoint=endpoints.binary_search_by(|t|t.total_cmp(&time)).map_err(|_| "measurement endpoint missing")?;
        let start=start.map(|start|endpoints.binary_search_by(|t|t.total_cmp(&start))
            .map_err(|_|"measurement window start missing")).transpose()?;
        contributions=contributions.checked_add(start.map_or(1,|first|row.endpoint-first+1))
            .filter(|&n|n<=MAX_SAMPLES).ok_or("time-window contributions exceed 65536 samples")?;
        starts.push(start);
    }
    cx.checkpoint().map_err(|_| "measurement parsing cancelled")?;
    Ok(Measurements {endpoints,rows,starts,contributions})
}

struct Fit {
    powers: Vec<f64>,
    measurements: Vec<(f64, ScanObservable, f64, f64, f64)>,
    exposure_starts_s: Vec<Option<f64>>,
    objective: f64,
    evaluations: usize,
    converged: bool,
    violation: f64,
}
fn fit(cx: &Cx<'_>, text: &str, initial: [f64;2], loss: ScanResidualLoss, substeps: usize) -> Result<Fit> {
    loss.validate()?;
    if initial.iter().any(|p|!p.is_finite() || *p<0.) { return Err("initial watts must be finite and nonnegative".into()); }
    let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let data=parse(cx,&mesh,text,substeps)?; let n=data.endpoints.len();
    let material=ConductivityModel::isotropic_declared(2.)?;
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([76;32]),
        [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)].into_iter()
            .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),vec![],0.)?;
    let length=Dims([1,0,0,0,0,0]);
    let source=MovingHeatSource::new([
        TimeSignal::Table {times:vec![0.,1.,2.],values:vec![0.25,0.75,0.75],dims:length,interp:Interp::Hold},
        TimeSignal::Constant(QtyAny::new(0.5,length)),TimeSignal::Constant(QtyAny::new(0.5,length))],
        [1.4,1.2,1.1],TimeSignal::Table {times:vec![0.,1.,2.],values:vec![initial[0],initial[1],0.],
            dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},0.,
        MovingHeatBudget {max_vertices:8,max_elements:6,min_active_vertices:1,power_relative_tolerance:1e-11})?;
    let factory=|cx: &Cx<'_>,sources: Vec<MovingHeatSource>|
        -> std::result::Result<ScanningEnthalpyRun<'_, '_, '_>,ScanSensitivityError> {
        let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&boundary,material:&material,
            element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
            EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}],
                of_element:&[0;6],phase_conductivity:None},&[150.;8],
            EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:8,max_elements:6},max_steps:n})?;
        Ok(ScanningEnthalpyRun::new(cx,base,sources,1,None)?)
    };
    let segment=n.min(8);
    let observations=data.rows.iter().cloned().zip(&data.starts).map(|(row,&start)| {
        let temporal=match start {
            None=>ScanTemporalObservation::instantaneous(row),
            Some(first)=>ScanTemporalObservation::trapezoidal_mean(cx,row,&data.endpoints,first,MAX_STEPS)?,
        };
        ScanObservationWithLoss::new(temporal,loss)
    })
        .collect::<std::result::Result<Vec<_>,_>>()?;
    let objective=ScanPowerObjective {endpoints_s:data.endpoints.clone(),observations,
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
            linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
            max_newton_iterations:64,energy_tolerance_j:1e-8},
        transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
        budget:ScanObjectiveBudget {max_steps:n,max_observations:data.contributions,max_retained_vertices:8*segment,
            max_retained_source_vertices:8*segment,max_schedule_controls:16,max_retained_radiation_entries:0},
        source_reference_time_s:0.};
    let controls=vec![ScanPowerControl {source:0,value:0,reference_w:1.},ScanPowerControl {source:0,value:1,reference_w:1.}];
    let mut study=ScanPowerStudy::new_joint_with_power_coordinates(cx,vec![source],controls,vec![],
        ScanPowerCoordinates::NonnegativeLinear,objective,factory,127)?;
    study.set_checkpointing(cx,ScanCheckpointBudget {segment_steps:segment,
        max_checkpoint_vertices:8*(n.div_ceil(segment)+1),max_forward_solves:2*n})?;
    let mut study=ConstrainedScanStudy::new(cx,study,vec![],vec![],8)?;
    let report=study.advance_optimizer(cx,1e-8,100)?;
    let spent=study.evaluations(); let exported=study.sources(cx)?;
    if study.evaluations()!=spent { return Err("source export spent a physical evaluation".into()); }
    let powers=exported[0].absorbed_power_values(3)?;
    if powers[2]!=0. { return Err("fitting changed the declared shutdown".into()); }
    study.extend_evaluation_budget(128)?;
    let (replay,limits)=study.replay_with_constraints(cx)?;
    if (report.solution.f-replay.physical.value).abs()>1e-10*(1.+report.solution.f.abs()) {
        return Err("replayed sources changed the accepted fitting objective".into());
    }
    let measurements=data.rows.iter().zip(&replay.physical.predictions)
        .zip(&replay.physical.standardized_residuals).map(|((row,&prediction),&residual)|
            (data.endpoints[row.endpoint],row.observable,row.target,prediction,residual)).collect();
    let exposure_starts_s=data.starts.iter().map(|start|start.map(|i|data.endpoints[i])).collect();
    Ok(Fit {powers,measurements,exposure_starts_s,objective:replay.physical.value,evaluations:study.evaluations(),
        converged:report.solution.converged,violation:limits.max_scaled_violation()})
}
fn main() -> Result<()> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    let options=options(&args)?;
    let mut text=String::new();
    std::fs::File::open(&options.path)?.take(MAX_BYTES as u64+1).read_to_string(&mut text)?;
    let gate=CancelGate::new(); let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result=pool.scope(|arena|fit(&Cx::new(&gate,arena,
        StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic),&text,options.initial,options.loss,options.substeps))?;
    println!("# Synthetic unit-cube reference model; not automatic experimental calibration.");
    println!("# residual_loss={:?}; residuals below are raw, not clipped",options.loss);
    println!("# substeps_per_event_interval={}; time means use normalized trapezoidal readouts",options.substeps);
    println!("# objective={:.17e} converged={} max_scaled_violation={:.17e} evaluations={}",
        result.objective,result.converged,result.violation,result.evaluations);
    println!("time_s,absorbed_power_w");
    for (time,power) in [0.,1.,2.].into_iter().zip(&result.powers) {println!("{time:.17e},{power:.17e}");}
    for (i,(time,quantity,target,prediction,residual)) in result.measurements.iter().enumerate() {
        let (term,slope)=options.loss.value_derivative(*residual)?;
        if let Some(start)=result.exposure_starts_s[i] {
            println!("# measurement={i} interpretation=time_mean exposure_start_s={start:.17e} exposure_end_s={time:.17e}");
        }
        println!("# measurement={i} time_s={time:.17e} observable={quantity:?} target={target:.17e} prediction={prediction:.17e} standardized_residual={residual:.17e} loss_term={term:.17e} dloss_dstandardized_residual={slope:.17e}");
    }
    if !result.converged || result.violation>1e-8 { return Err("fit did not pass convergence and nonnegative-domain checks".into()); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const DATA: &str=include_str!("spatial_measurements.csv");
    fn with_cx(f: impl FnOnce(&Cx<'_>)) {
        let gate=CancelGate::new(); let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena|f(&Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic)));
    }
    #[test]
    fn g1_csv_measurements_fit_their_actual_spatial_mean_from_off() {
        with_cx(|cx| {
            let solved=fit(cx,DATA,[0.,0.],ScanResidualLoss::Quadratic,1).unwrap();
            assert!(solved.converged && solved.violation<=1e-8 && solved.objective<1e-12);
            assert!((solved.powers[0]-2.).abs()<1e-5 && (solved.powers[1]-5.).abs()<1e-5);
            assert_eq!(solved.powers[2],0.);
            assert_eq!(solved.measurements.len(),3);
        });
    }
    #[test]
    fn g0_csv_rejects_bad_units_readouts_and_percentage_values_before_solving() {
        with_cx(|cx| {
            let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
            let mesh=ConductionMesh::new(complex,positions).unwrap();
            for line in ["1,temperature,310,1,vertex,0","1,liquid_mass_fraction,60,1,vertex,0",
                "1,specific_enthalpy_j_kg,152,0,vertex,0","NaN,specific_enthalpy_j_kg,152,1,vertex,0",
                "3,specific_enthalpy_j_kg,152,1,vertex,0","1,specific_enthalpy_j_kg,152,1,tet,0:-1:1:1:1",
                "1,specific_enthalpy_j_kg,152,1,volume,1;0","1,specific_enthalpy_j_kg,152,1,weighted,0:0;7:0"] {
                let error=parse(cx,&mesh,&format!("{HEADER}\n{line}\n"),1).err().unwrap().to_string();
                assert!(error.contains("line 2"),"{error}");
            }
            assert!(parse(cx,&mesh,HEADER,1).is_err());
            assert!(parse(cx,&mesh,&format!("{HEADER}\n{}","1,specific_enthalpy_j_kg,152,1,vertex,0\n".repeat(MAX_ROWS+1)),1).is_err());
        });
    }
    #[test]
    fn g3_csv_preserves_input_order_and_inserts_each_source_switch() {
        with_cx(|cx| {
            let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
            let mesh=ConductionMesh::new(complex,positions).unwrap();
            let text=format!("{HEADER}\n1.5,temperature_k,310,1,tet,0:1:2:3:4\n0.5,specific_enthalpy_j_kg,151,10,weighted,0:1;7:1\n");
            let data=parse(cx,&mesh,&text,1).unwrap();
            assert_eq!(data.endpoints,vec![0.5,1.,1.5,2.]);
            assert_eq!(data.rows.iter().map(|r|r.endpoint).collect::<Vec<_>>(),vec![2,0]);
            assert_eq!(data.rows[0].readout.evaluate(cx,&[310.;8]).unwrap(),310.);
        });
    }

    #[test]
    fn g1_csv_exposure_fit_recovers_original_heating_not_an_endpoint_surrogate() {
        const EXPOSURES: &str=include_str!("temporal_measurements.csv");
        with_cx(|cx| {
            for substeps in [1,4] {
                let solved=fit(cx,EXPOSURES,[0.,0.],ScanResidualLoss::Quadratic,substeps).unwrap();
                assert!(solved.converged && solved.violation<=1e-8 && solved.objective<1e-12);
                assert!((solved.powers[0]-2.).abs()<1e-4 && (solved.powers[1]-5.).abs()<1e-4);
                assert_eq!(solved.powers[2],0.);
                assert_eq!(solved.exposure_starts_s,vec![None,Some(1.),Some(1.)]);
                let wrong=fit(cx,&EXPOSURES.replace("1..2","2"),[0.,0.],ScanResidualLoss::Quadratic,substeps).unwrap();
                assert!(wrong.converged && wrong.objective<1e-12);
                assert!((wrong.powers[0]-2.).abs()<1e-4 && (wrong.powers[1]-2.5).abs()<1e-4);
                assert!(wrong.exposure_starts_s.iter().all(Option::is_none));
            }
        });
    }

    #[test]
    fn g3_csv_window_grid_retains_exact_events_samples_and_input_order() {
        with_cx(|cx| {
            let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
            let mesh=ConductionMesh::new(complex,positions).unwrap();
            let data=parse(cx,&mesh,include_str!("temporal_measurements.csv"),4).unwrap();
            assert_eq!(data.endpoints,vec![0.25,0.5,0.75,1.,1.25,1.5,1.75,2.]);
            assert_eq!(data.starts,vec![None,Some(3),Some(3)]);
            assert_eq!(data.rows.iter().map(|r|r.endpoint).collect::<Vec<_>>(),vec![3,7,7]);
            assert_eq!(data.contributions,11);
            let text=format!("{HEADER}\n0.5..1.5,temperature_k,310,1,vertex,0\n0.75,specific_enthalpy_j_kg,151.5,10,vertex,0\n");
            let data=parse(cx,&mesh,&text,1).unwrap();
            assert_eq!(data.endpoints,vec![0.5,0.75,1.,1.5,2.]);
            assert_eq!(data.starts,vec![Some(0),None]);
            assert_eq!(data.rows.iter().map(|r|r.endpoint).collect::<Vec<_>>(),vec![3,1]);
            assert_eq!(data.contributions,5);
        });
    }

    #[test]
    fn g0_csv_exposures_refuse_invalid_windows_excess_work_and_unrepresentable_substeps() {
        let opts=|a:&[&str]|options(&a.iter().map(|v|v.to_string()).collect::<Vec<_>>());
        assert_eq!(opts(&["data.csv"]).unwrap().substeps,1);
        let options=opts(&["--substeps","4","data.csv","0","0","--huber","1"]).unwrap();
        assert_eq!(options.substeps,4);assert_eq!(options.initial,[0.,0.]);
        assert_eq!(options.loss,ScanResidualLoss::Huber {threshold:1.});
        for args in [vec!["data.csv","--substeps"],vec!["data.csv","--substeps","0"],
            vec!["data.csv","--substeps","129"],vec!["data.csv","--substeps","2","--substeps","3"]] {
            assert!(opts(&args).is_err());
        }
        with_cx(|cx| {
            let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
            let mesh=ConductionMesh::new(complex,positions).unwrap();
            for time in ["0..1","1..1","1.5..1","1..3","NaN..2","1..NaN","1..2..3"] {
                let text=format!("{HEADER}\n{time},temperature_k,310,1,vertex,0\n");
                let error=parse(cx,&mesh,&text,1).err().unwrap().to_string();
                assert!(error.contains("line 2"),"{error}");
            }
            let mut many=format!("{HEADER}\n");
            for i in 1..=33 {many.push_str(&format!("{},temperature_k,310,1,vertex,0\n",f64::from(i)/64.));}
            assert!(parse(cx,&mesh,&many,128).is_err());
            let narrow=format!("{HEADER}\n1..{},temperature_k,310,1,vertex,0\n",1.0_f64.next_up());
            assert!(parse(cx,&mesh,&narrow,2).is_err());
            let mut overlaps=format!("{HEADER}\n");
            for i in 0..256 {overlaps.push_str(&format!("{}..2,temperature_k,310,1,vertex,0\n",f64::from(i+1)/512.));}
            assert!(parse(cx,&mesh,&overlaps,4).is_err());
        });
    }
}

#[cfg(test)]
mod robust_tests {
    use super::*;

    #[test]
    fn g0_huber_options_are_explicit_and_invalid_thresholds_refuse_before_io() {
        let parse=|items:&[&str]|options(&items.iter().map(|s|s.to_string()).collect::<Vec<_>>());
        assert_eq!(parse(&["data.csv"]).unwrap().loss,ScanResidualLoss::Quadratic);
        let chosen=parse(&["data.csv","0","0","--huber","1.5"]).unwrap();
        assert_eq!(chosen.initial,[0.,0.]);
        assert_eq!(chosen.loss,ScanResidualLoss::Huber {threshold:1.5});
        assert_eq!(parse(&["data.csv","--huber","1"]).unwrap().initial,[1.,1.]);
        for args in [vec![],vec!["--huber","1"],vec!["data.csv","--huber"],
            vec!["data.csv","--huber","0"],vec!["data.csv","--huber","NaN"],
            vec!["data.csv","--huber","-1"],vec!["data.csv","0"],
            vec!["data.csv","--huber","1","--huber","2"]] {
            assert!(parse(&args).is_err());
        }
    }

    #[test]
    fn g1_contaminated_mean_fit_limits_influence_without_hiding_raw_residuals() {
        let mut data=format!("{HEADER}\n");
        for _ in 0..4 {data.push_str("1,specific_enthalpy_j_kg,152,1,volume,0;1;2;3;4;5\n");}
        for _ in 0..4 {data.push_str("2,specific_enthalpy_j_kg,157,1,volume,0;1;2;3;4;5\n");}
        data.push_str("2,specific_enthalpy_j_kg,197,1,volume,0;1;2;3;4;5\n");
        let gate=CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx=Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
                Budget::INFINITE,ExecMode::Deterministic);
            for (loss,expected) in [(ScanResidualLoss::Quadratic,[2.,13.]),
                (ScanResidualLoss::Huber {threshold:1.},[2.,5.25])] {
                let result=fit(&cx,&data,[0.,0.],loss,1).unwrap();
                assert!(result.converged && result.violation<=1e-8);
                for (power,target) in result.powers.iter().zip(expected) {assert!((power-target).abs()<1e-4);}
                assert_eq!(result.powers[2],0.);
                assert_eq!(result.measurements.len(),9);
                assert!(result.measurements[8].4.abs()>30.);
                let manual:f64=result.measurements.iter().map(|m|loss.value_derivative(m.4).unwrap().0).sum();
                assert!((result.objective-manual).abs()<1e-10*(1.+manual));
            }
        });
    }
}

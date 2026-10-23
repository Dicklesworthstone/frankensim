//! CSV-driven inverse fitting on an explicit caller-supplied or reference specimen.
//! Numerical material/mesh/solver owners are shared with the library; input
//! parsing never infers physical conditions from the measurement values.
use std::io::Read;
use scan_power_fit::specimen::{ScanSpecimen, SpecimenBudget};
use fs_ascent::scan_design::{ScanPowerControl, ScanPowerCoordinates, ScanPowerObjective, ScanPowerStudy};
use fs_ascent::scan_design::constraints::ConstrainedScanStudy;
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
use fs_conduction::bc::history::march::moving::MovingHeatSource;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{ConductionMesh, ScalarField, ThermalBoundaryBuilder, LinearConfig};
use fs_exec::{Cx, Budget, CancelGate, ExecMode, StreamKey};
use fs_scenario::scenario::Environment;
use fs_solver::NewtonKrylovConfig;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const MAX_BYTES: usize = 128 * 1024;
const MAX_ROWS: usize = 256;
const MAX_STEPS: usize = 4096;
const MAX_SAMPLES: usize = 65_536;
const SPECIMEN_BUDGET: SpecimenBudget = SpecimenBudget {
    max_bytes: 2 << 20, max_vertices: 4096, max_elements: 16_384, max_phase_knots: 256,
};
const HEADER: &str = "time_s,observable,target,scale,readout,selection";
const USAGE: &str = "spatial_measurement_fit measurements.csv [initial_first_W initial_second_W] [--huber positive_threshold] [--substeps 1..128] [--specimen FILE]; time_s may be start..end for a mean or lag:start:end:tau:initial for a sensor reading";

struct Options {
    path: String,
    initial: [f64;2],
    loss: ScanResidualLoss,
    substeps: usize,
    specimen: Option<String>,
}
fn options(args: &[String]) -> Result<Options> {
    let mut positional=Vec::new();let mut threshold=None;let mut substeps=None;let mut specimen=None;
    let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--huber" if threshold.is_none() => {
                i+=1;threshold=Some(finite(args.get(i).ok_or(USAGE)?)?);
            }
            "--substeps" if substeps.is_none() => {
                i+=1;substeps=Some(args.get(i).ok_or(USAGE)?.parse::<usize>()?);
            }
            "--specimen" if specimen.is_none() => {
                i+=1;
                let path=args.get(i).filter(|p| !p.is_empty() && !p.starts_with("--")).ok_or(USAGE)?;
                specimen=Some(path.clone());
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
    Ok(Options {path:positional[0].to_string(),initial,loss,substeps,specimen})
}

struct Measurements<'m> {
    endpoints: Vec<f64>,
    rows: Vec<ScanSpatialObservation<'m>>,
    starts: Vec<Option<usize>>,
    sensor_lags: Vec<Option<SensorLag>>,
    contributions: usize,
}
#[derive(Debug, Clone, Copy, PartialEq)]
struct SensorLag { time_constant_s: f64, initial_reading: f64 }

fn measurement_time(text: &str) -> Result<(Option<f64>, f64, Option<SensorLag>)> {
    if let Some(spec) = text.strip_prefix("lag:") {
        let values: Vec<_> = spec.split(':').take(5).map(str::trim).collect();
        if values.len() != 4 { return Err("lag syntax is lag:start_s:end_s:tau_s:initial_reading".into()); }
        let start = finite(values[0])?;
        let end = finite(values[1])?;
        let time_constant_s = finite(values[2])?;
        let initial_reading = finite(values[3])?;
        if time_constant_s <= 0.0 { return Err("sensor time constant must be positive seconds".into()); }
        return Ok((Some(start), end, Some(SensorLag { time_constant_s, initial_reading })));
    }
    Ok(match text.split_once("..") {
        Some((start,end)) => (Some(finite(start.trim())?),finite(end.trim())?,None),
        None => (None,finite(text)?,None),
    })
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
#[cfg(test)]
fn parse<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh, text: &str, substeps: usize) -> Result<Measurements<'m>> {
    parse_for_clock(cx, mesh, text, substeps, [1.0, 2.0])
}
fn parse_for_clock<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh, text: &str,
    substeps: usize, source_times: [f64; 2]) -> Result<Measurements<'m>> {
    if source_times.iter().any(|t| !t.is_finite()) || source_times[0] <= 0.0
        || source_times[1] <= source_times[0] {
        return Err("source clock requires finite 0 < switch < shutdown".into());
    }
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
            let (start,time,sensor_lag)=measurement_time(fields[0])?;
            if time<=0. || time>source_times[1] { return Err("measurements require 0 < time_s <= the specimen shutdown".into()); }
            if start.is_some_and(|start|start<=0. || start>=time) {
                return Err("temporal readout needs 0 < start < end <= shutdown; initial-state windows are unsupported".into());
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
            if observable==ScanObservable::LiquidMassFraction
                && sensor_lag.is_some_and(|lag| !(0.0..=1.0).contains(&lag.initial_reading)) {
                return Err("initial sensor fraction must be in [0,1], not percent".into());
            }
            Ok(((start,time,sensor_lag),ScanSpatialObservation {endpoint:0,readout:readout(cx,mesh,fields[4],fields[5])?,observable,target,scale}))
        })().map_err(|e|format!("measurement line {}: {e}",line_index+2))?;
        times.push(row.0); rows.push(row.1);
    }
    if rows.is_empty() { return Err("at least one measurement is required".into()); }
    // Preserve exact representable times and every source switch. Inserting
    // measurements cannot make the grid step over the declared power change.
    let mut endpoints=Vec::new();
    for &(start,end,_) in &times {endpoints.extend(start);endpoints.push(end);}
    endpoints.extend(source_times);
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
    let mut sensor_lags=Vec::with_capacity(rows.len());
    for (row,(start,time,sensor_lag)) in rows.iter_mut().zip(times) {
        row.endpoint=endpoints.binary_search_by(|t|t.total_cmp(&time)).map_err(|_| "measurement endpoint missing")?;
        let start=start.map(|start|endpoints.binary_search_by(|t|t.total_cmp(&start))
            .map_err(|_|"measurement window start missing")).transpose()?;
        contributions=contributions.checked_add(start.map_or(1,|first|row.endpoint-first+1))
            .filter(|&n|n<=MAX_SAMPLES).ok_or("time-window contributions exceed 65536 samples")?;
        starts.push(start);
        sensor_lags.push(sensor_lag);
    }
    cx.checkpoint().map_err(|_| "measurement parsing cancelled")?;
    Ok(Measurements {endpoints,rows,starts,sensor_lags,contributions})
}

struct Fit {
    powers: Vec<f64>,
    source_times_s: [f64; 3],
    mesh_vertices: usize,
    mesh_elements: usize,
    reference_mass_kg: f64,
    measurements: Vec<(f64, ScanObservable, f64, f64, f64)>,
    exposure_starts_s: Vec<Option<f64>>,
    sensor_lags: Vec<Option<SensorLag>>,
    objective: f64,
    evaluations: usize,
    converged: bool,
    violation: f64,
}
#[cfg(test)]
fn fit(cx: &Cx<'_>, text: &str, initial: [f64;2], loss: ScanResidualLoss, substeps: usize) -> Result<Fit> {
    let model=ScanSpecimen::synthetic_reference(cx)?;
    fit_on_specimen(cx, &model, text, initial, loss, substeps)
}
fn fit_on_specimen(cx: &Cx<'_>, model: &ScanSpecimen, text: &str,
    initial: [f64;2], loss: ScanResidualLoss, substeps: usize) -> Result<Fit> {
    loss.validate()?;
    if initial.iter().any(|p|!p.is_finite() || *p<0.) { return Err("initial watts must be finite and nonnegative".into()); }
    let mesh=&model.mesh;
    let vertices=mesh.vertex_count(); let elements=mesh.element_count();
    if vertices>SPECIMEN_BUDGET.max_vertices || elements>SPECIMEN_BUDGET.max_elements {
        return Err("physical specimen exceeds the executable mesh allowance".into());
    }
    let data=parse_for_clock(cx,mesh,text,substeps,[model.switch_time_s,model.end_time_s])?;
    let n=data.endpoints.len();
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(mesh)
        .remainder("surface",model.exterior.clone())?,&Environment::earth_lab(),vec![],0.)?;
    let source=model.source(initial)?;
    let initial_h=vec![model.initial_specific_enthalpy_j_kg; vertices];
    let element_materials=vec![0; elements];
    let reference_mass_kg=(0..elements).map(|i| mesh.element_volume(i)*model.reference_density_kg_m3).sum::<f64>();
    if !reference_mass_kg.is_finite() || reference_mass_kg<=0.0 {
        return Err("specimen reference mass is not representable".into());
    }
    let factory=|cx: &Cx<'_>,sources: Vec<MovingHeatSource>|
        -> std::result::Result<ScanningEnthalpyRun<'_, '_, '_>,ScanSensitivityError> {
        let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&boundary,material:&model.material,
            element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
            EnthalpyHistoryMaterials {table:&[ReferenceEnthalpyMaterial {curve:&model.curve,reference_density_kg_m3:model.reference_density_kg_m3}],
                of_element:&element_materials,phase_conductivity:None},&initial_h,
            EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:vertices,max_elements:elements},max_steps:n})?;
        Ok(ScanningEnthalpyRun::new(cx,base,sources,1,None)?)
    };
    let segment=n.min(8);
    let observations=data.rows.iter().cloned().zip(&data.starts).zip(&data.sensor_lags).map(|((row,&start),&lag)| {
        let temporal=match start {
            None=>ScanTemporalObservation::instantaneous(row),
            Some(first)=>match lag {
                Some(lag)=>ScanTemporalObservation::first_order_lag(cx,row,&data.endpoints,first,
                    lag.time_constant_s,lag.initial_reading,MAX_STEPS)?,
                None=>ScanTemporalObservation::trapezoidal_mean(cx,row,&data.endpoints,first,MAX_STEPS)?,
            },
        };
        ScanObservationWithLoss::new(temporal,loss)
    })
        .collect::<std::result::Result<Vec<_>,_>>()?;
    let objective=ScanPowerObjective {endpoints_s:data.endpoints.clone(),observations,
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
            linear_restart:vertices.clamp(8,64),max_linear_cycles:if vertices<=8 {16} else {32},forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
            max_newton_iterations:64,energy_tolerance_j:1e-8},
        transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:vertices.clamp(8,64)},
        budget:ScanObjectiveBudget {max_steps:n,max_observations:data.contributions,max_retained_vertices:vertices*segment,
            max_retained_source_vertices:vertices*segment,max_schedule_controls:16,max_retained_radiation_entries:0},
        source_reference_time_s:0.};
    let controls=vec![ScanPowerControl {source:0,value:0,reference_w:1.},ScanPowerControl {source:0,value:1,reference_w:1.}];
    let mut study=ScanPowerStudy::new_joint_with_power_coordinates(cx,vec![source],controls,vec![],
        ScanPowerCoordinates::NonnegativeLinear,objective,factory,127)?;
    study.set_checkpointing(cx,ScanCheckpointBudget {segment_steps:segment,
        max_checkpoint_vertices:vertices*(n.div_ceil(segment)+1),max_forward_solves:2*n})?;
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
    Ok(Fit {powers,source_times_s:[0.0,model.switch_time_s,model.end_time_s],
        mesh_vertices:vertices,mesh_elements:elements,reference_mass_kg,
        measurements,exposure_starts_s,sensor_lags:data.sensor_lags,
        objective:replay.physical.value,evaluations:study.evaluations(),
        converged:report.solution.converged,violation:limits.max_scaled_violation()})
}
fn main() -> Result<()> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    let options=options(&args)?;
    let mut text=String::new();
    std::fs::File::open(&options.path)?.take(MAX_BYTES as u64+1).read_to_string(&mut text)?;
    let specimen_text=options.specimen.as_deref().map(|path| -> Result<String> {
        let mut contents=String::new();
        std::fs::File::open(path)?.take(SPECIMEN_BUDGET.max_bytes as u64+1).read_to_string(&mut contents)?;
        if contents.len()>SPECIMEN_BUDGET.max_bytes {return Err("specimen exceeds 2 MiB".into());}
        Ok(contents)
    }).transpose()?;
    let gate=CancelGate::new(); let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result=pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic);
        let model=match specimen_text.as_deref() {
            Some(input)=>ScanSpecimen::from_text(&cx,input,SPECIMEN_BUDGET)?,
            None=>ScanSpecimen::synthetic_reference(&cx)?,
        };
        fit_on_specimen(&cx,&model,&text,options.initial,options.loss,options.substeps)
    })?;
    println!("# specimen_mode={}; user declarations are not experimental validation",
        if options.specimen.is_some() {"explicit_file"} else {"synthetic_reference"});
    println!("# mesh_vertices={}; mesh_elements={}; reference_mass_kg={:.17e}",
        result.mesh_vertices,result.mesh_elements,result.reference_mass_kg);
    println!("# residual_loss={:?}; residuals below are raw, not clipped",options.loss);
    println!("# substeps_per_event_interval={}; time means use normalized trapezoidal readouts",options.substeps);
    println!("# objective={:.17e} converged={} max_scaled_violation={:.17e} evaluations={}",
        result.objective,result.converged,result.violation,result.evaluations);
    println!("time_s,absorbed_power_w");
    for (time,power) in result.source_times_s.into_iter().zip(&result.powers) {println!("{time:.17e},{power:.17e}");}
    for (i,(time,quantity,target,prediction,residual)) in result.measurements.iter().enumerate() {
        let (term,slope)=options.loss.value_derivative(*residual)?;
        if let Some(start)=result.exposure_starts_s[i] {
            if let Some(lag)=result.sensor_lags[i] {
                println!("# measurement={i} interpretation=first_order_lag response_start_s={start:.17e} response_end_s={time:.17e} tau_s={:.17e} initial_reading={:.17e}",
                    lag.time_constant_s,lag.initial_reading);
            } else {
                println!("# measurement={i} interpretation=time_mean exposure_start_s={start:.17e} exposure_end_s={time:.17e}");
            }
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

#[cfg(test)]
mod sensor_tests {
    use super::*;
    const DATA: &str=include_str!("sensor_measurements.csv");
    fn with_cx(f: impl FnOnce(&CancelGate,&Cx<'_>)) {
        let gate=CancelGate::new(); let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena|f(&gate,&Cx::new(&gate,arena,StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
            Budget::INFINITE,ExecMode::Deterministic)));
    }

    #[test]
    fn g1_lagged_csv_recovers_powers_and_exposes_instantaneous_misinterpretation() {
        with_cx(|_,cx| {
            for substeps in [1,4] {
                let actual=fit(cx,DATA,[0.,0.],ScanResidualLoss::Quadratic,substeps).unwrap();
                assert!(actual.converged && actual.violation<=1e-8 && actual.objective<1e-12);
                assert!((actual.powers[0]-2.).abs()<1e-4 && (actual.powers[1]-5.).abs()<1e-4);
                assert_eq!(actual.powers[2],0.);
                assert_eq!(actual.sensor_lags[0],None);
                assert_eq!(actual.sensor_lags[1],Some(SensorLag {time_constant_s:0.5,initial_reading:152.}));
                let wrong=DATA.replace("lag:1:2:0.5:152","2").replace("lag:1:2:0.5:0.52","2");
                let wrong=fit(cx,&wrong,[0.,0.],ScanResidualLoss::Quadratic,substeps).unwrap();
                assert!(wrong.converged && wrong.objective<1e-12);
                assert!((wrong.powers[0]-2.).abs()<1e-4);
                assert!((wrong.powers[1]-2.83833820809153).abs()<1e-4);
            }
        });
    }

    #[test]
    fn g3_sensor_csv_keeps_declared_clock_initial_reading_and_sample_accounting() {
        with_cx(|_,cx| {
            let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
            let mesh=ConductionMesh::new(complex,positions).unwrap();
            let text=format!("{HEADER}\nlag:0.25:1.75:0.5:151,specific_enthalpy_j_kg,154,10,volume,0;1;2;3;4;5\n1,temperature_k,310,1,vertex,0\n");
            let data=parse(cx,&mesh,&text,2).unwrap();
            for time in [0.25,1.0,1.75,2.0] {assert!(data.endpoints.contains(&time));}
            let first=data.starts[0].unwrap();
            assert_eq!(data.endpoints[first],0.25);
            assert_eq!(data.endpoints[data.rows[0].endpoint],1.75);
            assert_eq!(data.contributions,data.rows[0].endpoint-first+2);
            assert_eq!(data.sensor_lags,vec![Some(SensorLag {time_constant_s:0.5,initial_reading:151.}),None]);
        });
    }

    #[test]
    fn g0_sensor_csv_refuses_bad_units_policies_and_cancelled_input() {
        with_cx(|gate,cx| {
            let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
            let mesh=ConductionMesh::new(complex,positions).unwrap();
            for spec in ["lag:1:2:0:0.5","lag:1:2:-1:0.5","lag:1:2:NaN:0.5",
                "lag:1:2:0.5:60","lag:1:2:0.5:NaN","lag:0:2:0.5:0.5",
                "lag:2:1:0.5:0.5","lag:1:2:0.5","lag:1:2:0.5:0.5:9"] {
                let text=format!("{HEADER}\n{spec},liquid_mass_fraction,0.55,0.1,vertex,0\n");
                let error=parse(cx,&mesh,&text,1).err().unwrap().to_string();
                assert!(error.contains("line 2"),"{error}");
            }
            gate.request();
            assert!(parse(cx,&mesh,DATA,1).is_err());
        });
    }
}

#[cfg(test)]
#[path = "spatial_measurement_fit/specimen_tests.rs"]
mod specimen_tests;

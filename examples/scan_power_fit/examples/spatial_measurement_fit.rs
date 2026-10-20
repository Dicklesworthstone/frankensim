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
const HEADER: &str = "time_s,observable,target,scale,readout,selection";
const USAGE: &str = "spatial_measurement_fit measurements.csv [initial_first_W initial_second_W]";

struct Measurements<'m> {
    endpoints: Vec<f64>,
    rows: Vec<ScanSpatialObservation<'m>>,
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
fn parse<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh, text: &str) -> Result<Measurements<'m>> {
    if text.len()>MAX_BYTES { return Err("measurement file exceeds 128 KiB".into()); }
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
            let time=finite(fields[0])?;
            if time<=0. || time>2. { return Err("reference specimen measurements require 0 < time_s <= 2".into()); }
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
            Ok((time,ScanSpatialObservation {endpoint:0,readout:readout(cx,mesh,fields[4],fields[5])?,observable,target,scale}))
        })().map_err(|e|format!("measurement line {}: {e}",line_index+2))?;
        times.push(row.0); rows.push(row.1);
    }
    if rows.is_empty() { return Err("at least one measurement is required".into()); }
    // Preserve exact representable times and every source switch. Inserting
    // measurements cannot make the grid step over the one-second power change.
    let mut endpoints=times.clone(); endpoints.extend([1.0,2.0]);
    endpoints.sort_by(f64::total_cmp); endpoints.dedup();
    for (row,time) in rows.iter_mut().zip(times) {
        row.endpoint=endpoints.binary_search_by(|t|t.total_cmp(&time)).map_err(|_| "measurement endpoint missing")?;
    }
    cx.checkpoint().map_err(|_| "measurement parsing cancelled")?;
    Ok(Measurements {endpoints,rows})
}

struct Fit {
    powers: Vec<f64>,
    measurements: Vec<(f64, ScanObservable, f64, f64, f64)>,
    objective: f64,
    evaluations: usize,
    converged: bool,
    violation: f64,
}
fn fit(cx: &Cx<'_>, text: &str, initial: [f64;2]) -> Result<Fit> {
    if initial.iter().any(|p|!p.is_finite() || *p<0.) { return Err("initial watts must be finite and nonnegative".into()); }
    let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
    let mesh=ConductionMesh::new(complex,positions)?;
    let data=parse(cx,&mesh,text)?; let n=data.endpoints.len();
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
    let objective=ScanPowerObjective {endpoints_s:data.endpoints.clone(),observations:data.rows.clone(),
        sampling:ScanObjectiveSampling::EventLeft,
        step:EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
            linear_restart:8,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
            max_newton_iterations:64,energy_tolerance_j:1e-8},
        transpose:LinearConfig {tolerance:1e-10,max_iterations:1000,restart:8},
        budget:ScanObjectiveBudget {max_steps:n,max_observations:MAX_ROWS,max_retained_vertices:8*segment,
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
    Ok(Fit {powers,measurements,objective:replay.physical.value,evaluations:study.evaluations(),
        converged:report.solution.converged,violation:limits.max_scaled_violation()})
}
fn main() -> Result<()> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()!=1 && args.len()!=3 { return Err(USAGE.into()); }
    let mut text=String::new();
    std::fs::File::open(&args[0])?.take(MAX_BYTES as u64+1).read_to_string(&mut text)?;
    let initial=if args.len()==3 {[finite(&args[1])?,finite(&args[2])?]} else {[1.,1.]};
    let gate=CancelGate::new(); let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result=pool.scope(|arena|fit(&Cx::new(&gate,arena,
        StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic),&text,initial))?;
    println!("# Synthetic unit-cube reference model; not automatic experimental calibration.");
    println!("# objective={:.17e} converged={} max_scaled_violation={:.17e} evaluations={}",
        result.objective,result.converged,result.violation,result.evaluations);
    println!("time_s,absorbed_power_w");
    for (time,power) in [0.,1.,2.].into_iter().zip(&result.powers) {println!("{time:.17e},{power:.17e}");}
    for (i,(time,quantity,target,prediction,residual)) in result.measurements.iter().enumerate() {
        println!("# measurement={i} time_s={time:.17e} observable={quantity:?} target={target:.17e} prediction={prediction:.17e} standardized_residual={residual:.17e}");
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
            let solved=fit(cx,DATA,[0.,0.]).unwrap();
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
                let error=parse(cx,&mesh,&format!("{HEADER}\n{line}\n")).err().unwrap().to_string();
                assert!(error.contains("line 2"),"{error}");
            }
            assert!(parse(cx,&mesh,HEADER).is_err());
            assert!(parse(cx,&mesh,&format!("{HEADER}\n{}","1,specific_enthalpy_j_kg,152,1,vertex,0\n".repeat(MAX_ROWS+1))).is_err());
        });
    }
    #[test]
    fn g3_csv_preserves_input_order_and_inserts_each_source_switch() {
        with_cx(|cx| {
            let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[1.;3]);
            let mesh=ConductionMesh::new(complex,positions).unwrap();
            let text=format!("{HEADER}\n1.5,temperature_k,310,1,tet,0:1:2:3:4\n0.5,specific_enthalpy_j_kg,151,10,weighted,0:1;7:1\n");
            let data=parse(cx,&mesh,&text).unwrap();
            assert_eq!(data.endpoints,vec![0.5,1.,1.5,2.]);
            assert_eq!(data.rows.iter().map(|r|r.endpoint).collect::<Vec<_>>(),vec![2,0]);
            assert_eq!(data.rows[0].readout.evaluate(cx,&[310.;8]).unwrap(),310.);
        });
    }
}

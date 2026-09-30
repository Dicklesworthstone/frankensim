//! Fit spatial temperatures, a heater density and sensor bias from CSV readings.
//! Uses the production tetrahedral conduction and variational/L-BFGS owners.
//! The built-in 0.12 x 0.04 x 0.04 m slab has k=2 W/(m K), capacity
//! 1000 J/(m3 K), fixed 300 K at x=0 and insulated remaining faces. Its source
//! density is amplitude*(x/0.12), not an inferred arbitrary heat-source field.
//! `--nodes` prints the exact mesh indices used by `time_s,node,temperature_k,sigma_k`.
//! Without a CSV this runs labeled, noiseless synthetic observations. It does
//! not import arbitrary meshes or certify physical validity/identifiability.
//! `--substeps N` refines solver time only; observation knots and priors stay fixed.
use fs_ascent::conduction_assimilation::{ConductionWindowConfig, ConductionWindowModel, ConductionWindowPolicy};
use fs_ascent::transient::variational::{WeakConstraintWindow, WindowControl, WindowObjective};
use fs_ascent::transient::variational::joint::{JointEvaluation, JointWindow, JointWindowStudy, ParameterFamily};
use fs_ascent::transient::variational::study::StudySettings;
use fs_ascent::{LbfgsReport, StopReason};
use fs_ascent::transient::variational::intervals::IntervalScheme;
use fs_ascent::transient::variational::intervals::substeps::{CheckpointedIntervals, SubstepBudget, SubstepGrid};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig, StepLinearization};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use std::io::Read;
type Error = Box<dyn std::error::Error>;

struct Domain { mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel, profile: Vec<f64> }
impl Domain {
    fn new() -> Result<Self, ConductionError> {
        let (complex, positions) = box_grid([4,2,2],[0.12,0.04,0.04]);
        let mesh = ConductionMesh::new(complex,positions)?;
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("fixed",|f|on_box_face(f.centroid[0],0.0),ThermalBc::dirichlet(300.0)?)?
            .adiabatic_remainder().finish()?;
        let material = ConductivityModel::isotropic_declared(2.0)?;
        let profile = mesh.positions().iter().map(|p|p[0]/0.12).collect();
        Ok(Self {mesh,boundary,material,profile})
    }
}
fn config() -> ConductionWindowConfig {
    ConductionWindowConfig { step: StepConfig { linear: LinearConfig {
        tolerance:1e-11,max_iterations:2000,restart:24 }, energy_tolerance_j:1e-8 },
        nonlinear:None,max_vertices:128,max_elements:1024,max_intervals:8,max_parameters:2 }
}
#[derive(Debug,Clone,Copy)]
struct Reading { time:f64,node:usize,value:f64,sigma:f64 }
fn parse(text:&str,nodes:usize) -> Result<Vec<Reading>,String> {
    if text.len()>65_536 {return Err("CSV exceeds 64 KiB".into());}
    let mut lines=text.lines();
    if lines.next().map(str::trim)!=Some("time_s,node,temperature_k,sigma_k") {return Err("invalid CSV header".into());}
    let mut rows:Vec<Reading>=Vec::new();
    for line in lines.filter(|s|!s.trim().is_empty()) {
        if rows.len()==256 {return Err("more than 256 observations".into());}
        let c:Vec<_>=line.split(',').map(str::trim).collect();
        if c.len()!=4 {return Err("expected four CSV columns".into());}
        let row=Reading {time:c[0].parse().map_err(|_|"invalid time")?,node:c[1].parse().map_err(|_|"invalid node")?,
            value:c[2].parse().map_err(|_|"invalid temperature")?,sigma:c[3].parse().map_err(|_|"invalid sigma")?};
        if !row.time.is_finite() || row.time<0.0 || row.time>10.0 || row.node>=nodes
            || !row.value.is_finite() || row.value<0.0 || !row.sigma.is_finite() || row.sigma<=0.0
            || rows.last().is_some_and(|r|r.time>row.time)
            || rows.iter().any(|r|r.time==row.time && r.node==row.node)
        {return Err("invalid, repeated, or unordered observation (time limit 10 s)".into());}
        rows.push(row);
    }
    if rows.is_empty() || !rows.iter().any(|r|r.time>0.0) {return Err("at least one positive observation time required".into());}
    Ok(rows)
}
#[derive(Clone)]
struct Sample { frame:usize, slot:Option<usize>, value:f64, sigma:f64 }
struct Model<'a> { domain:&'a Domain, engine:&'a BackwardEuler<'a>, source:ScalarField,
    bias:f64, samples:Vec<Sample>, times:Vec<f64>, dimension:usize }
impl ConductionWindowModel for Model<'_> {
    fn engine(&self)->&BackwardEuler<'_> {self.engine}
    fn problem(&self,_:usize)->Result<ConductionProblem<'_>,ConductionError> {
        Ok(ConductionProblem {mesh:&self.domain.mesh,boundary:&self.domain.boundary,
            material:&self.domain.material,element_materials:None,source:&self.source})
    }
    fn parameter_count(&self)->usize {2}
    fn parameter_pullback(&self,_:usize,cx:&Cx<'_>,step:&StepLinearization<'_>,lambda:&[f64],out:&mut[f64])
        ->Result<(),ConductionError> {
        let density=step.source_density_pullback(cx,lambda)?;
        out[0]=density.iter().zip(&self.domain.profile).map(|(g,p)|g*p).sum();out[1]=0.0;Ok(())
    }
}
impl Model<'_> {
    fn residual(&self,s:&Sample,x:&[f64])->f64 {
        let temperature=s.slot.map_or(300.0,|slot|x[s.frame*self.dimension+slot]);
        (temperature+self.bias-s.value)/s.sigma
    }
}
impl WindowObjective for Model<'_> {
    fn evaluate(&self,times:&[f64],dimension:usize,x:&[f64],out:&mut[f64],check:&mut dyn FnMut()->bool)->Result<f64,String> {
        if times!=self.times || dimension!=self.dimension || x.len()!=dimension*times.len() || out.len()!=x.len() {
            return Err("observation frame/state shape changed".into());
        }
        out.fill(0.0);let mut cost=0.0;
        for s in &self.samples {
            if check() {return Err("observation cancelled".into());}
            let r=self.residual(s,x);cost+=0.5*r*r;
            if let Some(slot)=s.slot {out[s.frame*dimension+slot]+=r/s.sigma;}
        }
        Ok(cost)
    }
    fn parameter_partials(&self,_:&[f64],_:usize,x:&[f64],out:&mut[f64],check:&mut dyn FnMut()->bool)->Result<(),String> {
        out.fill(0.0);
        for s in &self.samples {if check() {return Err("sensor derivative cancelled".into());} out[1]+=self.residual(s,x)/s.sigma;}
        Ok(())
    }
}
struct Family<'a> { domain:&'a Domain,engine:&'a BackwardEuler<'a>,samples:Vec<Sample>,times:Vec<f64>,dimension:usize }
impl<'a> ParameterFamily for Family<'a> {
    type Model=Model<'a>;
    fn instantiate(&self,p:&[f64],check:&mut dyn FnMut()->bool)->Result<Self::Model,String> {
        if check() {return Err("factory cancelled".into());}
        if p.len()!=2 || p.iter().any(|v|!v.is_finite()) {return Err("two finite parameter coordinates required".into());}
        let source=ScalarField::nodal("heater profile",self.domain.mesh.vertex_count(),
            self.domain.profile.iter().map(|v|p[0]*v).collect()).map_err(|e|e.to_string())?;
        Ok(Model {domain:self.domain,engine:self.engine,source,bias:p[1],samples:self.samples.clone(),
            times:self.times.clone(),dimension:self.dimension})
    }
}
fn synthetic(cx:&Cx<'_>,domain:&Domain,engine:&BackwardEuler<'_>,substeps:usize)->Result<Vec<Reading>,Error> {
    let times=[0.0,0.25,0.5,1.0];let mut state=vec![300.0;domain.mesh.vertex_count()];let mut rows=Vec::new();
    let source=ScalarField::nodal("synthetic heater",state.len(),domain.profile.iter().map(|v|2000.0*v).collect())?;
    let grid=SubstepGrid::uniform(&times,substeps,256)?;
    for (k,&time) in times.iter().enumerate() {
        if k>0 {
            for step in grid.knot_indices()[k-1]..grid.knot_indices()[k] {
                let dt=grid.fine_times()[step+1]-grid.fine_times()[step];
                state=engine.advance(cx,ConductionProblem {mesh:&domain.mesh,boundary:&domain.boundary,
                    material:&domain.material,element_materials:None,source:&source},None,&state,dt,config().step)?.temperature;
            }
        }
        // Includes a fixed-boundary sensor: it informs bias but not a free state.
        for node in [4,13,22,40] {rows.push(Reading {time,node,value:state[node]+0.08,sigma:0.02});}
    }
    Ok(rows)
}
struct Fit { before:f64,evaluation:JointEvaluation,fields:Vec<Vec<f64>>,times:Vec<f64>,report:LbfgsReport }
fn fit(cx:&Cx<'_>,domain:&Domain,engine:&BackwardEuler<'_>,rows:&[Reading],model_sigma:f64,substeps:usize)->Result<Fit,Error> {
    if rows.is_empty() || rows.len()>256 || !model_sigma.is_finite() || model_sigma<=0.0 || !(1..=32).contains(&substeps) {return Err("invalid fit inputs".into());}
    let mut times=vec![0.0];
    for r in rows {if r.time>*times.last().unwrap() {times.push(r.time);}}
    let mut numerical=config();
    if times.len()<2 || times.len()-1>numerical.max_intervals {return Err("too many observation intervals".into());}
    let grid=SubstepGrid::uniform(&times,substeps,256)?;
    numerical.max_intervals=256; // Fine solver intervals, not optimization knots.
    let policy=ConductionWindowPolicy::new(cx,&domain.mesh,&domain.boundary,grid.fine_times(),numerical,&mut||false)?;
    let n=policy.dimension();let mut samples=Vec::new();
    for r in rows {
        let frame=times.binary_search_by(|t|t.total_cmp(&r.time)).map_err(|_|"observation not on the admitted clock")?;
        if r.node>=domain.mesh.vertex_count() {return Err("unknown mesh node".into());}
        samples.push(Sample {frame,slot:policy.slot_of(r.node),value:r.value,sigma:r.sigma});
    }
    let family=Family {domain,engine,samples,times:times.clone(),dimension:n};
    let window=WeakConstraintWindow::new(&times,&vec![300.0;n*times.len()],&vec![1.0;n],&vec![0.1;n],
        &vec![model_sigma;n*(times.len()-1)],10_000)?;
    let joint=JointWindow::new(&window,&family,&[1000.0,0.0],&[1000.0,0.1],&[2000.0,0.2],6)?;
    // The default retains the original one-step numerical/gradient path.
    let (before,evaluation,report)=if substeps==1 {
        run_study(&joint,policy.clone())?
    } else {
        run_study(&joint,CheckpointedIntervals::new(policy.clone(),grid,SubstepBudget {
            max_state_components:128,max_parameters:2,checkpoints:6,replayed_substeps:512,
        }))?
    };
    let fields=evaluation.window.states.chunks(n).map(|x|policy.expand_field(x)).collect::<Result<Vec<_>,_>>()?;
    Ok(Fit {before,evaluation,fields,times,report})
}
fn run_study<F:ParameterFamily,S:IntervalScheme<F::Model>>(joint:&JointWindow<'_,F>,policy:S)
    ->Result<(f64,JointEvaluation,LbfgsReport),Error>
{
    let settings=StudySettings {memory:17,gradient_tolerance:1e-5,max_evaluations:2500,max_optimizer_components:100_000};
    let mut control=WindowControl::new(2500,20_000,100_000);
    let mut study=JointWindowStudy::new(joint,&vec![0.0;joint.control_dimension()],policy,settings,&mut control,&mut||false)?;
    let before=study.accepted().value;let report=study.run(500,&mut control,&mut||false)?;
    Ok((before,study.accepted().clone(),report))
}
fn options(args:&[String])->Result<(Option<&str>,f64,usize),String> {
    let mut path=None;let mut sigma=None;let mut substeps=None;let mut i=0;
    while i<args.len() {
        if args[i]=="--model-sigma" && sigma.is_none() {
            let value:f64=args.get(i+1).ok_or("--model-sigma needs a scale in kelvin")?.parse().map_err(|_|"invalid model sigma")?;
            if !value.is_finite() || value<=0.0 {return Err("model sigma must be positive and finite".into());}
            sigma=Some(value);i+=2;
        } else if args[i]=="--substeps" && substeps.is_none() {
            let value:usize=args.get(i+1).ok_or("--substeps needs an integer")?.parse().map_err(|_|"invalid substep count")?;
            if !(1..=32).contains(&value) {return Err("substeps must be in 1..=32".into());}
            substeps=Some(value);i+=2;
        } else if !args[i].starts_with("--") && path.is_none() {path=Some(args[i].as_str());i+=1;}
        else {return Err("usage: conduction_assimilate [readings.csv] [--model-sigma kelvin] [--substeps 1..32] or --nodes".into());}
    }
    Ok((path,sigma.unwrap_or(0.03),substeps.unwrap_or(1)))
}
fn main()->Result<(),Error> {
    let args:Vec<_>=std::env::args().skip(1).collect();let domain=Domain::new()?;
    if args.len()==1 && args[0]=="--nodes" {
        println!("node,x_m,y_m,z_m,prescribed");
        for (node,p) in domain.mesh.positions().iter().enumerate() {println!("{node},{},{},{},{}",p[0],p[1],p[2],on_box_face(p[0],0.0));}
        return Ok(());
    }
    let (path,sigma,substeps)=options(&args)?;
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena|->Result<(),Error> {
        let cx=Cx::new(&gate,arena,StreamKey {seed:61,kernel_id:820,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let engine=BackwardEuler::uniform(&cx,&domain.mesh,VolumetricHeatCapacity::declared(1000.0)?)?;
        let (rows,source)=if let Some(path)=path {let mut text=String::new();std::fs::File::open(path)?.take(65_537).read_to_string(&mut text)?;
            (parse(&text,domain.mesh.vertex_count())?,"csv")
        } else {(synthetic(&cx,&domain,&engine,substeps)?,"synthetic-noiseless")};
        let result=fit(&cx,&domain,&engine,&rows,sigma,substeps)?;let e=&result.evaluation;
        println!("source={source} solver=P1-backward-Euler stop={:?} iterations={} evaluations={} model_sigma_k={sigma}",
            result.report.reason,result.report.iters,result.report.evals);
        println!("source_amplitude_w_m3={:.9} sensor_bias_k={:.9} objective_before={:.9e} objective_after={:.9e}",
            e.parameters[0],e.parameters[1],result.before,e.value);
        println!("observation_loss={:.9e} background_loss={:.9e} model_loss={:.9e} parameter_loss={:.9e}",
            e.window.observation_value,e.window.background_value,e.window.model_value,e.parameter_penalty);
        println!("substeps_per_interval={substeps} native_forward_steps={} native_replayed_steps={} controls={}",
            e.window.accepted_steps,e.window.replayed_steps,e.controls.len());
        println!("time_s,node,x_m,y_m,z_m,reconstructed_temperature_k");
        for (&time,field) in result.times.iter().zip(&result.fields) {for (node,p) in domain.mesh.positions().iter().enumerate() {
            println!("{time},{node},{},{},{},{:.9}",p[0],p[1],p[2],field[node]);
        }}
        println!("scope=numerical MAP on a declared slab; no physical-validation, calibrated-uncertainty or source-identifiability claim");
        if result.report.reason!=StopReason::GradNorm {return Err("reconstruction stopped without reaching its gradient tolerance".into());}
        Ok(())
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn csv_rejects_invalid_identity_units_and_order() {
        let h="time_s,node,temperature_k,sigma_k\n";
        assert!(parse(&format!("{h}0,4,300.08,0.02\n1,40,301,0.02\n"),45).is_ok());
        for data in ["1,45,300,0.1", "1,4,NaN,0.1", "1,4,300,0", "1,4,300,0.1\n0,13,300,0.1",
            "1,4,300,0.1\n1,4,300,0.2", "11,4,300,0.1", "0,4,300,0.1"] {
            assert!(parse(&format!("{h}{data}"),45).is_err());
        }
    }
    #[test]
    fn synthetic_fit_uses_pde_and_preserves_prescribed_temperatures() {
        let d=Domain::new().unwrap();let gate=CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx=Cx::new(&gate,arena,StreamKey {seed:61,kernel_id:820,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
            let engine=BackwardEuler::uniform(&cx,&d.mesh,VolumetricHeatCapacity::declared(1000.0).unwrap()).unwrap();
            let rows=synthetic(&cx,&d,&engine,1).unwrap();let got=fit(&cx,&d,&engine,&rows,0.03,1).unwrap();
            assert_eq!(got.report.reason,StopReason::GradNorm);
            assert!(got.evaluation.value<got.before*0.001);
            assert!((got.evaluation.parameters[0]-2000.0).abs()<10.0);
            assert!((got.evaluation.parameters[1]-0.08).abs()<0.01);
            for field in got.fields {for &(node,value) in d.boundary.dirichlet() {assert_eq!(field[node],value);}}
        });
    }
    #[test]
    fn model_error_options_are_explicit_and_finite() {
        for args in [vec!["--model-sigma"],vec!["--model-sigma","0"],vec!["--model-sigma","NaN"],vec!["--unknown"]] {
            assert!(options(&args.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
        }
    }
}

#[cfg(test)]
#[path = "conduction_assimilate/substeps_tests.rs"]
mod substeps_tests;

//! Multi-experiment mode of transient_fit, using the reusable campaign API.
//! CSV: experiment,amplitude,time_s,channel,value,sigma
//! Amplitude is a known experiment initial condition, not a fitted parameter.
//! Channel 0 observes x+b; channel 1 observes 2*x+b. Shared decisions are decay
//! rate [0.1,2]/s and sensor bias [-1,1] in the common signal unit. Sigma uses
//! that signal unit; Huber threshold is dimensionless. This is an illustrative
//! reduced-model fit, not a validated sensor or cooling model.
//! --shared-sigma adds one independent reference-error source per experiment,
//! shared by all its channels/times. CSV sigma then means INDEPENDENT noise.
//! This fixed-covariance GLS mode is mutually exclusive with --huber.
use fs_ascent::{SqpStop, transient::{TransientConfig,
    campaign::{CampaignControl, CampaignStudy, Experiment, TransientCampaign},
    observations::{ObservedFamily, SensorData, SensorFamily, SensorLoss, SensorModel, SensorReading}}};
use fs_time::{PiController, adaptive::adjoint::{OdeVjp, trajectory::{RecordingConfig, ReplayBudget}}};
use std::{collections::BTreeMap, io::Read};

// This file is itself included through #[path], so a bare `mod shared;`
// would resolve beside it (transient_fit/shared.rs), not under campaign/.
#[path = "campaign/shared.rs"]
mod shared;

const BOUNDS: [[f64;2];2] = [[0.1,2.0],[-1.0,1.0]];
const HEADER: &str = "experiment,amplitude,time_s,channel,value,sigma";
const MAX_ROWS: usize = 512;
const MAX_EXPERIMENTS: usize = 32;
const MAX_BYTES: usize = 131072;
struct Decay { amplitude: f64 }
struct Model { point: [f64;2], initial: [f64;1] }
impl SensorFamily for Decay {
    type Model = Model;
    fn bounds(&self) -> &[[f64;2]] { &BOUNDS }
    fn instantiate(&self, p: &[f64]) -> Result<Model,String> {
        Ok(Model { point: [p[0],p[1]], initial: [self.amplitude] })
    }
}
impl OdeVjp for Model {
    fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {2}
    fn rhs(&self,_:f64,x:&[f64],out:&mut[f64]) {out[0]=-self.point[0]*x[0];}
    fn rhs_vjp(&self,_:f64,x:&[f64],b:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        xb[0]=-self.point[0]*b[0];pb[0]=-x[0]*b[0];pb[1]=0.0;Ok(())
    }
}
impl SensorModel for Model {
    fn initial_values(&self)->&[f64] {&self.initial}
    fn initial_vjp(&self,_:&[f64],pb:&mut[f64])->Result<(),String> {pb.fill(0.0);Ok(())}
    fn predict(&self,channel:u64,_:f64,x:&[f64])->Result<f64,String> {
        if channel>1 {return Err("only channels 0 and 1 are supported".into());}
        Ok((channel+1) as f64*x[0]+self.point[1])
    }
    fn prediction_vjp(&self,channel:u64,_:f64,_:&[f64],b:f64,xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        xb[0]=(channel+1) as f64*b;pb[0]=0.0;pb[1]=b;Ok(())
    }
}
struct Data { id: u64, family: ObservedFamily<Decay> }
fn parse(text:&str,loss:SensorLoss)->Result<Vec<Data>,String> {
    if text.len()>MAX_BYTES {return Err("campaign CSV exceeds 128 KiB".into());}
    let mut lines=text.lines();
    if lines.next().map(str::trim)!=Some(HEADER) {return Err(format!("expected CSV header {HEADER}"));}
    let mut groups:BTreeMap<u64,(f64,Vec<SensorReading>)>=BTreeMap::new();
    let mut count=0;
    for (index,line) in lines.enumerate() {
        let line=line.trim();if line.is_empty() {continue;}
        if count==MAX_ROWS {return Err("campaign CSV exceeds 512 readings".into());}
        let fields=line.split(',').map(str::trim).collect::<Vec<_>>();
        if fields.len()!=6 {return Err(format!("row {} requires six columns",index+2));}
        let id:u64=fields[0].parse().map_err(|_|format!("invalid experiment ID at row {}",index+2))?;
        let amplitude:f64=fields[1].parse().map_err(|_|"invalid amplitude")?;
        let time:f64=fields[2].parse().map_err(|_|"invalid time")?;
        let channel:u64=fields[3].parse().map_err(|_|"invalid channel")?;
        let value:f64=fields[4].parse().map_err(|_|"invalid reading")?;
        let sigma:f64=fields[5].parse().map_err(|_|"invalid sigma")?;
        if !amplitude.is_finite() || amplitude<=0.0 || time<0.0 || channel>1 {
            return Err(format!("invalid amplitude, time or channel at row {}",index+2));
        }
        let reading=SensorReading::new(time,channel,value,sigma,loss)?;
        if !groups.contains_key(&id) && groups.len()==MAX_EXPERIMENTS {
            return Err("campaign CSV exceeds 32 experiments".into());
        }
        let group=groups.entry(id).or_insert_with(||(amplitude,Vec::new()));
        if group.0.to_bits()!=amplitude.to_bits() {return Err(format!("experiment {id} changes its initial amplitude"));}
        group.1.push(reading);count+=1;
    }
    if groups.is_empty() {return Err("campaign has no readings".into());}
    groups.into_iter().map(|(id,(amplitude,mut rows))| {
        // Stable sorting preserves declaration order within repeated times.
        rows.sort_by(|a,b|a.time().total_cmp(&b.time()));
        Ok(Data {id,family:ObservedFamily::new(Decay{amplitude},SensorData::new(&rows,MAX_ROWS)?)})
    }).collect()
}
fn config(end:f64)->TransientConfig {TransientConfig {start:0.0,initial_step:0.1,
    recording:RecordingConfig{end,rtol:1e-10,atol:1e-12,controller:PiController::default(),max_workspace_components:4096},
    max_state_components:1,max_samples:MAX_ROWS,max_attempts:10000,max_records:10000,
    replay:ReplayBudget{checkpoints:64,replayed_steps:100000},max_kkt_dimension:6}}
fn campaign(data:&[Data])->Result<TransientCampaign<'_,ObservedFamily<Decay>>,Box<dyn std::error::Error>> {
    let cases=data.iter().map(|row|Experiment {id:row.id,family:&row.family,weight:1.0,
        config:config(*row.family.data().times().last().expect("nonempty checked data"))}).collect();
    Ok(TransientCampaign::new(cases,MAX_EXPERIMENTS,4+5*MAX_EXPERIMENTS)?)
}
fn synthetic()->String {
    let mut text=format!("{HEADER}\n");
    for (id,amplitude,times) in [(7,1.2,[0.0,0.2,1.0,2.0]),(42,2.3,[0.1,0.5,1.4,3.0])] {
        for (i,time) in times.into_iter().enumerate() {
            let channel=i%2;let value=(channel+1) as f64*amplitude*fs_math::det::exp(-0.7*time)+0.15;
            text.push_str(&format!("{id},{amplitude},{time},{channel},{value:.17},0.1\n"));
        }
    }
    text
}
pub fn run(args:&[String])->Result<(),Box<dyn std::error::Error>> {
    let mut path=None;let mut loss=SensorLoss::Quadratic;let mut shared_sigma=None;let mut i=0;
    while i<args.len() {
        if args[i]=="--huber" {
            let threshold:f64=args.get(i+1).ok_or("--huber requires a threshold")?.parse()?;
            loss=SensorLoss::Huber{threshold};i+=2;
            // Check even when no CSV rows were yet parsed.
            SensorReading::new(0.0,0,0.0,1.0,loss)?;
        } else if args[i]=="--shared-sigma" {
            let sigma:f64=args.get(i+1).ok_or("--shared-sigma requires a scale")?.parse()?;
            if !sigma.is_finite() || sigma<0.0 || shared_sigma.replace(sigma).is_some() {
                return Err("shared sigma must be finite, nonnegative and specified once".into());
            }
            i+=2;
        } else if path.is_none() && !args[i].starts_with("--") {path=Some(&args[i]);i+=1;}
        else {return Err("usage: transient_fit --campaign [readings.csv] [--huber threshold | --shared-sigma scale]".into());}
    }
    if shared_sigma.is_some() && loss!=SensorLoss::Quadratic {
        return Err("--shared-sigma and --huber cannot be combined".into());
    }
    let (text,source)=if let Some(path)=path {
        let mut text=String::new();std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;(text,"csv")
    } else {(synthetic(),"synthetic-noiseless")};
    let data=parse(&text,loss)?;
    if let Some(sigma)=shared_sigma {return shared::run(data,sigma,source);}
    let campaign=campaign(&data)?;
    let mut control=CampaignControl::new(1500,1500*data.len());
    let mut study=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut control,&mut||false)?;
    let initial=study.accepted().value;let report=study.run(1e-6,100,1500,&mut||false)?;
    println!("source={source} experiments={} loss={loss:?} stop={:?} work={:?}",data.len(),report.stop,study.work());
    println!("rate_per_s={:.10} bias={:.10} objective_before={:.12e} objective_after={:.12e}",
        study.optimizer().point()[0],study.optimizer().point()[1],initial,study.accepted().value);
    for e in &study.accepted().experiments {println!("experiment={} observations={} objective={:.12e}",e.id,e.result.observations,e.result.value);}
    println!("scope=numerical-shared-parameter-fit; no physical-validation or identifiability claim");
    if report.stop!=SqpStop::Converged {return Err("campaign stopped without satisfying its local KKT tolerance".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn grouped_csv_fits_shared_parameters_through_the_production_campaign() {
        let data=parse(&synthetic(),SensorLoss::Quadratic).unwrap();let campaign=campaign(&data).unwrap();
        let mut work=CampaignControl::new(1500,3000);
        let mut study=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut work,&mut||false).unwrap();
        assert_eq!(study.run(1e-6,100,1500,&mut||false).unwrap().stop,SqpStop::Converged);
        for (p,target) in study.optimizer().point().iter().zip([0.7,0.15]) {assert!((p-target).abs()<2e-6);}
        assert_eq!(study.accepted().experiments.len(),2);
    }
    #[test]
    fn grouped_csv_refuses_changed_initial_conditions_and_invalid_noise() {
        for rows in ["1,1.2,0,0,1.35,0.1\n1,2.3,1,0,1.0,0.1\n",
            "1,1.2,0,0,1.35,0\n", "1,1.2,0,4,1.35,0.1\n", "1,1.2,NaN,0,1.35,0.1\n"] {
            assert!(parse(&format!("{HEADER}\n{rows}"),SensorLoss::Quadratic).is_err());
        }
    }
}

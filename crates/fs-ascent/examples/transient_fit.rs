//! Fit y(t) = amplitude * exp(-rate*t) + bias from irregular CSV observations.
//! This reduced model is illustrative, not a validated cooling/material law.
//! Run with `--features transient-design --example transient_fit -- readings.csv`.
//! CSV header: time_s,value. Times are ordered; values use the caller's common
//! signal unit. The explicit box is rate in [0.1,2]/s, amplitude in [0.2,3],
//! bias in [-0.5,0.5]. No file means labeled, noiseless synthetic data.
//! Use `--lag [readings.csv]` for joint physical-rate/sensor-response fitting.
use fs_ascent::transient::{TransientConfig, TransientFamily, TransientModel, TransientStudy};
use fs_ascent::SqpStop;
use fs_time::adaptive::adjoint::{OdeVjp, trajectory::{RecordingConfig, ReplayBudget, samples::SampleObjective}};
use fs_time::PiController;
use std::io::Read;
use std::sync::Arc;

#[path = "transient_fit/campaign.rs"]
mod campaign;
#[path = "transient_fit/lag.rs"]
mod lag;

const BOUNDS: [[f64;2];3] = [[0.1,2.0],[0.2,3.0],[-0.5,0.5]];
struct Data { times: Vec<f64>, values: Arc<Vec<f64>> }
impl Data {
    fn parse(text: &str) -> Result<Self, String> {
        if text.len()>65_536 {return Err("CSV exceeds the 64 KiB input cap".into());}
        let mut lines=text.lines();
        if lines.next().map(str::trim)!=Some("time_s,value") {return Err("expected CSV header time_s,value".into());}
        let (mut times,mut values)=(Vec::new(),Vec::new());
        for (i,line) in lines.enumerate() {
            let line=line.trim(); if line.is_empty() {continue;}
            if times.len()>=128 {return Err("CSV exceeds the 128 observation cap".into());}
            let fields:Vec<_>=line.split(',').collect();
            if fields.len()!=2 {return Err(format!("row {} requires two columns",i+2));}
            let time:f64=fields[0].trim().parse().map_err(|_|format!("invalid time at row {}",i+2))?;
            let value:f64=fields[1].trim().parse().map_err(|_|format!("invalid value at row {}",i+2))?;
            if !time.is_finite() || time<0.0 || !value.is_finite()
                || times.last().is_some_and(|previous| *previous>time) {
                return Err(format!("non-finite value or unordered/nonnegative-time violation at row {}",i+2));
            }
            times.push(time);values.push(value);
        }
        if times.len()<3 || times.windows(2).filter(|p|p[0]!=p[1]).count()<2 {
            return Err("at least three distinct observation times are required".into());
        }
        Ok(Self {times,values:Arc::new(values)})
    }
    fn synthetic() -> Self {
        let times=vec![0.0,0.13,0.4,0.8,1.3,2.0,3.0,4.0];
        let values=Arc::new(times.iter().map(|t|1.2*fs_math::det::exp(-0.7*t)+0.15).collect());
        Self {times,values}
    }
}
struct Model { point:[f64;3], initial:[f64;1], values:Arc<Vec<f64>> }
impl OdeVjp for Model {
    fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {3}
    fn rhs(&self,_:f64,state:&[f64],out:&mut[f64]) {out[0]=-self.point[0]*state[0];}
    fn rhs_vjp(&self,_:f64,state:&[f64],seed:&[f64],x:&mut[f64],p:&mut[f64])->Result<(),String> {
        x[0]=-self.point[0]*seed[0];p[0]=-state[0]*seed[0];p[1]=0.0;p[2]=0.0;Ok(())
    }
}
impl SampleObjective for Model {
    fn evaluate(&self,i:usize,_:f64,state:&[f64],x:&mut[f64],p:&mut[f64])->Result<f64,String> {
        let r=state[0]+self.point[2]-self.values[i];x[0]=r;p[0]=0.0;p[1]=0.0;p[2]=r;Ok(0.5*r*r)
    }
}
impl TransientModel for Model {
    fn initial_values(&self)->&[f64] {&self.initial}
    fn initial_vjp(&self,bar:&[f64],out:&mut[f64])->Result<(),String> {out[0]=0.0;out[1]=bar[0];out[2]=0.0;Ok(())}
}
impl TransientFamily for Data {
    type Model=Model;
    fn bounds(&self)->&[[f64;2]] {&BOUNDS}
    fn sample_times(&self)->&[f64] {&self.times}
    fn instantiate(&self,point:&[f64])->Result<Model,String> {
        Ok(Model{point:[point[0],point[1],point[2]],initial:[point[1]],values:self.values.clone()})
    }
}
fn config(end:f64)->TransientConfig {TransientConfig {
    start:0.0,initial_step:0.1,
    recording:RecordingConfig{end,rtol:1e-10,atol:1e-12,controller:PiController::default(),max_workspace_components:4096},
    max_state_components:1,max_samples:128,max_attempts:10000,max_records:10000,
    replay:ReplayBudget{checkpoints:64,replayed_steps:100000},max_kkt_dimension:9,
}}
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.first().is_some_and(|arg| arg == "--campaign") { return campaign::run(&args[1..]); }
    if args.first().is_some_and(|arg| arg == "--lag") { return lag::run(&args[1..]); }
    if args.len()>1 {return Err("usage: transient_fit [readings.csv]".into());}
    let (data,source)=if let Some(path)=args.first() {
        let mut text=String::new();std::fs::File::open(path)?.take(65_537).read_to_string(&mut text)?;
        (Data::parse(&text)?,"csv")
    } else {(Data::synthetic(),"synthetic-noiseless")};
    let end=*data.times.last().ok_or("missing observation times")?;
    let mut study=TransientStudy::new(&data,&[1.4,2.2,-0.2],config(end),&mut||false)?;
    let initial=study.accepted().value;let report=study.run(1e-7,100,1500,&mut||false)?;
    let p=study.optimizer().point();
    println!("source={source} model=single-exponential stop={:?} iterations={} evaluations={}",
        report.stop,study.optimizer().iterations(),study.optimizer().evaluations());
    println!("rate_per_s={:.10} amplitude={:.10} bias={:.10} objective_before={:.12e} objective_after={:.12e}",
        p[0],p[1],p[2],initial,study.accepted().value);
    println!("scope=numerical-fit; no physical-validation or parameter-identifiability claim");
    if report.stop!=SqpStop::Converged {return Err("fit stopped without satisfying its local KKT tolerance".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strict_csv_admission() {
        assert!(Data::parse("time_s,value\n0,1.35\n0.4,1.05694049\n2,0.44591636\n").is_ok());
        for text in ["time,value\n0,1\n1,2\n2,3", "time_s,value\n0,1\n2,2\n1,3",
            "time_s,value\n0,1\n1,NaN\n2,3", "time_s,value\n0,1\n0,2\n1,3"] {
            assert!(Data::parse(text).is_err());
        }
    }
    #[test]
    fn csv_and_synthetic_observations_use_the_same_production_fit() {
        let data=Data::synthetic();let mut text=String::from("time_s,value\n");
        for (t,v) in data.times.iter().zip(data.values.iter()) {text.push_str(&format!("{t:.17},{v:.17}\n"));}
        let parsed=Data::parse(&text).unwrap();let mut study=TransientStudy::new(&parsed,&[1.4,2.2,-0.2],config(4.0),&mut||false).unwrap();
        assert_eq!(study.run(1e-7,100,1500,&mut||false).unwrap().stop,SqpStop::Converged);
        for (actual,target) in study.optimizer().point().iter().zip([0.7,1.2,0.15]) {assert!((actual-target).abs()<2e-5);}
    }
}

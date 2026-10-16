//! CSV-driven joint fitting of physical decay and instrument response time.
//!
//! Model: x'=-k*x, x(0)=1.2; tau*z'=x-z, z(0)=0. Channel 0 is an
//! instantaneous reference x; channel 1 reads z+b. Decisions: k [0.1,2]/s,
//! tau [0.05,1] s and b [-0.5,0.5] in the common signal unit. These known
//! initial conditions and reference-channel assumptions are part of this
//! example, not inferred facts about arbitrary CSV files. Use the library
//! adapter to declare a different instrument/model. No validation claim.
use fs_ascent::{SqpStop, transient::{TransientConfig, TransientStudy,
    observations::{ObservedFamily, SensorData, SensorFamily, SensorLoss, SensorModel, SensorReading,
        lag::{LagInitial, LagSensor, LagValue, LaggedFamily}}}};
use fs_time::{PiController, adaptive::adjoint::{OdeVjp, trajectory::{RecordingConfig, ReplayBudget}}};
use std::io::Read;

const HEADER: &str = "time_s,channel,value,sigma";
const MAX_BYTES: usize = 65_536;
const MAX_ROWS: usize = 256;
const BOUNDS: [[f64;2];3] = [[0.1,2.0],[0.05,1.0],[-0.5,0.5]];
struct Physical;
struct Model { rate: f64 }
impl SensorFamily for Physical {
    type Model = Model;
    fn bounds(&self) -> &[[f64;2]] { &BOUNDS }
    fn instantiate(&self, p: &[f64]) -> Result<Model,String> {
        if p.len()!=3 {return Err("three parameters required".into());}
        Ok(Model{rate:p[0]})
    }
}
impl OdeVjp for Model {
    fn dimension(&self)->usize {1}
    fn parameter_count(&self)->usize {3}
    fn rhs(&self,_:f64,x:&[f64],out:&mut[f64]) {out[0]=-self.rate*x[0];}
    fn rhs_vjp(&self,_:f64,x:&[f64],b:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        xb[0]=-self.rate*b[0];pb.fill(0.0);pb[0]=-x[0]*b[0];Ok(())
    }
}
impl SensorModel for Model {
    fn initial_values(&self)->&[f64] {&[1.2]}
    fn initial_vjp(&self,_:&[f64],pb:&mut[f64])->Result<(),String> {pb.fill(0.0);Ok(())}
    fn predict(&self,c:u64,_:f64,x:&[f64])->Result<f64,String> {
        if c!=0 {return Err("only physical channel 0 exists before lag augmentation".into());} Ok(x[0])
    }
    fn prediction_vjp(&self,c:u64,_:f64,_:&[f64],b:f64,xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        if c!=0 {return Err("unknown physical channel".into());} xb[0]=b;pb.fill(0.0);Ok(())
    }
}
fn family(data:SensorData)->Result<ObservedFamily<LaggedFamily<Physical>>,String> {
    let sensor=LagSensor{channel:1,source:0,time_constant:LagValue::Parameter(1),
        initial:LagInitial::Value(LagValue::Fixed(0.0)),bias:LagValue::Parameter(2)};
    Ok(ObservedFamily::new(LaggedFamily::new(Physical,&[sensor],0.0,1,2)?,data))
}
fn parse(text:&str)->Result<SensorData,String> {
    if text.len()>MAX_BYTES {return Err("lag CSV exceeds 64 KiB".into());}
    let mut lines=text.lines();
    if lines.next().map(str::trim)!=Some(HEADER) {return Err(format!("expected header {HEADER}"));}
    let mut rows=Vec::new();let mut has_reference=false;let mut has_lag=false;
    for (i,line) in lines.enumerate() {
        if line.trim().is_empty() {continue;}
        if rows.len()==MAX_ROWS {return Err("lag CSV exceeds 256 readings".into());}
        let mut fields=line.split(',').map(str::trim);
        let time:f64=fields.next().ok_or("missing time")?.parse().map_err(|_|"invalid time")?;
        let channel:u64=fields.next().ok_or("missing channel")?.parse().map_err(|_|"invalid channel")?;
        let value:f64=fields.next().ok_or("missing value")?.parse().map_err(|_|"invalid value")?;
        let sigma:f64=fields.next().ok_or("missing sigma")?.parse().map_err(|_|"invalid sigma")?;
        if fields.next().is_some() || time<0.0 || channel>1 {return Err(format!("invalid row {}",i+2));}
        rows.push(SensorReading::new(time,channel,value,sigma,SensorLoss::Quadratic)?);
        has_reference|=channel==0;has_lag|=channel==1;
    }
    if !has_reference || !has_lag {return Err("both reference (0) and lagged (1) readings are required".into());}
    let data=SensorData::new(&rows,MAX_ROWS)?;
    if data.times().windows(2).filter(|p|p[0]!=p[1]).count()<2 {
        return Err("at least three distinct times are required; this is not an identifiability proof".into());
    }
    Ok(data)
}
fn synthetic()->String {
    let mut text=format!("{HEADER}\n");let (k,tau,bias)=(0.7,0.35,0.12);
    for t in [0.0,0.08,0.2,0.5,0.9,1.5,2.5,4.0] {
        let x=1.2*fs_math::det::exp(-k*t);
        // Independent convolution oracle, not the adapter's ODE implementation.
        let z=1.2*(fs_math::det::exp(-k*t)-fs_math::det::exp(-t/tau))/(1.0-k*tau);
        text.push_str(&format!("{t},0,{x:.17},0.05\n{t},1,{:.17},0.05\n",z+bias));
    }
    text
}
fn config(end:f64)->TransientConfig {TransientConfig{start:0.0,initial_step:0.02,
    recording:RecordingConfig{end,rtol:1e-10,atol:1e-12,controller:PiController::default(),max_workspace_components:4096},
    max_state_components:2,max_samples:MAX_ROWS,max_attempts:20000,max_records:20000,
    replay:ReplayBudget{checkpoints:64,replayed_steps:300000},max_kkt_dimension:9}}

pub fn run(args:&[String])->Result<(),Box<dyn std::error::Error>> {
    if args.len()>1 || args.first().is_some_and(|a|a.starts_with("--")) {
        return Err("usage: transient_fit --lag [readings.csv]".into());
    }
    let (text,source)=if let Some(path)=args.first() {
        let mut text=String::new();std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        (text,"csv")
    } else {(synthetic(),"synthetic-noiseless")};
    let data=parse(&text)?;let end=*data.times().last().ok_or("missing observations")?;let f=family(data)?;
    let mut study=TransientStudy::new(&f,&[1.2,0.7,-0.15],config(end),&mut||false)?;
    let before=study.accepted().value;let report=study.run(1e-6,120,1500,&mut||false)?;
    let p=study.optimizer().point();
    println!("source={source} model=decay-plus-first-order-sensor stop={:?} iterations={} evaluations={}",
        report.stop,study.optimizer().iterations(),study.optimizer().evaluations());
    println!("rate_per_s={:.10} sensor_tau_s={:.10} sensor_bias={:.10} objective_before={:.12e} objective_after={:.12e}",p[0],p[1],p[2],before,study.accepted().value);
    println!("declared_initial_physical=1.2 declared_initial_sensor=0; response-lag, not transport-delay; no physical-validation or identifiability claim");
    if report.stop!=SqpStop::Converged {return Err("fit did not meet its local KKT tolerance".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn joint_rate_response_and_bias_fit_uses_the_existing_sqp() {
        let f=family(parse(&synthetic()).unwrap()).unwrap();
        let mut study=TransientStudy::new(&f,&[1.2,0.7,-0.15],config(4.0),&mut||false).unwrap();
        let before=study.accepted().value;
        assert_eq!(study.run(1e-6,120,1500,&mut||false).unwrap().stop,SqpStop::Converged);
        for (got,want) in study.optimizer().point().iter().zip([0.7,0.35,0.12]) {assert!((got-want).abs()<2e-5,"{got} != {want}");}
        assert!(study.accepted().value<before*1e-8);
    }
    #[test]
    fn csv_requires_valid_noise_order_channels_and_a_reference() {
        assert!(parse(&synthetic()).is_ok());
        for rows in ["0,1,0.1,0.1\n1,1,0.3,0.1\n2,1,0.2,0.1", "0,0,1.2,0\n1,1,0.3,0.1\n2,1,0.2,0.1",
            "0,0,1.2,0.1\n2,1,0.2,0.1\n1,1,0.3,0.1", "0,0,NaN,0.1\n1,1,0.3,0.1\n2,1,0.2,0.1"] {
            assert!(parse(&format!("{HEADER}\n{rows}\n")).is_err());
        }
    }
    #[test]
    fn lagged_model_composes_with_fixed_shared_sensor_noise() {
        use fs_ascent::transient::{evaluate_transient, observations::correlated::{CorrelatedFamily,SharedNoiseLimits}};
        let data=parse(&synthetic()).unwrap();let n=data.readings().len();
        let f=CorrelatedFamily::new(family(data).unwrap(),1,&vec![0.02;n],
            SharedNoiseLimits{max_components:4096,max_work:100000,relative_tolerance:1e-10},&mut||false).unwrap();
        let p=[0.9,0.5,-0.1];let got=evaluate_transient(&f,&config(4.0),&p,&mut||false).unwrap().unwrap();
        for i in 0..3 {let (mut a,mut b)=(p,p);a[i]+=1e-5;b[i]-=1e-5;
            let fa=evaluate_transient(&f,&config(4.0),&a,&mut||false).unwrap().unwrap().value;
            let fb=evaluate_transient(&f,&config(4.0),&b,&mut||false).unwrap().unwrap().value;
            assert!((got.gradient[i]-(fa-fb)/2e-5).abs()<2e-5);
        }
    }
}

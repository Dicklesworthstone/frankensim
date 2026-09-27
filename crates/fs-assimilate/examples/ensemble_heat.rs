//! Estimate a 513-node temperature field from sparse thermocouple readings.
//!
//! A manufactured 1-D rod, 1 m long, has fixed 300 K endpoints and diffusivity
//! 0.001 m^2/s. The example's bounded explicit diffusion stencil is a forecast
//! producer, not a new production conduction solver or physical certificate.
//! Eight explicitly specified initial fields span three smooth uncertainty
//! modes; unrepresented uncertainties cannot be recovered or quantified here.
//!
//! Optional CSV: time_s,node,temperature_k,sigma_k (ordered times, zero-based
//! interior node). Sigma denotes independent measurement standard deviation.
//! No input file uses labeled noiseless synthetic readings, not experimental
//! data. Run: cargo run -p fs-assimilate --features ensemble --example ensemble_heat

use fs_assimilate::nonlinear::ensemble::{Ensemble, EnsembleControl, EnsembleError, EnsembleObservation};
use fs_assimilate::nonlinear::ensemble::forecast::ForecastStatus;
use std::io::Read;

const N: usize = 513;
const MEMBERS: usize = 8;
const AMBIENT: f64 = 300.0;
const DIFFUSIVITY: f64 = 0.001;
const MAX_ROWS: usize = 256;

fn field(coefficients: [f64;3]) -> Vec<f64> {
    (0..N).map(|i| AMBIENT + coefficients.iter().enumerate().map(|(k,a)|
        a * (((k+1)*(i+1)) as f64 * std::f64::consts::PI/(N+1) as f64).sin()).sum::<f64>()).collect()
}
fn prior() -> Result<Ensemble, EnsembleError> {
    let mut values=Vec::with_capacity(N*MEMBERS);
    for member in 0..MEMBERS {
        let mut c=[8.0,-3.0,0.0];
        for (k,amplitude) in [4.0,3.0,2.0].into_iter().enumerate() {
            c[k] += if (member>>k)&1 == 0 {amplitude} else {-amplitude};
        }
        values.extend(field(c));
    }
    Ensemble::new(0.0,N,&values,N*MEMBERS)
}

fn heat_step(
    _:usize, start:f64, end:f64, state:&[f64], out:&mut[f64], cancelled:&mut dyn FnMut()->bool,
) -> Result<(),String> {
    if state.len()!=N || out.len()!=N || end<=start {return Err("invalid rod forecast shape/interval".into());}
    let dx=1.0/(N+1) as f64;
    let count=((end-start)*DIFFUSIVITY/(0.45*dx*dx)).ceil();
    if !count.is_finite() || count>100_000.0 {return Err("rod forecast exceeds 100000 substeps".into());}
    let steps=(count as usize).max(1);let q=DIFFUSIVITY*((end-start)/steps as f64)/(dx*dx);
    let mut current=state.to_vec();let mut next=vec![0.0;N];
    for _ in 0..steps {
        if cancelled() {return Err("rod forecast cancelled".into());}
        for i in 0..N {
            let left=if i==0 {AMBIENT} else {current[i-1]};
            let right=if i+1==N {AMBIENT} else {current[i+1]};
            next[i]=q.mul_add((left-current[i])+(right-current[i]),current[i]);
        }
        std::mem::swap(&mut current,&mut next);
    }
    out.copy_from_slice(&current);Ok(())
}

#[derive(Clone,Copy)]
struct Reading { time:f64, node:usize, value:f64, sigma:f64 }
fn parse(text:&str)->Result<Vec<Reading>,String> {
    if text.len()>65_536 {return Err("CSV exceeds 64 KiB".into());}
    let mut lines=text.lines();
    if lines.next().map(str::trim)!=Some("time_s,node,temperature_k,sigma_k") {return Err("invalid CSV header".into());}
    let mut rows:Vec<Reading>=Vec::new();
    for line in lines.filter(|line| !line.trim().is_empty()) {
        if rows.len()==MAX_ROWS {return Err("more than 256 readings".into());}
        let columns:Vec<_>=line.split(',').map(str::trim).collect();
        if columns.len()!=4 {return Err("expected four CSV columns".into());}
        let row=Reading {time:columns[0].parse().map_err(|_|"invalid time")?,
            node:columns[1].parse().map_err(|_|"invalid node")?,value:columns[2].parse().map_err(|_|"invalid temperature")?,
            sigma:columns[3].parse().map_err(|_|"invalid sigma")?};
        if !row.time.is_finite() || row.time<0.0 || row.node>=N || !row.value.is_finite()
            || !row.sigma.is_finite() || row.sigma<=0.0 || rows.last().is_some_and(|last|last.time>row.time)
        {return Err("invalid or out-of-order thermocouple reading".into());}
        rows.push(row);
    }
    if rows.is_empty() {return Err("no thermocouple readings".into());}
    Ok(rows)
}
fn synthetic()->Result<(Vec<Reading>,Vec<f64>),String> {
    let mut truth=field([11.0,-1.0,2.0]);let mut time=0.0;let mut rows=Vec::new();
    for end in [0.25,0.5,1.0] {
        let mut next=vec![0.0;N];heat_step(0,time,end,&truth,&mut next,&mut || false)?;truth=next;time=end;
        for node in [102,257,410] {rows.push(Reading{time,node,value:truth[node],sigma:0.2});}
    }
    Ok((rows,truth))
}
fn estimate(rows:&[Reading], chunk:usize)->Result<(Ensemble,usize),EnsembleError> {
    if chunk == 0 { return Err(EnsembleError::Invalid("forecast chunk must be positive")); }
    let mut state=prior()?;let mut control=EnsembleControl::new(2*MEMBERS*MAX_ROWS,2*N*MEMBERS+N);
    for (id,row) in rows.iter().enumerate() {
        if row.time>state.time() {
            let mut job=state.forecast(row.time,&control,&mut || false)?;
            loop {
                if job.advance(chunk,&mut heat_step,&mut control,&mut || false)?.status==ForecastStatus::Complete {break;}
            }
            state=job.finish()?;
        }
        state.assimilate_scalar(EnsembleObservation{id:id as u64,time:row.time,value:row.value,sigma:row.sigma},
            None,&mut |_,_,x| Ok(x[row.node]),&mut control,&mut || false)?;
    }
    Ok((state,control.model_calls()))
}
fn rmse(a:&[f64],b:&[f64])->f64 { (a.iter().zip(b).map(|(a,b)|(a-b)*(a-b)).sum::<f64>()/a.len() as f64).sqrt() }
fn main()->Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>1 {return Err("usage: ensemble_heat [readings.csv]".into());}
    let (rows,truth,source)=if let Some(path)=args.first() {
        let mut text=String::new();std::fs::File::open(path)?.take(65_537).read_to_string(&mut text)?;
        (parse(&text)?,None,"csv")
    } else {let (rows,truth)=synthetic()?;(rows,Some(truth),"synthetic-noiseless")};
    let (result,calls)=estimate(&rows,2)?;let control=EnsembleControl::new(0,2*N);
    let moments=result.moments(&control,&mut || false)?;
    println!("source={source} state_nodes={N} members={MEMBERS} observations={} time_s={} model_calls={calls}",rows.len(),result.time());
    if let Some(truth)=truth {println!("synthetic_field_rmse_k={:.8}",rmse(&moments.mean,&truth));}
    println!("node,estimated_temperature_k,sample_std_k");
    for node in [0,102,257,410,N-1] {println!("{node},{:.8},{:.8}",moments.mean[node],moments.std[node]);}
    println!("scope=manufactured-rod ensemble estimate; no coverage, identifiability or physical-validation certificate");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_readings_recover_the_manufactured_field_and_split_forecasts_replay() {
        let (rows,truth)=synthetic().unwrap();let (a,calls)=estimate(&rows,MEMBERS).unwrap();let (b,split_calls)=estimate(&rows,2).unwrap();
        assert_eq!(a,b);assert_eq!(calls,split_calls);assert_eq!(calls,96);
        let estimate=a.moments(&EnsembleControl::new(0,2*N),&mut || false).unwrap();
        let mut free=field([8.0,-3.0,0.0]);let mut time=0.0;
        for end in [0.25,0.5,1.0] {let mut next=vec![0.0;N];heat_step(0,time,end,&free,&mut next,&mut || false).unwrap();free=next;time=end;}
        assert!(rmse(&estimate.mean,&truth)<0.05);
        assert!(rmse(&estimate.mean,&truth)<0.02*rmse(&free,&truth));
        assert!(estimate.std.iter().any(|v|*v>0.0));
    }
    #[test]
    fn csv_admission_refuses_missing_invalid_and_reordered_readings() {
        let header="time_s,node,temperature_k,sigma_k\n";
        assert!(parse(&format!("{header}0,102,310,0.2\n0.5,257,309,0.2\n")).is_ok());
        for data in ["", "0,513,310,0.2\n", "0,102,NaN,0.2\n", "0,102,310,0\n", "1,102,310,0.2\n0,102,310,0.2\n"] {
            assert!(parse(&format!("{header}{data}")).is_err());
        }
    }
}

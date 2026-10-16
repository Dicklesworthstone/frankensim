//! Late-delivery history mode of ensemble_heat. Same rod, prior and forecast.
//! CSV arrival_s,time_s,node,temperature_k,sigma_k separates ingestion from
//! acquisition time. A common reference may couple all rows in ONE delivery
//! batch across acquisition times; different batches remain independent.
use super::{N, MEMBERS, MAX_ROWS, field, prior, heat_step, rmse};
use fs_assimilate::nonlinear::ensemble::{EnsembleControl, EnsembleError};
use fs_assimilate::nonlinear::ensemble::forecast::{ForecastStatus, smoother::EnsembleSmoother};
use fs_assimilate::nonlinear::ensemble::correlated::ObservationBlock;
use std::io::Read;

const HEADER: &str = "arrival_s,time_s,node,temperature_k,sigma_k";
const MAX_FRAMES: usize = 32;
const MAX_BLOCK: usize = 32;
#[derive(Clone, Copy)]
struct DelayedReading { arrival: f64, time: f64, node: usize, value: f64, sigma: f64 }

fn validate(rows: &[DelayedReading]) -> Result<(), String> {
    if rows.is_empty() || rows.len() > MAX_ROWS { return Err("expected 1..256 delayed readings".into()); }
    let mut batch = 0;
    for (i,row) in rows.iter().enumerate() {
        if !row.time.is_finite() || !row.arrival.is_finite() || row.time < 0.0 || row.time > row.arrival
            || row.node >= N || !row.value.is_finite() || !row.sigma.is_finite() || row.sigma <= 0.0
            || (i > 0 && rows[i-1].arrival > row.arrival)
        { return Err("finite readings require 0 <= acquisition <= ordered arrival and positive noise".into()); }
        if rows[..i].iter().any(|p| p.time == row.time && p.node == row.node) {
            return Err("a thermocouple sample was delivered more than once".into());
        }
        batch = if i > 0 && rows[i-1].arrival == row.arrival { batch+1 } else { 1 };
        if batch > MAX_BLOCK { return Err("a delivery batch exceeds 32 readings".into()); }
    }
    Ok(())
}
fn parse(text: &str) -> Result<Vec<DelayedReading>, String> {
    if text.len() > 65_536 { return Err("CSV exceeds 64 KiB".into()); }
    let mut lines=text.lines();
    if lines.next().map(str::trim) != Some(HEADER) { return Err(format!("expected CSV header {HEADER}")); }
    let mut rows=Vec::new();
    for line in lines.filter(|l| !l.trim().is_empty()) {
        if rows.len()==MAX_ROWS { return Err("more than 256 delayed readings".into()); }
        let f:Vec<_>=line.split(',').map(str::trim).collect();
        if f.len()!=5 { return Err("expected five delayed-reading columns".into()); }
        rows.push(DelayedReading {arrival:f[0].parse().map_err(|_|"invalid arrival")?,
            time:f[1].parse().map_err(|_|"invalid acquisition time")?,node:f[2].parse().map_err(|_|"invalid node")?,
            value:f[3].parse().map_err(|_|"invalid temperature")?,sigma:f[4].parse().map_err(|_|"invalid sigma")?});
    }
    validate(&rows)?;Ok(rows)
}

fn estimate(rows: &[DelayedReading], frames: usize, chunk: usize, shared: f64)
    -> Result<(EnsembleSmoother, usize), EnsembleError>
{
    validate(rows).map_err(|_|EnsembleError::Invalid("invalid delayed-reading stream"))?;
    if frames==0 || frames>MAX_FRAMES || chunk==0 || !shared.is_finite() || shared<0.0 {
        return Err(EnsembleError::Invalid("frames must be 1..32, positive chunk and nonnegative finite shared noise"));
    }
    // The bounded offline input provides its recording timetable. This does not
    // invent missing past states for a truly live stream: a live caller must
    // have forecast/retained those acquisition times before a late arrival.
    let mut timeline=vec![0.0];
    for row in rows { timeline.push(row.time);timeline.push(row.arrival); }
    timeline.sort_by(|a,b|a.partial_cmp(b).expect("finite admitted time"));
    timeline.dedup_by(|a,b| *a==*b);
    for row in rows {
        let acquisition=timeline.binary_search_by(|t|t.partial_cmp(&row.time).unwrap()).unwrap();
        let arrival=timeline.binary_search_by(|t|t.partial_cmp(&row.arrival).unwrap()).unwrap();
        if arrival-acquisition>=frames { return Err(EnsembleError::Invalid("late sample exceeds --frames retention; increase the window")); }
    }
    let width=N*frames;
    let workspace=2*MEMBERS*(width+MAX_BLOCK)+MEMBERS+3*MAX_BLOCK+width+frames;
    let batches=1+rows.windows(2).filter(|p|p[0].arrival!=p[1].arrival).count();
    let calls=MEMBERS*(timeline.len()-1+batches);
    let mut control=EnsembleControl::new(calls,workspace);
    let mut state=EnsembleSmoother::new(prior()?,frames)?;
    let mut first=0;
    for frontier in timeline {
        if frontier>state.time() {
            let mut job=state.forecast(frontier,&control,&mut||false)?;
            while job.advance(chunk,&mut heat_step,&mut control,&mut||false)?.status!=ForecastStatus::Complete {}
            state=job.finish()?;
        }
        let mut end=first;
        while end<rows.len() && rows[end].arrival==frontier { end+=1; }
        if end==first { continue; }
        let group=&rows[first..end];
        let ids=(first..end).map(|i|i as u64).collect::<Vec<_>>();
        let values=group.iter().map(|r|r.value).collect::<Vec<_>>();
        let sigma=group.iter().map(|r|r.sigma).collect::<Vec<_>>();
        let times=group.iter().map(|r|r.time).collect::<Vec<_>>();
        let block=ObservationBlock::shared_reference(frontier,&ids,&values,&sigma,shared,
            MAX_BLOCK,2*MAX_BLOCK*(MAX_BLOCK+1),&mut||false)?;
        state.assimilate_correlated(&block,&times,&mut|_,states,out,check| {
            for (i,row) in group.iter().enumerate() {
                if check() { return Err("history observation cancelled".into()); }
                out[i]=states.state(i).ok_or("missing acquisition frame")?[row.node];
            }
            Ok(())
        },&mut control,&mut||false)?;
        first=end;
    }
    Ok((state,control.model_calls()))
}
fn synthetic() -> Result<(Vec<DelayedReading>, Vec<f64>), String> {
    let (rows,truth)=super::synthetic()?;
    let rows=rows.into_iter().map(|r|DelayedReading {arrival:if r.time==0.25 {0.5} else {1.0},
        time:r.time,node:r.node,value:r.value,sigma:r.sigma}).collect();
    Ok((rows,truth))
}
fn options(args: &[String]) -> Result<(Option<&str>,usize,f64),String> {
    let (mut path,mut frames,mut shared)=(None,None,None);let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--frames" if frames.is_none() => {
                let value:usize=args.get(i+1).ok_or("--frames requires a count")?.parse().map_err(|_|"invalid frame count")?;
                if value==0 || value>MAX_FRAMES {return Err("--frames must be in 1..32".into());}frames=Some(value);i+=2;
            }
            "--shared-sigma" if shared.is_none() => {
                let value:f64=args.get(i+1).ok_or("--shared-sigma requires kelvin")?.parse().map_err(|_|"invalid shared sigma")?;
                if !value.is_finite() || value<0.0 {return Err("shared sigma must be finite and nonnegative".into());}shared=Some(value);i+=2;
            }
            arg if path.is_none() && !arg.starts_with("--") => {path=Some(arg);i+=1;}
            _ => return Err("usage: ensemble_heat --smooth [readings.csv] [--frames 1..32] [--shared-sigma kelvin]".into()),
        }
    }
    Ok((path,frames.unwrap_or(8),shared.unwrap_or(0.0)))
}
pub fn run(args: &[String]) -> Result<(),Box<dyn std::error::Error>> {
    let (path,frames,shared)=options(args)?;
    let (rows,truth,source)=if let Some(path)=path {
        let mut text=String::new();std::fs::File::open(path)?.take(65_537).read_to_string(&mut text)?;
        (parse(&text)?,None,"csv")
    } else {let (rows,truth)=synthetic()?;(rows,Some(truth),"synthetic-noiseless-late-delivery")};
    let (state,calls)=estimate(&rows,frames,2,shared)?;
    let control=EnsembleControl::new(0,2*N);
    println!("source={source} mode=history-smoothing nodes={N} members={MEMBERS} retained_frames={} observations={} model_calls={calls}",state.times().len(),rows.len());
    println!("noise=independent-plus-reference-per-delivery shared_sigma_k={shared}");
    if let Some(truth)=truth {
        let latest=state.moments_at(state.times().len()-1,&control,&mut||false)?;
        println!("synthetic_latest_rmse_k={:.8}",rmse(&latest.mean,&truth));
        if state.times()[0]==0.0 {
            let initial=state.moments_at(0,&control,&mut||false)?;
            println!("synthetic_initial_rmse_k={:.8}",rmse(&initial.mean,&field([11.0,-1.0,2.0])));
        }
    }
    println!("time_s,node,smoothed_temperature_k,sample_std_k");
    for (frame,time) in state.times().iter().enumerate() {
        let moments=state.moments_at(frame,&control,&mut||false)?;
        for node in [102,257,410] {println!("{time},{node},{:.8},{:.8}",moments.mean[node],moments.std[node]);}
    }
    println!("scope=manufactured-rod numerical smoother; no coverage or physical-validation certificate");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn delayed_stream_revises_initial_field_and_replays_forecast_chunks() {
        let (rows,truth)=synthetic().unwrap();let (a,calls)=estimate(&rows,4,MEMBERS,0.0).unwrap();let (b,split)=estimate(&rows,4,1,0.0).unwrap();
        assert_eq!(a,b);assert_eq!(calls,40);assert_eq!(split,calls);assert_eq!(a.times(),[0.0,0.25,0.5,1.0]);
        let c=EnsembleControl::new(0,2*N);let initial=a.moments_at(0,&c,&mut||false).unwrap();
        let latest=a.moments_at(3,&c,&mut||false).unwrap();
        assert!(rmse(&initial.mean,&field([11.0,-1.0,2.0]))<0.02);assert!(rmse(&latest.mean,&truth)<0.02);
        let (on_time,_)=super::super::synthetic().unwrap();let (filtered,_)=super::super::estimate(&on_time,2).unwrap();
        let filtered=filtered.moments(&c,&mut||false).unwrap();assert!(rmse(&filtered.mean,&latest.mean)<1e-8);
    }
    #[test]
    fn cross_time_shared_noise_retains_larger_history_uncertainty() {
        let (rows,truth)=synthetic().unwrap();let (shared,_)=estimate(&rows,4,2,0.8).unwrap();let (diagonal,_)=estimate(&rows,4,2,0.0).unwrap();
        let c=EnsembleControl::new(0,2*N);let a=shared.moments_at(0,&c,&mut||false).unwrap();let b=diagonal.moments_at(0,&c,&mut||false).unwrap();
        assert!(a.std.iter().sum::<f64>()>b.std.iter().sum::<f64>());
        assert!(rmse(&shared.moments_at(3,&c,&mut||false).unwrap().mean,&truth)<0.3);
    }
    #[test]
    fn delivery_identity_retention_and_option_errors_refuse() {
        for data in ["", "0,1,1,310,0.2", "1,0,1,310,0", "2,1,1,310,0.2\n1,0,2,310,0.2",
            "1,0,1,310,0.2\n2,0,1,310,0.2", "1,0,513,310,0.2"] {
            assert!(parse(&format!("{HEADER}\n{data}\n")).is_err());
        }
        assert!(parse(&format!("{HEADER}\n1,0,1,310,0.2\n1,0.5,2,309,0.2\n")).is_ok());
        let (rows,_)=synthetic().unwrap();assert!(estimate(&rows,1,2,0.0).is_err());
        for args in [vec!["--frames","0"],vec!["--frames","33"],vec!["--shared-sigma","NaN"],vec!["--frames"],vec!["--bad"]] {
            assert!(options(&args.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
        }
    }
}

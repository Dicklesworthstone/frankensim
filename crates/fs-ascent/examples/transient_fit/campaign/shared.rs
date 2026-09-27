//! Shared-reference-error mode; no alternate model, optimizer, or CSV parser.
use super::*;
use fs_ascent::transient::observations::correlated::{CorrelatedFamily, SharedNoiseLimits};

type Prepared = Vec<(u64, CorrelatedFamily<Decay>)>;
fn prepare(data: Vec<Data>, sigma: f64) -> Result<Prepared, Box<dyn std::error::Error>> {
    if !sigma.is_finite() || sigma<0.0 {return Err("invalid shared reference scale".into());}
    let mut prepared=Vec::with_capacity(data.len());
    for row in data {
        let loadings=vec![sigma;row.family.data().readings().len()];
        let limits=SharedNoiseLimits {max_components:100_000,max_work:10_000_000,relative_tolerance:1e-9};
        prepared.push((row.id,CorrelatedFamily::new(row.family,1,&loadings,limits,&mut||false)?));
    }
    Ok(prepared)
}
fn build(data: &[(u64, CorrelatedFamily<Decay>)])
    -> Result<TransientCampaign<'_, CorrelatedFamily<Decay>>, Box<dyn std::error::Error>>
{
    let cases=data.iter().map(|(id,family)|Experiment {id:*id,family,weight:1.0,
        config:config(*family.data().times().last().expect("nonempty checked data"))}).collect();
    Ok(TransientCampaign::new(cases,MAX_EXPERIMENTS,4+5*MAX_EXPERIMENTS)?)
}
pub(super) fn run(data: Vec<Data>, sigma: f64, source: &str) -> Result<(), Box<dyn std::error::Error>> {
    let data=prepare(data,sigma)?;let campaign=build(&data)?;
    let mut control=CampaignControl::new(1500,1500*data.len());
    let mut study=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut control,&mut||false)?;
    let initial=study.accepted().value;let report=study.run(1e-6,100,1500,&mut||false)?;
    println!("source={source} experiments={} covariance=per-experiment-shared-reference shared_sigma={sigma} stop={:?} work={:?}",
        data.len(),report.stop,study.work());
    println!("rate_per_s={:.10} bias={:.10} objective_before={:.12e} objective_after={:.12e}",
        study.optimizer().point()[0],study.optimizer().point()[1],initial,study.accepted().value);
    for e in &study.accepted().experiments {
        println!("experiment={} observations={} objective={:.12e} replayed_steps={}",
            e.id,e.result.observations,e.result.value,e.result.replayed_steps);
    }
    println!("scope=fixed-covariance-numerical-fit; separate experiments have independent reference errors; no posterior or physical-validation claim");
    if report.stop!=SqpStop::Converged {return Err("correlated campaign stopped without satisfying its local KKT tolerance".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shared_csv_campaign_uses_the_existing_sqp_and_joint_pullback() {
        let data=prepare(parse(&synthetic(),SensorLoss::Quadratic).unwrap(),0.25).unwrap();
        let campaign=build(&data).unwrap();let mut work=CampaignControl::new(1500,3000);
        let mut study=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut work,&mut||false).unwrap();
        let report=study.run(1e-6,100,1500,&mut||false).unwrap();
        assert_eq!(report.stop,SqpStop::Converged,"{report:?}");
        for (p,target) in study.optimizer().point().iter().zip([0.7,0.15]) {assert!((p-target).abs()<2e-5);}
        for e in &study.accepted().experiments {assert!(e.result.replayed_steps>=2*e.result.forward.accepted);}
    }
    #[test]
    fn invalid_or_ambiguous_noise_options_refuse_before_fitting() {
        for args in [vec!["--shared-sigma","-1"],vec!["--shared-sigma","NaN"],
            vec!["--shared-sigma","0.2","--huber","1.5"],vec!["--shared-sigma","0.2","--shared-sigma","0.3"]] {
            let args:Vec<_>=args.into_iter().map(str::to_owned).collect();
            assert!(super::super::run(&args).is_err());
        }
    }
}

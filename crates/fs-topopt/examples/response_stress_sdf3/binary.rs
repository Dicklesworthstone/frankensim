//! Binary postprocessing keeps original limits and samples actual physical stress.
use super::*;
use fs_topopt::sdf3::response::limits::binary::BinaryEndpointOptions3;
use fs_topopt::sdf3::response::limits::binary::strength::{BinaryResponseStressOptions3, BinaryStressScanOptions3};
type Error = Box<dyn std::error::Error>;

pub(super) struct Config { thresholds: Vec<f64>, void_stiffness: f64 }
pub(super) fn parse(args: &mut Vec<String>) -> Result<Option<Config>, Error> {
    let Some(index)=args.iter().position(|v|v=="--binary") else {return Ok(None);};
    if index+3!=args.len() {return Err("--binary requires a final comma-separated threshold schedule and explicit void stiffness".into());}
    let count=args[index+1].split(',').count();
    if count>BinaryStressScanOptions3::default().max_candidates {return Err("too many binary thresholds".into());}
    let thresholds:Vec<f64>=args[index+1].split(',').map(str::parse).collect::<Result<_,_>>()?;
    let void_stiffness:f64=args[index+2].parse()?;
    if thresholds.is_empty() || thresholds.iter().any(|v|!v.is_finite()||*v<=0.0||*v>1.0)
        || !void_stiffness.is_finite() || void_stiffness<=0.0 || void_stiffness>=1.0 {
        return Err("binary thresholds must be in (0,1], void stiffness in (0,1)".into());
    }
    args.truncate(index);
    Ok(Some(Config{thresholds,void_stiffness}))
}

pub(super) fn run(study:&mut CutDensityStudy3<AdaptiveSolveSpace3>,rho:&[f64],
    cases:&[ResponseLimitCase3<'_>],options:ResponseStressDesignOptions3,config:Config,
    control:&mut SolveControl<'_>) -> Result<(),Error> {
    let policy=BinaryStressScanOptions3 {
        stress:BinaryResponseStressOptions3 {
            endpoint:BinaryEndpointOptions3 {void_stiffness:config.void_stiffness,limits:options.limits.response,..Default::default()},
            stress_limit:options.limits.stress_limit,aggregation_power:options.limits.stress.aggregation_power,
            max_points:options.limits.stress.max_points,
        },
        feasibility_tolerance:options.optimizer.tolerance,..Default::default()
    };
    eprintln!("binary_thresholds={:?}; binary_void_stiffness={:.17e}; physical_stress_cap={:.17e}; physical_stress_power={}; numerical_tolerance={:.9e}; binary_volume_cap={}; source_material_preserved=true",
        config.thresholds,config.void_stiffness,policy.stress.stress_limit,policy.stress.aggregation_power,
        policy.feasibility_tolerance,policy.volume_cap);
    let report=study.scan_binary_response_stress_limits(rho,cases,&config.thresholds,policy,control)?;
    println!("binary_candidate,threshold,solid_fraction,joint_constraint,worst_response_residual,physical_stress_norm,physical_max,stress_residual,response_feasible,joint_feasible,volume_feasible");
    for (i,r) in report.candidates.iter().enumerate() {
        println!("{i},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{},{}",r.threshold,r.volume_fraction,
            r.constraint,r.worst_response_residual,r.stress_aggregate,r.sampled_physical_max,r.stress_residual,
            r.response_feasible,r.joint_feasible,r.volume_feasible);
    }
    // Output completed evidence before propagating any late failure. An accepted
    // prefix is useful but is not relabeled a successfully completed full scan.
    if let Some(e)=&report.best {
        println!("binary_case,limit,sense,value,bound,scale,residual");
        for (i,c) in cases.iter().enumerate() {for (j,l) in c.limits.iter().enumerate() {
            println!("{i},{j},{:?},{:.17e},{:.17e},{:.17e},{:.17e}",l.sense,e.response.measurements[i][j],l.bound,l.scale,e.response.residuals[i][j]);
        }}
        println!("binary_case,points,physical_max,solid_physical_max,void_physical_max");
        for (i,c) in e.cases.iter().enumerate() {
            println!("{i},{},{:.17e},{:.17e},{:.17e}",c.points,c.sampled_physical_max,c.solid_max,c.void_max);
        }
        println!("binary_cell,solid,stiffness_scale");
        for (i,(&solid,scale)) in e.response.solid.iter().zip(e.response.stiffness_scales()).enumerate() {
            println!("{i},{solid},{scale:.17e}");
        }
        eprintln!("binary_best_index={:?}; binary_threshold={:.17e}; source_gray_volume={:.17e}; binary_solid_volume={:.17e}; points={}; joint_feasible={}",
            report.best_index,e.response.threshold,e.response.source_volume_fraction,e.response.volume_fraction,
            e.point_count,e.feasible(policy.feasibility_tolerance));
    }
    eprintln!("binary_complete={}; binary_candidates_completed={}; work={:?}; continuum_certified=false; zero_void_certified=false; manufacturing_certified=false",
        report.outcome.is_ok(),report.candidates.len(),report.work);
    report.outcome?;
    if report.best.is_none() {return Err("no scheduled binary endpoint satisfies response and physical-stress requirements".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(s:&str)->Vec<String>{s.split_whitespace().map(str::to_owned).collect()}
    #[test]
    fn g0_explicit_binary_schedule_preserves_order_ties_and_legacy_arguments(){
        let mut a=args("80 500000 0.0012 0.008 0.1 --binary 0.7,0.2,0.7 0.001");
        let c=parse(&mut a).unwrap().unwrap();assert_eq!(c.thresholds,[0.7,0.2,0.7]);assert_eq!(c.void_stiffness,0.001);
        assert_eq!(a,args("80 500000 0.0012 0.008 0.1"));
        let before=a.clone();assert!(parse(&mut a).unwrap().is_none());assert_eq!(a,before);
    }
    #[test]
    fn g0_invalid_binary_options_are_refused_without_consuming_arguments(){
        for s in ["--binary","--binary 0.5","--binary 0.5 0","--binary 0.5 1","--binary 0.5 NaN",
            "--binary 0.5,NaN 0.001","--binary 0,0.5 0.001","--binary 1.1 0.001",
            "--binary 0.5, 0.001","--binary 0.5 0.001 --binary 0.6 0.001"] {
            let mut a=args(s);let before=a.clone();assert!(parse(&mut a).is_err(),"{s}");assert_eq!(a,before);
        }
        let mut a=vec!["--binary".to_owned(),vec!["0.5";65].join(","),"0.001".to_owned()];
        assert!(parse(&mut a).is_err());
    }
}

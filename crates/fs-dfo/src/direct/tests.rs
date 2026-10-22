use super::*;
use std::{cell::Cell, convert::Infallible};
fn config(work:usize)->DirectConfig {
    DirectConfig {max_evaluations:work,max_rectangles:2048,max_rectangle_entries:4096,
        max_selection_comparisons:2048*2048,epsilon:1e-4,target:None}
}
fn state(work:usize)->DirectState {DirectState::new_with(&[[0.,1.];2],config(work),&mut ||true).unwrap()}
fn score(x:&[f64])->Result<f64,Infallible>{Ok((x[0]-0.2).powi(2)+(x[1]-0.8).powi(2))}

#[test]
fn g0_positive_slope_intersection_selects_hull_and_all_equal_size_ties() {
    let rectangles:Vec<_>=[(0.1,0.),(0.2,0.3),(0.4,0.2),(0.4,0.5)].into_iter()
        .map(|(radius,value)|Rectangle {center:vec![0.5],half:vec![radius],radius,value}).collect();
    assert_eq!(select(&rectangles,0.,1e-4,&mut ||true).unwrap().into_iter().collect::<Vec<_>>(),vec![0,2]);
    let tied:Vec<_>=[(0.1,0.),(0.2,0.3),(0.4,0.2),(0.4,0.2),(0.4,0.5)].into_iter()
        .map(|(radius,value)|Rectangle {center:vec![0.5],half:vec![radius],radius,value}).collect();
    assert_eq!(select(&tied,0.,1e-4,&mut ||true).unwrap().into_iter().collect::<Vec<_>>(),vec![0,2,3]);
    let flat:Vec<_>=[0.1,0.2,0.4,0.4].into_iter()
        .map(|radius|Rectangle {center:vec![0.5],half:vec![radius],radius,value:4.}).collect();
    assert_eq!(select(&flat,4.,1e-4,&mut ||true).unwrap().into_iter().collect::<Vec<_>>(),vec![2,3]);
}
#[test]
fn g0_subdivision_retains_volume_orders_axes_by_samples_and_publishes_atomically() {
    let mut s=state(3);
    let mut f=|x:&[f64]|Ok::<_,Infallible>(x[0]+10.*x[1]);
    let a=s.advance_with(99,&mut f,&mut ||true).unwrap();
    assert_eq!(a.evaluations,3);assert_eq!(a.rectangles,1);assert_eq!(a.pending_samples,2);
    assert_eq!(s.rectangles[0].half,vec![0.5;2]);
    s.extend_evaluation_budget(5).unwrap();
    let b=s.advance_with(99,&mut f,&mut ||true).unwrap();
    assert_eq!(b.rectangles,5);assert_eq!(b.iterations,1);assert_eq!(b.pending_samples,0);
    assert_eq!(s.rectangles[1].half,vec![0.5,0.5/3.]);
    assert_eq!(s.rectangles[3].half,vec![0.5/3.;2]);
    let volume:f64=s.rectangles.iter().map(|r|r.half.iter().map(|x|2.*x).product::<f64>()).sum();
    assert!((volume-1.).abs()<1e-14);
    for i in 0..s.rectangles.len() {for j in i+1..s.rectangles.len() {
        assert!((0..2).any(|axis|(s.rectangles[i].center[axis]-s.rectangles[j].center[axis]).abs()
            >=s.rectangles[i].half[axis]+s.rectangles[j].half[axis]-1e-14));
    }}
}
#[test]
fn g5_every_budget_boundary_and_pause_partition_preserves_candidate_trace() {
    let mut reference=state(41);let mut expected=Vec::new();
    let final_report=reference.advance_with(41,&mut |x|{expected.push(x.to_vec());score(x)},&mut ||true).unwrap();
    for width in [1,3,8,19] {
        let mut s=state(41);let mut actual=Vec::new();
        let mut report=s.advance_with(0,&mut score,&mut ||true).unwrap();
        while report.stop==DirectStop::Paused {
            report=s.advance_with(width,&mut |x|{actual.push(x.to_vec());score(x)},&mut ||true).unwrap();
        }
        assert_eq!(actual,expected);assert_eq!(report,final_report);
    }
    let mut s=state(0);let mut actual=Vec::new();
    for cap in 0..=41 {
        if cap>0 {s.extend_evaluation_budget(cap).unwrap();}
        let report=s.advance_with(99,&mut |x|{actual.push(x.to_vec());score(x)},&mut ||true).unwrap();
        assert_eq!(report.evaluations,cap);assert_eq!(report.stop,DirectStop::EvaluationLimit);
        assert_eq!(actual,expected[..cap]);
    }
}
#[test]
fn g4_success_survives_cancellation_and_failures_retry_the_original_candidate() {
    let permitted=Cell::new(true);let mut s=state(8);
    assert!(matches!(s.advance_with(1,&mut |x|{permitted.set(false);score(x)},&mut ||permitted.get()),
        Err(DirectRunError::Solver(DirectError::Cancelled))));
    assert_eq!(s.evaluations(),1);
    let report=s.advance_with(0,&mut |_|->Result<f64,Infallible>{panic!("cached root repeated")},&mut ||true).unwrap();
    assert_eq!(report.rectangles,1);assert_eq!(report.evaluations,1);
    let mut attempts=Vec::new();
    assert!(matches!(s.advance_with(1,&mut |x|{attempts.push(x.to_vec());Err::<f64,_>("physics")},&mut ||true),
        Err(DirectRunError::Evaluation("physics"))));
    assert!(s.advance_with(1,&mut |x|{attempts.push(x.to_vec());Ok::<_,Infallible>(f64::NAN)},&mut ||true).is_err());
    s.advance_with(1,&mut |x|{attempts.push(x.to_vec());score(x)},&mut ||true).unwrap();
    assert_eq!(s.evaluations(),4);assert_eq!(s.rectangles.len(),1);
    assert!(attempts.windows(2).all(|pair|pair[0]==pair[1]));
}
#[test]
fn g0_partition_limits_can_extend_and_resolution_is_not_convergence() {
    let mut cfg=config(12);cfg.max_rectangles=1;
    let mut s=DirectState::new_with(&[[0.,1.];2],cfg,&mut ||true).unwrap();
    assert_eq!(s.advance_with(99,&mut score,&mut ||true).unwrap().stop,DirectStop::RectangleLimit);
    assert_eq!(s.evaluations(),1);
    s.extend_partition_budget(32,4096,2048*2048).unwrap();
    assert_eq!(s.advance_with(4,&mut score,&mut ||true).unwrap().evaluations,5);
    let mut cfg=config(12);cfg.max_selection_comparisons=1;
    let mut s=DirectState::new_with(&[[0.,1.];2],cfg,&mut ||true).unwrap();
    assert_eq!(s.advance_with(99,&mut score,&mut ||true).unwrap().stop,DirectStop::SelectionLimit);
    assert_eq!(s.evaluations(),5);
    let hi=f64::from_bits(1.0_f64.to_bits()+1);
    let mut tiny=DirectState::new_with(&[[1.,hi]],config(8),&mut ||true).unwrap();
    assert_eq!(tiny.advance_with(8,&mut |_|Ok::<_,Infallible>(1.),&mut ||true).unwrap().stop,DirectStop::ResolutionLimit);
    assert_eq!(tiny.evaluations(),1);
    assert!(DirectState::new_with(&[],config(8),&mut ||true).is_err());
    assert!(DirectState::new_with(&[[1.,1.]],config(8),&mut ||true).is_err());
    assert!(DirectState::new_with(&[[-f64::MAX,f64::MAX]],config(8),&mut ||true).is_err());
}
#[test]
fn g1_global_exploration_crosses_a_flat_melting_response_without_gradients() {
    let mut cfg=config(800);cfg.target=Some(1e-5);
    let mut s=DirectState::new_with(&[[0.,50.];2],cfg,&mut ||true).unwrap();
    let report=s.advance_with(800,&mut |x|{
        let f1=((80.+x[0]-100.)/100.).clamp(0.,1.);
        let f2=((80.+x[0]+x[1]-100.)/100.).clamp(0.,1.);
        Ok::<_,Infallible>(0.5*((f1-0.1)/0.1).powi(2)+0.5*((f2-0.3)/0.1).powi(2))
    },&mut ||true).unwrap();
    assert_eq!(report.stop,DirectStop::TargetReached);
    let best=report.best.unwrap();assert!((best.coordinates[0]-30.).abs()<0.1);assert!((best.coordinates[1]-20.).abs()<0.1);
    let mut flat=state(63);
    let report=flat.advance_with(63,&mut |_|Ok::<_,Infallible>(5.),&mut ||true).unwrap();
    assert_eq!(report.stop,DirectStop::EvaluationLimit);assert_eq!(report.evaluations,63);
}
#[test]
fn g1_branin_finds_a_global_basin_instead_of_only_polishing_the_midpoint() {
    let mut s=DirectState::new_with(&[[-5.,10.],[0.,15.]],config(1200),&mut ||true).unwrap();
    let pi=std::f64::consts::PI;
    let report=s.advance_with(1200,&mut |x|Ok::<_,Infallible>(
        (x[1]-5.1/(4.*pi*pi)*x[0]*x[0]+5./pi*x[0]-6.).powi(2)
            +10.*(1.-1./(8.*pi))*x[0].cos()+10.),&mut ||true).unwrap();
    assert!(report.best.unwrap().value<0.3981);
}

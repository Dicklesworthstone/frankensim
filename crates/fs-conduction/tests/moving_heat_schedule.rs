//! Authored schedule gradients through the original normalized P1 projection.
mod support;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::moving::sensitivity::MovingHeatGradient;
use fs_conduction::{ConductionMesh, COMPONENT_POWER_DIMS};
use fs_conduction::fixtures::box_grid;
use fs_qty::{Dims, QtyAny};
use fs_scenario::signal::{Interp, TimeSignal};
use support::{with_cx, with_cancelled_cx};
const LENGTH: Dims = Dims([1,0,0,0,0,0]);
fn constant(v:f64,d:Dims)->TimeSignal { TimeSignal::Constant(QtyAny::new(v,d)) }
fn source(p:[f64;11], hold:bool)->MovingHeatSource {
    MovingHeatSource::new([
        TimeSignal::Ramp {t_start:0.0,t_end:2.0,from:QtyAny::new(p[0],LENGTH),to:QtyAny::new(p[1],LENGTH)},
        TimeSignal::Table {times:vec![0.0,1.0,2.0],values:p[2..5].to_vec(),dims:LENGTH,interp:Interp::Linear},
        constant(p[5],LENGTH)], [p[6],p[7],p[8]],
        TimeSignal::Table {times:vec![0.0,1.0],values:p[9..11].to_vec(),dims:COMPONENT_POWER_DIMS,
            interp:if hold {Interp::Hold} else {Interp::Linear}},0.0,
        MovingHeatBudget {max_vertices:1000,max_elements:1000,min_active_vertices:2,power_relative_tolerance:1e-12}).unwrap()
}
fn parameters()->[f64;11] { [0.1,0.7,0.2,0.6,0.4,0.3,0.85,0.8,0.9,2.0,4.0] }
#[test]
fn g3_schedule_value_gradients_match_primal_source_perturbations() {
    let (complex,positions)=box_grid([3;3],[1.0;3]);
    let mesh=ConductionMesh::new(complex,positions).unwrap();
    with_cx(|cx| {
        let p=parameters();let source=source(p,false);
        assert_eq!(source.schedule_control_count().unwrap(),11);
        let seed:Vec<_>=mesh.positions().iter().map(|x|1.0+x[0]+x[1]*x[2]).collect();
        for time in [-1.0,0.25,1.0,1.5,3.0] {
            let sampled=source.linearize_at(cx,&mesh,time).unwrap().pullback(cx,&seed).unwrap();
            let g=source.pullback_schedule_at(cx,time,sampled,11).unwrap();
            let flat:Vec<_>=g.center_values_m.iter().flatten().copied().chain(g.half_width_m)
                .chain(g.absorbed_power_values_w.iter().copied()).collect();
            let objective=|p| {
                let s=self::source(p,false).at(cx,&mesh,time).unwrap();
                seed.iter().enumerate().map(|(i,v)|v*s.source.at(i)).sum::<f64>()
            };
            for i in 0..11 {
                let(mut plus,mut minus)=(p,p);let eps=1e-6;plus[i]+=eps;minus[i]-=eps;
                let fd=(objective(plus)-objective(minus))/(2.0*eps);
                assert!((flat[i]-fd).abs()<3e-6*(1.0+fd.abs()),"time={time}, i={i}, {} != {fd}",flat[i]);
            }
        }
    });
}
#[test]
fn g1_held_knots_route_gradients_to_the_same_side_as_absorbed_power() {
    let source=source(parameters(),true);
    with_cx(|cx| {
        let seed=MovingHeatGradient {center_m:[2.0,3.0,5.0],half_width_m:[7.0,11.0,13.0],absorbed_power_w:17.0};
        let right=source.pullback_schedule_at(cx,1.0,seed,11).unwrap();
        let left=source.pullback_schedule_before(cx,1.0,seed,11).unwrap();
        assert_eq!(right.absorbed_power_values_w,vec![0.0,17.0]);
        assert_eq!(left.absorbed_power_values_w,vec![17.0,0.0]);
        assert_eq!(left.center_values_m,right.center_values_m);
        assert_eq!(left.center_values_m,[vec![1.0,1.0],vec![0.0,3.0,0.0],vec![5.0]]);
        assert_eq!(left.half_width_m,[7.0,11.0,13.0]);
        assert_eq!(source.pullback_schedule_before(cx,1.0_f64.next_up(),seed,11).unwrap().absorbed_power_values_w,
            vec![0.0,17.0]);
    });
}
#[test]
fn g0_schedule_budget_bad_seed_and_cancellation_return_no_partial_gradient() {
    let source=source(parameters(),false);let original=source.clone();
    with_cx(|cx| {
        assert!(source.pullback_schedule_at(cx,0.25,MovingHeatGradient::default(),10).is_err());
        assert!(source.pullback_schedule_at(cx,f64::NAN,MovingHeatGradient::default(),11).is_err());
        let invalid=MovingHeatGradient {absorbed_power_w:f64::NAN,..Default::default()};
        assert!(source.pullback_schedule_before(cx,0.25,invalid,11).is_err());
        assert_eq!(source,original);
    });
    with_cancelled_cx(|cx| assert!(source.pullback_schedule_at(cx,0.25,MovingHeatGradient::default(),11).is_err()));
    with_cx(|cx| assert_eq!(source.pullback_schedule_at(cx,0.25,MovingHeatGradient::default(),11).unwrap()
        .absorbed_power_values_w,vec![0.0,0.0]));
}

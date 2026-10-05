//! Independent member-speed perturbations re-run the production nominal roots.
use fs_airflow::{FanArrangement as A, FanBank, FanCurve, FanPoint, LossElement,
    LossNetwork, LossResistance, LeakageElement, EnclosureNetwork,
    SourceProvenance, ToleranceBasis, solve_operating_point};
use fs_airflow::composite::{compose_parallel, compose_series};
use fs_airflow::composite::speed::{operating_flow_speed_gradient as gradient, SpeedDerivativeError};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Pressure, VolumetricFlowRate as Q};

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 7, kernel_id: 837, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn bank(name: &str, points: &[(f64,f64)], count: usize, arrangement: A, speed: f64) -> FanBank {
    FanBank::new(FanCurve::new(name, points.iter().map(|&(q,p)| FanPoint::new(Q::new(q),Pressure::new(p))).collect(),
        SourceProvenance::new("synthetic piecewise fan fixture",name),0.0,ToleranceBasis::Analytic,
        Q::new(0.0),(0.25,2.0)).unwrap(), count, arrangement, speed).unwrap()
}
fn banks() -> Vec<FanBank> {
    vec![bank("a",&[(0.0,180.0),(0.2,150.0),(0.6,80.0),(1.4,0.0)],2,A::Series,0.9),
        bank("b",&[(0.0,210.0),(0.3,155.0),(0.7,65.0),(1.6,0.0)],3,A::Parallel,1.1),
        bank("c",&[(0.0,150.0),(0.25,120.0),(0.8,45.0),(1.5,0.0)],1,A::Series,1.0)]
}
fn at_speed(bank: &FanBank, speed: f64) -> FanBank {
    FanBank::new(bank.curve().clone(),bank.count(),bank.arrangement(),speed).unwrap()
}
fn network(r: f64) -> EnclosureNetwork {
    let element = |name| LossElement::new(name,LossResistance::new(4.0*r),0.0,
        SourceProvenance::new("analytic quadratic fixture",name),ToleranceBasis::Analytic).unwrap();
    EnclosureNetwork::new(LossNetwork::Element(element("vent")),LeakageElement::new(element("leak")))
}
fn physical(banks: &[FanBank], topology: A, r: f64) -> f64 {
    let fan = match topology {
        A::Series => compose_series(banks), A::Parallel => compose_parallel(banks),
    }.unwrap();
    solve_operating_point(&fan,&network(r)).unwrap().flow.value.value()
}

#[test]
fn independent_bank_speed_derivatives_match_production_series_and_parallel_resolves() {
    context(&CancelGate::new_clock_free(), |cx| {
        for (topology,r) in [(A::Series,400.0),(A::Parallel,20.0)] {
            let bs = banks();
            let q = physical(&bs,topology,r);
            let actual = gradient(cx,&bs,topology,Q::new(q),network(r).equivalent_resistance(),1e-8,12).unwrap();
            assert!(actual.relative_residual < 1e-8);
            assert!((actual.log_flow_per_log_speed.iter().sum::<f64>()-1.0).abs()<1e-8);
            for i in 0..bs.len() {
                let step=2e-5_f64;
                let mut plus=bs.clone(); let mut minus=bs.clone();
                plus[i]=at_speed(&bs[i],bs[i].speed_ratio()*step.exp());
                minus[i]=at_speed(&bs[i],bs[i].speed_ratio()*(-step).exp());
                let expected=(physical(&plus,topology,r).ln()-physical(&minus,topology,r).ln())/(2.0*step);
                let value=actual.log_flow_per_log_speed[i];
                assert!((value-expected).abs()<1e-5*expected.abs().max(1e-3),"{topology:?} bank {i}: {value:e} vs {expected:e}");
                assert!(value>0.0 && value<1.0,"one independent bank is not the whole fan system");
            }
            let reversed: Vec<_>=bs.iter().rev().cloned().collect();
            let other=gradient(cx,&reversed,topology,Q::new(q),LossResistance::new(r),1e-8,12).unwrap();
            assert_eq!(actual.log_flow_per_log_speed,other.log_flow_per_log_speed.into_iter().rev().collect::<Vec<_>>());
            for scale in [0.8,1.2] {
                let scaled:Vec<_>=bs.iter().map(|b| at_speed(b,b.speed_ratio()*scale)).collect();
                assert!((physical(&scaled,topology,r)/q-scale).abs()<1e-8);
            }
        }
    });
}

#[test]
fn fan_speed_derivative_refuses_kinks_endpoints_wrong_primals_and_budget() {
    context(&CancelGate::new_clock_free(), |cx| {
        let bs=banks();let q=physical(&bs,A::Series,400.0);
        assert!(matches!(gradient(cx,&bs,A::Series,Q::new(q*1.1),LossResistance::new(400.0),1e-8,12),
            Err(SpeedDerivativeError::Residual { .. })));
        assert!(matches!(gradient(cx,&bs,A::Series,Q::new(q),LossResistance::new(400.0),1e-8,11),
            Err(SpeedDerivativeError::PointBudget { required:12, allowed:11 })));
        let line=bank("line",&[(0.0,100.0),(0.5,50.0),(1.0,0.0)],1,A::Series,1.0);
        let same=gradient(cx,&[line.clone(),line.clone()],A::Series,Q::new(0.5),LossResistance::new(400.0),1e-10,6).unwrap();
        assert_eq!(same.log_flow_per_log_speed,[0.5,0.5]);
        let kink=bank("kink",&[(0.0,120.0),(0.5,50.0),(1.0,0.0)],1,A::Series,1.0);
        assert!(matches!(gradient(cx,&[kink,line.clone()],A::Series,Q::new(0.5),LossResistance::new(400.0),1e-10,6),
            Err(SpeedDerivativeError::NonSmooth { .. })));
        let endpoint=at_speed(&line,0.25);
        assert!(matches!(gradient(cx,&[endpoint],A::Series,Q::new(0.125),LossResistance::new(200.0),1e-10,3),
            Err(SpeedDerivativeError::NonSmooth { .. })));
        let plateau=bank("plateau",&[(0.0,100.0),(0.4,50.0),(0.6,50.0),(1.0,0.0)],1,A::Series,1.0);
        assert!(matches!(gradient(cx,&[plateau,line],A::Parallel,Q::new(1.0),LossResistance::new(50.0),1e-10,7),
            Err(SpeedDerivativeError::NonSmooth { .. })));
        assert!(gradient(cx,&[],A::Series,Q::new(1.0),LossResistance::new(1.0),1e-10,12).is_err());
        for invalid in [0.0,-1.0,f64::NAN,f64::INFINITY] {
            assert!(gradient(cx,&bs,A::Series,Q::new(invalid),LossResistance::new(400.0),1e-8,12).is_err());
        }
    });
    let gate=CancelGate::new_clock_free();gate.request();
    context(&gate,|cx| assert!(matches!(gradient(cx,&banks(),A::Series,Q::new(1.0),LossResistance::new(400.0),1e-8,12),
        Err(SpeedDerivativeError::Cancelled))));
}

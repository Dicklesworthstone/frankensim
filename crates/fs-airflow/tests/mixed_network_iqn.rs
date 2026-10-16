//! Real hydraulic split/merge topology with independently powered lumped solids.
//! The oracle follows branch energy balances and exponential effectiveness,
//! not a second fixed-point iteration. The legacy relaxation path is a control.
use fs_airflow::{AirflowError, LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_airflow::conjugate::{AirSegment, ConjugateConfig, SolidRegionState};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolution, GraphSolveConfig, LossGraph};
use fs_airflow::graph::thermal::coupled_transport::{
    solve_coupled_transport, solve_coupled_transport_iqn, solve_coupled_transport_iqn_from,
};
use fs_airflow::graph::thermal::transport::{
    BranchThermalModel, TransportAir, TransportConfig, TransportError, TransportInlet, TransportNetwork,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::iqn_ils::IqnIlsConfig;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(
        gate, arena, StreamKey { seed: 19, kernel_id: 731, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic,
    )))
}
fn flow(cx: &Cx<'_>) -> GraphSolution {
    let edge = |name, from, to, resistance| GraphBranch {
        from, to, loss: LossElement::new(name, LossResistance::new(resistance), 0.0,
            SourceProvenance::new("analytic fixture", "network-iqn-v1"), ToleranceBasis::Analytic).unwrap(),
    };
    LossGraph::new(3, vec![edge("left",0,1,1.0/9.0), edge("right",0,1,1.0), edge("tail",1,2,1.0/16.0)])
        .unwrap().solve(&[
            FixedPressure { node:0, pressure:Pressure::new(2.0) },
            FixedPressure { node:2, pressure:Pressure::new(0.0) },
        ], GraphSolveConfig { max_sweeps:4096, max_node_iterations:80,
            absolute_flow_tolerance:VolumetricFlowRate::new(1e-12), relative_flow_tolerance:1e-12 },cx).unwrap()
}
const G: [f64;3] = [150.0,25.0,100.0];
const P: [f64;3] = [10.0,20.0,30.0];
fn network<'a>(cx: &Cx<'_>, flow: &'a GraphSolution) -> TransportNetwork<'a> {
    TransportNetwork::new(cx, flow,
        TransportAir { density:Density::new(1.0), specific_heat_j_kg_k:1.0 },
        G.iter().enumerate().map(|(i,&g)| BranchThermalModel::Exchange(vec![
            AirSegment::new(&format!("wall-{i}"),1.0,g).unwrap(),
        ])).collect(),
        &[TransportInlet { node:0, temperature:Temperature::new(300.0) }],
        TransportConfig { absolute_flow_tolerance:VolumetricFlowRate::new(1e-10), relative_flow_tolerance:1e-10,
            absolute_heat_tolerance_w:1e-7, relative_heat_tolerance:1e-8 }).unwrap()
}
fn config() -> ConjugateConfig {
    ConjugateConfig { max_iterations:12, ..ConjugateConfig::default() }
}
fn solid(refs: &[f64]) -> Vec<SolidRegionState> {
    refs.iter().enumerate().map(|(i,&r)| SolidRegionState {
        region:format!("wall-{i}"), area_m2:1.0,
        mean_wall_temperature_k:r+P[i]/G[i], heat_rate_w:P[i], mean_reference_temperature_k:Some(r),
    }).collect()
}
fn close(a: f64, b: f64) { assert!((a-b).abs()<3e-7, "{a:.15e} != {b:.15e}"); }

#[test]
fn stiff_split_merge_converges_without_increasing_the_solid_solve_budget() {
    with_gate(&CancelGate::new_clock_free(), |cx| {
        let flow=flow(cx); let net=network(cx,&flow);
        assert!(matches!(solve_coupled_transport(cx,&net,&config(),|_,r| Ok(solid(r))),
            Err(TransportError::Airflow(AirflowError::ConjugateNotConverged { .. }))));
        let mut calls=0;
        let result=solve_coupled_transport_iqn(cx,&net,&config(),IqnIlsConfig::default(),|_,r| {
            calls+=1; Ok(solid(r))
        }).unwrap();
        assert_eq!(calls,result.iterations);
        assert!(calls<=config().max_iterations);
        let mixed=300.0+(P[0]+P[1])/4.0;
        let inlets=[300.0,300.0,mixed]; let capacity=[3.0,1.0,4.0];
        for i in 0..3 {
            let expected=inlets[i]+P[i]/(capacity[i]*-(-G[i]/capacity[i]).exp_m1());
            close(result.solid[i].mean_wall_temperature_k,expected);
            assert_eq!(result.solid[i].mean_reference_temperature_k.unwrap().to_bits(),
                result.reference_temperatures_k[i].to_bits());
        }
        close(result.transport.node_temperatures_k[1].unwrap(),mixed);
        close(result.transport.node_temperatures_k[2].unwrap(),315.0);
        close(result.transport.external_heat_gain_w,60.0);
    });
}

#[test]
fn cancellation_retains_all_actual_references_and_warm_start_recloses_physics() {
    let gate=CancelGate::new_clock_free();
    let mut last=Vec::new(); let mut calls=0;
    let error=with_gate(&gate,|cx| {
        let flow=flow(cx); let net=network(cx,&flow);
        solve_coupled_transport_iqn(cx,&net,&config(),IqnIlsConfig::default(),|_,r| {
            last=r.to_vec(); calls+=1; if calls==2 { gate.request(); } Ok(solid(r))
        }).unwrap_err()
    });
    let TransportError::Airflow(AirflowError::Cancelled { iteration,references_k })=error
        else { panic!("complete cancelled interface required") };
    assert_eq!(iteration,1); assert_eq!(references_k,last); assert_eq!(references_k.len(),3);
    with_gate(&CancelGate::new_clock_free(),|cx| {
        let flow=flow(cx); let net=network(cx,&flow);
        let result=solve_coupled_transport_iqn_from(cx,&net,&config(),&references_k,IqnIlsConfig::default(),
            |_,r| Ok(solid(r))).unwrap();
        close(result.transport.external_heat_gain_w,60.0);
        close(result.transport.node_temperatures_k[2].unwrap(),315.0);
    });
}

#[test]
fn malformed_policy_and_corrupt_heat_cannot_bypass_existing_admission() {
    with_gate(&CancelGate::new_clock_free(),|cx| {
        let flow=flow(cx); let net=network(cx,&flow);
        let invalid=IqnIlsConfig { max_history:0, ..IqnIlsConfig::default() };
        assert!(solve_coupled_transport_iqn(cx,&net,&config(),invalid,|_,_| panic!("invalid policy ran solid")).is_err());
        let result=solve_coupled_transport_iqn(cx,&net,&config(),IqnIlsConfig::default(),|_,r| {
            let mut response=solid(r); response[1].heat_rate_w+=0.01; Ok(response)
        });
        assert!(matches!(result,Err(TransportError::Airflow(AirflowError::ConjugateBalanceUnclosed { .. }))));
        let one=ConjugateConfig { max_iterations:1, ..config() };
        assert!(matches!(solve_coupled_transport_iqn(cx,&net,&one,IqnIlsConfig::default(),|_,r| Ok(solid(r))),
            Err(TransportError::Airflow(AirflowError::ConjugateNotConverged { iterations:1, .. }))));
    });
}

#[test]
fn nonphysical_extrapolation_never_reaches_the_solid_callback() {
    // Deliberately unstable synthetic solid response: its positive map samples
    // extrapolate toward a negative fixed point. This tests the safeguard, not
    // the validity of that synthetic constitutive law.
    with_gate(&CancelGate::new_clock_free(),|cx| {
        let flow=flow(cx); let net=network(cx,&flow);
        let config=ConjugateConfig { max_iterations:5, ..config() };
        let mut calls=0;
        let result=solve_coupled_transport_iqn(cx,&net,&config,IqnIlsConfig::default(),|_,r| {
            calls+=1;
            assert!(r.iter().all(|&v| v>0.0 && v.is_finite()));
            Ok(r.iter().enumerate().map(|(i,&reference)| SolidRegionState {
                region:format!("wall-{i}"), area_m2:1.0, mean_wall_temperature_k:2.0*reference,
                heat_rate_w:G[i]*reference, mean_reference_temperature_k:Some(reference),
            }).collect())
        });
        assert_eq!(calls,config.max_iterations);
        assert!(matches!(result,Err(TransportError::Airflow(AirflowError::ConjugateNotConverged { .. }))));
    });
}

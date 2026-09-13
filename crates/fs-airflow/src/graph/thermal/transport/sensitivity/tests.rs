//! Independent closed forms and perturbed production marches, not only J/J^T agreement.
use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_qty::{Density, Temperature, VolumetricFlowRate};

use super::*;
use crate::conjugate::AirSegment;
use crate::graph::{GraphSolution, LossGraph};
use crate::graph::tests::{boundary, bridge, config as flow_config, edge, with_cx};
use crate::graph::thermal::transport::{TransportAir, TransportConfig, TransportInlet};

fn config() -> TransportConfig {
    TransportConfig { absolute_flow_tolerance: VolumetricFlowRate::new(1.0e-8),
        relative_flow_tolerance: 1.0e-8, absolute_heat_tolerance_w: 1.0e-5, relative_heat_tolerance: 1.0e-8 }
}
fn exchange(rows: &[(&str, f64)]) -> BranchThermalModel {
    BranchThermalModel::Exchange(rows.iter().map(|&(name, u)| AirSegment::new(name, 1.0, u).unwrap()).collect())
}
fn inlet(node: usize, temperature: f64) -> TransportInlet {
    TransportInlet { node, temperature: Temperature::new(temperature) }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}");
}
struct Fixture {
    flow: GraphSolution,
    models: Vec<BranchThermalModel>,
    walls: Vec<f64>,
    inlets: Vec<TransportInlet>,
    config: TransportConfig,
}
impl Fixture {
    fn network<'a>(&'a self, cx: &Cx<'_>) -> TransportNetwork<'a> {
        TransportNetwork::new(cx, &self.flow,
            TransportAir { density: Density::new(1.0), specific_heat_j_kg_k: 2.0 },
            self.models.clone(), &self.inlets, self.config).unwrap()
    }
    fn shifted(&self, cx: &Cx<'_>, d: &TransportDirection, step: f64) -> TransportMarch {
        let base = self.network(cx);
        let indices: BTreeMap<&str, usize> = base.regions().into_iter().enumerate().map(|(i, name)| (name, i)).collect();
        let models = self.models.iter().map(|model| match model {
            BranchThermalModel::Adiabatic => BranchThermalModel::Adiabatic,
            BranchThermalModel::Exchange(rows) => BranchThermalModel::Exchange(rows.iter().map(|row| {
                let i = indices[row.region()];
                AirSegment::new(row.region(), row.area_m2(),
                    row.htc_w_per_m2_k() * (step * d.log_conductances[i]).exp()).unwrap()
            }).collect()),
        }).collect();
        let inlets: Vec<_> = self.inlets.iter().map(|source| inlet(source.node,
            source.temperature.value() + step * d.inlets_k[source.node])).collect();
        let walls: Vec<_> = self.walls.iter().zip(&d.walls_k).map(|(wall, delta)| wall + step * delta).collect();
        TransportNetwork::new(cx, &self.flow,
            TransportAir { density: Density::new(1.0), specific_heat_j_kg_k: 2.0 },
            models, &inlets, self.config).unwrap().march(cx, &walls).unwrap()
    }
}
fn split() -> Fixture {
    let graph = LossGraph::new(4, vec![edge("hot", 0, 2, 1.0), edge("bypass", 1, 2, 4.0),
        edge("mixed", 2, 3, 5.0/9.0)]).unwrap();
    Fixture { flow: with_cx(|cx| graph.solve(&[boundary(0, 9.0), boundary(1, 9.0), boundary(3, 0.0)], flow_config(), cx)).unwrap(),
        models: vec![exchange(&[("up-a", 1.2), ("up-b", 2.3)]), BranchThermalModel::Adiabatic, exchange(&[("down", 3.4)])],
        walls: vec![345.0, 320.0, 330.0], inlets: vec![inlet(0, 300.0), inlet(1, 285.0)], config: config() }
}
fn direction(lin: &TransportLinearization<'_, '_>) -> TransportDirection {
    let mut d = lin.zero_direction();
    for (i, delta) in d.walls_k.iter_mut().enumerate() { *delta = 0.7 - 0.3 * i as f64; }
    for (i, delta) in d.log_conductances.iter_mut().enumerate() { *delta = 0.2 + 0.15 * i as f64; }
    for (i, delta) in d.inlets_k.iter_mut().enumerate() {
        if lin.network.inlet_temperatures[i].is_some() { *delta = 0.3 - 0.2 * i as f64; }
    }
    d
}
fn objective(lin: &TransportLinearization<'_, '_>) -> TransportObjective {
    let mut w = lin.zero_objective();
    for (i, value) in w.node_temperatures.iter_mut().enumerate() {
        if lin.primal.node_temperatures_k[i].is_some() { *value = 0.4 - 0.3 * i as f64; }
    }
    for (i, value) in w.branch_outlets.iter_mut().enumerate() {
        if lin.primal.branches[i].outlet_temperature_k.is_some() { *value = -0.3 + 0.2 * i as f64; }
    }
    for (i, value) in w.references.iter_mut().enumerate() { *value = 0.5 - 0.3 * i as f64; }
    w.wall_heat_rate = 0.3; w.external_heat_gain = -0.2;
    w.heat_imbalance = 0.7; w.hydraulic_energy_defect = -0.5;
    w
}
fn value(m: &TransportMarch, w: &TransportObjective) -> f64 {
    m.node_temperatures_k.iter().zip(&w.node_temperatures).map(|(t, w)| t.unwrap_or(0.0) * w).sum::<f64>()
        + m.branches.iter().zip(&w.branch_outlets).map(|(b, w)| b.outlet_temperature_k.unwrap_or(0.0) * w).sum::<f64>()
        + m.reference_temperatures_k.iter().zip(&w.references).map(|(t, w)| t * w).sum::<f64>()
        + m.wall_heat_rate_w * w.wall_heat_rate + m.external_heat_gain_w * w.external_heat_gain
        + m.heat_imbalance_w * w.heat_imbalance + m.hydraulic_energy_defect_w * w.hydraulic_energy_defect
}
fn differential_value(d: &TransportDifferential, w: &TransportObjective) -> f64 {
    d.node_temperatures_k.iter().zip(&w.node_temperatures).map(|(t, w)| t.unwrap_or(0.0) * w).sum::<f64>()
        + d.branch_outlets_k.iter().zip(&w.branch_outlets).map(|(t, w)| t.unwrap_or(0.0) * w).sum::<f64>()
        + d.reference_temperatures_k.iter().zip(&w.references).map(|(t, w)| t * w).sum::<f64>()
        + d.wall_heat_rate_w * w.wall_heat_rate + d.external_heat_gain_w * w.external_heat_gain
        + d.heat_imbalance_w * w.heat_imbalance + d.hydraulic_energy_defect_w * w.hydraulic_energy_defect
}
fn check_finite_difference(f: &Fixture, cx: &Cx<'_>) {
    let net = f.network(cx);
    let lin = net.linearize(cx, &f.walls).unwrap();
    let d = direction(&lin);
    let tangent = lin.apply(cx, &d).unwrap();
    let step = 1.0e-4;
    let plus = f.shifted(cx, &d, step);
    let minus = f.shifted(cx, &d, -step);
    for ((got, a), b) in tangent.node_temperatures_k.iter().zip(&plus.node_temperatures_k).zip(&minus.node_temperatures_k) {
        match (got, a, b) { (Some(g), Some(a), Some(b)) => close(*g, (a-b)/(2.0*step), 2.0e-6),
            (None, None, None) => {}, other => panic!("shape changed: {other:?}") }
    }
    for ((got, a), b) in tangent.branch_outlets_k.iter().zip(&plus.branches).zip(&minus.branches) {
        match (got, a.outlet_temperature_k, b.outlet_temperature_k) {
            (Some(g), Some(a), Some(b)) => close(*g, (a-b)/(2.0*step), 2.0e-6),
            (None, None, None) => {}, other => panic!("shape changed: {other:?}") }
    }
    for ((got, a), b) in tangent.reference_temperatures_k.iter().zip(&plus.reference_temperatures_k).zip(&minus.reference_temperatures_k) {
        close(*got, (a-b)/(2.0*step), 2.0e-6);
    }
    for (got, a, b) in [(tangent.wall_heat_rate_w, plus.wall_heat_rate_w, minus.wall_heat_rate_w),
        (tangent.external_heat_gain_w, plus.external_heat_gain_w, minus.external_heat_gain_w),
        (tangent.heat_imbalance_w, plus.heat_imbalance_w, minus.heat_imbalance_w),
        (tangent.hydraulic_energy_defect_w, plus.hydraulic_energy_defect_w, minus.hydraulic_energy_defect_w)] {
        close(got, (a-b)/(2.0*step), 2.0e-6);
    }
    let w = objective(&lin);
    let grad = lin.pullback(cx, &w).unwrap();
    let dual = grad.walls.iter().zip(&d.walls_k).map(|(a,b)| a*b).sum::<f64>()
        + grad.inlets.iter().zip(&d.inlets_k).map(|(a,b)| a*b).sum::<f64>()
        + grad.log_conductances.iter().zip(&d.log_conductances).map(|(a,b)| a*b).sum::<f64>();
    close(dual, differential_value(&tangent, &w), 2.0e-12);
    // A fresh perturbed production solve checks EVERY adjoint control column.
    for (kind, values) in [&grad.walls, &grad.inlets, &grad.log_conductances].iter().enumerate() {
        for (index, &expected) in values.iter().enumerate() {
            if kind == 1 && lin.network.inlet_temperatures[index].is_none() { assert_eq!(expected, 0.0); continue; }
            let mut basis = lin.zero_direction();
            match kind { 0 => basis.walls_k[index] = 1.0, 1 => basis.inlets_k[index] = 1.0,
                _ => basis.log_conductances[index] = 1.0 }
            let a = value(&f.shifted(cx, &basis, step), &w);
            let b = value(&f.shifted(cx, &basis, -step), &w);
            close(expected, (a-b)/(2.0*step), 3.0e-6);
        }
    }
}

#[test]
fn single_channel_matches_closed_form_derivatives() {
    let graph = LossGraph::new(2, vec![edge("channel", 0, 1, 1.0)]).unwrap();
    let f = Fixture { flow: with_cx(|cx| graph.solve(&[boundary(0, 4.0), boundary(1, 0.0)], flow_config(), cx)).unwrap(),
        models: vec![exchange(&[("wall", 4.0)])], walls: vec![340.0], inlets: vec![inlet(0, 300.0)], config: config() };
    with_cx(|cx| {
        let net = f.network(cx); let lin = net.linearize(cx, &f.walls).unwrap();
        let a = (-1.0_f64).exp();
        close(lin.primal().node_temperatures_k[1].unwrap(), 340.0 - 40.0*a, 1.0e-12);
        let mut d = lin.zero_direction(); d.log_conductances[0] = 1.0;
        let got = lin.apply(cx, &d).unwrap();
        close(got.branch_outlets_k[0].unwrap(), 40.0*a, 1.0e-12);
        close(got.reference_temperatures_k[0], 40.0*(1.0-2.0*a), 1.0e-12);
        close(got.wall_heat_rate_w, 160.0*a, 1.0e-12);
        d.log_conductances[0] = 0.0; d.walls_k[0] = 1.0;
        let got = lin.apply(cx, &d).unwrap();
        close(got.branch_outlets_k[0].unwrap(), 1.0-a, 1.0e-12);
        close(got.wall_heat_rate_w, 4.0*(1.0-a), 1.0e-12);
    });
}

#[test]
fn split_merge_tangent_and_adjoint_match_all_perturbed_production_outputs() {
    let f = split(); with_cx(|cx| check_finite_difference(&f, cx));
}

#[test]
fn reverse_cross_branch_preserves_actual_segment_and_gradient_order() {
    let f = Fixture { flow: with_cx(|cx| bridge().solve(&[boundary(0,25.0), boundary(3,0.0)], flow_config(), cx)).unwrap(),
        models: vec![exchange(&[("feed", 1.0)]), BranchThermalModel::Adiabatic,
            exchange(&[("declared-start", 1.2), ("declared-end", 0.7)]), BranchThermalModel::Adiabatic,
            exchange(&[("out", 2.0)])], walls: vec![345.0, 320.0, 335.0, 330.0],
        inlets: vec![inlet(0, 290.0)], config: config() };
    with_cx(|cx| {
        assert_eq!(f.network(cx).regions(), vec!["feed", "declared-end", "declared-start", "out"]);
        check_finite_difference(&f, cx);
    });
}

#[test]
fn intermediate_reservoir_injection_and_extraction_use_the_mixed_temperature() {
    for middle in [4.0, 8.0] {
        let graph = LossGraph::new(3, vec![edge("a", 0, 1, 1.0), edge("b", 1, 2, 1.0)]).unwrap();
        let flow = with_cx(|cx| graph.solve(&[boundary(0,9.0),boundary(1,middle),boundary(2,0.0)], flow_config(), cx)).unwrap();
        let mut inlets = vec![inlet(0, 310.0)];
        if middle == 8.0 { inlets.push(inlet(1, 285.0)); }
        let f = Fixture { flow, models: vec![exchange(&[("a", 1.0)]),exchange(&[("b",2.0)])],
            walls: vec![345.0, 330.0], inlets, config: config() };
        with_cx(|cx| check_finite_difference(&f, cx));
    }
}

#[test]
fn uniform_temperature_shift_preserves_heat_and_moves_every_reference_once() {
    let f = split();
    with_cx(|cx| {
        let net = f.network(cx); let lin = net.linearize(cx, &f.walls).unwrap();
        let mut d = lin.zero_direction(); d.walls_k.fill(1.0); d.inlets_k[0] = 1.0; d.inlets_k[1] = 1.0;
        let got = lin.apply(cx, &d).unwrap();
        for value in got.node_temperatures_k.iter().chain(&got.branch_outlets_k).flatten() { close(*value, 1.0, 1.0e-13); }
        for value in got.reference_temperatures_k { close(value, 1.0, 1.0e-13); }
        close(got.wall_heat_rate_w, 0.0, 1.0e-12); close(got.external_heat_gain_w, 0.0, 1.0e-12);
    });
}

#[test]
fn small_ntu_reference_sensitivity_does_not_cancel_to_zero() {
    for x in [1.0e-16_f64, 1.0e-12, 1.0e-6] {
        let k = segment_jacobian(x, -(-x).exp_m1(), x, 1.0, 20.0).unwrap();
        assert!(k.reference[1] > 0.0 && k.reference[2] > 0.0);
        close(k.reference[1] / x, 0.5 - x/6.0 + x*x/24.0, 1.0e-13);
        close(k.reference[2] / (20.0*x), 0.5 - x/3.0 + x*x/8.0, 1.0e-13);
    }
    let saturated = segment_jacobian(1000.0, 1.0, 1000.0, 1.0, 20.0).unwrap();
    assert_eq!(saturated.outlet[2], 0.0); close(saturated.reference[2], 0.02, 1.0e-14);
}

#[test]
fn hydraulic_residual_and_moving_enthalpy_reference_are_not_corrected_away() {
    let graph = LossGraph::new(3, vec![edge("a",0,1,1.0),edge("b",1,2,1.0)]).unwrap();
    let mut flow = with_cx(|cx| graph.solve(&[boundary(0,2.0),boundary(2,0.0)],flow_config(),cx)).unwrap();
    flow.branches[1].flow = VolumetricFlowRate::new(1.01);
    let f = Fixture { flow, models: vec![exchange(&[("up",1.0)]), exchange(&[("down",2.0)])],
        walls: vec![340.0,330.0], inlets: vec![inlet(0,300.0)],
        config: TransportConfig { absolute_flow_tolerance: VolumetricFlowRate::new(0.1), absolute_heat_tolerance_w: 100.0, ..config() } };
    with_cx(|cx| {
        let net = f.network(cx); let lin = net.linearize(cx,&f.walls).unwrap();
        assert!(lin.primal().heat_imbalance_w.abs() > 0.1);
        check_finite_difference(&f,cx);
    });
}

#[test]
fn malformed_derivative_requests_and_unknown_stagnant_temperatures_refuse() {
    let graph = LossGraph::new(4,vec![edge("a",0,1,1.0),edge("stagnant",2,3,1.0)]).unwrap();
    let f = Fixture { flow: with_cx(|cx| graph.solve(&[boundary(0,1.0),boundary(1,0.0),boundary(2,0.0)],flow_config(),cx)).unwrap(),
        models:vec![exchange(&[("wall",1.0)]),BranchThermalModel::Adiabatic], walls:vec![330.0],
        inlets:vec![inlet(0,300.0)], config:config() };
    with_cx(|cx| {
        let net=f.network(cx); let lin=net.linearize(cx,&f.walls).unwrap();
        let mut d=lin.zero_direction(); d.walls_k.clear(); assert!(lin.apply(cx,&d).is_err());
        d=lin.zero_direction(); d.inlets_k[1]=1.0; assert!(lin.apply(cx,&d).is_err());
        for bad in [f64::NAN,f64::INFINITY,f64::NEG_INFINITY] {
            d=lin.zero_direction(); d.log_conductances[0]=bad; assert!(lin.apply(cx,&d).is_err());
            let mut w=lin.zero_objective(); w.wall_heat_rate=bad; assert!(lin.pullback(cx,&w).is_err());
        }
        let mut w=lin.zero_objective(); w.node_temperatures[2]=1.0; assert!(lin.pullback(cx,&w).is_err());
        w=lin.zero_objective(); w.branch_outlets[1]=1.0; assert!(lin.pullback(cx,&w).is_err());
        let zero=lin.apply(cx,&lin.zero_direction()).unwrap();
        assert_eq!(zero.node_temperatures_k[2],None); assert_eq!(zero.branch_outlets_k[1],None);
        assert!(net.linearize(cx,&[f64::NAN]).is_err());
    });
}

#[test]
fn cancellation_refuses_primal_tangent_and_adjoint_without_partial_output() {
    let f=split(); let gate=CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:7,kernel_id:712,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let net=f.network(&cx); let lin=net.linearize(&cx,&f.walls).unwrap(); gate.request();
        assert!(matches!(net.linearize(&cx,&f.walls),Err(TransportError::Interrupted)));
        assert!(matches!(lin.apply(&cx,&lin.zero_direction()),Err(TransportError::Interrupted)));
        assert!(matches!(lin.pullback(&cx,&lin.zero_objective()),Err(TransportError::Interrupted)));
    });
}

#[test]
fn linearization_and_sweeps_are_repeatable_and_equal_temperature_controls_vanish() {
    let mut f=split(); f.walls.fill(300.0); for source in &mut f.inlets {source.temperature=Temperature::new(300.0);}
    with_cx(|cx| {
        let net=f.network(cx); let a=net.linearize(cx,&f.walls).unwrap(); let b=net.linearize(cx,&f.walls).unwrap();
        let d=direction(&a); let w=objective(&a);
        assert_eq!(a.primal(),b.primal()); assert_eq!(a.apply(cx,&d).unwrap(),b.apply(cx,&d).unwrap());
        assert_eq!(a.pullback(cx,&w).unwrap(),b.pullback(cx,&w).unwrap());
        assert!(a.pullback(cx,&w).unwrap().log_conductances.iter().all(|g| *g==0.0));
    });
}

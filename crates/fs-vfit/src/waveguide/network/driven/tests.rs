//! G0/G3: independently supplied physical ports retain one wave-state update.
use super::*;
use crate::waveguide::network::{NetworkSegment, multiport::MimoTerminalSpec};
use crate::mimo::MimoModel;

fn graph(last: NetworkNode) -> WaveguideNetwork {
    WaveguideNetwork::new(&[NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Driven, last], &[
        NetworkSegment { nodes: [0, 1], one_way_samples: 2, impedance_pa_s_m3: 2.0 },
        NetworkSegment { nodes: [1, 2], one_way_samples: 3, impedance_pa_s_m3: 3.0 },
        NetworkSegment { nodes: [1, 3], one_way_samples: 5, impedance_pa_s_m3: 7.0 },
    ], 0.01, 1 << 20).unwrap()
}

#[test]
fn independent_ports_exchange_work_through_one_retained_wave_state() {
    let mut net = graph(NetworkNode::Driven);
    assert_eq!(net.driven_nodes().collect::<Vec<_>>(), [2, 3]);
    for i in 0..80 {
        let before = net.stored_energy_j();
        let waves = [DrivenWave { node: 2, outgoing_pressure_pa: if i < 5 { 1.3 } else { 0.0 } },
            DrivenWave { node: 3, outgoing_pressure_pa: if i < 3 { -0.7 } else { 0.0 } }];
        let frame = net.step_with_driven(if i < 4 { 0.8 } else { 0.0 }, &waves).unwrap();
        let actual: f64 = [2, 3].into_iter().map(|n| {
            let port = net.node_frame(n).unwrap();
            -port.pressure_pa * (port.net_flow_into_node_m3_s * 0.01)
        }).sum();
        assert_eq!(frame.driven_work_j, actual);
        assert_eq!(frame.storage_change_j, net.stored_energy_j() - before);
        assert!(frame.balance_residual_j().abs() < 1e-14);
        assert_eq!(net.accepted_steps(), i + 1);
    }
}

#[test]
fn refused_or_previewed_assignments_preserve_every_accepted_boundary() {
    let mut net = graph(NetworkNode::Driven);
    let waves = [DrivenWave { node: 2, outgoing_pressure_pa: 1.0 },
        DrivenWave { node: 3, outgoing_pressure_pa: -2.0 }];
    net.step_with_driven(3.0, &waves).unwrap();
    let state = |n: &WaveguideNetwork| (n.accepted_steps(), n.stored_energy_j(),
        n.incoming_pressure_pa(), *n.node_frame(2).unwrap(), *n.node_frame(3).unwrap(),
        n.driven_boundary(2).unwrap(), n.driven_boundary(3).unwrap());
    let original = state(&net);
    net.preview_step_with_driven(1.0, &waves).unwrap();
    assert_eq!(state(&net), original);
    assert!(net.step(0.0).is_err());
    assert!(net.step_with_driven(0.0, &waves[..1]).is_err());
    assert!(net.step_with_driven(0.0, &[waves[1], waves[0]]).is_err());
    assert!(net.step_with_driven(0.0, &[waves[0], waves[0]]).is_err());
    assert!(net.step_with_driven(0.0, &[waves[0], DrivenWave { node: 3,
        outgoing_pressure_pa: f64::NAN }]).is_err());
    assert!(net.set_terminal_reflection(2, 1.0).is_err());
    assert_eq!(state(&net), original);
}

#[test]
fn driven_ports_and_mutual_radiation_keep_their_separate_work_supplies() {
    let mut net = graph(NetworkNode::Termination { reflection: -1.0 })
        .with_mimo_terminals(MimoTerminalSpec { nodes: vec![3],
            model: MimoModel { ports: 1, poles: vec![], residues: vec![], d: vec![4.0],
                e: vec![0.0], reciprocal: true },
            band_rad_s: [1.0, core::f64::consts::PI / 0.01], sampled_passivity_points: 65,
            max_states: 4, max_memory_bytes: 1 << 20 }).unwrap();
    let mut radiation_work = 0.0;
    for i in 0..60 {
        let f = net.step_with_driven(0.0, &[DrivenWave { node: 2,
            outgoing_pressure_pa: if i < 4 { 2.0 } else { 0.0 } }]).unwrap();
        assert!(f.balance_residual_j().abs() < 1e-14);
        radiation_work += f.coupled_terminal_work_j;
    }
    assert!(radiation_work > 0.0);
}

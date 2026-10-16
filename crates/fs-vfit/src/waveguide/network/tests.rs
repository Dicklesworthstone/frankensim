use super::*;

fn segment(a: usize, b: usize, delay: usize, z: f64) -> NetworkSegment {
    NetworkSegment { nodes: [a, b], one_way_samples: delay, impedance_pa_s_m3: z }
}
fn construct(nodes: &[NetworkNode], edges: &[NetworkSegment]) -> WaveguideNetwork {
    WaveguideNetwork::new(nodes, edges, 1.0 / 48000.0, 1 << 20).unwrap()
}
fn snapshot(net: &WaveguideNetwork) -> Vec<u64> {
    let mut bytes = Vec::new();
    for line in &net.lines {
        bytes.push(line.head as u64);
        bytes.extend(line.waves.iter().chain(&line.energy).map(|v| v.to_bits()));
    }
    for n in &net.observed {
        bytes.extend([n.pressure_pa, n.net_flow_into_node_m3_s, n.absorbed_energy_j].map(f64::to_bits));
    }
    bytes
}
fn check_storage_and_power(net: &WaveguideNetwork, before: f64, frame: NetworkFrame) {
    let direct: f64 = net.lines.iter().map(|line| {
        line.waves.iter().map(|p| p * p).sum::<f64>()
            * line.spec.time_step_s / line.spec.impedance_pa_s_m3
    }).sum();
    let scale = before.abs() + frame.stored_energy_j.abs() + frame.inlet_work_j.abs()
        + frame.terminal_loss_j + f64::MIN_POSITIVE;
    assert!((direct - frame.stored_energy_j).abs() <= 2e-13 * scale);
    assert!(frame.balance_residual_j().abs() <= 2e-13 * scale, "{frame:?}");
    assert!(frame.junction_residual_j.abs() <= 2e-13 * scale);
    assert!(frame.terminal_loss_j >= 0.0);
    for (kind, node) in net.nodes.iter().zip(&net.observed) {
        if matches!(kind, NetworkNode::Termination { .. }) {
            let work = node.pressure_pa * node.net_flow_into_node_m3_s * net.time_step_s;
            assert!((work - node.absorbed_energy_j).abs() <= 2e-13 * scale);
        }
    }
}

#[test]
fn one_section_preserves_the_existing_waveguide_without_direction_bias() {
    for endpoints in [[0, 1], [1, 0]] {
        for r in [-1.0, -0.7, 0.0, 0.4, 1.0] {
            let mut graph = construct(&[NetworkNode::Inlet, NetworkNode::Termination { reflection: r }],
                &[NetworkSegment { nodes: endpoints, one_way_samples: 3, impedance_pa_s_m3: 1e6 }]);
            let mut line = PassiveWaveguide::new(WaveguideSpec {
                one_way_samples: 3, impedance_pa_s_m3: 1e6, time_step_s: 1.0 / 48000.0,
                reflection: r, max_memory_bytes: 1 << 20,
            }).unwrap();
            for n in 0..100 {
                let input = f64::from(n % 17) - 8.0;
                let a = graph.step(input).unwrap();
                let b = line.step(input).unwrap();
                for (x, y) in [
                    (a.incoming_pressure_pa, b.incoming_pressure_pa),
                    (a.inlet_pressure_pa, b.inlet_pressure_pa),
                    (a.inlet_flow_m3_s, b.inlet_flow_m3_s),
                    (a.stored_energy_j, b.stored_energy_j),
                    (a.storage_change_j, b.storage_change_j),
                    (a.inlet_work_j, b.inlet_work_j),
                    (a.terminal_loss_j, b.terminal_loss_j),
                    (graph.node_frame(1).unwrap().pressure_pa, b.terminal_pressure_pa),
                ] { assert_eq!(x.to_bits(), y.to_bits()); }
            }
        }
    }
}

#[test]
fn matched_sections_do_not_insert_an_extra_junction_sample() {
    let mut net = construct(&[NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: 0.5 }],
        &[segment(0, 1, 2, 1e6), segment(1, 2, 5, 1e6)]);
    for n in 0..20 {
        let f = net.step(if n == 0 { 20.0 } else { 0.0 }).unwrap();
        assert_eq!(f.incoming_pressure_pa, if n == 14 { 10.0 } else { 0.0 });
        assert_eq!(net.node_frame(2).unwrap().pressure_pa, if n == 7 { 30.0 } else { 0.0 });
    }
    assert_eq!(net.stored_energy_j(), 0.0);
}

#[test]
fn area_step_matches_analytical_reflection_and_transmission() {
    let (z1, z2) = (1e6, 4e6);
    let mut net = construct(&[NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: 0.0 }],
        &[segment(0, 1, 3, z1), segment(1, 2, 5, z2)]);
    let reflection = (z2 - z1) / (z2 + z1);
    for n in 0..20 {
        let before = net.stored_energy_j();
        let f = net.step(if n == 0 { 10.0 } else { 0.0 }).unwrap();
        check_storage_and_power(&net, before, f);
        if n == 3 {
            assert!((net.node_frame(1).unwrap().pressure_pa - 10.0 * (1.0 + reflection)).abs() < 1e-13);
            assert!(net.node_frame(1).unwrap().net_flow_into_node_m3_s.abs() < 1e-20);
        }
        assert!((f.incoming_pressure_pa - if n == 6 { 10.0 * reflection } else { 0.0 }).abs() < 1e-13);
        assert!((net.node_frame(2).unwrap().pressure_pa
            - if n == 8 { 10.0 * (1.0 + reflection) } else { 0.0 }).abs() < 1e-13);
    }
}

#[test]
fn branches_scatter_by_admittance_and_terminal_controls_retain_waves() {
    let mut net = construct(&[NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: 0.0 }, NetworkNode::Termination { reflection: 1.0 }],
        &[segment(0, 1, 2, 1e6), segment(1, 2, 3, 2e6), segment(1, 3, 5, 3e6)]);
    for n in 0..512 {
        if n == 17 {
            let before = snapshot(&net);
            net.set_terminal_reflection(3, -1.0).unwrap();
            assert_eq!(snapshot(&net), before);
        }
        let before = net.stored_energy_j();
        let f = net.step(if n == 0 { 10.0 } else { 0.0 }).unwrap();
        check_storage_and_power(&net, before, f);
        if n == 2 {
            let expected = 2.0 * (10.0 / 1e6) / (1.0 / 1e6 + 1.0 / 2e6 + 1.0 / 3e6);
            let junction = net.node_frame(1).unwrap();
            assert!((junction.pressure_pa - expected).abs() < 1e-13);
            assert!(junction.net_flow_into_node_m3_s.abs() < 1e-20);
        }
    }
}

#[test]
fn cyclic_network_uses_simultaneous_scattering_and_actual_energy() {
    let mut net = construct(&[NetworkNode::Inlet, NetworkNode::Junction, NetworkNode::Junction,
        NetworkNode::Junction, NetworkNode::Termination { reflection: -0.6 }],
        &[segment(0, 1, 2, 1e6), segment(1, 2, 3, 2e6), segment(2, 3, 4, 3e6),
          segment(3, 1, 5, 1.5e6), segment(3, 4, 6, 2.5e6)]);
    for n in 0..2000 {
        let before = net.stored_energy_j();
        let f = net.step(if n < 1000 { f64::from(n % 31) - 15.0 } else { 0.0 }).unwrap();
        check_storage_and_power(&net, before, f);
        if n >= 1000 { assert!(f.stored_energy_j <= before * (1.0 + 2e-13)); }
    }
}

#[test]
fn failed_preview_or_control_cannot_advance_any_section() {
    let mut net = construct(&[NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: 0.5 }],
        &[segment(0, 1, 2, 1e6), segment(1, 2, 3, 2e6)]);
    net.step(30.0).unwrap();
    let before = snapshot(&net);
    let expected = net.preview_step(7.0).unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(net.step(bad).is_err());
        assert!(net.set_terminal_reflection(2, bad).is_err());
    }
    assert!(net.set_terminal_reflection(0, 0.5).is_err());
    assert!(net.set_terminal_reflection(1, 0.5).is_err());
    assert!(net.set_terminal_reflection(99, 0.5).is_err());
    assert_eq!(snapshot(&net), before);
    assert_eq!(net.step(7.0).unwrap(), expected);
}

#[test]
fn topology_and_whole_network_memory_are_admitted() {
    let nodes = [NetworkNode::Inlet, NetworkNode::Termination { reflection: 0.5 }];
    let edges = [segment(0, 1, 3, 1e6)];
    let bytes = WaveguideNetwork::required_memory_bytes(2, &edges).unwrap();
    assert!(WaveguideNetwork::new(&nodes, &edges, 1e-5, bytes).is_ok());
    assert!(WaveguideNetwork::new(&nodes, &edges, 1e-5, bytes - 1).is_err());
    for edge in [segment(0, 0, 3, 1e6), segment(0, 2, 3, 1e6), segment(0, 1, 0, 1e6),
        segment(0, 1, usize::MAX, 1e6), segment(0, 1, 3, f64::NAN)] {
        assert!(WaveguideNetwork::new(&nodes, &[edge], 1e-5, usize::MAX).is_err());
    }
    for second in [NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: 1.01 }] {
        assert!(WaveguideNetwork::new(&[nodes[0], second], &edges, 1e-5, bytes).is_err());
    }
    let disconnected = [nodes[0], nodes[1], nodes[1], nodes[1]];
    assert!(WaveguideNetwork::new(&disconnected,
        &[edges[0], segment(2, 3, 4, 1e6)], 1e-5, 1 << 20).is_err());
}

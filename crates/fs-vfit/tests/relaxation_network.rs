use fs_vfit::impedance::{ImpedanceState, SeriesImpedanceSpec};
use fs_vfit::relaxation::{RelaxationImpedance, RelaxationImpedanceSpec, RelaxationTerm};
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};
fn spec() -> RelaxationImpedanceSpec {
    RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
        resistance_pa_s_m3: 0.4, inertance_pa_s2_m3: 0.3, compliance_m3_pa: Some(0.7),
    }, &[RelaxationTerm { resistance_pa_s_m3: 0.6, rate_per_s: 2.0 },
         RelaxationTerm { resistance_pa_s_m3: 0.3, rate_per_s: 40.0 }]).unwrap()
}
fn network(nodes: [usize; 2]) -> WaveguideNetwork {
    WaveguideNetwork::new(&[NetworkNode::Inlet, NetworkNode::Relaxation { load: spec() }],
        &[NetworkSegment { nodes, one_way_samples: 3, impedance_pa_s_m3: 1.0 }], 0.02, 1 << 20).unwrap()
}
#[test]
fn relaxation_terminal_matches_direct_load_and_two_actual_delay_buffers() {
    for orientation in [[0, 1], [1, 0]] {
        let mut net = network(orientation);
        let mut load = RelaxationImpedance::new(spec(), 1.0, 0.02).unwrap();
        let mut outward = [0.0_f64; 3];
        let mut returning = [0.0_f64; 3];
        for n in 0..2000 {
            let input = if n < 500 { (n as f64).sin() } else { 0.0 };
            let incoming = returning[n % 3];
            let direct = load.step(outward[n % 3]).unwrap();
            outward[n % 3] = input;
            returning[n % 3] = direct.port.reflected_pressure_pa;
            let before = net.stored_energy_j();
            let f = net.step(input).unwrap();
            let independent_wave = outward.iter().chain(&returning).map(|x| x * x * 0.02).sum::<f64>();
            let energy = independent_wave + direct.port.stored_energy_j;
            let scale = (before + energy + f.inlet_work_j.abs() + f.terminal_loss_j).max(f64::MIN_POSITIVE);
            assert_eq!(f.incoming_pressure_pa.to_bits(), incoming.to_bits());
            assert_eq!(net.node_frame(1).unwrap().pressure_pa.to_bits(), direct.port.pressure_pa.to_bits());
            assert_eq!(net.terminal_relaxation_flows(1).unwrap(), direct.branch_flows_m3_s());
            assert!((f.stored_energy_j - energy).abs() < 1e-12 * scale);
            assert!(f.balance_residual_j().abs() < 1e-12 * scale);
            if n >= 500 { assert!(f.stored_energy_j <= before + 1e-12 * scale); }
        }
    }
}
#[test]
fn terminal_previews_and_refused_control_never_publish_partial_history() {
    let mut a = network([0, 1]);
    let mut b = network([0, 1]);
    for _ in 0..12 { a.step(20.0).unwrap(); b.step(20.0).unwrap(); }
    let old = a.terminal_relaxation_flows(1).unwrap().to_vec();
    let base = a.terminal_state(1).unwrap();
    let energy = a.stored_energy_j();
    a.preview_step(7.0).unwrap();
    for invalid in [f64::NAN, f64::INFINITY, f64::MAX] { assert!(a.step(invalid).is_err()); }
    assert!(a.set_terminal_reflection(1, 0.0).is_err());
    assert_eq!(a.terminal_relaxation_flows(1).unwrap(), old);
    assert_eq!(a.terminal_state(1).unwrap(), base);
    assert_eq!(a.stored_energy_j().to_bits(), energy.to_bits());
    assert_eq!(a.step(7.0).unwrap(), b.step(7.0).unwrap());
}
#[test]
fn relaxation_node_memory_and_degree_are_admitted_before_stepping() {
    let nodes = [NetworkNode::Inlet, NetworkNode::Relaxation { load: spec() }];
    let edges = [NetworkSegment { nodes: [0, 1], one_way_samples: 3, impedance_pa_s_m3: 1.0 }];
    let bytes = WaveguideNetwork::required_memory_bytes(2, &edges).unwrap();
    assert!(WaveguideNetwork::new(&nodes, &edges, 0.02, bytes).is_ok());
    assert!(WaveguideNetwork::new(&nodes, &edges, 0.02, bytes - 1).is_err());
    let extra = NetworkSegment { nodes: [1, 2], ..edges[0] };
    assert!(WaveguideNetwork::new(&[nodes[0], nodes[1], NetworkNode::Termination { reflection: 0.0 }],
        &[edges[0], extra], 0.02, 1 << 20).is_err());
    let net = network([0, 1]);
    assert_eq!(net.terminal_state(1), Some(ImpedanceState::default()));
    assert!(net.terminal_relaxation_flows(0).is_none());
}

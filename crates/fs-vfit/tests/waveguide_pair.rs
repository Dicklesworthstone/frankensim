use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};

fn network(delay: usize, reflection: f64, dt: f64) -> WaveguideNetwork {
    WaveguideNetwork::new(&[NetworkNode::Inlet, NetworkNode::Termination { reflection }],
        &[NetworkSegment { nodes: [0, 1], one_way_samples: delay, impedance_pa_s_m3: 2e6 }],
        dt, 100_000).unwrap()
}

#[test]
fn paired_commit_matches_the_original_single_network_arithmetic() {
    let dt = 1.0 / 48_000.0;
    let (mut a, mut b) = (network(3, 0.7, dt), network(7, -0.8, dt));
    let (mut ra, mut rb) = (network(3, 0.7, dt), network(7, -0.8, dt));
    for n in 0..512 {
        let x = [f64::from(n % 13) - 6.0, 3.0 - f64::from(n % 7)];
        let actual = a.step_pair(&mut b, x).unwrap();
        assert_eq!(actual, [ra.step(x[0]).unwrap(), rb.step(x[1]).unwrap()]);
        assert_eq!(a.node_frame(1), ra.node_frame(1));
        assert_eq!(b.node_frame(1), rb.node_frame(1));
    }
}

#[test]
fn failed_second_network_does_not_publish_first_preview() {
    let dt = 1.0 / 48_000.0;
    let (mut a, mut b) = (network(3, 0.7, dt), network(7, -0.8, dt));
    let (mut ra, mut rb) = (network(3, 0.7, dt), network(7, -0.8, dt));
    a.step(2.0).unwrap(); b.step(-3.0).unwrap();
    ra.step(2.0).unwrap(); rb.step(-3.0).unwrap();
    for x in [[7.0, f64::MAX], [f64::NAN, 2.0], [7.0, f64::INFINITY]] {
        assert!(a.step_pair(&mut b, x).is_err());
        assert_eq!(a.stored_energy_j(), ra.stored_energy_j());
        assert_eq!(b.stored_energy_j(), rb.stored_energy_j());
        assert_eq!(a.node_frame(0), ra.node_frame(0));
        assert_eq!(b.node_frame(0), rb.node_frame(0));
    }
    // Future returned waves test the actual retained histories, not only energy.
    for _ in 0..32 {
        assert_eq!(a.step_pair(&mut b, [0.0; 2]).unwrap(),
            [ra.step(0.0).unwrap(), rb.step(0.0).unwrap()]);
    }
    let mut other_clock = network(3, 0.7, 2.0*dt);
    assert!(a.step_pair(&mut other_clock, [0.0; 2]).is_err());
}

//! Reactive endpoint integration: propagation, boundary state and total storage.
use fs_vfit::impedance::{SeriesImpedance, SeriesImpedanceSpec};
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};

fn load(r: f64, l: f64, c: Option<f64>) -> SeriesImpedanceSpec {
    SeriesImpedanceSpec { resistance_pa_s_m3: r, inertance_pa_s2_m3: l, compliance_m3_pa: c }
}
fn edge(a: usize, b: usize, n: usize, z: f64) -> NetworkSegment {
    NetworkSegment { nodes: [a,b], one_way_samples: n, impedance_pa_s_m3: z }
}

#[test]
fn reactive_reflection_returns_after_exact_transit_in_both_orientations() {
    for ends in [[0,1], [1,0]] {
        let spec = load(0.4, 0.3, Some(0.7));
        let mut net = WaveguideNetwork::new(&[NetworkNode::Inlet, NetworkNode::Impedance { load: spec }],
            &[NetworkSegment { nodes: ends, one_way_samples: 3, impedance_pa_s_m3: 1.0 }], 0.02, 1<<20).unwrap();
        let mut boundary = SeriesImpedance::new(spec, 1.0, 0.02).unwrap();
        let (mut forward, mut backward) = ([0.0_f64;3], [0.0_f64;3]);
        for n in 0..128 {
            let a = if n == 0 { 10.0 } else { 0.0 };
            let terminal = boundary.step(forward[n%3]).unwrap();
            let incoming = backward[n%3];
            let f = net.step(a).unwrap();
            assert_eq!(f.incoming_pressure_pa.to_bits(), incoming.to_bits());
            assert_eq!(net.node_frame(1).unwrap().pressure_pa.to_bits(), terminal.pressure_pa.to_bits());
            assert_eq!(net.terminal_state(1), Some(terminal.state));
            if n < 6 { assert_eq!(incoming, 0.0); }
            if n == 6 { assert!(incoming.abs() > 1.0); }
            forward[n%3] = a;
            backward[n%3] = terminal.reflected_pressure_pa;
            let wave = forward.iter().chain(&backward).map(|p| p*p).sum::<f64>() * 0.02;
            let scale = (wave + terminal.stored_energy_j).max(f64::MIN_POSITIVE);
            assert!((f.stored_energy_j - wave - terminal.stored_energy_j).abs() < 1e-12 * scale);
        }
    }
}

#[test]
fn multiple_reactive_loads_close_network_storage_and_return_energy() {
    let specs = [load(0.0, 0.3, Some(0.7)), load(0.4, 0.0, Some(0.6))];
    let mut net = WaveguideNetwork::new(&[NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.8 },
        NetworkNode::Impedance { load: specs[0] }, NetworkNode::Impedance { load: specs[1] }],
        &[edge(0,1,3,1.0), edge(1,2,5,2.0), edge(1,3,7,1.5), edge(1,4,4,0.5)],
        0.02, 1<<20).unwrap();
    let (mut saw_return, mut missing_storage_error) = (false, 0.0_f64);
    for n in 0..3000 {
        let before = net.stored_energy_j();
        let old_load = net.load_stored_energy_j();
        let f = net.step(if n<500 { f64::from(n%17)-8.0 } else { 0.0 }).unwrap();
        let scale = (before+f.stored_energy_j+f.inlet_work_j.abs()+f.terminal_loss_j).max(f64::MIN_POSITIVE);
        assert!(f.balance_residual_j().abs() <= 1e-12*scale);
        assert_eq!(net.stored_energy_j(), f.wave_stored_energy_j+f.load_stored_energy_j);
        let mut independent_load = 0.0;
        for (i, spec) in specs.iter().enumerate() {
            let s = net.terminal_state(i+3).unwrap();
            independent_load += 0.5*spec.inertance_pa_s2_m3*s.inertive_flow_m3_s.powi(2)
                + 0.5*spec.compliance_m3_pa.unwrap()*s.compliance_pressure_pa.powi(2);
            let node = net.node_frame(i+3).unwrap();
            let work = node.pressure_pa*node.net_flow_into_node_m3_s*0.02;
            assert!((node.storage_change_j+node.absorbed_energy_j-work).abs() <= 1e-12*scale);
            assert!(node.absorbed_energy_j >= 0.0);
            saw_return |= work < 0.0;
        }
        assert!((independent_load-f.load_stored_energy_j).abs() <= 1e-12*scale);
        missing_storage_error = missing_storage_error.max((f.storage_change_j
            - (f.load_stored_energy_j-old_load)+f.terminal_loss_j-f.inlet_work_j).abs());
        if n>=500 { assert!(f.stored_energy_j <= before+1e-12*scale); }
    }
    assert!(saw_return);
    assert!(missing_storage_error > 1e-5, "omitting boundary storage must fail the balance");
}

#[test]
fn network_reflection_matches_independent_impedance_and_travel_phase() {
    for bin in [7,23,83,173] {
        let theta = 2.0*core::f64::consts::PI*f64::from(bin)/2048.0;
        let mut net = WaveguideNetwork::new(&[NetworkNode::Inlet,
            NetworkNode::Impedance { load: load(0.4,0.3,Some(0.7)) }], &[edge(0,1,7,1.0)],0.02,1<<20).unwrap();
        let (mut re,mut im)=(0.0,0.0);
        for n in 0..12288 {
            let phase=theta*f64::from(n);
            let f=net.step(phase.cos()).unwrap();
            if n>=10240 {
                re+=f.incoming_pressure_pa*phase.cos()/1024.0;
                im-=f.incoming_pressure_pa*phase.sin()/1024.0;
            }
        }
        let omega=100.0*(0.5*theta).tan();
        let x=omega*0.3-1.0/(omega*0.7);
        let den=1.4*1.4+x*x;
        let zr=(0.4*0.4-1.0+x*x)/den;
        let zi=2.0*x/den;
        let phase=-14.0*theta;
        let expected_re=zr*phase.cos()-zi*phase.sin();
        let expected_im=zr*phase.sin()+zi*phase.cos();
        assert!((re-expected_re).hypot(im-expected_im)<1e-10);
    }
}

#[test]
fn reactive_state_survives_preview_failed_scattering_and_incompatible_control() {
    let nodes=[NetworkNode::Inlet, NetworkNode::Impedance { load:load(0.4,0.3,Some(0.7)) }];
    let mut a=WaveguideNetwork::new(&nodes,&[edge(0,1,3,1.0)],0.02,1<<20).unwrap();
    let mut b=WaveguideNetwork::new(&nodes,&[edge(0,1,3,1.0)],0.02,1<<20).unwrap();
    for n in 0..16 { a.step(f64::from(n)).unwrap(); b.step(f64::from(n)).unwrap(); }
    let state=a.terminal_state(1);
    let observation=*a.node_frame(1).unwrap();
    let energy=a.stored_energy_j();
    a.preview_step(7.0).unwrap();
    assert!(a.set_terminal_reflection(1,0.0).is_err());
    for input in [f64::NAN,f64::INFINITY,f64::MAX] { assert!(a.step(input).is_err()); }
    assert_eq!(a.terminal_state(1),state);
    assert_eq!(*a.node_frame(1).unwrap(),observation);
    assert_eq!(a.stored_energy_j(),energy);
    for n in 0..32 {
        let fa=a.step(f64::from(n%7)-3.0).unwrap();
        let fb=b.step(f64::from(n%7)-3.0).unwrap();
        assert_eq!(fa,fb);
        assert_eq!(a.terminal_state(1),b.terminal_state(1));
    }
}

#[test]
fn reactive_records_are_memory_admitted_and_invalid_loads_are_rejected() {
    let nodes=[NetworkNode::Inlet,NetworkNode::Impedance { load:load(0.4,0.3,Some(0.7)) }];
    let edges=[edge(0,1,3,1.0)];
    let bytes=WaveguideNetwork::required_memory_bytes(2,&edges).unwrap();
    assert!(WaveguideNetwork::new(&nodes,&edges,0.02,bytes).is_ok());
    assert!(WaveguideNetwork::new(&nodes,&edges,0.02,bytes-1).is_err());
    let invalid=[NetworkNode::Inlet,NetworkNode::Impedance { load:load(-1.0,0.3,Some(0.7)) }];
    assert!(WaveguideNetwork::new(&invalid,&edges,0.02,bytes).is_err());
}

//! Public-API circuit and energy checks for interior acoustic elements.
use fs_math::c64::C64;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationTerm};
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment, WaveguideNetwork};

fn load(r: f64, l: f64, c: Option<f64>, relaxation: bool) -> RelaxationImpedanceSpec {
    let terms = [RelaxationTerm { resistance_pa_s_m3: 0.6, rate_per_s: 2.0 },
        RelaxationTerm { resistance_pa_s_m3: 0.3, rate_per_s: 40.0 }];
    RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
        resistance_pa_s_m3: r, inertance_pa_s2_m3: l, compliance_m3_pa: c,
    }, if relaxation { &terms } else { &[] }).unwrap()
}
fn edges(reverse: bool) -> [NetworkSegment; 2] {
    [NetworkSegment { nodes: if reverse { [1, 0] } else { [0, 1] },
        one_way_samples: 3, impedance_pa_s_m3: 1.0 },
     NetworkSegment { nodes: if reverse { [2, 1] } else { [1, 2] },
        one_way_samples: 5, impedance_pa_s_m3: 2.0 }]
}
fn model(node: NetworkNode, reverse: bool) -> WaveguideNetwork {
    WaveguideNetwork::new(&[NetworkNode::Inlet, node,
        NetworkNode::Termination { reflection: 0.0 }], &edges(reverse), 0.02, 1 << 20).unwrap()
}

#[test]
fn series_resistance_has_the_analytical_pressure_jump_and_no_extra_delay() {
    for reverse in [false, true] {
        let mut net = model(NetworkNode::Series { load: load(0.5, 0.0, None, false) }, reverse);
        let u = 20.0 / 3.5;
        for n in 0..20 {
            let f = net.step(if n == 0 { 10.0 } else { 0.0 }).unwrap();
            assert!((f.incoming_pressure_pa - if n == 6 { 10.0 - u } else { 0.0 }).abs() < 1e-13);
            assert!((net.node_frame(2).unwrap().pressure_pa
                - if n == 8 { 2.0 * u } else { 0.0 }).abs() < 1e-13);
            if n == 3 {
                let node = net.node_frame(1).unwrap();
                assert!((node.load_flow_m3_s - u).abs() < 1e-13);
                assert!((node.load_pressure_pa - 0.5 * u).abs() < 1e-13);
                assert!((node.pressure_pa - node.series_other_pressure_pa.unwrap()
                    - node.load_pressure_pa).abs() < 1e-13);
                assert!(node.net_flow_into_node_m3_s.abs() < 1e-13);
                assert!((f.interior_loss_j - 0.5 * u * u * 0.02).abs() < 1e-13);
            }
        }
    }
}

#[test]
fn three_way_shunt_conserves_flow_and_counts_its_loss_once() {
    let mut sections = edges(false).to_vec();
    sections.push(NetworkSegment { nodes: [1, 3], one_way_samples: 7, impedance_pa_s_m3: 0.5 });
    let mut net = WaveguideNetwork::new(&[NetworkNode::Inlet,
        NetworkNode::Shunt { load: load(0.4, 0.0, None, false) },
        NetworkNode::Termination { reflection: 0.0 }, NetworkNode::Termination { reflection: 0.0 }],
        &sections, 0.02, 1 << 20).unwrap();
    for n in 0..4 { net.step(if n == 0 { 10.0 } else { 0.0 }).unwrap(); }
    let node = *net.node_frame(1).unwrap();
    let p = 20.0 / (1.0 + 0.5 + 2.0 + 1.0 / 0.4);
    assert!((node.pressure_pa - p).abs() < 1e-13);
    assert!((node.net_flow_into_node_m3_s - p / 0.4).abs() < 1e-13);
    assert!((node.load_flow_m3_s - p / 0.4).abs() < 1e-13);
    assert!((node.absorbed_energy_j - p * p / 0.4 * 0.02).abs() < 1e-13);
    assert_eq!(node.series_other_pressure_pa, None);
}

#[test]
fn interior_histories_obey_constitutive_laws_and_close_total_storage() {
    let spec = load(0.4, 0.3, Some(0.7), true);
    let mut loss_control = 0.0_f64;
    let mut storage_control = 0.0_f64;
    for kind in [NetworkNode::Series { load: spec }, NetworkNode::Shunt { load: spec }] {
        for reverse in [false, true] {
            let mut net = model(kind, reverse);
            for n in 0..2000 {
                let before = net.stored_energy_j();
                let state0 = net.terminal_state(1).unwrap();
                let flows0 = net.terminal_relaxation_flows(1).unwrap().to_vec();
                let f = net.step(if n < 500 { f64::from(n % 17) - 8.0 } else { 0.0 }).unwrap();
                let node = *net.node_frame(1).unwrap();
                let state = net.terminal_state(1).unwrap();
                let branches = net.terminal_relaxation_flows(1).unwrap();
                let scale = before + f.stored_energy_j + f.inlet_work_j.abs()
                    + f.terminal_loss_j + f.interior_loss_j + f64::MIN_POSITIVE;
                assert!(f.balance_residual_j().abs() <= 2e-12 * scale);
                assert!(f.interior_loss_j >= 0.0);
                let flow = node.load_flow_m3_s;
                let mut expected = 0.4 * flow + 0.3 * (state.inertive_flow_m3_s - state0.inertive_flow_m3_s) / 0.02
                    + 0.5 * (state.compliance_pressure_pa + state0.compliance_pressure_pa);
                let mut storage = 0.15 * state.inertive_flow_m3_s.powi(2)
                    + 0.35 * state.compliance_pressure_pa.powi(2);
                for (i, term) in spec.terms().iter().enumerate() {
                    let mid = 0.5 * (flows0[i] + branches[i]);
                    expected += term.resistance_pa_s_m3 * (flow - mid);
                    storage += 0.5 * term.resistance_pa_s_m3 / term.rate_per_s * branches[i].powi(2);
                    assert!(((branches[i] - flows0[i]) / 0.02
                        - term.rate_per_s * (flow - mid)).abs() < 2e-11 * (1.0 + flow.abs()));
                }
                assert!((expected - node.load_pressure_pa).abs() < 2e-11 * (1.0 + expected.abs()));
                assert!((storage - node.stored_energy_j).abs() < 2e-12 * scale);
                assert!((node.storage_change_j + node.absorbed_energy_j
                    - node.load_pressure_pa * flow * 0.02).abs() < 2e-12 * scale);
                assert_eq!(f.interior_loss_j, node.absorbed_energy_j);
                loss_control = loss_control.max((f.balance_residual_j() - f.interior_loss_j).abs());
                storage_control = storage_control.max((f.balance_residual_j() - node.storage_change_j).abs());
                if n >= 500 { assert!(f.stored_energy_j <= before + 2e-12 * scale); }
            }
        }
    }
    assert!(loss_control > 1e-3 && storage_control > 1e-3, "omitted-state/loss controls must fail");
}

#[test]
fn complex_reflection_matches_independent_series_and_parallel_impedance() {
    let spec = load(0.4, 0.3, Some(0.7), true);
    for series in [false, true] {
        for bin in [7, 23, 83, 173] {
            let kind = if series { NetworkNode::Series { load: spec } } else { NetworkNode::Shunt { load: spec } };
            let mut net = model(kind, false);
            let theta = 2.0 * core::f64::consts::PI * f64::from(bin) / 1024.0;
            let mut response = C64::ZERO;
            for n in 0..8192 {
                let angle = theta * f64::from(n);
                let f = net.step(angle.cos()).unwrap();
                if n >= 7168 {
                    response = response + C64::new(angle.cos(), -angle.sin())
                        .scale(f.incoming_pressure_pa / 512.0);
                }
            }
            let s = C64::new(0.0, 100.0 * (0.5 * theta).tan());
            let mut z = C64::from_re(0.4) + s.scale(0.3) + s.scale(0.7).recip();
            for term in spec.terms() {
                z = z + s.scale(term.resistance_pa_s_m3)
                    * (s + C64::from_re(term.rate_per_s)).recip();
            }
            let zin = if series { z + C64::from_re(2.0) }
                else { (z.recip() + C64::from_re(0.5)).recip() };
            let expected = (zin - C64::ONE) * (zin + C64::ONE).recip()
                * C64::new((6.0 * theta).cos(), -(6.0 * theta).sin());
            assert!((response - expected).abs() < 2e-10, "series={series}, bin={bin}");
        }
    }
}

#[test]
fn preview_and_late_refusal_preserve_interior_and_wave_state() {
    for kind in [NetworkNode::Series { load: load(0.4, 0.3, Some(0.7), true) },
        NetworkNode::Shunt { load: load(0.4, 0.3, Some(0.7), true) }] {
        let mut net = model(kind, false);
        let mut twin = model(kind, false);
        for n in 0..16 { net.step(f64::from(n)).unwrap(); twin.step(f64::from(n)).unwrap(); }
        let before = net.stored_energy_j();
        let state = net.terminal_state(1);
        let flows = net.terminal_relaxation_flows(1).unwrap().to_vec();
        let node = *net.node_frame(1).unwrap();
        net.preview_step(7.0).unwrap();
        for input in [f64::NAN, f64::INFINITY, f64::MAX] { assert!(net.step(input).is_err()); }
        assert!(net.set_terminal_reflection(1, 0.0).is_err());
        assert_eq!(net.stored_energy_j(), before);
        assert_eq!(net.terminal_state(1), state);
        assert_eq!(net.terminal_relaxation_flows(1).unwrap(), flows);
        assert_eq!(*net.node_frame(1).unwrap(), node);
        for _ in 0..32 { assert_eq!(net.step(0.0).unwrap(), twin.step(0.0).unwrap()); }
    }
}

#[test]
fn interior_degrees_equivalent_impedances_and_memory_are_admitted() {
    let shunt = NetworkNode::Shunt { load: load(0.4, 0.3, Some(0.7), false) };
    let series = NetworkNode::Series { load: load(0.4, 0.3, Some(0.7), false) };
    for kind in [shunt, series] {
        let nodes = [NetworkNode::Inlet, kind, NetworkNode::Termination { reflection: 0.0 }];
        let bytes = WaveguideNetwork::required_memory_bytes(nodes.len(), &edges(false)).unwrap();
        assert!(WaveguideNetwork::new(&nodes, &edges(false), 0.02, bytes).is_ok());
        assert!(WaveguideNetwork::new(&nodes, &edges(false), 0.02, bytes - 1).is_err());
        assert!(WaveguideNetwork::new(&[NetworkNode::Inlet, kind], &edges(false)[..1], 0.02, bytes).is_err());
    }
    let nodes = [NetworkNode::Inlet, series, NetworkNode::Termination { reflection: 0.0 },
        NetworkNode::Termination { reflection: 0.0 }];
    let mut three = edges(false).to_vec();
    three.push(NetworkSegment { nodes: [1, 3], one_way_samples: 1, impedance_pa_s_m3: 1.0 });
    assert!(WaveguideNetwork::new(&nodes, &three, 0.02, 1 << 20).is_err());
    let mut huge = edges(false);
    for edge in &mut huge { edge.impedance_pa_s_m3 = f64::MAX; }
    assert!(WaveguideNetwork::new(&nodes[..3], &huge, 0.02, 1 << 20).is_err());
}

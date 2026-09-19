//! Synthetic physics and public coupling tests; no measured-instrument claims.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, ApertureTerminal, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{
    ApertureNetwork, ApertureNetworkFrame, NetworkNode, TubeNetworkSpec, TubeSection,
};
use fs_couple::bernoulli_aperture::tube::{ApertureTube, TubeDrive, UniformTubeSpec};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

const DT: f64 = 1e-5;
fn section(a: usize, b: usize, delay: usize, radius: f64) -> TubeSection {
    TubeSection { nodes: [a, b], length_m: delay as f64 * (343.0 * DT), radius_m: radius,
        max_length_error_m: 1e-15 }
}
fn spec(reflection: f64) -> TubeNetworkSpec {
    TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
            NetworkNode::Termination { reflection: -0.8 }, NetworkNode::Termination { reflection }],
        sections: vec![section(0, 1, 8, 0.007), section(1, 2, 12, 0.009), section(1, 3, 5, 0.003)],
        sound_speed_m_s: 343.0, max_wave_memory_bytes: 1 << 20,
    }
}
fn aperture(z: f64, max_steps: u64, state: ApertureState) -> DynamicAperture {
    let mechanics = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35,
        density_kg_m3: 1.2, impedance_pa_s_m3: z, time_step_s: DT, max_steps,
    };
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic network contact fixture".into()).unwrap().with_internal_loss(5.0).unwrap();
    DynamicAperture::new(mechanics, state, contact).unwrap()
}
fn initial() -> ApertureState { ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 } }
fn model(s: TubeNetworkSpec, steps: u64, state: ApertureState) -> ApertureNetwork {
    ApertureNetwork::new(aperture(s.inlet_impedance(1.2).unwrap(), steps, state), s).unwrap()
}
fn bits(f: ApertureNetworkFrame) -> Vec<u64> {
    let a = f.aperture;
    let n = f.network;
    let mut values = vec![a.step];
    values.extend([a.time_s, a.state.opening_m, a.state.opening_velocity_m_s, a.midpoint_opening_m,
        a.outgoing_pressure_pa, a.bore_pressure_pa, a.jet_flow_m3_s, a.swept_flow_m3_s,
        a.bore_flow_m3_s, a.flow_residual_m3_s, a.stored_energy_j, a.storage_change_j,
        a.dissipated_energy_j, a.pressure_work_j, n.incoming_pressure_pa, n.inlet_pressure_pa,
        n.inlet_flow_m3_s, n.stored_energy_j, n.storage_change_j, n.inlet_work_j,
        n.terminal_loss_j, n.junction_residual_j, f.stored_energy_j, f.storage_change_j,
        f.dissipated_energy_j, f.upstream_work_j, f.body_work_j].map(f64::to_bits));
    values
}

#[test]
fn one_section_network_reproduces_the_existing_coupled_tube() {
    let uniform = UniformTubeSpec {
        length_m: 8.0 * (343.0 * DT), radius_m: 0.007, sound_speed_m_s: 343.0,
        terminal_reflection: -0.8, max_length_error_m: 1e-15, max_wave_memory_bytes: 1 << 20,
    };
    let z = uniform.characteristic_impedance(1.2).unwrap();
    let mut tube = ApertureTube::new(aperture(z, 128, initial()), uniform).unwrap();
    let mut s = spec(0.0);
    s.nodes = vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.8 }];
    s.sections = vec![section(0, 1, 8, 0.007)];
    let mut graph = model(s, 128, initial());
    for i in 0..128 {
        let drive = TubeDrive { upstream_pressure_pa: if i < 50 { 800.0 } else { 0.0 }, body_flow_m3_s: 2e-7 };
        let a = graph.step(drive).unwrap();
        let b = tube.step(drive).unwrap();
        assert_eq!(a.aperture, b.aperture);
        for (x, y) in [(a.stored_energy_j, b.stored_energy_j), (a.storage_change_j, b.storage_change_j),
            (a.dissipated_energy_j, b.dissipated_energy_j), (a.upstream_work_j, b.upstream_work_j),
            (a.body_work_j, b.body_work_j),
            (graph.node_frame(1).unwrap().pressure_pa, b.waveguide.terminal_pressure_pa)] {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }
}

#[test]
fn branched_feedback_closes_total_energy_with_a_separate_body_source() {
    let mut omitted_body_work: f64 = 0.0;
    for reflection in [-1.0, -0.7, 0.0, 0.6, 1.0] {
        for state in [initial(), ApertureState { opening_m: -1e-4, opening_velocity_m_s: -0.2 }] {
            let mut m = model(spec(reflection), 512, state);
            for n in 0..512 {
                let before = m.stored_energy_j();
                let f = m.step(TubeDrive {
                    upstream_pressure_pa: if n < 128 { 1200.0 } else { 0.0 },
                    body_flow_m3_s: if n < 256 { 2e-7 * (f64::from(n) * 0.11).sin() } else { 0.0 },
                }).unwrap();
                let scale = before + f.stored_energy_j + f.dissipated_energy_j
                    + f.upstream_work_j.abs() + f.body_work_j.abs();
                assert!(f.balance_residual_j().abs() <= 3e-10 * scale.max(f64::MIN_POSITIVE));
                assert!(f.dissipated_energy_j >= 0.0);
                assert_eq!(m.stored_energy_j().to_bits(), f.stored_energy_j.to_bits());
                omitted_body_work = omitted_body_work.max((f.storage_change_j + f.dissipated_energy_j - f.upstream_work_j).abs());
                if n >= 256 { assert!(f.stored_energy_j <= before + 3e-10 * scale); }
            }
        }
    }
    assert!(omitted_body_work > 1e-12, "uncredited body source must fail the balance");
}

#[test]
fn terminal_switch_preserves_motion_and_changes_feedback_after_physical_transit() {
    let mut fixed = model(spec(1.0), 256, initial());
    let mut switched = model(spec(1.0), 256, initial());
    let drive = TubeDrive { upstream_pressure_pa: 800.0, body_flow_m3_s: 0.0 };
    for _ in 0..128 { assert_eq!(bits(fixed.step(drive).unwrap()), bits(switched.step(drive).unwrap())); }
    let before = switched.stored_energy_j();
    let state = switched.aperture().state();
    switched.set_terminal_reflection(3, -1.0).unwrap();
    assert_eq!(switched.stored_energy_j().to_bits(), before.to_bits());
    assert_eq!(switched.aperture().state(), state);
    assert_eq!(switched.aperture().accepted_steps(), 128);
    assert_eq!(switched.spec().nodes[3], NetworkNode::Termination { reflection: -1.0 });
    let mut pressure_difference: f64 = 0.0;
    let mut motion_difference: f64 = 0.0;
    for n in 128..256 {
        let a = fixed.step(drive).unwrap();
        let b = switched.step(drive).unwrap();
        // Branch return (5) plus inlet section (8): no instant feedback.
        if n < 141 { assert_eq!(a.aperture, b.aperture); }
        pressure_difference = pressure_difference.max((a.aperture.bore_pressure_pa - b.aperture.bore_pressure_pa).abs());
        motion_difference = motion_difference.max((a.aperture.state.opening_m - b.aperture.state.opening_m).abs());
    }
    assert!(pressure_difference > 1.0);
    assert!(motion_difference > 1e-9);
}

#[test]
fn branch_radius_changes_actual_returning_pressure_and_valve_motion() {
    let mut wider = spec(-1.0);
    wider.sections[2].radius_m = 0.006;
    let mut a = model(spec(-1.0), 256, initial());
    let mut b = model(wider, 256, initial());
    let mut difference: f64 = 0.0;
    let mut motion: f64 = 0.0;
    for n in 0..256 {
        let drive = TubeDrive { upstream_pressure_pa: 800.0, body_flow_m3_s: 0.0 };
        let x = a.step(drive).unwrap();
        let y = b.step(drive).unwrap();
        if n < 16 { assert_eq!(x.aperture, y.aperture); }
        difference = difference.max((x.aperture.bore_pressure_pa - y.aperture.bore_pressure_pa).abs());
        motion = motion.max((x.aperture.state.opening_m - y.aperture.state.opening_m).abs());
    }
    assert!(difference > 1.0);
    assert!(motion > 1e-9);
}

#[test]
fn callback_cancel_and_budget_resume_keep_the_entire_network_trajectory() {
    let inputs: Vec<_> = (0..128).map(|i| TubeDrive {
        upstream_pressure_pa: if i < 50 { 800.0 } else { 0.0 },
        body_flow_m3_s: if i % 3 == 0 { -1e-7 } else { 2e-7 },
    }).collect();
    let gate = CancelGate::new_clock_free();
    let mut reference = model(spec(-0.7), 128, initial());
    let mut expected = vec![ApertureNetworkFrame::default(); 128];
    reference.advance_block(&inputs, &mut expected, &gate).unwrap();
    let sentinel = ApertureNetworkFrame { stored_energy_j: -1.0, ..ApertureNetworkFrame::default() };
    let mut actual = vec![sentinel; 128];
    let mut resumed = model(spec(-0.7), 12, initial());
    resumed.advance_block(&inputs[..8], &mut actual[..8], &gate).unwrap();
    let cancel = CancelGate::new_clock_free(); cancel.request();
    let paused = resumed.advance_block(&inputs[8..], &mut actual[8..], &cancel).unwrap();
    assert_eq!(paused.completed, 0); assert_eq!(paused.terminal, ApertureTerminal::Cancelled);
    assert!(actual[8..].iter().all(|f| *f == sentinel));
    let exhausted = resumed.advance_block(&inputs[8..], &mut actual[8..], &gate).unwrap();
    assert_eq!(exhausted.completed, 4); assert_eq!(exhausted.terminal, ApertureTerminal::BudgetExhausted);
    assert!(actual[12..].iter().all(|f| *f == sentinel));
    resumed.extend_step_budget(128).unwrap();
    assert_eq!(resumed.advance_block(&inputs[12..], &mut actual[12..], &gate).unwrap().terminal,
        ApertureTerminal::Complete);
    for (a, b) in actual.into_iter().zip(expected) { assert_eq!(bits(a), bits(b)); }
}

#[test]
fn refused_samples_and_controls_do_not_desynchronize_participants() {
    let mut a = model(spec(-0.7), 256, initial());
    let mut b = model(spec(-0.7), 256, initial());
    let drive = TubeDrive { upstream_pressure_pa: 800.0, body_flow_m3_s: 0.0 };
    for _ in 0..32 { a.step(drive).unwrap(); b.step(drive).unwrap(); }
    let energy = a.stored_energy_j(); let state = a.aperture().state();
    let nodes: Vec<_> = (0..4).map(|n| *a.node_frame(n).unwrap()).collect();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(a.step(TubeDrive { body_flow_m3_s: bad, ..drive }).is_err());
        assert!(a.set_terminal_reflection(3, bad).is_err());
    }
    assert!(a.set_terminal_reflection(1, 0.5).is_err());
    assert!(a.set_terminal_reflection(99, 0.5).is_err());
    assert_eq!(a.aperture().accepted_steps(), 32);
    assert_eq!(a.stored_energy_j().to_bits(), energy.to_bits()); assert_eq!(a.aperture().state(), state);
    for (n, frame) in nodes.iter().enumerate() { assert_eq!(a.node_frame(n).unwrap(), frame); }
    for _ in 0..128 { assert_eq!(bits(a.step(drive).unwrap()), bits(b.step(drive).unwrap())); }
}

#[test]
fn declared_geometry_clock_load_and_payload_limits_are_not_silently_changed() {
    let s = spec(0.5); let z = s.inlet_impedance(1.2).unwrap();
    let m = model(s.clone(), 8, initial());
    for (n, realized) in m.represented_sections().iter().enumerate() {
        assert!((realized.represented_length_m - s.sections[n].length_m).abs() <= 1e-15);
    }
    assert_eq!(m.represented_sections()[0].one_way_samples, 8);
    assert_eq!(m.represented_sections()[1].one_way_samples, 12);
    assert_eq!(m.represented_sections()[2].one_way_samples, 5);
    let mut bad_length = s.clone(); bad_length.sections[1].length_m += 0.3 * 343.0 * DT;
    bad_length.sections[1].max_length_error_m = 0.0;
    assert!(ApertureNetwork::new(aperture(z, 8, initial()), bad_length).is_err());
    let mut low_memory = s.clone(); low_memory.max_wave_memory_bytes = 1;
    assert!(ApertureNetwork::new(aperture(z, 8, initial()), low_memory).is_err());
    assert!(ApertureNetwork::new(aperture(z * 2.0, 8, initial()), s.clone()).is_err());
    let mut running = aperture(z, 8, initial());
    running.step(Default::default()).unwrap();
    assert!(ApertureNetwork::new(running, s).is_err());
}

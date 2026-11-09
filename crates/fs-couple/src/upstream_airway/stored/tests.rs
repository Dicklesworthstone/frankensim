use super::*;
use fs_vfit::impedance::{ImpedanceState, SeriesImpedanceSpec};
use fs_vfit::waveguide::network::{NetworkNode, NetworkSegment};

fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2.0e-12 * scale.max(a.abs()).max(b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}
fn source(c: f64, p: f64, z: f64, dt: f64) -> StoredPressureSource {
    StoredPressureSource::new(
        SeriesImpedance::new_with_state(
            SeriesImpedanceSpec {
                resistance_pa_s_m3: 0.0,
                inertance_pa_s2_m3: 0.0,
                compliance_m3_pa: Some(c),
            },
            z,
            dt,
            ImpedanceState {
                inertive_flow_m3_s: 0.0,
                compliance_pressure_pa: p,
            },
        )
        .unwrap(),
    )
}
fn line(z: f64, dt: f64, delay: usize) -> WaveguideNetwork {
    WaveguideNetwork::new(
        &[NetworkNode::Inlet, NetworkNode::Driven],
        &[NetworkSegment {
            nodes: [0, 1],
            one_way_samples: delay,
            impedance_pa_s_m3: z,
        }],
        dt,
        1_000_000,
    )
    .unwrap()
}

#[test]
fn g0_initial_storage_and_same_sample_flow_close_original_rlc_work() {
    let spec = SeriesImpedanceSpec {
        resistance_pa_s_m3: 3.0,
        inertance_pa_s2_m3: 0.02,
        compliance_m3_pa: Some(0.04),
    };
    let dt = 0.003;
    let mut initial = ImpedanceState {
        inertive_flow_m3_s: -0.03,
        compliance_pressure_pa: 5.0,
    };
    let load = SeriesImpedance::new_with_state(spec, 9.0, dt, initial).unwrap();
    let mut boundary = StoredPressureSource::new(load);
    close(
        boundary.impedance().stored_energy_j(),
        0.5 * (0.02 * 0.03 * 0.03 + 0.04 * 25.0),
        1.0,
    );
    let mut storage = boundary.impedance().stored_energy_j();
    let mut supplied = 0.0;
    let mut delivered = 0.0;
    let mut heat = 0.0;
    for n in 0..512 {
        let source_pa = (n % 7) as f64 - 3.0;
        let outgoing = 0.01 * ((n % 11) as f64 - 5.0);
        let port = boundary.port(source_pa).unwrap();
        let expected_pressure = source_pa + initial.compliance_pressure_pa
            - (2.0 * spec.inertance_pa_s2_m3 / dt) * initial.inertive_flow_m3_s
            - (spec.resistance_pa_s_m3
                + 2.0 * spec.inertance_pa_s2_m3 / dt
                + dt / (2.0 * spec.compliance_m3_pa.unwrap()))
                * outgoing;
        close(port.pressure_pa(outgoing), expected_pressure, 1.0);
        let before = boundary.impedance().state();
        let preview = boundary.preview_commit(outgoing).unwrap();
        assert_eq!(boundary.impedance().state(), before);
        let frame = boundary.commit(outgoing).unwrap();
        assert_eq!(frame, preview);
        close(frame.pressure_pa, expected_pressure, 1.0);
        initial = ImpedanceState {
            inertive_flow_m3_s: -2.0 * outgoing - initial.inertive_flow_m3_s,
            compliance_pressure_pa: initial.compliance_pressure_pa
                - dt * outgoing / spec.compliance_m3_pa.unwrap(),
        };
        close(
            frame.load.state.inertive_flow_m3_s,
            initial.inertive_flow_m3_s,
            1.0,
        );
        close(
            frame.load.state.compliance_pressure_pa,
            initial.compliance_pressure_pa,
            1.0,
        );
        close(
            frame.balance_residual_j(),
            0.0,
            frame.load.stored_energy_j + frame.outgoing_work_j.abs() + frame.source_work_j.abs(),
        );
        close(
            frame.load.dissipated_energy_j,
            spec.resistance_pa_s_m3 * outgoing * outgoing * dt,
            1.0,
        );
        supplied += frame.source_work_j;
        delivered += frame.outgoing_work_j;
        heat += frame.load.dissipated_energy_j;
    }
    storage = boundary.impedance().stored_energy_j() - storage;
    close(storage + heat + delivered, supplied, 1.0);
}

#[test]
fn g1_charged_compliance_reflection_and_time_refinement() {
    // Independent trapezoidal solution of C p' = -p/Z with no incident wave.
    let c = 0.002;
    let z = 2.0;
    let p0 = 7.0;
    let tau = c * z;
    let mut errors = [0.0; 3];
    for (rung, steps) in [8, 16, 32].into_iter().enumerate() {
        let dt = tau / f64::from(steps);
        let mut boundary = source(c, p0, z, dt);
        let mut network = line(z, dt, 1);
        let retention = (2.0 * tau - dt) / (2.0 * tau + dt);
        let mut expected = p0;
        for _ in 0..steps {
            let frame = boundary
                .step_network(&mut network, 1, 0.0, 0.0, || false)
                .unwrap();
            expected *= retention;
            close(frame.source.load.state.compliance_pressure_pa, expected, p0);
            close(frame.source.load.dissipated_energy_j, 0.0, 1.0);
            close(frame.source.source_work_j, 0.0, 1.0);
        }
        errors[rung] = (expected - p0 * fs_math::det::exp(-1.0)).abs();
    }
    assert!(errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9);
}

#[test]
fn g1_native_duct_transit_and_initial_reservoir_energy_are_not_source_work() {
    let (c, p0, z, dt, delay) = (0.002, 7.0, 2.0, 0.001, 3);
    let mut boundary = source(c, p0, z, dt);
    let mut network = line(z, dt, delay);
    let initial_energy = 0.5 * c * p0 * p0;
    let mut emitted = 0.0;
    let mut first = None;
    for n in 0..128 {
        let frame = boundary
            .step_network(&mut network, 1, 0.0, 0.0, || false)
            .unwrap();
        if frame.network.incoming_pressure_pa != 0.0 {
            first.get_or_insert(n);
        }
        let at_source = network.node_frame(1).unwrap();
        close(at_source.pressure_pa, frame.source.pressure_pa, p0);
        close(
            at_source.net_flow_into_node_m3_s,
            -frame.source.flow_out_m3_s,
            p0 / z,
        );
        close(
            frame.network.driven_work_j,
            frame.source.outgoing_work_j,
            initial_energy,
        );
        close(frame.balance_residual_j(), 0.0, initial_energy);
        emitted -= frame.network.inlet_work_j;
        close(
            frame.stored_energy_j() + emitted,
            initial_energy,
            initial_energy,
        );
        assert_eq!(frame.source.source_work_j, 0.0);
    }
    assert_eq!(first, Some(delay));
    assert!(emitted > 0.999 * initial_energy);
}

#[test]
fn g0_closed_duct_returns_energy_to_compliance_without_false_heat() {
    let mut boundary = source(0.002, 7.0, 2.0, 0.001);
    let mut network = line(2.0, 0.001, 4);
    let initial = boundary.impedance().stored_energy_j();
    let mut saw_absorption = false;
    let mut saw_release = false;
    for _ in 0..512 {
        // +1 pressure reflection at the inlet: exactly zero boundary volume flow.
        let incoming = network.incoming_pressure_pa();
        let frame = boundary
            .step_network(&mut network, 1, incoming, 0.0, || false)
            .unwrap();
        saw_absorption |= frame.source.outgoing_work_j < -1.0e-8;
        saw_release |= frame.source.outgoing_work_j > 1.0e-8;
        close(frame.network.inlet_work_j, 0.0, initial);
        close(frame.stored_energy_j(), initial, initial);
        close(frame.balance_residual_j(), 0.0, initial);
    }
    assert!(saw_absorption && saw_release);
}

#[test]
fn g4_failed_flow_network_and_cancellation_preserve_all_accepted_state() {
    let mut boundary = source(0.002, 7.0, 2.0, 0.001);
    let mut network = line(2.0, 0.001, 3);
    assert!(boundary.commit(0.0).is_err());
    let port = boundary.port(4.0).unwrap();
    let state = boundary.impedance().state();
    assert!(boundary.commit(f64::NAN).is_err());
    assert_eq!(boundary.impedance().state(), state);
    assert_eq!(boundary.port(4.0).unwrap(), port);
    let good_preview = boundary.preview_commit(0.1).unwrap();
    assert!(boundary.port(f64::INFINITY).is_err());
    assert_eq!(boundary.preview_commit(0.1).unwrap(), good_preview);
    for cancellation_call in [1, 2] {
        let mut polls = 0;
        assert!(
            boundary
                .step_network(&mut network, 1, 0.0, 0.0, || {
                    polls += 1;
                    polls == cancellation_call
                })
                .is_err()
        );
        assert_eq!(network.accepted_steps(), 0);
        assert_eq!(network.stored_energy_j(), 0.0);
        assert_eq!(boundary.impedance().state(), state);
        assert_eq!(boundary.preview_commit(0.1).unwrap(), good_preview);
    }
    assert!(
        boundary
            .step_network(&mut network, 1, f64::NAN, 0.0, || false)
            .is_err()
    );
    assert_eq!(network.accepted_steps(), 0);
    assert_eq!(boundary.impedance().state(), state);
    let mut wrong_clock = line(2.0, 0.002, 3);
    assert!(
        boundary
            .step_network(&mut wrong_clock, 1, 0.0, 0.0, || false)
            .is_err()
    );
    let mut wrong_reference = line(3.0, 0.001, 3);
    assert!(
        boundary
            .step_network(&mut wrong_reference, 1, 0.0, 0.0, || false)
            .is_err()
    );
    assert!(
        boundary
            .step_network(&mut network, 0, 0.0, 0.0, || false)
            .is_err()
    );
    let accepted = boundary
        .step_network(&mut network, 1, 0.0, 0.0, || false)
        .unwrap();
    assert_eq!(network.accepted_steps(), 1);
    close(accepted.balance_residual_j(), 0.0, 1.0);
    assert!(boundary.commit(0.0).is_err());
}

#[test]
fn g0_initialized_load_admission_and_native_wave_flow_parity() {
    let spec = SeriesImpedanceSpec {
        resistance_pa_s_m3: 4.0,
        inertance_pa_s2_m3: 0.05,
        compliance_m3_pa: Some(0.02),
    };
    let state = ImpedanceState {
        inertive_flow_m3_s: 0.01,
        compliance_pressure_pa: 3.0,
    };
    let mut wave = SeriesImpedance::new_with_state(spec, 8.0, 0.004, state).unwrap();
    let mut flow = wave;
    for i in 0..128 {
        let frame = wave.step((i % 9) as f64 - 4.0).unwrap();
        let driven = flow.step_at_flow(frame.flow_m3_s).unwrap();
        assert_eq!(frame.state, driven.state);
        close(frame.pressure_pa, driven.pressure_pa, 10.0);
        close(frame.stored_energy_j, driven.stored_energy_j, 1.0);
        close(driven.balance_residual_j(), 0.0, 1.0);
    }
    let absent = SeriesImpedanceSpec {
        resistance_pa_s_m3: 0.0,
        inertance_pa_s2_m3: 0.0,
        compliance_m3_pa: None,
    };
    for bad in [
        state,
        ImpedanceState {
            inertive_flow_m3_s: f64::NAN,
            compliance_pressure_pa: 0.0,
        },
        ImpedanceState {
            inertive_flow_m3_s: 0.0,
            compliance_pressure_pa: 1.0,
        },
    ] {
        assert!(SeriesImpedance::new_with_state(absent, 1.0, 0.001, bad).is_err());
    }
    assert!(
        SeriesImpedance::new_with_state(
            spec,
            1.0,
            0.001,
            ImpedanceState {
                inertive_flow_m3_s: 0.0,
                compliance_pressure_pa: f64::MAX
            }
        )
        .is_err()
    );
    let zero = SeriesImpedance::new(absent, 1.0, 0.001).unwrap();
    let initialized =
        SeriesImpedance::new_with_state(absent, 1.0, 0.001, ImpedanceState::default()).unwrap();
    assert_eq!(
        zero.preview_step(3.0).unwrap(),
        initialized.preview_step(3.0).unwrap()
    );
}

#[test]
fn g0_driven_resistive_source_closes_joint_storage_and_heat() {
    let mut source = StoredPressureSource::new(
        SeriesImpedance::new_with_state(
            SeriesImpedanceSpec {
                resistance_pa_s_m3: 0.4,
                inertance_pa_s2_m3: 0.01,
                compliance_m3_pa: Some(0.002),
            },
            2.0,
            0.001,
            ImpedanceState {
                inertive_flow_m3_s: 0.1,
                compliance_pressure_pa: 1.0,
            },
        )
        .unwrap(),
    );
    let mut network = line(2.0, 0.001, 3);
    let initial = source.impedance().stored_energy_j();
    let mut external = 0.0;
    let mut heat = 0.0;
    for n in 0..1024 {
        let drive = if n < 256 { (n % 9) as f64 - 4.0 } else { 0.0 };
        let frame = source
            .step_network(&mut network, 1, 0.0, drive, || false)
            .unwrap();
        close(
            frame.network.driven_work_j,
            frame.source.outgoing_work_j,
            1.0,
        );
        close(frame.source.balance_residual_j(), 0.0, 1.0);
        close(frame.balance_residual_j(), 0.0, 1.0);
        external += frame.network.inlet_work_j + frame.source.source_work_j;
        heat += frame.source.load.dissipated_energy_j;
        close(frame.stored_energy_j() - initial + heat, external, 1.0);
    }
    assert!(heat > 0.0);
}

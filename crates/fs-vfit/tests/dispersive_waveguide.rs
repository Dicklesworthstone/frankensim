//! G0/G1/G3: physical phase dispersion, actual storage and native network use.
use fs_math::c64::C64;
use fs_vfit::{
    impedance::SeriesImpedanceSpec,
    waveguide::{
        PassiveWaveguide, WaveguideSpec,
        dispersion::AllpassSection,
        fractional::FractionalTransit,
        network::{NetworkNode, NetworkSegment, WaveguideNetwork},
    },
};

const DT: f64 = 1.0 / 48_000.0;
fn spec(reflection: f64) -> WaveguideSpec {
    WaveguideSpec { one_way_samples: 7, impedance_pa_s_m3: 1e6,
        time_step_s: DT, reflection, max_memory_bytes: 1 << 20 }
}
fn sections<const N: usize>(a: [f64; N]) -> [AllpassSection; N] {
    a.map(|a| AllpassSection::new(a).unwrap())
}
fn fraction(d: f64) -> FractionalTransit { FractionalTransit::new(d).unwrap() }

#[test]
fn g3_empty_cascade_preserves_integer_and_fractional_frames() {
    for d in [0.0, 0.23] {
        let mut original = PassiveWaveguide::with_fractional_transit(spec(-0.7), fraction(d)).unwrap();
        let mut dispersive = PassiveWaveguide::with_dispersive_transit(spec(-0.7), fraction(d), &[]).unwrap();
        for n in 0..1000 {
            let x = f64::from(n % 19) - 9.0;
            assert_eq!(original.step(x).unwrap(), dispersive.step(x).unwrap());
        }
    }
}

#[test]
fn g1_cascade_matches_independent_direct_form_in_both_directions() {
    let coefficients = [-0.8, 0.0, 0.3, -0.5];
    let cascade = sections(coefficients);
    for reflection in [-1.0, 0.4, 1.0] {
        let mut line = PassiveWaveguide::with_dispersive_transit(spec(reflection), fraction(0.0), &cascade).unwrap();
        let mut buffers = [[0.0; 7]; 2];
        let mut x1 = [[0.0; 4]; 2];
        let mut y1 = x1;
        let mut work = 0.0;
        let mut loss = 0.0;
        for n in 0..4096 {
            let input = if n < 192 { (n as f64 * 0.37).sin() } else { 0.0 };
            let head = n % 7;
            let incident = buffers[0][head];
            let incoming = buffers[1][head];
            let before = line.stored_energy_j();
            let frame = line.step(input).unwrap();
            assert!((frame.incoming_pressure_pa - incoming).abs() < 2e-13);
            assert!((frame.terminal_pressure_pa - (1.0 + reflection) * incident).abs() < 2e-13);
            for (direction, mut wave) in [input, reflection * incident].into_iter().enumerate() {
                for (k, a) in coefficients.into_iter().enumerate() {
                    let next = a * wave + x1[direction][k] - a * y1[direction][k];
                    x1[direction][k] = wave;
                    y1[direction][k] = next;
                    wave = next;
                }
                buffers[direction][head] = wave;
            }
            // Independent buffer sum plus the publicly observed actual memories.
            let direct = (buffers.iter().flatten().map(|v| v * v).sum::<f64>()
                + line.dispersion_stages().iter().flat_map(|s| s.state_pa()).map(|v| v * v).sum::<f64>())
                * DT / 1e6;
            let scale = before + frame.stored_energy_j + frame.inlet_work_j.abs() + frame.terminal_loss_j;
            assert!((direct - frame.stored_energy_j).abs() < 2e-12 * scale.max(1e-300));
            assert!(frame.balance_residual_j().abs() < 2e-12 * scale.max(1e-300));
            assert_eq!(line.stored_energy_j(), frame.stored_energy_j);
            work += frame.inlet_work_j;
            loss += frame.terminal_loss_j;
        }
        assert!((work - loss - line.stored_energy_j()).abs() < 4e-23);
    }
}

#[test]
fn g1_measured_wave_has_unit_gain_and_original_frequency_dependent_phase() {
    let cascade = sections([-0.7, -0.35, 0.2]);
    let d = 0.29;
    for bin in [1, 11, 43, 113, 251] {
        let omega = core::f64::consts::TAU * f64::from(bin) / 512.0;
        let mut line = PassiveWaveguide::with_dispersive_transit(spec(0.0), fraction(d), &cascade).unwrap();
        // Independent z-transform; production diagnostics use real half-angle laws.
        let z = C64::new(omega.cos(), -omega.sin());
        let mut expected = C64::new((7.0 * omega).cos(), -(7.0 * omega).sin());
        for a in [(1.0 - d) / (1.0 + d), -0.7, -0.35, 0.2] {
            expected = expected * (C64::from_re(a) + z) / (C64::ONE + z.scale(a));
        }
        let mut measured = C64::ZERO;
        for n in 0..2560 {
            let angle = f64::from(n) * omega;
            let y = line.step(angle.cos()).unwrap().terminal_pressure_pa;
            if n >= 2048 { measured = measured + C64::new(y * angle.cos(), -y * angle.sin()); }
        }
        measured = measured.scale(1.0 / 256.0);
        assert!((measured - expected).abs() < 2e-11);
        assert!((measured.abs() - 1.0).abs() < 2e-11);
        let response = line.propagation_response(omega).unwrap();
        assert!((C64::new(response.phase_lag_rad.cos(), -response.phase_lag_rad.sin()) - expected).abs() < 2e-13);
        let h = 1e-6;
        let slope = (line.propagation_response(omega + h).unwrap().phase_lag_rad
            - line.propagation_response(omega - h).unwrap().phase_lag_rad) / (2.0 * h);
        assert!((slope - response.group_delay_samples).abs() < 3e-8);
    }
    // Authored stiff-string-direction dispersion delays low frequencies more.
    let stiff = AllpassSection::new(-0.75).unwrap();
    assert!(stiff.group_delay_samples(0.1).unwrap() > stiff.group_delay_samples(1.0).unwrap());
    assert!(stiff.phase_lag_rad(0.8).unwrap() < 8.0 * stiff.phase_lag_rad(0.1).unwrap());
}

#[test]
fn g0_limits_keep_exact_sample_delay_and_resolve_near_unit_circle_group_delay() {
    let zero = AllpassSection::new(0.0).unwrap();
    for omega in [0.0, 0.001, 1.0, core::f64::consts::PI.next_down()] {
        assert!((zero.group_delay_samples(omega).unwrap() - 1.0).abs() < 5e-16);
        assert!((zero.phase_lag_rad(omega).unwrap() - omega).abs() < 8e-16);
    }
    for a in [(-1.0_f64).next_up(), 1.0_f64.next_down()] {
        let section = AllpassSection::new(a).unwrap();
        let expected = (1.0 - a) / (1.0 + a);
        assert!((section.group_delay_samples(0.0).unwrap() / expected - 1.0).abs() < 5e-16);
        assert!(section.group_delay_samples(core::f64::consts::PI.next_down()).unwrap().is_finite());
    }
    for a in [-1.0, 1.0, f64::NAN, f64::INFINITY] { assert!(AllpassSection::new(a).is_err()); }
    for omega in [-0.1, core::f64::consts::PI, f64::NAN, f64::INFINITY] {
        assert!(zero.phase_lag_rad(omega).is_err());
        assert!(zero.group_delay_samples(omega).is_err());
    }
    let cascade = sections([0.0; 3]);
    let mut line = PassiveWaveguide::with_dispersive_transit(spec(0.0), fraction(0.0), &cascade).unwrap();
    for n in 0..20 {
        let y = line.step(if n == 0 { 1.0 } else { 0.0 }).unwrap().terminal_pressure_pa;
        assert_eq!(y, if n == 10 { 1.0 } else { 0.0 });
    }
}

#[test]
fn g0_budget_and_refused_step_preserve_every_retained_state() {
    let cascade = sections([-0.9, 0.0, 0.8]);
    let mut s = spec(-0.7);
    s.max_memory_bytes = PassiveWaveguide::required_memory_bytes(7, cascade.len()).unwrap();
    let mut line = PassiveWaveguide::with_dispersive_transit(s, fraction(0.3), &cascade).unwrap();
    s.max_memory_bytes -= 1;
    assert!(PassiveWaveguide::with_dispersive_transit(s, fraction(0.3), &cascade).is_err());
    assert!(PassiveWaveguide::required_memory_bytes(usize::MAX, 1).is_err());
    assert!(PassiveWaveguide::required_memory_bytes(7, usize::MAX).is_err());
    for n in 0..30 { line.step(f64::from(n % 7) - 3.0).unwrap(); }
    let mut reference = line.clone();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(line.step(bad).is_err());
        assert_eq!(line.dispersion_stages(), reference.dispersion_stages());
        assert_eq!(line.fractional_state_pa(), reference.fractional_state_pa());
        assert_eq!(line.stored_energy_j(), reference.stored_energy_j());
    }
    for n in 0..100 { assert_eq!(line.step(f64::from(n % 3)).unwrap(), reference.step(f64::from(n % 3)).unwrap()); }
}

fn branched() -> WaveguideNetwork {
    let nodes = [NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.7 }, NetworkNode::Impedance { load: SeriesImpedanceSpec {
            resistance_pa_s_m3: 2e5, inertance_pa_s2_m3: 20.0, compliance_m3_pa: Some(1e-9) } }];
    let segments = [NetworkSegment { nodes: [0, 1], one_way_samples: 3, impedance_pa_s_m3: 1e6 },
        NetworkSegment { nodes: [1, 2], one_way_samples: 5, impedance_pa_s_m3: 2e6 },
        NetworkSegment { nodes: [1, 3], one_way_samples: 7, impedance_pa_s_m3: 3e6 }];
    let a = sections([-0.7, 0.3]); let b = sections([-0.4]);
    WaveguideNetwork::with_dispersive_transits(&nodes, &segments, &[fraction(0.2), fraction(0.4), fraction(0.0)],
        &[&a, &b, &[]], DT, 1 << 20).unwrap()
}

#[test]
fn g1_network_uses_original_line_and_accounts_for_both_directions() {
    let cascade = sections([-0.7, -0.2]);
    let nodes = [NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.8 }];
    let segments = [NetworkSegment { nodes: [0, 1], one_way_samples: 7, impedance_pa_s_m3: 1e6 }];
    let needed = WaveguideNetwork::required_dispersive_memory_bytes(2, &segments, &[&cascade]).unwrap();
    let mut net = WaveguideNetwork::with_dispersive_transits(&nodes, &segments, &[fraction(0.37)], &[&cascade], DT, needed).unwrap();
    assert!(WaveguideNetwork::with_dispersive_transits(&nodes, &segments, &[fraction(0.37)], &[&cascade], DT, needed - 1).is_err());
    let mut line = PassiveWaveguide::with_dispersive_transit(spec(-0.8), fraction(0.37), &cascade).unwrap();
    for n in 0..1024 {
        let x = if n < 100 { (n as f64 * 0.3).sin() } else { 0.0 };
        let a = line.step(x).unwrap(); let b = net.step(x).unwrap();
        assert_eq!(a.incoming_pressure_pa, b.incoming_pressure_pa);
        assert_eq!(a.terminal_pressure_pa, net.node_frame(1).unwrap().pressure_pa);
        assert_eq!(a.stored_energy_j, b.stored_energy_j);
        assert_eq!(line.dispersion_stages(), net.section_dispersion_stages(0).unwrap());
    }
    assert_eq!(line.propagation_response(0.7).unwrap(), net.section_propagation_response(0, 0.7).unwrap());
    assert!(net.section_dispersion_stages(1).is_none());
    assert!(net.section_propagation_response(1, 0.7).is_err());
    let (mut a, mut b, mut ra, mut rb) = (branched(), branched(), branched(), branched());
    for n in 0..1024 {
        let before = a.stored_energy_j();
        let x = [if n < 128 { f64::from(n % 7) - 3.0 } else { 0.0 }, 0.0];
        let f = a.step_pair(&mut b, x).unwrap();
        assert_eq!(f, ra.step_pair(&mut rb, x).unwrap());
        let scale = before + f[0].stored_energy_j + f[0].inlet_work_j.abs() + f[0].terminal_loss_j;
        assert!(f[0].balance_residual_j().abs() < 3e-12 * scale.max(1e-300));
        assert_eq!(a.stored_energy_j(), f[0].stored_energy_j);
    }
    assert!(a.step_pair(&mut b, [1.0, f64::MAX]).is_err());
    for i in 0..3 {
        assert_eq!(a.section_dispersion_stages(i), ra.section_dispersion_stages(i));
        assert_eq!(b.section_dispersion_stages(i), rb.section_dispersion_stages(i));
    }
    for _ in 0..100 { assert_eq!(a.step_pair(&mut b, [0.0; 2]).unwrap(), ra.step_pair(&mut rb, [0.0; 2]).unwrap()); }
}

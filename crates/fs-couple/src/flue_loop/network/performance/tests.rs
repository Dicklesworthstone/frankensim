use super::*;
use crate::flue_loop::network::tests::{gas, source, voice};
use crate::bernoulli_aperture::network::TubeSection;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::gesture::{GestureEvent, GestureTarget, GestureTrack, GestureValue};

fn schedule(pressure: f64) -> GestureSchedule {
    GestureSchedule::try_new(700, vec![GestureTrack { id: "blow".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(pressure),
        events: vec![] }]).unwrap()
}
fn config(node: usize) -> FlueOutletConfig {
    FlueOutletConfig {
        sample_rate_hz: 48_000, samples: 513, max_block: 513, max_controls: 64,
        max_compile_work: 10_000, outlet_node: node,
        receiver: CircularOutletReceiver { position_m: [0.0, 0.0, 1.0], radial_rings: 4,
            angular_points: 16, maximum_frequency_hz: 1000.0 },
        maximum_abs_pressure_pa: 1.0e6,
    }
}

#[test]
fn exterior_output_matches_actual_terminal_flow_through_the_shared_receiver() {
    let g = gas();
    let spec = source(g, -0.9);
    let mut original = voice(spec.clone());
    original.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let cfg = config(1);
    let mut receiver = BaffledPressure::circular_outlet(0.006, cfg.sample_rate_hz, cfg.receiver,
        RayleighMedium { density: g.density, sound_speed: g.sound_speed }).unwrap();
    let tube_delay = original.network().unwrap().represented_sections()[0].one_way_samples;
    let first = tube_delay + receiver.delay_samples.0;
    let mut launched = vec![0.0; 513];
    let mut expected = Vec::new();
    original.step_block_observed(&mut launched, |line| {
        let FlueResonator::Network(n) = line else { panic!("network") };
        expected.push(receiver.step(&[n.node_frame(1).unwrap().net_flow_into_node_m3_s]).unwrap());
        Ok(())
    }).unwrap();
    let mut performance = FlueOutletPerformance::new(voice(spec), schedule(300.0), cfg).unwrap();
    let mut actual = vec![0.0; 513];
    performance.block(&mut actual).unwrap();
    assert_eq!(actual, expected);
    assert!(actual[..first].iter().all(|p| *p == 0.0));
    assert!(actual[first..].iter().any(|p| p.abs() > 1.0e-12));
    assert_ne!(actual, launched); // internal waves must not masquerade as microphone pressure
    assert_eq!(performance.voice().diagnostics().len(), 1);
}

#[test]
fn nondivisor_pressure_gestures_and_release_preserve_physics_across_callbacks() {
    let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
        id: "blow".into(), target: GestureTarget::BlowingPressure,
        initial: GestureValue::PressurePa(300.0), events: vec![
            GestureEvent { time_s: 1.0/700.0, transition_s: 4.0/700.0, value: GestureValue::PressurePa(450.0) },
            GestureEvent { time_s: 3.0/700.0, transition_s: 2.0/700.0, value: GestureValue::PressurePa(200.0) },
            GestureEvent { time_s: 6.0/700.0, transition_s: 0.0, value: GestureValue::PressurePa(0.0) },
        ],
    }]).unwrap();
    let mut expected = Vec::new();
    for partition in [513, 1, 37, 256] {
        let mut performance = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule.clone(), config(1)).unwrap();
        let mut out = vec![0.0; 513];
        for block in out.chunks_mut(partition) { performance.block(block).unwrap(); }
        let bits: Vec<_> = out.iter().map(|p| p.to_bits()).collect();
        if expected.is_empty() { expected = bits; } else { assert_eq!(bits, expected); }
        assert_eq!(performance.samples_rendered(), 513);
        assert_eq!(performance.remaining_samples(), 0);
        assert_eq!(performance.voice().jet_speed_m_s(), 0.0);
        assert!(performance.voice().network().unwrap().stored_energy_j() > 0.0);
        assert!(out[480..].iter().any(|p| p.abs() > 1.0e-12));
    }
}

#[test]
fn microphone_uses_the_outlets_local_gas_not_the_jets_inlet_gas() {
    let inlet = gas();
    let outlet = GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, inlet.pressure).unwrap();
    let mut spec = source(inlet, -0.9);
    spec.graph.nodes = vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.9 }];
    spec.graph.sections = vec![
        TubeSection { nodes: [0,1], length_m: 0.15, radius_m: 0.006, max_length_error_m: 0.005 },
        TubeSection { nodes: [1,2], length_m: 0.15, radius_m: 0.006, max_length_error_m: 0.005 },
    ];
    spec.gases = vec![inlet, outlet];
    let cfg = config(2);
    let mut direct = voice(spec.clone());
    direct.apply(FlueControl::SetBlowingPressure(300.0)).unwrap();
    let mut correct = BaffledPressure::circular_outlet(0.006, 48_000, cfg.receiver,
        RayleighMedium { density: outlet.density, sound_speed: outlet.sound_speed }).unwrap();
    let mut wrong = BaffledPressure::circular_outlet(0.006, 48_000, cfg.receiver,
        RayleighMedium { density: inlet.density, sound_speed: inlet.sound_speed }).unwrap();
    let mut expected = Vec::new();
    let mut incorrect = Vec::new();
    direct.step_block_observed(&mut [0.0; 513], |line| {
        let FlueResonator::Network(n) = line else { panic!("network") };
        let q = n.node_frame(2).unwrap().net_flow_into_node_m3_s;
        expected.push(correct.step(&[q]).unwrap());
        incorrect.push(wrong.step(&[q]).unwrap());
        Ok(())
    }).unwrap();
    let mut performance = FlueOutletPerformance::new(voice(spec), schedule(300.0), cfg).unwrap();
    let mut actual = vec![0.0; 513];
    performance.block(&mut actual).unwrap();
    assert_eq!(actual, expected);
    assert!(actual.iter().zip(&incorrect).any(|(a,b)| (a-b).abs() > 1.0e-9));
}

#[test]
fn zero_pressure_silence_and_failures_do_not_become_fake_audio() {
    let mut quiet = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(0.0), config(1)).unwrap();
    let mut out = [7.0; 513];
    quiet.block(&mut out).unwrap();
    assert!(out.iter().all(|p| *p == 0.0));
    let mut cfg = config(1);
    cfg.max_block = 64;
    let mut model = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(300.0), cfg).unwrap();
    let mut oversized = [7.0; 65];
    assert!(model.block(&mut oversized).is_err());
    assert_eq!(oversized, [7.0; 65]);
    assert_eq!(model.samples_rendered(), 0);
    assert!(model.validate_sample_rate(44_100).is_err());
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(1.0e9), config(1)).is_err());
    assert!(FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(300.0), config(0)).is_err());
    cfg = config(1);
    cfg.maximum_abs_pressure_pa = f64::MIN_POSITIVE;
    let mut limited = FlueOutletPerformance::new(voice(source(gas(), -0.9)), schedule(300.0), cfg).unwrap();
    assert!(limited.block(&mut out).is_err());
    assert!(limited.samples_rendered() > 0);
    assert!(matches!(limited.block(&mut [0.0]), Err(RenderError::Poisoned)));
}

mod embouchure;

fn material_observer_input(selection: &str) -> String {
    let source = include_str!("../../../../examples/flue-network.performance")
        .split_once("\ngesture\n")
        .unwrap()
        .0
        .replace("output 24000", "output 1024")
        .replace("exterior same", "exterior 310 101325 0.4")
        .replace(
            "receiver 0 0 0.5 4 16 1000",
            &format!("receiver 0.1 0 0.5 2 8 1000{selection}"),
        )
        .replace(
            "controls 4096",
            "mouth_aperture rectangle 0.012 0.008\nmouth_radiation -0.3 0 4 2\ncontrols 4096",
        );
    format!(
        "{source}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t3e2\nevents\t0\n"
    )
}

#[test]
fn material_exterior_reaches_every_flue_opening_without_changing_source_motion() {
    use file::FlueFilePerformance;
    use fs_exec::CancelGate;
    let gate = CancelGate::new_clock_free();
    let original = material_observer_input("");
    let selected = material_observer_input("\nobserver_gas stokes-kirchhoff 4096 1000000 50000000");
    let mut bare = FlueFilePerformance::from_bytes(original.as_bytes(), 1024, &gate)
        .unwrap()
        .into_renderer();
    let mut physical = FlueFilePerformance::from_bytes(selected.as_bytes(), 1024, &gate)
        .unwrap()
        .into_renderer();
    let mut partitioned = FlueFilePerformance::from_bytes(selected.as_bytes(), 37, &gate)
        .unwrap()
        .into_renderer();
    assert!(bare.observer_gas().is_none());
    assert!(bare.mouth_observer_gas().is_none());
    let terminal = physical.observer_gas().unwrap();
    let mouth = physical.mouth_observer_gas().unwrap();
    assert_eq!(terminal.gas, physical.radiation_input().unwrap().exterior);
    assert_ne!(terminal.gas, *physical.voice().gas());
    assert_eq!(terminal.gas, mouth.gas);
    assert_eq!(terminal.surface_filter_indices.len(), 2 * 2 * 8);
    assert_eq!(mouth.surface_filter_indices.len(), 4 * 2);
    assert_eq!(
        terminal.additional_latency_samples,
        mouth.additional_latency_samples
    );
    assert_eq!(
        physical.receiver().filter_delay_samples(),
        physical.mouth_receiver().unwrap().filter_delay_samples()
    );
    let paths: Vec<_> = terminal
        .filters
        .iter()
        .chain(&mouth.filters)
        .map(|r| r.range_m)
        .collect();
    assert!(
        paths.iter().copied().fold(0.0, f64::max)
            - paths.iter().copied().fold(f64::INFINITY, f64::min)
            > 0.1
    );
    let mut expected = vec![0.0; 1024];
    let mut plain = vec![0.0; 1024];
    let mut blocked = vec![0.0; 1024];
    for n in 0..1024 {
        bare.block(&mut plain[n..n + 1]).unwrap();
        physical.block(&mut expected[n..n + 1]).unwrap();
        assert_eq!(
            bare.voice().mouth_flow_m3_s(),
            physical.voice().mouth_flow_m3_s()
        );
        assert_eq!(
            bare.voice().network().unwrap().stored_energy_j(),
            physical.voice().network().unwrap().stored_energy_j()
        );
        for node in [1, 3] {
            assert_eq!(
                bare.voice()
                    .network()
                    .unwrap()
                    .node_frame(node)
                    .unwrap()
                    .net_flow_into_node_m3_s,
                physical
                    .voice()
                    .network()
                    .unwrap()
                    .node_frame(node)
                    .unwrap()
                    .net_flow_into_node_m3_s
            );
        }
    }
    for chunk in blocked.chunks_mut(37) {
        partitioned.block(chunk).unwrap();
    }
    assert_eq!(expected, blocked);
    assert!(
        expected
            .iter()
            .zip(&plain)
            .any(|(a, b)| (a - b).abs() > 1e-14)
    );
    assert_eq!(physical.samples_rendered(), 1024);
    assert_eq!(physical.voice().samples_rendered(), 2048);
}

#[test]
fn molecular_flue_observation_uses_the_explicit_exterior_species_and_refuses_missing_sources() {
    use crate::air_path::GasAbsorptionModel;
    use file::FlueFilePerformance;
    use fs_exec::CancelGate;
    let gate = CancelGate::new_clock_free();
    let input = material_observer_input(
        "\nobserver_gas molecular 4096 1000000 50000000 iso9613-terrestrial-air",
    );
    // The inlet has a CO2-rich mixture. The exterior supplies its own species,
    // so neither inlet composition nor a legacy air label can stand in for it.
    assert!(FlueFilePerformance::from_bytes(input.as_bytes(), 37, &gate).is_err());
    let input = input.replace(
        "exterior 310 101325 0.4",
        concat!(
            "exterior mixture 293.15 101325 mole 4 ",
            "N2 0.78084 O2 0.20946 Ar 0.00934 CO2 0.00036"
        ),
    );
    let parsed = FlueFilePerformance::from_bytes(input.as_bytes(), 37, &gate).unwrap();
    assert_eq!(
        parsed.renderer().observer_gas().unwrap().law,
        GasAbsorptionModel::MolecularRelaxation
    );
    assert_eq!(
        parsed.renderer().observer_gas().unwrap().gas,
        parsed.renderer().mouth_observer_gas().unwrap().gas
    );
    let nonterrestrial = input.replace(
        "N2 0.78084 O2 0.20946 Ar 0.00934 CO2 0.00036",
        "N2 0.5 O2 0.2 Ar 0.1 CO2 0.2",
    );
    assert!(FlueFilePerformance::from_bytes(nonterrestrial.as_bytes(), 37, &gate).is_err());
    for selected in [
        "\nobserver_gas stokes-kirchhoff 0 1000000 50000000",
        "\nobserver_gas lossless 4096 1000000 50000000 ignored",
    ] {
        assert!(
            FlueFilePerformance::from_bytes(
                material_observer_input(selected).as_bytes(),
                37,
                &gate
            )
            .is_err()
        );
    }
    let law = GasPropagation::classical(
        crate::pcm_wav::baffled::gas::GasPropagationLimits::surface_reference(),
    );
    let p = parsed.into_renderer();
    assert!(
        p.with_mouth_radiation(FlueMouthRadiation {
            center_m: [-0.3, 0.0],
            aperture_size_m: None,
            cells: [4, 2],
        })
        .is_err()
    );
    let cold = || {
        FlueFilePerformance::from_bytes(
            material_observer_input("").as_bytes(),
            37,
            &CancelGate::new_clock_free(),
        )
        .unwrap()
        .into_renderer()
    };
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(cold().with_observer_gas(&law, &cancelled).is_err());
    let mut warm = cold();
    warm.block(&mut [0.0]).unwrap();
    assert!(warm.with_observer_gas(&law, &gate).is_err());
    assert!(
        cold()
            .with_observer_gas(&law, &gate)
            .unwrap()
            .with_observer_gas(&law, &gate)
            .is_err()
    );
}

#[test]
fn measured_exit_negative_current_feedback_keeps_the_jet_mouth_residual_increasing() {
    // The committed D4 flute card with the opt-in measured jet exit (coarse kernel row). Its forced
    // kernel opens with a negative precursor at the labium, so g_0 < 0 and the same-sample jet
    // feedback is negative. At every source sample of a block, rebuild the residual the jet-mouth
    // solve brackets (this sample's boundary, radiation, previous flows and prepared kernel) and
    // require it strictly increasing over [lo, hi]; the solve itself must then accept the sample.
    use crate::flue_loop::{JetStepCoefficients, RetardedVelocity, solve_mouth_port};
    use crate::flue_loop::mouth_radiation::MouthRadiation;
    use crate::flue_loop::window_radiation::WindowRadiation;
    use file::FlueFilePerformance;
    use fs_exec::CancelGate;
    const CARD: &str =
        include_str!("../../../../../../data/instruments/descriptions/flute/flute-d4.performance");
    const JET: &str = "jet 90 1e-6 3 0\n";
    const SAMPLES: usize = 2048;
    const GRID: usize = 129;
    assert_eq!(CARD.matches(JET).count(), 1, "the D4 jet row changed; update this test");
    // In band: 6 Pa gives U0 = 3.19 m/s and Re_h = 236 on this card (the gate admits 100 to 400).
    const HOLD: &str = "pressure-pa\t3.32e2";
    assert_eq!(CARD.matches(HOLD).count(), 1, "the D4 hold pressure changed; update this test");
    let text = CARD.replace(HOLD, "pressure-pa\t6e0")
        .replace(JET, &format!("{JET}jet_response forced-bickley 0.25 40 exit measured\n"));
    let mut renderer = FlueFilePerformance::from_bytes(text.as_bytes(), 256,
        &CancelGate::new_clock_free()).unwrap().into_renderer();
    let v = &mut renderer.voice;
    // The residual below is the solve's branch for this card: no published transfer, no unsteady
    // windway, no mouth contraction.
    assert!(v.forced.is_some() && v.nolle.is_none() && v.unsteady_windway.is_none()
        && v.mouth_contraction.is_none());
    v.apply(FlueControl::SetBlowingPressure(6.0)).unwrap();
    let (dt, b, area) = (v.dt, 0.5*v.geometry.flue_height_m, v.geometry.mouth_area_m2());
    let (mut lowest, mut highest) = (f64::INFINITY, f64::NEG_INFINITY);
    for n in 0..SAMPLES {
        // This sample's inputs as advance_block forms them for this card.
        let u = v.jet_speed_m_s();
        let card_u = v.card_speed_at(v.blow_pa, u);
        let q_max = u*v.geometry.flue_height_m*v.geometry.flue_width_m;
        let seed_scale = v.island.seed_amplitude_rel*card_u;
        // The stream is Copy: the copy makes the voice's next draw and leaves the voice's own.
        let mut stream = v.seed;
        let seed = if seed_scale > 0.0 { seed_scale*(2.0*stream.next_f64()-1.0) } else { 0.0 };
        let omega_c = card_u/(4.0*v.momentum_thickness_m);
        let c = JetStepCoefficients { q_max, half_width_m: b,
            offset_m: v.geometry.effective_offset_m(), band_alpha: omega_c*dt/(1.0+omega_c*dt),
            hp_alpha: 1.0/(1.0+0.5*omega_c*dt), leak_rate: omega_c/100.0, eta_override: None,
            at_rest: omega_c == 0.0 };
        // prepare resets and recomputes the current pair from the retained intervals alone, so
        // advance_block's own prepare of this sample repeats it exactly.
        assert!(v.forced.as_mut().unwrap().prepare(u).unwrap().is_some());
        let g0 = v.forced.as_ref().unwrap().current().1;
        let delayed = RetardedVelocity { past_m_s: seed, current_weight: 1.0 };
        let boundary = v.line.mouth_boundary().unwrap();
        let (md, mh) = v.mouth_radiation.as_ref().map_or((0.0, 0.0), MouthRadiation::affine_pressure);
        let (wd, wh) = v.window_radiation.as_ref().map_or((0.0, 0.0), WindowRadiation::affine_pressure);
        let (direct_pa, history_pa) = (md+wd, mh+wh);
        let (l, r, zc) = (v.inertance, v.resistance+direct_pa, v.zc);
        let direct = boundary.direct_reflection;
        let compliance = (1.0-direct)/((l/dt+r)*(1.0-direct)+zc*(1.0+direct));
        let feedback = compliance*(v.source_coeff/dt)*(0.5*q_max)*(g0*delayed.current_weight/area)/b;
        assert!(g0 < 0.0 && feedback < 0.0, "sample {n}: g_0 {g0:e}, feedback {feedback:e}");
        lowest = lowest.min(feedback);
        highest = highest.max(feedback);
        let port = |source_pa: f64| solve_mouth_port(boundary, source_pa-history_pa, v.q_prev,
            l, r, zc, dt).unwrap().flow_m3_s;
        let lo = port(v.source_coeff*(0.0-v.flow_prev)/dt).min(0.0);
        let hi = port(v.source_coeff*(q_max-v.flow_prev)/dt).max(0.0);
        let residual = |q: f64| q-port(v.jet_trial(delayed.past_m_s+delayed.current_weight*q/area,
            0.0, c).unwrap().source_pa);
        let mut last = residual(lo);
        for k in 1..GRID {
            let q = lo+(hi-lo)*(k as f64)/((GRID-1) as f64);
            let next = residual(q);
            assert!(next > last, "sample {n}: residual falls at q {q:e}: {last:e} -> {next:e}");
            last = next;
        }
        v.step_block(&mut [0.0]).unwrap();
    }
    println!("D4 measured exit, 6 Pa: {SAMPLES} samples, current feedback {lowest:e} to {highest:e}");
}

// The example card with a short output, a held 300 Pa gesture (its own events fall outside a short
// output, as in `material_observer_input`) and the given rows after `jet`. The Auvray row states its
// own exit speed U = sqrt(2P/rho), so `bare` drops the example's slit windway.
fn auvray_input(rows: &str, bare: bool) -> String {
    let source = include_str!("../../../../examples/flue-network.performance")
        .split_once("\ngesture\n")
        .unwrap()
        .0
        .replace("output 24000", "output 1024")
        .replace("jet 90 1e-6 3 0\n", &format!("jet 90 1e-6 3 0\n{rows}"));
    let source = if bare { source.replace("windway slit 0.015\n", "") } else { source };
    format!(
        "{source}\ngesture\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t1000\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t3e2\nevents\t0\n"
    )
}

#[test]
fn auvray_jet_transfer_is_opt_in_by_name_changes_the_render_and_refuses_beside_another_jet_row() {
    use file::FlueFilePerformance;
    use fs_exec::CancelGate;
    let gate = CancelGate::new_clock_free();
    let parse = |text: &str| FlueFilePerformance::from_bytes(text.as_bytes(), 256, &gate);
    let row = "jet_transfer auvray-2014 0.3 0.4 numerical-fixture\n";
    let mut bare = parse(&auvray_input("", true)).unwrap().into_renderer();
    let mut named = parse(&auvray_input(row, true)).unwrap().into_renderer();
    assert!(bare.voice().auvray_jet().is_none());
    let jet = named.voice().auvray_jet().expect("the named row attaches");
    assert_eq!((jet.spec().beta, jet.spec().gamma, jet.spec().source.as_str()),
        (0.3, 0.4, "numerical-fixture"));
    assert!(named.voice().forced_jet().is_none() && named.voice().nolle_jet().is_none());
    let mut plain = vec![0.0; 1024];
    let mut driven = vec![0.0; 1024];
    for chunk in plain.chunks_mut(256) { bare.block(chunk).unwrap(); }
    for chunk in driven.chunks_mut(256) { named.block(chunk).unwrap(); }
    let energy: f64 = plain.iter().map(|p| p * p).sum();
    let difference: f64 = driven.iter().zip(&plain).map(|(a, b)| (a - b).powi(2)).sum();
    assert!(energy > 0.0);
    assert!(difference > 1e-10 * energy, "naming the row must change the render");
    // Refused, each with its own text: behind the slit windway, beside the forced jet or another transfer
    // in either order, and as a malformed or non-positive row.
    let refusal = |text: &str| match parse(text) {
        Ok(_) => panic!("accepted:\n{text}"),
        Err(e) => e.to_string(),
    };
    let windway = refusal(&auvray_input(row, false));
    assert!(windway.contains("refuses every windway law"), "{windway}");
    let one_row = "one jet_transfer or jet_response row at most: the jet transfers exclude each other";
    let beta = "Auvray jet growth coefficient beta must be finite and positive";
    let gamma = "Auvray jet convection ratio gamma must lie in (0, 1]";
    for (refused, text) in [
        (format!("{row}jet_response forced-bickley 0.1 90\n"), "jet_transfer and jet_response exclude each other"),
        (format!("jet_response forced-bickley 0.1 90\n{row}"), one_row),
        (format!("{row}jet_transfer nolle-price-2015 long-channel 0\n"), one_row),
        (format!("jet_transfer nolle-price-2015 long-channel 0\n{row}"), one_row),
        (format!("{row}{row}"), one_row),
        ("jet_transfer auvray-2014 0.3 0.4\n".to_owned(),
            "jet_transfer auvray-2014 needs an explicit source tag after BETA GAMMA"),
        ("jet_transfer auvray-2014 0.3\n".to_owned(), "missing or invalid numeric field"),
        ("jet_transfer auvray-2014 0.3 0.4 tag extra\n".to_owned(), "unexpected trailing fields"),
        ("jet_transfer auvray-2014 0 0.4 tag\n".to_owned(), beta),
        ("jet_transfer auvray-2014 -0.3 0.4 tag\n".to_owned(), beta),
        ("jet_transfer auvray-2014 0.3 0 tag\n".to_owned(), gamma),
        ("jet_transfer auvray-2014 0.3 -0.4 tag\n".to_owned(), gamma),
        ("jet_transfer auvray-2014 0.3 1.5 tag\n".to_owned(), gamma),
        ("jet_transfer auvray-2014 NaN 0.4 tag\n".to_owned(), "physical scalars must be finite"),
        ("jet_transfer auvray 0.3 0.4 tag\n".to_owned(),
            "jet_transfer requires nolle-price-2015 CHANNEL CHAMFER_DEPTH_M or auvray-2014 BETA GAMMA SOURCE_TAG"),
    ] {
        let got = refusal(&auvray_input(&refused, true));
        assert!(got.contains(text), "{refused}: refused as `{got}`, expected `{text}`");
    }
}

mod gas;

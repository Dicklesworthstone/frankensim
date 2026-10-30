//! G0/G3 native composition of accepted resolved escape flow and exterior Pa.
//! Synthetic geometry checks ownership and signs, not measured flute fidelity.
use fs_couple::nsf_dg::performance::ReservoirObservation;
use fs_couple::nsf_dg::performance::file::ResolvedDgFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
use std::sync::Arc;

const ROW: &str = "reservoir_observation baffled -0.08 0 0\n";

fn source(gas: bool) -> String {
    let mut text = include_str!("../examples/resolved-dg-side-aperture.performance")
        .replace(ROW, "")
        .replace("output 2048 20", "output 256 20")
        .replace("hold 128 20", "hold 64 20")
        .replace("hold 256 -20", "hold 128 -20")
        .replace("hold 2048 0", "hold 256 0")
        .replace(
            "reference 1\nreservoir_cell 5 2 1",
            "reference 3\nreservoir_cell 5 2 1\nreservoir_cell 3 2 1\nreservoir_cell 4 2 1",
        );
    if !gas {
        text = text.replace("observer_gas stokes-kirchhoff 4096 1000000 50000000\n", "");
    }
    text
}

fn observed(text: &str) -> String {
    text.replace(
        "observation baffled 0 0 0.2 8 32 2000\n",
        &format!("observation baffled 0 0 0.2 8 32 2000\n{ROW}"),
    )
}

fn read(text: &str, block: usize) -> ResolvedDgFilePerformance {
    ResolvedDgFilePerformance::from_bytes(
        text.as_bytes(),
        block,
        Arc::new(CancelGate::new_clock_free()),
    )
    .unwrap()
}

#[test]
fn every_accepted_escape_tile_reaches_the_shared_receiver_without_changing_source_state() {
    let text = source(false);
    let mut baseline = read(&text, 31).into_renderer();
    let mut observed = read(&observed(&text), 31).into_renderer();
    let mut bore = observed.baffled_receiver().unwrap().clone();
    let mut escape = observed.reservoir_receiver().unwrap().clone();
    let mut actual = [0.0];
    let mut old = [0.0];
    let mut escaped_peak = 0.0_f64;
    let mut distinct_tile_flow = false;
    for _ in 0..256 {
        baseline.block(&mut old).unwrap();
        observed.block(&mut actual).unwrap();
        let frame = observed.last_frame().unwrap();
        assert_eq!(frame.reservoir_flows.len(), 3);
        assert_eq!(
            frame
                .reservoir_flows
                .iter()
                .map(|flow| flow.cell)
                .collect::<Vec<_>>(),
            [[3, 2, 1], [4, 2, 1], [5, 2, 1]]
        );
        let q: Vec<_> = frame
            .terminals
            .iter()
            .map(|flow| flow.outward_flow_m3_s)
            .collect();
        let escaped_q: Vec<_> = frame
            .reservoir_flows
            .iter()
            .map(|flow| flow.outward_volume_flow_m3_s)
            .collect();
        distinct_tile_flow |= escaped_q.windows(2).any(|pair| pair[0] != pair[1]);
        let bore_pa = bore.step_trace(&q).unwrap();
        let escape_pa = escape.step_trace(&escaped_q).unwrap();
        assert_eq!(actual[0], bore_pa + escape_pa);
        assert_eq!(old[0], bore_pa);
        escaped_peak = escaped_peak.max(escape_pa.abs());
        assert_eq!(
            baseline.system().zone().state_hash(),
            observed.system().zone().state_hash()
        );
        assert_eq!(
            baseline.system().load().state(),
            observed.system().load().state()
        );
        assert_eq!(
            baseline.last_frame().unwrap().reservoir_flows,
            frame.reservoir_flows
        );
    }
    assert!(
        distinct_tile_flow,
        "the receiver must not collapse actual spatially varying flow"
    );
    assert!(
        escaped_peak > 1e-9,
        "the physical escape source must reach emitted pressure"
    );
}

#[test]
fn gas_paths_and_projected_clocks_preserve_all_tile_histories_across_callbacks() {
    // The physical grid stays fixed. At 48 kHz each output integrates two
    // original 96 kHz intervals before the same shared differentiator.
    let text = observed(&source(true))
        .replace(
            "wall adiabatic\n",
            "wall adiabatic\nsgs exploration wale 0.0001 1 100 0.9\n",
        )
        .replace("clock 96000 1", "clock 48000 1")
        .replace("output 256 20", "output 128 20")
        .replace("hold 64 20", "hold 32 20")
        .replace("hold 128 -20", "hold 64 -20")
        .replace("hold 256 0", "hold 128 0");
    let first = read(&text, 37);
    let sgs = first.renderer().system().zone().sgs_config().unwrap();
    assert_eq!(first.info().sgs, Some(sgs));
    assert!(matches!(
        sgs.closure.model,
        fs_lbm::sgs3::SgsModel::Wale { cw } if cw == 0.0001
    ));
    let bore = first.renderer().observer_gas().unwrap();
    let escape = first.renderer().reservoir_observer_gas().unwrap();
    assert_eq!(bore.gas, escape.gas);
    assert_eq!(bore.law, escape.law);
    assert_eq!(
        bore.additional_latency_samples,
        escape.additional_latency_samples
    );
    assert_eq!(escape.filters.len(), 3);
    let run = |file: ResolvedDgFilePerformance, block| {
        let count = file.info().samples as usize;
        let mut renderer = file.into_renderer();
        let mut pressure = vec![0.0; count];
        for chunk in pressure.chunks_mut(block) {
            renderer.block(chunk).unwrap();
        }
        (pressure, renderer)
    };
    let (whole, a) = run(first, 37);
    let (parts, b) = run(read(&text, 7), 7);
    assert_eq!(whole, parts);
    assert!(whole.iter().any(|value| value.abs() > 1e-10));
    assert_eq!(
        a.system().zone().state_hash(),
        b.system().zone().state_hash()
    );
    assert_eq!(a.system().load().state(), b.system().load().state());
    assert_eq!(a.system().accepted_steps(), 256);
}

#[test]
fn incompatible_or_late_reservoir_observation_refuses_before_an_output() {
    let text = source(false);
    for invalid in [
        observed(&text).replace(
            "observation baffled 0 0 0.2 8 32 2000",
            "observation inlet-pressure",
        ),
        observed(&text).replace(
            ROW,
            "reservoir_observation baffled -0.0079418628 -0.0238255884 0\n",
        ),
        observed(&text).replace("dg_exterior 300 100000 0", "dg_exterior 310 100000 0"),
        observed(&text).replace(ROW, "reservoir_observation baffled 0 0 NaN\n"),
        observed(&text).replace(ROW, "reservoir_observation baffled 0 0 0 extra\n"),
    ] {
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                invalid.as_bytes(),
                1,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err()
        );
    }
    let placement = ReservoirObservation {
        origin_m: [-0.08, 0.0],
        angle_rad: 0.0,
    };
    let with_gas = read(&source(true), 1).into_renderer();
    assert!(with_gas.with_reservoir_observation(placement).is_err());
    let duplicate = read(&observed(&text), 1).into_renderer();
    assert!(duplicate.with_reservoir_observation(placement).is_err());
    let mut late = read(&text, 1).into_renderer();
    late.block(&mut [0.0]).unwrap();
    assert!(late.with_reservoir_observation(placement).is_err());
}

// Same physical grid as the observation fixture, with twelve original NSF
// steps per 8 kHz acoustic interval. The complete mixed radiation source is
// then resolved by the original terminal surface quadrature over its 4 kHz band.
const FEEDBACK: &str = "reservoir_feedback -20 20 16 16\n";

fn feedback_source(gas: bool) -> String {
    observed(&source(gas))
        .replace(ROW, &format!("{ROW}{FEEDBACK}"))
        .replace("clock 96000 1", "clock 8000 12")
        .replace("output 256 20", "output 64 20")
        .replace("hold 64 20", "hold 16 2")
        .replace("hold 128 -20", "hold 32 -2")
        .replace("hold 256 0", "hold 64 0")
        .replace("coupling 64 9216", "coupling 256 442368")
}

#[test]
fn resolved_mouth_reaction_changes_original_flow_and_audits_both_acoustic_ports() {
    let text = feedback_source(false);
    let mut baseline = read(&text.replace(FEEDBACK, ""), 1).into_renderer();
    let mut coupled = read(&text, 1).into_renderer();
    let relation = coupled
        .system()
        .load()
        .independent_mouth_response()
        .unwrap();
    assert_ne!(relation.impedance_pa_s_m3[0][1], 0.0);
    assert_ne!(relation.impedance_pa_s_m3[1][0], 0.0);
    let mut independent_flow = false;
    let mut changed_flow = false;
    let mut reacted = false;
    let mut distinct_fv_work = false;
    for _ in 0..64 {
        let mass = coupled.system().zone().total_mass_kg();
        let energy = coupled.system().zone().total_energy_j();
        baseline.block(&mut [0.0]).unwrap();
        coupled.block(&mut [0.0]).unwrap();
        let frame = coupled.last_frame().unwrap();
        let mouth = frame.reservoir_radiation.as_ref().unwrap();
        let dt = coupled.system().load().propagator().dt();
        let resolved = &mouth.resolved_exchange;
        let tile_q: f64 = frame.reservoir_flows.iter()
            .map(|tile| tile.outward_volume_flow_m3_s).sum();
        let tile_q_abs: f64 = frame.reservoir_flows.iter()
            .map(|tile| tile.outward_volume_flow_m3_s.abs()).sum();
        let q = resolved.normal_volume_flow_m3_s;
        let resolved_work = resolved.acoustic_pressure_work_out_j;
        let trace_work: f64 = frame.reservoir_flows.iter()
            .map(|tile| tile.acoustic_pressure_work_out_j).sum();
        let trace_absolute_work: f64 = frame.reservoir_flows.iter()
            .map(|tile| tile.pressure_work_out_j.abs()).sum();
        let flow_scale = tile_q_abs + q.abs()
            + coupled.system().reference_gas().sound_speed * resolved.area_m2;
        assert!((tile_q - q).abs() <= 64.0 * f64::EPSILON * flow_scale);
        assert_eq!(mouth.outward_flow_m3_s, q);
        assert_eq!(mouth.mean_acoustic_pressure_pa, resolved.mean_acoustic_pressure_pa);
        assert_eq!(resolved.reference_pressure_pa,
            coupled.system().reference_gas().pressure);
        // Absolute-pressure scale includes centering and cancellation roundoff;
        // the witness must distinguish physical quadratures beyond those errors.
        let work_scale = resolved.pressure_work_out_j.abs() + trace_absolute_work
            + resolved.reference_pressure_pa * (flow_scale * dt);
        distinct_fv_work |= (resolved_work - trace_work).abs()
            > 256.0 * f64::EPSILON * work_scale;
        assert!(frame.pressure_residual_pa.abs() < 1e-5);
        assert!(mouth.pressure_residual_pa.abs() < 1e-5);
        assert_eq!(
            mouth.acoustic_work_residual_j,
            resolved_work - mouth.radiation_work_j
        );
        let projected_work = mouth.mean_acoustic_pressure_pa * (q * dt);
        assert_eq!(
            mouth.projection_covariance_j,
            resolved_work - projected_work
        );
        assert!(
            (mouth.radiation_work_j - mouth.radiation_pressure_pa * (q * dt)).abs()
                <= 64.0 * f64::EPSILON * mouth.radiation_work_j.abs().max(1e-30)
        );
        // Mass and total energy still use the original conservative accounts,
        // never the acoustic projection or reference-density volume flow.
        let mass_residual = frame.flow_mass_change_kg
            + frame.reservoir_mass_out_kg
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!(mass_residual.abs() < 5e-12 * mass);
        let supplied = frame.flow.open_energy_j
            + frame.flow.open_boundary_energy_j
            + frame.flow.wall_heat_j
            + frame.flow.source_energy_j;
        assert!((frame.flow_energy_change_j - supplied).abs() < 5e-11 * energy.abs());
        independent_flow |=
            (q + frame.inlet_flow_m3_s).abs() > 1e-6 * (q.abs() + frame.inlet_flow_m3_s.abs());
        changed_flow |= frame.inlet_flow_m3_s != baseline.last_frame().unwrap().inlet_flow_m3_s;
        reacted |= mouth.radiation_pressure_pa.abs() > 1e-9;
    }
    assert!(distinct_fv_work,
        "the fixture must distinguish original FV work from endpoint trace products");
    assert!(
        independent_flow,
        "resolved storage must not impose Q_mouth = -Q_bore"
    );
    assert!(
        changed_flow,
        "computed exterior reaction must reach the original source"
    );
    assert!(reacted);
}

#[test]
fn feedback_with_finite_upstream_and_material_observation_replays_across_callbacks() {
    let text = feedback_source(true)
        .replace("clock 8000 12", "clock 12000 12")
        .replace("output 64 20", "output 96 20")
        .replace("hold 16 2", "hold 24 2")
        .replace("hold 32 -2", "hold 48 -2")
        .replace("hold 64 0", "hold 96 0")
        .replace("observation baffled 0 0 0.2 8 32 2000",
            "upstream_dg dg_tube 0.03 0.00896143253443043 closed 1 2 32 50 2000 0.01 1e-11 512\nobservation baffled 0 0 0.2 8 32 2000")
        .replace("boundary pressure 3", "boundary upstream-pressure 3");
    let first = read(&text, 13);
    let bore = first.renderer().observer_gas().unwrap();
    let escape = first.renderer().reservoir_observer_gas().unwrap();
    assert_eq!(bore.gas, escape.gas);
    assert_eq!(bore.law, escape.law);
    assert_eq!(
        bore.additional_latency_samples,
        escape.additional_latency_samples
    );
    let run = |file: ResolvedDgFilePerformance, block| {
        let mut pressure = vec![0.0; file.info().samples as usize];
        let mut renderer = file.into_renderer();
        for chunk in pressure.chunks_mut(block) {
            renderer.block(chunk).unwrap();
        }
        (pressure, renderer)
    };
    let (a, whole) = run(first, 13);
    let (b, parts) = run(read(&text, 5), 5);
    assert_eq!(a, b);
    assert!(a.iter().any(|value| value.abs() > 1e-12));
    assert_eq!(
        whole.system().zone().state_hash(),
        parts.system().zone().state_hash()
    );
    assert_eq!(whole.system().load().state(), parts.system().load().state());
    assert_eq!(
        whole.system().upstream_load().unwrap().state(),
        parts.system().upstream_load().unwrap().state()
    );
    assert_eq!(whole.system().accepted_steps(), 64);
    let frame = whole.last_frame().unwrap();
    assert!(
        frame
            .reservoir_radiation
            .as_ref()
            .unwrap()
            .pressure_residual_pa
            .abs()
            < 1e-5
    );
    assert!(frame.upstream.as_ref().unwrap().pressure_residual_pa.abs() < 1e-5);
    assert!(frame.pressure_residual_pa.abs() < 1e-5);
}

#[test]
fn feedback_geometry_domains_and_cancelled_or_unbracketed_solves_preserve_owners() {
    let text = feedback_source(false);
    for invalid in [
        text.replace(FEEDBACK, "reservoir_feedback 1 1 16 16\n"),
        text.replace(FEEDBACK, "reservoir_feedback -20 20 0 16\n"),
        text.replace(ROW, "reservoir_observation baffled -0.08 0 0.2\n"),
        // The original aperture remains a valid pair of fluid tiles, but its
        // missing middle tile cannot silently become a filled rectangular port.
        text.replace("reference 3\nreservoir_cell", "reference 2\nreservoir_cell")
            .replace("reservoir_cell 4 2 1\n", ""),
        text.replace(ROW, ""),
        text.replace("boundary pressure 3", "boundary reservoir 3")
            .replace("hold 16 2", "hold 16 100002 300")
            .replace("hold 32 -2", "hold 32 99998 300")
            .replace("hold 64 0", "hold 64 100000 300"),
    ] {
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                invalid.as_bytes(),
                1,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err()
        );
    }
    let gate = Arc::new(CancelGate::new_clock_free());
    let file =
        ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 2, Arc::clone(&gate)).unwrap();
    let mut renderer = file.into_renderer();
    renderer.block(&mut [0.0; 2]).unwrap();
    let before = (
        renderer.system().zone().state_hash(),
        renderer.system().zone().cumulative_ledger(),
        renderer.system().zone().last_reservoir_flows().to_vec(),
        renderer.system().zone().last_reservoir_exchange().copied(),
        renderer.system().load().state().clone(),
        renderer.system().accepted_steps(),
    );
    gate.request();
    assert!(renderer.block(&mut [0.0]).is_err());
    assert_eq!(
        before,
        (
            renderer.system().zone().state_hash(),
            renderer.system().zone().cumulative_ledger(),
            renderer.system().zone().last_reservoir_flows().to_vec(),
            renderer.system().zone().last_reservoir_exchange().copied(),
            renderer.system().load().state().clone(),
            renderer.system().accepted_steps()
        )
    );

    let unreachable = text.replace(FEEDBACK, "reservoir_feedback 500 600 16 16\n");
    let mut refused = read(&unreachable, 1).into_renderer();
    let before = (
        refused.system().zone().state_hash(),
        refused.system().zone().last_reservoir_exchange().copied(),
        refused.system().load().state().clone(),
    );
    assert!(refused.block(&mut [0.0]).is_err());
    assert_eq!(refused.system().accepted_steps(), 0);
    assert_eq!(
        before,
        (
            refused.system().zone().state_hash(),
            refused.system().zone().last_reservoir_exchange().copied(),
            refused.system().load().state().clone()
        )
    );
}


fn separate_opening_source(gas: bool, upstream: bool) -> String {
    let mut text = feedback_source(gas)
        .replace("reservoir_aperture 1 high reference 3",
            "reservoir_apertures 2\nreservoir_aperture 1 low reference 1\nreservoir_cell 5 0 1\nreservoir_aperture 1 high reference 3")
        .replace(ROW, "reservoir_observation aperture 1 baffled -0.08 0 0\n")
        .replace("boundary pressure 3", "boundary aperture 0 3")
        .replace("hold 16 2", "hold 16 100002 300.2")
        .replace("hold 32 -2", "hold 32 99998 300.1")
        .replace("hold 64 0", "hold 64 100000 300");
    if upstream {
        text = text.replace("observation baffled 0 0 0.2 8 32 2000",
            "upstream_dg dg_tube 0.03 0.00896143253443043 closed 1 2 32 50 2000 0.01 1e-11 512\nobservation baffled 0 0 0.2 8 32 2000");
    }
    text
}

#[test]
fn separate_breath_and_escape_retain_conservative_accounts_and_selected_sound() {
    let text = separate_opening_source(false, false);
    let mut renderer = read(&text, 1).into_renderer();
    assert_eq!(renderer.system().zone().reservoir_aperture_count(), 2);
    assert_eq!(renderer.system().reservoir_radiation_aperture(), Some(1));
    assert_eq!(renderer.reservoir_observation_aperture(), Some(1));
    let mut bore = renderer.baffled_receiver().unwrap().clone();
    let mut escape = renderer.reservoir_receiver().unwrap().clone();
    let mut reacted = false;
    let mut independent_mass = false;
    for index in 0..64 {
        let mass = renderer.system().zone().total_mass_kg();
        let energy = renderer.system().zone().total_energy_j();
        let mut pressure = [0.0];
        renderer.block(&mut pressure).unwrap();
        let frame = renderer.last_frame().unwrap();
        let mouth = frame.reservoir_radiation.as_ref().unwrap();
        let accounts = &frame.reservoir_exchanges;
        assert_eq!(accounts.len(), 2);
        assert_eq!(accounts[0].aperture, 0);
        assert_eq!(accounts[1].aperture, 1);
        assert_eq!(mouth.resolved_exchange, accounts[1]);
        assert!(matches!(frame.drive,
            fs_couple::nsf_dg::NsfDgDrive::ReservoirAperture { aperture: 0, .. }));
        let tiles: Vec<_> = frame.reservoir_flows.iter().filter(|tile| tile.aperture == 1)
            .map(|tile| tile.outward_volume_flow_m3_s).collect();
        assert_eq!(tiles.len(), 3);
        assert_eq!(frame.reservoir_flows.iter().filter(|tile| tile.aperture == 0).count(), 1);
        let terminals: Vec<_> = frame.terminals.iter().map(|port| port.outward_flow_m3_s).collect();
        assert_eq!(pressure[0], bore.step_trace(&terminals).unwrap() + escape.step_trace(&tiles).unwrap());
        let dt = renderer.system().load().propagator().dt();
        assert_eq!(mouth.mass_projection_difference_kg,
            accounts[1].mass_out_kg - renderer.system().reference_gas().density
                * (accounts[1].normal_volume_flow_m3_s * dt));
        assert_eq!(mouth.acoustic_work_residual_j,
            accounts[1].acoustic_pressure_work_out_j - mouth.radiation_work_j);
        let reservoir_mass: f64 = accounts.iter().map(|owner| owner.mass_out_kg).sum();
        let reservoir_energy: f64 = accounts.iter().map(|owner| owner.total_energy_out_j).sum();
        assert!((reservoir_mass - frame.reservoir_mass_out_kg).abs() < 1e-13 * mass);
        assert!((reservoir_energy + frame.flow.open_energy_j).abs() < 1e-13 * energy.abs());
        let mass_residual = frame.flow_mass_change_kg + reservoir_mass
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!(mass_residual.abs() < 5e-12 * mass);
        let supplied = frame.flow.open_energy_j + frame.flow.open_boundary_energy_j
            + frame.flow.wall_heat_j + frame.flow.source_energy_j;
        assert!((frame.flow_energy_change_j - supplied).abs() < 5e-11 * energy.abs());
        assert!(frame.pressure_residual_pa.abs() < 1e-5);
        assert!(mouth.pressure_residual_pa.abs() < 1e-5);
        let expected_temperature = if index < 16 { 300.2 } else if index < 32 { 300.1 } else { 300.0 };
        assert_eq!(renderer.system().zone().reservoir_aperture_at(0).unwrap().state.temperature_k,
            expected_temperature);
        assert_eq!(renderer.system().zone().reservoir_aperture_at(1).unwrap().state.temperature_k, 300.0);
        reacted |= mouth.radiation_pressure_pa.abs() > 1e-9;
        independent_mass |= accounts[0].mass_out_kg != accounts[1].mass_out_kg;
    }
    assert!(reacted);
    assert!(independent_mass);
}

#[test]
fn separate_breath_escape_and_cork_duct_replay_and_cancel_as_one_system() {
    let text = separate_opening_source(true, true);
    let run = |block| {
        let mut renderer = read(&text, block).into_renderer();
        let mut output = vec![0.0; 64];
        for chunk in output.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        (output, renderer)
    };
    let (whole, a) = run(16);
    let (parts, b) = run(3);
    assert_eq!(whole, parts);
    assert!(whole.iter().any(|value| value.abs() > 1e-12));
    assert_eq!(a.system().zone().state_hash(), b.system().zone().state_hash());
    assert_eq!(a.system().load().state(), b.system().load().state());
    assert_eq!(a.system().upstream_load().unwrap().state(), b.system().upstream_load().unwrap().state());
    assert_eq!(a.last_frame().unwrap().reservoir_exchanges, b.last_frame().unwrap().reservoir_exchanges);
    let frame = a.last_frame().unwrap();
    assert!(frame.upstream.as_ref().unwrap().pressure_residual_pa.abs() < 1e-5);
    assert_eq!(frame.upstream.as_ref().unwrap().source_work_j, 0.0);
    assert_eq!(a.observer_gas().unwrap().additional_latency_samples,
        a.reservoir_observer_gas().unwrap().additional_latency_samples);

    let gate = Arc::new(CancelGate::new_clock_free());
    let mut renderer = ResolvedDgFilePerformance::from_bytes(
        text.as_bytes(), 1, Arc::clone(&gate)).unwrap().into_renderer();
    renderer.block(&mut [0.0]).unwrap();
    let before = (renderer.system().zone().state_hash(),
        renderer.system().zone().last_reservoir_exchanges().to_vec(),
        renderer.system().load().state().clone(),
        renderer.system().upstream_load().unwrap().state().clone());
    gate.request();
    assert!(renderer.block(&mut [0.0]).is_err());
    assert_eq!(before, (renderer.system().zone().state_hash(),
        renderer.system().zone().last_reservoir_exchanges().to_vec(),
        renderer.system().load().state().clone(),
        renderer.system().upstream_load().unwrap().state().clone()));
}

#[test]
fn unbracketed_escape_never_publishes_its_separate_breath_control() {
    let text = separate_opening_source(false, true)
        .replace(FEEDBACK, "reservoir_feedback 500 600 16 16\n");
    let mut renderer = read(&text, 1).into_renderer();
    let before = (renderer.system().zone().state_hash(),
        renderer.system().zone().reservoir_aperture_at(0).unwrap().state,
        renderer.system().zone().reservoir_aperture_at(1).unwrap().state,
        renderer.system().load().state().clone(),
        renderer.system().upstream_load().unwrap().state().clone());
    assert!(renderer.block(&mut [0.0]).is_err());
    assert_eq!(renderer.system().accepted_steps(), 0);
    assert_eq!(before, (renderer.system().zone().state_hash(),
        renderer.system().zone().reservoir_aperture_at(0).unwrap().state,
        renderer.system().zone().reservoir_aperture_at(1).unwrap().state,
        renderer.system().load().state().clone(),
        renderer.system().upstream_load().unwrap().state().clone()));
}


#[test]
fn both_characteristic_patches_keep_separate_breath_and_escape_in_the_same_trial() {
    let text = separate_opening_source(false, true).replace("reservoir_apertures 2",
        "coupled_patch 0 3 0 3 pressure 0\nsource_patch 0 3 0 3 pressure 0\nreservoir_apertures 2");
    let mut renderer = read(&text, 1).into_renderer();
    assert!(renderer.system().coupled_patch().is_some());
    assert!(renderer.system().source_patch().is_some());
    for _ in 0..8 {
        let mass = renderer.system().zone().total_mass_kg();
        renderer.block(&mut [0.0]).unwrap();
        let frame = renderer.last_frame().unwrap();
        let main = frame.coupled_patch.as_ref().unwrap();
        let source = frame.source_patch.as_ref().unwrap();
        let upstream = frame.upstream.as_ref().unwrap();
        let escape = frame.reservoir_radiation.as_ref().unwrap();
        assert_eq!(main.normal_volume_flow_m3_s, frame.inlet_flow_m3_s);
        assert_eq!(source.normal_volume_flow_m3_s, upstream.inlet_flow_m3_s);
        assert_eq!(frame.reservoir_exchanges.len(), 2);
        assert_eq!(escape.resolved_exchange, frame.reservoir_exchanges[1]);
        assert_eq!(upstream.source_work_j, 0.0);
        assert!(frame.pressure_residual_pa.abs() < 1e-5);
        assert!(upstream.pressure_residual_pa.abs() < 1e-5);
        assert!(escape.pressure_residual_pa.abs() < 1e-5);
        // Selected patch accounts are part of the full axial faces, not extra
        // transport to add a second time to the two transverse openings.
        let mass_residual = frame.flow_mass_change_kg + frame.reservoir_mass_out_kg
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!(mass_residual.abs() < 5e-12 * mass);
    }
}


const FIRST_ESCAPE: &str = "reservoir_observation aperture 0 baffled 0.06 0 0\n";
const SECOND_ESCAPE: &str = "reservoir_observation aperture 1 baffled -0.08 0 0\n";

fn multiple_escape_source(gas: bool) -> String {
    source(gas)
        .replace("reservoir_aperture 1 high reference 3",
            "reservoir_apertures 2\nreservoir_aperture 1 low reference 1\nreservoir_cell 6 0 1\nreservoir_aperture 1 high reference 3")
        .replace("observation baffled 0 0 0.2 8 32 2000\n",
            &format!("observation baffled 0 0 0.2 8 32 2000\n{FIRST_ESCAPE}{SECOND_ESCAPE}"))
}

#[test]
fn multiple_escape_surfaces_sum_each_signed_tile_once_at_its_original_distance() {
    let text = multiple_escape_source(false);
    let mut observed = read(&text, 1).into_renderer();
    let mut baseline = read(&text.replace(FIRST_ESCAPE, "").replace(SECOND_ESCAPE, ""), 1)
        .into_renderer();
    assert_eq!(observed.reservoir_observation_count(), 2);
    assert!(observed.reservoir_observation().is_none());
    assert!(observed.reservoir_observation_aperture().is_none());
    assert!(observed.reservoir_receiver().is_none());
    let mut bore = observed.baffled_receiver().unwrap().clone();
    let mut first = observed.reservoir_receiver_at(0).unwrap().clone();
    let mut second = observed.reservoir_receiver_at(1).unwrap().clone();
    let mut averaged = second.clone();
    assert_eq!(first.input_count(), 1);
    assert_eq!(second.input_count(), 3);
    let mut peaks = [0.0_f64; 2];
    let mut averaging_error = 0.0_f64;
    let mut signed = [false; 2];
    for _ in 0..256 {
        let mut output = [0.0];
        let mut old = [0.0];
        observed.block(&mut output).unwrap();
        baseline.block(&mut old).unwrap();
        let frame = observed.last_frame().unwrap();
        let terminals: Vec<_> = frame.terminals.iter().map(|port| port.outward_flow_m3_s).collect();
        let flows = |aperture| frame.reservoir_flows.iter()
            .filter(|tile| tile.aperture == aperture)
            .map(|tile| tile.outward_volume_flow_m3_s).collect::<Vec<_>>();
        let q0 = flows(0);
        let q1 = flows(1);
        assert_eq!(q0.len(), 1);
        assert_eq!(q1.len(), 3);
        for q in q0.iter().chain(&q1) {
            signed[0] |= *q < 0.0;
            signed[1] |= *q > 0.0;
        }
        let p0 = first.step_trace(&q0).unwrap();
        let p1 = second.step_trace(&q1).unwrap();
        let bore_pa = bore.step_trace(&terminals).unwrap();
        assert_eq!(output[0], (bore_pa + p0) + p1);
        assert_eq!(old[0], bore_pa);
        let merged = q1.iter().sum::<f64>() / q1.len() as f64;
        let wrong = averaged.step_trace(&vec![merged; q1.len()]).unwrap();
        averaging_error = averaging_error.max((p1 - wrong).abs());
        peaks[0] = peaks[0].max(p0.abs());
        peaks[1] = peaks[1].max(p1.abs());
        assert_eq!(observed.system().zone().state_hash(), baseline.system().zone().state_hash());
        assert_eq!(observed.system().load().state(), baseline.system().load().state());
    }
    assert!(peaks.into_iter().all(|peak| peak > 1e-9));
    assert!(signed.into_iter().all(|present| present), "both outward and inward flow must remain signed");
    assert!(averaging_error > 1e-12, "spatially resolved flow cannot be replaced by an equal-flow piston");

    let farther = text.replace(SECOND_ESCAPE,
        "reservoir_observation aperture 1 baffled -0.16 0 0\n");
    let near = read(&text, 16).into_renderer();
    let far = read(&farther, 16).into_renderer();
    assert!(far.reservoir_receiver_at(1).unwrap().delay_samples.0
        > near.reservoir_receiver_at(1).unwrap().delay_samples.0);
    assert_eq!(far.reservoir_receiver_at(0).unwrap().delay_samples,
        near.reservoir_receiver_at(0).unwrap().delay_samples);
}

#[test]
fn multiple_escape_material_paths_keep_one_latency_and_projected_callback_history() {
    let text = multiple_escape_source(true)
        .replace("clock 96000 1", "clock 48000 1")
        .replace("output 256 20", "output 128 20")
        .replace("hold 64 20", "hold 32 20")
        .replace("hold 128 -20", "hold 64 -20")
        .replace("hold 256 0", "hold 128 0");
    let first = read(&text, 37);
    let bore = first.renderer().observer_gas().unwrap();
    for (aperture, paths) in [(0, 1), (1, 3)] {
        let report = first.renderer().reservoir_observer_gas_at(aperture).unwrap();
        assert_eq!(report.gas, bore.gas);
        assert_eq!(report.law, bore.law);
        assert_eq!(report.additional_latency_samples, bore.additional_latency_samples);
        assert_eq!(report.filters.len(), paths);
    }
    assert!(first.renderer().reservoir_observer_gas().is_none());
    let run = |file: ResolvedDgFilePerformance, block| {
        let mut pressure = vec![0.0; file.info().samples as usize];
        let mut renderer = file.into_renderer();
        for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        (pressure, renderer)
    };
    let (a, whole) = run(first, 37);
    let (b, parts) = run(read(&text, 7), 7);
    assert_eq!(a, b);
    assert!(a.iter().any(|value| value.abs() > 1e-10));
    assert_eq!(whole.system().accepted_steps(), 256);
    assert_eq!(whole.system().zone().state_hash(), parts.system().zone().state_hash());
    assert_eq!(whole.system().load().state(), parts.system().load().state());

    let changed = text.replace("gas dry-air 300 100000", "gas dry-air 310 100000")
        .replace("dg_exterior 300 100000 0", "dg_exterior 310 100000 0");
    let warmer = read(&changed, 7);
    for aperture in 0..2 {
        assert_ne!(warmer.renderer().reservoir_observer_gas_at(aperture).unwrap().gas,
            whole.reservoir_observer_gas_at(aperture).unwrap().gas);
        assert_eq!(warmer.renderer().reservoir_observer_gas_at(aperture).unwrap().gas,
            warmer.renderer().observer_gas().unwrap().gas);
    }
    let (warm_pressure, warm) = run(warmer, 7);
    assert_ne!(a, warm_pressure);
    assert_eq!(warm.system().zone().law().dx_m, whole.system().zone().law().dx_m);
}

#[test]
fn multiple_escape_selection_refuses_duplicates_overlap_and_ambiguous_feedback() {
    let text = multiple_escape_source(false);
    let dx = read(&text, 1).renderer().system().zone().law().dx_m;
    // Place the low-face tile onto the first high-face tile, including a
    // non-axis-aligned case. This tests actual tile area, not source labels.
    let overlap = |angle: f64| {
        let x = -0.08 + 1.5 * dx;
        let y = 3.5 * dx;
        let origin_x = x - (angle.cos() * 1.5 + angle.sin() * 6.5) * dx;
        let origin_y = y - (angle.sin() * 1.5 - angle.cos() * 6.5) * dx;
        text.replace(FIRST_ESCAPE, &format!(
            "reservoir_observation aperture 0 baffled {origin_x:.17e} {origin_y:.17e} {angle:.17e}\n"))
    };
    for invalid in [
        text.replace(SECOND_ESCAPE, FIRST_ESCAPE),
        text.replace(SECOND_ESCAPE, "reservoir_observation aperture 2 baffled -0.08 0 0\n"),
        text.replace(SECOND_ESCAPE, &format!("{SECOND_ESCAPE}reservoir_feedback -20 20 16 16\n")),
        text.replace(SECOND_ESCAPE, &format!("{SECOND_ESCAPE}reservoir_feedback aperture 2 -20 20 16 16\n")),
        overlap(0.0), overlap(0.47),
    ] {
        assert!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 1,
            Arc::new(CancelGate::new_clock_free())).is_err());
    }

    let feedback = text
        .replace("clock 96000 1", "clock 8000 12")
        .replace("output 256 20", "output 64 20")
        .replace("hold 64 20", "hold 16 2")
        .replace("hold 128 -20", "hold 32 -2")
        .replace("hold 256 0", "hold 64 0")
        .replace("coupling 64 9216", "coupling 256 442368")
        .replace(SECOND_ESCAPE, &format!("{SECOND_ESCAPE}reservoir_feedback aperture 1 -20 20 16 16\n"));
    let mut renderer = read(&feedback, 1).into_renderer();
    assert_eq!(renderer.reservoir_observation_count(), 2);
    assert_eq!(renderer.system().reservoir_radiation_aperture(), Some(1));
    renderer.block(&mut [0.0]).unwrap();
    assert_eq!(renderer.last_frame().unwrap().reservoir_exchanges.len(), 2);
    assert!(renderer.last_frame().unwrap().reservoir_radiation.is_some());
}

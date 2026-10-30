//! G0/G3: warm breath, reciprocal escape and another observed escape coexist.
//! The published numerical fixture proves ownership, not instrument fidelity.
use super::*;

const ESCAPE: &str = "reservoir_observation aperture 1 baffled -0.08 0 0\n";
const EXTRA: &str = "reservoir_observation aperture 2 baffled 0.06 0 0\n";
const FEEDBACK: &str = "reservoir_feedback aperture 1 -20 20 16 16\n";

fn input() -> String {
    include_str!("../../../../../examples/resolved-dg-separate-openings.performance")
        .replace("reservoir_apertures 2", "reservoir_apertures 3")
        .replace("initial rest\n",
            "reservoir_aperture 1 high reference 1\nreservoir_cell 10 2 1\ninitial rest\n")
        // Observation order is deliberately different from fluid owner order;
        // the first selected surface is not the reciprocal radiation owner.
        .replace(ESCAPE, &format!("{EXTRA}{ESCAPE}"))
        .replace("reservoir_feedback -20 20 16 16\n", FEEDBACK)
}

#[test]
fn g3_native_extra_escape_preserves_warm_breath_feedback_and_complete_pressure_replay() {
    let text = input();
    let mut observed = read(&text, 8).into_renderer();
    let mut single = read(&text.replace(EXTRA, ""), 8).into_renderer();
    assert_eq!(observed.system().zone().reservoir_aperture_count(), 3);
    assert_eq!(observed.reservoir_observation_count(), 2);
    assert_eq!(observed.system().reservoir_radiation_aperture(), Some(1));
    assert!(observed.reservoir_receiver_at(0).is_none());
    assert!(observed.reservoir_observation_aperture().is_none());
    let mut bore = observed.baffled_receiver().unwrap().clone();
    let mut extra = observed.reservoir_receiver_at(2).unwrap().clone();
    let mut escape = observed.reservoir_receiver_at(1).unwrap().clone();
    assert_eq!(extra.input_count(), 1);
    assert_eq!(escape.input_count(), 1);
    for aperture in [1, 2] {
        let gas = observed.reservoir_observer_gas_at(aperture).unwrap();
        assert_eq!(gas.gas, observed.observer_gas().unwrap().gas);
        assert_eq!(gas.additional_latency_samples,
            observed.observer_gas().unwrap().additional_latency_samples);
    }
    let mut pressure = vec![0.0; 64];
    let mut extra_peak = 0.0_f64;
    for (index, output) in pressure.iter_mut().enumerate() {
        observed.block(core::slice::from_mut(output)).unwrap();
        single.block(&mut [0.0]).unwrap();
        let frame = observed.last_frame().unwrap();
        assert_eq!(frame.reservoir_exchanges.len(), 3);
        assert_eq!(frame.reservoir_flows.len(), 3);
        assert!(matches!(frame.drive,
            NsfDgDrive::ReservoirAperture { aperture: 0, .. }));
        let reaction = frame.reservoir_radiation.as_ref().unwrap();
        assert_eq!(reaction.resolved_exchange, frame.reservoir_exchanges[1]);
        assert_eq!(reaction.resolved_exchange.aperture, 1);
        let temperatures = observed.system().zone();
        assert_eq!(temperatures.reservoir_aperture_at(0).unwrap().state.temperature_k,
            if index < 32 { 300.15 } else { 300.0 });
        for aperture in [1, 2] {
            assert_eq!(temperatures.reservoir_aperture_at(aperture).unwrap().state.temperature_k, 300.0);
        }
        let q = |aperture| frame.reservoir_flows.iter()
            .filter(|flow| flow.aperture == aperture)
            .map(|flow| flow.outward_volume_flow_m3_s).collect::<Vec<_>>();
        let terminals: Vec<_> = frame.terminals.iter().map(|port| port.outward_flow_m3_s).collect();
        let bore_pa = bore.step_trace(&terminals).unwrap();
        let extra_pa = extra.step_trace(&q(2)).unwrap();
        let escape_pa = escape.step_trace(&q(1)).unwrap();
        assert_eq!(*output, (bore_pa + extra_pa) + escape_pa);
        extra_peak = extra_peak.max(extra_pa.abs());
        // Adding a microphone path must not alter either physical acoustic
        // load, the independent heated breath, or the existing escape reaction.
        assert_eq!(observed.system().zone().state_hash(), single.system().zone().state_hash());
        assert_eq!(observed.system().load().state(), single.system().load().state());
        assert_eq!(observed.system().upstream_load().unwrap().state(),
            single.system().upstream_load().unwrap().state());
        assert_eq!(frame.reservoir_exchanges, single.last_frame().unwrap().reservoir_exchanges);
        assert_eq!(frame.reservoir_radiation, single.last_frame().unwrap().reservoir_radiation);
    }
    assert!(extra_peak > 1e-12, "the additional original escape must reach pressure");
    let mut replay = read(&text, 3).into_renderer();
    let mut partitioned = vec![0.0; pressure.len()];
    for chunk in partitioned.chunks_mut(3) { replay.block(chunk).unwrap(); }
    assert_eq!(pressure, partitioned);
    assert_eq!(observed.system().zone().state_hash(), replay.system().zone().state_hash());
    assert_eq!(observed.system().load().state(), replay.system().load().state());
    assert_eq!(observed.system().upstream_load().unwrap().state(),
        replay.system().upstream_load().unwrap().state());
    assert_eq!(observed.system().accepted_steps(), 64);
}

#[test]
fn g0_native_multiple_observers_require_a_distinct_explicit_feedback_owner() {
    let text = input();
    for invalid in [
        text.replace(FEEDBACK, "reservoir_feedback -20 20 16 16\n"),
        text.replace(FEEDBACK, "reservoir_feedback aperture 0 -20 20 16 16\n"),
        text.replace("boundary aperture 0 3", "boundary aperture 1 3"),
        text.replace(EXTRA, ESCAPE),
        text.replace(ESCAPE, ""),
        text.replace(ESCAPE, "reservoir_observation aperture 1 baffled -0.08 0 0.2\n"),
    ] {
        assert!(matches!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free())), Err(ResolvedDgInputError::Input { .. })));
    }
}

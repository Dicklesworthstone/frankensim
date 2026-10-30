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
            "300 1\nreservoir_cell 5 2 1",
            "300 3\nreservoir_cell 5 2 1\nreservoir_cell 3 2 1\nreservoir_cell 4 2 1",
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

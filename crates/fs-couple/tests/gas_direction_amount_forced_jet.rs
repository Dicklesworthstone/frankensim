//! The flue gate of `gas_direction_amount.rs` once more under the opt-in forced jet (ruling 7518):
//! the same shipped examples/flue.performance, read by the same entry point, with the production
//! row `jet_response forced-bickley 0.1 90` after its `jet` row. The owner's test is unchanged and
//! grades the default chain. This twin fails today and stays ignored until it passes unchanged;
//! the default switches only then (bead pexfp.135.1).
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;

/// Autocorrelation f0 over a steady window, 60 Hz to `highest` Hz. A periodic
/// signal also correlates at two and three periods, so the period is the
/// shortest lag whose normalized score is within 5 % of the best (stated rule).
fn f0(signal: &[f64], rate: f64, highest: f64) -> f64 {
    let lags = (rate / highest) as usize..=(rate / 60.0) as usize;
    let score = |lag: usize| {
        let (a, b) = (&signal[..signal.len() - lag], &signal[lag..]);
        a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>()
            / (a.iter().map(|x| x * x).sum::<f64>() * b.iter().map(|y| y * y).sum::<f64>()).sqrt()
    };
    let scores: Vec<(usize, f64)> = lags.map(|lag| (lag, score(lag))).collect();
    let best = scores.iter().map(|(_, s)| *s).fold(f64::NEG_INFINITY, f64::max);
    // The first local maximum within 5 % of the best, not the best itself.
    let period = scores.windows(3).find(|w| w[1].1 >= w[0].1 && w[1].1 >= w[2].1 && w[1].1 >= 0.95 * best)
        .map_or(scores.iter().max_by(|p, q| p.1.total_cmp(&q.1)).unwrap().0, |w| w[1].0);
    rate / period as f64
}

const FLUE: &str = include_str!("../examples/flue.performance");
const FLUE_JET: &str = "jet 90 1e-6 3 0\n";
const FORCED_JET: &str = "jet 90 1e-6 3 0\njet_response forced-bickley 0.1 90\n";
const FLUE_AIR: &str = "gas 293.15 101325 0\n";
const FLUE_HELIOX: &str = "gas mixture 293.15 101325 mole 2 He 0.8 O2 0.2\n";
const FLUE_CO2: &str = "gas mixture 293.15 101325 mole 1 CO2 1\n";

fn flue(gas: &str) -> FlueFilePerformance {
    assert_eq!(FLUE.matches(FLUE_AIR).count(), 1);
    assert_eq!(FLUE.matches(FLUE_JET).count(), 1);
    let text = FLUE.replace(FLUE_AIR, gas).replace(FLUE_JET, FORCED_JET);
    let performance = FlueFilePerformance::from_bytes(text.as_bytes(), 512, &CancelGate::new()).unwrap();
    assert!(performance.renderer().voice().forced_jet().is_some(), "the twin renders the forced jet");
    performance
}

fn render_flue(record: &str) -> Vec<f64> {
    let mut renderer = flue(record).into_renderer();
    let mut pressure = vec![0.0; 24_000];
    for block in pressure.chunks_mut(480) { renderer.block(block).unwrap(); }
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure
}

#[test]
#[ignore = "frankensim-irp-v5-instrument-realism-pexfp.135.1: forced jet lands high: air and heliox speak mode 2 (1143 and 2087 Hz), CO2 mode 3 (1297 Hz), so CO2 sounds above air"]
fn the_forced_jet_flue_sounds_higher_in_heliox_and_lower_in_carbon_dioxide() {
    // 0.5 s at 48 kHz, as shipped; steady window 0.20-0.34 s, before the
    // shipped release at 0.35 s. Direction asserted; amount printed against
    // the c ratio with no pass band (jet-tube pulling, research-grade).
    let rate = 48_000.0;
    let steady = 9_600..16_320;
    let rms = |p: &[f64]| (p.iter().map(|x| x * x).sum::<f64>() / p.len() as f64).sqrt();
    let air = render_flue(FLUE_AIR);
    assert!(rms(&air[steady.clone()]) > 0.0);
    let f0_air = f0(&air[steady.clone()], rate, 2500.0);
    let c_air = flue(FLUE_AIR).renderer().voice().network().unwrap().section_gases()[0].sound_speed;
    for (name, record, up) in [("heliox", FLUE_HELIOX, true), ("co2", FLUE_CO2, false)] {
        let pressure = render_flue(record);
        let level = rms(&pressure[steady.clone()]) / rms(&air[steady.clone()]);
        let f0_gas = f0(&pressure[steady.clone()], rate, 2500.0);
        let c_gas = flue(record).renderer().voice().network().unwrap().section_gases()[0].sound_speed;
        println!("{{\"suite\":\"fs-couple\",\"case\":\"flue-rendered-direction-forced-jet\",\"gas\":\"{name}\",\"f0_air_hz\":{f0_air:.2},\"f0_gas_hz\":{f0_gas:.2},\"f0_ratio\":{:.4},\"c_ratio\":{:.4},\"rms_over_air\":{level:.3e},\"amount\":\"printed, no band\"}}",
            f0_gas / f0_air, c_gas / c_air);
        assert!(level > 1e-3, "{name}: the flue fell silent (rms {level:e} of air's)");
        assert_eq!(f0_gas > f0_air, up, "{name}: sounded {f0_gas} Hz against {f0_air} Hz");
    }
}

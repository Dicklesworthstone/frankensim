use super::*;
use crate::conjugate::AirSegment;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 57, kernel_id: 921, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn paths(capacity_scale: [f64; 2], h_scale: [f64; 4]) -> Vec<AirPath> {
    vec![AirPath::new(293.0, 1.7*capacity_scale[0], 1.0, vec![
        AirSegment::new("upstream", 0.03, 20.0*h_scale[0]).unwrap(),
        AirSegment::new("downstream", 0.04, 35.0*h_scale[1]).unwrap(),
        AirSegment::new("last", 0.02, 12.0*h_scale[2]).unwrap(),
    ]).unwrap(), AirPath::new(312.0, 2.3*capacity_scale[1], 1.0, vec![
        AirSegment::new("independent", 0.05, 17.0*h_scale[3]).unwrap(),
    ]).unwrap()]
}
const NAMES: [&str; 4] = ["upstream", "downstream", "last", "independent"];
const WALLS: [f64; 4] = [340.0, 320.0, 331.0, 299.0];
const BARS: [f64; 4] = [0.2, -0.7, 1.1, 0.4];
const DIRECT: [f64; 4] = [-2.0, 1.3, -0.8, 0.5];
fn value(paths: &[AirPath]) -> f64 {
    let mut at = 0;
    let mut value = 0.0;
    for path in paths {
        let end = at + path.segments().len();
        let march = path.march(&WALLS[at..end]).unwrap();
        for (j, (state, segment)) in march.segments.iter().zip(path.segments()).enumerate() {
            value += BARS[at+j]*state.reference_temperature_k
                + DIRECT[at+j]*segment.htc_w_per_m2_k().ln();
        }
        at = end;
    }
    value
}

#[test]
fn air_transport_controls_match_actual_multisegment_marches() {
    context(&CancelGate::new_clock_free(), |cx| {
        let paths = paths([1.0;2], [1.0;4]);
        let g = pullback_transport_controls(cx, &paths, &NAMES, &WALLS, &BARS, &DIRECT, 4).unwrap();
        let step = 1e-5_f64;
        for p in 0..6 {
            let mut values = Vec::new();
            for sign in [-1.0, 1.0] {
                let mut c = [1.0;2]; let mut h = [1.0;4];
                if p < 2 { c[p] = (sign*step).exp(); }
                else { h[p-2] = (sign*step).exp(); }
                values.push(value(&self::paths(c,h)));
            }
            let fd = (values[1]-values[0])/(2.0*step);
            let actual = if p < 2 { g.log_capacity_rates[p] } else { g.log_htc[p-2] };
            assert!((actual-fd).abs() < 1e-7*actual.abs().max(1.0), "control {p}: {actual} vs {fd}");
        }
        // Scaling C and every h equally preserves every NTU/reference.
        let common: f64 = g.log_capacity_rates.iter().chain(&g.log_htc).sum();
        assert!((common-DIRECT.iter().sum::<f64>()).abs() < 1e-12);
        assert!((g.log_htc[0]-DIRECT[0]).abs() > 0.1,
            "upstream h must include all downstream reference changes");
    });
}

#[test]
fn air_transport_controls_keep_weak_and_saturated_exchange_defined() {
    context(&CancelGate::new_clock_free(), |cx| {
        for z in [1e-14, 1e-8, 1e-3, 1000.0, 1e6] {
            let path = AirPath::new(300.0, 1.0, 1.0,
                vec![AirSegment::new("wall", 1.0, z).unwrap()]).unwrap();
            let g = pullback_transport_controls(cx, &[path], &["wall"], &[350.0], &[1.0], &[0.0], 1).unwrap();
            assert!(g.log_htc[0] > 0.0 && g.log_capacity_rates[0] < 0.0);
            assert_eq!(g.log_htc[0], -g.log_capacity_rates[0]);
            if z < 1e-8 { assert!((g.log_htc[0]/(25.0*z)-1.0).abs() < 1e-8); }
            if z >= 1000.0 { assert!((g.log_htc[0]-50.0/z).abs() < 1e-14); }
        }
    });
}

#[test]
fn air_transport_controls_refuse_bad_bindings_and_cancel_without_partial_results() {
    let paths = paths([1.0;2], [1.0;4]);
    context(&CancelGate::new_clock_free(), |cx| {
        assert!(pullback_transport_controls(cx, &paths, &NAMES, &WALLS, &BARS, &DIRECT, 3).is_err());
        assert!(pullback_transport_controls(cx, &paths, &NAMES[..3], &WALLS, &BARS, &DIRECT, 4).is_err());
        let mut names = NAMES; names.swap(0,1);
        assert!(pullback_transport_controls(cx, &paths, &names, &WALLS, &BARS, &DIRECT, 4).is_err());
        let mut bars = BARS; bars[2] = f64::NAN;
        assert!(pullback_transport_controls(cx, &paths, &NAMES, &WALLS, &bars, &DIRECT, 4).is_err());
        let mut walls = WALLS; walls[3] = 0.0;
        assert!(pullback_transport_controls(cx, &paths, &NAMES, &walls, &BARS, &DIRECT, 4).is_err());
        assert!(pullback_transport_controls(cx, &[], &[], &[], &[], &[], 0).is_err());
    });
    let gate = CancelGate::new_clock_free(); gate.request();
    context(&gate, |cx| assert!(matches!(pullback_transport_controls(cx, &paths, &NAMES,
        &WALLS, &BARS, &DIRECT, 4), Err(super::super::super::CoupledGoalError::Interrupted))));
}

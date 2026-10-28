//! G1/G3: independent scalar adjoints on a nonuniform pulse grid, including
//! zero amplitude and an earlier selected state. Malformed receipts refuse.
use super::*;
use fs_project::{TransientPowerStep, spec::dims};
use fs_qty::QtyAny;

const COARSE: [f64; 3] = [0.25, 0.75, 2.0];
const ENDS: [f64; 6] = [0.125, 0.25, 0.5, 0.75, 1.375, 2.0];
const POWER: [f64; 6] = [8.0, 8.0, 0.0, 0.0, 2.0, 2.0];
fn schedule(amplitude: f64) -> TransientRegionPower {
    TransientRegionPower { region: "chip".into(), source: "independent synthetic pulse".into(),
        steps: [(0.25, 8.0), (0.75, 0.0), (2.0, 2.0)].into_iter().map(|(t, p)|
            TransientPowerStep { until: QtyAny::new(t, dims::TIME),
                watts: QtyAny::new(p * amplitude, dims::POWER) }).collect() }
}
fn gradients(prefix: usize) -> Vec<f64> {
    let mut values = vec![0.0; ENDS.len()];
    let mut carry = 1.0;
    for i in (0..prefix).rev() {
        let dt = ENDS[i] - if i == 0 { 0.0 } else { ENDS[i - 1] };
        values[i] = carry * dt / (10.0 + 2.0 * dt);
        carry *= 10.0 / (10.0 + 2.0 * dt);
    }
    values
}
fn scalar(amplitude: f64, prefix: usize) -> f64 {
    let mut t = 300.0;
    let mut start = 0.0;
    for i in 0..prefix {
        let dt = ENDS[i] - start;
        t = (10.0 * t + dt * (POWER[i] * amplitude + 2.0 * 300.0)) / (10.0 + 2.0 * dt);
        start = ENDS[i];
    }
    t
}
fn text(values: &[f64]) -> String {
    let mut start = 0.0;
    let rows = ENDS.iter().zip(values).enumerate().map(|(i, (&end, &g))| {
        let row = format!(r#"{{"step":{},"start_s":{},"end_s":{},"regional_power":[{{"region":"chip","dtemperature_ddelivered_watts_k_per_w":{}}}]}}"#,
            i + 1, start, end, g);
        start = end;
        row
    }).collect::<Vec<_>>().join(",");
    format!(r#"{{"intervals":[{rows}]}}"#)
}

#[test]
fn g1_original_watts_contract_without_a_second_dt_or_duty_at_zero_and_nonunit_amplitudes() {
    for amplitude in [0.0, 0.5, 1.3, 2.0] {
        for prefix in [0, 2, 6] {
            let g = gradients(prefix);
            let report = J::parse(&text(&g)).unwrap();
            let derivative = contract(&report, &schedule(1.0), &schedule(amplitude), &COARSE, amplitude, prefix).unwrap();
            let reference: f64 = POWER.iter().zip(&g).map(|(p, g)| p * g).sum();
            assert!((derivative - reference).abs() <= 1e-14);
            let epsilon = 0.001;
            let difference = if amplitude == 0.0 {
                (scalar(epsilon, prefix) - scalar(0.0, prefix)) / epsilon
            } else {
                (scalar(amplitude + epsilon, prefix) - scalar(amplitude - epsilon, prefix)) / (2.0 * epsilon)
            };
            assert!((derivative - difference).abs() <= 1e-8, "{amplitude} {prefix}: {derivative} versus {difference}");
            if prefix == 0 { assert_eq!(derivative, 0.0); }
            else { assert!(derivative > 0.0, "a=0 is not a zero physical sensitivity"); }
        }
    }
}

#[test]
fn g3_step_identity_region_and_future_derivatives_are_checked_before_contraction() {
    let source = text(&gradients(2));
    for bad in [
        source.replace("\"step\":2", "\"step\":1"),
        source.replace("\"end_s\":0.25", "\"end_s\":0.24"),
        source.replace("\"start_s\":0.25", "\"start_s\":0.26"),
        source.replace("\"chip\"", "\"other\""),
        source.replace("\"regional_power\":[", "\"regional_power\":[{\"region\":\"chip\",\"dtemperature_ddelivered_watts_k_per_w\":0},"),
        "{\"intervals\":[]}".to_owned(),
        source.replace("\"dtemperature_ddelivered_watts_k_per_w\":0}", "\"dtemperature_ddelivered_watts_k_per_w\":1}"),
    ] {
        let report = J::parse(&bad).unwrap();
        assert!(contract(&report, &schedule(1.0), &schedule(1.3), &COARSE, 1.3, 2).is_err(), "{bad}");
    }
}

#[test]
fn g3_original_source_and_event_alignment_cannot_be_replaced_by_scaled_or_static_power() {
    let report = J::parse(&text(&gradients(6))).unwrap();
    let original = schedule(1.0);
    for mode in 0..5 {
        let mut sampled = schedule(1.3);
        match mode {
            0 => sampled.steps[0].watts.value *= 0.25,
            1 => sampled.steps[0].until.value = 0.2,
            2 => sampled.source = "unrelated source".into(),
            3 => sampled.steps[0].watts.dims = dims::ENERGY,
            _ => { sampled.steps.pop(); }
        }
        assert!(contract(&report, &original, &sampled, &COARSE, 1.3, 6).is_err());
    }
    let mut crossing = original.clone();
    crossing.steps[0].until.value = 0.2;
    assert!(contract(&report, &crossing, &crossing, &COARSE, 1.0, 6).is_err());
    assert!(contract(&report, &original, &original, &COARSE, 1.0, 7).is_err());
    assert!(contract(&report, &original, &original, &[], 1.0, 0).is_err());
}

#[test]
fn g0_nonfinite_and_unrepresentable_amplitude_contractions_refuse() {
    let report = J::parse(&text(&gradients(6))).unwrap();
    for amplitude in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(contract(&report, &schedule(1.0), &schedule(0.0), &COARSE, amplitude, 6).is_err());
    }
    let overflowing = J::parse(&text(&[f64::MAX; 6])).unwrap();
    assert!(contract(&overflowing, &schedule(1.0), &schedule(1.0), &COARSE, 1.0, 6).is_err());
    let mut tiny = schedule(1.0);
    tiny.steps[0].watts.value = f64::from_bits(1);
    let mut erased = tiny.clone();
    for step in &mut erased.steps { step.watts.value *= 0.25; }
    assert!(contract(&report, &tiny, &erased, &COARSE, 0.25, 6).is_err());
}

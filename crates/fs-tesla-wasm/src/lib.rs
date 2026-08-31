use fs_flux::lc::{TeslaCoilParams, step_tesla_coil};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

pub const MAX_RESONANT_FREQ_KHZ: f64 = 10_000.0;
pub const MAX_INPUT_KV: f64 = 1_000.0;
pub const MAX_SPARK_GAP_MM: f64 = 1_000.0;
pub const MAX_Q_FACTOR: f64 = 1_000_000.0;

fn refusal_json(code: &str, message: &str, first_repair: &str, second_repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\
         \"ranked_repairs\":[\"{first_repair}\",\"{second_repair}\"]}}}}"
    )
}

fn admit_inputs(
    resonant_freq_khz: f64,
    input_kv: f64,
    spark_gap_mm: f64,
    q_factor: f64,
) -> Result<(), String> {
    if ![resonant_freq_khz, input_kv, spark_gap_mm, q_factor]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(refusal_json(
            "non-finite-input",
            "every Tesla boundary input must be finite",
            "replace NaN or infinity with a finite value",
            "check host-side unit conversion before calling the kernel",
        ));
    }

    let admitted = (0.0..=MAX_RESONANT_FREQ_KHZ).contains(&resonant_freq_khz)
        && resonant_freq_khz > 0.0
        && (0.0..=MAX_INPUT_KV).contains(&input_kv)
        && input_kv > 0.0
        && (0.0..=MAX_SPARK_GAP_MM).contains(&spark_gap_mm)
        && spark_gap_mm > 0.0
        && (0.0..=MAX_Q_FACTOR).contains(&q_factor)
        && q_factor > 0.0;
    if !admitted {
        return Err(refusal_json(
            "input-outside-domain",
            "Tesla boundary input is outside the admitted positive finite domain",
            "use positive frequency voltage spark gap and quality factor values",
            "split larger studies into values below the documented boundary caps",
        ));
    }
    Ok(())
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn tesla_coil_step(
    resonant_freq_khz: f64,
    input_kv: f64,
    spark_gap_mm: f64,
    q_factor: f64,
) -> String {
    if let Err(refusal) = admit_inputs(resonant_freq_khz, input_kv, spark_gap_mm, q_factor) {
        return refusal;
    }

    let params = TeslaCoilParams {
        resonant_freq_khz,
        input_kv,
        spark_gap_mm,
        q_factor,
    };

    // House canonical emission: the Franken-only dependency policy
    // (Decalogue P1) forbids a serde edge here, and one flat record needs
    // nothing smarter than explicit field formatting.
    let result = step_tesla_coil(&params);
    if ![
        result.resonant_freq_khz,
        result.secondary_potential_mv,
        result.streamer_length_inches,
        result.streamer_length_meters,
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return refusal_json(
            "non-finite-output",
            "Tesla kernel produced a non-finite result",
            "reduce the admitted input magnitudes",
            "inspect the owning fs-flux kernel before retrying",
        );
    }
    if [
        result.resonant_freq_khz,
        result.secondary_potential_mv,
        result.streamer_length_inches,
        result.streamer_length_meters,
    ]
    .iter()
    .any(|value| *value <= 0.0)
    {
        return refusal_json(
            "output-outside-domain",
            "Tesla kernel produced a non-positive result",
            "increase inputs above floating-point underflow scale",
            "inspect the owning fs-flux kernel before retrying",
        );
    }
    format!(
        "{{\"ok\":{{\"resonant_freq_khz\":{},\"secondary_potential_mv\":{},\
         \"streamer_length_inches\":{},\"streamer_length_meters\":{}}}}}",
        result.resonant_freq_khz,
        result.secondary_potential_mv,
        result.streamer_length_inches,
        result.streamer_length_meters,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_step_returns_ok_envelope() {
        let result = tesla_coil_step(100.0, 15.0, 12.0, 145.0);
        assert!(result.starts_with("{\"ok\":"), "{result}");
        assert!(result.contains("\"secondary_potential_mv\":"), "{result}");
    }

    #[test]
    fn non_finite_input_is_a_typed_refusal() {
        let result = tesla_coil_step(100.0, f64::INFINITY, 12.0, 145.0);
        assert!(result.contains("\"code\":\"non-finite-input\""), "{result}");
        assert!(result.contains("\"ranked_repairs\""), "{result}");
    }

    #[test]
    fn non_positive_inputs_are_typed_refusals() {
        for result in [
            tesla_coil_step(0.0, 15.0, 12.0, 145.0),
            tesla_coil_step(100.0, 0.0, 12.0, 145.0),
            tesla_coil_step(100.0, 15.0, 0.0, 145.0),
            tesla_coil_step(100.0, 15.0, 12.0, 0.0),
        ] {
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }

    #[test]
    fn underflowed_output_refuses_instead_of_returning_zero_as_ok() {
        let result = tesla_coil_step(100.0, f64::from_bits(1), 12.0, 145.0);
        assert!(
            result.contains("\"code\":\"output-outside-domain\""),
            "{result}"
        );
    }

    #[test]
    fn admission_caps_hold_at_cap_and_refuse_cap_plus_one() {
        assert!(
            tesla_coil_step(
                MAX_RESONANT_FREQ_KHZ,
                MAX_INPUT_KV,
                MAX_SPARK_GAP_MM,
                MAX_Q_FACTOR,
            )
            .starts_with("{\"ok\":")
        );
        for result in [
            tesla_coil_step(
                MAX_RESONANT_FREQ_KHZ + 1.0,
                MAX_INPUT_KV,
                MAX_SPARK_GAP_MM,
                MAX_Q_FACTOR,
            ),
            tesla_coil_step(
                MAX_RESONANT_FREQ_KHZ,
                MAX_INPUT_KV + 1.0,
                MAX_SPARK_GAP_MM,
                MAX_Q_FACTOR,
            ),
            tesla_coil_step(
                MAX_RESONANT_FREQ_KHZ,
                MAX_INPUT_KV,
                MAX_SPARK_GAP_MM + 1.0,
                MAX_Q_FACTOR,
            ),
            tesla_coil_step(
                MAX_RESONANT_FREQ_KHZ,
                MAX_INPUT_KV,
                MAX_SPARK_GAP_MM,
                MAX_Q_FACTOR + 1.0,
            ),
        ] {
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }
}

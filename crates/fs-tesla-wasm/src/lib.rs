use fs_flux::lc::{QuarterWaveParams, step_quarter_wave};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

pub const MAX_FREQUENCY_HZ: f64 = 1_000_000_000.0;
pub const MAX_PROPAGATION_SPEED_MPS: f64 = 400_000_000.0;
pub const MAX_CONDUCTOR_LENGTH_M: f64 = 1_000_000_000.0;

fn refusal_json(code: &str, message: &str, first_repair: &str, second_repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\
         \"ranked_repairs\":[\"{first_repair}\",\"{second_repair}\"]}}}}"
    )
}

fn admit_inputs(
    frequency_hz: f64,
    propagation_speed_mps: f64,
    conductor_length_m: f64,
) -> Result<(), String> {
    if ![frequency_hz, propagation_speed_mps, conductor_length_m]
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

    let admitted = (0.0..=MAX_FREQUENCY_HZ).contains(&frequency_hz)
        && frequency_hz > 0.0
        && (0.0..=MAX_PROPAGATION_SPEED_MPS).contains(&propagation_speed_mps)
        && propagation_speed_mps > 0.0
        && (0.0..=MAX_CONDUCTOR_LENGTH_M).contains(&conductor_length_m)
        && conductor_length_m > 0.0;
    if !admitted {
        return Err(refusal_json(
            "input-outside-domain",
            "Tesla boundary input is outside the admitted positive finite domain",
            "use positive frequency propagation-speed and conductor-length values",
            "split larger studies into values below the documented boundary caps",
        ));
    }
    Ok(())
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn tesla_transformer_step(
    frequency_hz: f64,
    propagation_speed_mps: f64,
    conductor_length_m: f64,
) -> String {
    if let Err(refusal) = admit_inputs(frequency_hz, propagation_speed_mps, conductor_length_m) {
        return refusal;
    }

    let params = QuarterWaveParams {
        frequency_hz,
        propagation_speed_mps,
        conductor_length_m,
    };

    // House canonical emission: the Franken-only dependency policy
    // (Decalogue P1) forbids a serde edge here, and one flat record needs
    // nothing smarter than explicit field formatting.
    let result = step_quarter_wave(&params);
    if ![
        result.wavelength_m,
        result.quarter_wave_length_m,
        result.electrical_length_rad,
        result.quarter_wave_error_rad,
        result.length_error_m,
        result.length_ratio,
        result.remote_terminal_profile_fraction,
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
    if result.wavelength_m <= 0.0
        || result.quarter_wave_length_m <= 0.0
        || result.electrical_length_rad <= 0.0
        || result.length_ratio <= 0.0
        || !(0.0..=1.0).contains(&result.remote_terminal_profile_fraction)
    {
        return refusal_json(
            "output-outside-domain",
            "Tesla kernel produced a non-positive result",
            "increase inputs above floating-point underflow scale",
            "inspect the owning fs-flux kernel before retrying",
        );
    }
    format!(
        "{{\"ok\":{{\"wavelength_m\":{},\"quarter_wave_length_m\":{},\
         \"electrical_length_rad\":{},\"quarter_wave_error_rad\":{},\
         \"length_error_m\":{},\"length_ratio\":{},\
         \"remote_terminal_profile_fraction\":{}}}}}",
        result.wavelength_m,
        result.quarter_wave_length_m,
        result.electrical_length_rad,
        result.quarter_wave_error_rad,
        result.length_error_m,
        result.length_ratio,
        result.remote_terminal_profile_fraction,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_step_returns_ok_envelope() {
        let result = tesla_transformer_step(925.0, 185_000.0 * 1_609.344, 50.0 * 1_609.344);
        assert!(result.starts_with("{\"ok\":"), "{result}");
        assert!(result.contains("\"quarter_wave_length_m\":"), "{result}");
    }

    #[test]
    fn non_finite_input_is_a_typed_refusal() {
        let result = tesla_transformer_step(925.0, f64::INFINITY, 50.0 * 1_609.344);
        assert!(result.contains("\"code\":\"non-finite-input\""), "{result}");
        assert!(result.contains("\"ranked_repairs\""), "{result}");
    }

    #[test]
    fn non_positive_inputs_are_typed_refusals() {
        for result in [
            tesla_transformer_step(0.0, 185_000.0 * 1_609.344, 50.0 * 1_609.344),
            tesla_transformer_step(925.0, 0.0, 50.0 * 1_609.344),
            tesla_transformer_step(925.0, 185_000.0 * 1_609.344, 0.0),
        ] {
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }

    #[test]
    fn underflowed_output_refuses_instead_of_returning_zero_as_ok() {
        let result = tesla_transformer_step(f64::from_bits(1), 185_000.0, f64::from_bits(1));
        assert!(
            result.contains("\"code\":\"non-finite-output\""),
            "{result}"
        );
    }

    #[test]
    fn admission_caps_hold_at_cap_and_refuse_cap_plus_one() {
        assert!(
            tesla_transformer_step(
                MAX_FREQUENCY_HZ,
                MAX_PROPAGATION_SPEED_MPS,
                MAX_CONDUCTOR_LENGTH_M,
            )
            .starts_with("{\"ok\":")
        );
        for result in [
            tesla_transformer_step(
                MAX_FREQUENCY_HZ + 1.0,
                MAX_PROPAGATION_SPEED_MPS,
                MAX_CONDUCTOR_LENGTH_M,
            ),
            tesla_transformer_step(
                MAX_FREQUENCY_HZ,
                MAX_PROPAGATION_SPEED_MPS + 1.0,
                MAX_CONDUCTOR_LENGTH_M,
            ),
            tesla_transformer_step(
                MAX_FREQUENCY_HZ,
                MAX_PROPAGATION_SPEED_MPS,
                MAX_CONDUCTOR_LENGTH_M + 1.0,
            ),
        ] {
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }
}

use fs_mbd::goddard::{GoddardParams, step_goddard_rocket};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

pub const MAX_CHAMBER_PRESSURE_PSI: f64 = 10_000.0;
pub const MAX_FUEL_FLOW_KG_PER_SEC: f64 = 1_000.0;
pub const MAX_THROAT_AREA_CM2: f64 = 100_000.0;
pub const MIN_EXPANSION_RATIO: f64 = 1.4;
pub const MAX_EXPANSION_RATIO: f64 = 1_000.0;

fn refusal_json(code: &str, message: &str, first_repair: &str, second_repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\
         \"ranked_repairs\":[\"{first_repair}\",\"{second_repair}\"]}}}}"
    )
}

fn admit_inputs(
    chamber_pressure_psi: f64,
    fuel_flow_kg_per_sec: f64,
    throat_area_cm2: f64,
    expansion_ratio: f64,
) -> Result<(), String> {
    if ![
        chamber_pressure_psi,
        fuel_flow_kg_per_sec,
        throat_area_cm2,
        expansion_ratio,
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return Err(refusal_json(
            "non-finite-input",
            "every Goddard boundary input must be finite",
            "replace NaN or infinity with a finite value",
            "check host-side unit conversion before calling the kernel",
        ));
    }

    let admitted = (0.0..=MAX_CHAMBER_PRESSURE_PSI).contains(&chamber_pressure_psi)
        && chamber_pressure_psi > 0.0
        && (0.0..=MAX_FUEL_FLOW_KG_PER_SEC).contains(&fuel_flow_kg_per_sec)
        && fuel_flow_kg_per_sec > 0.0
        && (0.0..=MAX_THROAT_AREA_CM2).contains(&throat_area_cm2)
        && throat_area_cm2 > 0.0
        && (MIN_EXPANSION_RATIO..=MAX_EXPANSION_RATIO).contains(&expansion_ratio);
    if !admitted {
        return Err(refusal_json(
            "input-outside-domain",
            "Goddard boundary input is outside the admitted positive finite domain",
            "use positive pressure flow and throat area with expansion ratio at least 1.4",
            "split larger studies into values below the documented boundary caps",
        ));
    }
    Ok(())
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn goddard_rocket_step(
    chamber_pressure_psi: f64,
    fuel_flow_kg_per_sec: f64,
    throat_area_cm2: f64,
    expansion_ratio: f64,
) -> String {
    if let Err(refusal) = admit_inputs(
        chamber_pressure_psi,
        fuel_flow_kg_per_sec,
        throat_area_cm2,
        expansion_ratio,
    ) {
        return refusal;
    }

    let params = GoddardParams {
        chamber_pressure_psi,
        fuel_flow_kg_per_sec,
        throat_area_cm2,
        expansion_ratio,
    };

    // House canonical emission: the Franken-only dependency policy
    // (Decalogue P1) forbids a serde edge here, and one flat record needs
    // nothing smarter than explicit field formatting.
    let result = step_goddard_rocket(&params);
    if ![
        result.chamber_pressure_psi,
        result.chamber_pressure_pa,
        result.exhaust_velocity_mps,
        result.thrust_newtons,
        result.specific_impulse_sec,
        result.mach_exit,
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return refusal_json(
            "non-finite-output",
            "Goddard kernel produced a non-finite result",
            "reduce the admitted input magnitudes",
            "inspect the owning fs-mbd kernel before retrying",
        );
    }
    if [
        result.chamber_pressure_psi,
        result.chamber_pressure_pa,
        result.exhaust_velocity_mps,
        result.thrust_newtons,
        result.specific_impulse_sec,
        result.mach_exit,
    ]
    .iter()
    .any(|value| *value <= 0.0)
    {
        return refusal_json(
            "output-outside-domain",
            "Goddard kernel produced a non-positive result",
            "increase inputs above floating-point underflow scale",
            "inspect the owning fs-mbd kernel before retrying",
        );
    }
    format!(
        "{{\"ok\":{{\"chamber_pressure_psi\":{},\"chamber_pressure_pa\":{},\
         \"exhaust_velocity_mps\":{},\"thrust_newtons\":{},\
         \"specific_impulse_sec\":{},\"mach_exit\":{}}}}}",
        result.chamber_pressure_psi,
        result.chamber_pressure_pa,
        result.exhaust_velocity_mps,
        result.thrust_newtons,
        result.specific_impulse_sec,
        result.mach_exit,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_step_returns_ok_envelope() {
        let result = goddard_rocket_step(350.0, 1.8, 4.2, 3.5);
        assert!(result.starts_with("{\"ok\":"), "{result}");
        assert!(result.contains("\"thrust_newtons\":"), "{result}");
    }

    #[test]
    fn non_finite_input_is_a_typed_refusal() {
        let result = goddard_rocket_step(f64::NAN, 1.8, 4.2, 3.5);
        assert!(result.contains("\"code\":\"non-finite-input\""), "{result}");
        assert!(result.contains("\"ranked_repairs\""), "{result}");
    }

    #[test]
    fn non_positive_and_below_minimum_inputs_are_typed_refusals() {
        for result in [
            goddard_rocket_step(0.0, 1.8, 4.2, 3.5),
            goddard_rocket_step(350.0, 0.0, 4.2, 3.5),
            goddard_rocket_step(350.0, 1.8, 0.0, 3.5),
            goddard_rocket_step(350.0, 1.8, 4.2, MIN_EXPANSION_RATIO - 0.01),
        ] {
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }

    #[test]
    fn underflowed_output_refuses_instead_of_returning_zero_as_ok() {
        let result = goddard_rocket_step(350.0, f64::from_bits(1), 4.2, 3.5);
        assert!(
            result.contains("\"code\":\"output-outside-domain\""),
            "{result}"
        );
    }

    #[test]
    fn admission_caps_hold_at_cap_and_refuse_cap_plus_one() {
        assert!(
            goddard_rocket_step(
                MAX_CHAMBER_PRESSURE_PSI,
                MAX_FUEL_FLOW_KG_PER_SEC,
                MAX_THROAT_AREA_CM2,
                MAX_EXPANSION_RATIO,
            )
            .starts_with("{\"ok\":")
        );
        for result in [
            goddard_rocket_step(
                MAX_CHAMBER_PRESSURE_PSI + 1.0,
                MAX_FUEL_FLOW_KG_PER_SEC,
                MAX_THROAT_AREA_CM2,
                MAX_EXPANSION_RATIO,
            ),
            goddard_rocket_step(
                MAX_CHAMBER_PRESSURE_PSI,
                MAX_FUEL_FLOW_KG_PER_SEC + 1.0,
                MAX_THROAT_AREA_CM2,
                MAX_EXPANSION_RATIO,
            ),
            goddard_rocket_step(
                MAX_CHAMBER_PRESSURE_PSI,
                MAX_FUEL_FLOW_KG_PER_SEC,
                MAX_THROAT_AREA_CM2 + 1.0,
                MAX_EXPANSION_RATIO,
            ),
            goddard_rocket_step(
                MAX_CHAMBER_PRESSURE_PSI,
                MAX_FUEL_FLOW_KG_PER_SEC,
                MAX_THROAT_AREA_CM2,
                MAX_EXPANSION_RATIO + 1.0,
            ),
        ] {
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }
}

use fs_mbd::goddard::{
    GODDARD_MAX_ELAPSED_SECONDS, GODDARD_MAX_EXPANSION_RATIO, GODDARD_MAX_GYRO_SPIN_RPM,
    GODDARD_MAX_PRIMARY_SPIN_RPM, GODDARD_MAX_TUBE_LENGTH_RATIO, GODDARD_MIN_TUBE_LENGTH_RATIO,
    GoddardApparatusParams, GoddardParams, GoddardRocketError, step_goddard_apparatus,
    step_goddard_rocket,
};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

pub const MAX_CHAMBER_PRESSURE_PSI: f64 = 10_000.0;
pub const MAX_CHAMBER_TEMPERATURE_K: f64 = 10_000.0;
pub const MAX_SPECIFIC_GAS_CONSTANT_J_KG_K: f64 = 10_000.0;
pub const MAX_AMBIENT_PRESSURE_PA: f64 = 100_000_000.0;
pub const MAX_THROAT_AREA_CM2: f64 = 100_000.0;
pub const MIN_EXPANSION_RATIO: f64 = 1.0;
pub const MAX_EXPANSION_RATIO: f64 = GODDARD_MAX_EXPANSION_RATIO;

fn refusal_json(code: &str, message: &str, first_repair: &str, second_repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\
         \"ranked_repairs\":[\"{first_repair}\",\"{second_repair}\"]}}}}"
    )
}

fn admit_inputs(
    chamber_pressure_psi: f64,
    chamber_temperature_k: f64,
    specific_gas_constant_j_kg_k: f64,
    heat_capacity_ratio: f64,
    throat_area_cm2: f64,
    expansion_ratio: f64,
    ambient_pressure_pa: f64,
) -> Result<(), String> {
    if ![
        chamber_pressure_psi,
        chamber_temperature_k,
        specific_gas_constant_j_kg_k,
        heat_capacity_ratio,
        throat_area_cm2,
        expansion_ratio,
        ambient_pressure_pa,
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
        && (0.0..=MAX_CHAMBER_TEMPERATURE_K).contains(&chamber_temperature_k)
        && chamber_temperature_k > 0.0
        && (0.0..=MAX_SPECIFIC_GAS_CONSTANT_J_KG_K).contains(&specific_gas_constant_j_kg_k)
        && specific_gas_constant_j_kg_k > 0.0
        && heat_capacity_ratio > 1.0
        && heat_capacity_ratio <= 5.0 / 3.0
        && (0.0..=MAX_THROAT_AREA_CM2).contains(&throat_area_cm2)
        && throat_area_cm2 > 0.0
        && (MIN_EXPANSION_RATIO..=MAX_EXPANSION_RATIO).contains(&expansion_ratio)
        && (0.0..=MAX_AMBIENT_PRESSURE_PA).contains(&ambient_pressure_pa);
    if !admitted {
        return Err(refusal_json(
            "input-outside-domain",
            "ideal-nozzle input is outside the bounded calorically-perfect domain",
            "declare positive pressure temperature gas constant and throat area; use gamma in (1,5/3] and area ratio at least one",
            "split larger studies into values below the documented boundary caps",
        ));
    }
    Ok(())
}

fn admit_apparatus_inputs(
    elapsed_seconds: f64,
    primary_spin_rpm: f64,
    gyro_spin_rpm: f64,
    tube_length_ratio: f64,
    auxiliary_release_fraction: f64,
) -> Result<(), String> {
    if ![
        elapsed_seconds,
        primary_spin_rpm,
        gyro_spin_rpm,
        tube_length_ratio,
        auxiliary_release_fraction,
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return Err(refusal_json(
            "non-finite-input",
            "every Goddard apparatus input must be finite",
            "replace NaN or infinity with a finite value",
            "check host-side unit conversion before calling the kernel",
        ));
    }

    let admitted = (0.0..=GODDARD_MAX_ELAPSED_SECONDS).contains(&elapsed_seconds)
        && (0.0..=GODDARD_MAX_PRIMARY_SPIN_RPM).contains(&primary_spin_rpm)
        && (0.0..=GODDARD_MAX_GYRO_SPIN_RPM).contains(&gyro_spin_rpm)
        && (GODDARD_MIN_TUBE_LENGTH_RATIO..=GODDARD_MAX_TUBE_LENGTH_RATIO)
            .contains(&tube_length_ratio)
        && (0.0..=1.0).contains(&auxiliary_release_fraction);
    if !admitted {
        return Err(refusal_json(
            "input-outside-domain",
            "Goddard apparatus input is outside the bounded teaching domain",
            "use nonnegative time and spin speeds with release fraction between zero and one",
            "keep the tapered-tube ratio between one and twelve",
        ));
    }
    Ok(())
}

/// Source-bounded US 1,102,653 apparatus step composed from `fs-mbd`.
///
/// The two speed inputs are declared visitor controls because the facsimile
/// prints no numerical spin rates. The output deliberately contains no liquid
/// propellant, de Laval, Mach, thrust, or trajectory field.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn goddard_apparatus_step(
    elapsed_seconds: f64,
    primary_spin_rpm: f64,
    gyro_spin_rpm: f64,
    tube_length_ratio: f64,
    auxiliary_release_fraction: f64,
    primary_charge_substantially_consumed: bool,
    gyro_enabled: bool,
) -> String {
    if let Err(refusal) = admit_apparatus_inputs(
        elapsed_seconds,
        primary_spin_rpm,
        gyro_spin_rpm,
        tube_length_ratio,
        auxiliary_release_fraction,
    ) {
        return refusal;
    }

    let result = match step_goddard_apparatus(&GoddardApparatusParams {
        elapsed_seconds,
        primary_spin_rpm,
        gyro_spin_rpm,
        tube_length_ratio,
        auxiliary_release_fraction,
        primary_charge_substantially_consumed,
        gyro_enabled,
    }) {
        Ok(result) => result,
        Err(_) => {
            return refusal_json(
                "rigid-body-refusal",
                "the generic fs-mbd owner refused the normalized apparatus pose",
                "reduce elapsed time or the declared spin speeds",
                "inspect the owning fs-mbd Goddard apparatus kernel before retrying",
            );
        }
    };

    let scalar_outputs = [
        result.primary_angular_velocity_rad_per_sec,
        result.gyro_angular_velocity_rad_per_sec,
        result.camera_support_angular_velocity_rad_per_sec,
        result.primary_rim_speed_per_radius_mps_per_m,
        result.tube_length_ratio,
        result.claim_2_ratio_margin,
    ];
    if !scalar_outputs
        .iter()
        .chain(result.primary_quaternion.iter())
        .chain(result.gyro_quaternion.iter())
        .all(|value| value.is_finite())
    {
        return refusal_json(
            "non-finite-output",
            "Goddard apparatus kernel produced a non-finite result",
            "reduce elapsed time or the declared spin speeds",
            "inspect the owning fs-mbd Goddard apparatus kernel before retrying",
        );
    }

    format!(
        "{{\"ok\":{{\"primary_quaternion\":[{},{},{},{}],\
         \"gyro_quaternion\":[{},{},{},{}],\
         \"primary_angular_velocity_rad_per_sec\":{},\
         \"gyro_angular_velocity_rad_per_sec\":{},\
         \"camera_support_angular_velocity_rad_per_sec\":{},\
         \"primary_rim_speed_per_radius_mps_per_m\":{},\
         \"tube_length_ratio\":{},\"claim_2_ratio_margin\":{},\
         \"claim_2_satisfied\":{},\"claim_1_sequence_satisfied\":{},\
         \"auxiliary_nested\":{},\"gyro_enabled\":{}}}}}",
        result.primary_quaternion[0],
        result.primary_quaternion[1],
        result.primary_quaternion[2],
        result.primary_quaternion[3],
        result.gyro_quaternion[0],
        result.gyro_quaternion[1],
        result.gyro_quaternion[2],
        result.gyro_quaternion[3],
        result.primary_angular_velocity_rad_per_sec,
        result.gyro_angular_velocity_rad_per_sec,
        result.camera_support_angular_velocity_rad_per_sec,
        result.primary_rim_speed_per_radius_mps_per_m,
        result.tube_length_ratio,
        result.claim_2_ratio_margin,
        result.claim_2_satisfied,
        result.claim_1_sequence_satisfied,
        result.auxiliary_nested,
        result.gyro_enabled,
    )
}

/// Ideal calorically-perfect nozzle with explicit gas properties and ambient.
///
/// Mass flow is computed from the sonic throat. No propellant temperature or
/// heat-capacity ratio is inferred from chamber pressure. A four-argument
/// invocation of the former interface is invalid; every input is now explicit.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn goddard_rocket_step(
    chamber_pressure_psi: f64,
    chamber_temperature_k: f64,
    specific_gas_constant_j_kg_k: f64,
    heat_capacity_ratio: f64,
    throat_area_cm2: f64,
    expansion_ratio: f64,
    ambient_pressure_pa: f64,
) -> String {
    if let Err(refusal) = admit_inputs(
        chamber_pressure_psi,
        chamber_temperature_k,
        specific_gas_constant_j_kg_k,
        heat_capacity_ratio,
        throat_area_cm2,
        expansion_ratio,
        ambient_pressure_pa,
    ) {
        return refusal;
    }

    let params = GoddardParams {
        chamber_pressure_psi,
        chamber_temperature_k,
        specific_gas_constant_j_kg_k,
        heat_capacity_ratio,
        throat_area_cm2,
        expansion_ratio,
        ambient_pressure_pa,
    };

    // House canonical emission: the Franken-only dependency policy
    // (Decalogue P1) forbids a serde edge here, and one flat record needs
    // nothing smarter than explicit field formatting.
    let result = match step_goddard_rocket(&params) {
        Ok(result) => result,
        Err(GoddardRocketError::UnsupportedBackPressure { .. }) => {
            return refusal_json(
                "unsupported-back-pressure",
                "ambient pressure exceeds the isentropic exit pressure; this model does not resolve shocks or separation",
                "use a smaller area ratio or an explicitly lower ambient pressure",
                "use a shock and boundary-layer model for overexpanded operation",
            );
        }
        Err(GoddardRocketError::UnrepresentableOutput(_)) => {
            return refusal_json(
                "output-outside-domain",
                "a positive finite ideal-nozzle quantity is not representable",
                "keep pressure temperature gas constant and area away from floating-point underflow and overflow",
                "check the declared SI units and the pressure and area conversions",
            );
        }
        Err(_) => {
            return refusal_json(
                "nozzle-refusal",
                "the owning fs-mbd ideal-nozzle kernel refused this state",
                "check the explicit gas model and geometry against its documented domain",
                "inspect the owning fs-mbd Goddard rocket kernel before retrying",
            );
        }
    };
    if ![
        result.chamber_pressure_psi,
        result.chamber_pressure_pa,
        result.exit_temperature_k,
        result.exit_pressure_pa,
        result.mass_flow_kg_per_sec,
        result.exhaust_velocity_mps,
        result.momentum_thrust_newtons,
        result.pressure_thrust_newtons,
        result.thrust_newtons,
        result.specific_impulse_sec,
        result.characteristic_velocity_mps,
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
    format!(
        "{{\"ok\":{{\"chamber_pressure_psi\":{},\"chamber_pressure_pa\":{},\
         \"exit_temperature_k\":{},\"exit_pressure_pa\":{},\
         \"mass_flow_kg_per_sec\":{},\"exhaust_velocity_mps\":{},\
         \"momentum_thrust_newtons\":{},\"pressure_thrust_newtons\":{},\
         \"thrust_newtons\":{},\"specific_impulse_sec\":{},\
         \"characteristic_velocity_mps\":{},\"mach_exit\":{}}}}}",
        result.chamber_pressure_psi,
        result.chamber_pressure_pa,
        result.exit_temperature_k,
        result.exit_pressure_pa,
        result.mass_flow_kg_per_sec,
        result.exhaust_velocity_mps,
        result.momentum_thrust_newtons,
        result.pressure_thrust_newtons,
        result.thrust_newtons,
        result.specific_impulse_sec,
        result.characteristic_velocity_mps,
        result.mach_exit,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rocket(values: [f64; 7]) -> String {
        goddard_rocket_step(
            values[0], values[1], values[2], values[3], values[4], values[5], values[6],
        )
    }

    fn inputs() -> [f64; 7] {
        [350.0, 2800.0, 365.0, 1.24, 4.2, 3.5, 0.0]
    }

    #[test]
    fn valid_step_returns_ok_envelope() {
        let result = goddard_rocket_step(350.0, 2800.0, 365.0, 1.24, 4.2, 3.5, 0.0);
        assert!(result.starts_with("{\"ok\":"), "{result}");
        assert!(result.contains("\"thrust_newtons\":"), "{result}");
    }

    #[test]
    fn source_apparatus_step_returns_rigid_body_poses_and_claim_states() {
        let result = goddard_apparatus_step(0.25, 120.0, 6_000.0, 4.5, 0.0, false, true);
        assert!(result.starts_with("{\"ok\":"), "{result}");
        assert!(result.contains("\"primary_quaternion\":"), "{result}");
        assert!(result.contains("\"claim_2_satisfied\":true"), "{result}");
        assert!(result.contains("\"auxiliary_nested\":true"), "{result}");
        assert!(
            !result.contains("thrust") && !result.contains("mach") && !result.contains("liquid"),
            "{result}"
        );
    }

    #[test]
    fn source_apparatus_break_probes_remain_admitted_and_explicit() {
        let result = goddard_apparatus_step(0.25, 120.0, 6_000.0, 2.5, 0.5, false, false);
        assert!(result.contains("\"claim_2_satisfied\":false"), "{result}");
        assert!(
            result.contains("\"claim_1_sequence_satisfied\":false"),
            "{result}"
        );
        assert!(result.contains("\"gyro_enabled\":false"), "{result}");
    }

    #[test]
    fn source_apparatus_refuses_unbounded_or_non_finite_inputs() {
        for result in [
            goddard_apparatus_step(f64::NAN, 120.0, 6_000.0, 4.5, 0.0, false, true),
            goddard_apparatus_step(0.0, 120.0, 6_000.0, 0.5, 0.0, false, true),
            goddard_apparatus_step(0.0, 120.0, 6_000.0, 4.5, 1.1, false, true),
        ] {
            assert!(result.contains("\"refusal\":"), "{result}");
            assert!(result.contains("\"ranked_repairs\":"), "{result}");
        }
    }

    #[test]
    fn non_finite_input_is_a_typed_refusal() {
        for i in 0..7 {
            let mut values = inputs();
            values[i] = f64::NAN;
            let result = rocket(values);
            assert!(result.contains("\"code\":\"non-finite-input\""), "{result}");
            assert!(result.contains("\"ranked_repairs\""), "{result}");
        }
    }

    #[test]
    fn non_positive_and_below_minimum_inputs_are_typed_refusals() {
        for (index, value) in [
            (0, 0.0),
            (1, 0.0),
            (2, 0.0),
            (3, 1.0),
            (4, 0.0),
            (5, MIN_EXPANSION_RATIO - 0.01),
            (6, -1.0),
        ] {
            let mut values = inputs();
            values[index] = value;
            let result = rocket(values);
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }

    #[test]
    fn underflowed_output_refuses_instead_of_returning_zero_as_ok() {
        let mut values = inputs();
        values[4] = f64::from_bits(1);
        let result = rocket(values);
        assert!(
            result.contains("\"code\":\"output-outside-domain\""),
            "{result}"
        );
    }

    #[test]
    fn admission_caps_hold_at_cap_and_refuse_cap_plus_one() {
        let caps = [
            MAX_CHAMBER_PRESSURE_PSI,
            MAX_CHAMBER_TEMPERATURE_K,
            MAX_SPECIFIC_GAS_CONSTANT_J_KG_K,
            5.0 / 3.0,
            MAX_THROAT_AREA_CM2,
            MAX_EXPANSION_RATIO,
            0.0,
        ];
        assert!(rocket(caps).starts_with("{\"ok\":"));
        for i in 0..6 {
            let mut values = caps;
            values[i] += 1.0;
            let result = rocket(values);
            assert!(
                result.contains("\"code\":\"input-outside-domain\""),
                "{result}"
            );
        }
    }

    #[test]
    fn sonic_vacuum_and_overexpanded_states_have_distinct_envelopes() {
        let mut values = inputs();
        values[5] = 1.0;
        let sonic = rocket(values);
        assert!(sonic.starts_with("{\"ok\":"), "{sonic}");
        assert!(sonic.contains("\"mach_exit\":1}"), "{sonic}");
        assert!(sonic.contains("\"mass_flow_kg_per_sec\":"), "{sonic}");
        assert!(sonic.contains("\"pressure_thrust_newtons\":"), "{sonic}");
        values[6] = 2.0e6;
        let refused = rocket(values);
        assert!(
            refused.contains("\"code\":\"unsupported-back-pressure\""),
            "{refused}"
        );
        assert_eq!(rocket(inputs()), rocket(inputs()));
    }
}

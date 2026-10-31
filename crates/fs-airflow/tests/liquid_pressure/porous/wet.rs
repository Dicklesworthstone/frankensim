//! Prescribed liquid saturation through the actual pressure/heat owner.
//! Synthetic calibrated coefficients test equations, not a measured soil.
use super::*;
use fs_material::porous::wet::liquid::{
    LiquidCapillaryPressure, WetPoreLiquidLaw, WetPoreLiquidState, WetPoreLiquidSupport,
};

fn wet_support() -> WetPoreLiquidSupport {
    WetPoreLiquidSupport {
        liquid_saturation: (0.0, 1.0),
        temperature_k: (250.0, 450.0),
        liquid_pressure_pa: (90000.0, 1e9),
        maximum_superficial_speed_m_s: 1.0,
    }
}
fn wet_law(scale: f64, domain: WetPoreLiquidSupport) -> WetPoreLiquidLaw {
    WetPoreLiquidLaw::try_new(
        0.1,
        0.05,
        0.5,
        0.5,
        scale,
        domain,
        "synthetic wetting-liquid branch; no measured-material claim",
    )
    .unwrap()
}
fn liquid_kr(saturation: f64) -> f64 {
    // Independent closed form for this fixture's m=ell=1/2 calibration.
    let effective = (saturation - 0.1) / 0.85;
    effective.sqrt() * (1.0 - (1.0 - effective * effective).sqrt()).powi(2)
}
fn wet_budget() -> PressureDrivenBudget {
    PressureDrivenBudget {
        mass_flow_bracket_kg_per_s: (1e-7, 0.008),
        maximum_neglected_hydraulic_power_w: 1.0,
        ..pressure_budget()
    }
}
fn wet_sections<'a>(
    law: PorousLiquidResistance<'a>,
) -> [LiquidResistanceSection<'static, LiquidResistanceLaw<'a>>; 2] {
    [
        LiquidResistanceSection {
            region: "first",
            resistance: LiquidResistanceLaw::Porous(law),
        },
        LiquidResistanceSection {
            region: "last",
            resistance: LiquidResistanceLaw::Geometry(resistances()[1].resistance),
        },
    ]
}

#[test]
fn liquid_mobility_changes_native_flow_and_retains_original_layers_and_states() {
    let barrier = skeleton(5e9, support());
    let saturations = [0.7, 0.9];
    let factor = 5e9 * (0.3 / liquid_kr(saturations[0]) + 0.7 / liquid_kr(saturations[1]));
    // Changing the capillary calibration at fixed saturation changes Pc, but
    // cannot invent another hydraulic drag or thermal source.
    let mut first_capillary = None;
    for scale in [10000.0, 90000.0] {
        let law = wet_law(scale, wet_support());
        let states = saturations.map(|s| law.sample(s).unwrap());
        let wet = porous(&barrier, 100000.0)
            .with_liquid_states(&states)
            .unwrap();
        with_path(
            false,
            1.0,
            LiquidPropertyTemperature::BulkMean,
            |cx, path| {
                let result = path
                    .solve_pressure_driven(
                        cx,
                        &[300.0; 2],
                        &wet_sections(wet),
                        Pressure::new(50.0),
                        local(),
                        wet_budget(),
                    )
                    .unwrap();
                let expected_mass = 50.0 * 1000.0 / (0.002 * (factor + 5e9));
                close(result.mass_flow_kg_per_s, expected_mass);
                close(result.hydraulics.total_pressure_drop_pa, 50.0);
                close(result.thermal.total_heat_rate_w, 0.0);
                let report = result.hydraulics.regions[0].porous.as_ref().unwrap();
                assert_eq!(report.barrier, barrier);
                assert_eq!(report.liquid_states.as_deref(), Some(states.as_slice()));
                close(report.inlet_pressure_pa, 100050.0);
                let expected_loss = 0.002 * factor * expected_mass / 1000.0;
                close(
                    report.inlet_pressure_pa - report.outlet_pressure_pa,
                    expected_loss,
                );
                close(
                    result.hydraulics.regions[0].loss.dissipated_power_w,
                    expected_loss * expected_mass / 1000.0,
                );
            },
        );
        let LiquidCapillaryPressure::Finite(pc) = states[0].capillary_pressure() else {
            panic!("interior pressure must be finite");
        };
        if let Some(first) = first_capillary {
            close(pc, 9.0 * first);
        } else {
            first_capillary = Some(pc);
        }
    }
}

#[test]
fn changing_saturation_rebinds_original_k_without_compounding_previous_mobility() {
    let barrier = skeleton(5e9, support());
    let law = wet_law(10000.0, wet_support());
    let previous = [law.sample(0.7).unwrap(), law.sample(0.8).unwrap()];
    let next = [law.sample(0.85).unwrap(), law.sample(0.9).unwrap()];
    let saturated = [law.sample(0.95).unwrap(), law.sample(0.95).unwrap()];
    let original = porous(&barrier, 100000.0);
    let first = original.with_liquid_states(&previous).unwrap();
    let rebound = first.with_liquid_states(&next).unwrap();
    let fresh = original.with_liquid_states(&next).unwrap();
    let full = rebound.with_liquid_states(&saturated).unwrap();
    assert_eq!(
        rebound.source().geometric_resistance_per_m(),
        fresh.source().geometric_resistance_per_m()
    );
    close(
        full.source().geometric_resistance_per_m().unwrap(),
        original.source().geometric_resistance_per_m().unwrap(),
    );
    assert_eq!(first.source().liquid_states(), Some(previous.as_slice()));
    assert!(original.source().liquid_states().is_none());
    assert_eq!(rebound.source().barrier(), &barrier);
    with_path(
        false,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| {
            for wet in [rebound, fresh, full] {
                let result = path
                    .solve_pressure_driven(
                        cx,
                        &[300.0; 2],
                        &wet_sections(wet),
                        Pressure::new(50.0),
                        local(),
                        wet_budget(),
                    )
                    .unwrap();
                let factor = wet.source().geometric_resistance_per_m().unwrap() / area().value();
                close(
                    result.mass_flow_kg_per_s,
                    50.0 * 1000.0 / (0.002 * (factor + 5e9)),
                );
            }
        },
    );
}

#[test]
fn wet_mobility_and_original_temperature_viscosity_close_the_same_heat_operating_point() {
    let barrier = skeleton(5e9, support());
    let law = wet_law(10000.0, wet_support());
    let states = [law.sample(0.85).unwrap(), law.sample(0.9).unwrap()];
    let wet = porous(&barrier, 100000.0)
        .with_liquid_states(&states)
        .unwrap();
    let factor = 5e9 * (0.3 / liquid_kr(0.85) + 0.7 / liquid_kr(0.9));
    let mut lower = 1e-7;
    let mut upper = 0.008;
    for _ in 0..80 {
        let mass = 0.5 * (lower + upper);
        if hot_reference(mass, factor).0 < 50.0 {
            lower = mass;
        } else {
            upper = mass;
        }
    }
    let expected = 0.5 * (lower + upper);
    with_path(
        true,
        1.0,
        LiquidPropertyTemperature::FilmMean,
        |cx, path| {
            let result = path
                .solve_pressure_driven(
                    cx,
                    &[400.0, 370.0],
                    &wet_sections(wet),
                    Pressure::new(50.0),
                    local(),
                    wet_budget(),
                )
                .unwrap();
            close(result.mass_flow_kg_per_s, expected);
            close(
                result.thermal.regions[1].outlet.temperature_k(),
                hot_reference(expected, factor).1,
            );
            close(
                result.thermal.total_heat_rate_w,
                expected * 2000.0 * (hot_reference(expected, factor).1 - 300.0),
            );
            assert_eq!(
                result.hydraulics.regions[0]
                    .porous
                    .as_ref()
                    .unwrap()
                    .liquid_states
                    .as_deref(),
                Some(states.as_slice())
            );
        },
    );
}

#[test]
fn dry_liquid_seals_and_wet_source_limits_refuse_without_changing_original_law() {
    let barrier = skeleton(5e9, support());
    let law = wet_law(10000.0, wet_support());
    let dry = [law.sample(0.1).unwrap(), law.sample(0.9).unwrap()];
    let original = porous(&barrier, 100000.0);
    assert!(original.with_liquid_states(&dry).is_err());
    assert!(original.source().liquid_states().is_none());
    assert_eq!(
        dry[0].capillary_pressure(),
        LiquidCapillaryPressure::DryLimit
    );
    assert_eq!(dry[0].liquid_relative_permeability(), 0.0);
    for domain in [
        WetPoreLiquidSupport {
            maximum_superficial_speed_m_s: 0.05,
            ..wet_support()
        },
        WetPoreLiquidSupport {
            liquid_pressure_pa: (90000.0, 100080.0),
            ..wet_support()
        },
    ] {
        let bounded = wet_law(10000.0, domain);
        let states: [WetPoreLiquidState; 2] =
            [bounded.sample(0.95).unwrap(), bounded.sample(0.95).unwrap()];
        let wet = original.with_liquid_states(&states).unwrap();
        with_path(
            false,
            1.0,
            LiquidPropertyTemperature::BulkMean,
            |cx, path| {
                // The root lies in support, but an unsupported bracket trial is
                // rejected before it can lend authority to that root.
                assert!(
                    path.solve_pressure_driven(
                        cx,
                        &[300.0; 2],
                        &wet_sections(wet),
                        Pressure::new(50.0),
                        local(),
                        pressure_budget()
                    )
                    .is_err()
                );
            },
        );
    }
    assert!(original.source().liquid_states().is_none());
}

//! Prescribed liquid saturation through the actual pressure/heat owner.
//! Synthetic calibrated coefficients test equations, not a measured soil.
use super::*;
use fs_material::porous::wet::liquid::{
    LiquidCapillaryPressure, WetPoreLiquidLaw, WetPoreLiquidState, WetPoreLiquidSupport,
    transient::{LiquidSaturationCell, LiquidSaturationColumn, LiquidSaturationStepInput},
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
        // Every flow trial, the lower bracket end first, compiles the laminar thermal card,
        // which admits Re in [1, 2300] (fs-convection CircularDuctLaminarCwt). At 5e-5 kg/s,
        // Re = 5e-5 * 0.01 / (pi 0.01^2 / 4 * 0.002) = 3.18 at the largest viscosity these
        // fixtures reach (0.002 Pa s at 300 K); 1e-7 kg/s was Re 0.0064. The smallest root
        // here, saturation 0.6 in both layers, is 1.37e-4 kg/s.
        mass_flow_bracket_kg_per_s: (5e-5, 0.008),
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


#[test]
fn evolving_saturation_rebinds_the_live_pressure_heat_operating_point() {
    let barrier = skeleton(5e9, support());
    let base = porous(&barrier, 100000.0);
    let law = wet_law(10000.0, wet_support());
    let k = barrier.flow_constraint()
        .permeability(barrier.layers()[0].permeability(), barrier.normal())
        .unwrap();
    let mut column = LiquidSaturationColumn::new(
        barrier.layers().iter().map(|layer| LiquidSaturationCell {
            thickness_m: layer.thickness_m(),
            porosity: 0.4,
            intrinsic_permeability_m2: k,
            law: law.clone(),
        }).collect(),
        vec![0.6, 0.6],
    ).unwrap();

    let initial_states = column.states().to_vec();
    let initial = base.with_saturation_column(&column).unwrap();
    let initial_mass = with_path(
        false,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| path.solve_pressure_driven(
            cx,
            &[300.0; 2],
            &wet_sections(initial),
            Pressure::new(50.0),
            local(),
            wet_budget(),
        ).unwrap().mass_flow_kg_per_s,
    );

    let pc = match column.states()[0].capillary_pressure() {
        LiquidCapillaryPressure::Finite(pc) => pc,
        LiquidCapillaryPressure::DryLimit => panic!("fixture must be mobile"),
    };
    let pl = 150000.0 - pc;
    let candidate = column.step_candidate(
        LiquidSaturationStepInput {
            // Explicit step: only the left face carries flux (both cells and the right boundary
            // sit at pl), k kr / mu * 2000 Pa / (half of 3 mm) = 2.546e-8 * 0.02807 / 0.002 *
            // 1.333e6 = 0.4765 m/s, into cell 0 of capacity 0.4 * 0.003 = 1.2e-3 m. 0.02 s moved
            // its saturation by 7.94, out of the (0, 1) calibration. 5e-5 s moves it by 0.0199
            // (0.600 to 0.620); with dPc/dS = -42 kPa the cell's liquid pressure then rises 0.41
            // of the 2000 Pa push, so the explicit step stays below the boundary that fills it.
            // The later operating point carries 5.4 % more flow (1.365e-4 to 1.440e-4 kg/s).
            duration_s: 5e-5,
            temperature_k: 300.0,
            gas_pressure_pa: 150000.0,
            dynamic_viscosity_pa_s: 0.002,
            body_force_pressure_gradient_pa_m: 0.0,
            boundary_liquid_pressure_pa: [pl + 2000.0, pl],
        },
        || false,
    ).unwrap();
    assert!(candidate.state().saturations()[0] > column.saturations()[0]);
    assert!(candidate.step().balance_residual_m.abs()
        <= 64.0 * f64::EPSILON * candidate.step().liquid_volume_per_area_after_m);
    column = candidate.into_parts().0;

    let rebound = base.with_saturation_column(&column).unwrap();
    let evolved_mass = with_path(
        false,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| path.solve_pressure_driven(
            cx,
            &[300.0; 2],
            &wet_sections(rebound),
            Pressure::new(50.0),
            local(),
            wet_budget(),
        ).unwrap().mass_flow_kg_per_s,
    );
    assert!(evolved_mass > initial_mass);
    assert_ne!(column.states(), initial_states.as_slice());
    assert_eq!(rebound.source().barrier(), &barrier);
    assert_eq!(rebound.source().liquid_states(), Some(column.states()));

    let wrong = LiquidSaturationColumn::new(
        vec![
            LiquidSaturationCell {
                thickness_m: barrier.layers()[0].thickness_m() * 2.0,
                porosity: 0.4,
                intrinsic_permeability_m2: k,
                law: law.clone(),
            },
            LiquidSaturationCell {
                thickness_m: barrier.layers()[1].thickness_m(),
                porosity: 0.4,
                intrinsic_permeability_m2: k,
                law,
            },
        ],
        vec![0.6, 0.6],
    ).unwrap();
    assert!(base.with_saturation_column(&wrong).is_err());
}

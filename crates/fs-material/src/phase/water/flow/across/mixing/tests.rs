use super::*;
use crate::phase::water::pressure::{DenseIsobaricControls, DenseIsobaricWater};
use crate::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::{If97Region, If97SaturationPair, If97State};

fn pressure_controls() -> DensePressureSolve {
    DensePressureSolve {
        pressure_tolerance_pa: 1e-4,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    }
}
fn regional(pressure: f64, minimum: f64) -> RegionalWater {
    RegionalWater::new(
        Pressure::new(pressure),
        RegionalWaterControls {
            minimum_temperature_k: minimum,
            maximum_temperature_k: 873.15,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_saturation_entropy_defect_j_kg_k: 0.1,
            pressure: pressure_controls(),
            junctions: RegionalJunctionControls {
                maximum_enthalpy_jump_j_kg: 200.0,
                maximum_entropy_jump_j_kg_k: 0.25,
                maximum_volume_jump_m3_kg: 1e-5,
                maximum_temperature_gap_k: 1e-8,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap()
}
fn inlet(chart: &RegionalWater, t: f64, mass: f64) -> WaterMixingInlet {
    WaterMixingInlet {
        state: chart
            .state_at_temperature(Temperature::new(t), || false)
            .unwrap()
            .into(),
        mass_flow_kg_per_s: mass,
    }
}
fn controls() -> WaterFlowSolve {
    WaterFlowSolve {
        energy_tolerance_j_kg: 2e-5,
        entropy_tolerance_j_kg_k: 2e-8,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "{a:.15e} != {b:.15e}, tolerance {tolerance:.4e}"
    );
}

#[test]
fn cold_water_and_steam_mix_to_original_latent_quality_with_flow_work_once() {
    let chart = regional(1e6, 273.16);
    let law = WaterFlowChart::from(chart);
    let feeds = [inlet(&chart, 300.0, 2.0), inlet(&chart, 600.0, 1.0)];
    let mixed = law.mix_streams(&feeds, controls(), || false).unwrap();
    // Independent original Gibbs states and saturation lever rule, without
    // calling the mixture's enthalpy inverse or using its returned temperature.
    let p = Pressure::new(1e6);
    let cold = If97State::at_pressure_temperature(p, Temperature::new(300.0), If97Region::Region1)
        .unwrap();
    let steam = If97State::at_pressure_temperature(p, Temperature::new(600.0), If97Region::Region2)
        .unwrap();
    let pair = If97SaturationPair::at_pressure(p).unwrap();
    let h = (2.0 * cold.specific_enthalpy().value() + steam.specific_enthalpy().value()) / 3.0;
    let quality = (h - pair.liquid().specific_enthalpy().value())
        / (pair.vapor().specific_enthalpy().value() - pair.liquid().specific_enthalpy().value());
    assert!(quality > 0.0 && quality < 1.0);
    assert_eq!(mixed.mass_flow_kg_per_s, 3.0);
    near(mixed.outlet.vapor_quality().unwrap(), quality, 2e-14);
    near(
        mixed.outlet.temperature_k(),
        pair.temperature().value(),
        0.0,
    );
    near(mixed.outlet_enthalpy_rate_w, 3.0 * h, 1e-7);
    near(mixed.enthalpy_residual_w, 0.0, 6e-5);
    near(mixed.first_law_residual_w, 0.0, 6e-5);
    assert!(mixed.pressure_flow_work_change_rate_w.abs() > 1000.0);
    near(
        mixed.internal_energy_change_rate_w,
        -mixed.pressure_flow_work_change_rate_w,
        6e-5,
    );
    let s_out = pair.liquid().specific_entropy().value()
        + quality
            * (pair.vapor().specific_entropy().value() - pair.liquid().specific_entropy().value());
    let entropy =
        3.0 * s_out - 2.0 * cold.specific_entropy().value() - steam.specific_entropy().value();
    assert!(entropy > 0.0);
    near(mixed.entropy_generation_w_per_k, entropy, 1e-9);
}

#[test]
fn nonlinear_liquid_mixing_uses_enthalpy_and_scales_as_a_stream_junction() {
    let chart = regional(1e6, 273.16);
    let law = WaterFlowChart::from(chart);
    let feeds = [inlet(&chart, 300.0, 2.0), inlet(&chart, 420.0, 1.0)];
    let mixed = law.mix_streams(&feeds, controls(), || false).unwrap();
    let source = |t| {
        If97State::at_pressure_temperature(
            Pressure::new(1e6),
            Temperature::new(t),
            If97Region::Region1,
        )
        .unwrap()
    };
    let target = (2.0 * source(300.0).specific_enthalpy().value()
        + source(420.0).specific_enthalpy().value())
        / 3.0;
    let (mut low, mut high) = (300.0, 420.0);
    for _ in 0..64 {
        let mid = 0.5 * (low + high);
        if source(mid).specific_enthalpy().value() < target {
            low = mid;
        } else {
            high = mid;
        }
    }
    near(mixed.outlet.temperature_k(), 0.5 * (low + high), 1e-8);
    assert!((mixed.outlet.temperature_k() - 340.0).abs() > 0.01);
    assert_eq!(mixed.outlet.vapor_mass_fraction(), Some(0.0));
    let scaled = law
        .mix_streams(
            &feeds.map(|mut f| {
                f.mass_flow_kg_per_s *= 7.0;
                f
            }),
            controls(),
            || false,
        )
        .unwrap();
    near(
        scaled.outlet.temperature_k(),
        mixed.outlet.temperature_k(),
        1e-8,
    );
    near(
        scaled.inlet_enthalpy_rate_w,
        7.0 * mixed.inlet_enthalpy_rate_w,
        1e-7,
    );
    near(
        scaled.entropy_generation_w_per_k,
        7.0 * mixed.entropy_generation_w_per_k,
        1e-8,
    );
    let allowance = scaled.mass_flow_kg_per_s * controls().energy_tolerance_j_kg;
    near(scaled.enthalpy_residual_w, 0.0, allowance);
    near(scaled.first_law_residual_w, 0.0, allowance);
    let split = [
        feeds[0],
        WaterMixingInlet {
            mass_flow_kg_per_s: 0.25,
            ..feeds[1]
        },
        WaterMixingInlet {
            mass_flow_kg_per_s: 0.75,
            ..feeds[1]
        },
    ];
    let split = law.mix_streams(&split, controls(), || false).unwrap();
    near(
        split.outlet.specific_enthalpy_j_kg(),
        mixed.outlet.specific_enthalpy_j_kg(),
        2e-6,
    );
    let bypass = law.mix_streams(&feeds[..1], controls(), || false).unwrap();
    assert_eq!(bypass.outlet, feeds[0].state);
    assert_eq!(bypass.enthalpy_iterations, 0);
    assert_eq!(bypass.enthalpy_residual_w, 0.0);
    assert_eq!(bypass.entropy_generation_w_per_k, 0.0);
}

#[test]
fn dense_supercritical_mixing_retains_original_source_without_a_vapor_inventory() {
    let chart = DenseIsobaricWater::new(
        Pressure::new(25e6),
        DenseIsobaricControls {
            minimum_temperature_k: 650.0,
            maximum_temperature_k: 670.0,
            enthalpy_tolerance_j_kg: 1e-5,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            pressure: pressure_controls(),
        },
        || false,
    )
    .unwrap();
    let feeds = [(652.0, 1.0), (668.0, 2.0)].map(|(t, mass_flow_kg_per_s)| WaterMixingInlet {
        state: chart
            .state_at_temperature(Temperature::new(t), || false)
            .unwrap()
            .into(),
        mass_flow_kg_per_s,
    });
    let law = WaterFlowChart::from(chart);
    let mixed = law.mix_streams(&feeds, controls(), || false).unwrap();
    assert_eq!(mixed.outlet.law(), law);
    assert_eq!(mixed.outlet.vapor_quality(), None);
    assert_eq!(mixed.outlet.vapor_mass_fraction(), None);
    assert!(mixed.outlet.temperature_k() > 652.0 && mixed.outlet.temperature_k() < 668.0);
    near(
        mixed.outlet_enthalpy_rate_w,
        feeds
            .iter()
            .map(|f| f.mass_flow_kg_per_s * f.state.specific_enthalpy_j_kg())
            .sum(),
        6e-5,
    );
    near(mixed.first_law_residual_w, 0.0, 6e-5);
    assert!(mixed.entropy_generation_w_per_k > 0.0);
}

#[test]
fn mixing_refuses_bad_flows_foreign_charts_and_early_or_late_cancellation() {
    let chart = regional(1e6, 273.16);
    let law = WaterFlowChart::from(chart);
    let feeds = [inlet(&chart, 300.0, 2.0), inlet(&chart, 600.0, 1.0)];
    let original = feeds;
    assert!(law.mix_streams(&[], controls(), || false).is_err());
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.mix_streams(
                &[WaterMixingInlet {
                    mass_flow_kg_per_s: mass,
                    ..feeds[0]
                }],
                controls(),
                || false
            )
            .is_err()
        );
    }
    for foreign in [regional(2e6, 273.16), regional(1e6, 290.0)] {
        assert!(matches!(
            law.mix_streams(&[inlet(&foreign, 350.0, 1.0)], controls(), || false),
            Err(AcrossWaterFlowError::Flow(WaterFlowError::Invalid(_)))
        ));
    }
    assert!(
        law.mix_streams(
            &feeds,
            WaterFlowSolve {
                energy_tolerance_j_kg: -1.0,
                ..controls()
            },
            || false
        )
        .is_err()
    );
    let mut calls = 0;
    let expected = law
        .mix_streams(&feeds, controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    for cancel_at in [1, calls] {
        let mut current = 0;
        assert!(matches!(
            law.mix_streams(&feeds, controls(), || {
                current += 1;
                current == cancel_at
            }),
            Err(AcrossWaterFlowError::Flow(WaterFlowError::Cancelled))
        ));
        assert_eq!(feeds, original);
    }
    assert_eq!(
        law.mix_streams(&feeds, controls(), || false).unwrap(),
        expected
    );
}

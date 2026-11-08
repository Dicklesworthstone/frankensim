//! Original phase calorics for a flowing mixed chamber against a fixed wall.
#![cfg(feature = "material")]

#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{MIXTURE_MASS, model, source};
use fs_gas::chemistry::ThermalConductance;
use fs_gas::material::phase_heat_exchanger::{
    MixtureHeatExchangerError as E, MixtureHeatExchangerInlet as Inlet,
    wall_stream::{
        MixtureWallStreamControl as Control, MixtureWallStreamStep as Step,
        PreparedMixtureWallStream, exchange_mixture_stream_with_wall,
    },
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError, equilibrium::PhaseEnthalpyControls,
};
use fs_qty::{MassFlowRate, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}", actual - expected);
}
fn inlet(model: &ThermochemicalPengRobinsonModel, t: f64, flow: f64, p: f64) -> Inlet<'_> {
    Inlet { model, temperature: Temperature::new(t),
        pressure: Pressure::new(p), mass_flow: MassFlowRate::new(flow) }
}
fn controls(power: [f64; 2], temperature: [f64; 2]) -> Control {
    Control::new(power, PhaseEnthalpyControls::new(temperature.map(Temperature::new)))
}
fn ordinary_controls() -> Control {
    controls([-10_000.0, 10_000.0], [200.0, 240.0])
}
fn ordinary(model: &ThermochemicalPengRobinsonModel, control: Control,
    cancel: impl FnMut() -> bool) -> Result<Step, E> {
    exchange_mixture_stream_with_wall(
        inlet(model, 210.0, 0.5, 3e6), Temperature::new(228.945_170_130_525_04),
        ThermalConductance::new(500.0), control, cancel)
}

#[test]
fn g1_mixed_flow_matches_an_independent_original_two_phase_reference() {
    // Direct PR fugacity roots with analytic departure h/s give this state.
    // The fixed wall equals the independently computed hot mixed-chamber
    // outlet, so the same cold-stream root follows without using its Cp.
    let original = source();
    let result = ordinary(&original, ordinary_controls(), || false).unwrap();
    close(result.enthalpy_flow_change_w, 8401.552_172_252_039, 1e-4);
    close(result.outlet_equilibrium.state().temperature_k(), 212.142_065_786_020_96, 3e-7);
    close(result.outlet_equilibrium.state().vapor_mass_fraction().unwrap(),
        0.120_681_519_207_042_33, 2e-8);
    close(result.entropy_flow_change_w_k, 39.806_095_914_368_5, 1e-6);
    close(result.entropy_production_w_k, 3.109_308_786_419_092, 1e-6);
    assert_eq!(result.outlet_equilibrium.state().constitutive_law(), original.mechanics());
    assert_eq!(result.wall_entropy_flow_change_w_k,
        -result.enthalpy_flow_change_w / 228.945_170_130_525_04);
    assert!(result.enthalpy_target_residual_w.abs() < 2e-5);
    assert!(result.heat_transfer_residual_w.abs() < 2e-5);
}

#[test]
fn g1_mixed_flow_crosses_the_bubble_boundary_without_a_latent_heat_surrogate() {
    let original = source();
    let control = controls([0.0, 28_000.0], [200.0, 240.0]);
    let prepared = PreparedMixtureWallStream::new(
        inlet(&original, 205.0, 0.5, 3e6), control, || false).unwrap();
    let initial = prepared.evaluate(Temperature::new(205.0),
        ThermalConductance::new(1000.0), 1, 0, || false).unwrap();
    assert!(initial.outlet_equilibrium.state().phase().two_phase().is_none());
    let result = prepared.evaluate(Temperature::new(232.486_126_287_924_2),
        ThermalConductance::new(1000.0),
        control.maximum_duty_evaluations, control.maximum_pt_evaluations, || false).unwrap();
    close(result.enthalpy_flow_change_w, 21834.485_817_085_068, 1e-3);
    close(result.outlet_equilibrium.state().temperature_k(), 210.651_640_470_839_14, 3e-7);
    close(result.outlet_equilibrium.state().vapor_mass_fraction().unwrap(),
        0.095_018_183_508_969_51, 3e-8);
    assert!(result.entropy_production_w_k > 0.0);
}

#[test]
fn g1_dilute_heating_and_cooling_recover_the_analytic_mixed_chamber_limit() {
    let original = source();
    let cp = 4.0 * R / MIXTURE_MASS;
    let flow = 0.5;
    let conductance = 500.0;
    for (tin, wall) in [(300.0, 400.0), (400.0, 300.0)] {
        let expected_heat = conductance * (wall - tin) / (1.0 + conductance / (flow * cp));
        let expected_temperature = tin + expected_heat / (flow * cp);
        let expected_entropy = flow * cp * (expected_temperature / tin).ln() - expected_heat / wall;
        let result = exchange_mixture_stream_with_wall(
            inlet(&original, tin, flow, 1.0), Temperature::new(wall),
            ThermalConductance::new(conductance),
            controls([-50_000.0, 50_000.0], [295.0, 405.0]), || false).unwrap();
        // At finite 1-Pa pressure the original PR departure is retained.
        close(result.enthalpy_flow_change_w, expected_heat, 0.01);
        close(result.outlet_equilibrium.state().temperature_k(), expected_temperature, 1e-5);
        close(result.entropy_production_w_k, expected_entropy, 1e-5);
        assert!(result.entropy_production_w_k > 0.0);
    }
}

#[test]
fn g3_species_reference_datums_cancel_before_rounding() {
    let original = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let plain = ordinary(&original, ordinary_controls(), || false).unwrap();
    let gauged = ordinary(&shifted, ordinary_controls(), || false).unwrap();
    assert_eq!(plain.target_duty, gauged.target_duty);
    assert_eq!(plain.enthalpy_flow_change_w, gauged.enthalpy_flow_change_w);
    assert_eq!(plain.entropy_flow_change_w_k, gauged.entropy_flow_change_w_k);
    assert_eq!(plain.entropy_production_w_k, gauged.entropy_production_w_k);
    assert_eq!(plain.outlet_equilibrium.state().temperature_k(),
        gauged.outlet_equilibrium.state().temperature_k());
    assert_eq!(plain.outlet_equilibrium.state().phase(), gauged.outlet_equilibrium.state().phase());
    assert_eq!(plain.pt_evaluations, gauged.pt_evaluations);
}

#[test]
fn g0_prepared_rest_uses_no_endpoint_pt_and_preparation_is_charged_once() {
    let original = source();
    let prepared = PreparedMixtureWallStream::new(
        inlet(&original, 210.0, 0.5, 3e6), ordinary_controls(), || false).unwrap();
    assert_eq!(prepared.preparation_pt_evaluations(), 1);
    for (wall, conductance) in [(210.0, 500.0), (230.0, 0.0)] {
        let rest = prepared.evaluate(Temperature::new(wall), ThermalConductance::new(conductance),
            1, 0, || false).unwrap();
        assert_eq!(rest.enthalpy_flow_change_w, 0.0);
        assert_eq!(rest.entropy_flow_change_w_k, 0.0);
        assert_eq!(rest.entropy_production_w_k, 0.0);
        assert_eq!(rest.pt_evaluations, 0);
        assert_eq!(rest.outlet_equilibrium.pt_evaluations(), 0);
        assert_eq!(rest.duty_evaluations, 1);
    }
    let fresh = ordinary(&original, ordinary_controls(), || false).unwrap();
    let reused = prepared.evaluate(Temperature::new(228.945_170_130_525_04),
        ThermalConductance::new(500.0), ordinary_controls().maximum_duty_evaluations,
        ordinary_controls().maximum_pt_evaluations, || false).unwrap();
    assert_eq!(fresh.enthalpy_flow_change_w, reused.enthalpy_flow_change_w);
    assert_eq!(fresh.pt_evaluations, reused.pt_evaluations + 1);
    assert_eq!(fresh.eos_evaluation_upper_bound,
        reused.eos_evaluation_upper_bound + prepared.preparation_eos_evaluation_upper_bound());
    let mut control = ordinary_controls();
    control.maximum_pt_evaluations = 1;
    control.maximum_duty_evaluations = 1;
    let rest = exchange_mixture_stream_with_wall(
        inlet(&original, 210.0, 0.5, 3e6), Temperature::new(210.0),
        ThermalConductance::new(500.0), control, || false).unwrap();
    assert_eq!(rest.pt_evaluations, 1);
    assert_eq!(rest.enthalpy_flow_change_w, 0.0);
}

#[test]
fn g0_source_limits_shared_work_and_late_cancellation_are_refusals() {
    let original = source();
    let accepted = ordinary(&original, ordinary_controls(), || false).unwrap();
    let mut control = ordinary_controls();
    control.maximum_pt_evaluations = accepted.pt_evaluations;
    control.maximum_duty_evaluations = accepted.duty_evaluations;
    let tight = ordinary(&original, control, || false).unwrap();
    assert_eq!(tight.enthalpy_flow_change_w, accepted.enthalpy_flow_change_w);
    assert_eq!(tight.pt_evaluations, accepted.pt_evaluations);
    control.maximum_pt_evaluations -= 1;
    match ordinary(&original, control, || false).unwrap_err() {
        E::Equilibrium { completed_pt_evaluations,
            source: ThermochemicalFlashError::WorkLimit { flash_evaluations }, .. } => {
            assert_eq!(completed_pt_evaluations + flash_evaluations, control.maximum_pt_evaluations);
        }
        E::WorkLimit { pt_evaluations, .. } => assert_eq!(pt_evaluations, control.maximum_pt_evaluations),
        other => panic!("unexpected PT refusal: {other:?}"),
    }
    control = ordinary_controls();
    control.maximum_duty_evaluations = accepted.duty_evaluations - 1;
    assert!(matches!(ordinary(&original, control, || false),
        Err(E::WorkLimit { duty_evaluations, .. }) if duty_evaluations == control.maximum_duty_evaluations));
    control = ordinary_controls();
    control.equilibrium.enthalpy.temperature_bracket = [Temperature::new(210.0), Temperature::new(211.0)];
    assert!(matches!(ordinary(&original, control, || false),
        Err(E::Equilibrium { source: ThermochemicalFlashError::OutsideEnthalpyBracket, .. })));
    let mut polls = 0;
    ordinary(&original, ordinary_controls(), || { polls += 1; false }).unwrap();
    let mut replay = 0;
    assert_eq!(ordinary(&original, ordinary_controls(), || {
        replay += 1; replay == polls
    }).unwrap_err(), E::Cancelled);
    assert_eq!(ordinary(&original, ordinary_controls(), || false).unwrap(), accepted);
}

//! Steady stream enthalpy, phase change and entropy through two mixed chambers.
#![cfg(feature = "material")]

// Source input data only: the PR fugacity/caloric oracle values below were
// independently calculated without the production flash or scalar root.
#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{MIXTURE_MASS, model, source};
use fs_gas::chemistry::ThermalConductance;
use fs_gas::material::phase_heat_exchanger::{
    MixtureHeatExchangerControl as Control, MixtureHeatExchangerError as E,
    MixtureHeatExchangerInlet as Inlet, MixtureHeatExchangerStep as Step,
    exchange_steady_mixture_heat,
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError, equilibrium::PhaseEnthalpyControls,
};
use fs_qty::{MassFlowRate, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} != {expected}");
}
fn inlet(model: &ThermochemicalPengRobinsonModel, temperature: f64, flow: f64) -> Inlet<'_> {
    Inlet {
        model,
        temperature: Temperature::new(temperature),
        pressure: Pressure::new(3e6),
        mass_flow: MassFlowRate::new(flow),
    }
}
fn controls() -> Control {
    let phase = PhaseEnthalpyControls::new([Temperature::new(200.0), Temperature::new(240.0)]);
    Control::new([-10_000.0, 10_000.0], phase, phase)
}
fn exchange(
    model: &ThermochemicalPengRobinsonModel, control: Control, cancelled: impl FnMut() -> bool,
) -> Result<Step, E> {
    exchange_steady_mixture_heat(
        inlet(model, 210.0, 0.5), inlet(model, 230.0, 1.0),
        ThermalConductance::new(500.0), control, cancelled,
    )
}

#[test]
fn g1_two_phase_outlets_match_independent_fugacity_and_caloric_reference() {
    let original = source();
    let result = exchange(&original, controls(), || false).unwrap();
    let left = result.left_equilibrium.state();
    let right = result.right_equilibrium.state();
    // Independent binary PR flash: both cubic roots, equal species fugacities,
    // material balance, analytic departure h/s, and both flowing H balances.
    // PR constants 0.45724/0.07780; the fixture uses synthetic constant NASA Cp.
    close(result.duty_to_left.value(), 8401.552_172_252_039, 2e-5);
    close(left.temperature_k(), 212.142_065_786_020_96, 2e-7);
    close(right.temperature_k(), 228.945_170_130_525_04, 2e-7);
    close(left.phase().two_phase().unwrap().vapor_fraction_mol(), 0.154_949_049_468_175_75, 2e-8);
    close(right.phase().two_phase().unwrap().vapor_fraction_mol(), 0.456_113_317_904_763_4, 2e-8);
    close(left.vapor_mass_fraction().unwrap(), 0.120_681_519_207_042_33, 2e-8);
    close(right.vapor_mass_fraction().unwrap(), 0.389_653_085_002_195_37, 2e-8);
    close(result.entropy_production_w_k, 3.193_838_500_685_160_4, 2e-7);
    // Total generation also includes mixing of each incoming stream into its
    // chamber, so it exceeds just Q*(1/T_left_out - 1/T_right_out).
    let interface = result.duty_to_left.value() * (1.0 / left.temperature_k() - 1.0 / right.temperature_k());
    assert!(result.entropy_production_w_k > interface + 0.28);
    assert_eq!(left.constitutive_law(), original.mechanics());
    assert_eq!(right.constitutive_law(), original.mechanics());
    close(result.enthalpy_flow_change_w[0], result.duty_to_left.value(), 2e-6);
    close(result.enthalpy_flow_change_w[1], -result.duty_to_left.value(), 2e-6);
    assert!(result.energy_residual_w.abs() <= 2e-5);
    assert!(result.thermal_residual_w.iter().all(|r| r.abs() <= 2e-5));
    // The native enthalpy includes phase conversion, beyond the dilute NASA Cp.
    let cp = 4.0 * R / MIXTURE_MASS;
    let constant_cp_duty = 500.0 * 20.0 / (1.0 + 500.0 * (1.0 / (0.5 * cp) + 1.0 / cp));
    assert!(result.duty_to_left.value() > constant_cp_duty + 3000.0);
}

#[test]
fn g1_dilute_gas_limit_matches_analytic_constant_cp_exchanger() {
    let original = source();
    let phase = PhaseEnthalpyControls::new([Temperature::new(295.0), Temperature::new(405.0)]);
    let control = Control::new([0.0, 50_000.0], phase, phase);
    let mut left = inlet(&original, 300.0, 0.5);
    let mut right = inlet(&original, 400.0, 1.0);
    left.pressure = Pressure::new(1.0);
    right.pressure = Pressure::new(1.0);
    let result = exchange_steady_mixture_heat(left, right, ThermalConductance::new(500.0), control, || false).unwrap();
    let cp = 4.0 * R / MIXTURE_MASS;
    let q = 500.0 * 100.0 / (1.0 + 500.0 * (1.0 / (0.5 * cp) + 1.0 / cp));
    // At finite pressure, the original PR departures remain present. These
    // bounds allow their independently evaluated dilute-limit correction.
    close(result.duty_to_left.value(), 24511.444_464_505_9, 2e-5);
    close(result.duty_to_left.value(), q, 7e-4);
    close(result.left_equilibrium.state().temperature_k(), 300.0 + q / (0.5 * cp), 2e-6);
    close(result.right_equilibrium.state().temperature_k(), 400.0 - q / cp, 2e-6);
    assert!(result.left_equilibrium.state().phase().two_phase().is_none());
    assert!(result.right_equilibrium.state().phase().two_phase().is_none());
    assert!(result.entropy_production_w_k > 0.0);
}

#[test]
fn g1_flowing_heat_crosses_the_original_liquid_to_coexistence_boundary() {
    let original = source();
    let mut control = controls();
    control.duty_bracket_w = [0.0, 30_000.0];
    let left = inlet(&original, 205.0, 0.5);
    let right = inlet(&original, 235.0, 1.0);
    let initial = exchange_steady_mixture_heat(left, right, ThermalConductance::new(0.0), control, || false).unwrap();
    assert!(initial.left_equilibrium.state().phase().two_phase().is_none());
    let result = exchange_steady_mixture_heat(left, right, ThermalConductance::new(1000.0), control, || false).unwrap();
    // Independent PR flash and heat balance. Cold feed is homogeneous liquid;
    // delivered enthalpy produces vapor instead of a constant-Cp temperature rise.
    close(result.duty_to_left.value(), 21834.485_817_085_068, 4e-5);
    close(result.left_equilibrium.state().temperature_k(), 210.651_640_470_839_14, 2e-7);
    close(result.right_equilibrium.state().temperature_k(), 232.486_126_287_924_2, 2e-7);
    close(result.left_equilibrium.state().vapor_mass_fraction().unwrap(), 0.095_018_183_508_969_51, 2e-8);
    close(result.right_equilibrium.state().vapor_mass_fraction().unwrap(), 0.454_715_806_930_457_6, 2e-8);
    close(result.entropy_production_w_k, 11.575_208_290_513_729, 2e-7);
    assert!(result.energy_residual_w.abs() <= 4e-5);
    assert!(result.thermal_residual_w.iter().all(|r| r.abs() <= 4e-5));
}

#[test]
fn g3_independent_species_datums_on_both_streams_cancel_before_rounding() {
    let original = source();
    let a = exchange(&original, controls(), || false).unwrap();
    let shifted_left = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let shifted_right = model([-7e23, 2e22], [-5e18, 3e17], [100.0, 450.0]);
    let b = exchange_steady_mixture_heat(
        inlet(&shifted_left, 210.0, 0.5), inlet(&shifted_right, 230.0, 1.0),
        ThermalConductance::new(500.0), controls(), || false,
    ).unwrap();
    assert_eq!(a.duty_to_left, b.duty_to_left);
    assert_eq!(a.enthalpy_flow_change_w, b.enthalpy_flow_change_w);
    assert_eq!(a.entropy_flow_change_w_k, b.entropy_flow_change_w_k);
    assert_eq!(a.energy_residual_w, b.energy_residual_w);
    assert_eq!(a.thermal_residual_w, b.thermal_residual_w);
    assert_eq!(a.pt_evaluations, b.pt_evaluations);
    assert_eq!(a.eos_evaluation_upper_bound, b.eos_evaluation_upper_bound);
    for (plain, shifted) in [(&a.left_equilibrium, &b.left_equilibrium), (&a.right_equilibrium, &b.right_equilibrium)] {
        assert_eq!(plain.state().temperature_k(), shifted.state().temperature_k());
        assert_eq!(plain.state().phase(), shifted.state().phase());
    }
}

#[test]
fn g1_swapping_streams_reverses_duty_and_preserves_physical_outlets() {
    let original = source();
    let forward = exchange(&original, controls(), || false).unwrap();
    let reverse = exchange_steady_mixture_heat(
        inlet(&original, 230.0, 1.0), inlet(&original, 210.0, 0.5),
        ThermalConductance::new(500.0), controls(), || false,
    ).unwrap();
    close(reverse.duty_to_left.value(), -forward.duty_to_left.value(), 2e-5);
    close(reverse.left_equilibrium.state().temperature_k(), forward.right_equilibrium.state().temperature_k(), 2e-7);
    close(reverse.right_equilibrium.state().temperature_k(), forward.left_equilibrium.state().temperature_k(), 2e-7);
    close(reverse.entropy_production_w_k, forward.entropy_production_w_k, 2e-7);
}

#[test]
fn g0_rest_requires_only_two_original_pt_selections_and_exact_zero_duty() {
    let original = source();
    let mut control = controls();
    control.maximum_pt_evaluations = 2;
    control.maximum_duty_evaluations = 1;
    for (temperature, ua) in [(230.0, 0.0), (210.0, 500.0)] {
        let result = exchange_steady_mixture_heat(
            inlet(&original, 210.0, 0.5), inlet(&original, temperature, 1.0),
            ThermalConductance::new(ua), control, || false,
        ).unwrap();
        assert_eq!(result.duty_to_left.value(), 0.0);
        assert_eq!(result.enthalpy_flow_change_w, [0.0; 2]);
        assert_eq!(result.entropy_flow_change_w_k, [0.0; 2]);
        assert_eq!(result.thermal_residual_w, [0.0; 2]);
        assert_eq!(result.energy_residual_w, 0.0);
        assert_eq!(result.pt_evaluations, 2);
        assert_eq!(result.duty_evaluations, 1);
        assert_eq!(result.left_equilibrium.state().temperature_k(), 210.0);
        assert_eq!(result.right_equilibrium.state().temperature_k(), temperature);
    }
}

#[test]
fn g0_work_limits_cover_both_streams_every_trial_and_repeated_references() {
    let original = source();
    let accepted = exchange(&original, controls(), || false).unwrap();
    let mut control = controls();
    control.maximum_pt_evaluations = accepted.pt_evaluations;
    control.maximum_duty_evaluations = accepted.duty_evaluations;
    assert_eq!(exchange(&original, control, || false).unwrap(), accepted);
    control.maximum_pt_evaluations -= 1;
    match exchange(&original, control, || false).unwrap_err() {
        E::Equilibrium { completed_pt_evaluations, source: ThermochemicalFlashError::WorkLimit { flash_evaluations }, .. } => {
            assert_eq!(completed_pt_evaluations + flash_evaluations, control.maximum_pt_evaluations);
        }
        E::WorkLimit { pt_evaluations, .. } => assert_eq!(pt_evaluations, control.maximum_pt_evaluations),
        other => panic!("unexpected shared PT refusal: {other:?}"),
    }
    control.maximum_pt_evaluations = accepted.pt_evaluations;
    control.maximum_duty_evaluations -= 1;
    assert!(matches!(exchange(&original, control, || false), Err(E::WorkLimit { duty_evaluations, .. })
        if duty_evaluations == control.maximum_duty_evaluations));
}

#[test]
fn g0_invalid_flow_unbracketed_duty_and_original_domain_refuse() {
    let original = source();
    for flow in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(exchange_steady_mixture_heat(
            inlet(&original, 210.0, flow), inlet(&original, 230.0, 1.0),
            ThermalConductance::new(500.0), controls(), || false,
        ), Err(E::Invalid(_))));
    }
    let mut control = controls();
    control.duty_bracket_w = [0.0, 1000.0];
    assert!(matches!(exchange(&original, control, || false), Err(E::NotBracketed { .. })));
    control = controls();
    control.left_equilibrium.enthalpy.temperature_bracket = [Temperature::new(210.0), Temperature::new(211.0)];
    assert!(matches!(exchange(&original, control, || false), Err(E::Equilibrium {
        side: 0, source: ThermochemicalFlashError::OutsideEnthalpyBracket, ..
    })));
}

#[test]
fn g4_final_cancellation_publishes_no_outlet_pair_and_replay_is_exact() {
    let original = source();
    let mut polls = 0;
    let accepted = exchange(&original, controls(), || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert_eq!(exchange(&original, controls(), || { count += 1; count == polls }).unwrap_err(), E::Cancelled);
    assert_eq!(count, polls);
    assert_eq!(exchange(&original, controls(), || false).unwrap(), accepted);
}

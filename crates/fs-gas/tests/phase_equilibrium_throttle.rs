//! G0/G1/G3/G4: original enthalpy, phase change and entropy through a throttle.
#![cfg(feature = "material")]

// Reuse only source input data. The cryogenic anchor below is the independent
// 50-digit simultaneous cubic/fugacity/enthalpy solution retained by fs-material.
#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{METHANE, MIXTURE_MASS, model, source};
use fs_gas::material::phase_throttle::{
    MixtureThrottleControl as Control, MixtureThrottleError as E, MixtureThrottleInlet as Inlet,
    MixtureThrottleStep as Step, throttle_adiabatic_mixture_stream,
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel as Model;
use fs_material::thermochemical_gas::nonideal::flash::ThermochemicalFlashError;
use fs_qty::{MassFlowRate, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}",
        actual - expected,
    );
}
fn controls() -> Control {
    Control::new([Temperature::new(200.0), Temperature::new(235.0)])
}
fn inlet(model: &Model, temperature: f64, pressure: f64, flow: f64) -> Inlet<'_> {
    Inlet {
        model,
        temperature: Temperature::new(temperature),
        pressure: Pressure::new(pressure),
        mass_flow: MassFlowRate::new(flow),
    }
}
fn throttle(model: &Model, control: Control, cancelled: impl FnMut() -> bool) -> Result<Step, E> {
    throttle_adiabatic_mixture_stream(
        inlet(model, 220.0, 3e6, 0.5),
        Pressure::new(2e6),
        control,
        cancelled,
    )
}

#[test]
fn g1_cryogenic_throttle_retains_enthalpy_and_changes_phase_at_independent_endpoint() {
    let original = source();
    let control = controls();
    let result = throttle(&original, control, || false).unwrap();
    let outlet = result.equilibrium.state();
    let initial = original
        .evaluate_phase_equilibrium(
            Temperature::new(220.0),
            Pressure::new(3e6),
            control.equilibrium.phase_options(),
            || false,
        )
        .unwrap();
    // Independent five-variable PR cubic/fugacity/enthalpy equations solved at
    // 50 decimal digits with pressure continuation from 3 MPa to 2 MPa.
    close(outlet.temperature_k(), 209.606_128_901_491_32, 2e-7);
    close(outlet.density_kg_m3(), 74.238_856_531_517_66, 2e-6);
    let split = outlet.phase().two_phase().unwrap();
    close(split.vapor_fraction_mol(), 0.392_287_810_793_803_86, 2e-9);
    close(
        split.liquid_mole_fractions()[METHANE],
        0.290_967_906_203_411_7,
        2e-9,
    );
    close(
        split.vapor_mole_fractions()[METHANE],
        0.823_821_816_126_350_2,
        2e-9,
    );
    close(
        result.entropy_production_w_k,
        0.5 * 47.320_830_152_714_31,
        1e-6,
    );
    assert!(split.vapor_fraction_mol() > initial.phase().two_phase().unwrap().vapor_fraction_mol());
    assert_eq!(outlet.constitutive_law(), original.mechanics());
    assert_eq!(outlet.pressure_pa(), 2e6);
    assert_eq!(result.mass_flow.value(), 0.5);
    assert!(result.energy_residual_w.abs() <= control.absolute_power_tolerance_w);
    // A throttle preserves H, not U. The change of pressure-flow work is
    // already in H and cannot be charged as another heat or shaft-work source.
    let du = outlet.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg();
    let dpv = outlet.pressure_pa() / outlet.density_kg_m3()
        - initial.pressure_pa() / initial.density_kg_m3();
    assert!(du.abs() > 1000.0);
    close(
        du + dpv,
        result.equilibrium.accepted_increment_j_kg().unwrap(),
        1e-8,
    );
    assert!(
        result.equilibrium.pt_evaluations()
            <= control.equilibrium.enthalpy.maximum_flash_evaluations
    );
    assert!(
        result.equilibrium.eos_evaluation_upper_bound()
            <= control.equilibrium.maximum_eos_evaluation_bound().unwrap()
    );
}

#[test]
fn g1_dilute_limit_approaches_isothermal_ideal_gas_throttling_with_positive_entropy() {
    let original = source();
    let control = Control::new([Temperature::new(290.0), Temperature::new(310.0)]);
    let evaluate = |scale| {
        throttle_adiabatic_mixture_stream(
            inlet(&original, 300.0, 10.0 * scale, 0.5),
            Pressure::new(scale),
            control,
            || false,
        )
        .unwrap()
    };
    let first = evaluate(1.0);
    let dilute = evaluate(0.01);
    let ideal_entropy = 0.5 * R / MIXTURE_MASS * 10.0_f64.ln();
    close(dilute.equilibrium.state().temperature_k(), 300.0, 5e-6);
    close(dilute.entropy_production_w_k, ideal_entropy, 1e-5);
    assert!(dilute.equilibrium.state().phase().two_phase().is_none());
    assert!(
        (dilute.equilibrium.state().temperature_k() - 300.0).abs()
            < 0.02 * (first.equilibrium.state().temperature_k() - 300.0).abs()
    );
    assert!(dilute.entropy_production_w_k > 0.0);
}

#[test]
fn g1_throttle_heats_when_the_original_eos_has_negative_joule_thomson_response() {
    // Synthetic calorics remain explicit; this tests the PR process law, not
    // measured high-temperature methane/ethane calorimetry.
    let original = model([0.0; 2], [0.0; 2], [100.0, 2000.0]);
    let result = throttle_adiabatic_mixture_stream(
        inlet(&original, 1600.0, 5e6, 0.5),
        Pressure::new(1e6),
        Control::new([Temperature::new(1550.0), Temperature::new(1650.0)]),
        || false,
    )
    .unwrap();
    // Independent NumPy cubic roots and Brent solve of
    // h=4RT+[RT(Z-1)+(T da/dT-a)/(2 sqrt(2)b) ln((Z+(1+sqrt(2))B)/(Z+(1-sqrt(2))B))].
    // Entropy follows the same independent analytic PR departure and ideal
    // mixture s=4R ln(T)-R ln(p/p_ref)-R sum(x ln(x)).
    close(
        result.equilibrium.state().temperature_k(),
        1602.176_870_748_165_1,
        2e-7,
    );
    close(
        result.equilibrium.state().density_kg_m3(),
        1.726_730_280_214_807_2,
        2e-9,
    );
    close(
        result.entropy_production_w_k,
        0.5 * 583.768_153_846_929_6,
        2e-7,
    );
    assert!(result.equilibrium.state().temperature_k() > 1600.0);
    assert!(result.equilibrium.state().phase().two_phase().is_none());
    assert!(result.energy_residual_w.abs() <= controls().absolute_power_tolerance_w);
}

#[test]
fn g3_reference_datums_cancel_and_throughput_scales_extensive_accounts() {
    let original = source();
    let shifted = model([1e22, -3e21], [2e18, -9e17], [100.0, 450.0]);
    let result = throttle(&original, controls(), || false).unwrap();
    let other = throttle(&shifted, controls(), || false).unwrap();
    assert_eq!(
        result.equilibrium.state().temperature_k(),
        other.equilibrium.state().temperature_k()
    );
    assert_eq!(
        result.equilibrium.state().phase(),
        other.equilibrium.state().phase()
    );
    assert_eq!(
        result.equilibrium.accepted_increment_j_kg(),
        other.equilibrium.accepted_increment_j_kg()
    );
    assert_eq!(result.energy_residual_w, other.energy_residual_w);
    assert_eq!(result.entropy_production_w_k, other.entropy_production_w_k);
    let scaled = throttle_adiabatic_mixture_stream(
        inlet(&original, 220.0, 3e6, 2.0),
        Pressure::new(2e6),
        Control {
            absolute_power_tolerance_w: 4.0 * controls().absolute_power_tolerance_w,
            ..controls()
        },
        || false,
    )
    .unwrap();
    assert_eq!(scaled.equilibrium, result.equilibrium);
    assert_eq!(scaled.mass_flow.value(), 4.0 * result.mass_flow.value());
    assert_eq!(scaled.energy_residual_w, 4.0 * result.energy_residual_w);
    assert_eq!(
        scaled.entropy_production_w_k,
        4.0 * result.entropy_production_w_k
    );
}

#[test]
fn g0_equal_pressure_is_exact_rest_with_one_original_source_query() {
    let original = source();
    let mut control = controls();
    control.equilibrium.enthalpy.maximum_flash_evaluations = 1;
    let result = throttle_adiabatic_mixture_stream(
        inlet(&original, 220.0, 3e6, 0.5),
        Pressure::new(3e6),
        control,
        || false,
    )
    .unwrap();
    assert_eq!(result.equilibrium.state().temperature_k(), 220.0);
    assert_eq!(result.equilibrium.accepted_increment_j_kg(), Some(0.0));
    assert_eq!(result.energy_residual_w, 0.0);
    assert_eq!(result.entropy_production_w_k, 0.0);
    assert_eq!(result.equilibrium.pt_evaluations(), 1);
    assert!(matches!(
        throttle(&original, control, || false),
        Err(E::Equilibrium(ThermochemicalFlashError::WorkLimit { .. }))
    ));
}

#[test]
fn g0_uphill_pressure_invalid_inputs_and_unsupported_caloric_paths_refuse() {
    let original = source();
    for (temperature, inlet_p, outlet_p, flow) in [
        (220.0, 3e6, 4e6, 0.5),
        (220.0, 3e6, 0.0, 0.5),
        (220.0, 0.0, 2e6, 0.5),
        (f64::NAN, 3e6, 2e6, 0.5),
        (220.0, 3e6, 2e6, 0.0),
        (220.0, 3e6, 2e6, -0.5),
        (220.0, 3e6, 2e6, f64::INFINITY),
    ] {
        assert!(matches!(
            throttle_adiabatic_mixture_stream(
                inlet(&original, temperature, inlet_p, flow),
                Pressure::new(outlet_p),
                controls(),
                || false,
            ),
            Err(E::Invalid(_))
        ));
    }
    let narrow = Control::new([Temperature::new(215.0), Temperature::new(225.0)]);
    assert!(matches!(
        throttle(&original, narrow, || false),
        Err(E::Equilibrium(_))
    ));
    for control in [
        Control {
            absolute_power_tolerance_w: 0.0,
            ..controls()
        },
        Control {
            entropy_tolerance_w_k: -1.0,
            ..controls()
        },
        Control {
            entropy_tolerance_w_k: f64::NAN,
            ..controls()
        },
    ] {
        assert!(matches!(
            throttle(&original, control, || false),
            Err(E::Invalid(_))
        ));
    }
}

#[test]
fn g4_cancellation_never_publishes_an_outlet_and_replay_is_identical() {
    let original = source();
    let before = original.clone();
    let mut calls = 0;
    let accepted = throttle(&original, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    for cutoff in [1, calls / 2, calls] {
        let mut current = 0;
        assert_eq!(
            throttle(&original, controls(), || {
                current += 1;
                current == cutoff
            }),
            Err(E::Cancelled)
        );
        assert_eq!(original, before);
    }
    assert_eq!(throttle(&original, controls(), || false).unwrap(), accepted);
}

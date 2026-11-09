//! Original-source throttle/separator/mixer composition and enthalpy balances.
#![cfg(feature = "material")]

#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{METHANE, model, source};
use fs_gas::material::{
    phase_mixing::{MixtureMixingControl, mix_adiabatic_mixture_streams},
    phase_separator::{
        MixtureSeparatorControl as Control, MixtureSeparatorError as E,
        MixtureSeparatorInlet as Inlet, separate_adiabatic_mixture_stream,
        separate_resolved_mixture_stream,
    },
    phase_throttle::{MixtureThrottleControl, throttle_adiabatic_mixture_stream},
};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel as Model,
    flash::equilibrium::{PhaseEnthalpyControls, separation::PhaseSeparationError},
};
use fs_qty::{Amount, MassFlowRate, Pressure, Temperature};
use fs_thermochem::{SpeciesId, eos::peng_robinson::PengRobinsonPhaseV1};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}",
        actual - expected
    );
}
fn calorics() -> PhaseEnthalpyControls {
    PhaseEnthalpyControls::new([Temperature::new(200.0), Temperature::new(235.0)])
}
fn controls() -> Control {
    Control::new(calorics().phase_options())
}
fn inlet(model: &Model, temperature: f64, pressure: f64, flow: f64) -> Inlet<'_> {
    Inlet {
        model,
        temperature: Temperature::new(temperature),
        pressure: Pressure::new(pressure),
        mass_flow: MassFlowRate::new(flow),
    }
}

#[test]
fn throttle_split_and_remix_retain_original_flow_enthalpy_and_species() {
    let model = source();
    let saved = model.clone();
    let throttle = throttle_adiabatic_mixture_stream(
        inlet(&model, 220.0, 3e6, 0.5),
        Pressure::new(2e6),
        MixtureThrottleControl::new(calorics().enthalpy.temperature_bracket),
        || false,
    )
    .unwrap();
    let accepted = throttle.equilibrium.state().clone();
    let result = separate_resolved_mixture_stream(
        &model,
        throttle.equilibrium.into_state(),
        throttle.mass_flow,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.separation.inlet, accepted);
    assert_eq!(
        result.inlet_pt_evaluations, 0,
        "reuse must not repeat the inlet flash"
    );
    assert_eq!(result.separation.additional_eos_evaluations, 2);
    // Prior independent 50-digit PR/fugacity/enthalpy inlet solution, followed
    // by the mass lever rule using methane/ethane molar masses, not beta as quality.
    close(
        result.separation.mole_fractions[1],
        0.392_287_810_793_803_86,
        2e-9,
    );
    close(
        result.separation.mass_fractions[1],
        0.315_004_902_538_153_5,
        2e-9,
    );
    close(
        result.outlet_mass_flows[0].value(),
        0.342_497_548_730_923_25,
        1e-9,
    );
    close(
        result.outlet_mass_flows[1].value(),
        0.157_502_451_269_076_75,
        1e-9,
    );
    let methane = result
        .species
        .iter()
        .find(|s| s.species == SpeciesId::new("methane").unwrap())
        .unwrap();
    close(
        methane.outlet_mass_flow_kg_s[0],
        0.061_518_360_297_995_91,
        1e-9,
    );
    close(
        methane.outlet_mass_flow_kg_s[1],
        0.112_434_558_171_366_4,
        1e-9,
    );
    let mut power = 0.0;
    let mut entropy = 0.0;
    for (index, phase) in result.separation.outlets.iter().enumerate() {
        assert_eq!(phase.temperature_k(), accepted.temperature_k());
        assert_eq!(phase.pressure_pa(), accepted.pressure_pa());
        assert_eq!(
            phase.phase(),
            [PengRobinsonPhaseV1::Liquid, PengRobinsonPhaseV1::Vapor][index]
        );
        close(
            phase.specific_enthalpy_j_kg() - phase.specific_internal_energy_j_kg(),
            phase.pressure_pa() / phase.density_kg_m3(),
            1e-8,
        );
        assert_eq!(
            phase.model().source().thermodynamics().components(),
            model.source().thermodynamics().components()
        );
        assert_eq!(
            phase.model().source().transport_materials(),
            model.source().transport_materials()
        );
        assert_eq!(
            phase.model().mechanics().components(),
            model.mechanics().components()
        );
        assert_eq!(
            phase.model().mechanics().binary_interactions(),
            model.mechanics().binary_interactions()
        );
        power += result.outlet_mass_flows[index].value() * phase.specific_enthalpy_j_kg();
        entropy += result.outlet_mass_flows[index].value() * phase.specific_entropy_j_kg_k();
    }
    close(power, 0.5 * accepted.specific_enthalpy_j_kg(), 1e-7);
    close(entropy, 0.5 * accepted.specific_entropy_j_kg_k(), 1e-8);
    assert!(result.species.iter().all(|s| s.residual_kg_s.abs() < 1e-13));
    assert!(result.mass_residual_kg_s.abs() < 1e-13);
    assert!(result.energy_residual_w.abs() < controls().power_tolerance_w);
    assert!(result.entropy_residual_w_k.abs() < controls().entropy_tolerance_w_k);
    assert!(result.outlet(2).is_none());
    let remix = mix_adiabatic_mixture_streams(
        [result.outlet(0).unwrap(), result.outlet(1).unwrap()],
        None,
        MixtureMixingControl::new(Temperature::new(accepted.temperature_k()), calorics()),
        || false,
    )
    .unwrap();
    close(remix.mass_flow.value(), 0.5, 1e-14);
    close(
        remix.mixture.outlet.state().temperature_k(),
        accepted.temperature_k(),
        3e-7,
    );
    close(
        remix.mixture.outlet.state().density_kg_m3(),
        accepted.density_kg_m3(),
        3e-6,
    );
    close(
        remix.mixture.model.mechanics().mole_fractions()[METHANE],
        0.5,
        1e-13,
    );
    assert!(remix.energy_residual_w.abs() < 2e-5);
    assert!(
        remix.entropy_production_w_k.abs() < 1e-6,
        "routing and recombining already equilibrated phases creates no mixing entropy"
    );
    assert_eq!(model, saved);
}

#[test]
fn direct_split_retains_independent_phase_quality_and_reusable_phase_calorics() {
    let model = source();
    let result =
        separate_adiabatic_mixture_stream(inlet(&model, 220.0, 3e6, 1.0), controls(), || false)
            .unwrap();
    assert_eq!(result.inlet_pt_evaluations, 1);
    close(
        result.separation.mass_fractions[1],
        0.245_214_412_632_239_24,
        2e-10,
    );
    close(
        result.separation.outlets[0]
            .model()
            .mechanics()
            .mole_fractions()[METHANE],
        0.364_782_477_343_904_27,
        2e-10,
    );
    close(
        result.separation.outlets[1]
            .model()
            .mechanics()
            .mole_fractions()[METHANE],
        0.811_655_631_973_835_5,
        2e-10,
    );
    assert!(
        result.separation.outlets[0].density_kg_m3()
            > 5.0 * result.separation.outlets[1].density_kg_m3()
    );
    // The outlet owns all source data needed for later thermal processing.
    let liquid = &result.separation.outlets[0];
    let heated = liquid
        .model()
        .phase_equilibrium_pressure_enthalpy_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            Pressure::new(3e6),
            1000.0,
            PhaseEnthalpyControls::new([Temperature::new(210.0), Temperature::new(250.0)]),
            || false,
        )
        .unwrap();
    close(heated.accepted_increment_j_kg().unwrap(), 1000.0, 1e-5);
    assert!(heated.state().temperature_k() > liquid.temperature_k());
}

#[test]
fn species_energy_entropy_datums_and_throughput_do_not_change_separation_physics() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -9e17], [100.0, 450.0]);
    let a = separate_adiabatic_mixture_stream(inlet(&plain, 220.0, 3e6, 0.5), controls(), || false)
        .unwrap();
    let b =
        separate_adiabatic_mixture_stream(inlet(&shifted, 220.0, 3e6, 0.5), controls(), || false)
            .unwrap();
    assert_eq!(a.outlet_mass_flows, b.outlet_mass_flows);
    assert_eq!(a.species, b.species);
    assert_eq!(a.energy_residual_w, b.energy_residual_w);
    assert_eq!(a.entropy_residual_w_k, b.entropy_residual_w_k);
    for index in 0..2 {
        assert_eq!(
            a.separation.outlets[index].root(),
            b.separation.outlets[index].root()
        );
        assert_eq!(
            a.separation.outlets[index].model().mechanics(),
            b.separation.outlets[index].model().mechanics()
        );
        assert_ne!(
            a.separation.outlets[index].specific_enthalpy_j_kg(),
            b.separation.outlets[index].specific_enthalpy_j_kg()
        );
    }
    let c = separate_adiabatic_mixture_stream(inlet(&plain, 220.0, 3e6, 2.0), controls(), || false)
        .unwrap();
    assert_eq!(c.separation, a.separation);
    assert_eq!(
        c.outlet_mass_flows.map(|x| x.value()),
        a.outlet_mass_flows.map(|x| 4.0 * x.value())
    );
    assert_eq!(c.energy_residual_w, 4.0 * a.energy_residual_w);
    assert_eq!(c.entropy_residual_w_k, 4.0 * a.entropy_residual_w_k);
}

#[test]
fn homogeneous_and_amount_indeterminate_phase_boundaries_refuse() {
    let model = source();
    for temperature in [180.0, 280.0] {
        assert!(matches!(
            separate_adiabatic_mixture_stream(
                inlet(&model, temperature, 3e6, 0.5),
                controls(),
                || false
            ),
            Err(E::Separation(PhaseSeparationError::NoSeparatedPhases))
        ));
    }
    // Pure methane coexistence has distinct roots but PT alone cannot determine
    // their amounts. The source must refuse; a separator cannot invent quality.
    let source = model
        .source()
        .with_mole_amounts(&[Amount::new(0.0), Amount::new(1.0)])
        .unwrap();
    let critical = model.mechanics().components()[METHANE];
    let pure = Model::new(
        source,
        vec![PengRobinsonSpeciesBinding {
            species: SpeciesId::new("methane").unwrap(),
            critical_temperature: critical.critical_temperature,
            critical_pressure: critical.critical_pressure,
            acentric_factor: critical.acentric_factor,
        }],
        vec![0.0],
    )
    .unwrap();
    // Independent Maxwell equal-area pressure also used by the original
    // material/PR phase-screen regression, not a second saturation solve here.
    assert!(
        separate_adiabatic_mixture_stream(
            inlet(&pure, 150.0, 1_047_062.532_316_561_1, 0.5),
            controls(),
            || false,
        )
        .is_err()
    );
}

#[test]
fn mismatched_source_policy_invalid_flow_and_cancellation_publish_neither_outlet() {
    let model = source();
    let saved = model.clone();
    let mut polls = 0;
    let first =
        separate_adiabatic_mixture_stream(inlet(&model, 220.0, 3e6, 0.5), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        assert!(matches!(
            separate_adiabatic_mixture_stream(inlet(&model, 220.0, 3e6, 0.5), controls(), || {
                calls += 1;
                calls == stop
            }),
            Err(E::Cancelled)
        ));
    }
    for flow in [0.0, -1.0, f64::INFINITY] {
        assert!(matches!(
            separate_adiabatic_mixture_stream(inlet(&model, 220.0, 3e6, flow), controls(), || {
                false
            }),
            Err(E::Invalid(_))
        ));
    }
    let other = common::model([1e22, -3e21], [0.0; 2], [100.0, 450.0]);
    assert!(matches!(
        separate_resolved_mixture_stream(
            &other,
            first.separation.inlet.clone(),
            MassFlowRate::new(0.5),
            controls(),
            || false
        ),
        Err(E::Separation(PhaseSeparationError::Invalid(_)))
    ));
    let mut policy = controls();
    policy.phase.flash.log_k_tolerance *= 2.0;
    assert!(matches!(
        separate_resolved_mixture_stream(
            &model,
            first.separation.inlet.clone(),
            MassFlowRate::new(0.5),
            policy,
            || false
        ),
        Err(E::Invalid(_))
    ));
    let mut invalid = controls();
    invalid.power_tolerance_w = 0.0;
    assert!(matches!(
        separate_adiabatic_mixture_stream(inlet(&model, 220.0, 3e6, 0.5), invalid, || false),
        Err(E::Invalid(_))
    ));
    assert_eq!(
        separate_adiabatic_mixture_stream(inlet(&model, 220.0, 3e6, 0.5), controls(), || false)
            .unwrap(),
        first
    );
    assert_eq!(model, saved);
}

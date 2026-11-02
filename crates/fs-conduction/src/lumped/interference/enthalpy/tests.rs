use super::*;
use crate::lumped::{LumpedThermalTransport, solve_enthalpy_step};
use crate::material::CONDUCTIVITY_DIMS;
use crate::radiation::{STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_fft::harmonic_window::{
    HarmonicProductInput, HarmonicProductLaw, HarmonicProductSpectrumLimits,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::Dims;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn phase_curve(identity: ContentHash) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        identity,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 11_300.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 30_000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 11_100.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 55_000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 10_600.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 95_000.0,
                temperature_k: 800.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 10_300.0,
            },
        ],
    )
    .unwrap()
}

fn body(curve: &EquilibriumEnthalpyPhaseCurve) -> LumpedEnthalpyBody<'_> {
    LumpedEnthalpyBody::try_new("specimen", 1.0, 0.01, 0.0, 0.0, 0.001, 35.0, curve).unwrap()
}

fn config() -> LumpedEnthalpyMarchConfig {
    LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 20_000.0,
        ambient_temperature_k: 300.0,
        radiation_temperature_k: 300.0,
        internal_power_w: 0.0,
        duration_s: 1.0,
        maximum_step_s: 1.0,
        maximum_steps: 10_000,
        enthalpy_tolerance_j_kg: 1e-9,
    }
}

fn controls() -> LumpedEnthalpySpectrumControls {
    LumpedEnthalpySpectrumControls {
        start_s: 0.0,
        minimum_step_s: 1e-9,
        source_ripple_tolerance_j: 1.0,
        maximum_step_energy_residual_j: 1e-7,
        maximum_attempts: 20_000,
        maximum_spectral_work: 1_000_000,
        maximum_output_samples: 10_001,
        maximum_workspace_bytes: 8 * 1024 * 1024,
    }
}

fn spectrum(products: Vec<HarmonicProductInput>) -> PreparedHarmonicProductSpectrum {
    let count = products.len();
    PreparedHarmonicProductSpectrum::new(
        products,
        HarmonicProductSpectrumLimits {
            maximum_products: count,
            maximum_generated_terms: 2 * count,
            maximum_workspace_bytes: 64 * 1024,
        },
    )
    .unwrap()
}

fn product(frequency: f64, amplitude: Dd) -> HarmonicProductInput {
    HarmonicProductInput {
        left_frequency_hz: frequency,
        right_frequency_hz: 0.0,
        left_amplitude: [amplitude, Dd::ZERO],
        right_amplitude: [Dd::ONE, Dd::ZERO],
        law: HarmonicProductLaw::Real,
    }
}

fn constant(power: f64) -> PreparedHarmonicProductSpectrum {
    spectrum(vec![product(0.0, Dd::from_f64(power))])
}

fn modulated(dc: f64, amplitude: f64, frequency: f64) -> PreparedHarmonicProductSpectrum {
    spectrum(vec![
        product(0.0, Dd::from_f64(dc)),
        product(frequency, Dd::from_f64(amplitude)),
    ])
}

fn binding<'a>(
    body: &'a LumpedEnthalpyBody<'a>,
    spectrum: &'a PreparedHarmonicProductSpectrum,
) -> LumpedEnthalpyPowerBinding<'a> {
    LumpedEnthalpyPowerBinding {
        body,
        spectrum,
        input_unit: "W",
        watts_per_spectrum_unit: 1.0,
        source_clock: "shared laboratory trigger",
    }
}

fn run<'a>(
    binding: LumpedEnthalpyPowerBinding<'a>,
    config: LumpedEnthalpyMarchConfig,
    controls: LumpedEnthalpySpectrumControls,
) -> LumpedEnthalpySpectrumMarch<'a> {
    solve_lumped_enthalpy_spectrum(binding, BiotGate::corpus_default(), config, controls).unwrap()
}

#[test]
fn actual_native_phase_curve_retains_latent_plateau_under_interference_heating() {
    let curve = phase_curve(ContentHash([0x91; 32]));
    let body =
        LumpedEnthalpyBody::try_new("specimen", 2.0, 0.01, 0.0, 0.0, 0.001, 35.0, &curve).unwrap();
    let source = modulated(500.0, 500.0, 0.03);
    let march = run(
        binding(&body, &source),
        LumpedEnthalpyMarchConfig {
            duration_s: 100.0,
            ..config()
        },
        controls(),
    );
    assert!(core::ptr::eq(march.binding().body, &body));
    assert_eq!(
        march.binding().body.phase_curve().identity(),
        curve.identity()
    );
    assert_eq!(march.gate(), BiotGate::corpus_default());
    assert_eq!(march.maximum_biot(), 0.0);
    let final_state = march.samples().last().unwrap().phase_state();
    near(final_state.specific_enthalpy_j_kg(), 45_000.0, 2e-6);
    assert_eq!(final_state.temperature_k(), 600.0);
    near(final_state.liquid_mass_fraction(), 0.6, 1e-9);
    assert!(march.samples().iter().any(|sample| {
        let state = sample.phase_state();
        state.liquid_mass_fraction() > 0.0 && state.liquid_mass_fraction() < 1.0
    }));
    near(march.energy().source_input_j_dd().to_f64(), 50_000.0, 1e-8);
    near(
        march.energy().stored_energy_change_j_dd().to_f64(),
        50_000.0,
        4e-6,
    );
    assert_eq!(march.energy().convection_into_body_j_dd(), Dd::ZERO);
    assert_eq!(march.energy().radiation_into_body_j_dd(), Dd::ZERO);
    for sample in march.samples().iter().skip(1) {
        assert!(sample.source_ripple_bound_j() <= 1.0);
        assert!(sample.energy().residual_j_dd().to_f64().abs() <= 1e-7);
    }
}

#[test]
fn ripple_admission_catches_carrier_cycles_that_full_and_half_steps_miss() {
    let curve = phase_curve(ContentHash([0x92; 32]));
    let body = body(&curve);
    let source = modulated(0.0, 40_000.0, 1.0);
    // Full and half-cycle endpoint input integrals are all zero: simple step
    // doubling would miss the physically real quarter-cycle latent excursion.
    for (start, duration) in [(0.0, 1.0), (0.0, 0.5), (0.5, 0.5)] {
        near(
            source
                .causal_integral(start, duration, 0.0)
                .unwrap()
                .integral_dd()
                .to_f64(),
            0.0,
            1e-10,
        );
    }
    let march = run(
        binding(&body, &source),
        LumpedEnthalpyMarchConfig {
            initial_specific_enthalpy_j_kg: 25_000.0,
            maximum_steps: 128,
            ..config()
        },
        LumpedEnthalpySpectrumControls {
            source_ripple_tolerance_j: 100.0,
            maximum_output_samples: 129,
            ..controls()
        },
    );
    assert!(march.samples().len() > 3);
    assert!(march.admitted_maximum_step_s() < 0.1);
    let peak = march
        .samples()
        .iter()
        .map(|sample| sample.phase_state().liquid_mass_fraction())
        .fold(0.0, f64::max);
    assert!(peak > 0.05);
    assert!(
        march
            .samples()
            .iter()
            .any(|sample| sample.energy().source_input_j_dd().to_f64() < 0.0)
    );
    near(
        march
            .samples()
            .last()
            .unwrap()
            .phase_state()
            .specific_enthalpy_j_kg(),
        25_000.0,
        1e-7,
    );
    near(march.energy().source_input_j_dd().to_f64(), 0.0, 1e-9);
}

#[test]
fn unresolved_fast_carrier_is_allowed_only_when_its_primitive_is_below_budget() {
    let curve = phase_curve(ContentHash([0x93; 32]));
    let body = body(&curve);
    let source = spectrum(vec![HarmonicProductInput {
        left_frequency_hz: 1e12,
        right_frequency_hz: 1e12,
        left_amplitude: [Dd::ONE, Dd::ZERO],
        right_amplitude: [Dd::ONE, Dd::ZERO],
        law: HarmonicProductLaw::Real,
    }]);
    let march = run(
        LumpedEnthalpyPowerBinding {
            watts_per_spectrum_unit: 1_000.0,
            ..binding(&body, &source)
        },
        LumpedEnthalpyMarchConfig {
            maximum_steps: 1,
            ..config()
        },
        LumpedEnthalpySpectrumControls {
            source_ripple_tolerance_j: 1e-8,
            maximum_attempts: 2,
            maximum_output_samples: 2,
            ..controls()
        },
    );
    assert_eq!(march.samples().len(), 2);
    assert_eq!(march.admitted_maximum_step_s(), 1.0);
    assert_eq!(march.attempts(), 2);
    assert_eq!(march.spectral_work(), 2 * (source.tones().len() + 1));
    near(march.energy().source_input_j_dd().to_f64(), 500.0, 1e-8);
}

#[test]
fn constant_source_uses_original_bisection_and_explicit_volume_scale() {
    let curve = phase_curve(ContentHash([0x94; 32]));
    let body =
        LumpedEnthalpyBody::try_new("specimen", 1.0, 0.01, 2.0, 0.0, 0.001, 35.0, &curve).unwrap();
    let source = constant(2.5e6);
    let cfg = LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 1_000.0,
        maximum_step_s: 0.25,
        ..config()
    };
    let march = run(
        LumpedEnthalpyPowerBinding {
            input_unit: "W/m^3",
            watts_per_spectrum_unit: 2e-6,
            ..binding(&body, &source)
        },
        cfg,
        controls(),
    );
    let mut native_h = cfg.initial_specific_enthalpy_j_kg;
    let native_cfg = LumpedEnthalpyMarchConfig {
        internal_power_w: 5.0,
        ..cfg
    };
    for sample in march.samples().iter().skip(1) {
        native_h =
            solve_enthalpy_step(&body, native_cfg, native_h, sample.step_duration_s()).unwrap();
        near(
            sample.phase_state().specific_enthalpy_j_kg(),
            native_h,
            2e-9,
        );
    }
    near(march.energy().source_input_j_dd().to_f64(), 5.0, 1e-14);
    assert_eq!(march.binding().input_unit, "W/m^3");
}

fn thermal_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, knots) in [
        (
            "thermal_conductivity",
            CONDUCTIVITY_DIMS,
            vec![(300.0, 35.0), (600.0, 20.0), (800.0, 50.0)],
        ),
        (
            SURFACE_EMISSIVITY_PROPERTY,
            Dims::NONE,
            vec![(300.0, 0.2), (600.0, 0.6), (800.0, 0.8)],
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".to_owned(),
                    abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 800.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("synthetic varying {name}"),
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic phase specimen".to_owned(),
            phase: "solid-liquid".to_owned(),
            process: "test".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

// Independent RK4 temporal reference, using the original admitted phase chart
// and native boundary fluxes rather than a surrogate heat capacity.
fn rk4_reference(body: &LumpedEnthalpyBody<'_>, cfg: LumpedEnthalpyMarchConfig) -> f64 {
    let n = 12_000;
    let dt = cfg.duration_s / n as f64;
    let rhs = |t: f64, h: f64| {
        let state = phase_state(body, h).unwrap();
        let (convection, radiation, _) = enthalpy_power(body, cfg, state.temperature_k()).unwrap();
        let power = 20.0 + 15.0 * fs_math::det::cos(core::f64::consts::TAU * 0.2 * t);
        (power + convection + radiation) / body.mass_kg()
    };
    let mut h = cfg.initial_specific_enthalpy_j_kg;
    for step in 0..n {
        let t = step as f64 * dt;
        let k1 = rhs(t, h);
        let k2 = rhs(t + 0.5 * dt, h + 0.5 * dt * k1);
        let k3 = rhs(t + 0.5 * dt, h + 0.5 * dt * k2);
        let k4 = rhs(t + dt, h + dt * k3);
        h += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    h
}

#[test]
fn original_card_transport_and_nonlinear_radiation_converge_under_thermal_refinement() {
    let card = thermal_card();
    let curve = phase_curve(card.content_hash());
    let transport = LumpedThermalTransport::from_material_card(
        &card,
        "thermal_conductivity",
        &[300.0, 800.0],
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(transport.conductivity_receipts().len(), 3);
    let body = LumpedEnthalpyBody::try_new_with_transport(
        "specimen", 0.2, 0.01, 8.0, 0.001, transport, &curve,
    )
    .unwrap();
    let source = modulated(20.0, 15.0, 0.2);
    let cfg = LumpedEnthalpyMarchConfig {
        ambient_temperature_k: 700.0,
        radiation_temperature_k: 900.0,
        duration_s: 30.0,
        ..config()
    };
    let ctl = LumpedEnthalpySpectrumControls {
        source_ripple_tolerance_j: 50.0,
        ..controls()
    };
    let coarse = run(binding(&body, &source), cfg, ctl);
    let fine = run(
        binding(&body, &source),
        LumpedEnthalpyMarchConfig {
            maximum_step_s: 0.125,
            ..cfg
        },
        ctl,
    );
    let reference = rk4_reference(&body, cfg);
    let coarse_h = coarse
        .samples()
        .last()
        .unwrap()
        .phase_state()
        .specific_enthalpy_j_kg();
    let fine_h = fine
        .samples()
        .last()
        .unwrap()
        .phase_state()
        .specific_enthalpy_j_kg();
    assert!((fine_h - reference).abs() < 0.3 * (coarse_h - reference).abs());
    assert_eq!(
        fine.binding().body.transport().material_card_identity(),
        Some(card.content_hash())
    );
    assert_eq!(fine.maximum_biot(), body.maximum_biot(900.0));
    for sample in fine.samples() {
        let temperature = sample.phase_state().temperature_k();
        let emissivity = if temperature <= 600.0 {
            0.2 + (temperature - 300.0) * (0.4 / 300.0)
        } else {
            0.6 + (temperature - 600.0) * (0.2 / 200.0)
        };
        let enclosure_square = 900.0 * 900.0;
        let body_square = temperature * temperature;
        let expected = emissivity
            * STEFAN_BOLTZMANN_W_M2_K4
            * 0.01
            * (enclosure_square * enclosure_square - body_square * body_square);
        near(sample.radiation_into_body_w(), expected, 2e-12);
    }
    assert!(fine.energy().radiation_into_body_j_dd().to_f64() > 0.0);
    assert!(fine.energy().convection_into_body_j_dd().to_f64() > 0.0);
    near(fine.energy().source_input_j_dd().to_f64(), 600.0, 1e-9);
    assert!(fine.cumulative_absolute_energy_residual_j_dd().to_f64() < 1e-4);
}

#[test]
fn signed_extraction_and_both_boundary_reservoirs_are_counted_once() {
    let curve = phase_curve(ContentHash([0x95; 32]));
    let body =
        LumpedEnthalpyBody::try_new("specimen", 2.0, 0.01, 10.0, 0.5, 0.001, 35.0, &curve).unwrap();
    let source = constant(-5.0);
    let march = run(
        binding(&body, &source),
        LumpedEnthalpyMarchConfig {
            initial_specific_enthalpy_j_kg: 10_000.0,
            ambient_temperature_k: 500.0,
            radiation_temperature_k: 700.0,
            duration_s: 2.0,
            maximum_step_s: 0.25,
            ..config()
        },
        controls(),
    );
    let energy = march.energy();
    near(energy.source_input_j_dd().to_f64(), -10.0, 1e-12);
    assert!(energy.convection_into_body_j_dd().to_f64() > 0.0);
    assert!(energy.radiation_into_body_j_dd().to_f64() > 0.0);
    assert!(energy.stored_energy_change_j_dd().to_f64() > 0.0);
    let summed = march.samples().iter().fold(Dd::ZERO, |sum, sample| {
        sum + sample.energy().source_input_j_dd()
    });
    assert_eq!(summed, energy.source_input_j_dd());
    near(energy.residual_j_dd().to_f64(), 0.0, 1e-7);
    let cooling = run(
        binding(&body, &source),
        LumpedEnthalpyMarchConfig {
            initial_specific_enthalpy_j_kg: 75_000.0,
            maximum_step_s: 0.25,
            ..config()
        },
        controls(),
    );
    assert!(cooling.energy().convection_into_body_j_dd().to_f64() < 0.0);
    assert!(cooling.energy().radiation_into_body_j_dd().to_f64() < 0.0);
    assert!(cooling.energy().stored_energy_change_j_dd().to_f64() < 0.0);
}

#[test]
fn source_energy_low_components_and_exact_dark_equilibrium_survive() {
    let curve = phase_curve(ContentHash([0x96; 32]));
    let body = body(&curve);
    let low_source = spectrum(vec![product(0.0, Dd::ONE + Dd::from_f64(1e-20))]);
    let march = run(
        binding(&body, &low_source),
        LumpedEnthalpyMarchConfig {
            maximum_step_s: 0.25,
            ..config()
        },
        controls(),
    );
    near(
        (march.energy().source_input_j_dd() - Dd::ONE).to_f64(),
        1e-20,
        1e-34,
    );
    let cancelled_amplitude = Dd::from_f64(1e16) + Dd::ONE - Dd::from_f64(1e16) - Dd::ONE;
    let dark = spectrum(vec![HarmonicProductInput {
        left_frequency_hz: 10.0,
        right_frequency_hz: 10.0,
        left_amplitude: [cancelled_amplitude, Dd::ZERO],
        right_amplitude: [cancelled_amplitude, Dd::ZERO],
        law: HarmonicProductLaw::Real,
    }]);
    let start_h = 20_000.123_456_789;
    let dark_march = run(
        binding(&body, &dark),
        LumpedEnthalpyMarchConfig {
            initial_specific_enthalpy_j_kg: start_h,
            maximum_step_s: 0.125,
            ..config()
        },
        LumpedEnthalpySpectrumControls {
            source_ripple_tolerance_j: 0.0,
            ..controls()
        },
    );
    for sample in dark_march.samples() {
        assert_eq!(sample.phase_state().specific_enthalpy_j_kg(), start_h);
        assert_eq!(sample.energy(), LumpedEnthalpySpectrumEnergy::ZERO);
    }
}

#[test]
fn late_substeps_keep_their_absolute_source_clock_and_partition_input_exactly() {
    let curve = phase_curve(ContentHash([0x97; 32]));
    let body = body(&curve);
    let source = modulated(1_000.0, 1_000.0, 0.5);
    let cfg = LumpedEnthalpyMarchConfig {
        duration_s: 0.5,
        maximum_step_s: 0.125,
        ..config()
    };
    let ctl = LumpedEnthalpySpectrumControls {
        start_s: 1e16,
        source_ripple_tolerance_j: 10.0,
        ..controls()
    };
    let march = run(binding(&body, &source), cfg, ctl);
    assert_eq!(march.samples().len(), 5);
    let expected = source
        .causal_integral(1e16, 0.5, 0.0)
        .unwrap()
        .integral_dd();
    near(
        (march.energy().source_input_j_dd() - expected).to_f64(),
        0.0,
        1e-9,
    );
    let first = march.samples()[1].energy().source_input_j_dd();
    assert!(
        (march.energy().source_input_j_dd() - Dd::from_f64(4.0) * first)
            .to_f64()
            .abs()
            > 1e-3
    );
    assert_eq!(
        march.samples().last().unwrap().elapsed_s_dd(),
        Dd::from_f64(0.5)
    );
}

#[test]
fn original_biot_phase_domains_and_nonduplicated_sources_remain_mandatory() {
    let curve = phase_curve(ContentHash([0x98; 32]));
    let body = body(&curve);
    let source = constant(1.0);
    assert!(
        solve_lumped_enthalpy_spectrum(
            binding(&body, &source),
            BiotGate::corpus_default(),
            LumpedEnthalpyMarchConfig {
                internal_power_w: 1.0,
                ..config()
            },
            controls(),
        )
        .is_err()
    );
    let high_biot =
        LumpedEnthalpyBody::try_new("specimen", 1.0, 0.01, 100.0, 0.5, 0.1, 0.01, &curve).unwrap();
    assert!(
        solve_lumped_enthalpy_spectrum(
            binding(&high_biot, &source),
            BiotGate::corpus_default(),
            config(),
            controls(),
        )
        .is_err()
    );
    let too_much = constant(1e9);
    assert!(
        solve_lumped_enthalpy_spectrum(
            binding(&body, &too_much),
            BiotGate::corpus_default(),
            config(),
            controls(),
        )
        .is_err()
    );
    let before = body.identity();
    assert!(
        solve_lumped_enthalpy_spectrum(
            binding(&body, &source),
            BiotGate::corpus_default(),
            LumpedEnthalpyMarchConfig {
                initial_specific_enthalpy_j_kg: -1.0,
                ..config()
            },
            controls(),
        )
        .is_err()
    );
    assert_eq!(body.identity(), before);
}

#[test]
fn alias_resolution_work_output_and_energy_tolerances_fail_atomically() {
    let curve = phase_curve(ContentHash([0x99; 32]));
    let body = body(&curve);
    let source = modulated(0.0, 40_000.0, 1.0);
    let original_identity = body.identity();
    for ctl in [
        LumpedEnthalpySpectrumControls {
            minimum_step_s: 0.5,
            source_ripple_tolerance_j: 1.0,
            ..controls()
        },
        LumpedEnthalpySpectrumControls {
            maximum_attempts: 1,
            ..controls()
        },
        LumpedEnthalpySpectrumControls {
            maximum_spectral_work: 1,
            ..controls()
        },
        LumpedEnthalpySpectrumControls {
            maximum_output_samples: 1,
            ..controls()
        },
        LumpedEnthalpySpectrumControls {
            maximum_workspace_bytes: 0,
            ..controls()
        },
    ] {
        assert!(
            solve_lumped_enthalpy_spectrum(
                binding(&body, &source),
                BiotGate::corpus_default(),
                config(),
                ctl,
            )
            .is_err()
        );
    }
    let one = constant(1.0);
    assert!(
        solve_lumped_enthalpy_spectrum(
            binding(&body, &one),
            BiotGate::corpus_default(),
            LumpedEnthalpyMarchConfig {
                duration_s: 2.0,
                maximum_steps: 1,
                ..config()
            },
            controls(),
        )
        .is_err()
    );
    assert!(
        solve_lumped_enthalpy_spectrum(
            binding(&body, &one),
            BiotGate::corpus_default(),
            LumpedEnthalpyMarchConfig {
                enthalpy_tolerance_j_kg: 100.0,
                ..config()
            },
            LumpedEnthalpySpectrumControls {
                maximum_step_energy_residual_j: 1e-12,
                ..controls()
            },
        )
        .is_err()
    );
    assert_eq!(body.identity(), original_identity);
    assert!(
        solve_lumped_enthalpy_spectrum(
            binding(&body, &one),
            BiotGate::corpus_default(),
            config(),
            controls(),
        )
        .is_ok()
    );
}

#[test]
fn metadata_physical_controls_and_zero_duration_are_explicit() {
    let curve = phase_curve(ContentHash([0x9a; 32]));
    let body = body(&curve);
    let source = constant(1.0);
    for invalid in [
        LumpedEnthalpyPowerBinding {
            input_unit: "",
            ..binding(&body, &source)
        },
        LumpedEnthalpyPowerBinding {
            source_clock: " ",
            ..binding(&body, &source)
        },
        LumpedEnthalpyPowerBinding {
            watts_per_spectrum_unit: 0.0,
            ..binding(&body, &source)
        },
        LumpedEnthalpyPowerBinding {
            watts_per_spectrum_unit: f64::INFINITY,
            ..binding(&body, &source)
        },
    ] {
        assert!(
            solve_lumped_enthalpy_spectrum(
                invalid,
                BiotGate::corpus_default(),
                config(),
                controls(),
            )
            .is_err()
        );
    }
    for ctl in [
        LumpedEnthalpySpectrumControls {
            start_s: f64::NAN,
            ..controls()
        },
        LumpedEnthalpySpectrumControls {
            minimum_step_s: 0.0,
            ..controls()
        },
        LumpedEnthalpySpectrumControls {
            source_ripple_tolerance_j: -1.0,
            ..controls()
        },
        LumpedEnthalpySpectrumControls {
            maximum_step_energy_residual_j: f64::NAN,
            ..controls()
        },
    ] {
        assert!(
            solve_lumped_enthalpy_spectrum(
                binding(&body, &source),
                BiotGate::corpus_default(),
                config(),
                ctl,
            )
            .is_err()
        );
    }
    let zero = run(
        binding(&body, &source),
        LumpedEnthalpyMarchConfig {
            duration_s: 0.0,
            ..config()
        },
        LumpedEnthalpySpectrumControls {
            maximum_output_samples: 1,
            ..controls()
        },
    );
    assert_eq!(zero.samples().len(), 1);
    assert_eq!(zero.energy(), LumpedEnthalpySpectrumEnergy::ZERO);
    assert_eq!(zero.spectral_work(), 0);
    assert_eq!(zero.attempts(), 0);
}

#[test]
fn cancellation_inside_native_root_solver_never_publishes_a_partial_march() {
    let curve = phase_curve(ContentHash([0x9b; 32]));
    let body = body(&curve);
    let source = constant(1.0);
    let initial = phase_state(&body, config().initial_specific_enthalpy_j_kg).unwrap();
    let mut calls = 0;
    assert!(matches!(
        solve_lumped_enthalpy_spectrum_cancellable(
            binding(&body, &source),
            BiotGate::corpus_default(),
            config(),
            controls(),
            || {
                calls += 1;
                calls == 12
            },
        ),
        Err(ConductionError::Cancelled { .. })
    ));
    assert_eq!(
        phase_state(&body, config().initial_specific_enthalpy_j_kg).unwrap(),
        initial
    );
    let retry = run(binding(&body, &source), config(), controls());
    assert_eq!(retry.samples()[0].phase_state(), initial);
    near(retry.energy().source_input_j_dd().to_f64(), 1.0, 1e-14);
}

#[test]
fn physical_source_joules_survive_unit_conversion_before_tiny_time_integration() {
    let curve = phase_curve(ContentHash([0x9c; 32]));
    let body = body(&curve);
    let dc = constant(1e-200);
    let harmonic = modulated(0.0, 1e-200, 3.0);
    for source in [&dc, &harmonic] {
        let march = run(
            LumpedEnthalpyPowerBinding {
                input_unit: "declared thermal-power unit",
                watts_per_spectrum_unit: 1e200,
                ..binding(&body, source)
            },
            LumpedEnthalpyMarchConfig {
                duration_s: 1e-200,
                maximum_step_s: 1e-200,
                maximum_steps: 1,
                enthalpy_tolerance_j_kg: 1e-250,
                ..config()
            },
            LumpedEnthalpySpectrumControls {
                minimum_step_s: 1e-201,
                maximum_step_energy_residual_j: 1e-190,
                maximum_output_samples: 2,
                ..controls()
            },
        );
        let energy = march.energy();
        assert!(energy.source_input_j_dd().hi > 0.0);
        near(energy.source_input_j_dd().to_f64(), 1e-200, 4e-215);
        // The original f64 phase coordinate cannot represent this tiny change.
        // Its admitted numerical residual must expose the nonzero input energy,
        // not relabel it as a dark source.
        assert_eq!(energy.stored_energy_change_j_dd(), Dd::ZERO);
        near(energy.residual_j_dd().to_f64(), -1e-200, 4e-215);
        assert_eq!(
            march.samples()[1].energy().source_input_j_dd(),
            energy.source_input_j_dd(),
        );
    }
}

#[test]
fn dd_clock_feedback_intervals_keep_fractional_epoch_and_native_accepted_history() {
    let curve = phase_curve(ContentHash([0xa1; 32]));
    let body = body(&curve);
    let source = modulated(1_000.0, 1_000.0, 0.5);
    let origin = Dd::from_f64(1e16) + Dd::from_f64(0.25);
    assert_eq!(origin.to_f64(), 1e16);
    assert_ne!(origin.lo, 0.0);
    let cfg = LumpedEnthalpyMarchConfig {
        duration_s: 0.5,
        maximum_step_s: 0.125,
        ..config()
    };
    // The DD entrypoint's source_start_s is authoritative; the compatibility
    // controls view must never replace its physically significant low part.
    let ctl = LumpedEnthalpySpectrumControls {
        start_s: -123.0,
        source_ripple_tolerance_j: 100.0,
        ..controls()
    };
    let full = solve_lumped_enthalpy_spectrum_from_clock_dd(
        binding(&body, &source),
        BiotGate::corpus_default(),
        cfg,
        ctl,
        origin,
    )
    .unwrap();
    assert_eq!(full.source_start_s_dd(), origin);
    assert_eq!(full.controls().start_s, origin.to_f64());
    let exact_source = source
        .causal_integral_from_clock_dd(origin, cfg.duration_s, 0.0)
        .unwrap()
        .integral_dd();
    near(
        (full.energy().source_input_j_dd() - exact_source).to_f64(),
        0.0,
        1e-9,
    );
    let rounded = run(
        binding(&body, &source),
        cfg,
        LumpedEnthalpySpectrumControls {
            start_s: origin.to_f64(),
            ..ctl
        },
    );
    assert!(
        (rounded.energy().source_input_j_dd() - full.energy().source_input_j_dd())
            .to_f64()
            .abs()
            > 1.0
    );

    let half_cfg = LumpedEnthalpyMarchConfig {
        duration_s: 0.25,
        ..cfg
    };
    let first = solve_lumped_enthalpy_spectrum_from_clock_dd(
        binding(&body, &source),
        BiotGate::corpus_default(),
        half_cfg,
        ctl,
        origin,
    )
    .unwrap();
    let first_energy = first.energy();
    let first_last = *first.samples().last().unwrap();
    let resumed_origin = first.source_start_s_dd() + first_last.elapsed_s_dd();
    let original_buffer = first.samples().as_ptr();
    let retained = first.into_samples();
    assert_eq!(retained.as_ptr(), original_buffer);
    assert_eq!(retained.as_slice(), &full.samples()[..retained.len()]);
    let resumed = solve_lumped_enthalpy_spectrum_from_clock_dd_cancellable(
        binding(&body, &source),
        BiotGate::corpus_default(),
        LumpedEnthalpyMarchConfig {
            initial_specific_enthalpy_j_kg: first_last.phase_state().specific_enthalpy_j_kg(),
            ..half_cfg
        },
        ctl,
        resumed_origin,
        || false,
    )
    .unwrap();
    assert_eq!(resumed.source_start_s_dd(), origin + Dd::from_f64(0.25));
    assert_eq!(
        resumed.samples().last().unwrap().phase_state(),
        full.samples().last().unwrap().phase_state()
    );
    near(
        (first_energy.source_input_j_dd() + resumed.energy().source_input_j_dd()
            - full.energy().source_input_j_dd())
        .to_f64(),
        0.0,
        1e-25,
    );
    for (resumed_sample, full_sample) in resumed.samples()[1..]
        .iter()
        .zip(&full.samples()[retained.len()..])
    {
        assert_eq!(resumed_sample.phase_state(), full_sample.phase_state());
        assert_eq!(resumed_sample.energy(), full_sample.energy());
        assert_eq!(
            resumed_sample.elapsed_s_dd() + first_last.elapsed_s_dd(),
            full_sample.elapsed_s_dd()
        );
    }
}

#[test]
fn dd_clock_origin_normalizes_and_rejects_nonfinite_or_unrepresentable_pairs() {
    let curve = phase_curve(ContentHash([0xa2; 32]));
    let body = body(&curve);
    let source = constant(0.0);
    let cfg = LumpedEnthalpyMarchConfig {
        duration_s: 0.0,
        ..config()
    };
    for (hi, lo) in [(1e16, 2.25), (-1e16, -2.25), (1e16, -1e16)] {
        let march = solve_lumped_enthalpy_spectrum_from_clock_dd(
            binding(&body, &source),
            BiotGate::corpus_default(),
            cfg,
            controls(),
            Dd { hi, lo },
        )
        .unwrap();
        let canonical = Dd::from_f64(hi) + Dd::from_f64(lo);
        assert_eq!(march.source_start_s_dd(), canonical);
        assert_eq!(march.controls().start_s, canonical.to_f64());
        assert_eq!(march.samples().len(), 1);
        assert_eq!(march.energy(), LumpedEnthalpySpectrumEnergy::ZERO);
    }
    for origin in [
        Dd {
            hi: f64::NAN,
            lo: 0.0,
        },
        Dd {
            hi: 0.0,
            lo: f64::INFINITY,
        },
        Dd {
            hi: f64::NEG_INFINITY,
            lo: 0.0,
        },
        Dd {
            hi: f64::MAX,
            lo: f64::MAX,
        },
    ] {
        assert!(
            solve_lumped_enthalpy_spectrum_from_clock_dd(
                binding(&body, &source),
                BiotGate::corpus_default(),
                cfg,
                controls(),
                origin,
            )
            .is_err()
        );
    }
    assert!(matches!(
        solve_lumped_enthalpy_spectrum_from_clock_dd_cancellable(
            binding(&body, &source),
            BiotGate::corpus_default(),
            cfg,
            controls(),
            Dd::ZERO,
            || true,
        ),
        Err(ConductionError::Cancelled { at: 0, .. })
    ));
}

#[test]
fn original_f64_clock_api_retains_exact_native_results() {
    let curve = phase_curve(ContentHash([0xa3; 32]));
    let body = body(&curve);
    let source = modulated(200.0, 150.0, 0.2);
    let cfg = LumpedEnthalpyMarchConfig {
        duration_s: 0.5,
        maximum_step_s: 0.125,
        ..config()
    };
    let ctl = LumpedEnthalpySpectrumControls {
        start_s: -17.25,
        source_ripple_tolerance_j: 100.0,
        ..controls()
    };
    let original = run(binding(&body, &source), cfg, ctl);
    let dd = solve_lumped_enthalpy_spectrum_from_clock_dd(
        binding(&body, &source),
        BiotGate::corpus_default(),
        cfg,
        ctl,
        Dd::from_f64(ctl.start_s),
    )
    .unwrap();
    assert_eq!(original.source_start_s_dd(), Dd::from_f64(ctl.start_s));
    assert_eq!(dd.samples(), original.samples());
    assert_eq!(dd.energy(), original.energy());
    assert_eq!(dd.attempts(), original.attempts());
    assert_eq!(dd.spectral_work(), original.spectral_work());
    assert_eq!(
        dd.cumulative_absolute_energy_residual_j_dd(),
        original.cumulative_absolute_energy_residual_j_dd()
    );
}

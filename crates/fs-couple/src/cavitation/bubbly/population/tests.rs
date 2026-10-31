use super::*;
use fs_duct::Segment;
use fs_duct::compose;
use fs_flux::cavitation::{BubbleGas, BubbleLimits, BubblePopulationClass, KellerMiksisLimits};

fn material() -> CavitationState {
    CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap()
}
fn limits() -> BubblyLiquidLimits {
    BubblyLiquidLimits {
        maximum_void_fraction: 0.001,
        maximum_wavenumber_radius: 0.1,
        maximum_wavenumber_spacing: 1.0,
        maximum_fractional_radius: 0.01,
        maximum_neighbor_scattered_pressure_ratio: 0.1,
        maximum_multiple_scattering_parameter: 0.1,
    }
}
fn bubble(radius: f64) -> KellerMiksisBubble {
    let liquid = material();
    KellerMiksisBubble::try_new(
        liquid,
        BubbleGas::equilibrium(liquid, radius, 100_000.0, 1.4).unwrap(),
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.5 * radius,
                maximum_radius_m: 2.0 * radius,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap()
}
fn class(radius: f64, fraction: f64) -> BubblePopulationClass {
    let number = fraction / (4.0 * core::f64::consts::PI * radius.powi(3) / 3.0);
    BubblePopulationClass::try_new(bubble(radius), number).unwrap()
}
fn distribution() -> BubblePopulation {
    BubblePopulation::try_new(material(), vec![class(5e-5, 1e-5), class(2e-4, 2e-6)], 2).unwrap()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?}"
    );
}
fn near_scalar(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-30));
}
fn options() -> ComposeOptions {
    ComposeOptions {
        preserve_well_conditioned_bits: false,
        ..ComposeOptions::default()
    }
}

// Reconstruct the complex radial susceptibility from gas/interface properties,
// independently of both the cached native coefficients and coupled helper.
fn independent_chi(model: &KellerMiksisBubble, omega: f64) -> C64 {
    let liquid = model.liquid();
    let r = model.gas().reference_radius_m();
    let rho = liquid.liquid_density_kg_m3();
    let mu = liquid.dynamic_viscosity_pa_s();
    let c = liquid.liquid_sound_speed_m_s();
    let stiffness = (3.0 * model.gas().polytropic_exponent() * model.gas().reference_pressure_pa()
        - 2.0 * liquid.surface_tension_n_m() / r)
        / (rho * r * r);
    let inertia = 1.0 + 4.0 * mu / (rho * c * r);
    let damping = 4.0 * mu / (rho * r * r) + stiffness * r / c;
    -(C64::new(1.0, -omega * r / c)
        / C64::new(stiffness - inertia * omega * omega, -damping * omega))
    .scale(1.0 / (rho * r * r))
}

#[test]
fn zero_distribution_is_exact_original_liquid_and_one_bin_recovers_mono_owner() {
    let source = material();
    let empty = PolydisperseBubblyLiquid::new(
        &source,
        BubblePopulation::try_new(source, vec![], 0).unwrap(),
        limits(),
    )
    .unwrap();
    let pure = empty.response(1e5, 100_000.0).unwrap();
    let native = fluid::segment_wave(&source, 0.005, 1e5).unwrap();
    let wave = pure.segment_wave(0.005).unwrap();
    assert_eq!(wave.wavenumber, native.wavenumber);
    assert_eq!(wave.specific_impedance, native.specific_impedance);
    assert_eq!(
        wave.characteristic_impedance,
        native.characteristic_impedance
    );
    assert_eq!(pure.diagnostics().void_fraction_per_pa, C64::ZERO);
    assert_eq!(pure.diagnostics().mean_bubble_spacing_m, None);

    let member = class(1e-4, 1e-4);
    let mono = DiluteBubblyLiquid::new(&source, *member.bubble(), member.void_fraction(), limits())
        .unwrap();
    let poly = PolydisperseBubblyLiquid::new(
        &source,
        BubblePopulation::try_new(source, vec![member], 1).unwrap(),
        limits(),
    )
    .unwrap();
    let omega = 0.3
        * member
            .linearization()
            .unwrap()
            .stiffness_coefficient_per_s2
            .sqrt();
    let a = mono.response(omega, 1.0).unwrap();
    let b = poly.response(omega, 1.0).unwrap();
    near(b.wavenumber_per_m(), a.wavenumber_per_m(), 4e-15);
    near(
        b.diagnostics().void_fraction_per_pa,
        a.diagnostics()
            .fractional_radius_per_pa
            .scale(3.0 * member.void_fraction()),
        4e-15,
    );
    near_scalar(
        b.diagnostics().neighbor_scattered_pressure_ratio,
        a.diagnostics().neighbor_scattered_pressure_ratio,
        8e-15,
    );
    near_scalar(
        b.diagnostics().multiple_scattering_parameter,
        a.diagnostics().multiple_scattering_parameter,
        8e-15,
    );
}

#[test]
fn splitting_a_class_does_not_change_waves_or_relax_scattering_admission() {
    let source = material();
    let member = class(1e-4, 1e-4);
    let partition =
        BubblePopulationClass::try_new(*member.bubble(), member.number_density_per_m3() / 8.0)
            .unwrap();
    let original = PolydisperseBubblyLiquid::new(
        &source,
        BubblePopulation::try_new(source, vec![member], 1).unwrap(),
        limits(),
    )
    .unwrap();
    let split = PolydisperseBubblyLiquid::new(
        &source,
        BubblePopulation::try_new(source, vec![partition; 8], 8).unwrap(),
        limits(),
    )
    .unwrap();
    let omega = 0.3
        * member
            .linearization()
            .unwrap()
            .stiffness_coefficient_per_s2
            .sqrt();
    let a = original.response(omega, 0.1).unwrap();
    let b = split.response(omega, 0.1).unwrap();
    assert_eq!(a.wavenumber_per_m(), b.wavenumber_per_m());
    assert_eq!(
        a.diagnostics().void_fraction_per_pa,
        b.diagnostics().void_fraction_per_pa
    );
    assert_eq!(
        a.diagnostics().number_density_per_m3,
        b.diagnostics().number_density_per_m3
    );
    assert_eq!(
        a.diagnostics().multiple_scattering_parameter,
        b.diagnostics().multiple_scattering_parameter
    );
    assert_eq!(
        a.diagnostics().neighbor_scattered_pressure_ratio,
        b.diagnostics().neighbor_scattered_pressure_ratio
    );
}

#[test]
fn two_radius_spectrum_matches_independent_foldy_and_sum_of_scattering_cross_sections() {
    let source = material();
    let population = distribution();
    let omega_big = population.classes()[1]
        .linearization()
        .unwrap()
        .stiffness_coefficient_per_s2
        .sqrt();
    let mixture = PolydisperseBubblyLiquid::new(&source, population.clone(), limits()).unwrap();
    for ratio in [0.2, 0.65, 1.8] {
        let omega = ratio * omega_big;
        let result = mixture.response(omega, 0.01).unwrap();
        let mut k_squared = C64::from_re((omega / 1500.0).powi(2));
        let mut volume_response = C64::ZERO;
        let mut total_cross_section = 0.0;
        let mut forward_amplitude_density = C64::ZERO;
        for member in population.classes() {
            let r = member.bubble().gas().reference_radius_m();
            let number = member.number_density_per_m3();
            let chi = independent_chi(member.bubble(), omega);
            let f = chi.scale(-1000.0 * r.powi(3) * omega * omega);
            k_squared = k_squared + f.scale(4.0 * core::f64::consts::PI * number);
            volume_response =
                volume_response + chi.scale(4.0 * core::f64::consts::PI * number * r.powi(3));
            total_cross_section += number * 4.0 * core::f64::consts::PI * f.abs().powi(2);
            forward_amplitude_density = forward_amplitude_density + f.scale(number);
        }
        near(result.wavenumber_per_m(), k_squared.sqrt(), 3e-13);
        near(
            result.diagnostics().void_fraction_per_pa,
            volume_response,
            3e-13,
        );
        near_scalar(
            result.diagnostics().multiple_scattering_parameter,
            total_cross_section / (omega / 1500.0),
            3e-13,
        );
        if ratio > 1.0 {
            assert!(independent_chi(population.classes()[0].bubble(), omega).re < 0.0);
            assert!(independent_chi(population.classes()[1].bubble(), omega).re > 0.0);
            let incorrect_mean_cross_section =
                4.0 * core::f64::consts::PI * forward_amplitude_density.abs().powi(2)
                    / population.number_density_per_m3();
            assert!(total_cross_section > 1.1 * incorrect_mean_cross_section);
        }
    }
}

#[test]
fn low_frequency_compliance_and_bimodal_response_cannot_be_replaced_by_mean_radius() {
    let source = material();
    let population = distribution();
    let big_frequency = population.classes()[1]
        .linearization()
        .unwrap()
        .stiffness_coefficient_per_s2
        .sqrt();
    let mixture = PolydisperseBubblyLiquid::new(&source, population.clone(), limits()).unwrap();
    let omega = 1e-4 * big_frequency;
    let response = mixture.response(omega, 0.0).unwrap();
    let mut static_slowness_squared = 1.0 / 1500.0_f64.powi(2);
    for member in population.classes() {
        let gas = member.bubble().gas();
        let stiffness = 3.0 * gas.polytropic_exponent() * gas.reference_pressure_pa()
            - 2.0 * 0.072 / gas.reference_radius_m();
        static_slowness_squared += 3.0 * 1000.0 * member.void_fraction() / stiffness;
    }
    let actual = (response.wavenumber_per_m() * response.wavenumber_per_m()).re / (omega * omega);
    near_scalar(actual, static_slowness_squared, 1e-8);

    // This one-radius surrogate preserves BOTH total number and total volume.
    // It still cannot preserve the two distinct radial resonance responses.
    let mean_radius = (3.0 * population.void_fraction()
        / (4.0 * core::f64::consts::PI * population.number_density_per_m3()))
    .cbrt();
    let surrogate_class =
        BubblePopulationClass::try_new(bubble(mean_radius), population.number_density_per_m3())
            .unwrap();
    let surrogate = PolydisperseBubblyLiquid::new(
        &source,
        BubblePopulation::try_new(source, vec![surrogate_class], 1).unwrap(),
        limits(),
    )
    .unwrap();
    let frequency = 1.8 * big_frequency;
    let exact = mixture
        .response(frequency, 0.01)
        .unwrap()
        .wavenumber_per_m();
    let collapsed = surrogate
        .response(frequency, 0.01)
        .unwrap()
        .wavenumber_per_m();
    assert!((exact - collapsed).abs() / exact.abs() > 1e-4);
}

#[test]
fn population_wave_uses_native_duct_momentum_phase_and_attenuation() {
    let source = material();
    let population = distribution();
    let frequency = 0.65
        * population.classes()[1]
            .linearization()
            .unwrap()
            .stiffness_coefficient_per_s2
            .sqrt();
    let mixture = PolydisperseBubblyLiquid::new(&source, population, limits()).unwrap();
    let response = mixture.response(frequency, 0.01).unwrap();
    let radius = 0.005;
    let length = 0.1;
    let wave = response.segment_wave(radius).unwrap();
    near(
        wave.wavenumber * wave.specific_impedance,
        C64::from_re(1000.0 * frequency),
        3e-15,
    );
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let chain = response
        .duct_chain(&duct, WallModel::Rigid, ConeSlices::Auto, options(), || {
            false
        })
        .unwrap();
    let solved = compose::compose_chain(
        &chain.blocks,
        Some(wave.characteristic_impedance),
        options(),
    )
    .unwrap();
    let phase = wave.wavenumber.re * length;
    let expected = C64::new(phase.cos(), phase.sin()).scale((-wave.wavenumber.im * length).exp());
    near(solved.z_in, wave.characteristic_impedance, 3e-12);
    near(solved.u_out_over_u_in, expected, 3e-12);
    assert!(expected.abs() < 1.0);
}

#[test]
fn total_void_fraction_and_rare_class_amplitude_cannot_be_hidden_by_weighting() {
    let source = material();
    let too_much =
        BubblePopulation::try_new(source, vec![class(1e-4, 0.0006), class(2e-4, 0.0006)], 2)
            .unwrap();
    assert!(matches!(
        PolydisperseBubblyLiquid::new(&source, too_much, limits()),
        Err(BubblyLiquidError::OutsideModel {
            quantity: "void fraction",
            ..
        })
    ));

    let small = class(5e-5, 1e-4);
    let rare = BubblePopulationClass::try_new(bubble(2e-4), 1e-100).unwrap();
    let omega = 0.98
        * rare
            .linearization()
            .unwrap()
            .stiffness_coefficient_per_s2
            .sqrt();
    let only_small = PolydisperseBubblyLiquid::new(
        &source,
        BubblePopulation::try_new(source, vec![small], 1).unwrap(),
        limits(),
    )
    .unwrap();
    assert!(only_small.response(omega, 2000.0).is_ok());
    let with_rare = PolydisperseBubblyLiquid::new(
        &source,
        BubblePopulation::try_new(source, vec![small, rare], 2).unwrap(),
        limits(),
    )
    .unwrap();
    assert!(matches!(
        with_rare.response(omega, 2000.0),
        Err(BubblyLiquidError::OutsideModel {
            quantity: "fractional radius",
            ..
        })
    ));
}

#[test]
fn cancellation_and_source_mismatch_do_not_mutate_or_reweight_the_distribution() {
    let source = material();
    let population = distribution();
    let changed = CavitationState::try_new(300.0, 999.0, 0.001, 0.072, 2000.0, 1500.0).unwrap();
    assert!(matches!(
        PolydisperseBubblyLiquid::new(&changed, population.clone(), limits()),
        Err(BubblyLiquidError::MaterialMismatch)
    ));
    let owner = PolydisperseBubblyLiquid::new(&source, population.clone(), limits()).unwrap();
    let omega = 0.2
        * population.classes()[1]
            .linearization()
            .unwrap()
            .stiffness_coefficient_per_s2
            .sqrt();
    let accepted = owner.response(omega, 0.1).unwrap();
    for stop in [1, 2, 3, 4] {
        let mut polls = 0;
        assert_eq!(
            owner.response_controlled(omega, 0.1, &mut || {
                polls += 1;
                polls == stop
            }),
            Err(BubblyLiquidError::Bubble(CavitationError::Cancelled))
        );
        assert_eq!(owner.population(), &population);
        assert_eq!(owner.response(omega, 0.1).unwrap(), accepted);
    }
}

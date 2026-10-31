//! Independent dilute-volume, native radial-dynamics and travelling-wave checks.
use super::*;
use fs_duct::compose;
use fs_duct::{DuctError, RadiusProfile, RiccatiSteps, Segment};
use fs_flux::cavitation::{
    BubbleGas, BubbleLimits, BubblePressureSample, BubbleState, KellerMiksisLimits,
};
use fs_material::cavitation::SaturatedWaterCavitation;
use fs_qty::Temperature;
use fs_thermochem::water::if97::If97SaturationPair;
use std::sync::Arc;

const RADIUS: f64 = 1e-4;

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
fn model(state: CavitationState) -> KellerMiksisBubble {
    KellerMiksisBubble::try_new(
        state,
        BubbleGas::equilibrium(state, RADIUS, 100_000.0, 1.4).unwrap(),
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.5 * RADIUS,
                maximum_radius_m: 2.0 * RADIUS,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap()
}
fn fixture(viscosity: f64) -> CavitationState {
    CavitationState::try_new(300.0, 1000.0, viscosity, 0.072, 2000.0, 1500.0).unwrap()
}
fn omega0(bubble: KellerMiksisBubble) -> f64 {
    bubble
        .linearization()
        .unwrap()
        .stiffness_coefficient_per_s2
        .sqrt()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?}"
    );
}
fn options() -> ComposeOptions {
    ComposeOptions {
        preserve_well_conditioned_bits: false,
        ..ComposeOptions::default()
    }
}

#[test]
fn exact_zero_population_reuses_original_water_waves_and_all_native_geometries() {
    let water = SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(Temperature::new(293.15)).unwrap(),
        0.072,
    )
    .unwrap();
    let bubble = model(water.state());
    let empty = DiluteBubblyLiquid::new(&water, bubble, 0.0, limits()).unwrap();
    // No bubbles exist, so bubble compactness and linear-radius limits do not
    // restrict this exact original-fluid endpoint.
    let omega = 20.0 * omega0(bubble);
    let response = empty.response(omega, 100_000.0).unwrap();
    assert_eq!(response.liquid(), water.acoustic_state().unwrap());
    assert_eq!(response.diagnostics().mean_bubble_spacing_m, None);
    let wave = response.segment_wave(0.01).unwrap();
    let original_wave = fluid::segment_wave(&water, 0.01, omega).unwrap();
    assert_eq!(wave.wavenumber, original_wave.wavenumber);
    assert_eq!(wave.specific_impedance, original_wave.specific_impedance);
    assert_eq!(
        wave.characteristic_impedance,
        original_wave.characteristic_impedance
    );
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.01,
                length: 0.1,
            },
            Segment::Cone {
                inlet_radius: 0.01,
                outlet_radius: 0.012,
                length: 0.05,
            },
            Segment::Profile {
                profile: Arc::new(RadiusProfile::exponential(0.012, 0.2, 0.05).unwrap()),
                steps: RiccatiSteps(8),
            },
        ],
    };
    let chain = response
        .duct_chain(&duct, WallModel::Rigid, ConeSlices::Auto, options(), || {
            false
        })
        .unwrap();
    let original = fluid::duct_chain(
        &duct,
        &water,
        omega,
        WallModel::Rigid,
        ConeSlices::Auto,
        options(),
        || false,
    )
    .unwrap();
    assert_eq!(chain.reference_impedance, original.reference_impedance);
    assert_eq!(chain.blocks.len(), original.blocks.len());
    for (actual, expected) in chain.blocks.iter().zip(&original.blocks) {
        assert_eq!(actual.scattering, expected.scattering);
    }
}

#[test]
fn low_frequency_compliance_contains_actual_gas_pressure_and_capillary_stiffness() {
    let source = fixture(0.001);
    let bubble = model(source);
    let alpha = 1e-4;
    let mixture = DiluteBubblyLiquid::new(&source, bubble, alpha, limits()).unwrap();
    let omega = 1e-4 * omega0(bubble);
    let response = mixture.response(omega, 1.0).unwrap();
    // Independent compression of unit bulk volume at leading dilute order:
    // bubble gas pressure changes by -3*kappa*p_g*dR/R; capillarity changes
    // the exterior pressure derivative by +2*sigma*dR/R^2.
    let bubble_modulus = bubble.gas().polytropic_exponent() * bubble.gas().reference_pressure_pa()
        - 2.0 * source.surface_tension_n_m() / (3.0 * RADIUS);
    let static_compliance =
        1.0 / source.acoustic_state().unwrap().adiabatic_bulk_modulus_pa() + alpha / bubble_modulus;
    let expected_k = omega * (source.liquid_density_kg_m3() * static_compliance).sqrt();
    assert!((response.wavenumber_per_m().re / expected_k - 1.0).abs() < 2e-8);
    let static_susceptibility = -1.0 / (3.0 * bubble_modulus);
    assert!(
        (response.diagnostics().fractional_radius_per_pa.re / static_susceptibility - 1.0).abs()
            < 2e-8
    );
    // A low-frequency sound-speed reduction and a stronger dispersive change
    // below resonance both reach the physical propagation output.
    let below = mixture.response(0.8 * omega0(bubble), 1.0).unwrap();
    let slow = omega / response.wavenumber_per_m().re;
    let dispersive = below.omega_rad_s() / below.wavenumber_per_m().re;
    assert!(dispersive < 0.8 * slow && slow < source.liquid_sound_speed_m_s());
}

#[test]
fn coupled_susceptibility_matches_native_compressible_acceleration_and_radiates_without_viscosity()
{
    for viscosity in [0.0, 0.001] {
        let source = fixture(viscosity);
        let bubble = model(source);
        let mixture = DiluteBubblyLiquid::new(&source, bubble, 1e-4, limits()).unwrap();
        let omega = 0.5 * omega0(bubble);
        let response = mixture.response(omega, 0.1).unwrap();
        let chi = response.diagnostics().fractional_radius_per_pa;
        assert!(response.wavenumber_per_m().im > 0.0);
        let equilibrium = bubble
            .linearization()
            .unwrap()
            .equilibrium_absolute_pressure_pa;
        for phase in [0.0_f64, 0.7, 1.9] {
            let eta = 0.1 * (chi.re * phase.cos() + chi.im * phase.sin());
            let velocity = RADIUS * omega * 0.1 * (-chi.re * phase.sin() + chi.im * phase.cos());
            let acceleration = |sign: f64| {
                bubble
                    .observe(
                        &BubbleState::try_new(0.0, RADIUS * (1.0 + sign * eta), sign * velocity)
                            .unwrap(),
                        BubblePressureSample::try_new(
                            equilibrium + sign * 0.1 * phase.cos(),
                            -sign * omega * 0.1 * phase.sin(),
                        )
                        .unwrap(),
                    )
                    .unwrap()
                    .wall_acceleration_m_s2
            };
            let actual = (acceleration(1.0) - acceleration(-1.0)) / 2.0;
            let expected = -RADIUS * omega * omega * eta;
            assert!((actual - expected).abs() < 2e-7 * expected.abs().max(1e-8));
        }
    }
}

#[test]
fn matched_bubbly_pipe_preserves_native_momentum_impedance_phase_and_attenuation() {
    let source = fixture(0.001);
    let bubble = model(source);
    let mixture = DiluteBubblyLiquid::new(&source, bubble, 1e-4, limits()).unwrap();
    let response = mixture.response(0.5 * omega0(bubble), 1.0).unwrap();
    let radius = 0.005;
    let length = 0.1;
    let wave = response.segment_wave(radius).unwrap();
    near(
        wave.wavenumber * wave.specific_impedance,
        C64::from_re(source.liquid_density_kg_m3() * response.omega_rad_s()),
        2e-15,
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
    near(
        solved.p_out_over_u_in,
        expected * wave.characteristic_impedance,
        3e-12,
    );
    assert!(expected.abs() < 1.0);
    assert!(matches!(
        response.duct_chain(&duct, WallModel::Rigid, ConeSlices::Auto, options(), || {
            true
        }),
        Err(BubblyLiquidError::Duct(FluidDuctError::Duct(
            DuctError::Cancelled { .. }
        )))
    ));
    // Public complex-wave entry refuses active propagation, not only this model.
    assert!(
        fluid::segment_wave_with_wavenumber(
            &source,
            radius,
            response.omega_rad_s(),
            C64::new(1.0, -0.1)
        )
        .is_err()
    );
}

#[test]
fn local_amplitude_homogenization_and_resonance_interactions_are_real_admission_limits() {
    let source = fixture(0.001);
    let bubble = model(source);
    let mixture = DiluteBubblyLiquid::new(&source, bubble, 1e-4, limits()).unwrap();
    let omega = omega0(bubble);
    assert!(mixture.response(0.8 * omega, 1.0).is_ok());
    assert!(matches!(
        mixture.response(0.8 * omega, 1e5),
        Err(BubblyLiquidError::OutsideModel {
            quantity: "fractional radius",
            ..
        })
    ));
    assert!(matches!(
        mixture.response(1.2 * omega, 1.0),
        Err(BubblyLiquidError::OutsideModel {
            quantity: "multiple scattering",
            ..
        })
    ));
    assert!(mixture.response(omega, 1e-6).is_err()); // Small drive does not fix strong bubble interactions.
    let sparse = DiluteBubblyLiquid::new(&source, bubble, 1e-12, limits()).unwrap();
    assert!(matches!(
        sparse.response(0.5 * omega, 1.0),
        Err(BubblyLiquidError::OutsideModel {
            quantity: "homogenization",
            ..
        })
    ));
    let different = CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 1490.0).unwrap();
    assert!(matches!(
        DiluteBubblyLiquid::new(&different, bubble, 1e-4, limits()),
        Err(BubblyLiquidError::MaterialMismatch)
    ));
    assert!(DiluteBubblyLiquid::new(&source, bubble, 0.1, limits()).is_err());
}

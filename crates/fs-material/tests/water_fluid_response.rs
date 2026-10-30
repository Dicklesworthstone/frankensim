//! Original pressure-dependent water/steam EOS and transport reach every shared
//! fluid role without independently changing density or imposing a gas law.

use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::transport::If97TransportState;

fn state(p: f64, t: f64) -> If97TransportState {
    If97TransportState::at_pressure_temperature(
        Pressure::new(p),
        Temperature::new(t),
        DensePressureSolve {
            pressure_tolerance_pa: 1e-4,
            density_tolerance_kg_m3: 1e-8,
            maximum_iterations: 128,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn liquid_steam_and_dense_states_share_original_sound_heat_and_momentum_response() {
    for (p, t) in [(1e5, 300.0), (80e6, 300.0), (0.3e6, 650.0), (25e6, 650.0)] {
        let water = state(p, t);
        let original = water.thermodynamics();
        let acoustic = water.acoustic_state().unwrap();
        let transport = water.transport_state().unwrap();
        let thermal = water.thermal_transport_state().unwrap();
        let sound_loss = water.thermoviscous_state().unwrap();
        assert_eq!(acoustic.density_kg_m3(), original.density().value());
        assert_eq!(transport.density_kg_m3(), acoustic.density_kg_m3());
        assert_eq!(acoustic.sound_speed_m_s(), original.sound_speed().value());
        assert_eq!(transport.specific_heat_j_kg_k(), original.cp().value());
        assert_eq!(sound_loss.specific_heat_cv_j_kg_k(), original.cv().value());
        assert_eq!(
            sound_loss.heat_capacity_ratio(),
            original.cp().value() / original.cv().value()
        );
        assert_eq!(sound_loss.transport(), transport);
        assert_eq!(sound_loss.acoustic(), acoustic);
        assert_eq!(thermal.transport(), transport);
        assert_eq!(thermal.temperature_k(), original.temperature().value());
        assert_eq!(
            thermal.volumetric_expansion_per_k(),
            Some(original.thermal_expansion().value())
        );
        assert_eq!(
            transport.thermal_conductivity_w_m_k(),
            water.thermal_conductivity().value()
        );
        assert_eq!(
            transport.dynamic_viscosity_pa_s(),
            water.dynamic_viscosity().value()
        );
        let beta_s =
            original.isothermal_compressibility().value() / sound_loss.heat_capacity_ratio();
        assert!((acoustic.adiabatic_compressibility_per_pa() / beta_s - 1.0).abs() < 3e-10);
    }
}

#[test]
fn pressure_refresh_changes_diffusion_and_sound_and_preserves_signed_liquid_expansion() {
    let low = state(1e5, 300.0);
    let high = state(80e6, 300.0);
    let a = low.thermoviscous_state().unwrap();
    let b = high.thermoviscous_state().unwrap();
    assert!(b.acoustic().adiabatic_bulk_modulus_pa() > a.acoustic().adiabatic_bulk_modulus_pa());
    assert_ne!(
        a.transport().thermal_diffusivity_m2_s(),
        b.transport().thermal_diffusivity_m2_s()
    );
    assert_ne!(
        a.transport().kinematic_viscosity_m2_s(),
        b.transport().kinematic_viscosity_m2_s()
    );
    let cold = state(1e5, 275.0).thermal_transport_state().unwrap();
    assert!(cold.volumetric_expansion_per_k().unwrap() < 0.0);
    let borrowed: &dyn ThermoviscousFluid = &low;
    assert_eq!((&borrowed).thermoviscous_state().unwrap(), a);
}

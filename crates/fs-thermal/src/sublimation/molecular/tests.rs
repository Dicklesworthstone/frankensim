use super::super::{
    MOLAR_GAS_CONSTANT, SublimationError, SublimationLaw, SublimationProperties, SurfaceKinetics,
    VaporReservoir, surface_flux,
};
use super::molecular_mass_flux;

const MOLAR_MASS: f64 = 0.020;

struct ReferenceLaw;

impl SublimationLaw for ReferenceLaw {
    fn temperature_range_k(&self) -> [f64; 2] {
        [150.0, 300.0]
    }

    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        if !t.is_finite() || !(150.0..=300.0).contains(&t) {
            return Err(SublimationError::TemperatureOutOfRange {
                temperature_k: t,
                minimum_k: 150.0,
                maximum_k: 300.0,
            });
        }
        let scale = 500_000.0 * MOLAR_MASS / MOLAR_GAS_CONSTANT;
        let pressure = 100.0 * (scale * (1.0 / 250.0 - 1.0 / t)).exp();
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * scale / t / t,
            latent_heat_j_kg: 500_000.0,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}

fn bath(pressure: f64, temperature: f64) -> VaporReservoir {
    VaporReservoir::Prescribed {
        partial_pressure_pa: pressure,
        temperature_k: temperature,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_detailed_balance_retains_finite_opposing_molecular_fluxes() {
    let alpha = 0.37;
    let t = 250.0;
    let flux = molecular_mass_flux(&ReferenceLaw, t, alpha, bath(100.0, t)).unwrap();
    let expected =
        alpha * 100.0 * (MOLAR_MASS / (2.0 * std::f64::consts::PI * MOLAR_GAS_CONSTANT * t)).sqrt();
    close(
        flux.outgoing_kg_m2_s,
        expected,
        16.0 * f64::EPSILON * expected,
    );
    assert_eq!(flux.outgoing_kg_m2_s, flux.incoming_kg_m2_s);
    assert_eq!(flux.net_kg_m2_s, 0.0);
    assert!(flux.outgoing_kg_m2_s > 0.0);
    let sealed = molecular_mass_flux(&ReferenceLaw, 220.0, 0.0, bath(500.0, 280.0)).unwrap();
    assert_eq!(sealed.outgoing_kg_m2_s, 0.0);
    assert_eq!(sealed.incoming_kg_m2_s, 0.0);
    assert_eq!(sealed.net_kg_m2_s, 0.0);
}

#[test]
fn g1_unequal_temperature_zero_net_mass_still_resolves_energy_exchange() {
    let ts = 220.0;
    let tg: f64 = 280.0;
    let pressure = ReferenceLaw.properties(ts).unwrap().saturation_pressure_pa * (tg / ts).sqrt();
    let flux = molecular_mass_flux(&ReferenceLaw, ts, 0.63, bath(pressure, tg)).unwrap();
    assert!(flux.outgoing_kg_m2_s > 0.0 && flux.incoming_kg_m2_s > 0.0);
    close(
        flux.net_kg_m2_s,
        0.0,
        16.0 * f64::EPSILON * flux.outgoing_kg_m2_s,
    );
    // Independent calorically-perfect example: each Maxwellian half-space
    // carries u + Rs*T/2. Net mass alone would incorrectly report zero energy.
    let effective_capacity = 1000.0 + 0.5 * MOLAR_GAS_CONSTANT / MOLAR_MASS;
    let energy_flux = flux.outgoing_kg_m2_s * effective_capacity * ts
        - flux.incoming_kg_m2_s * effective_capacity * tg;
    assert!(energy_flux < 0.0);
    close(
        energy_flux,
        flux.outgoing_kg_m2_s * effective_capacity * (ts - tg),
        1.0e-8,
    );
}

#[test]
fn g1_net_flux_agrees_with_original_surface_flux_for_both_directions() {
    let alpha = 0.41;
    for (ts, reservoir) in [
        (250.0, bath(0.0, 250.0)),
        (250.0, bath(50.0, 250.0)),
        (250.0, bath(200.0, 250.0)),
        (220.0, bath(30.0, 280.0)),
        (220.0, bath(300.0, 280.0)),
        (250.0, bath(100.0_f64.next_down(), 250.0)),
        (
            250.0,
            VaporReservoir::Finite {
                vapor_mass_kg: 1.0e-5,
                volume_m3: 0.1,
                temperature_k: 260.0,
            },
        ),
    ] {
        let gross = molecular_mass_flux(&ReferenceLaw, ts, alpha, reservoir).unwrap();
        let signed = surface_flux(
            &ReferenceLaw,
            ts,
            SurfaceKinetics {
                accommodation_coefficient: alpha,
                mass_transfer_coefficient_m_s: None,
            },
            reservoir,
        )
        .unwrap();
        assert_eq!(gross.net_kg_m2_s, signed.mass_flux_kg_m2_s);
        assert!(gross.outgoing_kg_m2_s >= 0.0);
        assert!(gross.incoming_kg_m2_s >= 0.0);
    }
}

#[test]
fn g0_molecular_flux_refuses_invalid_accommodation_pressure_and_domains() {
    for alpha in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(molecular_mass_flux(&ReferenceLaw, 250.0, alpha, bath(100.0, 250.0)).is_err());
    }
    for (ts, pressure, tg) in [
        (149.0, 100.0, 250.0),
        (301.0, 100.0, 250.0),
        (f64::NAN, 100.0, 250.0),
        (250.0, -1.0, 250.0),
        (250.0, f64::NAN, 250.0),
        (250.0, 100.0, 0.0),
        (250.0, 100.0, f64::INFINITY),
    ] {
        assert!(molecular_mass_flux(&ReferenceLaw, ts, 1.0, bath(pressure, tg)).is_err());
    }
}

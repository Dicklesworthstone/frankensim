use super::*;
use crate::sublimation::{
    MOLAR_GAS_CONSTANT, SublimationLaw, SublimationProperties, VaporReservoir,
    molecular::molecular_mass_flux,
};

const MOLAR_MASS: f64 = 0.018;
const R: f64 = MOLAR_GAS_CONSTANT / MOLAR_MASS;

struct Phase(f64);
impl SublimationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [1.0, 2000.0]
    }
    fn properties(&self, _: f64) -> Result<SublimationProperties, SublimationError> {
        Ok(SublimationProperties {
            saturation_pressure_pa: self.0,
            saturation_pressure_derivative_pa_k: 0.0,
            latent_heat_j_kg: 2e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 917.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}
fn flux(p_sat: f64, p_gas: f64, ts: f64, tg: f64, sticking: f64) -> MolecularMassFlux {
    molecular_mass_flux(
        &Phase(p_sat),
        ts,
        sticking,
        VaporReservoir::Prescribed {
            partial_pressure_pa: p_gas,
            temperature_k: tg,
        },
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

// Independent Simpson integration of the normalized one-dimensional Maxwellian.
// Crossing mass weights z*f(z); normal momentum weights z^2*f(z).
fn half_space_moments() -> (f64, f64) {
    let n = 8192;
    let step = 12.0 / n as f64;
    let (mut mass, mut momentum) = (0.0, 0.0);
    for i in 0..=n {
        let z = i as f64 * step;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let density = (-0.5 * z * z).exp() / (2.0 * core::f64::consts::PI).sqrt();
        mass += weight * z * density;
        momentum += weight * z * z * density;
    }
    (mass * step / 3.0, momentum * step / 3.0)
}

#[test]
fn original_gross_beams_match_independent_half_space_integrals() {
    let (i1, i2) = half_space_moments();
    near(i2, 0.5, 3e-14);
    for (ts, tg, p_sat, p_gas, sticking) in [
        (220.0, 350.0, 14.0, 37.0, 0.3),
        (600.0, 200.0, 900.0, 2.0, 1.0),
    ] {
        let actual =
            phase_normal_momentum(flux(p_sat, p_gas, ts, tg, sticking), R, ts, tg).unwrap();
        let source_density = p_sat / (R * ts);
        let gas_density = p_gas / (R * tg);
        near(
            actual.outgoing_pa,
            sticking * source_density * R * ts * i2,
            3e-11,
        );
        near(
            actual.incoming_pa,
            sticking * gas_density * R * tg * i2,
            3e-11,
        );
        near(
            flux_weighted_normal_speed(R, ts).unwrap(),
            (R * ts).sqrt() * i2 / i1,
            2e-10,
        );
    }
}

#[test]
fn equilibrium_and_sealed_specular_interfaces_retain_pressure() {
    let p = 137.0;
    for beta in [0.0, 0.01, 0.5, 1.0] {
        let mass = flux(p, p, 300.0, 300.0, beta);
        assert_eq!(mass.net_kg_m2_s, 0.0);
        let phase = phase_normal_momentum(mass, R, 300.0, 300.0).unwrap();
        for chi in [0.0, 0.2, 1.0] {
            let reflected = maxwell_reflected_normal_pressure(p, 300.0, 300.0, beta, chi).unwrap();
            near(phase.total_pa + reflected, p, 6e-14);
        }
    }
    let mass = flux(1000.0, p, 600.0, 200.0, 0.0);
    assert_eq!(
        phase_normal_momentum(mass, R, 600.0, 200.0)
            .unwrap()
            .total_pa,
        0.0
    );
    assert_eq!(
        maxwell_reflected_normal_pressure(p, 600.0, 200.0, 0.0, 0.0).unwrap(),
        p
    );
}

#[test]
fn vacuum_emission_and_zero_net_mass_have_nonzero_recoil() {
    let phase =
        phase_normal_momentum(flux(400.0, 0.0, 500.0, 300.0, 0.4), R, 500.0, 300.0).unwrap();
    near(phase.outgoing_pa, 80.0, 3e-14);
    assert_eq!(phase.incoming_pa, 0.0);
    assert_eq!(
        maxwell_reflected_normal_pressure(0.0, 500.0, 300.0, 0.4, 1.0).unwrap(),
        0.0
    );
    // p_sat/sqrt(Ts) = p_gas/sqrt(Tg): equal opposing mass is not zero momentum.
    let mass = flux(200.0, 100.0, 800.0, 200.0, 0.5);
    near(mass.net_kg_m2_s, 0.0, 1e-16);
    let phase = phase_normal_momentum(mass, R, 800.0, 200.0).unwrap();
    near(phase.total_pa, 75.0, 3e-14);
}

#[test]
fn explicit_maxwell_reflection_matches_incident_and_reemitted_integrals() {
    let (i1, i2) = half_space_moments();
    for (ts, tg, beta, chi) in [(500.0, 200.0, 0.3, 0.6), (200.0, 500.0, 0.0, 1.0)] {
        let p = 77.0;
        let incident_mass = (1.0 - beta) * (p / (R * tg)) * (R * tg).sqrt() * i1;
        let incident_momentum = (1.0 - beta) * p * i2;
        let diffuse_momentum = incident_mass * (R * ts).sqrt() * i2 / i1;
        let expected = incident_momentum + (1.0 - chi) * incident_momentum + chi * diffuse_momentum;
        near(
            maxwell_reflected_normal_pressure(p, ts, tg, beta, chi).unwrap(),
            expected,
            2e-11,
        );
    }
}

#[test]
fn finite_extremes_avoid_squaring_temperatures_and_preserve_empty_beams() {
    let coefficient = core::f64::consts::FRAC_PI_2.sqrt();
    near(
        flux_weighted_normal_speed(1e250, 1e250).unwrap() / 1e250,
        coefficient,
        5e-16,
    );
    near(
        flux_weighted_normal_speed(1e-250, 1e-250).unwrap() / 1e-250,
        coefficient,
        5e-16,
    );
    let zero = MolecularMassFlux {
        outgoing_kg_m2_s: 0.0,
        incoming_kg_m2_s: 0.0,
        net_kg_m2_s: 0.0,
    };
    assert_eq!(
        phase_normal_momentum(zero, f64::MAX, f64::MAX, f64::MAX)
            .unwrap()
            .total_pa,
        0.0
    );
    assert_eq!(
        maxwell_reflected_normal_pressure(17.0, f64::MAX, f64::MIN_POSITIVE, 0.0, 0.0).unwrap(),
        17.0
    );
    assert_eq!(
        maxwell_reflected_normal_pressure(f64::MAX, f64::MAX, f64::MIN_POSITIVE, 1.0, 1.0).unwrap(),
        0.0
    );
    // T_s/T_g overflows; the square-root ratio and final pressure are ordinary.
    near(
        maxwell_reflected_normal_pressure(1e-300, 1e300, 1e-300, 0.0, 1.0).unwrap(),
        0.5,
        2e-16,
    );
    // Even sqrt(T_s)/sqrt(T_g) overflows, but weighted pressure restores range.
    let tiny = f64::from_bits(1);
    let actual = maxwell_reflected_normal_pressure(1e-300, 1e300, tiny, 0.0, 1.0).unwrap();
    let expected = (0.5e-300 / tiny.sqrt()) * 1e150;
    near(actual / expected, 1.0, 3e-16);
}

#[test]
fn invalid_inputs_and_unrepresentable_nonzero_moments_refuse() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(flux_weighted_normal_speed(bad, 300.0).is_err());
        assert!(flux_weighted_normal_speed(R, bad).is_err());
        assert!(maxwell_reflected_normal_pressure(0.0, bad, 300.0, 1.0, 0.0).is_err());
        assert!(maxwell_reflected_normal_pressure(0.0, 300.0, bad, 1.0, 0.0).is_err());
    }
    for bad in [-1.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(maxwell_reflected_normal_pressure(0.0, 300.0, 300.0, bad, 0.0).is_err());
        assert!(maxwell_reflected_normal_pressure(0.0, 300.0, 300.0, 1.0, bad).is_err());
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(maxwell_reflected_normal_pressure(bad, 300.0, 300.0, 0.0, 0.0).is_err());
        let beam = MolecularMassFlux {
            outgoing_kg_m2_s: bad,
            incoming_kg_m2_s: 0.0,
            net_kg_m2_s: 0.0,
        };
        assert!(phase_normal_momentum(beam, R, 300.0, 300.0).is_err());
    }
    assert!(flux_weighted_normal_speed(f64::MAX, f64::MAX).is_err());
    let tiny = f64::from_bits(1);
    let beam = MolecularMassFlux {
        outgoing_kg_m2_s: tiny,
        incoming_kg_m2_s: 0.0,
        net_kg_m2_s: tiny,
    };
    assert!(phase_normal_momentum(beam, 1e-200, 1e-200, 1e-200).is_err());
    assert!(maxwell_reflected_normal_pressure(f64::MAX, 1200.0, 300.0, 0.0, 1.0).is_err());
}

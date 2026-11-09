use super::{ClassicalBulkAbsorption, dispersion, invalid};
use crate::MaterialError;
use fs_math::c64::C64;

mod energy;
pub use energy::{BulkWaveAbsorption, BulkWaveCellAbsorption, BulkWaveField};

/// One forward homogeneous Navier--Stokes--Fourier acoustic mode, using
/// `exp(i k x - i omega t)`. This excludes the independent thermal mode and
/// transverse/wall boundary layers; it does not close thermal end conditions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkThermoviscousWave {
    wavenumber: C64,
    traction_impedance: C64,
    pressure_impedance: C64,
    thermal_denominator: C64,
}

impl BulkThermoviscousWave {
    /// Forward acoustic wavenumber [rad/m], with nonnegative attenuation.
    pub const fn wavenumber_rad_m(self) -> C64 {
        self.wavenumber
    }

    /// Total compressive normal traction / longitudinal velocity [Pa s/m].
    /// With `eta = 4 mu/3 + zeta`, traction is `t = p - eta * du/dx`;
    /// momentum gives `t/u = rho * omega/k`. Mechanical power flux is
    /// `Re(t conj(u))/2`. The full thermal acoustic-energy flux also contains
    /// `Re(T' conj(q'))/(2 T0)`; it is not inferred from this mechanical port.
    pub const fn traction_impedance_pa_s_m(self) -> C64 {
        self.traction_impedance
    }

    /// Thermodynamic pressure / longitudinal velocity [Pa s/m], distinct
    /// from traction: `p/u = rho * omega/k + i eta k`. Both traveling
    /// directions share the same pressure/traction ratio. The response keeps
    /// pressure resolved even when the two terms above nearly cancel.
    pub const fn pressure_impedance_pa_s_m(self) -> C64 {
        self.pressure_impedance
    }
}

impl ClassicalBulkAbsorption {
    /// Evaluate the acoustic mode and its two physically distinct impedances
    /// from this frozen material. The existing acoustic branch/refusals and
    /// constant-coefficient continuum limits apply; no material is queried
    /// again and neither wall nor molecular loss is added. At DC return the
    /// continuous traveling-wave limit `k=0`, `p/u=t/u=rho*c`.
    pub fn thermoviscous_wave(
        &self,
        omega_rad_s: f64,
    ) -> Result<BulkThermoviscousWave, MaterialError> {
        let acoustic = self.state.acoustic();
        let transport = self.state.transport();
        let c = acoustic.sound_speed_m_s();
        let nu = (4.0 / 3.0) * transport.kinematic_viscosity_m2_s()
            + self.bulk_viscosity_pa_s / acoustic.density_kg_m3();
        let (wavenumber, traction, pressure, thermal_denominator) = dispersion::acoustic_response(
            omega_rad_s,
            c,
            (nu / c) / c,
            (transport.thermal_diffusivity_m2_s() / c) / c,
            self.state.heat_capacity_ratio(),
        )
        .map_err(invalid)?;
        let reference = acoustic.characteristic_impedance_pa_s_m();
        let traction = traction.scale(reference);
        let pressure = pressure.scale(reference);
        if [traction, pressure]
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite() || z.abs() == 0.0)
            || traction.re <= 0.0
        {
            return Err(invalid("bulk acoustic impedances are unrepresentable"));
        }
        Ok(BulkThermoviscousWave {
            wavenumber,
            traction_impedance: traction,
            pressure_impedance: pressure,
            thermal_denominator,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_acoustics::FluidAcousticState;
    use crate::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
    use crate::fluid_transport::FluidTransportState;

    fn material(gamma: f64, mu: f64, bulk: f64) -> ClassicalBulkAbsorption {
        let state = ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(2.0, 10.0).unwrap(),
            FluidTransportState::try_new(2.0, 8.0, 0.6, mu).unwrap(),
            8.0 / gamma,
        )
        .unwrap();
        ClassicalBulkAbsorption::from_material(&state, bulk).unwrap()
    }

    fn near(actual: C64, expected: C64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
            "actual={actual:?}, expected={expected:?}"
        );
    }

    #[test]
    fn g1_bulk_wave_closes_momentum_calorics_and_total_energy_flux() {
        // Analytic ideal-gas fixture: rho=2, c=10, Cp=8, Cv=4, R=4,
        // T0=12.5, p0=100. Solve the ORIGINAL heat equation for T'/u;
        // neither the production caloric factor nor a pressure fit is used.
        let law = material(2.0, 0.003, 0.002);
        let (rho, eta, kappa, cv, r, temperature, p0) = (2.0, 0.006, 0.6, 4.0, 4.0, 12.5, 100.0);
        for omega in [0.01, 1.0, 100.0, 10_000.0] {
            let wave = law.thermoviscous_wave(omega).unwrap();
            let k = wave.wavenumber_rad_m();
            assert_eq!(k, law.thermoviscous_wavenumber(omega).unwrap());
            let traction = wave.traction_impedance_pa_s_m();
            let pressure = wave.pressure_impedance_pa_s_m();
            let ik = C64::new(-k.im, k.re);
            near(traction * k, C64::from_re(rho * omega), 3e-15);
            near(pressure - ik.scale(eta), traction, 3e-14);
            let t = -ik.scale(p0) / ((k * k).scale(kappa) + C64::new(0.0, -omega * rho * cv));
            near(
                pressure,
                k.scale(rho * r * temperature / omega) + t.scale(rho * r),
                2e-13,
            );
            let flux = 0.5 * traction.re + kappa * k.im * t.norm_sq() / (2.0 * temperature);
            let heat = 0.5 * (eta * k.norm_sq() + kappa * (k * t).norm_sq() / temperature);
            assert!((2.0 * k.im * flux / heat - 1.0).abs() < 2e-13);
        }
    }

    #[test]
    fn g1_bulk_wave_retains_pressure_in_the_strong_viscosity_limit() {
        // Cp=Cv decouples heat exactly: k=omega/sqrt(c²-i omega eta/rho).
        // Same longitudinal viscosity, independently varied shear and bulk.
        for (mu, bulk) in [(0.003, 0.002), (0.0015, 0.004)] {
            let law = material(1.0, mu, bulk);
            for x in [1e-10, 1.0, 1e8, 1e20] {
                let omega = x * 100.0 / 0.003;
                let wave = law.thermoviscous_wave(omega).unwrap();
                let speed = C64::new(100.0, -omega * 0.003).sqrt();
                near(wave.wavenumber_rad_m(), C64::from_re(omega) / speed, 4e-15);
                near(wave.traction_impedance_pa_s_m(), speed.scale(2.0), 4e-15);
                near(
                    wave.pressure_impedance_pa_s_m(),
                    C64::from_re(200.0) / speed,
                    4e-15,
                );
                let k = wave.wavenumber_rad_m();
                let shear_heat = 0.5 * (4.0 * mu / 3.0) * k.norm_sq();
                let bulk_heat = 0.5 * bulk * k.norm_sq();
                let lost_flux = k.im * wave.traction_impedance_pa_s_m().re;
                assert!((lost_flux / (shear_heat + bulk_heat) - 1.0).abs() < 4e-15);
            }
        }
    }

    #[test]
    fn g1_bulk_wave_preserves_pressure_when_diffusion_scales_coincide() {
        // nu_L=gamma*chi=0.075. At large x=omega*nu_L/c² the exact
        // quadratic gives p/u -> rho*c*sqrt((gamma-1)/gamma). Neither
        // subtracting viscous stress nor dividing by the near-zero thermal
        // denominator resolves this limit from rounded k alone.
        let law = material(2.0, 0.075, 0.05);
        let wave = law.thermoviscous_wave(1e22 / 0.075).unwrap();
        near(
            wave.pressure_impedance_pa_s_m(),
            C64::from_re(20.0 * fs_math::det::sqrt(0.5)),
            2e-10,
        );
    }

    #[test]
    fn g1_bulk_wave_dc_weak_loss_and_invalid_frequency() {
        let law = material(2.0, 0.003, 0.002);
        let dc = law.thermoviscous_wave(0.0).unwrap();
        assert_eq!(dc.wavenumber_rad_m(), C64::ZERO);
        assert_eq!(dc.pressure_impedance_pa_s_m(), C64::from_re(20.0));
        assert_eq!(dc.traction_impedance_pa_s_m(), C64::from_re(20.0));
        let omega = 1e-5;
        let wave = law.thermoviscous_wave(omega).unwrap();
        let (x, y) = (omega * 0.003 / 100.0, omega * 0.0375 / 100.0);
        assert!((wave.traction_impedance_pa_s_m().im / (-10.0 * (x + y)) - 1.0).abs() < 2e-15);
        assert!((wave.pressure_impedance_pa_s_m().im / (10.0 * (x - y)) - 1.0).abs() < 2e-15);
        for bad in [-1.0, f64::NAN, f64::INFINITY, f64::MAX, f64::from_bits(1)] {
            assert!(law.thermoviscous_wave(bad).is_err());
        }
    }
}

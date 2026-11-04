//! Classical longitudinal absorption of a uniform Newtonian fluid.
//!
//! Stokes–Kirchhoff gives the pressure-amplitude decay coefficient
//! `alpha = omega² [4 mu/3 + zeta + kappa (1/Cv - 1/Cp)] / (2 rho c³)`.
//! Density, adiabatic c, Cp/Cv and transport all belong to one original
//! [`ThermoviscousFluid`] evaluation; no ideal-gas EOS is used. See the
//! generalized relation, Eq. (2), in Kustova et al. (2023),
//! <https://doi.org/10.3390/fluids8020048>.
//!
//! This is the weak-loss, small-signal, homogeneous-continuum limit. The
//! supplied bulk viscosity is a real, frequency-independent coefficient over
//! the caller's band. It may be zero only as an explicit omission/assumption;
//! shear viscosity does not determine it. If it represents internal
//! relaxation, adding that same relaxation again would double-count loss.
//! [`ClassicalBulkAbsorption::thermoviscous_wavenumber`] also exposes the
//! homogeneous linear Navier--Stokes--Fourier acoustic branch, including its
//! coupled attenuation and phase dispersion. Neither path infers wall loss,
//! molecular spectra or scattering. Original source support and frequency
//! validity stay with the caller; this snapshot grants no source authority.

use super::{ThermoviscousFluid, ThermoviscousFluidState};
use crate::MaterialError;
use fs_math::c64::C64;

mod dispersion;

/// Prepared classical bulk absorption, with no per-frequency material queries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClassicalBulkAbsorption {
    state: ThermoviscousFluidState,
    bulk_viscosity_pa_s: f64,
    relaxation_time_s: f64,
}

impl ClassicalBulkAbsorption {
    /// Freeze one coherent original material evaluation and an explicitly
    /// declared nonnegative bulk viscosity [Pa s]. Original material refusals
    /// propagate unchanged. A successful construction does not establish the
    /// weak-loss assumption at any particular frequency; inspect
    /// [`Self::weak_loss_parameter`] over the intended band.
    pub fn from_material<M: ThermoviscousFluid + ?Sized>(
        material: &M,
        bulk_viscosity_pa_s: f64,
    ) -> Result<Self, MaterialError> {
        if !bulk_viscosity_pa_s.is_finite() || bulk_viscosity_pa_s < 0.0 {
            return Err(invalid(
                "classical bulk viscosity must be finite and nonnegative",
            ));
        }
        let state = material.thermoviscous_state()?;
        let transport = state.transport();
        let acoustic = state.acoustic();
        // Divide each transport contribution by rho before combining them.
        // Avoid forming rho*c³, which can overflow for an admitted medium
        // even when its final attenuation and relaxation time are finite.
        let diffusivity = (4.0 / 3.0) * transport.kinematic_viscosity_m2_s()
            + bulk_viscosity_pa_s / acoustic.density_kg_m3()
            + (state.heat_capacity_ratio() - 1.0) * transport.thermal_diffusivity_m2_s();
        let relaxation_time_s =
            (diffusivity / acoustic.sound_speed_m_s()) / acoustic.sound_speed_m_s();
        if !diffusivity.is_finite() || !relaxation_time_s.is_finite() || relaxation_time_s <= 0.0 {
            return Err(invalid("classical bulk absorption is unrepresentable"));
        }
        Ok(Self {
            state,
            bulk_viscosity_pa_s,
            relaxation_time_s,
        })
    }

    /// The original frozen acoustic/caloric/transport response.
    pub const fn state(&self) -> ThermoviscousFluidState {
        self.state
    }

    /// Explicit compressional viscosity [Pa s], distinct from shear viscosity.
    pub const fn bulk_viscosity_pa_s(&self) -> f64 {
        self.bulk_viscosity_pa_s
    }

    /// Longitudinal diffusion time `[4 mu/3 + zeta + kappa(1/Cv-1/Cp)]/(rho c²)`.
    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }

    /// Dimensionless `omega*tau`, which must be much smaller than one for
    /// this leading-order attenuation law. This is an applicability diagnostic,
    /// not an error bound or a fitted molecular relaxation time.
    pub fn weak_loss_parameter(&self, omega_rad_s: f64) -> Result<f64, MaterialError> {
        if !omega_rad_s.is_finite() || omega_rad_s < 0.0 {
            return Err(invalid(
                "classical absorption needs finite nonnegative angular frequency",
            ));
        }
        let value = omega_rad_s * self.relaxation_time_s;
        if !value.is_finite() {
            return Err(invalid("classical weak-loss parameter is unrepresentable"));
        }
        Ok(value)
    }

    /// Pressure-amplitude attenuation [Np/m], with exact zero at DC.
    /// Power/intensity decays as `exp(-2 alpha r)`, not `exp(-alpha r)`.
    /// Frequencies outside the weak-loss regime are not silently clipped;
    /// the caller must apply its explicit applicability threshold.
    pub fn absorption_np_per_m(&self, omega_rad_s: f64) -> Result<f64, MaterialError> {
        let parameter = self.weak_loss_parameter(omega_rad_s)?;
        let alpha = (0.5 * omega_rad_s / self.state.acoustic().sound_speed_m_s()) * parameter;
        if !alpha.is_finite() || (omega_rad_s > 0.0 && alpha == 0.0) {
            return Err(invalid(
                "classical absorption coefficient is unrepresentable",
            ));
        }
        Ok(alpha)
    }

    /// Complex acoustic wavenumber [rad/m] of the full homogeneous linear
    /// Navier--Stokes--Fourier equations at real angular frequency [rad/s].
    ///
    /// For `exp(i k r - i omega t)`, `Re(k)` is the spatial phase rate and
    /// `Im(k) >= 0` is pressure-amplitude attenuation. Intensity decays as
    /// `exp(-2 Im(k) r)`. The acoustic branch continues from `k = omega/c`
    /// at zero frequency; the separate thermal diffusion mode is not returned.
    /// This is an allocation-free evaluation of the original frozen density,
    /// adiabatic sound speed, Cp/Cv, shear/bulk viscosity and conductivity,
    /// with no ideal-gas substitution or per-frequency material query.
    ///
    /// Unlike [`Self::absorption_np_per_m`], the attenuation is not truncated
    /// to its weak-loss limit. The same original continuum, linear-response
    /// and constant-coefficient validity limits still apply. Unresolved
    /// acoustic/thermal branch crossings and unrepresentable results refuse.
    /// Molecular relaxation or wall losses are not added; a bulk viscosity
    /// already representing relaxation must not be counted again. Fourier
    /// heat conduction and Newtonian viscosity are parabolic, so this model
    /// does not establish a finite wavefront at `r/c`.
    pub fn thermoviscous_wavenumber(&self, omega_rad_s: f64) -> Result<C64, MaterialError> {
        let acoustic = self.state.acoustic();
        let transport = self.state.transport();
        let c = acoustic.sound_speed_m_s();
        let longitudinal_diffusivity = (4.0 / 3.0) * transport.kinematic_viscosity_m2_s()
            + self.bulk_viscosity_pa_s / acoustic.density_kg_m3();
        dispersion::acoustic_wavenumber(
            omega_rad_s,
            c,
            (longitudinal_diffusivity / c) / c,
            (transport.thermal_diffusivity_m2_s() / c) / c,
            self.state.heat_capacity_ratio(),
        )
        .map_err(invalid)
    }
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_acoustics::FluidAcousticState;
    use crate::fluid_transport::FluidTransportState;
    use crate::gas::{GasSpec, GasState};

    fn fluid(cv: f64) -> ThermoviscousFluidState {
        ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(2.0, 10.0).unwrap(),
            FluidTransportState::try_new(2.0, 8.0, 0.6, 0.003).unwrap(),
            cv,
        )
        .unwrap()
    }

    #[test]
    fn g1_general_fluid_matches_independent_specific_heat_formula() {
        let law = ClassicalBulkAbsorption::from_material(&fluid(4.0), 0.002).unwrap();
        let omega = 3.0;
        let expected = omega * omega / (2.0 * 2.0 * 10.0_f64.powi(3))
            * (4.0 * 0.003 / 3.0 + 0.002 + 0.6 * (1.0 / 4.0 - 1.0 / 8.0));
        assert!((law.absorption_np_per_m(omega).unwrap() / expected - 1.0).abs() < 3e-15);
        assert_eq!(law.state(), fluid(4.0));
        assert_eq!(law.bulk_viscosity_pa_s(), 0.002);
        assert_eq!(law.absorption_np_per_m(0.0).unwrap(), 0.0);
    }

    #[test]
    fn g3_bulk_and_thermal_contributions_remain_independent() {
        let shear = ClassicalBulkAbsorption::from_material(&fluid(8.0), 0.0).unwrap();
        let bulk = ClassicalBulkAbsorption::from_material(&fluid(8.0), 0.004).unwrap();
        let thermal = ClassicalBulkAbsorption::from_material(&fluid(4.0), 0.0).unwrap();
        assert!(
            (bulk.absorption_np_per_m(10.0).unwrap() / shear.absorption_np_per_m(10.0).unwrap()
                - 2.0)
                .abs()
                < 1e-15
        );
        assert!(
            thermal.absorption_np_per_m(10.0).unwrap() > bulk.absorption_np_per_m(10.0).unwrap()
        );
        assert_eq!(
            shear.absorption_np_per_m(20.0).unwrap(),
            4.0 * shear.absorption_np_per_m(10.0).unwrap()
        );
        assert_eq!(
            bulk.weak_loss_parameter(20.0).unwrap(),
            2.0 * bulk.weak_loss_parameter(10.0).unwrap()
        );
    }


    #[test]
    fn g1_full_dispersion_retains_original_state_and_classical_limit() {
        for cv in [4.0, 8.0] {
            for bulk in [0.0, 0.002] {
                let original = fluid(cv);
                let law = ClassicalBulkAbsorption::from_material(&original, bulk).unwrap();
                let omega = 0.01;
                let wave = law.thermoviscous_wavenumber(omega).unwrap();
                assert!((wave.im / law.absorption_np_per_m(omega).unwrap() - 1.0).abs() < 1e-8);
                assert!((wave.re / (omega / 10.0) - 1.0).abs() < 1e-8);
                assert_eq!(law.state(), original);
                assert_eq!(law.thermoviscous_wavenumber(0.0).unwrap(), C64::ZERO);
            }
        }
    }

    #[test]
    fn g3_full_dispersion_does_not_extend_the_quadratic_loss_approximation() {
        let shear_bulk = ClassicalBulkAbsorption::from_material(&fluid(8.0), 0.002).unwrap();
        let thermal = ClassicalBulkAbsorption::from_material(&fluid(4.0), 0.002).unwrap();
        let omega = 1000.0;
        let viscous_wave = shear_bulk.thermoviscous_wavenumber(omega).unwrap();
        let thermal_wave = thermal.thermoviscous_wavenumber(omega).unwrap();
        assert!(thermal_wave.im > viscous_wave.im);
        assert!((thermal_wave.im / thermal.absorption_np_per_m(omega).unwrap() - 1.0).abs() > 0.01);
        assert!((thermal_wave.re / (omega / 10.0) - 1.0).abs() > 1e-3);
        for bad in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(thermal.thermoviscous_wavenumber(bad).is_err());
        }
    }

    #[test]
    fn g3_gas_uses_its_original_material_and_classical_limit() {
        for spec in [GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()] {
            let gas = GasState::try_new(&spec, 400.0, 120_000.0).unwrap();
            let law = ClassicalBulkAbsorption::from_material(&gas, 0.0).unwrap();
            for omega in [1.0, 1000.0, 50_000.0] {
                assert!(
                    (law.absorption_np_per_m(omega).unwrap()
                        / gas.stokes_kirchhoff_absorption(omega)
                        - 1.0)
                        .abs()
                        < 8e-15
                );
            }
            let mut stale = gas;
            stale.density *= 2.0;
            assert!(ClassicalBulkAbsorption::from_material(&stale, 0.0).is_err());
        }
    }

    #[test]
    fn g1_iapws_liquid_steam_and_dense_water_retain_their_original_response() {
        use fs_qty::{Pressure, Temperature};
        use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
        use fs_thermochem::water::transport::If97TransportState;
        for (temperature, pressure) in [(300.0, 3e6), (700.0, 0.5e6), (650.0, 25e6)] {
            let water = If97TransportState::at_pressure_temperature(
                Pressure::new(pressure),
                Temperature::new(temperature),
                DensePressureSolve {
                    pressure_tolerance_pa: 1e-3,
                    density_tolerance_kg_m3: 1e-8,
                    maximum_iterations: 256,
                },
                || false,
            )
            .unwrap();
            // Zero is an explicit bulk-loss omission, not an IAPWS bulk claim.
            let law = ClassicalBulkAbsorption::from_material(&water, 0.0).unwrap();
            let source = water.thermodynamics();
            let rho = source.density().value();
            let c = source.sound_speed().value();
            let cp = source.cp().value();
            let cv = source.cv().value();
            let omega = 10_000.0;
            let expected = omega * omega / (2.0 * rho * c.powi(3))
                * (4.0 * water.dynamic_viscosity().value() / 3.0
                    + water.thermal_conductivity().value() * (1.0 / cv - 1.0 / cp));
            assert_eq!(law.state().acoustic().density_kg_m3(), rho);
            assert_eq!(law.state().acoustic().sound_speed_m_s(), c);
            assert!((law.absorption_np_per_m(omega).unwrap() / expected - 1.0).abs() < 1e-13);
        }
    }

    #[test]
    fn g0_invalid_bulk_frequency_and_unrepresentable_loss_refuse() {
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(ClassicalBulkAbsorption::from_material(&fluid(4.0), bad).is_err());
        }
        let law = ClassicalBulkAbsorption::from_material(&fluid(4.0), 0.0).unwrap();
        for bad in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(law.absorption_np_per_m(bad).is_err());
        }
        // The original nonzero coefficient may not disappear as zero loss.
        assert!(law.absorption_np_per_m(f64::from_bits(1)).is_err());
    }
}

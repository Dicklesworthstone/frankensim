//! Molecular-kinetic contact-line friction with explicit calibration and domain.
//!
//! Blake-Haynes molecular kinetics uses
//! U = 2 K lambda sinh[sigma lambda^2 (cos(theta_0)-cos(theta))/(2 k_B T)].
//! This implementation inverts that relation at the actual signed line speed:
//! f = 2 k_B T/lambda^2 asinh[U/(2 K lambda)],
//! cos(theta) = cos(theta_0) - f/sigma.
//!
//! The force f is per unit contact-line length [N/m]; f*U is nonnegative
//! dissipation per unit length [W/m]. It is a diagnostic of this local closure,
//! not an evolved heat source or a global surface-energy balance.
//!
//! Reference: Farrokhi et al., JFM (2025), equation (3.2),
//! <https://doi.org/10.1017/jfm.2025.10407>.
//! Jump length and equilibrium jump frequency are caller-supplied calibrated
//! model parameters, never inferred from viscosity, a material name or a static
//! contact angle. The original constructor holds both constant. The explicit
//! Arrhenius constructor instead selects K(T)=K_ref exp[(E_a/k_B)(1/T_ref-1/T)]
//! within the same admitted temperature interval. E_a is joules per molecular
//! event, not J/mol or equilibrium adhesion work. This separately calibrated
//! constant-prefactor model follows <https://doi.org/10.1351/goldbook.A00446>;
//! it does not infer an Eyring prefactor or activation energy from static data.
//! A selected static advancing/receding endpoint or explicit equilibrium angle
//! supplies theta_0. Combining a hysteresis endpoint with this rate law does not
//! prove microscopic pinning; the equilibrium mode introduces no deadband.
//!
//! O(1) arithmetic, no allocation, fixed deterministic transcendental owner.
//! Hydrodynamic matching, gas-wedge dissipation, slip, precursor films,
//! entrainment, chemical ageing and surfactant transport are not represented.

/// Exact SI Boltzmann constant [J/K].
const BOLTZMANN_J_K: f64 = 1.380_649e-23;

/// Explicit model calibration and applicability. There are no universal defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularKineticControls {
    /// Effective distance between adsorption sites [m], positive.
    pub jump_length_m: f64,
    /// Positive equilibrium jump frequency [1/s], constant for `try_new` or
    /// evaluated at the explicit reference temperature for `try_new_arrhenius`.
    pub equilibrium_jump_frequency_per_s: f64,
    /// Lower temperature supported by this calibration [K], positive.
    pub minimum_temperature_k: f64,
    /// Upper supported temperature [K], at least the lower bound.
    pub maximum_temperature_k: f64,
    /// Maximum absolute contact-line speed supported by this calibration [m/s].
    pub maximum_speed_m_s: f64,
    /// Minimum static AND applied angle [rad], at least 1e-4.
    pub minimum_angle_radians: f64,
    /// Maximum static AND applied angle [rad], at most pi-1e-4.
    pub maximum_angle_radians: f64,
}

/// Refusal leaves the boundary owner free to retain its accepted state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MolecularKineticError {
    Invalid(&'static str),
    OutsideDomain(&'static str),
    NonRepresentable,
}
impl core::fmt::Display for MolecularKineticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid molecular-kinetic input: {reason}"),
            Self::OutsideDomain(reason) => write!(f, "molecular-kinetic domain: {reason}"),
            Self::NonRepresentable => f.write_str("molecular-kinetic response is unrepresentable"),
        }
    }
}
impl std::error::Error for MolecularKineticError {}

/// The local model result, kept separate from original material measurements.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularKineticSample {
    pub applied_angle_radians: f64,
    pub temperature_k: f64,
    pub surface_tension_n_m: f64,
    /// Actual model jump frequency at this sample's temperature [1/s].
    pub equilibrium_jump_frequency_per_s: f64,
    /// Actual speed consumed by this rate law. A pinned branch consumes zero.
    pub signed_speed_m_s: f64,
    /// Signed rate-dependent excess force relative to the selected static
    /// endpoint [N/m]. A hysteresis endpoint need not be the equilibrium angle.
    pub driving_force_per_length_n_m: f64,
    /// Nonnegative f*U [W/m]; not automatically added to liquid enthalpy.
    pub dissipation_per_length_w_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct ArrheniusJumpFrequency {
    reference_temperature_k: f64,
    activation_energy_j: f64,
    activation_temperature_k: f64,
    log_reference_frequency: f64,
}

/// Admitted jump kinetics; temperature and sigma are sampled at use.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularKineticLaw {
    controls: MolecularKineticControls,
    velocity_scale_m_s: f64,
    force_scale_per_k_n_m: f64,
    arrhenius: Option<ArrheniusJumpFrequency>,
}

impl MolecularKineticLaw {
    /// Select the original constant-frequency calibration over its temperature domain.
    pub fn try_new(controls: MolecularKineticControls) -> Result<Self, MolecularKineticError> {
        if !controls.jump_length_m.is_finite()
            || controls.jump_length_m <= 0.0
            || !controls.equilibrium_jump_frequency_per_s.is_finite()
            || controls.equilibrium_jump_frequency_per_s <= 0.0
            || !controls.minimum_temperature_k.is_finite()
            || controls.minimum_temperature_k <= 0.0
            || !controls.maximum_temperature_k.is_finite()
            || controls.maximum_temperature_k < controls.minimum_temperature_k
            || !controls.maximum_speed_m_s.is_finite()
            || controls.maximum_speed_m_s <= 0.0
            || !controls.minimum_angle_radians.is_finite()
            || controls.minimum_angle_radians < 1e-4
            || !controls.maximum_angle_radians.is_finite()
            || controls.maximum_angle_radians > core::f64::consts::PI - 1e-4
            || controls.minimum_angle_radians >= controls.maximum_angle_radians
        {
            return Err(MolecularKineticError::Invalid(
                "jump kinetics, temperature, speed or angular bounds",
            ));
        }
        let velocity_scale_m_s =
            (2.0 * controls.jump_length_m) * controls.equilibrium_jump_frequency_per_s;
        let force_scale_per_k_n_m =
            (2.0 * BOLTZMANN_J_K / controls.jump_length_m) / controls.jump_length_m;
        for value in [
            velocity_scale_m_s,
            force_scale_per_k_n_m,
            force_scale_per_k_n_m * controls.minimum_temperature_k,
            force_scale_per_k_n_m * controls.maximum_temperature_k,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(MolecularKineticError::NonRepresentable);
            }
        }
        Ok(Self {
            controls,
            velocity_scale_m_s,
            force_scale_per_k_n_m,
            arrhenius: None,
        })
    }

    /// Select an explicitly calibrated, nonnegative, temperature-independent Arrhenius
    /// activation energy [J per event]. The controls' frequency is K(T_ref).
    /// Zero activation exactly recovers the constant-frequency model. Both
    /// endpoint rates and velocity scales must be representable; no clipping,
    /// extrapolation, activation-heat source or source claim is introduced.
    pub fn try_new_arrhenius(
        controls: MolecularKineticControls,
        reference_temperature_k: f64,
        activation_energy_j: f64,
    ) -> Result<Self, MolecularKineticError> {
        let mut law = Self::try_new(controls)?;
        if !reference_temperature_k.is_finite()
            || reference_temperature_k < controls.minimum_temperature_k
            || reference_temperature_k > controls.maximum_temperature_k
            || !activation_energy_j.is_finite()
            || activation_energy_j < 0.0
        {
            return Err(MolecularKineticError::Invalid("Arrhenius calibration"));
        }
        let activation_temperature_k = activation_energy_j / BOLTZMANN_J_K;
        if !activation_temperature_k.is_finite() {
            return Err(MolecularKineticError::NonRepresentable);
        }
        law.arrhenius = Some(ArrheniusJumpFrequency {
            reference_temperature_k,
            activation_energy_j,
            activation_temperature_k,
            log_reference_frequency: fs_math::det::ln(controls.equilibrium_jump_frequency_per_s),
        });
        // For nonnegative activation the frequency is monotone in temperature.
        for temperature in [
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        ] {
            law.scales_at_temperature(temperature)?;
        }
        Ok(law)
    }

    /// Optional original (reference temperature [K], activation energy [J/event]).
    #[must_use]
    pub fn arrhenius_calibration(&self) -> Option<(f64, f64)> {
        self.arrhenius
            .map(|a| (a.reference_temperature_k, a.activation_energy_j))
    }

    /// Evaluate the selected jump frequency [1/s] inside its calibrated domain.
    pub fn jump_frequency_per_s_at(
        &self,
        temperature_k: f64,
    ) -> Result<f64, MolecularKineticError> {
        if !temperature_k.is_finite() {
            return Err(MolecularKineticError::Invalid("temperature"));
        }
        if temperature_k < self.controls.minimum_temperature_k
            || temperature_k > self.controls.maximum_temperature_k
        {
            return Err(MolecularKineticError::OutsideDomain("temperature"));
        }
        let reference = self.controls.equilibrium_jump_frequency_per_s;
        let Some(a) = self.arrhenius else {
            return Ok(reference);
        };
        if temperature_k == a.reference_temperature_k || a.activation_temperature_k == 0.0 {
            return Ok(reference);
        }
        // Factored reciprocal difference retains nearby temperatures. Adding
        // logarithms avoids overflow of exp(shift) when K_ref is very small.
        let shift = (a.activation_temperature_k / temperature_k.min(a.reference_temperature_k))
            * ((temperature_k - a.reference_temperature_k)
                / temperature_k.max(a.reference_temperature_k));
        let frequency = fs_math::det::exp(a.log_reference_frequency + shift);
        if !frequency.is_finite() || frequency <= 0.0 {
            return Err(MolecularKineticError::NonRepresentable);
        }
        Ok(frequency)
    }

    fn scales_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<(f64, f64), MolecularKineticError> {
        let frequency = self.jump_frequency_per_s_at(temperature_k)?;
        let velocity_scale = if frequency == self.controls.equilibrium_jump_frequency_per_s {
            self.velocity_scale_m_s
        } else {
            (2.0 * self.controls.jump_length_m) * frequency
        };
        if !velocity_scale.is_finite() || velocity_scale <= 0.0 {
            return Err(MolecularKineticError::NonRepresentable);
        }
        Ok((frequency, velocity_scale))
    }

    #[must_use]
    pub const fn controls(&self) -> &MolecularKineticControls {
        &self.controls
    }

    pub fn admit_angle(&self, angle_radians: f64) -> Result<(), MolecularKineticError> {
        if !angle_radians.is_finite() {
            return Err(MolecularKineticError::Invalid("contact angle"));
        }
        if angle_radians < self.controls.minimum_angle_radians
            || angle_radians > self.controls.maximum_angle_radians
        {
            return Err(MolecularKineticError::OutsideDomain("contact angle"));
        }
        Ok(())
    }

    /// Invert the full nonlinear rate law. Positive speed advances liquid.
    ///
    /// Requires the actual local absolute temperature and surface tension used
    /// by the fluid boundary. Static endpoints retain their original source
    /// receipts; this prediction does not create a measured dynamic-angle claim.
    /// An unsupported angle, temperature or speed refuses without clipping.
    pub fn apparent_angle(
        &self,
        endpoint_angle_radians: f64,
        signed_speed_m_s: f64,
        temperature_k: f64,
        surface_tension_n_m: f64,
    ) -> Result<MolecularKineticSample, MolecularKineticError> {
        self.admit_angle(endpoint_angle_radians)?;
        if !signed_speed_m_s.is_finite()
            || !temperature_k.is_finite()
            || !surface_tension_n_m.is_finite()
            || surface_tension_n_m <= 0.0
        {
            return Err(MolecularKineticError::Invalid(
                "speed, temperature or surface tension",
            ));
        }
        let (equilibrium_jump_frequency_per_s, velocity_scale_m_s) =
            self.scales_at_temperature(temperature_k)?;
        if signed_speed_m_s.abs() > self.controls.maximum_speed_m_s {
            return Err(MolecularKineticError::OutsideDomain("contact-line speed"));
        }
        if signed_speed_m_s == 0.0 {
            return Ok(MolecularKineticSample {
                applied_angle_radians: endpoint_angle_radians,
                temperature_k,
                surface_tension_n_m,
                equilibrium_jump_frequency_per_s,
                signed_speed_m_s,
                driving_force_per_length_n_m: 0.0,
                dissipation_per_length_w_m: 0.0,
            });
        }
        let magnitude = signed_speed_m_s.abs();
        let ratio = magnitude / velocity_scale_m_s;
        let inverse = if ratio.is_infinite() {
            // For an overflowing positive ratio the exponentially small
            // correction to ln(2*x) is below any representable f64 precision.
            fs_math::det::ln(magnitude) - fs_math::det::ln(velocity_scale_m_s)
                + core::f64::consts::LN_2
        } else if ratio > 0.0 {
            asinh_positive(ratio)
        } else {
            return Err(MolecularKineticError::NonRepresentable);
        };
        let force_magnitude = (self.force_scale_per_k_n_m * temperature_k) * inverse;
        let driving_force_per_length_n_m = force_magnitude.copysign(signed_speed_m_s);
        let cosine_shift = driving_force_per_length_n_m / surface_tension_n_m;
        let dissipation_per_length_w_m = force_magnitude * magnitude;
        if !force_magnitude.is_finite()
            || force_magnitude <= 0.0
            || !cosine_shift.is_finite()
            || cosine_shift == 0.0
            || !dissipation_per_length_w_m.is_finite()
        {
            return Err(MolecularKineticError::NonRepresentable);
        }
        let cosine = fs_math::det::cos(endpoint_angle_radians) - cosine_shift;
        if !(-1.0..=1.0).contains(&cosine) {
            return Err(MolecularKineticError::OutsideDomain(
                "requested speed has no finite contact angle",
            ));
        }
        let applied_angle_radians = fs_math::det::acos(cosine);
        self.admit_angle(applied_angle_radians)?;
        Ok(MolecularKineticSample {
            applied_angle_radians,
            temperature_k,
            surface_tension_n_m,
            equilibrium_jump_frequency_per_s,
            signed_speed_m_s,
            driving_force_per_length_n_m,
            dissipation_per_length_w_m,
        })
    }
}

// Stable odd inverse is evaluated on |U| and the sign is applied to force.
// Avoid cancellation at small speed and squaring overflow at large speed.
fn asinh_positive(x: f64) -> f64 {
    if x < 1e-4 {
        let square = x * x;
        return x * (1.0 + square * (-1.0 / 6.0 + square * (3.0 / 40.0)));
    }
    if x > 1e150 {
        return fs_math::det::ln(x) + core::f64::consts::LN_2;
    }
    let y = x + (x / (1.0 + (1.0 + x * x).sqrt())) * x;
    let one_plus_y = 1.0 + y;
    // Correct for rounding in 1+y; deterministic ln remains the only log owner.
    fs_math::det::ln(one_plus_y) * (y / (one_plus_y - 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn controls() -> MolecularKineticControls {
        MolecularKineticControls {
            jump_length_m: 5e-10,
            equilibrium_jump_frequency_per_s: 1e8,
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 400.0,
            maximum_speed_m_s: 1.0,
            minimum_angle_radians: 0.05,
            maximum_angle_radians: 3.0,
        }
    }

    #[test]
    fn g1_arrhenius_rates_recover_forward_kinetics_and_linear_friction() {
        let energy = 2000.0 * BOLTZMANN_J_K;
        let law = MolecularKineticLaw::try_new_arrhenius(controls(), 300.0, energy).unwrap();
        assert_eq!(law.arrhenius_calibration(), Some((300.0, energy)));
        let mut previous_friction = f64::INFINITY;
        for temperature in [250.0_f64, 300.0, 400.0] {
            let frequency = 1e8 * (2000.0 * (1.0 / 300.0 - 1.0 / temperature)).exp();
            assert!(
                (law.jump_frequency_per_s_at(temperature).unwrap() / frequency - 1.0).abs() < 2e-14
            );
            let friction = BOLTZMANN_J_K * temperature / (frequency * 5e-10_f64.powi(3));
            assert!(friction < previous_friction);
            previous_friction = friction;
            for speed in [-0.003_f64, -1e-9, 0.0, 1e-9, 0.003] {
                let sample = law.apparent_angle(1.2, speed, temperature, 0.07).unwrap();
                assert!((sample.equilibrium_jump_frequency_per_s / frequency - 1.0).abs() < 2e-14);
                let exponent =
                    0.07 * 5e-10_f64.powi(2) * (1.2_f64.cos() - sample.applied_angle_radians.cos())
                        / (2.0 * BOLTZMANN_J_K * temperature);
                let recovered = 2.0 * frequency * 5e-10 * exponent.sinh();
                assert!((recovered - speed).abs() < 5e-16);
                assert_eq!(
                    sample.dissipation_per_length_w_m,
                    sample.driving_force_per_length_n_m * speed
                );
                if speed.abs() == 1e-9 {
                    assert!(
                        (sample.driving_force_per_length_n_m / speed / friction - 1.0).abs()
                            < 2e-14
                    );
                }
            }
        }
    }

    #[test]
    fn g3_arrhenius_reference_and_zero_activation_preserve_original_bits() {
        let fixed = MolecularKineticLaw::try_new(controls()).unwrap();
        let zero = MolecularKineticLaw::try_new_arrhenius(controls(), 300.0, 0.0).unwrap();
        let activated =
            MolecularKineticLaw::try_new_arrhenius(controls(), 300.0, 2000.0 * BOLTZMANN_J_K)
                .unwrap();
        assert_eq!(fixed.arrhenius_calibration(), None);
        for speed in [-0.003, 0.0, 0.003] {
            assert_eq!(
                activated.apparent_angle(1.2, speed, 300.0, 0.07),
                fixed.apparent_angle(1.2, speed, 300.0, 0.07)
            );
            for temperature in [250.0, 300.0, 400.0] {
                assert_eq!(
                    zero.apparent_angle(1.2, speed, temperature, 0.07),
                    fixed.apparent_angle(1.2, speed, temperature, 0.07)
                );
            }
        }
    }

    #[test]
    fn g0_arrhenius_admits_only_calibrated_representable_rates() {
        for (reference, energy) in [
            (249.0, 1e-20),
            (401.0, 1e-20),
            (300.0, -1e-20),
            (f64::NAN, 0.0),
            (300.0, f64::NAN),
            (300.0, 1e7 * BOLTZMANN_J_K),
            (300.0, f64::MAX),
        ] {
            assert!(MolecularKineticLaw::try_new_arrhenius(controls(), reference, energy).is_err());
        }
        let law = MolecularKineticLaw::try_new_arrhenius(controls(), 300.0, 1e-20).unwrap();
        for temperature in [0.0, 249.0, 401.0, f64::NAN, f64::INFINITY] {
            assert!(law.jump_frequency_per_s_at(temperature).is_err());
            assert!(law.apparent_angle(1.2, 0.0, temperature, 0.07).is_err());
        }
        // exp(833...) overflows, but the actual frequency is representable.
        let wide = MolecularKineticLaw::try_new_arrhenius(
            MolecularKineticControls {
                equilibrium_jump_frequency_per_s: 1e-200,
                minimum_temperature_k: 300.0,
                maximum_temperature_k: 600.0,
                ..controls()
            },
            300.0,
            5e5 * BOLTZMANN_J_K,
        )
        .unwrap();
        let expected = (1e-200_f64.ln() + 5e5 * (1.0 / 300.0 - 1.0 / 600.0)).exp();
        assert!((wide.jump_frequency_per_s_at(600.0).unwrap() / expected - 1.0).abs() < 2e-13);
    }

    #[test]
    fn g1_full_inverse_matches_independent_forward_molecular_kinetics() {
        let law = MolecularKineticLaw::try_new(controls()).unwrap();
        for speed in [-0.03_f64, -1e-4, -1e-8, 1e-8, 1e-4, 0.03] {
            let sample = law.apparent_angle(1.2, speed, 300.0, 0.07).unwrap();
            let exponent = 0.07
                * controls().jump_length_m.powi(2)
                * (1.2_f64.cos() - sample.applied_angle_radians.cos())
                / (2.0 * BOLTZMANN_J_K * 300.0);
            let recovered = 2.0
                * controls().equilibrium_jump_frequency_per_s
                * controls().jump_length_m
                * exponent.sinh();
            assert!((recovered - speed).abs() < 3e-16, "{speed}: {recovered}");
            assert!(sample.dissipation_per_length_w_m > 0.0);
            assert_eq!(
                sample.dissipation_per_length_w_m,
                sample.driving_force_per_length_n_m * speed,
            );
        }
    }

    #[test]
    fn g1_small_speed_recovers_the_independent_linear_friction_limit() {
        let law = MolecularKineticLaw::try_new(controls()).unwrap();
        let sample = law.apparent_angle(1.2, 1e-9, 300.0, 0.07).unwrap();
        let friction = BOLTZMANN_J_K * 300.0
            / (controls().equilibrium_jump_frequency_per_s * controls().jump_length_m.powi(3));
        assert!((sample.driving_force_per_length_n_m / 1e-9 / friction - 1.0).abs() < 5e-15);
        let linearized_angle_cosine = 1.2_f64.cos() - friction * 1e-9 / 0.07;
        assert!((sample.applied_angle_radians.cos() - linearized_angle_cosine).abs() < 3e-16);
    }

    #[test]
    fn g3_temperature_tension_and_calibration_change_the_actual_angle() {
        let law = MolecularKineticLaw::try_new(controls()).unwrap();
        let cold = law.apparent_angle(1.2, 0.01, 280.0, 0.07).unwrap();
        let hot = law.apparent_angle(1.2, 0.01, 350.0, 0.07).unwrap();
        let lower_tension = law.apparent_angle(1.2, 0.01, 280.0, 0.035).unwrap();
        assert!(hot.applied_angle_radians > cold.applied_angle_radians);
        assert!(lower_tension.applied_angle_radians > cold.applied_angle_radians);
        assert!(
            (hot.driving_force_per_length_n_m / cold.driving_force_per_length_n_m - 350.0 / 280.0)
                .abs()
                < 3e-16
        );
        assert_eq!(
            cold.driving_force_per_length_n_m,
            lower_tension.driving_force_per_length_n_m
        );
        let faster = MolecularKineticLaw::try_new(MolecularKineticControls {
            equilibrium_jump_frequency_per_s: 2e8,
            ..controls()
        })
        .unwrap()
        .apparent_angle(1.2, 0.01, 280.0, 0.07)
        .unwrap();
        assert!(faster.applied_angle_radians < cold.applied_angle_radians);
        assert!(faster.dissipation_per_length_w_m < cold.dissipation_per_length_w_m);
    }

    #[test]
    fn g0_zero_speed_retains_exact_endpoints_and_has_zero_loss() {
        let law = MolecularKineticLaw::try_new(controls()).unwrap();
        for angle in [
            controls().minimum_angle_radians,
            1.2,
            controls().maximum_angle_radians,
        ] {
            for speed in [0.0, -0.0] {
                let sample = law.apparent_angle(angle, speed, 300.0, 0.07).unwrap();
                assert_eq!(sample.applied_angle_radians.to_bits(), angle.to_bits());
                assert_eq!(sample.driving_force_per_length_n_m, 0.0);
                assert_eq!(sample.dissipation_per_length_w_m, 0.0);
            }
        }
    }

    #[test]
    fn g3_angle_reflection_reverses_force_and_preserves_loss() {
        let law = MolecularKineticLaw::try_new(controls()).unwrap();
        let direct = law.apparent_angle(1.2, 0.02, 300.0, 0.07).unwrap();
        let reflected = law
            .apparent_angle(core::f64::consts::PI - 1.2, -0.02, 300.0, 0.07)
            .unwrap();
        assert!(
            (direct.applied_angle_radians + reflected.applied_angle_radians
                - core::f64::consts::PI)
                .abs()
                < 5e-16
        );
        assert_eq!(
            direct.driving_force_per_length_n_m,
            -reflected.driving_force_per_length_n_m
        );
        assert_eq!(
            direct.dissipation_per_length_w_m,
            reflected.dissipation_per_length_w_m
        );
    }

    #[test]
    fn g0_unsupported_domains_refuse_instead_of_clipping_the_angle() {
        let law = MolecularKineticLaw::try_new(controls()).unwrap();
        for (angle, speed, temperature, sigma) in [
            (1.2, 1.01, 300.0, 0.07),
            (1.2, 0.0, 249.0, 0.07),
            (1.2, 0.0, 401.0, 0.07),
            (0.04, 0.0, 300.0, 0.07),
            (3.01, 0.0, 300.0, 0.07),
            (1.2, f64::NAN, 300.0, 0.07),
            (1.2, 0.01, 300.0, 0.0),
            (1.2, 1.0, 300.0, 1e-6),
        ] {
            assert!(
                law.apparent_angle(angle, speed, temperature, sigma)
                    .is_err()
            );
        }
        for invalid in [
            MolecularKineticControls {
                jump_length_m: 0.0,
                ..controls()
            },
            MolecularKineticControls {
                equilibrium_jump_frequency_per_s: -1.0,
                ..controls()
            },
            MolecularKineticControls {
                maximum_temperature_k: 200.0,
                ..controls()
            },
            MolecularKineticControls {
                minimum_angle_radians: 0.0,
                ..controls()
            },
            MolecularKineticControls {
                maximum_speed_m_s: 0.0,
                ..controls()
            },
        ] {
            assert!(MolecularKineticLaw::try_new(invalid).is_err());
        }
    }

    #[test]
    fn g1_inverse_hyperbolic_sine_handles_near_equilibrium_and_large_ratios() {
        for x in [1e-200_f64, 1e-10, 9e-5, 1e-4, 0.1, 1.0, 1e4, 1e150, 1e250] {
            let expected = x.asinh();
            assert!((asinh_positive(x) / expected - 1.0).abs() < 5e-15, "{x}");
        }
    }
}

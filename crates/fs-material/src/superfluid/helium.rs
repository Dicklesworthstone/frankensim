//! A single measured helium-II reference at 1.8 K and saturated vapor pressure.
//!
//! Source: Donnelly and Barenghi, *The Observed Properties of Liquid Helium at
//! the Saturated Vapor Pressure*, J. Phys. Chem. Ref. Data 27, 1217 (1998),
//! <https://doi.org/10.1063/1.556028>; open paper:
//! <https://srd.nist.gov/jpcrdreprint/1.556028.pdf>.
//!
//! The source's density derivative follows the saturation curve; it is not
//! the isobaric expansion needed by the two-fluid equations. Its heat capacity
//! likewise follows saturation. This owner retains those original values and
//! converts both derivatives, also distinguishing the measured fast sound mode
//! from the thermodynamic adiabatic speed. Tables are independently fitted and
//! rounded: agreement with the independent second-sound value is about 0.8%,
//! not a numerical certificate of experimental accuracy.
//!
//! This is ONE equilibrium point, not an interpolation table or a helium
//! equation of state. No other temperature, pressure, phase transition,
//! critical velocity, mutual-friction coefficient, or material uncertainty is
//! inferred. Perturbation solvers must declare their own linearization domain.

use super::{TwoFluidProperties, TwoFluidState, invalid};
use crate::MaterialError;
use fs_math::det;

/// Source of all tabulated values and the ITS-90 pressure relation below.
pub const SOURCE_DOI: &str = "10.1063/1.556028";

/// Original reported row, with its original mass units where indicated.
/// Returned by value; editing a copy never changes the admitted reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeliumSvpSourceRow {
    pub temperature_k: f64,
    /// Table 1.2, total liquid density [g/cm³].
    pub density_g_cm3: f64,
    /// Table 2.3, superfluid density [g/cm³]. Normal density is obtained by
    /// subtracting this from the more precisely tabulated total density.
    pub superfluid_density_g_cm3: f64,
    /// Table 8.5, entropy from calorimetry [J/(g K)].
    pub entropy_j_g_k: f64,
    /// Table 7.4, heat capacity along saturation [J/(mol K)].
    pub saturation_heat_capacity_j_mol_k: f64,
    /// Table 1.2 / equation 1.6: `-(1/rho) d rho_sat/dT` [1/K].
    pub saturation_expansion_per_k: f64,
    /// Table 3.3, measured fast two-fluid sound mode [m/s].
    pub fast_sound_m_s: f64,
    /// Table 4.3, independent second-sound observation [m/s].
    pub second_sound_m_s: f64,
    /// Table 6.3, independent fourth-sound observation [m/s].
    pub fourth_sound_m_s: f64,
    /// Table 12.3, normal dynamic shear viscosity [Pa s].
    pub normal_viscosity_pa_s: f64,
    /// Table 17.3, rounded saturation pressure [Pa].
    pub reported_saturation_pressure_pa: f64,
}

const ROW: HeliumSvpSourceRow = HeliumSvpSourceRow {
    temperature_k: 1.8,
    density_g_cm3: 0.145_353_8,
    superfluid_density_g_cm3: 0.099_82,
    entropy_j_g_k: 0.539_9,
    saturation_heat_capacity_j_mol_k: 11.81,
    saturation_expansion_per_k: -0.006_79,
    fast_sound_m_s: 231.3,
    second_sound_m_s: 19.90,
    fourth_sound_m_s: 192.6,
    normal_viscosity_pa_s: 1.298e-6,
    reported_saturation_pressure_pa: 1638.0,
};

// Section 1 / section 7 unit conversion; the source uses 4.0026 g/mol.
const MOLAR_MASS_KG_MOL: f64 = 0.004_002_6;

/// Immutable admitted material and pressure at exactly the retained source row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeliumSvpReference {
    state: TwoFluidState,
    saturation_pressure_pa: f64,
    saturation_pressure_derivative_pa_k: f64,
}

impl HeliumSvpReference {
    /// Construct the 1.8 K SVP point; this API deliberately has no temperature
    /// or pressure argument and offers no interpolation or extrapolation.
    pub fn at_1_8_k() -> Result<Self, MaterialError> {
        let (pressure, slope) = svp_at_1_8_k()?;
        let t = ROW.temperature_k;
        let rho = 1000.0 * ROW.density_g_cm3;
        let rho_s = 1000.0 * ROW.superfluid_density_g_cm3;
        let s = 1000.0 * ROW.entropy_j_g_k;
        let c_sigma = ROW.saturation_heat_capacity_j_mol_k / MOLAR_MASS_KG_MOL;
        let x = ROW.fast_sound_m_s * ROW.fast_sound_m_s;
        let d_entropy = (rho_s / (rho - rho_s)) * t * s * s;
        let h = t * slope / rho;
        let a0 = ROW.saturation_expansion_per_k + slope / (rho * x);
        let q = c_sigma - d_entropy / x;
        let denominator = q - h * a0;
        if ![q, denominator].into_iter().all(super::positive) {
            return Err(invalid("helium SVP thermodynamic conversion is singular"));
        }
        // Solve these identities together, using the fast Landau root:
        // alpha_p = alpha_sigma + kappa_T * p'_sat;
        // cp = c_sigma + (T/rho) alpha_p * p'_sat (source section 18);
        // kappa_T = 1/(rho*c_ad²) + T*alpha_p²/(rho*cp).
        // Eliminating cp and c_ad² gives alpha_p = a0*q/(q-h*a0).
        let alpha = a0 * q / denominator;
        let cp = c_sigma + h * alpha;
        let b = d_entropy / cp;
        let expansion = t * alpha * alpha / cp;
        let sound_denominator = x - b + x * b * expansion;
        if ![cp, x - b, sound_denominator]
            .into_iter()
            .all(super::positive)
        {
            return Err(invalid("helium SVP fast-mode conversion is singular"));
        }
        let c_ad2 = x * (x - b) / sound_denominator;
        let state = TwoFluidState::try_new(TwoFluidProperties {
            temperature_k: t,
            density_kg_m3: rho,
            superfluid_fraction: rho_s / rho,
            entropy_j_kg_k: s,
            specific_heat_cp_j_kg_k: cp,
            adiabatic_sound_speed_m_s: c_ad2.sqrt(),
            thermal_expansion_per_k: alpha,
            normal_viscosity_pa_s: ROW.normal_viscosity_pa_s,
        })?;
        Ok(Self {
            state,
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: slope,
        })
    }

    pub const fn state(&self) -> &TwoFluidState {
        &self.state
    }
    pub const fn raw(&self) -> HeliumSvpSourceRow {
        ROW
    }
    /// ITS-90 evaluation; more digits than the independently rounded table row.
    pub const fn saturation_pressure_pa(&self) -> f64 {
        self.saturation_pressure_pa
    }
    pub const fn saturation_pressure_derivative_pa_k(&self) -> f64 {
        self.saturation_pressure_derivative_pa_k
    }
}

// Table 17.1: T = sum A_i z^i, z = (ln(p/Pa)-5.6)/2.9.
// Only used to recover the value and derivative at the one admitted point.
fn its90_temperature_and_derivative(z: f64) -> (f64, f64) {
    const A: [f64; 9] = [
        1.392_408, 0.527_153, 0.166_756, 0.050_988, 0.026_514, 0.001_975, -0.017_976, 0.005_409,
        0.013_259,
    ];
    let mut temperature = A[8];
    let mut derivative = 0.0;
    for coefficient in A[..8].iter().rev() {
        derivative = derivative * z + temperature;
        temperature = temperature * z + coefficient;
    }
    (temperature, derivative)
}

fn svp_at_1_8_k() -> Result<(f64, f64), MaterialError> {
    let mut z = (det::ln(ROW.reported_saturation_pressure_pa) - 5.6) / 2.9;
    for _ in 0..8 {
        let (temperature, derivative) = its90_temperature_and_derivative(z);
        if !super::positive(derivative) {
            return Err(invalid("helium SVP ITS-90 inversion is singular"));
        }
        z -= (temperature - ROW.temperature_k) / derivative;
    }
    let (temperature, derivative) = its90_temperature_and_derivative(z);
    let pressure = det::exp(2.9 * z + 5.6);
    let slope = 2.9 * pressure / derivative;
    if ![pressure, slope].into_iter().all(super::positive)
        || (temperature - ROW.temperature_k).abs() > 1e-13
    {
        return Err(invalid("helium SVP ITS-90 inversion did not converge"));
    }
    Ok((pressure, slope))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g2_pressure_agrees_with_source_and_its90_derivative() {
        let reference = HeliumSvpReference::at_1_8_k().unwrap();
        let pressure = reference.saturation_pressure_pa();
        assert!((pressure - reference.raw().reported_saturation_pressure_pa).abs() < 0.5);
        let temperature_at = |p| its90_temperature_and_derivative((det::ln(p) - 5.6) / 2.9).0;
        assert!((temperature_at(pressure) - 1.8).abs() < 1e-13);
        let dp_dt = 0.02 / (temperature_at(pressure + 0.01) - temperature_at(pressure - 0.01));
        assert!((dp_dt / reference.saturation_pressure_derivative_pa_k() - 1.0).abs() < 1e-8);
    }

    #[test]
    fn g0_isobaric_properties_close_both_saturation_chain_rules() {
        let reference = HeliumSvpReference::at_1_8_k().unwrap();
        let state = reference.state();
        let row = reference.raw();
        let slope = reference.saturation_pressure_derivative_pa_k();
        let alpha_sigma =
            state.thermal_expansion_per_k() - state.isothermal_compressibility_per_pa() * slope;
        let c_sigma = state.specific_heat_cp_j_kg_k()
            - state.temperature_k() * state.thermal_expansion_per_k() * slope
                / state.density_kg_m3();
        assert!((alpha_sigma - row.saturation_expansion_per_k).abs() < 1e-15);
        assert!((c_sigma * MOLAR_MASS_KG_MOL - row.saturation_heat_capacity_j_mol_k).abs() < 1e-12);
        assert!((state.thermal_expansion_per_k() - row.saturation_expansion_per_k).abs() > 0.0007);
        assert!(state.specific_heat_cp_j_kg_k() < c_sigma);
    }

    #[test]
    fn g2_fast_mode_retains_source_and_second_sound_predicts_independent_measurement() {
        let reference = HeliumSvpReference::at_1_8_k().unwrap();
        let modes = reference.state().sound_modes();
        assert!((modes.fast_m_s - reference.raw().fast_sound_m_s).abs() < 1e-10);
        // Independent fitted data are not mutually exact thermodynamic inputs.
        assert!((modes.slow_m_s / reference.raw().second_sound_m_s - 1.0).abs() < 0.01);
        assert_ne!(
            reference.state().adiabatic_sound_speed_m_s(),
            modes.fast_m_s
        );
    }

    #[test]
    fn g2_clamped_normal_mode_predicts_independent_fourth_sound() {
        let reference = HeliumSvpReference::at_1_8_k().unwrap();
        let state = reference.state();
        let ct2 = state.isothermal_sound_speed_m_s().powi(2);
        let entropy_term = state.entropy_j_kg_k() - state.thermal_expansion_per_k() * ct2;
        // v_n=0 fixes entropy per volume: c4²=rho_s*(d mu/d rho)_(rho*s).
        let fourth = (state.superfluid_fraction()
            * (ct2
                + state.temperature_k() / state.specific_heat_cv_j_kg_k() * entropy_term.powi(2)))
        .sqrt();
        assert!((fourth / reference.raw().fourth_sound_m_s - 1.0).abs() < 0.005);
        assert_eq!(
            state.normal_viscosity_pa_s(),
            reference.raw().normal_viscosity_pa_s
        );
    }
}

//! Discrete measured helium-II references at saturated vapor pressure.
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
//! rounded: the independent second-sound discrepancy is about 2.8% at 1.4 K,
//! below 0.8% at 1.6/1.8 K, and below 0.1% at 2.0 K. These are consistency
//! checks between source fits, not certificates of experimental accuracy.
//!
//! Only the four retained equilibrium temperatures are admitted. There is no
//! interpolation, extrapolation or helium equation of state. No other
//! temperature, pressure, phase transition, critical velocity, mutual-friction
//! coefficient, or material uncertainty is inferred. Perturbation solvers must
//! declare their own linearization domain.

use super::{TwoFluidProperties, TwoFluidState, invalid};
use crate::MaterialError;
use fs_math::det;

/// Source of all tabulated values and the ITS-90 pressure relation below.
pub const SOURCE_DOI: &str = "10.1063/1.556028";

/// Discrete ITS-90 temperatures retained from the source tables [K].
/// Selection is exact: nearby values are not rounded to a reference point.
pub const SUPPORTED_TEMPERATURES_K: [f64; 4] = [1.4, 1.6, 1.8, 2.0];

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

// Printed pages 1224 (table 1.2), 1227 (2.3), 1230 (3.3), 1233 (4.3),
// 1235 (6.3), 1241 (7.4), 1245 (8.5), 1257 (12.3), and 1268 (17.3).
// The 10^3 factor in table 1.2's expansion column is removed here; density,
// entropy and heat capacity retain the original source units until admission.
const ROWS: [HeliumSvpSourceRow; 4] = [
    HeliumSvpSourceRow {
        temperature_k: 1.4,
        density_g_cm3: 0.145_137_3,
        superfluid_density_g_cm3: 0.134_57,
        entropy_j_g_k: 0.132_3,
        saturation_heat_capacity_j_mol_k: 3.092,
        saturation_expansion_per_k: -0.001_45,
        fast_sound_m_s: 236.3,
        second_sound_m_s: 19.56,
        fourth_sound_m_s: 227.8,
        normal_viscosity_pa_s: 1.416e-6,
        reported_saturation_pressure_pa: 281.9,
    },
    HeliumSvpSourceRow {
        temperature_k: 1.6,
        density_g_cm3: 0.145_207_1,
        superfluid_density_g_cm3: 0.121_63,
        entropy_j_g_k: 0.282_4,
        saturation_heat_capacity_j_mol_k: 6.285,
        saturation_expansion_per_k: -0.003_53,
        fast_sound_m_s: 234.5,
        second_sound_m_s: 20.33,
        fourth_sound_m_s: 214.7,
        normal_viscosity_pa_s: 1.306e-6,
        reported_saturation_pressure_pa: 746.5,
    },
    HeliumSvpSourceRow {
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
    },
    HeliumSvpSourceRow {
        temperature_k: 2.0,
        density_g_cm3: 0.145_621_7,
        superfluid_density_g_cm3: 0.065_07,
        entropy_j_g_k: 0.962_1,
        saturation_heat_capacity_j_mol_k: 21.28,
        saturation_expansion_per_k: -0.012_2,
        fast_sound_m_s: 226.2,
        second_sound_m_s: 16.78,
        fourth_sound_m_s: 152.3,
        normal_viscosity_pa_s: 1.468e-6,
        reported_saturation_pressure_pa: 3130.0,
    },
];

// Section 1 / section 7 unit conversion; the source uses 4.0026 g/mol.
const MOLAR_MASS_KG_MOL: f64 = 0.004_002_6;

/// Immutable admitted material and pressure at exactly the retained source row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeliumSvpReference {
    source_row: HeliumSvpSourceRow,
    state: TwoFluidState,
    saturation_pressure_pa: f64,
    saturation_pressure_derivative_pa_k: f64,
}

impl HeliumSvpReference {
    /// Construct the original 1.8 K SVP point with unchanged source values
    /// and thermodynamic conversion.
    pub fn at_1_8_k() -> Result<Self, MaterialError> {
        Self::at_temperature_k(1.8)
    }

    /// Construct one of [`SUPPORTED_TEMPERATURES_K`] at its saturated vapor
    /// pressure. Every other temperature, including values between the
    /// retained points, refuses. Each result is a separate frozen equilibrium
    /// reference; selecting a point does not evolve an existing channel's
    /// temperature or provide a caloric law between points.
    pub fn at_temperature_k(temperature_k: f64) -> Result<Self, MaterialError> {
        let row = ROWS
            .into_iter()
            .find(|row| row.temperature_k == temperature_k)
            .ok_or_else(|| {
                invalid("helium SVP reference requires exactly 1.4, 1.6, 1.8 or 2.0 K")
            })?;
        let (pressure, slope) = svp_at_row(row)?;
        let t = row.temperature_k;
        let rho = 1000.0 * row.density_g_cm3;
        let rho_s = 1000.0 * row.superfluid_density_g_cm3;
        let s = 1000.0 * row.entropy_j_g_k;
        let c_sigma = row.saturation_heat_capacity_j_mol_k / MOLAR_MASS_KG_MOL;
        let x = row.fast_sound_m_s * row.fast_sound_m_s;
        let d_entropy = (rho_s / (rho - rho_s)) * t * s * s;
        let h = t * slope / rho;
        let a0 = row.saturation_expansion_per_k + slope / (rho * x);
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
            normal_viscosity_pa_s: row.normal_viscosity_pa_s,
        })?;
        Ok(Self {
            source_row: row,
            state,
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: slope,
        })
    }

    pub const fn state(&self) -> &TwoFluidState {
        &self.state
    }
    pub const fn raw(&self) -> HeliumSvpSourceRow {
        self.source_row
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
// The retained points all use the table's 0.65 K to 2.1768 K coefficient set.
// Only used to recover the value and derivative at those discrete points.
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

fn svp_at_row(row: HeliumSvpSourceRow) -> Result<(f64, f64), MaterialError> {
    let mut z = (det::ln(row.reported_saturation_pressure_pa) - 5.6) / 2.9;
    for _ in 0..8 {
        let (temperature, derivative) = its90_temperature_and_derivative(z);
        if !super::positive(derivative) {
            return Err(invalid("helium SVP ITS-90 inversion is singular"));
        }
        z -= (temperature - row.temperature_k) / derivative;
    }
    let (temperature, derivative) = its90_temperature_and_derivative(z);
    let pressure = det::exp(2.9 * z + 5.6);
    let slope = 2.9 * pressure / derivative;
    if ![pressure, slope].into_iter().all(super::positive)
        || (temperature - row.temperature_k).abs() > 1e-13
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
        for (temperature, rounding_half_unit_pa) in SUPPORTED_TEMPERATURES_K
            .into_iter()
            .zip([0.05, 0.05, 0.5, 0.5])
        {
            let reference = HeliumSvpReference::at_temperature_k(temperature).unwrap();
            let pressure = reference.saturation_pressure_pa();
            assert!(
                (pressure - reference.raw().reported_saturation_pressure_pa).abs()
                    < rounding_half_unit_pa
            );
            let temperature_at = |p| its90_temperature_and_derivative((det::ln(p) - 5.6) / 2.9).0;
            assert!((temperature_at(pressure) - temperature).abs() < 1e-13);
            let dp_dt = 0.02 / (temperature_at(pressure + 0.01) - temperature_at(pressure - 0.01));
            assert!((dp_dt / reference.saturation_pressure_derivative_pa_k() - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn g0_isobaric_properties_close_both_saturation_chain_rules() {
        for temperature in SUPPORTED_TEMPERATURES_K {
            let reference = HeliumSvpReference::at_temperature_k(temperature).unwrap();
            let state = reference.state();
            let row = reference.raw();
            let slope = reference.saturation_pressure_derivative_pa_k();
            let alpha_sigma =
                state.thermal_expansion_per_k() - state.isothermal_compressibility_per_pa() * slope;
            let c_sigma = state.specific_heat_cp_j_kg_k()
                - state.temperature_k() * state.thermal_expansion_per_k() * slope
                    / state.density_kg_m3();
            assert!((alpha_sigma - row.saturation_expansion_per_k).abs() < 1e-15);
            assert!(
                (c_sigma * MOLAR_MASS_KG_MOL - row.saturation_heat_capacity_j_mol_k).abs() < 1e-12
            );
            assert_ne!(
                state.thermal_expansion_per_k(),
                row.saturation_expansion_per_k
            );
            assert!(state.specific_heat_cp_j_kg_k() < c_sigma);
        }
    }

    #[test]
    fn g2_fast_mode_retains_source_and_second_sound_predicts_independent_measurement() {
        // Independent fits are not mutually exact thermodynamic inputs. At
        // 1.4 K these unmodified source inputs predict 20.0989 m/s compared
        // with table 4.3's 19.56 m/s (2.76%); do not hide it by fitting s or rs.
        for (temperature, relative_discrepancy) in SUPPORTED_TEMPERATURES_K
            .into_iter()
            .zip([0.03, 0.01, 0.01, 0.001])
        {
            let reference = HeliumSvpReference::at_temperature_k(temperature).unwrap();
            let modes = reference.state().sound_modes();
            assert!((modes.fast_m_s - reference.raw().fast_sound_m_s).abs() < 1e-10);
            assert!(
                (modes.slow_m_s / reference.raw().second_sound_m_s - 1.0).abs()
                    < relative_discrepancy
            );
            assert_ne!(
                reference.state().adiabatic_sound_speed_m_s(),
                modes.fast_m_s
            );
        }
    }

    #[test]
    fn g2_clamped_normal_mode_predicts_independent_fourth_sound() {
        for temperature in SUPPORTED_TEMPERATURES_K {
            let reference = HeliumSvpReference::at_temperature_k(temperature).unwrap();
            let state = reference.state();
            let ct2 = state.isothermal_sound_speed_m_s().powi(2);
            let entropy_term = state.entropy_j_kg_k() - state.thermal_expansion_per_k() * ct2;
            // v_n=0 fixes entropy per volume: c4²=rho_s*(d mu/d rho)_(rho*s).
            let fourth = (state.superfluid_fraction()
                * (ct2
                    + state.temperature_k() / state.specific_heat_cv_j_kg_k()
                        * entropy_term.powi(2)))
            .sqrt();
            assert!((fourth / reference.raw().fourth_sound_m_s - 1.0).abs() < 0.005);
            assert_eq!(
                state.normal_viscosity_pa_s(),
                reference.raw().normal_viscosity_pa_s
            );
        }
    }

    #[test]
    fn g0_discrete_selection_retains_the_selected_source_and_original_constructor() {
        for temperature in SUPPORTED_TEMPERATURES_K {
            let reference = HeliumSvpReference::at_temperature_k(temperature).unwrap();
            assert_eq!(reference.raw().temperature_k, temperature);
            assert_eq!(reference.state().temperature_k(), temperature);
            assert_eq!(
                reference.state().density_kg_m3(),
                1000.0 * reference.raw().density_g_cm3
            );
        }
        assert_eq!(
            HeliumSvpReference::at_1_8_k().unwrap(),
            HeliumSvpReference::at_temperature_k(1.8).unwrap()
        );
    }

    #[test]
    fn g0_unretained_temperature_or_phase_refuses_instead_of_interpolating() {
        for temperature in [
            -1.0,
            0.0,
            1.5,
            1.7,
            1.9,
            2.1,
            2.1768,
            f64::from_bits(1.8_f64.to_bits() - 1),
            f64::from_bits(1.8_f64.to_bits() + 1),
            f64::NAN,
            f64::INFINITY,
        ] {
            assert!(HeliumSvpReference::at_temperature_k(temperature).is_err());
        }
    }
}

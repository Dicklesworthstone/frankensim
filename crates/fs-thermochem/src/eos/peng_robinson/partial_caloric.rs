//! Source-bound partial molar calorics of one homogeneous PR mixture.
//!
//! Species-selective energy streams carry the upstream partial molar enthalpy,
//! not a share of the bulk mixture departure. These derivatives hold T and p
//! fixed while adding one species; the supplied density identifies the original
//! homogeneous state. Its volume response follows the same pressure equation.
//! No composition perturbation, phase selection or numerical differentiation is
//! performed. The caller adds each original species' ideal-reference enthalpy.

use super::{
    Density, MAX_PENG_ROBINSON_COMPONENTS_V1, MixingRecord, PengRobinsonDensityCaloricEvaluationV1,
    PengRobinsonErrorV1, PengRobinsonMixtureV1, Temperature,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

/// One component's fixed-T,p response, in its original mixture environment.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PengRobinsonPartialMolarCaloricV1 {
    residual_enthalpy_j_mol: f64,
    molar_volume_m3_mol: f64,
}
impl PengRobinsonPartialMolarCaloricV1 {
    /// Partial H^R [J/mol], relative to that species' ideal h_i(T).
    /// Add the original ideal enthalpy once; flow work is already included.
    #[must_use]
    pub const fn residual_enthalpy_j_mol(&self) -> f64 {
        self.residual_enthalpy_j_mol
    }
    /// (dV/dn_i) at fixed T,p and all other species amounts [m³/mol].
    /// Partial volumes can be negative and are not clipped or phase volumes.
    #[must_use]
    pub const fn molar_volume_m3_mol(&self) -> f64 {
        self.molar_volume_m3_mol
    }
}

/// Scalar state and partial calorics on the borrowed original component axis.
/// Storage is bounded to 64 records; evaluation and copying allocate nothing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PengRobinsonPartialCaloricEvaluationV1<'a> {
    mixture: &'a PengRobinsonMixtureV1,
    state: PengRobinsonDensityCaloricEvaluationV1,
    components: [PengRobinsonPartialMolarCaloricV1; MAX_PENG_ROBINSON_COMPONENTS_V1],
}
impl<'a> PengRobinsonPartialCaloricEvaluationV1<'a> {
    /// Original source constants, pair interactions and ordered mole fractions.
    #[must_use]
    pub const fn mixture(&self) -> &'a PengRobinsonMixtureV1 {
        self.mixture
    }
    /// The identical state returned by the ordinary density caloric evaluator.
    #[must_use]
    pub const fn state(&self) -> &PengRobinsonDensityCaloricEvaluationV1 {
        &self.state
    }
    /// One response for each original component, including exact-zero fractions.
    #[must_use]
    pub fn components(&self) -> &[PengRobinsonPartialMolarCaloricV1] {
        &self.components[..self.mixture.components.len()]
    }
}

#[derive(Clone, Copy, Default)]
struct Row {
    covolume: f64,
    attraction: f64,
    temperature_derivative: f64,
}

impl PengRobinsonMixtureV1 {
    /// Analytic partial enthalpies and volumes at an explicit homogeneous state.
    ///
    /// Uses the same original pair values and temperature derivatives in one
    /// mixing traversal. The shared density constructor enforces positive p,
    /// free volume and local mechanical stability. Complete Cv, species ideal
    /// calorics, source/transport domains and phase stability remain with their
    /// existing owners. In particular this does not admit a two-phase mixture.
    ///
    /// Euler identities are sum(x_i hbar_i^R)=H^R and sum(x_i vbar_i)=v;
    /// no post-hoc redistribution is used to force their floating-point sums.
    pub fn evaluate_density_partial_calorics(
        &self,
        temperature: Temperature,
        density: Density,
    ) -> Result<PengRobinsonPartialCaloricEvaluationV1<'_>, PengRobinsonErrorV1> {
        let rho = density.value();
        if !rho.is_finite() || rho <= 0.0 {
            return Err(PengRobinsonErrorV1::InvalidDensity {
                bits: rho.to_bits(),
            });
        }
        let mut rows = [Row::default(); MAX_PENG_ROBINSON_COMPONENTS_V1];
        let mixed = self.mixing_caloric(temperature.value(), |record| match record {
            MixingRecord::Component {
                index, covolume, ..
            } => {
                rows[index].covolume = covolume;
            }
            MixingRecord::Pair {
                row,
                column,
                attraction,
                derivative,
                ..
            } => {
                let x = self.mole_fractions[column];
                rows[row].attraction += x * attraction;
                rows[row].temperature_derivative += x * derivative;
            }
        })?;
        let state = mixed.density_caloric_evaluation(temperature, density, self.molar_mass)?;
        let t = temperature.value();
        let p = state.pressure().value();
        let v = state.molar_volume_m3_per_mol();
        let b = mixed.b_mix;
        let denominator = v * v + 2.0 * b * v - b * b;
        let pressure_b =
            R * t / ((v - b) * (v - b)) + 2.0 * mixed.a_mix * (v - b) / (denominator * denominator);
        let pressure_v = state.pressure_derivatives().1;
        let integral = state.caloric.residual_geometry.attraction_integral;
        let integral_b = integral_b_scaled(b / v, b * integral) / (v * v);
        let energy_coefficient = t * mixed.da_mix_dt - mixed.a_mix;
        let bulk_enthalpy = state.residual_enthalpy_j_mol();
        let enthalpy_volume_derivative = p - energy_coefficient / denominator;
        if ![
            pressure_b,
            integral_b,
            energy_coefficient,
            enthalpy_volume_derivative,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(PengRobinsonErrorV1::Unrepresentable {
                what: "partial molar caloric derivatives",
            });
        }
        let mut components =
            [PengRobinsonPartialMolarCaloricV1::default(); MAX_PENG_ROBINSON_COMPONENTS_V1];
        for (i, row) in rows[..self.components.len()].iter().enumerate() {
            // n*d(x_j)/d(n_i)=delta_ij-x_j. The quadratic mixing law
            // therefore gives da_i=2*(sum_j x_j a_ij-a), while db_i=b_i-b.
            let delta_b = row.covolume - b;
            let delta_a = 2.0 * (row.attraction - mixed.a_mix);
            let delta_da_dt = 2.0 * (row.temperature_derivative - mixed.da_mix_dt);
            // Keeping p fixed determines n*dv/dn_i=vbar_i-v. Retain this
            // correction BEFORE adding v, including in the dilute limit.
            let delta_v = (delta_a / denominator - pressure_b * delta_b) / pressure_v;
            let molar_volume_m3_mol = v + delta_v;
            // d(n H^R)/dn_i, with H^R=(T a'-a)I+p*v-R*T.
            let terms = [
                bulk_enthalpy,
                (t * delta_da_dt - delta_a) * integral,
                energy_coefficient * integral_b * delta_b,
                enthalpy_volume_derivative * delta_v,
            ];
            let mut residual_enthalpy_j_mol = 0.0;
            let mut correction = 0.0;
            for term in terms {
                let adjusted = term - correction;
                let next = residual_enthalpy_j_mol + adjusted;
                correction = (next - residual_enthalpy_j_mol) - adjusted;
                residual_enthalpy_j_mol = next;
            }
            if !molar_volume_m3_mol.is_finite() || !residual_enthalpy_j_mol.is_finite() {
                return Err(PengRobinsonErrorV1::Unrepresentable {
                    what: "partial molar enthalpy/volume",
                });
            }
            components[i] = PengRobinsonPartialMolarCaloricV1 {
                residual_enthalpy_j_mol,
                molar_volume_m3_mol,
            };
        }
        Ok(PengRobinsonPartialCaloricEvaluationV1 {
            mixture: self,
            state,
            components,
        })
    }
}

// v²*(dI/db)_v = [eta/(1+2eta-eta²)-b*I]/eta². Direct subtraction
// cancels in dilute states. Integrating the original rational derivative gives
// sum_{k>=1} [k/(k+1)] c_k eta^(k-1), c_0=1,c_1=-2,c_k=-2c_{k-1}+c_{k-2}.
// At eta<=1/32 the 18 retained terms leave <3e-20 absolute truncation error.
// This is an analytic evaluation of the original logarithm, not an EOS fit.
fn integral_b_scaled(eta: f64, scaled_integral: f64) -> f64 {
    if eta <= 1.0 / 32.0 {
        let mut previous = 1.0;
        let mut coefficient = -2.0;
        let mut power = 1.0;
        let mut sum = 0.0;
        for k in 1..=18 {
            sum += (k as f64 / (k + 1) as f64) * coefficient * power;
            let next = -2.0 * coefficient + previous;
            previous = coefficient;
            coefficient = next;
            power *= eta;
        }
        sum
    } else {
        (eta / (1.0 + 2.0 * eta - eta * eta) - scaled_integral) / (eta * eta)
    }
}

#[cfg(test)]
mod tests;

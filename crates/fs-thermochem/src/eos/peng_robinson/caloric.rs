//! Finite residual calorics without subtracting large endpoint departures.
//!
//! The reference borrows the original mixture and fixes both conserved density
//! coordinates. One original mixing traversal supplies each trial state and
//! the attraction differences; density integrals are prepared once. These are
//! homogeneous endpoint changes, not heat, boundary work or latent phase heat.

use super::{
    Density, MAX_PENG_ROBINSON_COMPONENTS_V1, MixingRecord, MolarMass,
    PengRobinsonDensityCaloricEvaluationV1, PengRobinsonErrorV1, PengRobinsonMixtureV1,
    Temperature, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, det, log_one_plus,
};

#[derive(Clone, Copy, Debug, Default)]
struct ComponentReference {
    // sqrt(a0) * [(1+kappa) - kappa*sqrt(T/Tc)], BEFORE its absolute value.
    constant: f64,
    root_coefficient: f64,
    attraction: f64,
    logarithmic_derivative: f64,
    negative_alpha: bool,
}

/// Original fixed-composition reference and a prescribed final density.
/// Preparation and trial evaluation allocate nothing. Source identity is kept
/// by borrowing the mixture; caller-supplied unrelated endpoint states cannot
/// be combined. At most 64 component records are retained on the stack.
#[derive(Debug)]
pub struct PengRobinsonCaloricReferenceV1<'a> {
    mixture: &'a PengRobinsonMixtureV1,
    initial: PengRobinsonDensityCaloricEvaluationV1,
    final_density: Density,
    components: [ComponentReference; MAX_PENG_ROBINSON_COMPONENTS_V1],
    density: DensityChanges,
}

/// One admitted scalar endpoint and finite molar changes in its original PR
/// departures. Ideal-reference calorics must be added by their source owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PengRobinsonCaloricChangeV1 {
    state: PengRobinsonDensityCaloricEvaluationV1,
    initial_temperature: Temperature,
    log_density_ratio: f64,
    internal_energy_j_mol: f64,
    enthalpy_j_mol: f64,
    entropy_at_density_j_mol_k: f64,
}

impl PengRobinsonCaloricChangeV1 {
    #[must_use]
    pub const fn state(&self) -> &PengRobinsonDensityCaloricEvaluationV1 {
        &self.state
    }
    /// Delta U^R [J/mol], relative to ideal gas at the same T and composition.
    #[must_use]
    pub const fn residual_internal_energy_change_j_mol(&self) -> f64 {
        self.internal_energy_j_mol
    }
    /// Delta H^R [J/mol], including the original residual flow work once.
    #[must_use]
    pub const fn residual_enthalpy_change_j_mol(&self) -> f64 {
        self.enthalpy_j_mol
    }
    /// Delta S^R [J/(mol K)] relative to ideal gas at each endpoint's T/rho.
    #[must_use]
    pub const fn residual_entropy_at_density_change_j_mol_k(&self) -> f64 {
        self.entropy_at_density_j_mol_k
    }
    /// Add this to an ideal-reference entropy change evaluated at UNCHANGED
    /// pressure: -R*[ln(T1/T0)+ln(rho1/rho0)]. Together with the density
    /// departure this gives total entropy, without cancellation of R ln Z
    /// against entropy evaluated at separately rounded real-gas pressures.
    #[must_use]
    pub fn ideal_entropy_density_correction_j_mol_k(&self) -> f64 {
        -R * (log_ratio(
            self.state.temperature().value(),
            self.initial_temperature.value(),
        ) + self.log_density_ratio)
    }
}

impl PengRobinsonMixtureV1 {
    /// Prepare stable finite calorics from an admitted initial homogeneous
    /// state to a prescribed final density. The final temperature remains free.
    /// Every trial retains the ordinary positive-pressure, covolume and local
    /// mechanical-stability checks. Complete Cv and material/source support
    /// remain with the ideal-reference and material owners. No phase is selected.
    pub fn prepare_density_caloric_change(
        &self,
        initial_temperature: Temperature,
        initial_density: Density,
        final_density: Density,
    ) -> Result<PengRobinsonCaloricReferenceV1<'_>, PengRobinsonErrorV1> {
        let rho = initial_density.value();
        if !rho.is_finite() || rho <= 0.0 {
            return Err(PengRobinsonErrorV1::InvalidDensity {
                bits: rho.to_bits(),
            });
        }
        let mut components = [ComponentReference::default(); MAX_PENG_ROBINSON_COMPONENTS_V1];
        let mixed = self.mixing_caloric(initial_temperature.value(), |record| {
            if let MixingRecord::Component {
                index,
                attraction,
                derivative,
                alpha_base,
                a0,
                kappa,
                ..
            } = record
            {
                let scale = det::sqrt(a0);
                components[index] = ComponentReference {
                    constant: scale * (1.0 + kappa),
                    root_coefficient: scale * kappa
                        / det::sqrt(self.components[index].critical_temperature.value()),
                    attraction,
                    logarithmic_derivative: derivative / attraction,
                    negative_alpha: alpha_base.is_sign_negative(),
                };
            }
        })?;
        let initial = mixed.density_caloric_evaluation(
            initial_temperature,
            initial_density,
            self.molar_mass,
        )?;
        if components[..self.components.len()].iter().any(|c| {
            !c.constant.is_finite()
                || !c.root_coefficient.is_finite()
                || !c.logarithmic_derivative.is_finite()
        }) {
            return Err(unrepresentable("caloric increment component coefficients"));
        }
        let density = DensityChanges::new(&initial, final_density, self.molar_mass)?;
        Ok(PengRobinsonCaloricReferenceV1 {
            mixture: self,
            initial,
            final_density,
            components,
            density,
        })
    }
}

impl PengRobinsonCaloricReferenceV1<'_> {
    #[must_use]
    pub const fn initial_state(&self) -> &PengRobinsonDensityCaloricEvaluationV1 {
        &self.initial
    }
    #[must_use]
    pub const fn final_density(&self) -> Density {
        self.final_density
    }

    /// Fused original scalar state and finite residual change. The pair loop
    /// has the same fixed bounds as ordinary evaluation, with no second pair
    /// traversal, nonlinear solve, fitted law, allocation or source clone.
    pub fn evaluate(
        &self,
        temperature: Temperature,
    ) -> Result<PengRobinsonCaloricChangeV1, PengRobinsonErrorV1> {
        let t0 = self.initial.temperature().value();
        let t1 = temperature.value();
        let root0 = det::sqrt(t0);
        let root1 = det::sqrt(t1);
        let root_sum = root0 + root1;
        let root_change = (t1 - t0) / root_sum;
        let root_product = root0 * root1;
        let mut negative = [false; MAX_PENG_ROBINSON_COMPONENTS_V1];
        let mut attraction = Sum::default();
        let mut energy_factor = Sum::default();
        let mut derivative = Sum::default();
        let mixed = self.mixture.mixing_caloric(t1, |record| match record {
            MixingRecord::Component {
                index, alpha_base, ..
            } => {
                negative[index] = alpha_base.is_sign_negative();
            }
            MixingRecord::Pair {
                row,
                column,
                attraction: a1,
                derivative: da1,
                weight,
            } => {
                let left = self.components[row];
                let right = self.components[column];
                let before = left.negative_alpha != right.negative_alpha;
                let after = negative[row] != negative[column];
                let interaction = 1.0
                    - self.mixture.binary_interactions
                        [row * self.mixture.components.len() + column];
                let (delta_a, delta_energy, delta_derivative) = if before == after {
                    // EXACT branch-local algebra of sqrt(a_i*a_j), not a
                    // global polynomial fit across absolute-alpha sign changes.
                    let sign = if before { -1.0 } else { 1.0 };
                    let half_b = (0.5 * left.constant * right.root_coefficient
                        + 0.5 * right.constant * left.root_coefficient)
                        * interaction
                        * sign;
                    let c = left.root_coefficient * right.root_coefficient * interaction * sign;
                    let delta_energy = half_b * root_change;
                    (
                        root_change * (c * root_sum - 2.0 * half_b),
                        delta_energy,
                        delta_energy / root_product,
                    )
                } else {
                    // The original absolute-alpha branch changes its derivative.
                    // Keep both ORIGINAL endpoint pair values and their jump.
                    let a0 = det::sqrt(left.attraction * right.attraction) * interaction;
                    let da0 =
                        a0 * (0.5 * (left.logarithmic_derivative + right.logarithmic_derivative));
                    (a1 - a0, (t1 * da1 - a1) - (t0 * da0 - a0), da1 - da0)
                };
                attraction.add(weight * delta_a);
                energy_factor.add(weight * delta_energy);
                derivative.add(weight * delta_derivative);
            }
        })?;
        let state = mixed.density_caloric_evaluation(
            temperature,
            self.final_density,
            self.mixture.molar_mass,
        )?;
        let initial = self.initial.caloric;
        let final_state = state.caloric;
        let mean_integral = mean(
            initial.residual_geometry.attraction_integral,
            final_state.residual_geometry.attraction_integral,
        );
        let mean_energy = mean(
            t0 * initial.parameters.da_mix_dt - initial.parameters.a_mix,
            t1 * final_state.parameters.da_mix_dt - final_state.parameters.a_mix,
        );
        let internal_energy_j_mol = sum(&[
            mean_integral * energy_factor.value(),
            mean_energy * self.density.integral,
        ]);
        let enthalpy_j_mol = sum(&[
            internal_energy_j_mol,
            R * self.density.mean_excluded_volume * (t1 - t0),
            R * mean(t0, t1) * self.density.excluded_volume,
            -self.density.mean_attractive_volume * attraction.value(),
            -mean(initial.parameters.a_mix, final_state.parameters.a_mix)
                * self.density.attractive_volume,
        ]);
        let entropy_at_density_j_mol_k = sum(&[
            R * self.density.log_free_volume,
            mean_integral * derivative.value(),
            mean(
                initial.parameters.da_mix_dt,
                final_state.parameters.da_mix_dt,
            ) * self.density.integral,
        ]);
        if [
            internal_energy_j_mol,
            enthalpy_j_mol,
            entropy_at_density_j_mol_k,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(unrepresentable("finite residual caloric changes"));
        }
        Ok(PengRobinsonCaloricChangeV1 {
            state,
            initial_temperature: self.initial.temperature(),
            log_density_ratio: self.density.log_density_ratio,
            internal_energy_j_mol,
            enthalpy_j_mol,
            entropy_at_density_j_mol_k,
        })
    }
}

#[derive(Clone, Copy, Debug)]
struct DensityChanges {
    integral: f64,
    excluded_volume: f64,
    attractive_volume: f64,
    log_free_volume: f64,
    mean_excluded_volume: f64,
    mean_attractive_volume: f64,
    log_density_ratio: f64,
}

impl DensityChanges {
    fn new(
        initial: &PengRobinsonDensityCaloricEvaluationV1,
        final_density: Density,
        molar_mass: MolarMass,
    ) -> Result<Self, PengRobinsonErrorV1> {
        let rho0 = initial.density().value();
        let rho1 = final_density.value();
        if !rho1.is_finite() || rho1 <= 0.0 {
            return Err(PengRobinsonErrorV1::InvalidDensity {
                bits: rho1.to_bits(),
            });
        }
        let m = molar_mass.value();
        let b = initial.caloric.parameters.b_mix;
        let v1 = m / rho1;
        if !v1.is_finite() || v1 <= b {
            return Err(PengRobinsonErrorV1::OutsideCovolume);
        }
        // Canonical density ordering retains exact endpoint reversal, including
        // the logarithms, while all small differences come from input rho.
        let (lo, hi, sign) = if rho1 >= rho0 {
            (rho0, rho1, 1.0)
        } else {
            (rho1, rho0, -1.0)
        };
        let eta0 = b / (m / lo);
        let eta1 = b / (m / hi);
        let free0 = 1.0 - eta0;
        let free1 = 1.0 - eta1;
        let d0 = 1.0 + 2.0 * eta0 - eta0 * eta0;
        let d1 = 1.0 + 2.0 * eta1 - eta1 * eta1;
        let dr = (hi - lo) / m;
        let deta = b * dr;
        let denominator = 1.0 + (eta0 + eta1) - eta0 * eta1;
        let q = det::sqrt(2.0) * deta / denominator;
        let atanh_over_q = if q == 0.0 {
            1.0
        } else {
            (log_one_plus(q) - log_one_plus(-q)) / (2.0 * q)
        };
        let offset = -deta / free0;
        let log_free_volume = if offset > -1.0 {
            log_one_plus(offset)
        } else {
            det::ln(free1) - det::ln(free0)
        };
        let value = Self {
            integral: sign * (dr / denominator) * atanh_over_q,
            excluded_volume: sign * deta / (free0 * free1),
            attractive_volume: sign * dr * (1.0 + eta0 * eta1) / (d0 * d1),
            log_free_volume: sign * log_free_volume,
            mean_excluded_volume: mean(eta0 / free0, eta1 / free1),
            mean_attractive_volume: mean((lo / m) / d0, (hi / m) / d1),
            log_density_ratio: log_ratio(rho1, rho0),
        };
        if [
            value.integral,
            value.excluded_volume,
            value.attractive_volume,
            value.log_free_volume,
            value.mean_excluded_volume,
            value.mean_attractive_volume,
            value.log_density_ratio,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(unrepresentable("finite caloric density geometry"));
        }
        Ok(value)
    }
}

fn mean(a: f64, b: f64) -> f64 {
    0.5 * a + 0.5 * b
}
fn unrepresentable(what: &'static str) -> PengRobinsonErrorV1 {
    PengRobinsonErrorV1::Unrepresentable { what }
}

fn log_ratio(after: f64, before: f64) -> f64 {
    if after == before {
        return 0.0;
    }
    let (lo, hi, sign) = if after > before {
        (before, after, 1.0)
    } else {
        (after, before, -1.0)
    };
    let offset = (hi - lo) / lo;
    sign * if offset.is_finite() {
        log_one_plus(offset)
    } else {
        det::ln(hi) - det::ln(lo)
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let increment = value - self.correction;
        let next = self.value + increment;
        self.correction = (next - self.value) - increment;
        self.value = next;
    }
    fn value(&self) -> f64 {
        self.value
    }
}
fn sum(values: &[f64]) -> f64 {
    let mut sum = Sum::default();
    for &value in values {
        sum.add(value);
    }
    sum.value()
}

#[cfg(test)]
mod tests;

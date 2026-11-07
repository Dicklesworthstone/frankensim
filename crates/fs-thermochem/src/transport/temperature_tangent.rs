//! Temperature response of the original, already evaluated Wilke mixture.

use super::{
    MAX_WILKE_GAS_COMPONENTS_V1, WilkeAllocationStageV1, WilkeGasViscosityErrorV1,
    WilkeGasViscosityEvaluationV1, interaction_terms, reserve_vec,
};

/// One analytic derivative and the original mixing-denominator derivatives.
/// Fractions and molar masses remain fixed; pure-species derivatives are inputs.
#[derive(Clone, Debug, PartialEq)]
pub struct WilkeGasViscosityTemperatureTangentV1 {
    /// d mu_mix / d T [Pa s/K].
    pub dynamic_viscosity_derivative_pa_s_k: f64,
    denominator_derivatives_per_k: Vec<f64>,
}

impl WilkeGasViscosityTemperatureTangentV1 {
    /// d(sum_j x_j Phi_ij)/dT [1/K], in the original canonical species order.
    /// WMS conductivity consumers reuse these same temperature responses.
    pub fn denominator_derivatives_per_k(&self) -> &[f64] {
        &self.denominator_derivatives_per_k
    }
}

/// Refusals while differentiating one unchanged original Wilke evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum WilkeTemperatureTangentErrorV1 {
    /// There must be one derivative for every original active species.
    Length {
        /// Original canonical component count.
        components: usize,
        /// Offered derivative count.
        derivatives: usize,
    },
    /// The original interaction or bounded allocation refused.
    Mixing(WilkeGasViscosityErrorV1),
    /// An input or signed derivative operation is not representable.
    Arithmetic {
        /// Canonical outer species index.
        component: usize,
        /// First unrepresentable derivative operation.
        quantity: &'static str,
    },
}
impl core::fmt::Display for WilkeTemperatureTangentErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Length {
                components,
                derivatives,
            } => write!(
                f,
                "Wilke tangent needs {components} canonical derivatives, received {derivatives}",
            ),
            Self::Mixing(error) => error.fmt(f),
            Self::Arithmetic {
                component,
                quantity,
            } => write!(
                f,
                "Wilke component {component} has an unrepresentable {quantity}",
            ),
        }
    }
}
impl core::error::Error for WilkeTemperatureTangentErrorV1 {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Mixing(error) => Some(error),
            _ => None,
        }
    }
}
impl From<WilkeGasViscosityErrorV1> for WilkeTemperatureTangentErrorV1 {
    fn from(error: WilkeGasViscosityErrorV1) -> Self {
        Self::Mixing(error)
    }
}

impl WilkeGasViscosityEvaluationV1 {
    /// Differentiate this ORIGINAL mixture using canonical pure dmu/dT [Pa s/K].
    ///
    /// The exact pure values, fractions, masses and denominators remain in the
    /// retained evaluation. This uses its original interaction owner and O(N²)
    /// fixed-order arithmetic, with O(N) storage; no pair matrix or perturbed
    /// value evaluation is built. Signed and exact-zero derivatives are allowed.
    /// The caller must bind derivatives to these same pure-species sources and
    /// temperature. A numerical tangent does not establish physical source
    /// support, dilute-gas validity or an uncertainty bound.
    pub fn temperature_tangent(
        &self,
        pure_viscosity_derivatives_pa_s_k: &[f64],
    ) -> Result<WilkeGasViscosityTemperatureTangentV1, WilkeTemperatureTangentErrorV1> {
        let receipt = self.receipt();
        let count = receipt.species().len();
        if pure_viscosity_derivatives_pa_s_k.len() != count {
            return Err(WilkeTemperatureTangentErrorV1::Length {
                components: count,
                derivatives: pure_viscosity_derivatives_pa_s_k.len(),
            });
        }
        if count == 1 {
            let derivative = finite(
                0,
                "pure temperature derivative",
                pure_viscosity_derivatives_pa_s_k[0],
            )?;
            let mut denominators = reserve_vec(1, WilkeAllocationStageV1::DenominatorReceipt)?;
            denominators.push(0.0);
            return Ok(WilkeGasViscosityTemperatureTangentV1 {
                dynamic_viscosity_derivative_pa_s_k: derivative,
                denominator_derivatives_per_k: denominators,
            });
        }
        let mut logarithmic = [0.0; MAX_WILKE_GAS_COMPONENTS_V1];
        for (i, &derivative) in pure_viscosity_derivatives_pa_s_k.iter().enumerate() {
            finite(i, "pure temperature derivative", derivative)?;
            let viscosity = f64::from_bits(receipt.pure_viscosity_bits()[i]);
            logarithmic[i] = finite(
                i,
                "logarithmic temperature derivative",
                derivative / viscosity,
            )?;
            if derivative != 0.0 && logarithmic[i] == 0.0 {
                return Err(WilkeTemperatureTangentErrorV1::Arithmetic {
                    component: i,
                    quantity: "vanished logarithmic temperature derivative",
                });
            }
        }
        let mut denominator_derivatives =
            reserve_vec(count, WilkeAllocationStageV1::DenominatorReceipt)?;
        let mut mixture_derivative = 0.0;
        for i in 0..count {
            let viscosity_i = f64::from_bits(receipt.pure_viscosity_bits()[i]);
            let mass_i = f64::from_bits(receipt.molar_mass_bits()[i]);
            let mut derivative = 0.0;
            for j in 0..count {
                // A component's self interaction is constant. This also makes
                // the exact pure-species limit independent of derivative size.
                if i == j {
                    continue;
                }
                let (phi, z) = interaction_terms(
                    viscosity_i,
                    mass_i,
                    f64::from_bits(receipt.pure_viscosity_bits()[j]),
                    f64::from_bits(receipt.molar_mass_bits()[j]),
                    i,
                    j,
                )?;
                // Phi=(1+z)^2/B, z'/z=(mu_i'/mu_i-mu_j'/mu_j)/2.
                let response = finite(
                    i,
                    "interaction logarithmic derivative",
                    (z / (1.0 + z)) * (logarithmic[i] - logarithmic[j]),
                )?;
                let term = finite(
                    i,
                    "weighted interaction derivative",
                    (f64::from_bits(receipt.mole_fraction_bits()[j]) * phi) * response,
                )?;
                derivative = finite(i, "denominator temperature derivative", derivative + term)?;
            }
            denominator_derivatives.push(derivative);
            let denominator = f64::from_bits(receipt.denominator_bits()[i]);
            let contribution =
                f64::from_bits(receipt.mole_fraction_bits()[i]) * (viscosity_i / denominator);
            let response = finite(
                i,
                "component temperature derivative",
                f64::from_bits(receipt.mole_fraction_bits()[i])
                    * (pure_viscosity_derivatives_pa_s_k[i] / denominator)
                    - contribution * (derivative / denominator),
            )?;
            mixture_derivative = finite(
                i,
                "mixture temperature derivative",
                mixture_derivative + response,
            )?;
        }
        Ok(WilkeGasViscosityTemperatureTangentV1 {
            dynamic_viscosity_derivative_pa_s_k: mixture_derivative,
            denominator_derivatives_per_k: denominator_derivatives,
        })
    }
}

fn finite(
    component: usize,
    quantity: &'static str,
    value: f64,
) -> Result<f64, WilkeTemperatureTangentErrorV1> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WilkeTemperatureTangentErrorV1::Arithmetic {
            component,
            quantity,
        })
    }
}

#[cfg(test)]
mod tests;

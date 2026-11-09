use super::*;
use crate::solution::caloric::transport::{divided_enthalpy, reciprocal_difference};
use crate::solution::multicomponent::stability::MulticomponentStabilityControlV1;
use fs_math::dd::Dd;

/// One original homogeneous composition/temperature, in canonical species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentLiquidTransportPointV1<'a> {
    pub temperature: Temperature,
    pub mole_fractions: &'a [f64],
}

/// Source-owned finite-face forces and energy carriage in the molar frame.
///
/// With beta=1/T and zero total molar flux, these satisfy
/// `Delta(mu_i/T) = chemical_effort_i + h_star_i*Delta(beta) + common_gauge`.
/// The common species-independent gauge does no work in that frame. Thus
/// energy carriage `sum_i h_star_i*J_i` removes thermal chemical forcing;
/// an authored Soret/Dufour model, if present, must be added reciprocally by
/// the spatial transport owner. Unequal T at uniform composition gives zero
/// chemical effort without assuming equal pure-component heat capacities.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidTransportEffortsV1 {
    pub chemical_efforts_j_per_mol_k: Vec<f64>,
    pub partial_enthalpies_j_per_mol: Vec<f64>,
    /// h_star_i-h_star_0 formed before rounding the absolute enthalpies. Use
    /// these in the zero-net-molar energy flux so a large common datum cannot
    /// erase finite original excess-enthalpy or sensible-enthalpy contrasts.
    pub partial_enthalpy_contrasts_j_per_mol: Vec<f64>,
    pub mean_mole_fractions: Vec<f64>,
}

impl MulticomponentLiquidCaloricModelV1 {
    /// Evaluate the original A+B/T potential and pure Cp/T divided differences.
    /// Uses arithmetic face composition and reciprocal-temperature averaging
    /// of chi. All present species remain independent; no dependent projection
    /// or logarithm of an absent species is introduced.
    pub fn transport_efforts(
        &self,
        points: [MulticomponentLiquidTransportPointV1<'_>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidTransportEffortsV1> {
        poll(&mut cancelled)?;
        for point in points {
            self.admit_temperature(point.temperature)?;
            positive_support(point.mole_fractions)?;
        }
        let left = self.thermodynamics().evaluate(
            points[0].temperature,
            points[0].mole_fractions,
            &mut cancelled,
        )?;
        let right = self.thermodynamics().evaluate(
            points[1].temperature,
            points[1].mole_fractions,
            &mut cancelled,
        )?;
        let x = [left.mole_fractions(), right.mole_fractions()];
        let n = x[0].len();
        let mean: Vec<_> = (0..n).map(|i| 0.5 * (x[0][i] + x[1][i])).collect();
        let temperatures = points.map(|p| p.temperature);
        let sampled =
            Temperature::new(0.5 * temperatures[0].value() + 0.5 * temperatures[1].value());
        let middle = self
            .thermodynamics()
            .evaluate(sampled, &mean, &mut cancelled)?;
        let delta_beta = reciprocal_difference(temperatures[0].value(), temperatures[1].value())?;
        let average_beta = 0.5 / temperatures[0].value() + 0.5 / temperatures[1].value();
        // The mean of endpoint pair energies differs from the energy at the
        // midpoint by one quarter of the pair quadratic in Delta(x). Retain
        // this common partial-enthalpy term even though net molar flux is zero.
        let mut excess_correction = Sum::default();
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in i + 1..n {
                let coefficient =
                    0.25 * R * self.thermodynamics().chi_inverse_temperature_k()[i * n + j];
                excess_correction.add((coefficient * (x[1][i] - x[0][i])) * (x[1][j] - x[0][j]));
            }
        }
        let excess_correction = checked(excess_correction.value(), "divided excess enthalpy")?;
        let mut chemical_efforts_j_per_mol_k = Vec::with_capacity(n);
        let mut partial_enthalpies_j_per_mol = Vec::with_capacity(n);
        let mut partial_enthalpy_contrasts_j_per_mol = Vec::with_capacity(n);
        let mut reference_enthalpy = Dd::from_f64(0.0);
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut effort = Sum::default();
            effort.add(ideal_jump(x[0], x[1], i));
            for j in 0..n {
                let chi = self.thermodynamics().chi_constants()[i * n + j]
                    + self.thermodynamics().chi_inverse_temperature_k()[i * n + j] * average_beta;
                effort.add(chi * (x[1][j] - x[0][j]));
            }
            chemical_efforts_j_per_mol_k
                .push(checked(R * effort.value(), "finite-face chemical effort")?);
            let pure = divided_enthalpy(
                &self.heat_capacities()[i],
                temperatures[0],
                temperatures[1],
                sampled,
                delta_beta,
            )?;
            let h = pure + Dd::from_f64(middle.partial_molar_mixing_enthalpies_j_per_mol()[i])
                - Dd::from_f64(excess_correction);
            if !h.hi.is_finite() || !h.lo.is_finite() {
                return Err(MulticomponentCaloricErrorV1::Arithmetic(
                    "finite-face transported partial enthalpy",
                ));
            }
            partial_enthalpies_j_per_mol.push(checked(
                h.to_f64(),
                "finite-face transported partial enthalpy",
            )?);
            if i == 0 {
                reference_enthalpy = h;
            }
            partial_enthalpy_contrasts_j_per_mol.push(checked(
                (h - reference_enthalpy).to_f64(),
                "finite-face transported partial enthalpy contrast",
            )?);
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentLiquidTransportEffortsV1 {
            chemical_efforts_j_per_mol_k,
            partial_enthalpies_j_per_mol,
            partial_enthalpy_contrasts_j_per_mol,
            mean_mole_fractions: middle.mole_fractions().to_vec(),
        })
    }

    /// Sufficient global convexity margin on the complete temperature bracket.
    ///
    /// For sum(v)=0, the ideal Hessian obeys sum(v_i²/x_i)>=2||v||².
    /// At uniform composition the original stability matrix is I+P chi P/N;
    /// its bounded eigensolve therefore supplies a composition-independent
    /// lower curvature bound. Since chi=A+B/T is affine in 1/T, both endpoint
    /// bounds imply the same bound throughout the interval. This is a
    /// sufficient numerical admission test, not a phase-equilibrium claim or
    /// an interval-arithmetic certificate. Work and scratch are bounded by the
    /// original maximum species count and 64 eigensolver sweeps per endpoint.
    pub fn admit_transport_temperature_interval(
        &self,
        bracket: [Temperature; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        poll(&mut cancelled)?;
        for t in bracket {
            self.admit_temperature(t)?;
        }
        if bracket[0].value() >= bracket[1].value() {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "increasing transport temperature bracket required",
            ));
        }
        let n = self.heat_capacities().len();
        if n == 1 {
            poll(&mut cancelled)?;
            return Ok(2.0);
        }
        let uniform = vec![1.0 / n as f64; n];
        let mut lower = f64::INFINITY;
        for t in bracket {
            let state = self
                .thermodynamics()
                .evaluate(t, &uniform, &mut cancelled)?;
            let stability = state
                .local_stability(
                    MulticomponentStabilityControlV1 {
                        max_sweeps: 64,
                        ..MulticomponentStabilityControlV1::default()
                    },
                    &mut cancelled,
                )
                .map_err(|error| match error {
                    super::super::super::stability::MulticomponentStabilityErrorV1::Cancelled => {
                        MulticomponentCaloricErrorV1::Cancelled
                    }
                    error => MulticomponentCaloricErrorV1::Stability(error),
                })?;
            let eigenvalue =
                stability
                    .smallest_eigenvalue()
                    .ok_or(MulticomponentCaloricErrorV1::Arithmetic(
                        "missing composition curvature",
                    ))?;
            let margin = checked(
                2.0 + n as f64 * (eigenvalue - 1.0 - stability.numerical_uncertainty()),
                "global composition curvature bound",
            )?;
            lower = lower.min(margin);
        }
        if lower <= 64.0 * f64::EPSILON * n as f64 {
            return Err(MulticomponentCaloricErrorV1::TransportConvexity { lower_bound: lower });
        }
        poll(&mut cancelled)?;
        Ok(lower)
    }
}

// Same abundant-carrier convention as the original potential. A retained
// trace can change while the dominant stored fraction remains rounded to one.
fn ideal_jump(left: &[f64], right: &[f64], species: usize) -> f64 {
    if left[species] <= 0.5 && right[species] <= 0.5 {
        return log_ratio(right[species], left[species]);
    }
    if left[species] <= 0.5 || right[species] <= 0.5 {
        return original_log(right, species) - original_log(left, species);
    }
    let mut left_other = Sum::default();
    let mut right_other = Sum::default();
    let mut difference = Sum::default();
    for i in 0..left.len() {
        if i != species {
            left_other.add(left[i]);
            right_other.add(right[i]);
            difference.add(left[i] - right[i]);
        }
    }
    let delta = difference.value();
    let relative = delta.abs() / (1.0 - left_other.value().max(right_other.value()));
    let shifted = 1.0 + relative;
    let result = if shifted == 1.0 {
        relative
    } else {
        ln(shifted) * relative / (shifted - 1.0)
    };
    result.copysign(delta)
}

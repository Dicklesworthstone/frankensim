//! Source-derived affine tensor free energy and material entropy.
//!
//! At fixed total and viscous strain, psi = (e:Ceq(T):e + sum q:Cj(T):q)/2.
//! Affine tensors give S = -psi_T and U = psi + T*S, with zero fixed-strain
//! mechanical heat-capacity contribution. Density, geometry and eigenstrain
//! remain fixed. Independent affine K/G can change the derived Poisson ratio;
//! other equilibrium representations require a constant Poisson ratio.
//! Source interpolation is admitted over every
//! original knot; merely querying two endpoints does not establish this model.

use super::*;
use crate::Tangent6;
use crate::state_point::{affine_temperature_slope, pinned_temperature_line};
use fs_qty::{
    Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::collections::BTreeMap;

/// Immutable derivatives of the original affine stored tensors, expressed in
/// the consuming strain frame. Constructed only by the source-law owner.
#[derive(Clone, Debug)]
pub struct TensorTemperatureDerivative {
    pub(crate) equilibrium_pa_k: Tangent6,
    pub(crate) branches_pa_k: Vec<Tangent6>,
    pub(crate) interval_k: [f64; 2],
}

/// Original source selection and complete bounded affine thermodynamic law.
#[derive(Clone, Debug)]
pub struct AffineTensorThermal {
    source: ResolvedTensorMaterial,
    temperature_axis: String,
    derivative: TensorTemperatureDerivative,
}
impl AffineTensorThermal {
    /// Admit independent affine K(T)/G(T), E(T) with fixed isotropic nu,
    /// or a common affine scale of isotropic E/G or all orthotropic E/G with
    /// fixed Poisson ratios. Each of the 21 coefficients
    /// of every Maxwell branch may have its OWN affine slope, including signed
    /// cross terms. Relaxation times retain positive piecewise-linear support.
    ///
    /// Every source axis/validity restriction remains owned. Changing singular
    /// branch factors retains the native explicit history-transfer refusal.
    /// This law supplies neither thermal expansion nor a finite-strain model.
    pub fn try_new(
        source: &ResolvedTensorMaterial,
        card: &MaterialCard,
        temperature_axis: &str,
        interval_k: [f64; 2],
    ) -> Result<Self> {
        let invalid = |s: &str| TensorSourceError::Invalid(s.into());
        let [low, high] = interval_k;
        let t0 = source
            .state
            .state_coordinate(temperature_axis)
            .ok_or_else(|| invalid("affine tensor thermal coordinate is absent"))?;
        let quantity = source.state.axis_quantities().get(temperature_axis);
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || !(low..=high).contains(&t0)
            || quantity.is_some_and(|q| {
                *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS)
            })
            || quantity.is_none() && temperature_axis != "T"
        {
            return Err(invalid(
                "affine tensor thermodynamics requires a bounded absolute-temperature interval containing the original state",
            ));
        }
        source.admit_kelvin_thermal_interval(card, temperature_axis, interval_k)?;
        let query = |t| {
            source
                .state
                .query()
                .with(temperature_axis, t)
                .map_err(|e| invalid(&e.to_string()))
        };
        let lower = query(low)?;
        let upper = query(high)?;
        // Also authenticates the original card and selected claims.
        let endpoints = [
            source.at_query_point(card, &lower)?,
            source.at_query_point(card, &upper)?,
        ];
        let mut slopes = BTreeMap::new();
        for property in source.state.properties() {
            if !property.has_only_fixed_external_conditions(Some(temperature_axis), true) {
                return Err(invalid(
                    "affine tensor thermodynamics cannot freeze unowned source coordinates",
                ));
            }
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or_else(|| invalid("affine tensor source claim is missing"))?;
            let path = pinned_temperature_line(
                card,
                &claim.key,
                temperature_axis,
                &lower,
                &upper,
                property.answer().receipt.selected,
            )?;
            let name = property.requirement().name();
            let relaxation = source
                .properties
                .branches
                .iter()
                .any(|b| b.relaxation_time == name);
            if relaxation {
                if path.iter().any(|(_, v)| !v.is_finite() || *v <= 0.0) {
                    return Err(invalid(
                        "relaxation must stay positive at every original thermal knot",
                    ));
                }
                continue;
            }
            // Signed cross coefficients may vanish at the reference state. Use
            // an actual path endpoint for the shared affine admission and then
            // independently match the exact original source sample.
            let first = *path
                .first()
                .ok_or_else(|| invalid("empty affine source path"))?;
            let slope = affine_temperature_slope(&path, first.0, first.1).map_err(invalid)?;
            let expected = slope.mul_add(t0 - first.0, first.1);
            let scale = path.iter().fold(
                expected.abs().max(property.value_si().abs()),
                |m, (_, v)| m.max(v.abs()),
            );
            if !expected.is_finite()
                || (expected - property.value_si()).abs() > 128.0 * f64::EPSILON * scale
            {
                return Err(invalid(
                    "affine stored tensor path differs from its original resolved sample",
                ));
            }
            slopes.insert(name.to_owned(), slope);
        }
        let fixed = |name: &str| slopes.get(name).is_some_and(|s| *s == 0.0);
        if !fixed("density") || !fixed(&source.properties.strain_limit) {
            return Err(invalid(
                "affine tensor heat feedback requires fixed reference density and strain validity",
            ));
        }
        let equilibrium_pa_k = match source.properties.elasticity {
            TensorElasticity::IsotropicBulkShear => {
                // The isotropic tensor is affine in K and G independently:
                // normal diagonal K+4G/3, off-diagonal K-2G/3, Mandel shear 2G.
                // The derived E and nu need not be affine or proportional.
                let bulk_slope = slopes["bulk_modulus"];
                let shear_slope = slopes["shear_modulus"];
                let diagonal = (4.0 / 3.0_f64).mul_add(shear_slope, bulk_slope);
                let cross = (-2.0 / 3.0_f64).mul_add(shear_slope, bulk_slope);
                core::array::from_fn(|i| {
                    core::array::from_fn(|j| {
                        if i < 3 && j < 3 {
                            if i == j { diagonal } else { cross }
                        } else if i == j {
                            2.0 * shear_slope
                        } else {
                            0.0
                        }
                    })
                })
            }
            _ => {
                let moduli: Vec<&str> = match source.properties.elasticity {
                    TensorElasticity::Isotropic => {
                        if !fixed("poisson_ratio") {
                            return Err(invalid(
                                "temperature-dependent Poisson ratio requires a different non-affine thermodynamic law",
                            ));
                        }
                        vec!["young_modulus"]
                    }
                    TensorElasticity::IsotropicYoungShear => vec!["young_modulus", "shear_modulus"],
                    TensorElasticity::IsotropicBulkShear => {
                        unreachable!("bulk/shear handled directly")
                    }
                    TensorElasticity::Orthotropic => {
                        if ORTHOTROPIC_POISSON_RATIO_PROPERTIES
                            .iter()
                            .any(|key| !fixed(key))
                        {
                            return Err(invalid(
                                "orthotropic Poisson ratios must stay constant over the complete thermal interval",
                            ));
                        }
                        ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES
                            .into_iter()
                            .chain(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES)
                            .collect()
                    }
                };
                let ratio =
                    slopes[moduli[0]] / source.state.property(moduli[0]).unwrap().value_si();
                for name in moduli {
                    let this = slopes[name] / source.state.property(name).unwrap().value_si();
                    if !this.is_finite()
                        || !ratio.is_finite()
                        || (this - ratio).abs() > 128.0 * f64::EPSILON * this.abs().max(ratio.abs())
                    {
                        return Err(invalid(
                            "selected equilibrium moduli require one common affine scale so the stiffness tensor itself is affine",
                        ));
                    }
                }
                source
                    .law
                    .equilibrium_mandel()
                    .map(|row| row.map(|v| v * ratio))
            }
        };
        let mut branches_pa_k = Vec::with_capacity(source.properties.branches.len());
        for branch in &source.properties.branches {
            let mut matrix = [[0.0; 6]; 6];
            let mut entry = 0;
            for row in 0..6 {
                for column in row..6 {
                    let slope = slopes[&branch.stiffness_upper[entry]];
                    matrix[row][column] = slope;
                    matrix[column][row] = slope;
                    entry += 1;
                }
            }
            branches_pa_k.push(matrix);
        }
        if equilibrium_pa_k
            .iter()
            .flatten()
            .chain(branches_pa_k.iter().flatten().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(invalid("affine tensor derivatives overflow"));
        }
        // Confirm the derived tensor model against both exact original source
        // endpoint operators. Their SPD/PSD admission remains with TensorProny.
        for (endpoint, t) in endpoints.iter().zip(interval_k) {
            check_affine(
                source.law.equilibrium_mandel(),
                endpoint.law.equilibrium_mandel(),
                &equilibrium_pa_k,
                t - t0,
            )?;
            for ((old, new), slope) in source
                .law
                .branches()
                .iter()
                .zip(endpoint.law.branches())
                .zip(&branches_pa_k)
            {
                check_affine(
                    old.stiffness_mandel(),
                    new.stiffness_mandel(),
                    slope,
                    t - t0,
                )?;
            }
        }
        Ok(Self {
            source: source.clone(),
            temperature_axis: temperature_axis.into(),
            derivative: TensorTemperatureDerivative {
                equilibrium_pa_k,
                branches_pa_k,
                interval_k,
            },
        })
    }
    pub const fn source(&self) -> &ResolvedTensorMaterial {
        &self.source
    }
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    pub const fn temperature_interval_k(&self) -> [f64; 2] {
        self.derivative.interval_k
    }
    /// Rotate derivatives through the SAME admitted proper material frame.
    pub fn in_frame(&self, rotation: [[f64; 3]; 3]) -> Result<TensorTemperatureDerivative> {
        self.source.rotated_law(rotation)?;
        let derivative = TensorTemperatureDerivative {
            equilibrium_pa_k: crate::tensor::rotate_mandel(
                self.derivative.equilibrium_pa_k,
                rotation,
            ),
            branches_pa_k: self
                .derivative
                .branches_pa_k
                .iter()
                .map(|m| crate::tensor::rotate_mandel(*m, rotation))
                .collect(),
            interval_k: self.derivative.interval_k,
        };
        if derivative
            .equilibrium_pa_k
            .iter()
            .flatten()
            .chain(derivative.branches_pa_k.iter().flatten().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(TensorSourceError::Invalid(
                "rotated thermal derivative overflow".into(),
            ));
        }
        Ok(derivative)
    }
}

fn check_affine(old: &Tangent6, new: &Tangent6, slope: &Tangent6, dt: f64) -> Result<()> {
    let scale = old
        .iter()
        .flatten()
        .chain(new.iter().flatten())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    for i in 0..6 {
        for j in 0..6 {
            let expected = slope[i][j].mul_add(dt, old[i][j]);
            if !expected.is_finite() || (expected - new[i][j]).abs() > 2048.0 * f64::EPSILON * scale
            {
                return Err(TensorSourceError::Invalid(
                    "source tensor is not the admitted affine thermodynamic operator".into(),
                ));
            }
        }
    }
    Ok(())
}

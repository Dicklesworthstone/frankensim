//! Source-derived tensor free energy, material entropy and thermal curvature.
//!
//! At fixed total and viscous strain, psi = (e:Ceq(T):e + sum q:Cj(T):q)/2.
//! Affine tensors give S = -psi_T and U = psi + T*S, with zero fixed-strain
//! mechanical heat-capacity contribution. General affine engineering constants
//! retain exact rational bulk or inverse-compliance response and thermal curvature.
//! Reference density and geometry remain fixed. Source interpolation is admitted over every
//! original knot; merely querying two endpoints does not establish this model.

use super::*;
use crate::Tangent6;
use crate::state_point::{affine_temperature_slope, pinned_temperature_line};
use fs_qty::{
    Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::collections::BTreeMap;
use std::sync::Arc;

mod compliance;
use compliance::ComplianceThermal;

/// Symmetric operator-norm bounds over the entire original temperature interval.
/// Symmetric operator norms remain invariant under the admitted material rotation.
#[derive(Clone, Copy, Debug)]
pub(crate) struct TensorTemperatureBounds {
    pub(crate) stiffness_norm_pa: f64,
    pub(crate) slope_norm_pa_k: f64,
    pub(crate) curvature_norm_pa_k2: f64,
}

/// Derivatives of the original stored tensors in the consuming strain frame.
/// A nonlinear snapshot is tied to one exact temperature; refresh it before
/// evaluating another endpoint. Immutable branch slopes are shared by trials.
#[derive(Clone, Debug)]
pub struct TensorTemperatureDerivative {
    pub(crate) equilibrium_pa_k: Tangent6,
    pub(crate) equilibrium_pa_k2: Tangent6,
    pub(crate) branches_pa_k: Arc<[Tangent6]>,
    pub(crate) interval_k: [f64; 2],
    pub(crate) resolved_temperature_k: Option<f64>,
    pub(crate) equilibrium_bounds: Option<TensorTemperatureBounds>,
    young_shear: Option<AffineYoungShear>,
    compliance: Option<ComplianceThermal>,
}

#[derive(Clone, Copy, Debug)]
struct AffineYoungShear {
    reference_temperature_k: f64,
    young_pa: f64,
    shear_pa: f64,
    young_slope_pa_k: f64,
    shear_slope_pa_k: f64,
    // (E'G-EG')/G_reference is constant along both admitted affine paths.
    cross_pa_k: f64,
}

#[derive(Clone, Copy)]
struct IsotropicThermalResponse {
    bulk_pa: f64,
    shear_pa: f64,
    bulk_slope_pa_k: f64,
    shear_slope_pa_k: f64,
    bulk_curvature_pa_k2: f64,
}

impl AffineYoungShear {
    fn at_temperature(self, temperature_k: f64) -> Result<IsotropicThermalResponse> {
        let dt = temperature_k - self.reference_temperature_k;
        let young = self.young_slope_pa_k.mul_add(dt, self.young_pa);
        let shear = self.shear_slope_pa_k.mul_add(dt, self.shear_pa);
        let ratio = young / shear;
        let gap = 3.0 - ratio;
        if !young.is_finite()
            || !shear.is_finite()
            || young <= 0.0
            || shear <= 0.0
            || !gap.is_finite()
            || gap <= 0.0
        {
            return Err(TensorSourceError::Invalid(
                "thermal Young/shear sources leave positive isotropic stability".into(),
            ));
        }
        // K=E G/[3(3G-E)]=E/[3(3-E/G)]. Normalize before multiplying
        // pressures, so a representable response never requires forming E*G.
        let bulk = (young / 3.0) / gap;
        let young_rate = self.young_slope_pa_k / young;
        let ratio_rate = (self.shear_pa / shear) * (self.cross_pa_k / shear) / gap;
        let bulk_slope = bulk * (young_rate + ratio_rate);
        // W=E'G-EG' is constant and K''=2W^2/(3G-E)^3 >= 0.
        // This factored form avoids both the squared pressures and W itself.
        let bulk_curvature = (shear * ratio_rate) * (2.0 * ratio_rate / gap);
        if !bulk.is_finite()
            || bulk <= 0.0
            || !bulk_slope.is_finite()
            || !bulk_curvature.is_finite()
            || bulk_curvature < 0.0
        {
            return Err(TensorSourceError::Invalid(
                "rational isotropic thermal response is unrepresentable".into(),
            ));
        }
        Ok(IsotropicThermalResponse {
            bulk_pa: bulk,
            shear_pa: shear,
            bulk_slope_pa_k: bulk_slope,
            shear_slope_pa_k: self.shear_slope_pa_k,
            bulk_curvature_pa_k2: bulk_curvature,
        })
    }
}

fn isotropic_mandel(bulk: f64, shear: f64) -> Tangent6 {
    let diagonal = (4.0 / 3.0_f64).mul_add(shear, bulk);
    let cross = (-2.0 / 3.0_f64).mul_add(shear, bulk);
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i < 3 && j < 3 {
                if i == j { diagonal } else { cross }
            } else if i == j {
                2.0 * shear
            } else {
                0.0
            }
        })
    })
}

impl TensorTemperatureDerivative {
    /// Whether the complete stored tensor is affine over its admitted interval.
    /// A projection can remove nonlinear directions, but cannot widen this
    /// source-level guarantee for consumers caching one constant derivative.
    pub const fn is_affine(&self) -> bool {
        self.resolved_temperature_k.is_none()
    }

    pub(crate) fn from_affine_parts(
        equilibrium: Tangent6,
        branches: Vec<Tangent6>,
        interval: [f64; 2],
    ) -> Self {
        Self {
            equilibrium_pa_k: equilibrium,
            equilibrium_pa_k2: [[0.0; 6]; 6],
            branches_pa_k: branches.into(),
            interval_k: interval,
            resolved_temperature_k: None,
            equilibrium_bounds: None,
            young_shear: None,
            compliance: None,
        }
    }

    /// Refresh the same admitted source law in its existing frame. The original
    /// card paths, branch slopes and interval bounds are immutable and shared;
    /// this operation performs no material queries, allocation or branch rotation.
    pub fn at_temperature(&self, temperature_k: f64) -> Result<Self> {
        if !temperature_k.is_finite()
            || !(self.interval_k[0]..=self.interval_k[1]).contains(&temperature_k)
        {
            return Err(TensorSourceError::Invalid(
                "thermal derivative temperature leaves its original source interval".into(),
            ));
        }
        let mut result = self.clone();
        if let Some(pair) = self.young_shear {
            let response = pair.at_temperature(temperature_k)?;
            result.equilibrium_pa_k =
                isotropic_mandel(response.bulk_slope_pa_k, response.shear_slope_pa_k);
            result.equilibrium_pa_k2 = isotropic_mandel(response.bulk_curvature_pa_k2, 0.0);
            result.resolved_temperature_k = Some(temperature_k);
        } else if let Some(model) = self.compliance {
            let jet = model.principal_jet(temperature_k)?;
            result.equilibrium_pa_k = model.rotate(jet[1]);
            result.equilibrium_pa_k2 = model.rotate(jet[2]);
            if result
                .equilibrium_pa_k
                .iter()
                .flatten()
                .chain(result.equilibrium_pa_k2.iter().flatten())
                .any(|v| !v.is_finite())
            {
                return Err(TensorSourceError::Invalid(
                    "rotated thermal compliance derivative overflow".into(),
                ));
            }
            result.resolved_temperature_k = Some(temperature_k);
        }
        Ok(result)
    }
}

/// Original source selection and complete bounded tensor thermodynamic law.
/// Stored Maxwell branches remain affine; nonlinear equilibrium sources retain
/// their exact rational bulk or inverse-compliance response.
#[derive(Clone, Debug)]
pub struct TensorThermalLaw {
    source: ResolvedTensorMaterial,
    temperature_axis: String,
    derivative: TensorTemperatureDerivative,
}
impl TensorThermalLaw {
    /// Admit independent affine K(T)/G(T), E(T) with fixed isotropic nu,
    /// independent affine E(T)/G(T), E(T)/nu(T), or independently affine
    /// orthotropic E/G/major Poisson ratios. Isotropic compliance uses exact
    /// eigenvalue interval bounds, including auxetic sources. Orthotropic paths
    /// require a positive scaled diagonal-dominance bound over the whole interval;
    /// this sufficient proof can refuse valid non-dominant sources. Existing
    /// affine tensor limits retain their original admission. Each of the 21 coefficients
    /// of every Maxwell branch may have its OWN affine slope, including signed
    /// cross terms. Relaxation times retain positive piecewise-linear support.
    ///
    /// Every source axis/validity restriction remains owned. Changing branch
    /// factors preserve the same observable strain range; changed rank/nullspace
    /// or unrecoverable factor rows retain the native history-transfer refusal.
    /// This law supplies neither thermal expansion nor a finite-strain model.
    pub fn try_new(
        source: &ResolvedTensorMaterial,
        card: &MaterialCard,
        temperature_axis: &str,
        interval_k: [f64; 2],
    ) -> Result<Self> {
        Self::admit(source, card, temperature_axis, interval_k, true)
    }

    fn admit(
        source: &ResolvedTensorMaterial,
        card: &MaterialCard,
        temperature_axis: &str,
        interval_k: [f64; 2],
        nonlinear_equilibrium: bool,
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
        let mut young_shear = None;
        let mut compliance = None;
        let mut equilibrium_pa_k2 = [[0.0; 6]; 6];
        let mut equilibrium_bounds = None;
        let equilibrium_pa_k = match source.properties.elasticity {
            TensorElasticity::Isotropic | TensorElasticity::Orthotropic
                if nonlinear_equilibrium && affine_equilibrium_scale(source, &slopes).is_err() =>
            {
                let (model, bounds) = ComplianceThermal::admit(source, &slopes, t0, interval_k)?;
                let reference = model.principal_jet(t0)?;
                check_affine(
                    &reference[0],
                    source.law.equilibrium_mandel(),
                    &[[0.0; 6]; 6],
                    0.0,
                )?;
                equilibrium_pa_k2 = reference[2];
                equilibrium_bounds = Some(bounds);
                compliance = Some(model);
                reference[1]
            }
            TensorElasticity::IsotropicYoungShear if nonlinear_equilibrium => {
                let young = source.state.property("young_modulus").unwrap().value_si();
                let shear = source.state.property("shear_modulus").unwrap().value_si();
                let young_slope = slopes["young_modulus"];
                let shear_slope = slopes["shear_modulus"];
                let young_rate = young_slope / young;
                let shear_rate = shear_slope / shear;
                if young_rate.is_finite() && young_rate == shear_rate {
                    // Preserve the original exactly proportional affine limit.
                    source
                        .law
                        .equilibrium_mandel()
                        .map(|row| row.map(|v| v * young_rate))
                } else {
                    let pair = AffineYoungShear {
                        reference_temperature_k: t0,
                        young_pa: young,
                        shear_pa: shear,
                        young_slope_pa_k: young_slope,
                        shear_slope_pa_k: shear_slope,
                        cross_pa_k: (-young / shear).mul_add(shear_slope, young_slope),
                    };
                    let first = pair.at_temperature(low)?;
                    let last = pair.at_temperature(high)?;
                    let reference = pair.at_temperature(t0)?;
                    // q=3G-E is positive affine, W=E'G-EG' is constant,
                    // and K''=2W^2/q^3. Curvature attains its maximum where
                    // q is smallest; convex K and monotone K' attain their
                    // absolute upper bounds at the interval endpoints.
                    let guard = 1.0 + 4096.0 * f64::EPSILON;
                    let bounds = TensorTemperatureBounds {
                        stiffness_norm_pa: (3.0 * first.bulk_pa.max(last.bulk_pa))
                            .max(2.0 * first.shear_pa.max(last.shear_pa))
                            * guard,
                        slope_norm_pa_k: (3.0
                            * first.bulk_slope_pa_k.abs().max(last.bulk_slope_pa_k.abs()))
                        .max(2.0 * shear_slope.abs())
                            * guard,
                        curvature_norm_pa_k2: 3.0
                            * first.bulk_curvature_pa_k2.max(last.bulk_curvature_pa_k2)
                            * guard,
                    };
                    if [
                        bounds.stiffness_norm_pa,
                        bounds.slope_norm_pa_k,
                        bounds.curvature_norm_pa_k2,
                    ]
                    .iter()
                    .any(|v| !v.is_finite() || *v < 0.0)
                    {
                        return Err(invalid("rational tensor thermal interval bounds overflow"));
                    }
                    equilibrium_pa_k2 = isotropic_mandel(reference.bulk_curvature_pa_k2, 0.0);
                    equilibrium_bounds = Some(bounds);
                    young_shear = Some(pair);
                    isotropic_mandel(reference.bulk_slope_pa_k, reference.shear_slope_pa_k)
                }
            }
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
                let ratio = affine_equilibrium_scale(source, &slopes)?;
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
            if let Some(pair) = young_shear {
                let response = pair.at_temperature(t)?;
                check_affine(
                    &isotropic_mandel(response.bulk_pa, response.shear_pa),
                    endpoint.law.equilibrium_mandel(),
                    &[[0.0; 6]; 6],
                    0.0,
                )?;
            } else if let Some(model) = compliance {
                check_affine(
                    &model.principal_jet(t)?[0],
                    endpoint.law.equilibrium_mandel(),
                    &[[0.0; 6]; 6],
                    0.0,
                )?;
            } else {
                check_affine(
                    source.law.equilibrium_mandel(),
                    endpoint.law.equilibrium_mandel(),
                    &equilibrium_pa_k,
                    t - t0,
                )?;
            }
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
        let mut derivative = TensorTemperatureDerivative::from_affine_parts(
            equilibrium_pa_k,
            branches_pa_k,
            interval_k,
        );
        derivative.equilibrium_pa_k2 = equilibrium_pa_k2;
        derivative.resolved_temperature_k =
            (young_shear.is_some() || compliance.is_some()).then_some(t0);
        derivative.equilibrium_bounds = equilibrium_bounds;
        derivative.young_shear = young_shear;
        derivative.compliance = compliance;
        Ok(Self {
            source: source.clone(),
            temperature_axis: temperature_axis.into(),
            derivative,
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
        let compliance = self
            .derivative
            .compliance
            .map(|model| model.in_frame(rotation));
        let derivative = TensorTemperatureDerivative {
            // Isotropic derivatives are frame invariant. Anisotropic compliance
            // snapshots retain the cached Bond map used by every later endpoint.
            equilibrium_pa_k: if self.derivative.young_shear.is_some() {
                self.derivative.equilibrium_pa_k
            } else if let Some(model) = compliance {
                model.rotate(self.derivative.equilibrium_pa_k)
            } else {
                crate::tensor::rotate_mandel(self.derivative.equilibrium_pa_k, rotation)
            },
            equilibrium_pa_k2: compliance.map_or(self.derivative.equilibrium_pa_k2, |model| {
                model.rotate(self.derivative.equilibrium_pa_k2)
            }),
            branches_pa_k: self
                .derivative
                .branches_pa_k
                .iter()
                .map(|m| crate::tensor::rotate_mandel(*m, rotation))
                .collect(),
            interval_k: self.derivative.interval_k,
            resolved_temperature_k: self.derivative.resolved_temperature_k,
            equilibrium_bounds: self.derivative.equilibrium_bounds,
            young_shear: self.derivative.young_shear,
            compliance,
        };
        if derivative
            .equilibrium_pa_k
            .iter()
            .flatten()
            .chain(derivative.branches_pa_k.iter().flatten().flatten())
            .chain(derivative.equilibrium_pa_k2.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(TensorSourceError::Invalid(
                "rotated thermal derivative overflow".into(),
            ));
        }
        Ok(derivative)
    }
}

/// Explicitly affine subset of the shared tensor thermal source owner.
/// Consumers retaining one constant stiffness derivative must use this
/// admission; non-affine engineering-constant sources still refuse here.
#[derive(Clone, Debug)]
pub struct AffineTensorThermal {
    law: TensorThermalLaw,
}

impl AffineTensorThermal {
    pub fn try_new(
        source: &ResolvedTensorMaterial,
        card: &MaterialCard,
        temperature_axis: &str,
        interval_k: [f64; 2],
    ) -> Result<Self> {
        TensorThermalLaw::admit(source, card, temperature_axis, interval_k, false)
            .map(|law| Self { law })
    }

    pub const fn source(&self) -> &ResolvedTensorMaterial {
        self.law.source()
    }

    pub fn temperature_axis(&self) -> &str {
        self.law.temperature_axis()
    }

    pub const fn temperature_interval_k(&self) -> [f64; 2] {
        self.law.temperature_interval_k()
    }

    pub fn in_frame(&self, rotation: [[f64; 3]; 3]) -> Result<TensorTemperatureDerivative> {
        self.law.in_frame(rotation)
    }
}

fn affine_equilibrium_scale(
    source: &ResolvedTensorMaterial,
    slopes: &BTreeMap<String, f64>,
) -> Result<f64> {
    let invalid = |s: &str| TensorSourceError::Invalid(s.into());
    let fixed = |name: &str| slopes.get(name).is_some_and(|s| *s == 0.0);
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
    let ratio = slopes[moduli[0]] / source.state.property(moduli[0]).unwrap().value_si();
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
    Ok(ratio)
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

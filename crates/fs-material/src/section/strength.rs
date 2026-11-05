//! Axial tensile capacity of a bonded section from per-layer sourced limits.
//!
//! A perfectly bonded section under axial force N alone has one common strain
//! N/EA, so layer i carries the stress `E_i N / EA`. Layer i reaches its own
//! limit at `N_i = limit_i EA / E_i`; the section capacity is the smallest
//! `N_i`, and that layer governs. The caller names the limit property (a yield
//! stress, proof stress or ultimate tensile strength): this module never picks
//! one, applies a safety factor or turns a typical strength into an allowable.
//!
//! Every limit enters with its fs-matdb answer, and the capacity keeps it: a
//! limit read from a layer's state bundle, or one resolved at another state
//! point (a strength-by-diameter claim at the drawn wire's own diameter),
//! carries its pack, claim and evaluation decision into the capacity's
//! receipt. No constructor takes a bare number.
//!
//! The original force-only capacity is complemented by signed axial/bending
//! assessments and an assessment of the existing thermal section's actual
//! layer stresses. Perfect bonding, concentric geometry and small-strain linear
//! elasticity are assumed. Residual drawing/winding stress, shear failure,
//! yielding redistribution, fracture evolution, buckling and fatigue are absent.

use super::{ConcentricElasticSection, MaterialStatePointError, invalid, positive};
use super::thermal::ConcentricThermoelasticSection;
use crate::state_point::ResolvedScalarProperty;
use fs_qty::{Pressure, QuantitySpec};

/// Sourced limits of every layer of one bonded section under axial force.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialTensileCapacity {
    limit_property: String,
    section: ConcentricElasticSection,
    layer_limits_pa: Vec<f64>,
    layer_limit_sources: Vec<ResolvedScalarProperty>,
    capacity_n: f64,
    governing_layer: usize,
}

/// Axial normal-stress extrema on one complete annulus, tension positive.
/// Compression is reported, but is not compared with a tensile limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerTensileAssessment {
    /// Minimum normal stress over this layer [Pa].
    pub minimum_axial_stress_pa: f64,
    /// Maximum normal stress over this layer [Pa].
    pub maximum_axial_stress_pa: f64,
    /// `max(0, maximum_axial_stress_pa) / layer_limit_pa`.
    pub tensile_utilization: f64,
}

/// Static tensile-limit assessment in the capacity's original layer order.
/// This result neither changes material state nor asserts a failure mechanism.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricTensileAssessment {
    layers: Vec<LayerTensileAssessment>,
    governing_layer: usize,
}

impl ConcentricTensileAssessment {
    /// Extrema and utilization of each original layer.
    #[must_use]
    pub fn layers(&self) -> &[LayerTensileAssessment] { &self.layers }
    /// Layer with the greatest tensile utilization; the first wins an exact tie.
    /// If all layers are in compression, this is zero and utilization is zero.
    #[must_use]
    pub const fn governing_layer(&self) -> usize { self.governing_layer }
    /// Largest positive tensile stress divided by its own layer's limit.
    #[must_use]
    pub fn utilization(&self) -> f64 { self.layers[self.governing_layer].tensile_utilization }
    /// Whether at least one layer exceeds its declared limit, without tolerance.
    #[must_use]
    pub fn exceeds_limit(&self) -> bool { self.utilization() > 1.0 }
}

impl ConcentricElasticSection {
    /// Bind the named tensile-limit claim [Pa] from EVERY layer's state bundle.
    ///
    /// The limit must be a positive constant at the bundle's fixed state, like
    /// the density and modulus it is compared with. A layer without the claim is
    /// missing data: it refuses rather than being treated as unbreakable.
    ///
    /// # Errors
    /// Missing, ambiguous, wrong-quantity, nonpositive or state-varying limits
    /// and unrepresentable capacities refuse; no partial capacity is returned.
    pub fn axial_tensile_capacity(
        &self,
        limit_property: &str,
    ) -> Result<AxialTensileCapacity, MaterialStatePointError> {
        let mut sources = Vec::with_capacity(self.layers().len());
        for layer in self.layers() {
            let limit = layer.material.property(limit_property)
                .ok_or_else(|| invalid("every tension-carrying layer needs its named tensile-limit claim; missing strength data refuses"))?;
            if limit.requirement().quantity() != QuantitySpec::dimensional(Pressure::DIMS) {
                return Err(invalid("tensile limit must be dimension-only Pa"));
            }
            if !limit.is_static_at_fixed_external_conditions() {
                return Err(invalid("tensile limit requires fixed external conditions; unowned rate, strain or frequency support refuses"));
            }
            sources.push(limit.clone());
        }
        self.capacity_from_sources(limit_property, sources)
    }

    /// Bind tensile limits resolved at state points other than the layers'
    /// fixed state bundles, one per layer in spatial order: for example a
    /// strength-by-diameter claim evaluated at the drawn wire's own diameter.
    /// Each limit is the complete fs-matdb answer for `limit_property`, and
    /// the capacity keeps it as that layer's limit receipt (pack, claim and
    /// interpolation decision). Which claim belongs to which layer is the
    /// caller's declaration; the capacity arithmetic is the same as
    /// [`Self::axial_tensile_capacity`].
    ///
    /// # Errors
    /// A limit count other than the layer count, an answer for another
    /// property or quantity, nonpositive or nonfinite limits and
    /// unrepresentable capacities refuse.
    pub fn axial_tensile_capacity_from_resolved(
        &self,
        limit_property: &str,
        layer_limits: &[ResolvedScalarProperty],
    ) -> Result<AxialTensileCapacity, MaterialStatePointError> {
        if layer_limits.len() != self.layers().len() {
            return Err(invalid("one resolved tensile limit per section layer"));
        }
        for limit in layer_limits {
            if limit.requirement().name() != limit_property {
                return Err(invalid("a resolved tensile limit answers another property than the declared limit"));
            }
            if limit.requirement().quantity() != QuantitySpec::dimensional(Pressure::DIMS) {
                return Err(invalid("tensile limit must be dimension-only Pa"));
            }
        }
        self.capacity_from_sources(limit_property, layer_limits.to_vec())
    }

    fn capacity_from_sources(
        &self,
        limit_property: &str,
        sources: Vec<ResolvedScalarProperty>,
    ) -> Result<AxialTensileCapacity, MaterialStatePointError> {
        let axial = self.axial_stiffness_n();
        let mut layer_limits_pa = Vec::with_capacity(self.layers().len());
        let mut capacity_n = f64::INFINITY;
        let mut governing_layer = 0;
        for (index, (geometry, source)) in self.layer_integrals().iter().zip(&sources).enumerate() {
            let limit = positive(source.value_si(), "positive finite tensile limit")?;
            let modulus = geometry.young_modulus_pa();
            let layer_capacity = positive(limit * axial / modulus, "layer axial-force capacity")?;
            if layer_capacity < capacity_n {
                capacity_n = layer_capacity;
                governing_layer = index;
            }
            layer_limits_pa.push(limit);
        }
        Ok(AxialTensileCapacity {
            limit_property: limit_property.into(), section: self.clone(),
            layer_limits_pa, layer_limit_sources: sources, capacity_n, governing_layer,
        })
    }
}

impl AxialTensileCapacity {
    /// Original geometry, constitutive parameters and material source bundles.
    #[must_use]
    pub const fn section(&self) -> &ConcentricElasticSection { &self.section }
    /// Exact property name whose claims supplied the limits.
    #[must_use]
    pub fn limit_property(&self) -> &str { &self.limit_property }
    /// Per-layer limits [Pa], in the section's spatial layer order.
    #[must_use]
    pub fn layer_limits_pa(&self) -> &[f64] { &self.layer_limits_pa }
    /// Per-layer fs-matdb answers behind [`Self::layer_limits_pa`], in the same
    /// order: each limit's pack, claim and evaluation decision.
    #[must_use]
    pub fn layer_limit_sources(&self) -> &[ResolvedScalarProperty] { &self.layer_limit_sources }
    /// Smallest axial force [N] at which some layer reaches its limit.
    #[must_use]
    pub const fn capacity_n(&self) -> f64 { self.capacity_n }
    /// Index of the first layer to reach its limit as the force grows.
    #[must_use]
    pub const fn governing_layer(&self) -> usize { self.governing_layer }

    /// Stress [Pa] in one layer under axial force `axial_force_n` [N].
    ///
    /// # Errors
    /// A layer index outside the section or a negative/nonfinite force refuses.
    pub fn layer_stress_pa(&self, layer: usize, axial_force_n: f64)
        -> Result<f64, MaterialStatePointError> {
        let modulus = self.section.layer_integrals().get(layer)
            .ok_or_else(|| invalid("tensile check layer index is outside the section"))?
            .young_modulus_pa();
        if !axial_force_n.is_finite() || axial_force_n < 0.0 {
            return Err(invalid("tensile check needs a finite nonnegative axial force"));
        }
        let stress = modulus * (axial_force_n / self.section.axial_stiffness_n());
        if stress.is_finite() { Ok(stress) } else { Err(invalid("layer axial stress overflow")) }
    }

    /// Force over capacity: above one, the governing layer exceeds its limit.
    ///
    /// # Errors
    /// A negative/nonfinite force or unrepresentable utilization refuses.
    pub fn utilization(&self, axial_force_n: f64) -> Result<f64, MaterialStatePointError> {
        if !axial_force_n.is_finite() || axial_force_n < 0.0 {
            return Err(invalid("tensile check needs a finite nonnegative axial force"));
        }
        checked_utilization(axial_force_n, self.capacity_n)
    }

    /// Assess signed axial force [N] and bending moments [N m] about two
    /// orthogonal centroidal transverse axes. For a concentric section their
    /// orientation does not change the annular extrema:
    /// `sigma_i = E_i N/EA +/- E_i r_outer_i |M|/EI`.
    ///
    /// Compression is admitted and can coexist with tensile bending stress.
    /// The reported normal stresses use elastic load sharing even above the
    /// named limit; plastic redistribution, shear/torsion and buckling are not
    /// assessed. No safety factor or fracture law is inferred.
    ///
    /// # Errors
    /// Nonfinite loads and unrepresentable stress or utilization refuse.
    pub fn assess_axial_bending(
        &self,
        axial_force_n: f64,
        bending_moment_n_m: [f64; 2],
    ) -> Result<ConcentricTensileAssessment, MaterialStatePointError> {
        if !axial_force_n.is_finite() {
            return Err(invalid("section tensile assessment requires finite signed axial force"));
        }
        let strain = axial_force_n / self.section.axial_stiffness_n();
        if !strain.is_finite() || axial_force_n != 0.0 && strain == 0.0 {
            return Err(invalid("section tensile assessment axial strain is unrepresentable"));
        }
        self.assess_stresses(bending_moment_n_m, |i, _| {
            let stress = self.section.layer_integrals()[i].young_modulus_pa() * strain;
            if strain != 0.0 && stress == 0.0 {
                return Err(invalid("section tensile assessment axial stress is unrepresentable"));
            }
            Ok(stress)
        })
    }

    /// Assess the existing thermoelastic section's actual per-layer axial
    /// stresses, plus the same concentric bending stress as above. The thermal
    /// section must exactly match the retained geometry and material sources.
    /// Its native integrated eigenstrains and present moduli are reused; a
    /// small or zero resultant cannot conceal self-equilibrated layer tension.
    ///
    /// `total_axial_strain` is the common centroidal engineering strain.
    /// Both total and elastic strain extrema, including bending, must lie
    /// inside the explicitly declared small-strain limit. Thermal strain is
    /// bounded by the original thermal response as well.
    ///
    /// # Errors
    /// Foreign/stale sections, invalid thermal response, nonfinite moments,
    /// excessive strain extrema or unrepresentable arithmetic refuse atomically.
    pub fn assess_thermoelastic_axial_bending(
        &self,
        thermal: &ConcentricThermoelasticSection,
        total_axial_strain: f64,
        bending_moment_n_m: [f64; 2],
        linear_strain_limit: f64,
    ) -> Result<ConcentricTensileAssessment, MaterialStatePointError> {
        if thermal.section() != &self.section {
            return Err(invalid("thermal tensile assessment must use the exact original section and material states"));
        }
        let response = thermal.response_at_axial_strain(total_axial_strain, linear_strain_limit)?;
        self.assess_stresses(bending_moment_n_m, |i, bending_strain| {
            let layer = &response.layers()[i];
            if total_axial_strain.abs() + bending_strain > linear_strain_limit
                || layer.elastic_strain().abs() + bending_strain > linear_strain_limit
            {
                return Err(invalid("thermal tensile assessment strain extrema exceed the declared linear regime"));
            }
            Ok(layer.stress_pa())
        })
    }

    fn assess_stresses(
        &self,
        bending_moment_n_m: [f64; 2],
        mut mean_stress: impl FnMut(usize, f64) -> Result<f64, MaterialStatePointError>,
    ) -> Result<ConcentricTensileAssessment, MaterialStatePointError> {
        if bending_moment_n_m.iter().any(|x| !x.is_finite()) {
            return Err(invalid("section tensile assessment requires finite bending moments"));
        }
        // Scale before squaring, retaining finite large/small moments.
        let scale = bending_moment_n_m[0].abs().max(bending_moment_n_m[1].abs());
        let curvature = if scale == 0.0 { 0.0 } else {
            let x = bending_moment_n_m[0] / scale;
            let y = bending_moment_n_m[1] / scale;
            (scale / self.section.bending_stiffness_n_m2()) * fs_math::det::sqrt(x * x + y * y)
        };
        if !curvature.is_finite() || scale > 0.0 && curvature == 0.0 {
            return Err(invalid("section tensile assessment curvature is unrepresentable"));
        }
        let mut layers = Vec::with_capacity(self.section.layers().len());
        let mut governing_layer = 0;
        let mut maximum_utilization = 0.0;
        for (i, layer) in self.section.layers().iter().enumerate() {
            let bending_strain = curvature * layer.outer_radius_m;
            if !bending_strain.is_finite() || curvature > 0.0 && bending_strain == 0.0 {
                return Err(invalid("section tensile assessment bending strain is unrepresentable"));
            }
            let mean = mean_stress(i, bending_strain)?;
            let amplitude = self.section.layer_integrals()[i].young_modulus_pa() * bending_strain;
            if bending_strain > 0.0 && amplitude == 0.0 {
                return Err(invalid("section tensile assessment bending stress is unrepresentable"));
            }
            let minimum = mean - amplitude;
            let maximum = mean + amplitude;
            if !minimum.is_finite() || !maximum.is_finite() {
                return Err(invalid("section tensile assessment normal-stress extrema are unrepresentable"));
            }
            let utilization = checked_utilization(maximum.max(0.0), self.layer_limits_pa[i])?;
            if utilization > maximum_utilization {
                maximum_utilization = utilization;
                governing_layer = i;
            }
            layers.push(LayerTensileAssessment { minimum_axial_stress_pa: minimum,
                maximum_axial_stress_pa: maximum, tensile_utilization: utilization });
        }
        Ok(ConcentricTensileAssessment { layers, governing_layer })
    }
}

fn checked_utilization(demand: f64, capacity: f64) -> Result<f64, MaterialStatePointError> {
    let ratio = demand / capacity;
    if !ratio.is_finite() || demand > 0.0 && ratio == 0.0 {
        Err(invalid("tensile utilization is unrepresentable"))
    } else {
        Ok(ratio)
    }
}

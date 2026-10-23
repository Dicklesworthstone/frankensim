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
//! Only the axial-force stress is checked. Residual drawing or winding stress,
//! thermal eigenstrain, the bending stress of vibration and fatigue are absent.

use super::{ConcentricElasticSection, MaterialStatePointError, invalid, positive};
use crate::state_point::{ResolvedScalarProperty, YOUNG_MODULUS_PROPERTY};
use fs_qty::{Pressure, QuantitySpec};

/// Sourced limits of every layer of one bonded section under axial force.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialTensileCapacity {
    limit_property: String,
    axial_stiffness_n: f64,
    layer_moduli_pa: Vec<f64>,
    layer_limits_pa: Vec<f64>,
    layer_limit_sources: Vec<ResolvedScalarProperty>,
    capacity_n: f64,
    governing_layer: usize,
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
            if !limit.is_constant_at_fixed_temperature() {
                return Err(invalid("tensile limit must be a constant or T-only curve at fixed positive absolute T"));
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
        let mut layer_moduli_pa = Vec::with_capacity(self.layers().len());
        let mut layer_limits_pa = Vec::with_capacity(self.layers().len());
        let mut capacity_n = f64::INFINITY;
        let mut governing_layer = 0;
        for (index, (layer, source)) in self.layers().iter().zip(&sources).enumerate() {
            let limit = positive(source.value_si(), "positive finite tensile limit")?;
            let modulus = layer.material.property(YOUNG_MODULUS_PROPERTY)
                .ok_or_else(|| invalid("tensile check lost a layer modulus"))?.value_si();
            let modulus = positive(modulus, "positive finite layer modulus")?;
            let layer_capacity = positive(limit * axial / modulus, "layer axial-force capacity")?;
            if layer_capacity < capacity_n {
                capacity_n = layer_capacity;
                governing_layer = index;
            }
            layer_moduli_pa.push(modulus);
            layer_limits_pa.push(limit);
        }
        Ok(AxialTensileCapacity {
            limit_property: limit_property.into(), axial_stiffness_n: axial,
            layer_moduli_pa, layer_limits_pa, layer_limit_sources: sources, capacity_n, governing_layer,
        })
    }
}

impl AxialTensileCapacity {
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
        let modulus = *self.layer_moduli_pa.get(layer)
            .ok_or_else(|| invalid("tensile check layer index is outside the section"))?;
        if !axial_force_n.is_finite() || axial_force_n < 0.0 {
            return Err(invalid("tensile check needs a finite nonnegative axial force"));
        }
        let stress = modulus * (axial_force_n / self.axial_stiffness_n);
        if stress.is_finite() { Ok(stress) } else { Err(invalid("layer axial stress overflow")) }
    }

    /// Force over capacity: above one, the governing layer exceeds its limit.
    ///
    /// # Errors
    /// A negative or nonfinite force refuses.
    pub fn utilization(&self, axial_force_n: f64) -> Result<f64, MaterialStatePointError> {
        if !axial_force_n.is_finite() || axial_force_n < 0.0 {
            return Err(invalid("tensile check needs a finite nonnegative axial force"));
        }
        Ok(axial_force_n / self.capacity_n)
    }
}

//! Source-bound sections for slender, perfectly bonded elastic members.
//!
//! The concentric geometry has a common centroid and neutral axis, so the
//! one-dimensional parallel-strain model gives `mass/length = sum(rho A)`,
//! `EA = sum(E A)` and `EI = sum(E I)`. Every layer keeps its own material
//! state; there is no fictitious homogeneous material or averaged density.
//!
//! These are Euler-Bernoulli section integrals, not a three-dimensional
//! composite constitutive law: shear deformation, Poisson constraint, slip,
//! winding geometry, damage and dissipation are not modeled. The explicit
//! [`thermal`] binding adds per-layer thermal eigenstrain and axial stress.
//! Source validity applies at each retained query point, not at future states.

pub mod thermal;

use fs_qty::{Density, Pressure, QuantitySpec};
use crate::state_point::{
    DENSITY_PROPERTY, MaterialStatePointError, ResolvedMaterialStatePoint,
    YOUNG_MODULUS_PROPERTY,
};

/// Hard cold-path bound; no work proportional to an unbounded layer count.
pub const MAX_CONCENTRIC_LAYERS: usize = 128;

/// One annular region. Its inner radius is the previous region's outer radius,
/// or the section's explicitly declared bore radius for the first region.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricMaterialLayer {
    /// Outer radius in metres, strictly increasing in declaration order.
    pub outer_radius_m: f64,
    /// Exact source bundle containing dimensioned density and Young's modulus.
    pub material: ResolvedMaterialStatePoint,
}

/// The same per-layer integrals used in every aggregate section coefficient.
/// Their source is the corresponding entry in `ConcentricElasticSection::layers`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticLayerIntegrals {
    area_m2: f64,
    second_moment_m4: f64,
    line_density_kg_m: f64,
    axial_stiffness_n: f64,
    bending_stiffness_n_m2: f64,
}
impl ElasticLayerIntegrals {
    /// Layer cross-sectional area [m2].
    #[must_use]
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Layer second moment about the shared centroid [m4].
    #[must_use]
    pub const fn second_moment_m4(&self) -> f64 { self.second_moment_m4 }
    /// Layer mass per length [kg/m].
    #[must_use]
    pub const fn line_density_kg_m(&self) -> f64 { self.line_density_kg_m }
    /// Layer axial stiffness [N].
    #[must_use]
    pub const fn axial_stiffness_n(&self) -> f64 { self.axial_stiffness_n }
    /// Layer bending stiffness [N m2].
    #[must_use]
    pub const fn bending_stiffness_n_m2(&self) -> f64 { self.bending_stiffness_n_m2 }
}

/// Immutable section integrals with all original source bundles and geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricElasticSection {
    inner_radius_m: f64,
    layers: Vec<ConcentricMaterialLayer>,
    layer_integrals: Vec<ElasticLayerIntegrals>,
    area_m2: f64,
    second_moment_m4: f64,
    line_density_kg_m: f64,
    axial_stiffness_n: f64,
    bending_stiffness_n_m2: f64,
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}
fn positive(value: f64, quantity: &'static str) -> Result<f64, MaterialStatePointError> {
    if value.is_finite() && value > 0.0 { Ok(value) } else { Err(invalid(quantity)) }
}

impl ConcentricElasticSection {
    /// Integrate a solid rod (`inner_radius_m = 0`) or a hollow layered tube.
    /// Layer order is spatial and MUST NOT be sorted by material identity.
    ///
    /// # Errors
    /// Refuses absent/wrong-quantity/nonpositive properties, empty or oversized
    /// sections, reversed/zero-thickness annuli and unrepresentable arithmetic.
    /// No partial section or substituted property is returned.
    pub fn try_new(
        inner_radius_m: f64,
        layers: Vec<ConcentricMaterialLayer>,
    ) -> Result<Self, MaterialStatePointError> {
        if !inner_radius_m.is_finite() || inner_radius_m < 0.0 {
            return Err(invalid("finite nonnegative section bore radius"));
        }
        if layers.is_empty() || layers.len() > MAX_CONCENTRIC_LAYERS {
            return Err(invalid("section requires 1..=128 layers"));
        }
        let mut inner = inner_radius_m;
        let mut area = 0.0;
        let mut inertia = 0.0;
        let mut line_density = 0.0;
        let mut axial = 0.0;
        let mut bending = 0.0;
        let mut layer_integrals = Vec::with_capacity(layers.len());
        for layer in &layers {
            let outer = layer.outer_radius_m;
            if !outer.is_finite() || outer <= inner {
                return Err(invalid("section radii must be finite and strictly increasing"));
            }
            let property = |name, quantity| {
                let resolved = layer.material.property(name)
                    .ok_or_else(|| invalid("section layer is missing a required property"))?;
                if resolved.requirement().quantity() != quantity {
                    return Err(invalid("section layer property has the wrong quantity"));
                }
                positive(resolved.value_si(), "positive finite section layer property")
            };
            let rho = property(DENSITY_PROPERTY, QuantitySpec::dimensional(Density::DIMS))?;
            let young = property(YOUNG_MODULUS_PROPERTY, QuantitySpec::dimensional(Pressure::DIMS))?;
            // Factored annular differences avoid cancellation of nearly equal
            // squared/fourth powers when a thin layer is still representable.
            let a = positive(core::f64::consts::PI * (outer - inner) * (outer + inner), "annular area")?;
            let i = positive(0.25 * a * (outer * outer + inner * inner), "annular second moment")?;
            let layer_mass = positive(rho * a, "layer line density")?;
            let layer_axial = positive(young * a, "layer axial stiffness")?;
            let layer_bending = positive(young * i, "layer bending stiffness")?;
            area = positive(area + a, "section area")?;
            inertia = positive(inertia + i, "section second moment")?;
            line_density = positive(line_density + layer_mass, "section line density")?;
            axial = positive(axial + layer_axial, "section axial stiffness")?;
            bending = positive(bending + layer_bending, "section bending stiffness")?;
            layer_integrals.push(ElasticLayerIntegrals {
                area_m2: a, second_moment_m4: i, line_density_kg_m: layer_mass,
                axial_stiffness_n: layer_axial, bending_stiffness_n_m2: layer_bending,
            });
            inner = outer;
        }
        Ok(Self { inner_radius_m, layers, layer_integrals, area_m2: area, second_moment_m4: inertia,
            line_density_kg_m: line_density, axial_stiffness_n: axial, bending_stiffness_n_m2: bending })
    }

    /// Original ordered material assignments; no averaged material receipt.
    #[must_use]
    pub fn layers(&self) -> &[ConcentricMaterialLayer] { &self.layers }
    /// Integrals in exactly the same spatial order as `layers()`.
    #[must_use]
    pub fn layer_integrals(&self) -> &[ElasticLayerIntegrals] { &self.layer_integrals }
    /// Explicit central void radius [m].
    #[must_use]
    pub const fn inner_radius_m(&self) -> f64 { self.inner_radius_m }
    /// Physical outer radius [m], not an equivalent stiffness radius.
    #[must_use]
    pub fn outer_radius_m(&self) -> f64 { self.layers[self.layers.len() - 1].outer_radius_m }
    /// Solid cross-sectional area [m2].
    #[must_use]
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Geometric second moment about either centroidal diameter [m4].
    #[must_use]
    pub const fn second_moment_m4(&self) -> f64 { self.second_moment_m4 }
    /// Mass per physical length [kg/m].
    #[must_use]
    pub const fn line_density_kg_m(&self) -> f64 { self.line_density_kg_m }
    /// Parallel-strain axial stiffness [N].
    #[must_use]
    pub const fn axial_stiffness_n(&self) -> f64 { self.axial_stiffness_n }
    /// Common-curvature flexural stiffness [N m2].
    #[must_use]
    pub const fn bending_stiffness_n_m2(&self) -> f64 { self.bending_stiffness_n_m2 }
}

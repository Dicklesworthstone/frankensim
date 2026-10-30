//! Source-bound sections for slender, perfectly bonded elastic members.
//!
//! The concentric geometry has a common centroid and neutral axis, so the
//! one-dimensional parallel-strain model gives `mass/length = sum(rho A)`,
//! `EA = sum(E A)` and `EI = sum(E I)`. Every layer keeps its own material
//! state; there is no fictitious homogeneous material or averaged density.
//!
//! These are section integrals, not a three-dimensional composite law. The
//! optional shear/torsion binding adds a uniform-shear rod reduction with an
//! explicit correction factor; Kelvin-Voigt bindings add sourced bending and
//! linear axial loss. [`winding`] and the mass-only [`wound`] add helical
//! wraps, and [`strength`] checks axial stress against sourced limits.
//! Winding slip and damage are not modeled.
//! The explicit [`thermal`] binding adds per-layer thermal eigenstrain and
//! axial stress. Source validity applies at retained states, not future states.

/// Explicit uniaxial viscosity and longitudinal modal loss.
pub mod axial;
/// Per-layer sourced loss tangents as hysteretic bending loss.
pub mod loss_angle;
/// Per-layer generalized-Maxwell bending spectra without averaged time constants.
pub mod relaxation;
/// Axial tensile capacity from per-layer sourced strength limits.
pub mod strength;
/// Computed Zener thermoelastic bending loss of homogeneous solid members.
pub mod ted;
pub mod thermal;
pub mod transport;
/// Close-coiled, sliding round-wire winding over an independently sourced core.
pub mod winding;
/// Helical wraps that add inertia, but no tension or stiffness, to a core.
pub mod wound;
/// Rotary inertia and explicitly sourced torsional shear loss.
pub mod torsion;

use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec};
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

    /// Bind an explicitly named shear-modulus claim [Pa] from EVERY layer's
    /// existing state bundle. `GA = correction * sum(G A)` is a uniform-shear
    /// rod reduction; `GJ = sum(G * 2I)` is circular Saint-Venant torsion.
    /// The correction is authored geometry/model data, not a material property.
    /// No Poisson ratio, shear correction, or material-name preset is inferred.
    ///
    /// # Errors
    /// Missing, ambiguous, nonpositive or wrong-quantity shear data, correction
    /// outside (0,1], and unrepresentable products/sums refuse atomically.
    pub fn with_shear_modulus(
        self,
        property_name: &str,
        shear_correction: f64,
    ) -> Result<ConcentricShearSection, MaterialStatePointError> {
        positive(shear_correction, "finite positive shear correction")?;
        if shear_correction > 1.0 { return Err(invalid("shear correction must not exceed one")); }
        let mut shear = 0.0;
        let mut torsion = 0.0;
        for (layer, geometry) in self.layers.iter().zip(&self.layer_integrals) {
            let property = layer.material.property(property_name)
                .ok_or_else(|| invalid("every layer needs an unambiguous shear modulus"))?;
            if property.requirement().quantity() != QuantitySpec::dimensional(Pressure::DIMS) {
                return Err(invalid("section shear modulus must be dimension-only Pa"));
            }
            let modulus = positive(property.value_si(), "positive finite layer shear modulus")?;
            let ga = positive(modulus * geometry.area_m2, "layer shear stiffness")?;
            let gj = positive(2.0 * (modulus * geometry.second_moment_m4), "layer torsional stiffness")?;
            shear = positive(shear + ga, "summed shear stiffness")?;
            torsion = positive(torsion + gj, "summed torsional stiffness")?;
        }
        let shear_stiffness_n = positive(shear_correction * shear, "corrected shear stiffness")?;
        Ok(ConcentricShearSection {
            elastic: self, shear_property: property_name.into(), shear_correction,
            shear_stiffness_n, torsional_stiffness_n_m2: torsion,
        })
    }

    /// Bind a named solid bending-viscosity claim `eta` [Pa s] in each layer.
    /// The source author must identify this as bending strain-rate viscosity;
    /// equal dimensions do not justify selecting fluid shear viscosity.
    /// `C_b = sum(eta I)` acts on curvature rate without averaging moduli or
    /// assigning one layer's decay to its neighbors. Explicit zero eta is allowed;
    /// a missing datum is NOT treated as a lossless layer.
    ///
    /// Density, E and eta must be constants at the held typed environment.
    /// Each eta declares a frequency band; all coefficient bands
    /// intersect. The consumer must admit every actual mode against this band.
    /// This is frozen-state linear bending loss, not axial creep or thermal history.
    ///
    /// # Errors
    /// Refuses missing/ambiguous/wrong-quantity data, negative eta, non-broadband
    /// curves, incompatible frequency support and unrepresentable integrals.
    pub fn with_kelvin_voigt_bending(
        self,
        viscosity_property: &str,
    ) -> Result<ConcentricKelvinVoigtSection, MaterialStatePointError> {
        let mut band = (0.0_f64, f64::MAX);
        let mut viscous_stiffness_n_m2_s = 0.0;
        let mut layer_viscous_stiffness_n_m2_s = Vec::with_capacity(self.layers.len());
        for (layer, geometry) in self.layers.iter().zip(&self.layer_integrals) {
            let viscosity = layer.material.property(viscosity_property)
                .ok_or_else(|| invalid("every layer needs an unambiguous bending viscosity"))?;
            let eta = viscosity.value_si();
            if viscosity.requirement().quantity() != QuantitySpec::dimensional(DynViscosity::DIMS)
                || !eta.is_finite() || eta < 0.0
            {
                return Err(invalid("layer bending viscosity must be nonnegative dimension-only Pa s"));
            }
            let eta_band = viscosity.angular_frequency_band_rad_s()
                .map_err(|_| invalid("invalid bending-viscosity frequency semantics"))?
                .ok_or_else(|| invalid("each bending viscosity needs an explicit frequency validity band"))?;
            band = (band.0.max(eta_band.0), band.1.min(eta_band.1));
            for key in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, viscosity_property] {
                let coefficient = layer.material.property(key)
                    .ok_or_else(|| invalid("layer bending law needs density, modulus and viscosity"))?;
                if !coefficient.is_constant_at_fixed_environment() {
                    return Err(invalid("layer bending coefficients need a held, explicitly typed environment and no unowned response dependence"));
                }
                if let Some((lo, hi)) = coefficient.angular_frequency_band_rad_s()
                    .map_err(|_| invalid("invalid layer coefficient frequency semantics"))?
                {
                    band = (band.0.max(lo), band.1.min(hi));
                }
            }
            let value = eta * geometry.second_moment_m4;
            if !value.is_finite() || (eta > 0.0 && value == 0.0) {
                return Err(invalid("layer viscous flexural stiffness is unrepresentable"));
            }
            viscous_stiffness_n_m2_s += value;
            if !viscous_stiffness_n_m2_s.is_finite() {
                return Err(invalid("summed viscous flexural stiffness is unrepresentable"));
            }
            layer_viscous_stiffness_n_m2_s.push(value);
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid("layer bending frequency validity bands do not intersect"));
        }
        Ok(ConcentricKelvinVoigtSection {
            elastic: self, viscosity_property: viscosity_property.into(),
            layer_viscous_stiffness_n_m2_s, viscous_stiffness_n_m2_s,
            omega_band_rad_s: band,
        })
    }
}

/// Source-retaining axial, bending, shear and torsional rod section.
/// This positive diagonal rod law is not a certification of a 3D elastic tensor.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricShearSection {
    elastic: ConcentricElasticSection,
    shear_property: String,
    shear_correction: f64,
    shear_stiffness_n: f64,
    torsional_stiffness_n_m2: f64,
}
impl ConcentricShearSection {
    /// Common geometry and every original material bundle.
    #[must_use]
    pub const fn elastic(&self) -> &ConcentricElasticSection { &self.elastic }
    /// Explicitly selected shear-modulus property name.
    #[must_use]
    pub fn shear_property(&self) -> &str { &self.shear_property }
    /// Authored uniform-shear correction.
    #[must_use]
    pub const fn shear_correction(&self) -> f64 { self.shear_correction }
    /// Corrected shear stiffness [N].
    #[must_use]
    pub const fn shear_stiffness_n(&self) -> f64 { self.shear_stiffness_n }
    /// Circular torsional stiffness [N m2].
    #[must_use]
    pub const fn torsional_stiffness_n_m2(&self) -> f64 { self.torsional_stiffness_n_m2 }
}

/// A Kelvin-Voigt bending law over the SAME source-bound elastic section.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricKelvinVoigtSection {
    elastic: ConcentricElasticSection,
    viscosity_property: String,
    layer_viscous_stiffness_n_m2_s: Vec<f64>,
    viscous_stiffness_n_m2_s: f64,
    omega_band_rad_s: (f64, f64),
}
impl ConcentricKelvinVoigtSection {
    /// Common geometry, elastic integrals and every original source bundle.
    #[must_use]
    pub const fn elastic(&self) -> &ConcentricElasticSection { &self.elastic }
    /// Caller-selected solid bending-viscosity property.
    #[must_use]
    pub fn viscosity_property(&self) -> &str { &self.viscosity_property }
    /// Per-layer `eta I` in spatial order [N m2 s].
    #[must_use]
    pub fn layer_viscous_stiffness_n_m2_s(&self) -> &[f64] { &self.layer_viscous_stiffness_n_m2_s }
    /// Summed `eta I` multiplying curvature rate [N m2 s].
    #[must_use]
    pub const fn viscous_stiffness_n_m2_s(&self) -> f64 { self.viscous_stiffness_n_m2_s }
    /// Intersection of all source-supported angular-frequency intervals [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }
}

#[cfg(test)]
mod constitutive_tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use crate::state_point::{MaterialPropertySelection, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point};

    fn source(eta: f64, shear: f64, band: (f64, f64), curve: bool) -> ResolvedMaterialStatePoint {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, value, dims) in [
            (DENSITY_PROPERTY, 8000.0, Density::DIMS),
            (YOUNG_MODULUS_PROPERTY, 2e11, Pressure::DIMS),
            ("shear_modulus", shear, Pressure::DIMS),
            ("bending_eta", eta, DynViscosity::DIMS),
        ] {
            let varying = curve && name == "bending_eta";
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varying { PropertyValue::Curve {
                    abscissa: "omega".into(), abscissa_dims: fs_qty::Dims([0, 0, -1, 0, 0, 0]),
                    knots: vec![(band.0, value), (band.1, value * 2.0)], dims,
                } } else { PropertyValue::Scalar { value, dims } },
                validity: ValidityDomain::unconstrained().with("omega", band.0, band.1),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if varying { InterpolationPolicy::LinearInside }
                    else { InterpolationPolicy::ConstantWithinValidity },
                observations: vec![],
                provenance: Provenance { source: "synthetic section law".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: "synthetic".into(), phase: "solid".into(), process: "test".into(), revision: 0,
        }, claims, vec![]).unwrap();
        resolve_material_state_point(&card, &QueryPoint::new().with("omega", (band.0 + band.1) / 2.0).unwrap(),
            &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
    }
    fn section(a: ResolvedMaterialStatePoint, b: ResolvedMaterialStatePoint) -> ConcentricElasticSection {
        ConcentricElasticSection::try_new(0.0, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.001, material: a },
            ConcentricMaterialLayer { outer_radius_m: 0.002, material: b },
        ]).unwrap()
    }
    fn near(a: f64, b: f64) { assert!((a-b).abs() < 2e-13 * a.abs().max(b.abs()), "{a} != {b}"); }

    fn moisture_source(at: f64, typed: bool, extra_strain: bool) -> ResolvedMaterialStatePoint {
        use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};
        let moisture = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static));
        let mut claims = ClaimSet::new();
        let mut required = Vec::new();
        for (name, dims, low, high) in [
            (DENSITY_PROPERTY, Density::DIMS, 8000.0, 9000.0),
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 2e11, 1e11),
            ("shear_modulus", Pressure::DIMS, 8e10, 4e10),
            ("bending_eta", DynViscosity::DIMS, 1e7, 2e7),
            ("delta_e", Pressure::DIMS, 1e9, 2e9),
            ("tau", fs_qty::Time::DIMS, 0.001, 0.002),
        ] {
            let mut validity = ValidityDomain::unconstrained().with("omega", 1.0, 20000.0);
            validity = if typed { validity.with_quantity("moisture_content", moisture, 0.0, 0.2) }
                else { validity.with("moisture_content", 0.0, 0.2) };
            if extra_strain { validity = validity.with("strain", 0.0, 0.1); }
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve { abscissa: "moisture_content".into(),
                    abscissa_dims: fs_qty::Dims::NONE, knots: vec![(0.0, low), (0.2, high)], dims },
                validity, uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside, observations: vec![],
                provenance: Provenance { source: "manufactured held-moisture section".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            required.push(ScalarPropertyRequirement::try_new(name, dims,
                ScalarAdmissibility::StrictlyPositive).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: "manufactured".into(),
            phase: "solid".into(), process: "held environment".into(), revision: 0 }, claims, vec![]).unwrap();
        let point = QueryPoint::new().with("omega", 1000.0).unwrap();
        let mut point = if typed { point.with_quantity("moisture_content", moisture, at).unwrap() }
            else { point.with("moisture_content", at).unwrap() };
        if extra_strain { point = point.with("strain", 0.01).unwrap(); }
        resolve_material_state_point(&card, &point, &required,
            MaterialPropertySelection::SingleClaimOnly).unwrap()
    }

    #[test]
    fn held_moisture_reaches_section_mass_elasticity_and_all_retained_loss_channels() {
        let dry = moisture_source(0.0, true, false);
        let wet = moisture_source(0.2, true, false);
        assert_eq!(dry.card_identity(), wet.card_identity());
        assert_ne!(dry.identity(), wet.identity());
        let make = |material| ConcentricElasticSection::try_new(0.0,
            vec![ConcentricMaterialLayer { outer_radius_m: 0.001, material }]).unwrap();
        let a = make(dry.clone());
        let b = make(wet.clone());
        near(b.line_density_kg_m() / a.line_density_kg_m(), 9.0 / 8.0);
        near(b.axial_stiffness_n() / a.axial_stiffness_n(), 0.5);
        near(b.bending_stiffness_n_m2() / a.bending_stiffness_n_m2(), 0.5);
        let loss_a = a.clone().with_kelvin_voigt_bending("bending_eta").unwrap();
        let loss_b = b.clone().with_kelvin_voigt_bending("bending_eta").unwrap();
        near(loss_b.viscous_stiffness_n_m2_s() / loss_a.viscous_stiffness_n_m2_s(), 2.0);
        assert_eq!(loss_b.omega_band_rad_s(), (1.0, 20000.0));
        assert_eq!(loss_b.elastic().layers()[0].material, wet);
        let axial_a = a.clone().with_axial_viscosity("bending_eta").unwrap();
        let axial_b = b.clone().with_axial_viscosity("bending_eta").unwrap();
        near(axial_b.viscous_rigidity_n_s() / axial_a.viscous_rigidity_n_s(), 2.0);
        let torsion_a = a.clone().with_shear_modulus("shear_modulus", 1.0).unwrap()
            .with_torsional_viscosity("bending_eta").unwrap();
        let torsion_b = b.clone().with_shear_modulus("shear_modulus", 1.0).unwrap()
            .with_torsional_viscosity("bending_eta").unwrap();
        near(torsion_b.stiffness_n_m2() / torsion_a.stiffness_n_m2(), 0.5);
        near(torsion_b.rotary_inertia_kg_m() / torsion_a.rotary_inertia_kg_m(), 9.0 / 8.0);
        let spectrum = [relaxation::LayerRelaxationProperties {
            equilibrium_modulus: YOUNG_MODULUS_PROPERTY.into(),
            branches: vec![("delta_e".into(), "tau".into())],
        }];
        let prony = b.with_prony_bending(&spectrum).unwrap();
        assert_eq!(prony.branches()[0].relaxation_time_s, 0.002);
        near(prony.branches()[0].relaxing_stiffness_n_m2,
            2e9 * prony.elastic().second_moment_m4());
        let property = wet.property(YOUNG_MODULUS_PROPERTY).unwrap();
        assert!(property.is_constant_at_fixed_environment());
        assert!(!property.is_static_at_fixed_external_conditions(), "frequency applicability stays explicit");
        assert!(!property.is_constant_at_fixed_temperature(), "the old predicate is unchanged");
        assert!(matches!(property.answer().receipt.decision, fs_matdb::EvaluationDecision::ExactTabulated { .. }));
    }

    #[test]
    fn held_environment_never_freezes_untyped_moisture_unowned_strain_or_frequency_response() {
        for material in [moisture_source(0.1, false, false), moisture_source(0.1, true, true),
            source(1e7, 8e10, (1.0, 20000.0), true)] {
            assert!(!material.property("bending_eta").unwrap().is_constant_at_fixed_environment());
            assert!(section(material.clone(), material).with_kelvin_voigt_bending("bending_eta").is_err());
        }
    }

    #[test]
    fn layer_integrals_drive_shear_torsion_and_bending_loss_without_averaging() {
        let a = source(1e7, 8e10, (1.0, 20000.0), false);
        let b = source(4e7, 4e10, (10.0, 10000.0), false);
        let elastic = section(a.clone(), b);
        let g = elastic.layer_integrals();
        let shear = elastic.clone().with_shear_modulus("shear_modulus", 0.8).unwrap();
        near(shear.shear_stiffness_n(), 0.8 * (8e10*g[0].area_m2() + 4e10*g[1].area_m2()));
        near(shear.torsional_stiffness_n_m2(), 2.0 * (8e10*g[0].second_moment_m4() + 4e10*g[1].second_moment_m4()));
        let loss = elastic.clone().with_kelvin_voigt_bending("bending_eta").unwrap();
        near(loss.viscous_stiffness_n_m2_s(), 1e7*g[0].second_moment_m4() + 4e7*g[1].second_moment_m4());
        assert_eq!(loss.omega_band_rad_s(), (10.0, 10000.0));
        assert_eq!(loss.elastic(), &elastic);
        assert_eq!(shear.elastic().layers()[0].material, a);
    }

    #[test]
    fn homogeneous_subdivision_is_invariant_for_hollow_and_solid_sections() {
        let material = source(1e7, 8e10, (1.0, 20000.0), false);
        for bore in [0.0, 0.0005] {
            let make = |radii: &[f64]| ConcentricElasticSection::try_new(bore,
                radii.iter().map(|&r| ConcentricMaterialLayer { outer_radius_m: r, material: material.clone() }).collect()).unwrap();
            let whole = make(&[0.002]);
            let split = make(&[0.001, 0.002]);
            let a = whole.clone().with_shear_modulus("shear_modulus", 0.9).unwrap();
            let b = split.clone().with_shear_modulus("shear_modulus", 0.9).unwrap();
            near(a.shear_stiffness_n(), b.shear_stiffness_n());
            near(a.torsional_stiffness_n_m2(), b.torsional_stiffness_n_m2());
            near(whole.with_kelvin_voigt_bending("bending_eta").unwrap().viscous_stiffness_n_m2_s(),
                split.with_kelvin_voigt_bending("bending_eta").unwrap().viscous_stiffness_n_m2_s());
        }
    }

    #[test]
    fn explicit_zero_loss_is_admitted_but_missing_and_negative_loss_refuse() {
        let zero = source(0.0, 8e10, (1.0, 20000.0), false);
        let elastic = section(zero.clone(), zero);
        assert_eq!(elastic.clone().with_kelvin_voigt_bending("bending_eta").unwrap().viscous_stiffness_n_m2_s(), 0.0);
        assert!(elastic.with_kelvin_voigt_bending("missing").is_err());
        let negative = source(-1.0, 8e10, (1.0, 20000.0), false);
        assert!(section(negative.clone(), negative).with_kelvin_voigt_bending("bending_eta").is_err());
    }

    #[test]
    fn sampled_frequency_laws_and_disjoint_layer_bands_are_not_broadband_models() {
        let a = source(1e7, 8e10, (1.0, 100.0), false);
        let b = source(1e7, 8e10, (101.0, 200.0), false);
        assert!(section(a, b).with_kelvin_voigt_bending("bending_eta").is_err());
        let sampled = source(1e7, 8e10, (1.0, 20000.0), true);
        assert!(section(sampled.clone(), sampled).with_kelvin_voigt_bending("bending_eta").is_err());
    }

    #[test]
    fn invalid_shear_inputs_do_not_fall_back_to_poisson_or_material_names() {
        let a = source(1e7, 8e10, (1.0, 20000.0), false);
        let elastic = section(a.clone(), a);
        for correction in [0.0, -1.0, 1.01, f64::NAN, f64::INFINITY] {
            assert!(elastic.clone().with_shear_modulus("shear_modulus", correction).is_err());
        }
        assert!(elastic.clone().with_shear_modulus("missing", 1.0).is_err());
        assert!(elastic.with_shear_modulus("density", 1.0).is_err());
        let negative = source(1e7, -1.0, (1.0, 20000.0), false);
        assert!(section(negative.clone(), negative).with_shear_modulus("shear_modulus", 1.0).is_err());
    }
}
